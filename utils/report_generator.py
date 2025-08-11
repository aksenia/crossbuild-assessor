#!/usr/bin/env python3
"""
HTML Report Generator - Clinical Evidence Focus

Creates a comprehensive HTML report from CrossBuild Assessor outputs.
Prefers structured JSON data, falls back to text parsing for backward compatibility.
"""

import argparse
import base64
import json
import pandas as pd
from datetime import datetime
from pathlib import Path
from jinja2 import Template
from data_utils import format_consequence_relationship


class ReportGenerator:
    """Generate HTML reports from CrossBuild Assessor outputs"""
    
    def __init__(self, input_dir):
        """Initialize with input directory containing analysis outputs"""
        self.input_dir = Path(input_dir)
        self.report_data = {}

    def _debug_summary_data(self):
        """Debug what's actually in the summary data"""
        print("\n=== DEBUG: Summary data structure ===")
        
        if 'summaries' in self.report_data:
            for section_name, section_data in self.report_data['summaries'].items():
                print(f"\nSection: {section_name}")
                print(f"Type: {type(section_data)}")
                
                if isinstance(section_data, dict):
                    for key, value in section_data.items():
                        print(f"  {key}: {type(value)} = {value if len(str(value)) < 100 else str(value)[:100] + '...'}")
                else:
                    print(f"  Content: {section_data}")
        else:
            print("No 'summaries' key found in report_data")
        
    def generate_report(self, output_file):
        """Generate complete HTML report"""
        print("Generating HTML report...")
        
        # Collect all data
        self._collect_metadata()
        self._collect_summary_data()
        self._debug_summary_data() # DEBUG
        self._collect_plot_images()
        self._collect_variant_data()
        
        # Generate HTML
        html_content = self._render_html()
        
        # Save report
        with open(output_file, 'w', encoding='utf-8') as f:
            f.write(html_content)
        
        print(f"✓ HTML report saved to: {output_file}")
        return output_file
    
    def _format_consequence_relationship(self, variant):
        """Format consequence relationship with unified display format"""
        relationship = variant.get('Consequence_Relationship', 'unknown')
        change = variant.get('Consequence_Change', 'no data')
        
        return f"{relationship}: {change}"
    
    def _collect_metadata(self):
        """Collect basic metadata"""
        self.report_data['timestamp'] = datetime.now().strftime('%Y-%m-%d %H:%M:%S')
        self.report_data['input_dir'] = str(self.input_dir)
        
    def _collect_summary_data(self):
        """Load structured JSON data if available, fallback to text parsing"""
        summaries = {}
        
        # Try to load JSON data first (preferred)
        liftover_json = self.input_dir / 'liftover_analysis.json'
        if liftover_json.exists():
            print("Loading liftover data from JSON...")
            with open(liftover_json) as f:
                summaries['liftover'] = json.load(f)
        else:
            # Fallback to text parsing
            print("JSON not found, parsing liftover text summary...")
            liftover_summary = self.input_dir / 'liftover_analysis_summary.txt'
            if liftover_summary.exists():
                summaries['liftover'] = self._parse_liftover_summary_text(liftover_summary)
        
        # Try to load prioritization JSON
        priority_json = self.input_dir / 'prioritization_results.json'
        if priority_json.exists():
            print("Loading prioritization data from JSON...")
            with open(priority_json) as f:
                summaries['prioritization'] = json.load(f)
        else:
            # Fallback to text parsing
            print("JSON not found, parsing prioritization text summary...")
            priority_summary = self.input_dir / 'variant_prioritization_summary.txt'
            if priority_summary.exists():
                summaries['prioritization'] = self._parse_prioritization_summary_text(priority_summary)
        
        self.report_data['summaries'] = summaries
    
    def _parse_liftover_summary_text(self, file_path):
        """Extract key metrics from liftover summary (fallback)"""
        with open(file_path, 'r') as f:
            content = f.read()
        
        # Extract key numbers (basic parsing)
        data = {}
        lines = content.split('\n')
        
        for line in lines:
            if 'Total variants analyzed:' in line:
                data['total_variants'] = line.split(':')[1].strip()
            elif 'Position match rate:' in line:
                data['position_match_rate'] = line.split(':')[1].strip()
            elif 'Genotype match rate:' in line:
                data['genotype_match_rate'] = line.split(':')[1].strip()
                
        return {'dataset_overview': data}
    
    def _parse_prioritization_summary_text(self, file_path):
        """Extract key metrics from prioritization summary (fallback)"""
        with open(file_path, 'r') as f:
            content = f.read()
        
        data = {}
        lines = content.split('\n')
        
        for line in lines:
            if 'Total discordant variants analyzed:' in line:
                data['total_discordant'] = line.split(':')[1].strip()
            elif 'Variants included in Excel output:' in line:
                data['excel_output'] = line.split(':')[1].strip()
            elif line.strip().startswith('CRITICAL:'):
                data['critical_count'] = line.split(':')[1].split('(')[0].strip()
            elif line.strip().startswith('HIGH:'):
                data['high_count'] = line.split(':')[1].split('(')[0].strip()
                
        return {'dataset_overview': data}
    
    def _collect_plot_images(self):
        """Convert plot images to base64 for embedding"""
        images = {}
        
        plot_files = {
            'liftover_analysis': 'liftover_analysis.png',
            'position_differences': 'position_differences_analysis.png', 
            'prioritization_plots': 'variant_prioritization_plots.png'
        }
        
        for key, filename in plot_files.items():
            plot_path = self.input_dir / filename
            if plot_path.exists():
                with open(plot_path, 'rb') as f:
                    encoded = base64.b64encode(f.read()).decode('utf-8')
                    images[key] = f"data:image/png;base64,{encoded}"
        
        self.report_data['images'] = images
    
    def _collect_variant_data(self):
        """Load variant data for clinical evidence table - FIXED: Show only variants with actual changes"""
        csv_file = self.input_dir / 'prioritized_variants.csv'
        if csv_file.exists():
            df = pd.read_csv(csv_file)
            
            # Debug: Print available columns to help identify the correct column name
            print(f"Available columns in CSV: {list(df.columns)}")
            
            # Top 10 variants with clinical evidence focus
            if len(df) > 0:
                # Filter to variants with meaningful clinical changes only
                has_changes_mask = pd.Series([False] * len(df))

                # Clinical significance changes
                if 'Has_Clinical_Change' in df.columns:
                    has_changes_mask |= (df['Has_Clinical_Change'] == 'YES')

                # Consequence relationship changes  
                if 'Has_Consequence_Change' in df.columns:
                    has_changes_mask |= (df['Has_Consequence_Change'] == 'YES')

                # HGVS discordance (replacing Has_Impact_Change)
                if 'HGVSc_MATCHED_discordant' in df.columns:
                    def has_hgvs_discordance(value):
                        """Check if there are discordant HGVS matches"""
                        if pd.isna(value) or str(value).strip() in ['', '-', 'nan']:
                            return False
                        try:
                            # Count non-empty items in comma-separated list
                            items = [item.strip() for item in str(value).split(',') if item.strip()]
                            return len(items) > 0
                        except:
                            return False
                    
                    has_changes_mask |= df['HGVSc_MATCHED_discordant'].apply(has_hgvs_discordance)
                
                # Filter to variants with changes
                df_with_changes = df[has_changes_mask]
                print(f"Filtered to {len(df_with_changes)} variants with actual changes (from {len(df)} total)")
                
                # Sort by Priority_Score descending to get highest priority first
                if 'Priority_Score' in df_with_changes.columns and len(df_with_changes) > 0:
                    df_sorted = df_with_changes.sort_values('Priority_Score', ascending=False)
                elif 'Rank' in df_with_changes.columns and len(df_with_changes) > 0:
                    df_sorted = df_with_changes.sort_values('Rank', ascending=True)  # Lower rank = higher priority
                else:
                    df_sorted = df_with_changes  # Use original order if no sorting column available
                
                # Select columns for clinical review
                clinical_columns = [
                    'Rank', 'Chromosome', 'Position_hg19', 'Gene_hg19', 'Gene_hg38',
                    'Clinical_Significance_hg19', 'Clinical_Significance_hg38',
                    'Impact_hg19', 'Impact_hg38',
                    'SIFT_hg19', 'SIFT_hg38', 'PolyPhen_hg19', 'PolyPhen_hg38',
                    'Priority_Score', 'Priority_Category',
                    'Consequence_Relationship', 'Consequence_Change',
                    'High_Impact_Consequences_hg19', 'High_Impact_Consequences_hg38',
                    # HGVS columns
                    'CANONICAL_HGVSc_Match', 'HGVSc_MATCHED_concordant', 'HGVSc_MATCHED_discordant'
                ]
                
                # Only include columns that exist
                available_columns = [col for col in clinical_columns if col in df_sorted.columns]
                print(f"Using columns: {available_columns}")
                
                # Take top 10 highest priority variants with changes
                if len(df_sorted) > 0:
                    top_variants = df_sorted[available_columns].head(10)
                    self.report_data['top_variants'] = top_variants.to_dict('records')
                else:
                    print("No variants with changes found")
                    self.report_data['top_variants'] = []
            else:
                self.report_data['top_variants'] = []
    
    def _get_summary_value(self, data_path, fallback='N/A'):
        """Safely extract values from nested summary data with detailed debugging"""
        try:
            current = self.report_data['summaries']
            for i, key in enumerate(data_path):
                if key not in current:
                    print(f"DEBUG: Key '{key}' not found at path {data_path[:i+1]}")
                    return fallback
                current = current[key]
            
            # Debug what we're actually returning
            if data_path == ['liftover', 'flip_swap_analysis'] or data_path == ['prioritization', 'priority_distribution']:
                print(f"DEBUG: Path {data_path} returned type {type(current)}: {current}")
            
            return current
        except (KeyError, TypeError) as e:
            print(f"DEBUG: Error accessing {data_path}: {e}")
            return fallback
    
    def _render_html(self):
        """Render HTML using Jinja2 template with clinical focus"""
        template_str = '''<!DOCTYPE html>
<html lang="en">
<head>
    <meta charset="UTF-8">
    <meta name="viewport" content="width=device-width, initial-scale=1.0">
    <title>CrossBuild Assessor Report</title>
    <style>
        body { font-family: 'Segoe UI', Tahoma, Geneva, Verdana, sans-serif; margin: 0; padding: 20px; background: #f5f5f5; }
        .container { max-width: 1200px; margin: 0 auto; background: white; padding: 30px; border-radius: 8px; box-shadow: 0 2px 10px rgba(0,0,0,0.1); }
        .header { border-bottom: 3px solid #2c5f8a; padding-bottom: 20px; margin-bottom: 30px; }
        .header h1 { color: #2c5f8a; margin: 0; font-size: 28px; }
        .header .meta { color: #666; margin-top: 5px; }
        .section { margin-bottom: 40px; }
        .section h2 { color: #2c5f8a; border-bottom: 2px solid #e0e0e0; padding-bottom: 10px; margin-bottom: 20px; }
        .metrics-grid { display: grid; grid-template-columns: repeat(auto-fit, minmax(200px, 1fr)); gap: 15px; margin-bottom: 20px; }
        .metric-card { background: #f8f9fa; padding: 15px; border-radius: 6px; border-left: 4px solid #2c5f8a; }
        .metric-value { font-size: 24px; font-weight: bold; color: #2c5f8a; }
        .metric-label { color: #666; font-size: 14px; }
        .plot-container { text-align: center; margin: 20px 0; }
        .plot-container img { max-width: 100%; height: auto; border: 1px solid #ddd; border-radius: 4px; }
        .table-container { overflow-x: auto; }
        table { width: 100%; border-collapse: collapse; margin-top: 15px; }
        th, td { padding: 8px 12px; text-align: left; border-bottom: 1px solid #ddd; font-size: 14px; }
        th { background: #f8f9fa; font-weight: 600; color: #2c5f8a; }
        .clinical-change { background: #ffebee; color: #c62828; font-weight: bold; }
        .clinical-stable { background: #e8f5e8; color: #2e7d32; }
        .no-data { color: #999; font-style: italic; text-align: center; padding: 20px; }
        .summary-text { background: #f8f9fa; padding: 15px; border-radius: 6px; margin: 15px 0; }
        @media print { body { background: white; } .container { box-shadow: none; } }
    </style>
</head>
<body>
    <div class="container">
        <div class="header">
            <h1>CrossBuild Assessor Report</h1>
            <div class="meta">Generated: {{ timestamp }} | Analysis: Genome build comparison (hg19 ↔ hg38)</div>
        </div>

        <!-- EXECUTIVE SUMMARY -->
        <div class="section">
            <h2>Executive Summary</h2>
            <div class="metrics-grid">
                <div class="metric-card">
                    <div class="metric-value">{{ get_summary_value(['liftover', 'dataset_overview', 'total_variants']) }}</div>
                    <div class="metric-label">Total variants analyzed</div>
                </div>
                <div class="metric-card">
                    <div class="metric-value">{{ get_summary_value(['liftover', 'dataset_overview', 'position_match_percentage']) }}%</div>
                    <div class="metric-label">Position match rate</div>
                </div>
                <div class="metric-card">
                    <div class="metric-value">{{ get_summary_value(['liftover', 'dataset_overview', 'genotype_match_percentage']) }}%</div>
                    <div class="metric-label">Genotype match rate</div>
                </div>
                <div class="metric-card">
                    <div class="metric-value">{{ get_summary_value(['prioritization', 'top_variants_summary', 'critical_variants_distinct']) }}</div>
                    <div class="metric-label">Critical variants</div>
                </div>
                <div class="metric-card">
                    <div class="metric-value">{{ get_summary_value(['prioritization', 'top_variants_summary', 'clinical_changes_count']) }}</div>
                    <div class="metric-label">Clinical significance changes</div>
                </div>
            </div>
            
            <div class="summary-text">
                <strong>Key findings:</strong> 
                {{ get_summary_value(['liftover', 'dataset_overview', 'position_match_percentage']) }}% of {{ get_summary_value(['liftover', 'dataset_overview', 'total_variants']) }} variants show coordinate concordance between liftover tools. 
                {{ get_summary_value(['prioritization', 'top_variants_summary', 'critical_variants_distinct']) }} variants require immediate review due to functional or clinical significance changes.
                {{ get_summary_value(['prioritization', 'top_variants_summary', 'clinical_changes_count']) }} variants show clinical significance transitions that may affect interpretation.
            </div>
        </div>

        <!-- LIFTOVER QUALITY CONTROL -->
        <div class="section">
            <h2>Liftover Quality Control</h2>
            
            <h3>Tool performance summary</h3>
            <div class="metrics-grid">
                <div class="metric-card">
                    <div class="metric-value">{{ get_summary_value(['liftover', 'dataset_overview', 'concordant_variants']) }}</div>
                    <div class="metric-label">Concordant variants</div>
                </div>
                <div class="metric-card">
                    <div class="metric-value">{{ get_summary_value(['liftover', 'dataset_overview', 'discordant_variants']) }}</div>
                    <div class="metric-label">Discordant variants</div>
                </div>
                <div class="metric-card">
                    <div class="metric-value">{{ get_summary_value(['liftover', 'quality_summary', 'position_mismatches']) }}</div>
                    <div class="metric-label">Position mismatches</div>
                </div>
                <div class="metric-card">
                    <div class="metric-value">{{ get_summary_value(['liftover', 'quality_summary', 'genotype_mismatches']) }}</div>
                    <div class="metric-label">Genotype mismatches</div>
                </div>
            </div>

            <!-- STRAND FLIP SWAP ANALYSIS -->
            <h3>Strand flip & allele swap analysis</h3>
            <div class="table-container">
                <table>
                    <thead>
                        <tr>
                            <th>Category</th>
                            <th>Count</th>
                            <th>Description</th>
                        </tr>
                    </thead>
                    <tbody>
                        {% for category, count in get_summary_value(['liftover', 'flip_swap_analysis']).items() %}
                        <tr>
                            <td>{{ category }}</td>
                            <td>{{ count }}</td>
                            <td>
                                {% if category == "No Changes Required" %}
                                    Variants lifted successfully without modifications
                                {% elif category == "Strand Flip Only" %}
                                    Strand orientation corrected during liftover
                                {% elif category == "Allele Swap Successful" %}
                                    REF and ALT alleles were swapped to match reference genome
                                {% elif category == "Strand Flip + Allele Swap Successful" %}
                                    Both strand flip and allele swap occurred successfully
                                {% elif category == "Allele Swap Failed" %}
                                    Allele swap attempted but failed due to ambiguous alleles
                                {% elif category == "Strand Flip + Allele Swap Failed" %}
                                    Strand flip occurred but allele swap failed due to ambiguous alleles
                                {% else %}
                                    {{ category }}
                                {% endif %}
                            </td>
                        </tr>
                        {% endfor %}
                    </tbody>
                </table>
            </div>

            <!-- LIFTOVER PLOTS -->
            {% if images.liftover_analysis %}
            <h3>Quality control visualizations</h3>
            <div class="plot-container">
                <img src="{{ images.liftover_analysis }}" alt="Liftover Analysis">
            </div>
            {% endif %}
            
            {% if images.position_differences %}
            <div class="plot-container">
                <img src="{{ images.position_differences }}" alt="Position Differences Analysis">
            </div>
            {% endif %}
        </div>

        <!-- VARIANT BUILD ANNOTATIONS COMPARISON -->
        <div class="section">
            <h2>Variant build annotations comparison</h2>
            
            <h3>Priority distribution</h3>
            <div class="summary-text">
                <strong>Priority definitions:</strong>
                <strong>CRITICAL:</strong> Clinical interpretation changes or high impact transitions |
                <strong>HIGH:</strong> Functionally significant changes |
                <strong>MODERATE:</strong> Pathogenicity prediction changes |
                <strong>LOW:</strong> Technical issues, annotation differences, and gene symbol changes
            </div>
            <div class="metrics-grid">
                {% for category, count in get_summary_value(['prioritization', 'priority_distribution']).items() %}
                <div class="metric-card">
                    <div class="metric-value">{{ count }}</div>
                    <div class="metric-label">{{ category }} priority</div>
                </div>
                {% endfor %}
            </div>

            <!-- HGVS NOMENCLATURE ANALYSIS (UPDATED) -->
            <h3>HGVS nomenclature analysis</h3>
            <div class="metrics-grid">
                <div class="metric-card">
                    <div class="metric-value">{{ get_summary_value(['prioritization', 'hgvs_analysis', 'canonical_hgvsc_match', 'match_rate_percentage']) }}%</div>
                    <div class="metric-label">Canonical HGVSc Match Rate</div>
                </div>
                <div class="metric-card">
                    <div class="metric-value">{{ get_summary_value(['prioritization', 'hgvs_analysis', 'canonical_hgvsp_match', 'match_rate_percentage']) }}%</div>
                    <div class="metric-label">Canonical HGVSp Match Rate</div>
                </div>
                <div class="metric-card">
                    <div class="metric-value">{{ get_summary_value(['prioritization', 'hgvs_analysis', 'hgvsc_concordance', 'concordant_rate_percentage']) }}%</div>
                    <div class="metric-label">HGVSc Concordance Rate</div>
                </div>
                <div class="metric-card">
                    <div class="metric-value">{{ get_summary_value(['prioritization', 'hgvs_analysis', 'hgvsp_concordance', 'concordant_rate_percentage']) }}%</div>
                    <div class="metric-label">HGVSp Concordance Rate</div>
                </div>
                <div class="metric-card">
                    <div class="metric-value">{{ get_summary_value(['prioritization', 'hgvs_analysis', 'matched_transcripts', 'hgvsc_average_per_variant']) }}</div>
                    <div class="metric-label">Avg HGVSc-Matched Transcripts</div>
                </div>
                <div class="metric-card">
                    <div class="metric-value">{{ get_summary_value(['prioritization', 'hgvs_analysis', 'matched_transcripts', 'hgvsp_average_per_variant']) }}</div>
                    <div class="metric-label">Avg HGVSp-Matched Transcripts</div>
                </div>
            </div>

            <!-- CLINICAL SIGNIFICANCE TRANSITIONS -->
            <h3>Clinical significance transitions (hg19→hg38)</h3>
            <div class="metrics-grid">
                <div class="metric-card">
                    <div class="metric-value">{{ get_summary_value(['prioritization', 'clinical_transitions', 'total_stable']) }}</div>
                    <div class="metric-label">Stable annotations</div>
                </div>
                <div class="metric-card">
                    <div class="metric-value">{{ get_summary_value(['prioritization', 'clinical_transitions', 'total_changing']) }}</div>
                    <div class="metric-label">Changing annotations</div>
                </div>
            </div>
            
            <div class="table-container">
                <table>
                    <thead>
                        <tr>
                            <th>Transition</th>
                            <th>Variant count</th>
                            <th>Impact level</th>
                        </tr>
                    </thead>
                    <tbody>
                        {% set ordered_changes = get_summary_value(['prioritization', 'clinical_transitions', 'directional_changes_ordered']) %}
                        {% set regular_changes = get_summary_value(['prioritization', 'clinical_transitions', 'directional_changes']) %}
                        
                        {% if ordered_changes and ordered_changes != 'N/A' %}
                            {% for item in ordered_changes %}
                            <tr>
                                <td>{{ item.transition }}</td>
                                <td>{{ item.count }}</td>
                                <td>
                                    {% if item.clinical_priority == 'CRITICAL' %}
                                        <span class="clinical-change">CRITICAL</span>
                                    {% elif item.clinical_priority == 'HIGH' %}
                                        <span class="clinical-change">HIGH</span>
                                    {% elif item.clinical_priority == 'MODERATE' %}
                                        <span class="clinical-stable">MODERATE</span>
                                    {% elif item.clinical_priority == 'LOW' %}
                                        <span class="clinical-stable">LOW</span>
                                    {% else %}
                                        <span class="clinical-stable">{{ item.clinical_priority }}</span>
                                    {% endif %}
                                </td>
                            </tr>
                            {% endfor %}
                        {% elif regular_changes and regular_changes != 'N/A' %}
                            {% for transition, data in regular_changes.items() %}
                            <tr>
                                <td>{{ transition }}</td>
                                <td>{{ data.count }}</td>
                                <td>
                                    {% if data.clinical_priority == 'CRITICAL' %}
                                        <span class="clinical-change">CRITICAL</span>
                                    {% elif data.clinical_priority == 'HIGH' %}
                                        <span class="clinical-change">HIGH</span>
                                    {% elif data.clinical_priority == 'MODERATE' %}
                                        <span class="clinical-stable">MODERATE</span>
                                    {% elif data.clinical_priority == 'LOW' %}
                                        <span class="clinical-stable">LOW</span>
                                    {% else %}
                                        <span class="clinical-stable">{{ data.clinical_priority }}</span>
                                    {% endif %}
                                </td>
                            </tr>
                            {% endfor %}
                        {% else %}
                            <tr>
                                <td colspan="3" class="no-data">No clinical significance transitions found</td>
                            </tr>
                        {% endif %}
                    </tbody>
                </table>
            </div>

            <!-- PRIORITIZATION PLOTS -->
            {% if images.prioritization_plots %}
            <h3>Prioritization visualizations</h3>
            <div class="plot-container">
                <img src="{{ images.prioritization_plots }}" alt="Variant Prioritization">
            </div>
            {% endif %}
        </div>

        <!-- TOP PRIORITY VARIANTS -->
        {% if top_variants %}
        <div class="section">
            <h2>Top priority variants</h2>
            <p><em>Showing example variants with discrepancies</em></p>
            <div class="table-container">
                <table>
                    <thead>
                        <tr>
                            <th>Priority</th>
                            <th>Location</th>
                            <th>Gene (hg19/hg38)</th>
                            <th>HGVSc Matches</th>
                            <th>Clinical significance</th>
                            <th>Impact level</th>
                            <th>Consequence Relationship</th>
                        </tr>
                    </thead>
                    <tbody>
                        {% for variant in top_variants %}
                        <tr>
                            <td>
                                {% if variant.get('Priority_Category') %}
                                    <span class="{% if variant.Priority_Category == 'CRITICAL' %}clinical-change{% else %}clinical-stable{% endif %}">
                                        {{ variant.Priority_Category }}
                                    </span>
                                    {% if variant.get('Priority_Score') %}
                                        <br><small>({{ variant.Priority_Score }})</small>
                                    {% endif %}
                                {% else %}
                                    {{ variant.get('Rank', 'N/A') }}
                                {% endif %}
                            </td>
                            <td>{{ variant.get('Chromosome', 'N/A') }}:{{ variant.get('Position_hg19', 'N/A') }}</td>
                            <!-- Gene hg19/hg38 concatenated  -->
                            <td>
                                {% set gene_hg19_raw = variant.get('Gene_hg19', 'N/A') %}
                                {% set gene_hg38_raw = variant.get('Gene_hg38', 'N/A') %}
                                
                                {# Normalize NaN and empty values #}
                                {% set gene_hg19 = 'N/A' if gene_hg19_raw | string == 'nan' or gene_hg19_raw == '' else gene_hg19_raw %}
                                {% set gene_hg38 = 'N/A' if gene_hg38_raw | string == 'nan' or gene_hg38_raw == '' else gene_hg38_raw %}
                                
                                {% if gene_hg19 == gene_hg38 %}
                                    <span class="clinical-stable">{{ gene_hg19 }}</span>
                                {% else %}
                                    <span class="clinical-change">{{ gene_hg19 }} → {{ gene_hg38 }}</span>
                                {% endif %}
                            </td>
                            <!-- HGVS Matches with HGVSp priority + grouping + clear labeling -->
                            <td style="font-size: 10px; max-width: 250px; vertical-align: top; line-height: 1.2;">
                                {% set hgvsc_concordant_raw = variant.get('HGVSc_MATCHED_concordant', '') %}
                                {% set hgvsc_discordant_raw = variant.get('HGVSc_MATCHED_discordant', '') %}
                                {% set hgvsp_concordant_raw = variant.get('HGVSp_MATCHED_concordant', '') %}
                                {% set hgvsp_discordant_raw = variant.get('HGVSp_MATCHED_discordant', '') %}
                                
                                {# Convert to string and handle NaN/float values #}
                                {% set hgvsc_concordant = hgvsc_concordant_raw | string if hgvsc_concordant_raw is not none else '' %}
                                {% set hgvsc_discordant = hgvsc_discordant_raw | string if hgvsc_discordant_raw is not none else '' %}
                                {% set hgvsp_concordant = hgvsp_concordant_raw | string if hgvsp_concordant_raw is not none else '' %}
                                {% set hgvsp_discordant = hgvsp_discordant_raw | string if hgvsp_discordant_raw is not none else '' %}
                                
                                {# Clean up values #}
                                {% set hgvsc_concordant = '' if hgvsc_concordant in ['', '-', 'nan'] else hgvsc_concordant %}
                                {% set hgvsc_discordant = '' if hgvsc_discordant in ['', '-', 'nan'] else hgvsc_discordant %}
                                {% set hgvsp_concordant = '' if hgvsp_concordant in ['', '-', 'nan'] else hgvsp_concordant %}
                                {% set hgvsp_discordant = '' if hgvsp_discordant in ['', '-', 'nan'] else hgvsp_discordant %}
                                
                                {# Parse HGVSc lists for summary #}
                                {% set hgvsc_concordant_items = [] %}
                                {% set hgvsc_discordant_items = [] %}
                                
                                {% if hgvsc_concordant %}
                                    {% for item in hgvsc_concordant.replace(';', ',').split(',') %}
                                        {% set clean_item = item.strip() %}
                                        {% if clean_item and clean_item != '' %}
                                            {% set _ = hgvsc_concordant_items.append(clean_item) %}
                                        {% endif %}
                                    {% endfor %}
                                {% endif %}
                                
                                {% if hgvsc_discordant %}
                                    {% for item in hgvsc_discordant.replace(';', ',').split(',') %}
                                        {% set clean_item = item.strip() %}
                                        {% if clean_item and clean_item != '' %}
                                            {% set _ = hgvsc_discordant_items.append(clean_item) %}
                                        {% endif %}
                                    {% endfor %}
                                {% endif %}
                                
                                {# Parse HGVSp lists for details #}
                                {% set hgvsp_concordant_items = [] %}
                                {% set hgvsp_discordant_items = [] %}
                                
                                {% if hgvsp_concordant %}
                                    {% for item in hgvsp_concordant.replace(';', ',').split(',') %}
                                        {% set clean_item = item.strip() %}
                                        {% if clean_item and clean_item != '' %}
                                            {% set _ = hgvsp_concordant_items.append(clean_item) %}
                                        {% endif %}
                                    {% endfor %}
                                {% endif %}
                                
                                {% if hgvsp_discordant %}
                                    {% for item in hgvsp_discordant.replace(';', ',').split(',') %}
                                        {% set clean_item = item.strip() %}
                                        {% if clean_item and clean_item != '' %}
                                            {% set _ = hgvsp_discordant_items.append(clean_item) %}
                                        {% endif %}
                                    {% endfor %}
                                {% endif %}
                                
                                {# Count canonical transcripts for HGVSc summary #}
                                {% set hgvsc_concordant_canonical_count = 0 %}
                                {% set hgvsc_discordant_canonical_count = 0 %}
                                
                                {% for item in hgvsc_concordant_items %}
                                    {% if 'NM_' in item %}
                                        {% set hgvsc_concordant_canonical_count = hgvsc_concordant_canonical_count + 1 %}
                                    {% endif %}
                                {% endfor %}
                                
                                {% for item in hgvsc_discordant_items %}
                                    {% if 'NM_' in item %}
                                        {% set hgvsc_discordant_canonical_count = hgvsc_discordant_canonical_count + 1 %}
                                    {% endif %}
                                {% endfor %}
                                
                                {# Separate canonical vs other for details display #}
                                {% set hgvsp_discordant_canonical = [] %}
                                {% set hgvsp_discordant_other = [] %}
                                {% set hgvsc_discordant_canonical = [] %}
                                {% set hgvsc_discordant_other = [] %}
                                
                                {% for item in hgvsp_discordant_items %}
                                    {% if 'NM_' in item %}
                                        {% set _ = hgvsp_discordant_canonical.append(item) %}
                                    {% else %}
                                        {% set _ = hgvsp_discordant_other.append(item) %}
                                    {% endif %}
                                {% endfor %}
                                
                                {% for item in hgvsc_discordant_items %}
                                    {% if 'NM_' in item %}
                                        {% set _ = hgvsc_discordant_canonical.append(item) %}
                                    {% else %}
                                        {% set _ = hgvsc_discordant_other.append(item) %}
                                    {% endif %}
                                {% endfor %}
                                
                                {# Group HGVSp canonical by identical HGVS changes #}
                                {% set grouped_hgvsp_canonical = {} %}
                                {% for item in hgvsp_discordant_canonical %}
                                    {% if '(' in item %}
                                        {% set tx_id = item.split('(')[0].strip() %}
                                        {% set hgvs_part = '(' + item.split('(')[1] %}
                                        {% if hgvs_part in grouped_hgvsp_canonical %}
                                            {% set current_txs = grouped_hgvsp_canonical[hgvs_part] %}
                                            {% set _ = grouped_hgvsp_canonical.update({hgvs_part: current_txs + ', ' + tx_id}) %}
                                        {% else %}
                                            {% set _ = grouped_hgvsp_canonical.update({hgvs_part: tx_id}) %}
                                        {% endif %}
                                    {% else %}
                                        {% set _ = grouped_hgvsp_canonical.update({item: item}) %}
                                    {% endif %}
                                {% endfor %}
                                
                                {# Group HGVSc canonical by identical HGVS changes #}
                                {% set grouped_hgvsc_canonical = {} %}
                                {% for item in hgvsc_discordant_canonical %}
                                    {% if '(' in item %}
                                        {% set tx_id = item.split('(')[0].strip() %}
                                        {% set hgvs_part = '(' + item.split('(')[1] %}
                                        {% if hgvs_part in grouped_hgvsc_canonical %}
                                            {% set current_txs = grouped_hgvsc_canonical[hgvs_part] %}
                                            {% set _ = grouped_hgvsc_canonical.update({hgvs_part: current_txs + ', ' + tx_id}) %}
                                        {% else %}
                                            {% set _ = grouped_hgvsc_canonical.update({hgvs_part: tx_id}) %}
                                        {% endif %}
                                    {% else %}
                                        {% set _ = grouped_hgvsc_canonical.update({item: item}) %}
                                    {% endif %}
                                {% endfor %}
                                
                                {% if hgvsc_concordant_items or hgvsc_discordant_items %}
                                    <div class="{% if hgvsc_discordant_items %}clinical-change{% else %}clinical-stable{% endif %}">
                                        {# HGVSc Summary (as requested) #}
                                        <strong>Summary:</strong> {{ hgvsc_concordant_items | length }} concordant ({{ hgvsc_concordant_canonical_count }} canonical), {{ hgvsc_discordant_items | length }} discordant ({{ hgvsc_discordant_canonical_count }} canonical)<br><br>
                                        
                                        {# Details: HGVSp first (protein-level priority) with grouping #}
                                        {% if grouped_hgvsp_canonical %}
                                            <strong style="color: #c62828;">Discordant HGVSp:</strong><br>
                                            {% for hgvs_part, tx_list in grouped_hgvsp_canonical.items() %}
                                                <div style="margin-left: 4px; margin-top: 2px;">
                                                    <span style="white-space: nowrap; font-weight: bold;">{{ tx_list }}</span><br>
                                                    <span style="margin-left: 8px; word-break: break-all;">{{ hgvs_part | replace('→', '<br>&nbsp;&nbsp;&nbsp;&nbsp;→') | safe }}</span>
                                                </div>
                                            {% endfor %}
                                            <br>
                                        {% elif grouped_hgvsc_canonical %}
                                            {# Fallback to HGVSc when no HGVSp available - with grouping #}
                                            <strong style="color: #c62828;">Discordant HGVSc:</strong><br>
                                            {% for hgvs_part, tx_list in grouped_hgvsc_canonical.items() %}
                                                <div style="margin-left: 4px; margin-top: 2px;">
                                                    <span style="white-space: nowrap; font-weight: bold;">{{ tx_list }}</span><br>
                                                    <span style="margin-left: 8px; word-break: break-all;">{{ hgvs_part | replace('→', '<br>&nbsp;&nbsp;&nbsp;&nbsp;→') | safe }}</span>
                                                </div>
                                            {% endfor %}
                                            <br>
                                        {% endif %}
                                        
                                        {# Show other discordant if no canonical shown #}
                                        {% if not grouped_hgvsp_canonical and not grouped_hgvsc_canonical %}
                                            {% if hgvsp_discordant_other %}
                                                <strong style="color: #c62828;">Discordant HGVSp (Other, top 3):</strong><br>
                                                {% for item in hgvsp_discordant_other[:3] %}
                                                    <div style="margin-left: 4px; margin-top: 2px;">
                                                        {% if item | length > 60 %}
                                                            <span style="word-break: break-all;">{{ item[:60] }}...</span>
                                                        {% else %}
                                                            <span style="word-break: break-all;">{{ item }}</span>
                                                        {% endif %}
                                                    </div>
                                                {% endfor %}
                                            {% elif hgvsc_discordant_other %}
                                                <strong style="color: #c62828;">Discordant HGVSc (Other, top 3):</strong><br>
                                                {% for item in hgvsc_discordant_other[:3] %}
                                                    <div style="margin-left: 4px; margin-top: 2px;">
                                                        {% if item | length > 60 %}
                                                            <span style="word-break: break-all;">{{ item[:60] }}...</span>
                                                        {% else %}
                                                            <span style="word-break: break-all;">{{ item }}</span>
                                                        {% endif %}
                                                    </div>
                                                {% endfor %}
                                            {% endif %}
                                        {% endif %}
                                    </div>
                                {% else %}
                                    N/A
                                {% endif %}
                            </td>
                            <!-- Clinical significance column -->
                            <td>
                                {% set hg19_clin = variant.get('Clinical_Significance_hg19', 'N/A') %}
                                {% set hg38_clin = variant.get('Clinical_Significance_hg38', 'N/A') %}
                                {% if hg19_clin != hg38_clin %}
                                    <span class="clinical-change">{{ hg19_clin }} → {{ hg38_clin }}</span>
                                {% else %}
                                    <span class="clinical-stable">{{ hg19_clin }}</span>
                                {% endif %}
                            </td>
                            <!-- Impact level column -->
                            <td>
                                {% set hg19_impact = variant.get('Impact_hg19', 'N/A') %}
                                {% set hg38_impact = variant.get('Impact_hg38', 'N/A') %}
                                {% if hg19_impact != hg38_impact %}
                                    <span class="clinical-change">{{ hg19_impact }} → {{ hg38_impact }}</span>
                                {% else %}
                                    <span class="clinical-stable">{{ hg19_impact }}</span>
                                {% endif %}
                            </td>
                            <td>
                                {{ format_consequence_relationship(variant) }}
                            </td>
                        </tr>
                        {% endfor %}
                    </tbody>
                </table>
            </div>
        </div>
        {% endif %}

        <!-- DATA COVERAGE ANALYSIS -->
        <div class="section">
            <h2>Data coverage analysis</h2>
            
            <h3>Evidence distribution by build</h3>
            <div class="metrics-grid">
                <div class="metric-card">
                    <div class="metric-value">{{ get_summary_value(['prioritization', 'clinical_coverage', 'clinical_annotations', 'percentage_with_annotations']) }}%</div>
                    <div class="metric-label">Clinical data coverage</div>
                </div>
                <div class="metric-card">
                    <div class="metric-value">{{ get_summary_value(['prioritization', 'clinical_coverage', 'pathogenicity_predictions', 'percentage_with_predictions']) }}%</div>
                    <div class="metric-label">Prediction coverage</div>
                </div>
            </div>
            
            <!-- Show basic clinical categories if detailed data not available -->
            <div class="table-container">
                <table>
                    <thead>
                        <tr>
                            <th>Category</th>
                            <th>hg19 count</th>
                            <th>hg38 count</th>
                        </tr>
                    </thead>
                    <tbody>
                        <tr>
                            <td>PATHOGENIC</td>
                            <td>{{ get_summary_value(['prioritization', 'clinical_coverage', 'clinical_annotations', 'hg19_distribution', 'PATHOGENIC'], 0) }}</td>
                            <td>{{ get_summary_value(['prioritization', 'clinical_coverage', 'clinical_annotations', 'hg38_distribution', 'PATHOGENIC'], 0) }}</td>
                        </tr>
                        <tr>
                            <td>BENIGN</td>
                            <td>{{ get_summary_value(['prioritization', 'clinical_coverage', 'clinical_annotations', 'hg19_distribution', 'BENIGN'], 0) }}</td>
                            <td>{{ get_summary_value(['prioritization', 'clinical_coverage', 'clinical_annotations', 'hg38_distribution', 'BENIGN'], 0) }}</td>
                        </tr>
                        <tr>
                            <td>VUS</td>
                            <td>{{ get_summary_value(['prioritization', 'clinical_coverage', 'clinical_annotations', 'hg19_distribution', 'VUS'], 0) }}</td>
                            <td>{{ get_summary_value(['prioritization', 'clinical_coverage', 'clinical_annotations', 'hg38_distribution', 'VUS'], 0) }}</td>
                        </tr>
                        <tr>
                            <td>NONE</td>
                            <td>{{ get_summary_value(['prioritization', 'clinical_coverage', 'clinical_annotations', 'hg19_distribution', 'NONE'], 0) }}</td>
                            <td>{{ get_summary_value(['prioritization', 'clinical_coverage', 'clinical_annotations', 'hg38_distribution', 'NONE'], 0) }}</td>
                        </tr>
                    </tbody>
                </table>
            </div>
        </div>

        <!-- TECHNICAL NOTES -->
        <div class="section">
            <h2>Technical Notes</h2>
            <div class="summary-text">
                <p><strong>Coordinate system:</strong> All coordinates use VEP normalization (SNVs: original position, Indels: original position + 1)</p>
                <p><strong>Evidence-first scoring:</strong> Prioritizes clinical significance changes over annotation differences between builds</p>
                <p><strong>Priority categories:</strong> CRITICAL (immediate review) → HIGH (priority review) → MODERATE (standard review) → LOW (secondary review)</p>
                <p><strong>Quality control:</strong> Liftover concordance analysis compares CrossMap and bcftools coordinate conversion results</p>
                <p><strong>For complete details:</strong> Refer to the accompanying summary text files and CSV output for comprehensive analysis results.</p>
            </div>
        </div>
    </div>
</body>
</html>'''
        
        template = Template(template_str)
        template.globals['get_summary_value'] = self._get_summary_value
        template.globals['format_consequence_relationship'] = self._format_consequence_relationship  
        return template.render(**self.report_data)


def main():
    parser = argparse.ArgumentParser(
        description='Generate HTML report from CrossBuild Assessor outputs'
    )
    parser.add_argument('--input-dir', '-i', required=True,
                       help='Directory containing analysis outputs')
    parser.add_argument('--output', '-o', default='crossbuild_report.html',
                       help='Output HTML file (default: crossbuild_report.html)')
    
    args = parser.parse_args()
    
    try:
        generator = ReportGenerator(args.input_dir)
        output_file = generator.generate_report(args.output)
        
        print(f"\n✓ Report generated successfully!")
        print(f"Open {output_file} in a web browser to view the report.")
        
    except Exception as e:
        print(f"Error generating report: {e}")
        return 1

if __name__ == "__main__":
    main()