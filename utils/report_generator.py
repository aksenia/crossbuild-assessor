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
    #    self._debug_summary_data() # DEBUG
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
        """Load variant data for clinical evidence table : Show only variants with actual changes"""
        csv_file = self.input_dir / 'prioritized_variants.csv'
        if csv_file.exists():
            df = pd.read_csv(csv_file)
            
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
                    'Rank', 'Chromosome_hg19', 'Position_hg19', 'Gene_hg19', 'Gene_hg38',
                    'Priority_Score', 'Priority_Category', 'Discordance_Summary',
                    'Priority_Transcript_CrossBuild', 'MANE_Flag_hg38', 'HGVS_c_hg19', 'HGVS_c_hg38', 
                    'HGVS_p_hg19', 'HGVS_p_hg38', 'HGVS_c_Concordance', 'HGVS_p_Concordance',
                    'Consequence_Relationship', 'Consequence_Change'
                ]
                
                # Only include columns that exist
                available_columns = [col for col in clinical_columns if col in df_sorted.columns]
                
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
                    <div class="metric-value">{{ get_summary_value(['prioritization', 'clinical_transitions', 'total_changing']) }}</div>
                    <div class="metric-label">Clinical significance changes</div>
                </div>
            </div>
            
            <div class="summary-text">
                <strong>Key findings:</strong> 
                {{ get_summary_value(['prioritization', 'dataset_overview', 'total_discordant_variants']) }} discordant variants analyzed with {{ get_summary_value(['prioritization', 'dataset_overview', 'variants_in_excel_output']) }} variants included in Excel output.
                {{ get_summary_value(['prioritization', 'top_variants_summary', 'critical_variants_distinct']) }} variants require immediate review.
                {{ get_summary_value(['prioritization', 'clinical_transitions', 'total_changing']) }} variants show clinical significance changes that may affect interpretation.
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
            <strong>Priority definitions (4-Category System):</strong>
            <strong>CRITICAL:</strong> HGVS mismatches on priority transcript OR major clinical significance changes (Pathogenic↔Benign) |
            <strong>MODERATE:</strong> Priority transcript unavailable OR moderate clinical changes OR serious functional differences |
            <strong>LOW:</strong> Minor changes (VUS-Benign, prediction changes, gene symbol differences) |
            <strong>CONCORDANT:</strong> Perfect priority transcript match with identical HGVS nomenclature
        </div>
        <div class="metrics-grid">
            <div class="metric-card">
                <div class="metric-value">{{ get_summary_value(['prioritization', 'priority_distribution', 'CRITICAL'], 0) }}</div>
                <div class="metric-label">CRITICAL priority</div>
            </div>
            <div class="metric-card">
                <div class="metric-value">{{ get_summary_value(['prioritization', 'priority_distribution', 'MODERATE'], 0) }}</div>
                <div class="metric-label">MODERATE priority</div>
            </div>
            <div class="metric-card">
                <div class="metric-value">{{ get_summary_value(['prioritization', 'priority_distribution', 'LOW'], 0) }}</div>
                <div class="metric-label">LOW priority</div>
            </div>
            <div class="metric-card">
                <div class="metric-value">{{ get_summary_value(['prioritization', 'priority_distribution', 'CONCORDANT'], 0) }}</div>
                <div class="metric-label">CONCORDANT</div>
            </div>
        </div>

            <!-- PRIORITY TRANSCRIPT HGVS ANALYSIS -->
            <h3>Priority transcript HGVS concordance</h3>
            <div class="metrics-grid">
                <div class="metric-card">
                    <div class="metric-value">{{ get_summary_value(['prioritization', 'priority_transcript_analysis', 'hgvsc_concordance_rate'], 'N/A') }}%</div>
                    <div class="metric-label">Priority HGVSc Concordance</div>
                </div>
                <div class="metric-card">
                    <div class="metric-value">{{ get_summary_value(['prioritization', 'priority_transcript_analysis', 'hgvsp_concordance_rate'], 'N/A') }}%</div>
                    <div class="metric-label">Priority HGVSp Concordance</div>
                </div>
                <div class="metric-card">
                    <div class="metric-value">{{ get_summary_value(['prioritization', 'priority_transcript_analysis', 'mane_select_both_builds'], 'N/A') }}</div>
                    <div class="metric-label">MANE Select Both Builds</div>
                </div>
                <div class="metric-card">
                    <div class="metric-value">{{ get_summary_value(['prioritization', 'priority_transcript_analysis', 'mane_hg38_only'], 'N/A') }}</div>
                    <div class="metric-label">MANE hg38 Only</div>
                </div>
            </div>

            <!-- Show basic clinical categories if detailed data not available -->
            <h3>Clinical significance annotations from ClinVar </h3>
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
            <h3>Discrepancies visualizations</h3>
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
                                <th>Priority Transcript HGVS</th>
                                <th>Score Breakdown</th>
                                <th>Consequence Relationship</th>
                            </tr>
                        </thead>
                        <tbody>
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
                                <td>{{ variant.get('Chromosome_hg19', 'N/A') }}:{{ variant.get('Position_hg19', 'N/A') }}</td>
                                <td>
                                    {% set gene_hg19_raw = variant.get('Gene_hg19', 'N/A') %}
                                    {% set gene_hg38_raw = variant.get('Gene_hg38', 'N/A') %}
                                    {% set gene_hg19 = 'N/A' if gene_hg19_raw | string == 'nan' or gene_hg19_raw == '' else gene_hg19_raw %}
                                    {% set gene_hg38 = 'N/A' if gene_hg38_raw | string == 'nan' or gene_hg38_raw == '' else gene_hg38_raw %}
                                    
                                    {% if gene_hg19 == gene_hg38 %}
                                        <span class="clinical-stable">{{ gene_hg19 }}</span>
                                    {% else %}
                                        <span class="clinical-change">{{ gene_hg19 }} → {{ gene_hg38 }}</span>
                                    {% endif %}
                                </td>
                                <td style="font-size: 11px;">
                                    {% set transcript_id = variant.get('Priority_Transcript_CrossBuild', 'N/A') %}
                                    {% set mane_flag = variant.get('MANE_Flag_hg38', '') %}
                                    
                                    {# Normalize values to handle nan, None, empty strings #}
                                    {% set transcript_clean = 'N/A' if transcript_id in ['NONE', 'nan', None, ''] or transcript_id|string == 'nan' else transcript_id %}
                                    {% set mane_clean = '' if mane_flag in ['nan', None, ''] or mane_flag|string == 'nan' else mane_flag %}
                                    
                                    <strong>Transcript:</strong> 
                                    {{ transcript_clean }}
                                    {% if mane_clean and mane_clean != 'N/A' and mane_clean != '' %}
                                        <span style="color: #2c5f8a; font-weight: bold;">({{ mane_clean }})</span>
                                    {% endif %}
                                    <br>
                                    
                                    {# Normalize HGVS values #}
                                    {% set hgvsc_hg19_raw = variant.get('HGVS_c_hg19', 'N/A') %}
                                    {% set hgvsc_hg38_raw = variant.get('HGVS_c_hg38', 'N/A') %}
                                    {% set hgvsc_hg19 = 'N/A' if hgvsc_hg19_raw|string == 'nan' else hgvsc_hg19_raw %}
                                    {% set hgvsc_hg38 = 'N/A' if hgvsc_hg38_raw|string == 'nan' else hgvsc_hg38_raw %}
                                    
                                    <strong>HGVSc:</strong> 
                                    {% if hgvsc_hg19 == hgvsc_hg38 %}
                                        {{ hgvsc_hg19 }}
                                    {% else %}
                                        <span class="clinical-change">{{ hgvsc_hg19 }} → {{ hgvsc_hg38 }}</span>
                                    {% endif %}
                                    <br>
                                    
                                    {# Normalize HGVSp values #}
                                    {% set hgvsp_hg19_raw = variant.get('HGVS_p_hg19', 'N/A') %}
                                    {% set hgvsp_hg38_raw = variant.get('HGVS_p_hg38', 'N/A') %}
                                    {% set hgvsp_hg19 = 'N/A' if hgvsp_hg19_raw|string == 'nan' else hgvsp_hg19_raw %}
                                    {% set hgvsp_hg38 = 'N/A' if hgvsp_hg38_raw|string == 'nan' else hgvsp_hg38_raw %}
                                    
                                    {% if hgvsp_hg19 != 'N/A' or hgvsp_hg38 != 'N/A' %}
                                        <strong>HGVSp:</strong> 
                                        {% if hgvsp_hg19 == hgvsp_hg38 %}
                                            {{ hgvsp_hg19 }}
                                        {% else %}
                                            <span class="clinical-change">{{ hgvsp_hg19 }} → {{ hgvsp_hg38 }}</span>
                                        {% endif %}
                                    {% endif %}
                                </td>
                                <td style="font-size: 10px; max-width: 250px;">
                                    {{ variant.get('score_breakdown', variant.get('Discordance_Summary', 'N/A')) }}
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
            

            <!-- GENE-LEVEL TECHNICAL ANALYSIS -->
            <div class="section">
                <h2>Gene-level technical issues</h2>              
                <div class="metrics-grid">
                    <div class="metric-card">
                        <div class="metric-value">{{ get_summary_value(['prioritization', 'gene_technical_analysis', 'dataset_mean_technical_rate']) }}</div>
                        <div class="metric-label">Dataset mean technical rate</div>
                    </div>
                    <div class="metric-card">
                        <div class="metric-value">{{ get_summary_value(['prioritization', 'gene_technical_analysis', 'total_genes_analyzed']) }}</div>
                        <div class="metric-label">Total genes analyzed</div>
                    </div>
                    <div class="metric-card">
                        <div class="metric-value">{{ get_summary_value(['prioritization', 'gene_technical_analysis', 'flagged_genes_count']) }}</div>
                        <div class="metric-label">Genes above threshold</div>
                    </div>
                </div>
                
                <h3>Statistical outliers (≥ mean + 2SD)</h3>
                <div class="table-container">
                    <table>
                        <thead>
                            <tr>
                                <th>Gene</th>
                                <th>Technical issue rate</th>
                                <th>Total issues</th>
                                <th>Total variants</th>
                            </tr>
                        </thead>
                        <tbody>
                            {% set flagged_genes = get_summary_value(['prioritization', 'gene_technical_analysis', 'flagged_genes_above_average']) %}
                            {% if flagged_genes and flagged_genes != 'N/A' and flagged_genes|length > 0 %}
                                {% for gene in flagged_genes %}
                                <tr>
                                    <td><strong>{{ gene.gene }}</strong></td>
                                    <td class="clinical-change">{{ gene.technical_issue_rate }}</td>
                                    <td>{{ gene.total_technical_issues }}</td>
                                    <td>{{ gene.total_variants }}</td>
                                </tr>
                                {% endfor %}
                            {% else %}
                                <tr>
                                    <td colspan="4" class="no-data">No genes flagged as statistical outliers</td>
                                </tr>
                            {% endif %}
                        </tbody>
                    </table>
                </div>
            </div>

        <!-- TECHNICAL NOTES -->
        <div class="section">
            <h2>Technical Notes</h2>
            <div class="summary-text">
                <p><strong>Priority transcript-based scoring:</strong> HGVS concordance on priority transcript drives prioritization, using MANE-first transcript selection</p>
                <p><strong>Priority categories:</strong> CRITICAL (immediate review) → MODERATE (standard review) → LOW (secondary review) → CONCORDANT (no review needed)</p>
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