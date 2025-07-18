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
        
    def generate_report(self, output_file):
        """Generate complete HTML report"""
        print("Generating HTML report...")
        
        # Collect all data
        self._collect_metadata()
        self._collect_summary_data()
        self._collect_plot_images()
        self._collect_variant_data()
        
        # Generate HTML
        html_content = self._render_html()
        
        # Save report
        with open(output_file, 'w', encoding='utf-8') as f:
            f.write(html_content)
        
        print(f"✓ HTML report saved to: {output_file}")
        return output_file
    
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
        """Load variant data for clinical evidence table"""
        csv_file = self.input_dir / 'prioritized_variants.csv'
        if csv_file.exists():
            df = pd.read_csv(csv_file)
            
            # Debug: Print available columns to help identify the correct column name
            print(f"Available columns in CSV: {list(df.columns)}")
            
            # Top 10 variants with clinical evidence focus
            if len(df) > 0:
                # Select columns for clinical review - let's check what's actually available
                clinical_columns = [
                    'Rank', 'Chromosome', 'Position_hg19', 'Gene_hg19', 'Gene_hg38',
                    'Clinical_Significance_hg19', 'Clinical_Significance_hg38',
                    'Impact_hg19', 'Impact_hg38',
                    'SIFT_hg19', 'SIFT_hg38', 'PolyPhen_hg19', 'PolyPhen_hg38',
                    'Consequence_hg19', 'Consequence_hg38'
                ]
                
                # Only include columns that exist
                available_columns = [col for col in clinical_columns if col in df.columns]
                print(f"Using columns: {available_columns}")
                
                top_variants = df[available_columns].head(10)
                self.report_data['top_variants'] = top_variants.to_dict('records')
            else:
                self.report_data['top_variants'] = []
    
    def _get_summary_value(self, data_path, fallback='N/A'):
        """Safely extract values from nested summary data"""
        try:
            current = self.report_data['summaries']
            for key in data_path:
                current = current[key]
            return current
        except (KeyError, TypeError):
            return fallback
    
    def _render_html(self):
        """Render HTML using Jinja2 template with clinical focus"""
        template_str = """
<!DOCTYPE html>
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
        .impact-high { background: #ffebee; color: #c62828; font-weight: bold; }
        .impact-moderate { background: #fff3e0; color: #ef6c00; }
        .impact-low { background: #f3e5f5; color: #7b1fa2; }
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
                {% set total_variants = get_summary_value(['liftover', 'dataset_overview', 'total_variants']) %}
                {% set position_match = get_summary_value(['liftover', 'dataset_overview', 'position_match_percentage']) %}
                {% set critical_variants = get_summary_value(['prioritization', 'top_variants_summary', 'critical_variants_distinct']) %}
                {% set clinical_changes = get_summary_value(['prioritization', 'top_variants_summary', 'clinical_changes_count']) %}
                
                {{ position_match }}% of {{ total_variants }} variants show coordinate concordance between liftover tools. 
                {{ critical_variants }} variants require immediate review due to functional or clinical significance changes.
                {{ clinical_changes }} variants show clinical significance transitions that may affect interpretation.
            </div>
        </div>

        <!-- LIFTOVER QUALITY CONTROL -->
        {% if get_summary_value(['liftover']) %}
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

            <!-- Tool Performance Table -->
            {% if get_summary_value(['liftover', 'mapping_status_breakdown']) %}
            <h3>Mapping status breakdown</h3>
            <div class="table-container">
                <table>
                    <thead>
                        <tr>
                            <th>Mapping status</th>
                            <th>Variant count</th>
                            <th>Position match rate</th>
                            <th>Genotype match rate</th>
                        </tr>
                    </thead>
                    <tbody>
                        {% for status, data in get_summary_value(['liftover', 'mapping_status_breakdown']).items() %}
                        <tr>
                            <td>{{ status }}</td>
                            <td>{{ data.count }}</td>
                            <td>{{ data.pos_match_percentage }}%</td>
                            <td>{{ data.gt_match_percentage }}%</td>
                        </tr>
                        {% endfor %}
                    </tbody>
                </table>
            </div>
            {% endif %}

            <!-- Match Categories -->
            {% if get_summary_value(['liftover', 'match_categories']) %}
            <h3>Variant match categories</h3>
            <div class="table-container">
                <table>
                    <thead>
                        <tr>
                            <th>Match category</th>
                            <th>Count</th>
                            <th>Percentage</th>
                        </tr>
                    </thead>
                    <tbody>
                        {% set match_categories = get_summary_value(['liftover', 'match_categories']) %}
                        <tr>
                            <td>Both match</td>
                            <td>{{ match_categories.both_match.count }}</td>
                            <td>{{ match_categories.both_match.percentage }}%</td>
                        </tr>
                        <tr>
                            <td>Position only</td>
                            <td>{{ match_categories.position_only.count }}</td>
                            <td>{{ match_categories.position_only.percentage }}%</td>
                        </tr>
                        <tr>
                            <td>Genotype only</td>
                            <td>{{ match_categories.genotype_only.count }}</td>
                            <td>{{ match_categories.genotype_only.percentage }}%</td>
                        </tr>
                        <tr>
                            <td>Both mismatch</td>
                            <td>{{ match_categories.both_mismatch.count }}</td>
                            <td>{{ match_categories.both_mismatch.percentage }}%</td>
                        </tr>
                    </tbody>
                </table>
            </div>
            {% endif %}

            <!-- Flip/Swap Analysis -->
            {% if get_summary_value(['liftover', 'flip_swap_analysis']) %}
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
                                {% elif category == "Allele Swap Only" %}
                                    REF/ALT alleles swapped to match reference
                                {% elif category == "Strand Flip + Allele Swap" %}
                                    Both strand and allele corrections applied
                                {% elif category == "Swap Failed (Ambiguous)" %}
                                    Allele swap attempted but failed due to ambiguous alleles
                                {% else %}
                                    {{ category }}
                                {% endif %}
                            </td>
                        </tr>
                        {% endfor %}
                    </tbody>
                </table>
            </div>
            {% endif %}

            <!-- Liftover QC Plots -->
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
        {% endif %}

        <!-- VARIANT PRIORITIZATION -->
        {% if get_summary_value(['prioritization']) %}
        <div class="section">
            <h2>Variant Prioritization Analysis</h2>
            
            <h3>Priority distribution</h3>
            <div class="metrics-grid">
                {% if get_summary_value(['prioritization', 'priority_distribution']) %}
                {% for category, count in get_summary_value(['prioritization', 'priority_distribution']).items() %}
                <div class="metric-card">
                    <div class="metric-value">{{ count }}</div>
                    <div class="metric-label">{{ category }} priority</div>
                </div>
                {% endfor %}
                {% endif %}
            </div>

            <!-- Clinical Significance Transitions -->
            {% if get_summary_value(['prioritization', 'clinical_transitions']) %}
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
            
            {% if get_summary_value(['prioritization', 'clinical_transitions', 'directional_changes']) %}
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
                        {% for transition, data in get_summary_value(['prioritization', 'clinical_transitions', 'directional_changes']).items() %}
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
                    </tbody>
                </table>
            </div>
            {% endif %}
            {% endif %}

            <!-- Prioritization Plots -->
            {% if images.prioritization_plots %}
            <h3>Prioritization visualizations</h3>
            <div class="plot-container">
                <img src="{{ images.prioritization_plots }}" alt="Variant Prioritization">
            </div>
            {% endif %}
        </div>
        {% endif %}

        <!-- TOP PRIORITY VARIANTS -->
        {% if top_variants %}
        <div class="section">
            <h2>Top priority variants</h2>
            <div class="table-container">
                <table>
                    <thead>
                        <tr>
                            <th>Location</th>
                            <th>Gene hg19</th>
                            <th>Gene hg38</th>
                            <th>Clinical significance</th>
                            <th>Impact level</th>
                            <th>Pathogenicity predictions</th>
                            <th>Consequence</th>
                        </tr>
                    </thead>
                    <tbody>
                        {% for variant in top_variants %}
                        <tr>
                            <td>{{ variant.get('Chromosome', 'N/A') }}:{{ variant.get('Position_hg19', 'N/A') }}</td>
                            <td>{{ variant.get('Gene_hg19', 'N/A') }}</td>
                            <td>{{ variant.get('Gene_hg38', 'N/A') }}</td>
                            <td>
                                {% set hg19_clin = variant.get('Clinical_Significance_hg19', 'N/A') %}
                                {% set hg38_clin = variant.get('Clinical_Significance_hg38', 'N/A') %}
                                {% if hg19_clin != hg38_clin %}
                                    <span class="clinical-change">{{ hg19_clin }} → {{ hg38_clin }}</span>
                                {% else %}
                                    <span class="clinical-stable">{{ hg19_clin }}</span>
                                {% endif %}
                            </td>
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
                                {% set hg19_sift = variant.get('SIFT_hg19', '') %}
                                {% set hg38_sift = variant.get('SIFT_hg38', '') %}
                                {% set hg19_polyphen = variant.get('PolyPhen_hg19', '') %}
                                {% set hg38_polyphen = variant.get('PolyPhen_hg38', '') %}
                                
                                {% if hg19_sift and hg38_sift and hg19_sift != hg38_sift %}
                                    SIFT: {{ hg19_sift }} → {{ hg38_sift }}<br>
                                {% elif hg19_sift %}
                                    SIFT: {{ hg19_sift }}<br>
                                {% endif %}
                                
                                {% if hg19_polyphen and hg38_polyphen and hg19_polyphen != hg38_polyphen %}
                                    PolyPhen: {{ hg19_polyphen }} → {{ hg38_polyphen }}
                                {% elif hg19_polyphen %}
                                    PolyPhen: {{ hg19_polyphen }}
                                {% endif %}
                                
                                {% if not (hg19_sift or hg19_polyphen) %}
                                    N/A
                                {% endif %}
                            </td>
                            <td>
                                {% set hg19_consequence = variant.get('Consequence_hg19', 'N/A') %}
                                {% set hg38_consequence = variant.get('Consequence_hg38', 'N/A') %}
                                {% if hg19_consequence != hg38_consequence %}
                                    <span class="clinical-change">{{ hg19_consequence }} → {{ hg38_consequence }}</span>
                                {% else %}
                                    {{ hg19_consequence }}
                                {% endif %}
                            </td>
                        </tr>
                        {% endfor %}
                    </tbody>
                </table>
            </div>
        </div>
        {% endif %}

        <!-- DATA COVERAGE ANALYSIS -->
        {% if get_summary_value(['prioritization', 'clinical_coverage']) %}
        <div class="section">
            <h2>Data coverage analysis</h2>
            
            {% if get_summary_value(['prioritization', 'clinical_coverage', 'clinical_annotations']) %}
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
            
            {% set build_comparison = get_summary_value(['prioritization', 'clinical_coverage', 'clinical_annotations', 'build_comparison']) %}
            {% if build_comparison %}
            <div class="table-container">
                <table>
                    <thead>
                        <tr>
                            <th>Category</th>
                            <th>hg19 count</th>
                            <th>hg38 count</th>
                            <th>Net change</th>
                        </tr>
                    </thead>
                    <tbody>
                        {% for category, data in build_comparison.items() %}
                        {% if data.hg19_count > 0 or data.hg38_count > 0 %}
                        <tr>
                            <td>{{ category }}</td>
                            <td>{{ data.hg19_count }}</td>
                            <td>{{ data.hg38_count }}</td>
                            <td>
                                {% if data.difference > 0 %}
                                    +{{ data.difference }}
                                {% elif data.difference < 0 %}
                                    {{ data.difference }}
                                {% else %}
                                    0
                                {% endif %}
                            </td>
                        </tr>
                        {% endif %}
                        {% endfor %}
                    </tbody>
                </table>
            </div>
            {% endif %}
            {% endif %}
        </div>
        {% endif %}

        <!-- TECHNICAL NOTES -->
        <div class="section">
            <h2>Technical Notes</h2>
            <div class="summary-text">
                <p><strong>Coordinate system:</strong> All coordinates use VEP normalization (SNVs: original position, Indels: original position + 1)</p>
                <p><strong>Evidence-first scoring:</strong> Prioritizes clinical significance changes over annotation differences between builds</p>
                <p><strong>Priority categories:</strong> CRITICAL (immediate review) → HIGH (priority review) → MODERATE (standard review) → LOW (secondary review)</p>
                <p><strong>Quality control:</strong> Liftover concordance analysis compares CrossMap and bcftools coordinate conversion results</p>
                
                {% if get_summary_value(['prioritization', 'metadata', 'timestamp']) %}
                <p><strong>Analysis timestamp:</strong> {{ get_summary_value(['prioritization', 'metadata', 'timestamp']) }}</p>
                {% endif %}
                {% if get_summary_value(['liftover', 'metadata', 'timestamp']) %}
                <p><strong>Liftover analysis:</strong> {{ get_summary_value(['liftover', 'metadata', 'timestamp']) }}</p>
                {% endif %}
                
                <p><strong>For complete details:</strong> Refer to the accompanying summary text files and CSV output for comprehensive analysis results.</p>
            </div>
        </div>
    </div>
</body>
</html>
        """
        
        template = Template(template_str)
        # Add only the data access helper function
        template.globals['get_summary_value'] = self._get_summary_value
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