"""
Prioritization Plot Generator

Creates comprehensive 4-plot visualization for variant prioritization analysis
with clinical evidence categorization and discordance type breakdown.
"""

import matplotlib.pyplot as plt
import pandas as pd
import numpy as np
from pathlib import Path


class PrioritizationPlotter:
    """Generate prioritization analysis plots with clinical evidence focus"""
    
    def __init__(self, plot_colors, figure_config):
        """
        Initialize plotter with color schemes and figure configuration
        
        Args:
            plot_colors: Dictionary with color palettes for different plot types
            figure_config: Dictionary with figure layout and styling configuration
        """
        self.colors_priority = plot_colors['priority']
        self.colors_clinical = plot_colors['clinical']
        self.figure_config = figure_config
    
    def create_all_plots(self, df, conn, output_dir):
        """Create enhanced visualization plots for variant prioritization analysis"""
        print("Creating enhanced prioritization visualizations...")
        
        # Create enhanced visualization with 4 plots
        fig, axes = plt.subplots(
            self.figure_config['main_figure']['rows'], 
            self.figure_config['main_figure']['cols'], 
            figsize=self.figure_config['main_figure']['figsize']
        )
        fig.suptitle('Enhanced Variant Prioritization Analysis', fontsize=16, fontweight='bold')
        
        # PLOT 1: Primary Discordance Types (moved to first position)
        self._plot_discordance_types(df, axes[0,0])
        
        # PLOT 2: Impact Level Transitions (hg19 vs hg38)
        self._plot_impact_transitions(df, axes[0,1])
        
        # PLOT 3: Clinical Significance Transitions (hg19 vs hg38)
        self._plot_clinical_significance_transitions(df, axes[1,0])
        
        # PLOT 4: Directional Clinical Significance Changes
        self._plot_directional_clinical_changes(df, axes[1,1])
        
        plt.tight_layout()
        
        # Save plot
        output_path = Path(output_dir) / 'variant_prioritization_plots.png'
        plt.savefig(output_path, dpi=300, bbox_inches='tight')
        plt.show()
        
        # Print summary statistics
        self._print_summary_statistics(df, conn)
    
    def _plot_discordance_types(self, df, ax):
        """Plot primary discordance types (moved from position 3 to 1)"""
        print("1. Creating primary discordance types visualization...")
        
        if len(df) > 0:
            df['discordance_primary'] = df['discordance_summary'].apply(self._categorize_discordance_primary)
            discordance_counts = df['discordance_primary'].value_counts()
            
            # Order categories by clinical importance
            discordance_order = [
                'Same Transcript\nConsequence Change',
                'Clinical Significance\nChange',
                'Impact Level\nChange',
                'Gene Annotation\nChange',
                'Pathogenicity\nPrediction Change',
                'Unmatched\nConsequences',
                'Position\nMismatch',
                'Genotype\nMismatch',
                'Other/Multiple'
            ]
            
            discordance_counts = discordance_counts.reindex([cat for cat in discordance_order if cat in discordance_counts.index])
            
            bars = ax.bar(range(len(discordance_counts)), discordance_counts.values,
                            color=self.colors_priority[:len(discordance_counts)], alpha=0.8)
            
            ax.set_xlabel('Primary Discordance Type', fontweight='bold')
            ax.set_ylabel('Variant Count', fontweight='bold')
            ax.set_title('Variants by Primary Discordance Type\n(Aggregated Categories)',
                           fontsize=12, fontweight='bold')
            ax.set_xticks(range(len(discordance_counts)))
            ax.set_xticklabels(discordance_counts.index, rotation=45, ha='right', fontsize=9)
            ax.grid(True, alpha=0.3, axis='y')
            
            # Add count labels on bars
            for i, (bar, count) in enumerate(zip(bars, discordance_counts.values)):
                ax.text(bar.get_x() + bar.get_width()/2, bar.get_height() + count*0.01,
                          str(count), ha='center', va='bottom', fontweight='bold', fontsize=10)
        else:
            ax.text(0.5, 0.5, 'No discordant variants\nfound', ha='center', va='center',
                      fontsize=12, transform=ax.transAxes)
            ax.set_title('Variants by Primary Discordance Type\n(Aggregated Categories)',
                           fontsize=12, fontweight='bold')
    
    def _plot_impact_transitions(self, df, ax):
        """Plot impact level transitions (hg19 vs hg38)"""
        print("2. Creating impact level transitions...")
        
        if len(df) > 0:
            # Create transition matrix
            impact_levels = ['HIGH', 'MODERATE', 'LOW', 'MODIFIER']
            transition_matrix = pd.crosstab(df['hg19_impact'], df['hg38_impact'], 
                                          rownames=['hg19 Impact'], colnames=['hg38 Impact'])
            
            # Reindex to ensure all impact levels are present
            transition_matrix = transition_matrix.reindex(index=impact_levels, columns=impact_levels, fill_value=0)
            
            # Create stacked bar chart
            bottom = np.zeros(len(impact_levels))
            for i, hg38_impact in enumerate(impact_levels):
                values = [transition_matrix.loc[hg19_impact, hg38_impact] for hg19_impact in impact_levels]
                ax.bar(impact_levels, values, bottom=bottom, 
                      label=f'hg38: {hg38_impact}', color=self.colors_clinical[i], alpha=0.8)
                bottom += values
            
            ax.set_xlabel('hg19 Impact Level', fontweight='bold')
            ax.set_ylabel('Variant Count', fontweight='bold')
            ax.set_title('Impact Level Transitions\n(hg19 → hg38)', fontsize=12, fontweight='bold')
            ax.legend(title='hg38 Impact', bbox_to_anchor=(1.05, 1), loc='upper left', fontsize=8)
            ax.grid(True, alpha=0.3, axis='y')
        else:
            ax.text(0.5, 0.5, 'No variants\nfound', ha='center', va='center',
                      fontsize=12, transform=ax.transAxes)
            ax.set_title('Impact Level Transitions\n(hg19 → hg38)', fontsize=12, fontweight='bold')
    
    def _plot_clinical_significance_transitions(self, df, ax):
        """Plot clinical significance transitions (hg19 vs hg38)"""
        print("3. Creating clinical significance transitions...")
        
        if len(df) > 0 and 'hg19_clin_sig_normalized' in df.columns:
            # Create transition matrix for normalized clinical significance
            clin_categories = ['PATHOGENIC', 'BENIGN', 'VUS', 'RISK', 'DRUG_RESPONSE', 'PROTECTIVE', 'OTHER', 'NONE']
            
            transition_matrix = pd.crosstab(df['hg19_clin_sig_normalized'], df['hg38_clin_sig_normalized'], 
                                          rownames=['hg19 Clinical Sig'], colnames=['hg38 Clinical Sig'])
            
            # Reindex to ensure all categories are present (only include those with data)
            present_hg19 = [cat for cat in clin_categories if cat in transition_matrix.index]
            present_hg38 = [cat for cat in clin_categories if cat in transition_matrix.columns]
            
            if present_hg19 and present_hg38:
                transition_matrix = transition_matrix.reindex(index=present_hg19, columns=present_hg38, fill_value=0)
                
                # Create stacked bar chart
                bottom = np.zeros(len(present_hg19))
                for i, hg38_clin in enumerate(present_hg38):
                    values = [transition_matrix.loc[hg19_clin, hg38_clin] for hg19_clin in present_hg19]
                    ax.bar(present_hg19, values, bottom=bottom, 
                          label=f'hg38: {hg38_clin}', color=self.colors_clinical[i % len(self.colors_clinical)], alpha=0.8)
                    bottom += values
                
                ax.set_xlabel('hg19 Clinical Significance', fontweight='bold')
                ax.set_ylabel('Variant Count', fontweight='bold')
                ax.set_title('Clinical Significance Transitions\n(hg19 → hg38)', fontsize=12, fontweight='bold')
                ax.legend(title='hg38 Clinical Sig', bbox_to_anchor=(1.05, 1), loc='upper left', fontsize=7)
                ax.tick_params(axis='x', rotation=45)
                ax.grid(True, alpha=0.3, axis='y')
            else:
                ax.text(0.5, 0.5, 'No clinical significance\ndata found', ha='center', va='center',
                          fontsize=12, transform=ax.transAxes)
                ax.set_title('Clinical Significance Transitions\n(hg19 → hg38)', fontsize=12, fontweight='bold')
        else:
            ax.text(0.5, 0.5, 'No clinical significance\ndata available', ha='center', va='center',
                      fontsize=12, transform=ax.transAxes)
            ax.set_title('Clinical Significance Transitions\n(hg19 → hg38)', fontsize=12, fontweight='bold')
    
    def _plot_directional_clinical_changes(self, df, ax):
        """Plot directional clinical significance changes"""
        print("4. Creating directional clinical significance changes...")
        
        if len(df) > 0 and 'clin_sig_change' in df.columns:
            # Count all directional changes
            change_counts = df['clin_sig_change'].value_counts()
            
            # Sort by importance: Stable categories first, then directional changes
            stable_changes = [change for change in change_counts.index if change.startswith('STABLE_')]
            directional_changes = [change for change in change_counts.index if '_TO_' in change and not change.startswith('STABLE_')]
            
            # Order by clinical importance
            ordered_changes = stable_changes + directional_changes
            change_counts = change_counts.reindex([change for change in ordered_changes if change in change_counts.index])
            
            if len(change_counts) > 0:
                bars = ax.bar(range(len(change_counts)), change_counts.values,
                                color=self.colors_clinical[:len(change_counts)], alpha=0.8)
                
                ax.set_xlabel('Clinical Significance Change', fontweight='bold')
                ax.set_ylabel('Variant Count', fontweight='bold')
                ax.set_title('Directional Clinical Significance Changes\n(All Transitions)', fontsize=12, fontweight='bold')
                ax.set_xticks(range(len(change_counts)))
                ax.set_xticklabels([change.replace('_', ' ') for change in change_counts.index], 
                                  rotation=45, ha='right', fontsize=8)
                ax.grid(True, alpha=0.3, axis='y')
                
                # Add count labels on bars
                for i, (bar, count) in enumerate(zip(bars, change_counts.values)):
                    ax.text(bar.get_x() + bar.get_width()/2, bar.get_height() + count*0.01,
                              str(count), ha='center', va='bottom', fontweight='bold', fontsize=9)
            else:
                ax.text(0.5, 0.5, 'No clinical significance\nchanges found', ha='center', va='center',
                          fontsize=12, transform=ax.transAxes)
                ax.set_title('Directional Clinical Significance Changes\n(All Transitions)', fontsize=12, fontweight='bold')
        else:
            ax.text(0.5, 0.5, 'No clinical significance\nchange data available', ha='center', va='center',
                      fontsize=12, transform=ax.transAxes)
            ax.set_title('Directional Clinical Significance Changes\n(All Transitions)', fontsize=12, fontweight='bold')
    
    def _categorize_discordance_primary(self, summary):
        """Categorize discordance without numbers for aggregation"""
        if pd.isna(summary) or summary == '':
            return 'Other'
        
        summary_lower = summary.lower()
        
        if 'same transcript consequence' in summary_lower:
            return 'Same Transcript\nConsequence Change'
        elif 'gene annotation' in summary_lower:
            return 'Gene Annotation\nChange'
        elif 'impact level' in summary_lower:
            return 'Impact Level\nChange'
        elif 'clinical significance' in summary_lower:
            return 'Clinical Significance\nChange'
        elif 'sift change' in summary_lower or 'polyphen change' in summary_lower:
            return 'Pathogenicity\nPrediction Change'
        elif 'unmatched consequences' in summary_lower:
            return 'Unmatched\nConsequences'
        elif 'position mismatch' in summary_lower:
            return 'Position\nMismatch'
        elif 'genotype mismatch' in summary_lower:
            return 'Genotype\nMismatch'
        else:
            return 'Other/Multiple'
    
    def _print_summary_statistics(self, df, conn):
        """Print enhanced summary statistics"""
        print("\n=== ENHANCED VISUALIZATION SUMMARY ===")
        print(f"Total discordant variants analyzed: {len(df):,}")
        
        if len(df) > 0:
            # Impact transition summary
            if 'hg19_impact' in df.columns and 'hg38_impact' in df.columns:
                print(f"\nImpact level transitions:")
                impact_changes = df[df['hg19_impact'] != df['hg38_impact']]
                print(f"  Variants with impact changes: {len(impact_changes):,}")
                
            # Clinical significance transition summary
            if 'hg19_clin_sig_normalized' in df.columns and 'hg38_clin_sig_normalized' in df.columns:
                print(f"\nClinical significance transitions:")
                clin_changes = df[df['hg19_clin_sig_normalized'] != df['hg38_clin_sig_normalized']]
                print(f"  Variants with clinical significance changes: {len(clin_changes):,}")
                
                # Show top clinical transitions
                if 'clin_sig_change' in df.columns:
                    top_changes = df['clin_sig_change'].value_counts().head(5)
                    print(f"  Top clinical significance changes:")
                    for change, count in top_changes.items():
                        print(f"    {change.replace('_', ' ')}: {count:,}")
        
        print("✓ Enhanced prioritization visualizations completed")