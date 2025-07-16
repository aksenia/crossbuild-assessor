"""
Prioritization Plot Generator

Creates comprehensive 4-plot visualization for variant prioritization analysis
with clinical evidence categorization and discordance type breakdown.
"""

import matplotlib.pyplot as plt
import pandas as pd
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
        
        # PLOT 1: Priority categories bar chart
        self._plot_priority_categories(df, conn, axes[0,0])
        
        # PLOT 2: Priority Categories × Clinical Evidence
        self._plot_clinical_evidence(df, axes[0,1])
        
        # PLOT 3: Primary Discordance Types
        self._plot_discordance_types(df, axes[1,0])
        
        # PLOT 4: Clinical Evidence Transitions
        self._plot_clinical_transitions(df, axes[1,1])
        
        plt.tight_layout()
        
        # Save plot
        output_path = Path(output_dir) / 'variant_prioritization_plots.png'
        plt.savefig(output_path, dpi=300, bbox_inches='tight')
        plt.show()
        
        # Print summary statistics
        self._print_summary_statistics(df, conn)
    
    def _plot_priority_categories(self, df, conn, ax):
        """Plot priority categories with concordant variants included"""
        print("1. Creating priority categories bar chart...")
        
        # Calculate total variants including concordant ones
        total_query = """SELECT COUNT(*) as total FROM comparison"""
        total_result = pd.read_sql_query(total_query, conn)
        total_variants = total_result.iloc[0]['total']
        
        category_counts = df['priority_category'].value_counts()
        
        # Add actual concordant variants count
        discordant_count = len(df)
        concordant_count = max(0, total_variants - discordant_count)
        if concordant_count > 0:
            category_counts['CONCORDANT'] = concordant_count
        
        # Order categories by priority
        category_order = ['CRITICAL', 'HIGH', 'MODERATE', 'INVESTIGATE', 'LOW', 'CONCORDANT']
        category_counts = category_counts.reindex([cat for cat in category_order if cat in category_counts.index])
        
        # Create bar chart with percentages
        bars = ax.bar(range(len(category_counts)), category_counts.values,
                      color=self.colors_priority[:len(category_counts)], alpha=0.8)
        
        # Add percentage labels on bars
        for i, (bar, count) in enumerate(zip(bars, category_counts.values)):
            pct = count / category_counts.sum() * 100
            ax.text(bar.get_x() + bar.get_width()/2, bar.get_height() + count*0.01,
                    f'{pct:.1f}%\n({count:,})', ha='center', va='bottom',
                    fontweight='bold', fontsize=9)
        
        ax.set_xticks(range(len(category_counts)))
        ax.set_xticklabels(category_counts.index, rotation=45, ha='right', fontsize=10)
        ax.set_ylabel('Variant Count', fontweight='bold')
        ax.set_title('Variant Priority Categories\n(Clinical Review Levels)', fontsize=12, fontweight='bold')
        ax.grid(True, alpha=0.3, axis='y')
    
    def _plot_clinical_evidence(self, df, ax):
        """Plot priority categories by clinical evidence with intelligent gap filling"""
        print("2. Creating priority categories × clinical evidence...")
        
        if len(df) > 0:
            df['clinical_evidence'] = df.apply(self._get_clinical_evidence_color, axis=1)
            
            # Create contingency table
            evidence_priority = pd.crosstab(df['priority_category'], df['clinical_evidence'])
            category_order = ['CRITICAL', 'HIGH', 'MODERATE', 'INVESTIGATE', 'LOW']
            evidence_priority = evidence_priority.reindex(
                index=[cat for cat in category_order if cat in evidence_priority.index],
                fill_value=0
            )
            
            # Create stacked bar chart
            evidence_priority.plot(kind='bar', ax=ax, stacked=True,
                                  color=self.colors_clinical[:len(evidence_priority.columns)], alpha=0.8)
            
            ax.set_title('Priority Categories by Clinical Evidence\n(Intelligent Gap Filling)',
                           fontsize=12, fontweight='bold')
            ax.set_xlabel('Priority Category', fontweight='bold')
            ax.set_ylabel('Variant Count', fontweight='bold')
            ax.legend(title='Clinical Evidence', bbox_to_anchor=(1.05, 1), loc='upper left', fontsize=8)
            ax.tick_params(axis='x', rotation=45)
            ax.grid(True, alpha=0.3, axis='y')
        else:
            ax.text(0.5, 0.5, 'No variants\nfound', ha='center', va='center',
                      fontsize=12, transform=ax.transAxes)
            ax.set_title('Priority Categories by Clinical Evidence\n(Intelligent Gap Filling)',
                           fontsize=12, fontweight='bold')
    
    def _plot_discordance_types(self, df, ax):
        """Plot primary discordance types (aggregated, no numbers)"""
        print("3. Creating primary discordance types visualization...")
        
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
    
    def _plot_clinical_transitions(self, df, ax):
        """Plot clinical evidence transitions"""
        print("4. Creating clinical evidence transitions...")
        
        if len(df) > 0:
            df['clinical_transitions'] = df.apply(self._categorize_clinical_transitions, axis=1)
            transition_counts = df['clinical_transitions'].value_counts()
            
            # Order categories by clinical concern level
            transition_order = [
                'Pathogenic → Benign',
                'Benign → Pathogenic',
                'Stable Pathogenic',
                'Prediction Changes',
                'High Impact Only',
                'Stable Benign',
                'Other/Unknown'
            ]
            
            transition_counts = transition_counts.reindex([cat for cat in transition_order if cat in transition_counts.index])
            
            bars = ax.bar(range(len(transition_counts)), transition_counts.values,
                            color=self.colors_clinical[:len(transition_counts)], alpha=0.8)
            
            ax.set_xlabel('Clinical Evidence Transition', fontweight='bold')
            ax.set_ylabel('Variant Count', fontweight='bold')
            ax.set_title('Clinical Evidence Transitions\n(Dynamic vs Static Evidence)',
                           fontsize=12, fontweight='bold')
            ax.set_xticks(range(len(transition_counts)))
            ax.set_xticklabels(transition_counts.index, rotation=45, ha='right', fontsize=9)
            ax.grid(True, alpha=0.3, axis='y')
            
            # Add count labels on bars
            for i, (bar, count) in enumerate(zip(bars, transition_counts.values)):
                ax.text(bar.get_x() + bar.get_width()/2, bar.get_height() + count*0.01,
                          str(count), ha='center', va='bottom', fontweight='bold', fontsize=10)
        else:
            ax.text(0.5, 0.5, 'No variants\nfound', ha='center', va='center',
                      fontsize=12, transform=ax.transAxes)
            ax.set_title('Clinical Evidence Transitions\n(Dynamic vs Static Evidence)',
                           fontsize=12, fontweight='bold')
    
    def _get_clinical_evidence_color(self, row):
        """Smart gap filling for clinical evidence"""
        # Priority 1: Clinical significance
        if row['clin_sig_change']:
            return 'Clinical Sig Change'
        elif row['hg19_is_pathogenic'] or row['hg38_is_pathogenic']:
            return 'Pathogenic'
        elif row['hg19_is_benign'] or row['hg38_is_benign']:
            return 'Benign'
        
        # Priority 2: SIFT + PolyPhen (combined)
        elif row['sift_change'] or row['polyphen_change']:
            return 'Prediction Change'
        elif ('deleterious' in str(row['hg19_sift']).lower() or 'deleterious' in str(row['hg38_sift']).lower() or
              'probably_damaging' in str(row['hg19_polyphen']).lower() or 'probably_damaging' in str(row['hg38_polyphen']).lower()):
            return 'Pred. Pathogenic'
        elif ('tolerated' in str(row['hg19_sift']).lower() or 'tolerated' in str(row['hg38_sift']).lower() or
              'benign' in str(row['hg19_polyphen']).lower() or 'benign' in str(row['hg38_polyphen']).lower()):
            return 'Pred. Benign'
        
        # Priority 3: Impact level fallback
        elif row['hg19_impact'] == 'HIGH' or row['hg38_impact'] == 'HIGH':
            return 'High Impact'
        else:
            return 'No Evidence'
    
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
    
    def _categorize_clinical_transitions(self, row):
        """Categorize variants by clinical evidence transitions"""
        # Check for clinical significance transitions
        if row['clin_sig_change'] == 'PATHOGENIC_TO_BENIGN':
            return 'Pathogenic → Benign'
        elif row['clin_sig_change'] == 'BENIGN_TO_PATHOGENIC':
            return 'Benign → Pathogenic'
        elif row['hg19_is_pathogenic'] and row['hg38_is_pathogenic']:
            return 'Stable Pathogenic'
        elif row['hg19_is_benign'] and row['hg38_is_benign']:
            return 'Stable Benign'
        elif row['sift_change'] or row['polyphen_change']:
            return 'Prediction Changes'
        elif row['hg19_impact'] == 'HIGH' or row['hg38_impact'] == 'HIGH':
            return 'High Impact Only'
        else:
            return 'Other/Unknown'
    
    def _print_summary_statistics(self, df, conn):
        """Print enhanced summary statistics"""
        print("\n=== ENHANCED VISUALIZATION SUMMARY ===")
        
        # Calculate total variants including concordant ones
        total_query = """SELECT COUNT(*) as total FROM comparison"""
        total_result = pd.read_sql_query(total_query, conn)
        total_variants = total_result.iloc[0]['total']
        
        print(f"Total variants in database: {total_variants:,}")
        print(f"Discordant variants analyzed: {len(df):,}")
        
        discordant_count = len(df)
        concordant_count = max(0, total_variants - discordant_count)
        print(f"Concordant variants: {concordant_count:,}")
        
        if len(df) > 0:
            category_counts = df['priority_category'].value_counts()
            if concordant_count > 0:
                category_counts['CONCORDANT'] = concordant_count
            
            print("\nPriority category breakdown:")
            for category, count in category_counts.items():
                pct = count / category_counts.sum() * 100
                print(f"  {category}: {count:,} ({pct:.1f}%)")
            
            if 'clinical_transitions' in df.columns:
                transition_counts = df['clinical_transitions'].value_counts()
                print(f"\nClinical evidence transitions:")
                for transition, count in transition_counts.items():
                    pct = count / len(df) * 100
                    print(f"  {transition}: {count:,} ({pct:.1f}%)")
            
            if 'discordance_primary' in df.columns:
                discordance_counts = df['discordance_primary'].value_counts()
                print(f"\nPrimary discordance types:")
                for disc_type, count in discordance_counts.items():
                    pct = count / len(df) * 100
                    print(f"  {disc_type.replace(chr(10), ' ')}: {count:,} ({pct:.1f}%)")
        
        print("✓ Enhanced prioritization visualizations completed")
