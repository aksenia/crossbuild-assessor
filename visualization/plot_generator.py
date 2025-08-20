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
        self.plot_colors = plot_colors
        self.figure_config = figure_config
        
        # Build-specific colors (colorblind friendly)
        self.colors_builds = ['#1f77b4', '#ff7f0e']  # Blue for hg19, Orange for hg38
    
    def create_all_plots(self, df, conn, output_dir):
        """Create enhanced visualization plots for variant prioritization analysis"""
        print("Creating enhanced prioritization visualizations...")
        
        # Create enhanced visualization with 4 plots
        fig, axes = plt.subplots(
            self.figure_config['main_figure']['rows'], 
            self.figure_config['main_figure']['cols'], 
            figsize=self.figure_config['main_figure']['figsize']
        )
        fig.suptitle('Discrepant Variant Prioritization Analysis', fontsize=16, fontweight='bold')
        
        # PLOT 1: Clinical Evidence Distribution by Build 
        self._plot_clinical_evidence_by_build(df, axes[0,0])

        # PLOT 2: HGVSc Concordance Overview
        self._plot_hgvs_concordance_overview(df, axes[0,1])
        
        # PLOT 3: Clinical Significance Transitions with Grouping 
        self._plot_clinical_significance_transitions_grouped(df, axes[1,0])
        
        # PLOT 4: Primary Discordance Types 
        self._plot_discordance_types(df, axes[1,1])
        
        plt.tight_layout()
        
        # Save plot
        output_path = Path(output_dir) / 'variant_prioritization_plots.png'
        plt.savefig(output_path, dpi=300, bbox_inches='tight')
        plt.show()
        
        # Print summary statistics
        self._print_summary_statistics(df, conn)
    
    def _plot_clinical_evidence_by_build(self, df, ax):
        """Plot 1: Clinical evidence distribution by build (NEW)"""
        print("1. Creating clinical evidence distribution by build...")
        
        if len(df) > 0 and 'hg19_clin_sig_normalized' in df.columns:
            # Get categories present in the data
            all_categories = set(df['hg19_clin_sig_normalized'].unique()) | set(df['hg38_clin_sig_normalized'].unique())
            categories = ['PATHOGENIC', 'BENIGN', 'VUS', 'RISK', 'DRUG_RESPONSE', 'PROTECTIVE', 'OTHER', 'NONE']
            present_categories = [cat for cat in categories if cat in all_categories]
            
            if present_categories:
                # Count occurrences in each build
                hg19_counts = df['hg19_clin_sig_normalized'].value_counts()
                hg38_counts = df['hg38_clin_sig_normalized'].value_counts()
                
                x_positions = np.arange(len(present_categories))
                width = 0.35
                
                hg19_values = [hg19_counts.get(cat, 0) for cat in present_categories]
                hg38_values = [hg38_counts.get(cat, 0) for cat in present_categories]
                
                bars1 = ax.bar(x_positions - width/2, hg19_values, width, 
                              label='hg19', color=self.colors_builds[0], alpha=0.8)
                bars2 = ax.bar(x_positions + width/2, hg38_values, width,
                              label='hg38', color=self.colors_builds[1], alpha=0.8)
                
                ax.set_xlabel('Clinical Significance Category', fontweight='bold')
                ax.set_ylabel('Variant Count', fontweight='bold')
                ax.set_title('Clinical Significance Distribution by Build', fontsize=12, fontweight='bold')
                ax.set_xticks(x_positions)
                ax.set_xticklabels(present_categories, rotation=45, ha='right', fontsize=9)
                ax.legend()
                ax.grid(True, alpha=0.3, axis='y')
                
                # Add count labels on bars with overflow protection
                for bars in [bars1, bars2]:
                    for bar in bars:
                        height = bar.get_height()
                        if height > 0:
                            # Check if there's space for text above the bar
                            max_height = ax.get_ylim()[1]
                            if height < max_height * 0.9:  # If bar is less than 90% of max, put text above
                                ax.text(bar.get_x() + bar.get_width()/2, height + max_height*0.01,
                                       f'{int(height)}', ha='center', va='bottom', fontweight='bold', fontsize=8)
                            else:  # Put text inside the bar
                                ax.text(bar.get_x() + bar.get_width()/2, height/2,
                                       f'{int(height)}', ha='center', va='center', fontweight='bold', 
                                       fontsize=8, color='white')
            else:
                ax.text(0.5, 0.5, 'No clinical significance\ndata available', ha='center', va='center',
                          fontsize=12, transform=ax.transAxes)
                ax.set_title('Clinical Significance Distribution by Build', fontsize=12, fontweight='bold')
        else:
            ax.text(0.5, 0.5, 'No clinical significance\ndata available', ha='center', va='center',
                      fontsize=12, transform=ax.transAxes)
            ax.set_title('Clinical Significance Distribution by Build', fontsize=12, fontweight='bold')

    def _plot_hgvs_concordance_overview(self, df, ax):
        """Plot 2: Simple HGVS Concordance Overview (HGVSc vs HGVSp)"""
        print("2. Creating HGVS Concordance Overview...")
        
        if len(df) > 0 and 'priority_hgvsc_concordance' in df.columns and 'priority_hgvsp_concordance' in df.columns:
            
            # Calculate concordance rates for HGVSc
            hgvsc_total = len(df)
            hgvsc_match = (df['priority_hgvsc_concordance'] == 'Match').sum()
            hgvsc_mismatch = (df['priority_hgvsc_concordance'] == 'Mismatch').sum()
            hgvsc_no_analysis = (df['priority_hgvsc_concordance'] == 'No_Analysis').sum()
            
            # Calculate concordance rates for HGVSp
            hgvsp_total = len(df)
            hgvsp_match = (df['priority_hgvsp_concordance'] == 'Match').sum()
            hgvsp_mismatch = (df['priority_hgvsp_concordance'] == 'Mismatch').sum()
            hgvsp_no_analysis = (df['priority_hgvsp_concordance'] == 'No_Analysis').sum()
            
            # Convert to percentages
            hgvsc_match_pct = (hgvsc_match / hgvsc_total) * 100
            hgvsc_mismatch_pct = (hgvsc_mismatch / hgvsc_total) * 100
            hgvsc_no_analysis_pct = (hgvsc_no_analysis / hgvsc_total) * 100
            
            hgvsp_match_pct = (hgvsp_match / hgvsp_total) * 100
            hgvsp_mismatch_pct = (hgvsp_mismatch / hgvsp_total) * 100
            hgvsp_no_analysis_pct = (hgvsp_no_analysis / hgvsp_total) * 100
            
            # Data for plotting
            categories = ['Match', 'Mismatch', 'No Analysis']
            hgvsc_values = [hgvsc_match_pct, hgvsc_mismatch_pct, hgvsc_no_analysis_pct]
            hgvsp_values = [hgvsp_match_pct, hgvsp_mismatch_pct, hgvsp_no_analysis_pct]
            
            # Use colors from config
            priority_colors = self.plot_colors.get('priority', ['#d62728', '#ff7f0e', '#2ca02c', '#1f77b4'])
            colors = [priority_colors[2], priority_colors[0], '#7f7f7f']  # Green, Red, Gray
            
            # Create grouped bar plot
            x = np.arange(len(categories))
            width = 0.35
            
            bars1 = ax.bar(x - width/2, hgvsc_values, width, 
                        label=f'HGVSc (n={hgvsc_total})', color=colors, alpha=0.8)
            bars2 = ax.bar(x + width/2, hgvsp_values, width,
                        label=f'HGVSp (n={hgvsp_total})', color=colors, alpha=0.6)
            
            # Customize plot using config
            ax.set_xlabel('HGVS Concordance Status', fontweight='bold')
            ax.set_ylabel('Percentage of Variants (%)', fontweight='bold')
            ax.set_title('HGVS Concordance Overview\n(Priority Transcript Analysis)', 
                        fontsize=12, fontweight='bold')
            ax.set_xticks(x)
            ax.set_xticklabels(categories)
            ax.legend()
            ax.grid(True, alpha=0.3, axis='y')
            ax.set_ylim(0, 100)
            
            # Add percentage and count labels on bars
            for i, (bar1, bar2) in enumerate(zip(bars1, bars2)):
                # HGVSc labels
                if hgvsc_values[i] > 3:  # Only show if bar is big enough
                    count_c = [hgvsc_match, hgvsc_mismatch, hgvsc_no_analysis][i]
                    ax.text(bar1.get_x() + bar1.get_width()/2, bar1.get_height() + 1,
                        f'{hgvsc_values[i]:.1f}%\n({count_c})', ha='center', va='bottom',
                        fontweight='bold', fontsize=8)
                
                # HGVSp labels  
                if hgvsp_values[i] > 3:  # Only show if bar is big enough
                    count_p = [hgvsp_match, hgvsp_mismatch, hgvsp_no_analysis][i]
                    ax.text(bar2.get_x() + bar2.get_width()/2, bar2.get_height() + 1,
                        f'{hgvsp_values[i]:.1f}%\n({count_p})', ha='center', va='bottom',
                        fontweight='bold', fontsize=8)
            
            # Add summary text box
            summary_text = f"Total Variants: {hgvsc_total:,}\nHGVSc Issues: {hgvsc_mismatch+hgvsc_no_analysis} ({((hgvsc_mismatch+hgvsc_no_analysis)/hgvsc_total*100):.1f}%)\nHGVSp Issues: {hgvsp_mismatch+hgvsp_no_analysis} ({((hgvsp_mismatch+hgvsp_no_analysis)/hgvsp_total*100):.1f}%)"
            ax.text(0.02, 0.98, summary_text, transform=ax.transAxes, fontsize=9,
                verticalalignment='top', bbox=dict(boxstyle='round', facecolor='wheat', alpha=0.8))
            
        else:
            ax.text(0.5, 0.5, 'HGVS concordance data\nnot available', 
                ha='center', va='center', fontsize=12, transform=ax.transAxes)
            ax.set_title('HGVS Concordance Overview\n(Priority Transcript Analysis)', 
                        fontsize=12, fontweight='bold')
                
    def _plot_clinical_significance_transitions_grouped(self, df, ax):
        """Plot 3: Clinical significance transitions with grouping (REDESIGNED)"""
        print("3. Creating clinical significance transitions with grouping...")
        
        if len(df) > 0 and 'hg19_clin_sig_normalized' in df.columns:
            # Categorize transitions
            def categorize_clin_sig_transition(row):
                hg19 = row['hg19_clin_sig_normalized']
                hg38 = row['hg38_clin_sig_normalized']
                
                if hg19 == hg38:
                    return f'Stable {hg19}'
                elif hg38 == 'PATHOGENIC':
                    return '→ PATHOGENIC'
                elif hg19 == 'PATHOGENIC':
                    return 'FROM PATHOGENIC →'
                elif hg38 == 'BENIGN':
                    return '→ BENIGN'
                elif hg19 == 'BENIGN':
                    return 'FROM BENIGN →'
                else:
                    return 'Other Transitions'
            
            df_copy = df.copy()
            df_copy['clin_sig_transition_category'] = df_copy.apply(categorize_clin_sig_transition, axis=1)
            
            # Get all stable categories that exist in the data
            stable_categories = [cat for cat in df_copy['clin_sig_transition_category'].unique() if cat.startswith('Stable')]
            stable_categories.sort()  # Sort alphabetically
            
            # Order categories: stable first, then critical transitions, then other
            category_order = stable_categories + [
                '→ PATHOGENIC', 'FROM PATHOGENIC →',
                '→ BENIGN', 'FROM BENIGN →',
                'Other Transitions'
            ]
            
            transition_counts = df_copy['clin_sig_transition_category'].value_counts()
            transition_counts = transition_counts.reindex([cat for cat in category_order if cat in transition_counts.index])
            
            if len(transition_counts) > 0:
                # Color coding: stable (blue), pathogenic (red), benign (green), other (gray)
                colors = []
                for cat in transition_counts.index:
                    if 'Stable' in cat:
                        colors.append('#1f77b4')  # Blue for stable
                    elif 'PATHOGENIC' in cat:
                        colors.append('#d62728')  # Red for pathogenic transitions
                    elif 'BENIGN' in cat:
                        colors.append('#2ca02c')  # Green for benign transitions
                    else:
                        colors.append('#7f7f7f')  # Gray for other
                
                bars = ax.bar(range(len(transition_counts)), transition_counts.values,
                             color=colors, alpha=0.8)
                
                ax.set_xlabel('Clinical Significance Transition', fontweight='bold')
                ax.set_ylabel('Variant Count', fontweight='bold')
                ax.set_title('Clinical Significance Transitions (hg19→hg38)', fontsize=12, fontweight='bold')
                ax.set_xticks(range(len(transition_counts)))
                ax.set_xticklabels(transition_counts.index, rotation=45, ha='right', fontsize=8)
                ax.grid(True, alpha=0.3, axis='y')
                
                # Add count labels
                for bar in bars:
                    height = bar.get_height()
                    if height > 0:
                        ax.text(bar.get_x() + bar.get_width()/2, height + transition_counts.max()*0.01,
                               f'{int(height)}', ha='center', va='bottom', fontweight='bold', fontsize=9)
            else:
                ax.text(0.5, 0.5, 'No clinical significance\ntransitions found', ha='center', va='center',
                          fontsize=12, transform=ax.transAxes)
                ax.set_title('Clinical Significance Transitions (hg19→hg38)', fontsize=12, fontweight='bold')
        else:
            ax.text(0.5, 0.5, 'No clinical significance\ndata available', ha='center', va='center',
                      fontsize=12, transform=ax.transAxes)
            ax.set_title('Clinical Significance Transitions (hg19→hg38)', fontsize=12, fontweight='bold')
    
    def _plot_discordance_types(self, df, ax):
        """Plot 4: Primary discordance types"""
        print("4. Creating primary discordance types visualization...")
        
        if len(df) > 0:
            df['discordance_primary'] = df['score_breakdown'].apply(self._categorize_discordance_primary)
            discordance_counts = df['discordance_primary'].value_counts()
            
            # Order exactly following our BASE_SCORES hierarchy (high to low points)
            discordance_order = [
                'HGVS Changes',              # 100/50 points - CRITICAL
                'Major Clinical\nChanges',   # 90 points - CRITICAL (P↔B only)
                'Transcript\nIssues',        # 60 points - MODERATE  
                'Other Clinical\nChanges',   # 40/25/20 points - MODERATE/LOW
                'Consequence\nDifference',   # 35/20 points - MODERATE/LOW
                'Technical\nLiftover Issues', # 20/15/10 points - LOW
                'Annotation\nChanges',       # 15 points - LOW
                'Prediction\nChanges',       # 10 points - LOW
                'Other Issues',
                'No Issues'
            ]
            discordance_counts = discordance_counts.reindex([cat for cat in discordance_order if cat in discordance_counts.index])
            
            bars = ax.bar(range(len(discordance_counts)), discordance_counts.values,
                            color=self.colors_priority[:len(discordance_counts)], alpha=0.8)
            
            ax.set_xlabel('Primary Discordance Type', fontweight='bold')
            ax.set_ylabel('Variant Count', fontweight='bold')
            ax.set_title('Variants by Primary Discordance Type\n(Functional Changes)', fontsize=12, fontweight='bold')
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
            ax.set_title('Variants by Primary Discordance Type\n(Functional Changes)', fontsize=12, fontweight='bold')
    
    def _categorize_discordance_primary(self, summary):
        """Categorize discordance based on score breakdown from new scoring system"""
        if pd.isna(summary) or summary == '':
            return 'No Issues'
        
        summary_lower = str(summary).lower()
        
        # PRIORITY ORDER: Exactly following our BASE_SCORES hierarchy
        if 'hgvsc mismatch' in summary_lower or 'hgvsp mismatch' in summary_lower:
            return 'HGVS Changes'                 # 100/50 points - CRITICAL
        elif 'pathogenic↔benign change' in summary_lower:
            return 'Major Clinical\nChanges'      # 90 points - CRITICAL (P↔B, LP↔B)
        elif any(x in summary_lower for x in ['pathogenic↔vus change', 'vus↔benign change', 'minor clinical change']):
            return 'Other Clinical\nChanges'      # 40/25/20 points - MODERATE/LOW
        elif 'transcript mismatch' in summary_lower:
            return 'Transcript\nIssues'           # 60 points - MODERATE
        elif 'serious consequence difference' in summary_lower or 'minor consequence difference' in summary_lower:
            return 'Consequence\nDifference'      # 35/20 points - MODERATE/LOW
        elif any(x in summary_lower for x in ['position mismatch', 'genotype mismatch', 'ref/alt swap']):
            return 'Technical\nLiftover Issues'   # 20/15/10 points - LOW
        elif 'sift change' in summary_lower or 'polyphen change' in summary_lower:
            return 'Prediction\nChanges'          # 10 points - LOW
        elif 'gene changes' in summary_lower or 'impact changes' in summary_lower:
            return 'Annotation\nChanges'          # 15 points - LOW
        else:
            return 'Other Issues'
    
    def _print_summary_statistics(self, df, conn):
        """Print enhanced summary statistics"""
        print("\n=== VISUALIZATION SUMMARY ===")
        print(f"Total discordant variants analyzed: {len(df):,}")
        
        if len(df) > 0:
            # Clinical evidence summary
            if 'hg19_clin_sig_normalized' in df.columns:
                print(f"\nClinical evidence distribution:")
                hg19_with_clin = (df['hg19_clin_sig_normalized'] != 'NONE').sum()
                hg38_with_clin = (df['hg38_clin_sig_normalized'] != 'NONE').sum()
                print(f"  hg19 variants with clinical data: {hg19_with_clin:,}")
                print(f"  hg38 variants with clinical data: {hg38_with_clin:,}")
                
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
                
                # Critical transitions
                pathogenic_transitions = df[
                    (df['hg19_clin_sig_normalized'] == 'PATHOGENIC') | 
                    (df['hg38_clin_sig_normalized'] == 'PATHOGENIC')
                ]
                pathogenic_changes = pathogenic_transitions[
                    pathogenic_transitions['hg19_clin_sig_normalized'] != pathogenic_transitions['hg38_clin_sig_normalized']
                ]
                print(f"  PATHOGENIC-related transitions: {len(pathogenic_changes):,}")
        
        print("✓ Prioritization visualizations completed")