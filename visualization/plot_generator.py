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
        
        # PLOT 2: CANONICAL HGVSc Match vs Clinical Changes 
        self._plot_canonical_hgvsc_match(df, axes[0,1])
        
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


    def _plot_canonical_hgvsc_match(self, df, ax):
        """Plot 2: Placeholder - MANE-based HGVS Analysis (To be implemented in Phase 2)"""
        print("2. Creating placeholder for MANE-based HGVS Analysis...")
        
        if len(df) > 0 and 'Clinical_Change_Direction' in df.columns:
            # Create simple clinical change distribution as placeholder
            clinical_counts = df['Clinical_Change_Direction'].value_counts()
            
            # Group all stable categories together
            stable_count = sum(count for cat, count in clinical_counts.items() if str(cat).startswith('STABLE_'))
            
            # Create simplified data for placeholder
            plot_data = {'MANE Analysis (Pending)': stable_count}
            
            # Add other categories
            other_categories = {cat: count for cat, count in clinical_counts.items() 
                              if not str(cat).startswith('STABLE_')}
            if other_categories:
                plot_data.update(other_categories)
            
            # Use same color scheme as original
            colors = []
            clinical_colors = self.plot_colors.get('clinical_transitions', {})
            
            for cat in plot_data.keys():
                if cat == 'MANE Analysis (Pending)':
                    colors.append(clinical_colors.get('stable', '#1f77b4'))
                elif '→ PATHOGENIC' in str(cat):
                    colors.append(clinical_colors.get('to_pathogenic', '#d62728'))
                elif 'FROM PATHOGENIC →' in str(cat):
                    colors.append(clinical_colors.get('from_pathogenic', '#ff7f0e'))
                elif '→ BENIGN' in str(cat):
                    colors.append(clinical_colors.get('to_benign', '#2ca02c'))
                elif 'FROM BENIGN →' in str(cat):
                    colors.append(clinical_colors.get('from_benign', '#17becf'))
                else:
                    colors.append(clinical_colors.get('other', '#9467bd'))
            
            # Create simple bar plot
            ax.bar([0], [sum(plot_data.values())], color=colors[0], width=0.7)
            
            ax.set_xlabel('MANE-Based HGVS Analysis', fontweight='bold')
            ax.set_ylabel('Variant Count', fontweight='bold')
            ax.set_title('MANE-Based HGVS Analysis vs Clinical Changes (Placeholder)', 
                        fontsize=12, fontweight='bold')
            ax.set_xticks([0])
            ax.set_xticklabels(['Pending MANE\nImplementation'], rotation=0, ha='center')
            ax.grid(True, alpha=0.3, axis='y')
            
            # Add count label on single bar
            total_count = sum(plot_data.values())
            ax.text(0, total_count + total_count*0.01, f'{int(total_count)}', 
                   ha='center', va='bottom', fontweight='bold')
        else:
            ax.text(0.5, 0.5, 'Placeholder for\nMANE-based HGVS Analysis', ha='center', va='center',
                    fontsize=12, transform=ax.transAxes)
            ax.set_title('MANE-Based HGVS Analysis vs Clinical Changes (Placeholder)', 
                        fontsize=12, fontweight='bold')
        
    def _plot_impact_transitions_grouped(self, df, ax):
        """Plot 2: Impact transitions with grouping (REDESIGNED)"""
        print("2. Creating impact transitions with grouping...")
        
        if len(df) > 0:
            # Group impacts: HIGH, MODERATE, LOW-IMPACT (LOW+MODIFIER)
            df_copy = df.copy()
            df_copy['hg19_impact_grouped'] = df_copy['hg19_impact'].apply(
                lambda x: 'LOW-IMPACT' if x in ['LOW', 'MODIFIER'] else x
            )
            df_copy['hg38_impact_grouped'] = df_copy['hg38_impact'].apply(
                lambda x: 'LOW-IMPACT' if x in ['LOW', 'MODIFIER'] else x
            )
            
            # Categorize transitions
            def categorize_impact_transition(row):
                hg19 = row['hg19_impact_grouped']
                hg38 = row['hg38_impact_grouped']
                
                if hg19 == hg38:
                    return f'Stable {hg19}'
                elif hg38 == 'HIGH':
                    return '→ HIGH'
                elif hg38 == 'MODERATE':
                    return '→ MODERATE'
                elif hg19 == 'HIGH':
                    return 'FROM HIGH →'
                elif hg19 == 'MODERATE':
                    return 'FROM MODERATE →'
                elif hg19 == 'LOW-IMPACT' and hg38 == 'LOW-IMPACT':
                    return 'LOW-IMPACT Transitions'
                else:
                    return 'Other Transitions'
            
            df_copy['impact_transition_category'] = df_copy.apply(categorize_impact_transition, axis=1)
            
            # Order categories
            category_order = [
                'Stable HIGH', 'Stable MODERATE', 'Stable LOW-IMPACT',
                '→ HIGH', '→ MODERATE',
                'FROM HIGH →', 'FROM MODERATE →',
                'LOW-IMPACT Transitions', 'Other Transitions'
            ]
            
            transition_counts = df_copy['impact_transition_category'].value_counts()
            transition_counts = transition_counts.reindex([cat for cat in category_order if cat in transition_counts.index])
            
            if len(transition_counts) > 0:
                # Color coding: stable (blue), gains (red), losses (orange), low-impact (gray)
                colors = []
                for cat in transition_counts.index:
                    if 'Stable' in cat:
                        colors.append('#1f77b4')  # Blue for stable
                    elif '→' in cat and 'FROM' not in cat:
                        colors.append('#d62728')  # Red for clinical gains
                    elif 'FROM' in cat:
                        colors.append('#ff7f0e')  # Orange for clinical losses
                    elif 'LOW-IMPACT' in cat:
                        colors.append('#7f7f7f')  # Gray for low impact
                    else:
                        colors.append('#9467bd')  # Purple for other
                
                bars = ax.bar(range(len(transition_counts)), transition_counts.values,
                             color=colors, alpha=0.8)
                
                ax.set_xlabel('Impact Transition Category', fontweight='bold')
                ax.set_ylabel('Variant Count', fontweight='bold')
                ax.set_title('Impact Transitions (hg19→hg38)\nGrouped by Clinical Significance', fontsize=12, fontweight='bold')
                ax.set_xticks(range(len(transition_counts)))
                ax.set_xticklabels(transition_counts.index, rotation=45, ha='right', fontsize=9)
                ax.grid(True, alpha=0.3, axis='y')
                
                # Add count labels
                for bar in bars:
                    height = bar.get_height()
                    if height > 0:
                        ax.text(bar.get_x() + bar.get_width()/2, height + transition_counts.max()*0.01,
                               f'{int(height)}', ha='center', va='bottom', fontweight='bold', fontsize=9)
            else:
                ax.text(0.5, 0.5, 'No impact transitions\nfound', ha='center', va='center',
                          fontsize=12, transform=ax.transAxes)
                ax.set_title('Impact Transitions (hg19→hg38)\nGrouped by Clinical Significance', fontsize=12, fontweight='bold')
        else:
            ax.text(0.5, 0.5, 'No variants\nfound', ha='center', va='center',
                      fontsize=12, transform=ax.transAxes)
            ax.set_title('Impact Transitions (hg19→hg38)\nGrouped by Clinical Significance', fontsize=12, fontweight='bold')
    
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
            df['discordance_primary'] = df['discordance_summary'].apply(self._categorize_discordance_primary)
            discordance_counts = df['discordance_primary'].value_counts()
            
            # Order categories by clinical importance
            discordance_order = [
                'CRITICAL Clinical\nSignificance Change',
                'Clinical\nSignificance Change', 
                'Same Transcript\nConsequence Change',
                'Disjoint Functional\nConsequences',
                'Impact Level\nChange',
                'Pathogenicity\nPrediction Change',
                'Technical/Annotation\nIssues',
                'Other'
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
        """Categorize discordance with combined low-priority categories for visual clarity"""
        if pd.isna(summary) or summary == '':
            return 'Other'
        
        summary_lower = summary.lower()
        
        # PRIORITY ORDER: Clinical impact first (highest priority)
        if 'critical clinical significance' in summary_lower:
            return 'CRITICAL Clinical\nSignificance Change'
        elif 'clinical significance' in summary_lower:
            return 'Clinical\nSignificance Change'
        elif 'same transcript consequence' in summary_lower:
            return 'Same Transcript\nConsequence Change'
        elif 'disjoint functional consequences' in summary_lower:
            return 'Disjoint Functional\nConsequences'
        elif 'impact level' in summary_lower:
            return 'Impact Level\nChange'
        elif 'sift prediction change' in summary_lower or 'polyphen prediction change' in summary_lower:
            return 'Pathogenicity\nPrediction Change'
        elif any(term in summary_lower for term in [
            'gene symbol changes', 'consequence subset', 'consequence subset relationship',
            'position mismatch', 'genotype mismatch', 'annotation noise',
            'partial consequence overlap', 'same consequence, different transcripts'
        ]):
            return 'Technical/Annotation\nIssues'
        else:
            return 'Technical/Annotation\nIssues' 
    
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