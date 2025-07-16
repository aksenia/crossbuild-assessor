"""
Variant Processor

Main orchestrator for the analysis pipeline, coordinating VEP analysis,
caching, and clinical scoring.
"""

from .vep_analyzer import VEPAnalyzer
from .scoring_engine import ClinicalScorer
from .cache_manager import CacheManager


class VariantProcessor:
    """Main orchestrator for variant analysis pipeline"""
    
    def __init__(self):
        """Initialize variant processor with analysis components"""
        self.vep_analyzer = VEPAnalyzer()
        self.clinical_scorer = ClinicalScorer()
    
    def process_all_variants(self, conn, cache_file=None, force_recalculate=False):
        """
        Calculate priority scores using clinical evidence-driven analysis with memory optimization and caching
        
        Args:
            conn: Database connection
            cache_file: Path to cache file for VEP analysis results
            force_recalculate: If True, ignore cache and recalculate
            
        Returns:
            pd.DataFrame: Scored variants with priority categories
        """
        print("Calculating variant priority scores with clinical evidence-driven analysis...")
        
        # Initialize cache manager
        cache_manager = CacheManager(cache_file)
        
        # Check for cached results (raw VEP analysis only, no scores)
        if cache_manager.should_use_cache(force_recalculate):
            try:
                cached_vep_analysis = cache_manager.load_cache()
                
                # Calculate scores on-the-fly (not cached)
                result_df = self.clinical_scorer.calculate_scores_from_analysis(cached_vep_analysis)
                return result_df
            except Exception as e:
                print(f"Warning: Could not load cache file ({e}), recalculating...")
        
        # Calculate fresh VEP analysis
        vep_analysis_df = self.vep_analyzer.analyze_all_variants(conn)
        
        # Save VEP analysis to cache (without scores)
        cache_manager.save_cache(vep_analysis_df)
        
        # Calculate scores on-the-fly
        result_df = self.clinical_scorer.calculate_scores_from_analysis(vep_analysis_df)
        return result_df
