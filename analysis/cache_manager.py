"""
Cache Manager

Handles caching of VEP analysis results for faster subsequent runs.
Stores raw VEP analysis only (no scores) to allow score recalibration.
"""

import pandas as pd
from pathlib import Path


class CacheManager:
    """Manages caching of VEP analysis results"""
    
    def __init__(self, cache_file=None):
        """
        Initialize cache manager
        
        Args:
            cache_file: Path to cache file (pickle format)
        """
        self.cache_file = Path(cache_file) if cache_file else None
    
    def should_use_cache(self, force_recalculate=False):
        """
        Determine if cache should be used
        
        Args:
            force_recalculate: If True, ignore cache and recalculate
            
        Returns:
            bool: True if cache should be used
        """
        if force_recalculate:
            return False
        
        if not self.cache_file:
            return False
            
        return self.cache_file.exists()
    
    def load_cache(self):
        """
        Load cached VEP analysis results
        
        Returns:
            pd.DataFrame: Cached VEP analysis data
            
        Raises:
            Exception: If cache cannot be loaded
        """
        if not self.cache_file or not self.cache_file.exists():
            raise ValueError("Cache file does not exist")
        
        print(f"Loading cached VEP analysis from: {self.cache_file}")
        
        try:
            cached_vep_analysis = pd.read_pickle(self.cache_file)
            print(f"Loaded {len(cached_vep_analysis):,} cached VEP analyses")
            return cached_vep_analysis
        except Exception as e:
            raise Exception(f"Could not load cache file: {e}")
    
    def save_cache(self, vep_analysis_df):
        """
        Save VEP analysis results to cache
        
        Args:
            vep_analysis_df: DataFrame with VEP analysis results to cache
        """
        if not self.cache_file:
            return
        
        try:
            vep_analysis_df.to_pickle(self.cache_file)
            print(f"✓ Cached VEP analysis saved to: {self.cache_file}")
        except Exception as e:
            print(f"Warning: Could not save cache file ({e})")
