"""
Configuration management for scMetabolism.
"""

import os
import json
import yaml
from typing import Dict, Any, Optional
from pathlib import Path
import warnings


class Config:
    """Configuration manager for scMetabolism."""
    
    def __init__(self, config_file: Optional[str] = None):
        """
        Initialize configuration.
        
        Parameters:
        -----------
        config_file : str, optional
            Path to configuration file
        """
        self.config_dir = Path.home() / ".scmetabolism"
        self.config_file = config_file or (self.config_dir / "config.yaml")
        self.cache_dir = self.config_dir / "cache"
        
        # Create directories
        self.config_dir.mkdir(exist_ok=True)
        self.cache_dir.mkdir(exist_ok=True)
        
        # Default configuration
        self.defaults = {
            "analysis": {
                "default_method": "aucell",
                "default_n_cores": 2,
                "default_imputation": False,
                "min_genes_per_pathway": 5,
                "max_genes_per_pathway": 500
            },
            "visualization": {
                "default_figsize": [8, 6],
                "default_dpi": 300,
                "default_colormap": "viridis",
                "save_format": "png"
            },
            "data": {
                "cache_dir": str(self.cache_dir),
                "auto_download": True,
                "download_timeout": 300
            },
            "go_analysis": {
                "default_organism": "human",
                "evidence_codes": ["EXP", "IDA", "IPI", "IMP", "IGI", "IEP"],
                "exclude_electronic": True,
                "min_genes": 5,
                "max_genes": 500
            },
            "performance": {
                "use_sparse_matrices": True,
                "batch_size": 5000,
                "memory_limit_gb": 8,
                "enable_parallel": True
            }
        }
        
        # Load configuration
        self.config = self.load_config()
    
    def load_config(self) -> Dict[str, Any]:
        """Load configuration from file."""
        if self.config_file.exists():
            try:
                with open(self.config_file, 'r') as f:
                    if self.config_file.suffix.lower() == '.json':
                        user_config = json.load(f)
                    else:
                        user_config = yaml.safe_load(f)
                
                # Merge with defaults
                config = self.defaults.copy()
                self._deep_update(config, user_config)
                return config
                
            except Exception as e:
                warnings.warn(f"Error loading config file: {e}. Using defaults.")
                return self.defaults.copy()
        else:
            return self.defaults.copy()
    
    def save_config(self):
        """Save current configuration to file."""
        try:
            with open(self.config_file, 'w') as f:
                if self.config_file.suffix.lower() == '.json':
                    json.dump(self.config, f, indent=2)
                else:
                    yaml.dump(self.config, f, default_flow_style=False)
        except Exception as e:
            warnings.warn(f"Error saving config file: {e}")
    
    def get(self, key: str, default: Any = None) -> Any:
        """Get configuration value using dot notation."""
        keys = key.split('.')
        value = self.config
        
        for k in keys:
            if isinstance(value, dict) and k in value:
                value = value[k]
            else:
                return default
        
        return value
    
    def set(self, key: str, value: Any):
        """Set configuration value using dot notation."""
        keys = key.split('.')
        config = self.config
        
        for k in keys[:-1]:
            if k not in config:
                config[k] = {}
            config = config[k]
        
        config[keys[-1]] = value
    
    def reset_to_defaults(self):
        """Reset configuration to defaults."""
        self.config = self.defaults.copy()
    
    def _deep_update(self, base_dict: Dict, update_dict: Dict):
        """Deep update dictionary."""
        for key, value in update_dict.items():
            if key in base_dict and isinstance(base_dict[key], dict) and isinstance(value, dict):
                self._deep_update(base_dict[key], value)
            else:
                base_dict[key] = value
    
    def create_user_config_template(self):
        """Create a user configuration template."""
        template_file = self.config_dir / "config_template.yaml"
        
        template = {
            "# scMetabolism Configuration": None,
            "analysis": {
                "default_method": "aucell  # Options: aucell, ssgsea, gsva",
                "default_n_cores": 2,
                "default_imputation": False
            },
            "visualization": {
                "default_figsize": [8, 6],
                "default_colormap": "viridis"
            },
            "data": {
                "cache_dir": str(self.cache_dir),
                "auto_download": True
            }
        }
        
        with open(template_file, 'w') as f:
            yaml.dump(template, f, default_flow_style=False)
        
        print(f"Configuration template created at: {template_file}")


# Global configuration instance
_config = None

def get_config() -> Config:
    """Get global configuration instance."""
    global _config
    if _config is None:
        _config = Config()
    return _config

def set_config_value(key: str, value: Any):
    """Set configuration value globally."""
    config = get_config()
    config.set(key, value)

def get_config_value(key: str, default: Any = None) -> Any:
    """Get configuration value globally."""
    config = get_config()
    return config.get(key, default)


class AnalysisSettings:
    """Analysis-specific settings manager."""
    
    def __init__(self):
        self.config = get_config()
    
    @property
    def default_method(self) -> str:
        return self.config.get('analysis.default_method', 'aucell')
    
    @property
    def default_n_cores(self) -> int:
        return self.config.get('analysis.default_n_cores', 2)
    
    @property
    def default_imputation(self) -> bool:
        return self.config.get('analysis.default_imputation', False)
    
    @property
    def min_genes_per_pathway(self) -> int:
        return self.config.get('analysis.min_genes_per_pathway', 5)
    
    @property
    def max_genes_per_pathway(self) -> int:
        return self.config.get('analysis.max_genes_per_pathway', 500)


class VisualizationSettings:
    """Visualization-specific settings manager."""
    
    def __init__(self):
        self.config = get_config()
    
    @property
    def default_figsize(self) -> list:
        return self.config.get('visualization.default_figsize', [8, 6])
    
    @property
    def default_dpi(self) -> int:
        return self.config.get('visualization.default_dpi', 300)
    
    @property
    def default_colormap(self) -> str:
        return self.config.get('visualization.default_colormap', 'viridis')
    
    @property
    def save_format(self) -> str:
        return self.config.get('visualization.save_format', 'png')


class DataSettings:
    """Data-specific settings manager."""
    
    def __init__(self):
        self.config = get_config()
    
    @property
    def cache_dir(self) -> str:
        return self.config.get('data.cache_dir', str(Path.home() / ".scmetabolism" / "cache"))
    
    @property
    def auto_download(self) -> bool:
        return self.config.get('data.auto_download', True)
    
    @property
    def download_timeout(self) -> int:
        return self.config.get('data.download_timeout', 300)