"""
lotus.config.settings — Configuration Management
===============================================

Pure Python configuration system for LOTUS.
No external dependencies beyond Python stdlib.
"""

import os
import json
from typing import Dict, Any, Optional, List
from decimal import Decimal

class LotusConfig:
    """
    Configuration manager for LOTUS framework.
    """

    # Default configuration
    DEFAULT_CONFIG = {
        # Numerical precision settings
        'precision': {
            'decimal_places': 50,  # For Decimal arithmetic
            'float_tolerance': 1e-15,  # For floating point comparisons
            'max_iterations': 1000,  # For iterative calculations
        },

        # Caching settings
        'cache': {
            'enabled': True,
            'max_size': 5000,  # Maximum cache entries
            'ttl_seconds': 3600,  # Time-to-live (1 hour)
            'auto_clear_threshold': 0.9,  # Clear when 90% full
        },

        # Spectral computation settings
        'spectral': {
            'zeta_terms': 1000,  # Terms in zeta function series
            'eta_terms': 500,  # Terms in eta invariant calculation
            'integration_points': 1000,  # Points for numerical integration
            'eigenvalue_precision': 1e-12,  # Eigenvalue convergence
        },

        # Testing settings
        'testing': {
            'tolerance_percent': 1.0,  # Default prediction tolerance
            'max_test_time': 300,  # Maximum test execution time (seconds)
            'verbose_output': False,  # Verbose test output
        },

        # Visualization settings
        'visualization': {
            'plot_width': 80,  # ASCII plot width
            'plot_height': 24,  # ASCII plot height
            'default_title': 'LOTUS Plot',  # Default plot title
        },

        # Documentation settings
        'documentation': {
            'auto_generate': True,  # Auto-generate docs
            'include_private': False,  # Include private members
            'output_format': 'markdown',  # Output format
        },

        # Framework settings
        'framework': {
            'debug_mode': False,  # Enable debug output
            'log_level': 'INFO',  # Logging level
            'thread_safe': True,  # Thread safety
        }
    }

    def __init__(self, config_file: str = None):
        """
        Initialize configuration.

        Args:
            config_file: Path to configuration file (optional)
        """
        self.config = self.DEFAULT_CONFIG.copy()
        self.config_file = config_file or self._get_default_config_path()

        # Load user configuration
        self.load_config()

    def _get_default_config_path(self) -> str:
        """Get default configuration file path."""
        # Try to find lotus package directory
        try:
            import lotus
            package_dir = os.path.dirname(lotus.__file__)
            return os.path.join(package_dir, 'config.json')
        except ImportError:
            # Fallback to current directory
            return 'lotus_config.json'

    def load_config(self) -> None:
        """
        Load configuration from file.
        """
        if os.path.exists(self.config_file):
            try:
                with open(self.config_file, 'r') as f:
                    user_config = json.load(f)
                self._merge_config(user_config)
            except (json.JSONDecodeError, IOError) as e:
                print(f"Warning: Could not load config file {self.config_file}: {e}")
                print("Using default configuration.")

    def save_config(self) -> None:
        """
        Save current configuration to file.
        """
        try:
            # Create directory if it doesn't exist
            os.makedirs(os.path.dirname(self.config_file), exist_ok=True)

            with open(self.config_file, 'w') as f:
                json.dump(self.config, f, indent=2)
        except IOError as e:
            print(f"Warning: Could not save config file {self.config_file}: {e}")

    def _merge_config(self, user_config: Dict[str, Any]) -> None:
        """
        Merge user configuration with defaults.

        Args:
            user_config: User configuration dictionary
        """
        def merge_dict(target: Dict[str, Any], source: Dict[str, Any]) -> None:
            for key, value in source.items():
                if key in target and isinstance(target[key], dict) and isinstance(value, dict):
                    merge_dict(target[key], value)
                else:
                    target[key] = value

        merge_dict(self.config, user_config)

    def get(self, key: str, default: Any = None) -> Any:
        """
        Get configuration value.

        Args:
            key: Dot-separated configuration key (e.g., 'cache.max_size')
            default: Default value if key not found

        Returns:
            Configuration value
        """
        keys = key.split('.')
        value = self.config

        try:
            for k in keys:
                value = value[k]
            return value
        except (KeyError, TypeError):
            return default

    def set(self, key: str, value: Any) -> None:
        """
        Set configuration value.

        Args:
            key: Dot-separated configuration key
            value: Value to set
        """
        keys = key.split('.')
        config = self.config

        # Navigate to the parent dictionary
        for k in keys[:-1]:
            if k not in config:
                config[k] = {}
            config = config[k]

        # Set the value
        config[keys[-1]] = value

    def reset_to_defaults(self) -> None:
        """Reset configuration to defaults."""
        self.config = self.DEFAULT_CONFIG.copy()

    def validate_config(self) -> List[str]:
        """
        Validate current configuration.

        Returns:
            List of validation error messages
        """
        errors = []

        # Validate numerical ranges
        if not (10 <= self.get('precision.decimal_places') <= 100):
            errors.append("precision.decimal_places must be between 10 and 100")

        if not (0 < self.get('precision.float_tolerance') < 1e-10):
            errors.append("precision.float_tolerance must be between 0 and 1e-10")

        if not (100 <= self.get('cache.max_size') <= 100000):
            errors.append("cache.max_size must be between 100 and 100000")

        if not (0 <= self.get('cache.ttl_seconds') <= 86400):  # Max 24 hours
            errors.append("cache.ttl_seconds must be between 0 and 86400")

        # Validate plot dimensions
        if not (40 <= self.get('visualization.plot_width') <= 200):
            errors.append("visualization.plot_width must be between 40 and 200")

        if not (10 <= self.get('visualization.plot_height') <= 50):
            errors.append("visualization.plot_height must be between 10 and 50")

        return errors

    def apply_environment_variables(self) -> None:
        """
        Apply configuration overrides from environment variables.

        Environment variables should be prefixed with LOTUS_ and use underscores
        instead of dots (e.g., LOTUS_CACHE_MAX_SIZE=10000)
        """
        prefix = 'LOTUS_'

        for env_var, value in os.environ.items():
            if env_var.startswith(prefix):
                # Convert to config key
                config_key = env_var[len(prefix):].lower().replace('_', '.')

                # Try to convert value to appropriate type
                if value.lower() in ('true', 'false'):
                    value = value.lower() == 'true'
                elif value.isdigit():
                    value = int(value)
                elif self._is_float(value):
                    value = float(value)

                self.set(config_key, value)

    def _is_float(self, s: str) -> bool:
        """Check if string represents a float."""
        try:
            float(s)
            return True
        except ValueError:
            return False

    def __str__(self) -> str:
        """String representation of configuration."""
        return json.dumps(self.config, indent=2)

    def get_summary(self) -> str:
        """
        Get configuration summary.

        Returns:
            Formatted configuration summary
        """
        summary = "LOTUS Configuration Summary\n"
        summary += "=" * 40 + "\n\n"

        sections = [
            ('Precision', 'precision'),
            ('Cache', 'cache'),
            ('Spectral', 'spectral'),
            ('Testing', 'testing'),
            ('Visualization', 'visualization'),
            ('Documentation', 'documentation'),
            ('Framework', 'framework')
        ]

        for section_name, section_key in sections:
            summary += f"{section_name}:\n"
            section = self.get(section_key, {})
            for key, value in section.items():
                summary += f"  {key}: {value}\n"
            summary += "\n"

        return summary

# Global configuration instance
config = LotusConfig()

# Apply environment variable overrides
config.apply_environment_variables()

# Convenience functions
def get_config(key: str, default: Any = None) -> Any:
    """Get configuration value."""
    return config.get(key, default)

def set_config(key: str, value: Any) -> None:
    """Set configuration value."""
    config.set(key, value)

def save_config() -> None:
    """Save configuration to file."""
    config.save_config()

def validate_config() -> List[str]:
    """Validate configuration."""
    return config.validate_config()

def config_summary() -> str:
    """Get configuration summary."""
    return config.get_summary()