"""
lotus.core.cache — Pure Python Caching System
==============================================

Self-contained LRU cache implementation for spectral computations.
No external dependencies beyond Python stdlib.
"""

import time
from typing import Any, Dict, Optional, Tuple
from collections import OrderedDict
import threading

class LRUCache:
    """
    Least Recently Used (LRU) cache implementation.

    Thread-safe and memory-efficient caching for spectral computations.
    """

    def __init__(self, capacity: int = 1000, ttl: Optional[float] = None):
        """
        Initialize the LRU cache.

        Args:
            capacity: Maximum number of items to store
            ttl: Time-to-live in seconds (None for no expiration)
        """
        self.capacity = capacity
        self.ttl = ttl
        self.cache: OrderedDict = OrderedDict()
        self.timestamps: Dict[Any, float] = {}
        self.lock = threading.RLock()

    def get(self, key: Any) -> Optional[Any]:
        """
        Retrieve an item from the cache.

        Args:
            key: Cache key

        Returns:
            Cached value or None if not found/expired
        """
        with self.lock:
            if key not in self.cache:
                return None

            # Check TTL
            if self.ttl is not None:
                if time.time() - self.timestamps[key] > self.ttl:
                    del self.cache[key]
                    del self.timestamps[key]
                    return None

            # Move to end (most recently used)
            self.cache.move_to_end(key)
            return self.cache[key]

    def put(self, key: Any, value: Any) -> None:
        """
        Store an item in the cache.

        Args:
            key: Cache key
            value: Value to store
        """
        with self.lock:
            if key in self.cache:
                # Update existing
                self.cache.move_to_end(key)
            else:
                # Add new
                if len(self.cache) >= self.capacity:
                    # Remove least recently used
                    oldest_key, _ = self.cache.popitem(last=False)
                    if oldest_key in self.timestamps:
                        del self.timestamps[oldest_key]

            self.cache[key] = value
            self.timestamps[key] = time.time()

    def clear(self) -> None:
        """Clear all items from the cache."""
        with self.lock:
            self.cache.clear()
            self.timestamps.clear()

    def size(self) -> int:
        """Return the current number of items in the cache."""
        with self.lock:
            return len(self.cache)

    def __contains__(self, key: Any) -> bool:
        """Check if key is in cache (and not expired)."""
        with self.lock:
            if key not in self.cache:
                return False

            if self.ttl is not None:
                if time.time() - self.timestamps[key] > self.ttl:
                    del self.cache[key]
                    del self.timestamps[key]
                    return False

            return True

    def __len__(self) -> int:
        """Return the number of items in the cache."""
        return self.size()

class SpectralCache:
    """
    Specialized cache for spectral geometry computations.

    Provides caching for expensive spectral calculations with
    automatic cache key generation based on parameters.
    """

    def __init__(self, capacity: int = 5000):
        self.cache = LRUCache(capacity=capacity, ttl=3600)  # 1 hour TTL

    def _make_key(self, func_name: str, *args, **kwargs) -> Tuple:
        """
        Generate a cache key from function name and arguments.

        Args:
            func_name: Name of the function
            *args: Positional arguments
            **kwargs: Keyword arguments

        Returns:
            Tuple that can be used as a cache key
        """
        # Convert args to hashable types
        hashable_args = []
        for arg in args:
            if isinstance(arg, (list, tuple)):
                hashable_args.append(tuple(arg))
            elif isinstance(arg, dict):
                hashable_args.append(tuple(sorted(arg.items())))
            else:
                hashable_args.append(arg)

        hashable_kwargs = tuple(sorted(kwargs.items()))
        return (func_name, tuple(hashable_args), hashable_kwargs)

    def cached(self, func_name: str = None):
        """
        Decorator to cache function results.

        Args:
            func_name: Optional function name override

        Returns:
            Decorated function
        """
        def decorator(func):
            name = func_name or func.__name__

            def wrapper(*args, **kwargs):
                key = self._make_key(name, *args, **kwargs)
                result = self.cache.get(key)
                if result is not None:
                    return result

                result = func(*args, **kwargs)
                self.cache.put(key, result)
                return result

            return wrapper
        return decorator

    def get(self, func_name: str, *args, **kwargs):
        """
        Manually retrieve a cached result.

        Args:
            func_name: Function name
            *args, **kwargs: Function arguments

        Returns:
            Cached result or None
        """
        key = self._make_key(func_name, *args, **kwargs)
        return self.cache.get(key)

    def put(self, func_name: str, result: Any, *args, **kwargs) -> None:
        """
        Manually store a result in the cache.

        Args:
            func_name: Function name
            result: Result to cache
            *args, **kwargs: Function arguments
        """
        key = self._make_key(func_name, *args, **kwargs)
        self.cache.put(key, result)

    def clear(self) -> None:
        """Clear all cached results."""
        self.cache.clear()

    def stats(self) -> Dict[str, int]:
        """
        Get cache statistics.

        Returns:
            Dictionary with cache statistics
        """
        return {
            'size': self.cache.size(),
            'capacity': self.cache.capacity
        }

# Global spectral cache instance
spectral_cache = SpectralCache()

# Convenience functions
def cached(func_name: str = None):
    """Decorator for caching spectral computations."""
    return spectral_cache.cached(func_name)

def get_cached_result(func_name: str, *args, **kwargs):
    """Get a cached result."""
    return spectral_cache.get(func_name, *args, **kwargs)

def put_cached_result(func_name: str, result: Any, *args, **kwargs) -> None:
    """Store a result in the cache."""
    spectral_cache.put(func_name, result, *args, **kwargs)

def clear_cache() -> None:
    """Clear the spectral cache."""
    spectral_cache.clear()

def cache_stats() -> Dict[str, int]:
    """Get cache statistics."""
    return spectral_cache.stats()