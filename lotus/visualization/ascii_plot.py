"""
lotus.visualization.ascii_plot — ASCII-Based Plotting System
===========================================================

Pure Python ASCII plotting for spectral geometry visualizations.
No external dependencies beyond Python stdlib.
"""

import math
from typing import List, Tuple, Dict, Any, Optional
from collections import defaultdict

class ASCIIPlot:
    """
    ASCII-based plotting system for spectral data visualization.
    """

    def __init__(self, width: int = 80, height: int = 24):
        """
        Initialize the ASCII plot.

        Args:
            width: Plot width in characters
            height: Plot height in characters
        """
        self.width = width
        self.height = height
        self.data_points = []
        self.x_range = None
        self.y_range = None

    def add_data(self, x: float, y: float, label: str = None) -> None:
        """
        Add a data point to the plot.

        Args:
            x: X coordinate
            y: Y coordinate
            label: Optional point label
        """
        self.data_points.append({'x': x, 'y': y, 'label': label})

    def add_series(self, x_data: List[float], y_data: List[float],
                  label: str = None) -> None:
        """
        Add a series of data points.

        Args:
            x_data: List of x coordinates
            y_data: List of y coordinates
            label: Series label
        """
        for x, y in zip(x_data, y_data):
            self.add_data(x, y, label)

    def set_ranges(self, x_min: float = None, x_max: float = None,
                  y_min: float = None, y_max: float = None) -> None:
        """
        Set the plot ranges manually.

        Args:
            x_min, x_max: X axis range
            y_min, y_max: Y axis range
        """
        self.x_range = (x_min, x_max) if x_min is not None and x_max is not None else None
        self.y_range = (y_min, y_max) if y_min is not None and y_max is not None else None

    def _calculate_ranges(self) -> Tuple[Tuple[float, float], Tuple[float, float]]:
        """Calculate automatic ranges from data."""
        if not self.data_points:
            return (0, 1), (0, 1)

        x_values = [p['x'] for p in self.data_points]
        y_values = [p['y'] for p in self.data_points]

        x_min, x_max = min(x_values), max(x_values)
        y_min, y_max = min(y_values), max(y_values)

        # Add some padding
        x_padding = (x_max - x_min) * 0.1 if x_max != x_min else 1
        y_padding = (y_max - y_min) * 0.1 if y_max != y_min else 1

        x_range = (x_min - x_padding, x_max + x_padding)
        y_range = (y_min - y_padding, y_max + y_padding)

        return x_range, y_range

    def _scale_point(self, x: float, y: float,
                    x_range: Tuple[float, float],
                    y_range: Tuple[float, float]) -> Tuple[int, int]:
        """
        Scale a data point to plot coordinates.

        Args:
            x, y: Data coordinates
            x_range, y_range: Data ranges

        Returns:
            Tuple of (plot_x, plot_y)
        """
        x_min, x_max = x_range
        y_min, y_max = y_range

        # Scale to 0-1 range
        if x_max != x_min:
            x_scaled = (x - x_min) / (x_max - x_min)
        else:
            x_scaled = 0.5

        if y_max != y_min:
            y_scaled = (y - y_min) / (y_max - y_min)
        else:
            y_scaled = 0.5

        # Scale to plot dimensions (leave room for axes)
        plot_x = int(x_scaled * (self.width - 8)) + 4
        plot_y = int((1 - y_scaled) * (self.height - 4)) + 2

        return plot_x, plot_y

    def plot(self, title: str = "ASCII Plot") -> str:
        """
        Generate the ASCII plot.

        Args:
            title: Plot title

        Returns:
            ASCII plot as string
        """
        if not self.data_points:
            return f"{title}\n(No data)"

        # Calculate ranges
        x_range = self.x_range or self._calculate_ranges()[0]
        y_range = self.y_range or self._calculate_ranges()[1]

        # Initialize plot grid
        grid = [[' ' for _ in range(self.width)] for _ in range(self.height)]

        # Draw axes
        self._draw_axes(grid, x_range, y_range)

        # Plot data points
        for point in self.data_points:
            plot_x, plot_y = self._scale_point(point['x'], point['y'], x_range, y_range)
            if 0 <= plot_x < self.width and 0 <= plot_y < self.height:
                grid[plot_y][plot_x] = '●'

        # Add labels
        self._add_labels(grid, x_range, y_range)

        # Convert to string
        plot_str = f"{title}\n"
        plot_str += "=" * self.width + "\n"

        for row in grid:
            plot_str += ''.join(row) + "\n"

        return plot_str

    def _draw_axes(self, grid: List[List[str]],
                  x_range: Tuple[float, float],
                  y_range: Tuple[float, float]) -> None:
        """Draw the plot axes."""
        # Y-axis
        for y in range(2, self.height - 2):
            grid[y][4] = '│'

        # X-axis
        for x in range(4, self.width - 4):
            grid[self.height - 3][x] = '─'

        # Axis intersection
        grid[self.height - 3][4] = '┼'

        # Add tick marks and labels
        x_min, x_max = x_range
        y_min, y_max = y_range

        # X-axis ticks
        for i in range(5):
            x_pos = 4 + i * (self.width - 8) // 4
            grid[self.height - 3][x_pos] = '┬'
            # Add tick label
            tick_value = x_min + i * (x_max - x_min) / 4
            label = f"{tick_value:.1f}"
            for j, char in enumerate(label):
                if x_pos + j < self.width:
                    grid[self.height - 2][x_pos + j] = char

        # Y-axis ticks
        for i in range(5):
            y_pos = 2 + i * (self.height - 4) // 4
            grid[y_pos][4] = '├'
            # Add tick label
            tick_value = y_max - i * (y_max - y_min) / 4
            label = f"{tick_value:.1f}"
            label_start = 0
            for j, char in enumerate(label):
                if label_start + j < 4:
                    grid[y_pos][label_start + j] = char

    def _add_labels(self, grid: List[List[str]],
                   x_range: Tuple[float, float],
                   y_range: Tuple[float, float]) -> None:
        """Add axis labels."""
        # X-axis label
        x_label = "X"
        x_label_start = self.width // 2 - len(x_label) // 2
        for i, char in enumerate(x_label):
            if x_label_start + i < self.width:
                grid[self.height - 1][x_label_start + i] = char

        # Y-axis label
        y_label = "Y"
        for i, char in enumerate(y_label):
            if i < self.height - 4:
                grid[i + 2][0] = char

class SpectralPlotter:
    """
    Specialized plotter for spectral geometry data.
    """

    def __init__(self):
        self.plotter = ASCIIPlot()

    def plot_eigenvalues(self, eigenvalues: List[float],
                        title: str = "Eigenvalue Spectrum") -> str:
        """
        Plot eigenvalue spectrum.

        Args:
            eigenvalues: List of eigenvalues
            title: Plot title

        Returns:
            ASCII plot string
        """
        self.plotter = ASCIIPlot()
        for i, eigenval in enumerate(eigenvalues):
            self.plotter.add_data(i, eigenval)

        return self.plotter.plot(title)

    def plot_mass_spectrum(self, masses: List[float], labels: List[str] = None,
                          title: str = "Mass Spectrum") -> str:
        """
        Plot hadron mass spectrum.

        Args:
            masses: List of masses
            labels: Optional mass labels
            title: Plot title

        Returns:
            ASCII plot string
        """
        self.plotter = ASCIIPlot()
        for i, mass in enumerate(masses):
            label = labels[i] if labels and i < len(labels) else None
            self.plotter.add_data(i, mass, label)

        return self.plotter.plot(title)

    def plot_prediction_vs_measured(self, predictions: List[float],
                                  measured: List[float],
                                  labels: List[str] = None,
                                  title: str = "Prediction vs Measured") -> str:
        """
        Plot predictions vs measured values.

        Args:
            predictions: List of predicted values
            measured: List of measured values
            labels: Optional point labels
            title: Plot title

        Returns:
            ASCII plot string
        """
        self.plotter = ASCIIPlot()

        # Add prediction points
        for i, (pred, meas) in enumerate(zip(predictions, measured)):
            label = labels[i] if labels and i < len(labels) else f"P{i+1}"
            self.plotter.add_data(pred, meas, label)

        # Add diagonal reference line
        min_val = min(min(predictions), min(measured))
        max_val = max(max(predictions), max(measured))
        self.plotter.add_data(min_val, min_val, "Reference")
        self.plotter.add_data(max_val, max_val, "Reference")

        return self.plotter.plot(title)

    def plot_error_distribution(self, errors: List[float],
                              title: str = "Error Distribution") -> str:
        """
        Plot distribution of prediction errors.

        Args:
            errors: List of relative errors
            title: Plot title

        Returns:
            ASCII plot string
        """
        self.plotter = ASCIIPlot()

        # Create histogram data
        if errors:
            min_err = min(errors)
            max_err = max(errors)
            bins = 10
            bin_width = (max_err - min_err) / bins if max_err != min_err else 1

            histogram = [0] * bins
            for error in errors:
                bin_idx = min(int((error - min_err) / bin_width), bins - 1)
                histogram[bin_idx] += 1

            for i, count in enumerate(histogram):
                bin_center = min_err + (i + 0.5) * bin_width
                self.plotter.add_data(bin_center, count)

        return self.plotter.plot(title)

# Convenience functions
def plot_eigenvalues(eigenvalues: List[float], title: str = "Eigenvalue Spectrum") -> str:
    """Plot eigenvalue spectrum."""
    plotter = SpectralPlotter()
    return plotter.plot_eigenvalues(eigenvalues, title)

def plot_mass_spectrum(masses: List[float], labels: List[str] = None,
                      title: str = "Mass Spectrum") -> str:
    """Plot hadron mass spectrum."""
    plotter = SpectralPlotter()
    return plotter.plot_mass_spectrum(masses, labels, title)

def plot_prediction_accuracy(predictions: List[float], measured: List[float],
                           labels: List[str] = None,
                           title: str = "Prediction vs Measured") -> str:
    """Plot prediction accuracy."""
    plotter = SpectralPlotter()
    return plotter.plot_prediction_vs_measured(predictions, measured, labels, title)

def plot_errors(errors: List[float], title: str = "Error Distribution") -> str:
    """Plot error distribution."""
    plotter = SpectralPlotter()
    return plotter.plot_error_distribution(errors, title)