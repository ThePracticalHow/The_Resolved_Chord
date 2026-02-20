"""
lotus.testing.framework — Self-Contained Testing Framework
==========================================================

Pure Python unit testing and validation framework for LOTUS.
No external dependencies beyond Python stdlib (unittest).
"""

import unittest
import time
import sys
from typing import List, Dict, Any, Callable, Optional
from collections import defaultdict
import traceback

class SpectralTestCase(unittest.TestCase):
    """
    Base test case for spectral geometry tests.

    Provides utilities for numerical comparisons and spectral validation.
    """

    def assertAlmostEqualSpectral(self, first: float, second: float,
                                 places: int = 7, msg: str = None) -> None:
        """
        Assert that two spectral values are almost equal.

        Uses relative error for spectral comparisons.

        Args:
            first: First value
            second: Second value
            places: Decimal places for comparison
            msg: Optional error message
        """
        if abs(second) < 1e-15:
            # Absolute comparison for very small values
            self.assertAlmostEqual(first, second, places=places, msg=msg)
        else:
            # Relative comparison
            relative_error = abs((first - second) / second)
            max_error = 10 ** (-places)
            if relative_error > max_error:
                standard_msg = f"Spectral values not almost equal: {first} vs {second} (relative error: {relative_error:.2e})"
                self.fail(self._formatMessage(msg, standard_msg))

    def assertPredictionAccuracy(self, predicted: float, measured: float,
                                tolerance_percent: float, msg: str = None) -> None:
        """
        Assert that a prediction is within tolerance of measured value.

        Args:
            predicted: Predicted value
            measured: Measured value
            tolerance_percent: Tolerance as percentage
            msg: Optional error message
        """
        if abs(measured) < 1e-15:
            self.fail(f"Measured value too small for percentage comparison: {measured}")

        error_percent = abs((predicted - measured) / measured) * 100
        if error_percent > tolerance_percent:
            standard_msg = f"Prediction accuracy {error_percent:.3f}% exceeds tolerance {tolerance_percent}%: predicted {predicted}, measured {measured}"
            self.fail(self._formatMessage(msg, standard_msg))

class SpectralTestSuite:
    """
    Test suite for running spectral geometry tests.

    Provides test discovery, execution, and reporting.
    """

    def __init__(self):
        self.tests = []
        self.results = defaultdict(list)
        self.start_time = None
        self.end_time = None

    def add_test(self, test_case: unittest.TestCase) -> None:
        """
        Add a test case to the suite.

        Args:
            test_case: Test case instance
        """
        self.tests.append(test_case)

    def discover_tests(self, module) -> None:
        """
        Discover and add all test cases from a module.

        Args:
            module: Module to search for test cases
        """
        for name in dir(module):
            obj = getattr(module, name)
            if (isinstance(obj, type) and
                issubclass(obj, unittest.TestCase) and
                obj != unittest.TestCase):
                # Create instance and add to suite
                try:
                    test_instance = obj()
                    self.add_test(test_instance)
                except Exception as e:
                    print(f"Warning: Could not instantiate test {name}: {e}")

    def run_tests(self, verbosity: int = 1) -> Dict[str, Any]:
        """
        Run all tests in the suite.

        Args:
            verbosity: Verbosity level (0-2)

        Returns:
            Dictionary with test results
        """
        self.start_time = time.time()

        loader = unittest.TestLoader()
        suite = unittest.TestSuite()

        for test_case in self.tests:
            suite.addTest(loader.loadTestsFromTestCase(type(test_case)))

        runner = unittest.TextTestRunner(
            verbosity=verbosity,
            stream=sys.stdout,
            resultclass=SpectralTestResult
        )

        result = runner.run(suite)
        self.end_time = time.time()

        return {
            'tests_run': result.testsRun,
            'failures': len(result.failures),
            'errors': len(result.errors),
            'skipped': len(result.skipped),
            'expected_failures': len(result.expectedFailures),
            'unexpected_successes': len(result.unexpectedSuccesses),
            'execution_time': self.end_time - self.start_time,
            'success': result.wasSuccessful()
        }

class SpectralTestResult(unittest.TextTestResult):
    """
    Custom test result class with spectral-specific reporting.
    """

    def __init__(self, stream, descriptions, verbosity):
        super().__init__(stream, descriptions, verbosity)
        self.spectral_errors = []

    def addFailure(self, test, err):
        super().addFailure(test, err)
        self._record_spectral_error(test, err, 'FAILURE')

    def addError(self, test, err):
        super().addError(test, err)
        self._record_spectral_error(test, err, 'ERROR')

    def _record_spectral_error(self, test, err, error_type: str) -> None:
        """Record spectral-specific error information."""
        error_info = {
            'test': str(test),
            'type': error_type,
            'exception': str(err[0]),
            'traceback': ''.join(traceback.format_exception(*err))
        }
        self.spectral_errors.append(error_info)

class PredictionValidator:
    """
    Validator for checking prediction accuracy against experimental data.
    """

    def __init__(self):
        self.predictions = []
        self.rms_errors = []

    def add_prediction(self, name: str, predicted: float, measured: float,
                      uncertainty: Optional[float] = None) -> None:
        """
        Add a prediction for validation.

        Args:
            name: Name of the prediction
            predicted: Predicted value
            measured: Measured value
            uncertainty: Measurement uncertainty (optional)
        """
        if abs(measured) < 1e-15:
            print(f"Warning: Measured value for {name} is too small: {measured}")
            return

        error = abs((predicted - measured) / measured)
        self.predictions.append({
            'name': name,
            'predicted': predicted,
            'measured': measured,
            'uncertainty': uncertainty,
            'relative_error': error
        })
        self.rms_errors.append(error)

    def calculate_rms_error(self) -> float:
        """
        Calculate RMS error across all predictions.

        Returns:
            RMS error as percentage
        """
        if not self.rms_errors:
            return 0.0

        squared_errors = [err ** 2 for err in self.rms_errors]
        mean_squared_error = sum(squared_errors) / len(squared_errors)
        rms_error = (mean_squared_error ** 0.5) * 100  # Convert to percentage

        return rms_error

    def get_worst_predictions(self, n: int = 5) -> List[Dict[str, Any]]:
        """
        Get the predictions with the worst accuracy.

        Args:
            n: Number of predictions to return

        Returns:
            List of worst predictions
        """
        sorted_predictions = sorted(
            self.predictions,
            key=lambda x: x['relative_error'],
            reverse=True
        )
        return sorted_predictions[:n]

    def get_best_predictions(self, n: int = 5) -> List[Dict[str, Any]]:
        """
        Get the predictions with the best accuracy.

        Args:
            n: Number of predictions to return

        Returns:
            List of best predictions
        """
        sorted_predictions = sorted(
            self.predictions,
            key=lambda x: x['relative_error']
        )
        return sorted_predictions[:n]

    def report(self) -> str:
        """
        Generate a validation report.

        Returns:
            Formatted report string
        """
        if not self.predictions:
            return "No predictions to validate."

        rms = self.calculate_rms_error()
        worst = self.get_worst_predictions(3)
        best = self.get_best_predictions(3)

        report = f"""
LOTUS Prediction Validation Report
==================================

Total Predictions: {len(self.predictions)}
RMS Error: {rms:.3f}%

Worst Predictions:
"""

        for pred in worst:
            report += f"  {pred['name']}: {pred['relative_error']*100:.3f}% error\n"
            report += f"    Predicted: {pred['predicted']}, Measured: {pred['measured']}\n"

        report += "\nBest Predictions:\n"
        for pred in best:
            report += f"  {pred['name']}: {pred['relative_error']*100:.6f}% error\n"
            report += f"    Predicted: {pred['predicted']}, Measured: {pred['measured']}\n"

        return report

# Global instances
test_suite = SpectralTestSuite()
validator = PredictionValidator()

# Convenience functions
def run_spectral_tests(verbosity: int = 1) -> Dict[str, Any]:
    """Run all spectral tests."""
    return test_suite.run_tests(verbosity=verbosity)

def validate_predictions(predictions: List[Dict[str, Any]]) -> str:
    """
    Validate a list of predictions.

    Args:
        predictions: List of prediction dictionaries with keys:
                    'name', 'predicted', 'measured', 'uncertainty' (optional)

    Returns:
        Validation report
    """
    for pred in predictions:
        validator.add_prediction(
            pred['name'],
            pred['predicted'],
            pred['measured'],
            pred.get('uncertainty')
        )

    return validator.report()

def get_rms_error() -> float:
    """Get current RMS error across all predictions."""
    return validator.calculate_rms_error()