# LOTUS Development Makefile
# Common tasks for development, testing, and deployment

.PHONY: help install test verify falsify clean gate autoclean docker-build docker-run docker-test pypi-build pypi-upload

# Default target
help:
	@echo "LOTUS Development Commands:"
	@echo ""
	@echo "Installation:"
	@echo "  make install          Install LOTUS in development mode"
	@echo "  make install-prod     Install LOTUS for production"
	@echo ""
	@echo "Testing:"
	@echo "  make test             Run all tests"
	@echo "  make verify           Run verification scripts"
	@echo "  make falsify          Run falsification test suite"
	@echo "  make compile          Compile the universe (quick test)"
	@echo ""
	@echo "Docker:"
	@echo "  make docker-build     Build Docker image"
	@echo "  make docker-run       Run LOTUS in Docker container"
	@echo "  make docker-test      Test Docker container functionality"
	@echo ""
	@echo "Publishing:"
	@echo "  make pypi-build       Build PyPI distribution"
	@echo "  make pypi-upload      Upload to PyPI (requires credentials)"
	@echo ""
	@echo "Maintenance:"
	@echo "  make clean            Remove build artifacts and cache files"
	@echo "  make gate             Run robustness release gate checks"
	@echo "  make autoclean        Auto-clean caches/build artifacts (apply mode)"
	@echo "  make lint             Check code style"
	@echo "  make format           Format code with black"

# Installation
install:
	pip install -e ".[dev]"

install-prod:
	pip install .

# Testing
test:
	python -m pytest tools/falsification/ -v

verify:
	python verification/run_verification.py

falsify:
	python -m pytest tools/falsification/ -v

compile:
	python verification/compile_universe.py

# Docker
docker-build:
	docker build -t lotus-bloom .

docker-run:
	docker run --rm lotus-bloom

docker-test:
	docker run --rm -v $(PWD):/app lotus-bloom python tools/tests/test_docker.py

# PyPI
pypi-build:
	python -m build

pypi-upload: pypi-build
	twine upload dist/*

# Maintenance
clean:
	python tools/auto_clean.py --apply

gate:
	python tools/release_gate.py

autoclean:
	python tools/auto_clean.py --apply

lint:
	flake8 lotus/ verification/ tools/falsification/ --max-line-length=100 --ignore=E203,W503
	mypy lotus/ --ignore-missing-imports

format:
	black lotus/ verification/ tools/falsification/ --line-length=100
	isort lotus/ verification/ tools/falsification/