# LOTUS Model - Containerized Environment
# The Spectral Geometry of Everything

FROM python:3.11-slim

# Metadata
LABEL maintainer="Jixiang Leng <jixiangleng1@gmail.com>"
LABEL description="LOTUS: Lens Orbifold Theory of the Unified Spectrum"
LABEL version="1.0.0"
LABEL org.opencontainers.image.source="https://github.com/ThePracticalHow/The_Resolved_Chord"

# Set working directory
WORKDIR /app

# Install system dependencies for scientific computing
RUN apt-get update && apt-get install -y \
    git \
    && rm -rf /var/lib/apt/lists/*

# Copy and install Python dependencies
COPY requirements.txt .
RUN pip install --no-cache-dir -r requirements.txt

# Copy the entire package
COPY . .

# Install LOTUS in development mode
RUN pip install -e .

# Make test script executable
RUN chmod +x tools/tests/test_docker.py

# Create non-root user for security
RUN useradd --create-home --shell /bin/bash lotus
USER lotus

# Default command: compile the universe
CMD ["python", "verification/compile_universe.py"]

# Expose port for potential web interface (future)
# EXPOSE 8000

# Health check
HEALTHCHECK --interval=30s --timeout=10s --start-period=5s --retries=3 \
    CMD python -c "import lotus; u = lotus.Universe(); print('LOTUS ready')" || exit 1