# Use latest Python image based on Linux
FROM python:3.11-slim

# Set working directory
WORKDIR /app

# Install system dependencies including SQLite
RUN apt-get update && apt-get install -y \
    gcc \
    sqlite3 \
    && rm -rf /var/lib/apt/lists/*

# Install Python packages
RUN pip install --no-cache-dir \
    pandas==2.1.4 \
    matplotlib==3.8.2 \
    seaborn==0.13.0 \
    numpy==1.24.3

# Create data and output directories
RUN mkdir -p data output

# Copy analysis scripts
COPY db_loader.py .
COPY db_analyzer.py .
COPY variant_prioritizer.py .

# Set environment variables
ENV PYTHONUNBUFFERED=1
ENV MPLCONFIGDIR=/tmp/matplotlib

# Default command (can be overridden)
CMD ["python", "--help"]
