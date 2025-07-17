FROM python:3.11-slim

# Set working directory
WORKDIR /app

# Install system dependencies
RUN apt-get update && apt-get install -y \
    gcc \
    g++ \
    sqlite3 \
    && rm -rf /var/lib/apt/lists/*

# Copy requirements and install Python dependencies
COPY requirements.txt .
RUN pip install --no-cache-dir -r requirements.txt

# Copy main scripts at root level
COPY db_loader.py .
COPY db_analyzer.py .
COPY variant_prioritizer.py .

# Copy configuration modules
COPY config/ ./config/
COPY analysis/ ./analysis/
COPY visualization/ ./visualization/
COPY utils/ ./utils/

# Copy documentation
COPY docs/ ./docs/

# Copy example configuration
COPY config.json.example ./config.json.example

# Create data and output directories
RUN mkdir -p /data /output

# Set environment variables
ENV PYTHONPATH="/app"
ENV PYTHONUNBUFFERED=1

# Create non-root user for security
RUN groupadd -r crossbuild && useradd -r -g crossbuild crossbuild
RUN chown -R crossbuild:crossbuild /app /data /output
USER crossbuild

# Default command shows help
CMD ["python", "--help"]

# Labels for metadata
LABEL maintainer="aksenia"
LABEL description="Genomic variant analysis tool for liftover quality assessment"
LABEL version="1.0"

# Usage examples in labels
LABEL usage.db_loader="docker run -v /path/to/data:/data crossbuild-assessor python db_loader.py --config /data/config.json"
LABEL usage.analyzer="docker run -v /path/to/data:/data crossbuild-assessor python db_analyzer.py --db-path /data/genomic_analysis.db --output-dir /data/results/"
LABEL usage.prioritizer="docker run -v /path/to/data:/data crossbuild-assessor python variant_prioritizer.py --db-path /data/genomic_analysis.db --output-dir /data/results/"