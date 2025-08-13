FROM python:3.11-slim

# Set working directory
WORKDIR /app

# Install system dependencies including build tools, perl, and others needed for VEP and bcftools
RUN apt-get update && apt-get install -y \
    curl \
    wget \
    tar \
    unzip \
    gzip \
    bzip2 \
    make \
    gcc \
    g++ \
    perl \
    cpanminus \
    libperl-dev \
    zlib1g-dev \
    libbz2-dev \
    liblzma-dev \
    libcurl4-openssl-dev \
    libwww-perl \
    && cpanm --notest DBI Archive::Extract Archive::Tar Archive::Zip LWP::Simple \
    && rm -rf /var/lib/apt/lists/*

# Install CrossMap python library and other python dependencies
RUN pip install --no-cache-dir CrossMap pysam snakemake

# Install htslib for bgzip and tabix
RUN curl -L https://github.com/samtools/htslib/releases/download/1.22/htslib-1.22.tar.bz2 -o htslib-1.22.tar.bz2 \
    && tar -xjf htslib-1.22.tar.bz2 \
    && cd htslib-1.22 \
    && ./configure --prefix=/usr/local \
    && make \
    && make install


# Install latest bcftools with liftover plugin
RUN wget https://github.com/samtools/bcftools/releases/download/1.22/bcftools-1.22.tar.bz2 && \
    tar xjf bcftools-1.22.tar.bz2 && \
    cd bcftools-1.22 && \
    # Remove all plugin source files to avoid building unwanted plugins
    /bin/rm -f plugins/*.{c,h,mk} && \
    # Download only liftover plugin source from freeseek/score repo
    wget -P plugins https://raw.githubusercontent.com/freeseek/score/master/liftover.c && \
    # Remove plugin makefile if exists (e.g. pgs.mk), just in case
    /bin/rm -f plugins/pgs.mk || true && \
    # Configure and build bcftools
    ./configure --prefix=/usr/local && \
    make -j$(nproc) && \
    # Build liftover plugin explicitly
    make plugins/liftover.so && \
    mkdir -p /usr/local/bin /usr/local/lib/bcftools_plugins && \
    # Copy bcftools executable
    cp bcftools /usr/local/bin/ && \
    # Copy only the liftover plugin shared library
    cp plugins/liftover.so /usr/local/lib/bcftools_plugins/ && \
    cd .. && \
    rm -rf bcftools-1.22 bcftools-1.22.tar.bz2

ENV BCFTOOLS_PLUGINS=/usr/local/lib/bcftools_plugins

# Install Ensembl VEP version 113.3 (without caches)
RUN wget https://github.com/Ensembl/ensembl-vep/archive/refs/tags/release/113.3.tar.gz -O vep.tar.gz && \
    tar -xzf vep.tar.gz && \
    cd ensembl-vep-release-113.3 && \
    perl INSTALL.pl -n -a a --NO_UPDATE --NO_HTSLIB --NO_PLUGIN --NO_TEST && \
    cd .. && \
    rm -f vep.tar.gz

# Copy CrossMap worker script
COPY liftover/ ./liftover/

# Copy the bash script for downloading VEP caches
COPY download_vep_caches.sh .

COPY snake/ ./snake/

# Make scripts executable
RUN chmod +x ./liftover/*.py
RUN chmod +x download_vep_caches.sh

# Create non-root user for security and set ownership
RUN groupadd -r crossbuild && useradd -r -g crossbuild crossbuild && \
    mkdir -p /home/crossbuild && \
    chown -R crossbuild:crossbuild /app /home/crossbuild

USER crossbuild

ENV HOME=/home/crossbuild
ENV PERL5LIB=/app/ensembl-vep-release-113.3/modules
ENV PATH="/app/ensembl-vep-release-113.3:/usr/local/bin:${PATH}"

# Default command to show help; run the download script or your annotation script explicitly
CMD ["python", "--help"]


