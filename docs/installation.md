
# Installation

Simple Docker-based setup for CrossBuild Assessor.

## Requirements

- Docker
- 8GB+ RAM for large datasets

## Setup

```bash
git clone https://github.com/your-org/crossbuild-assessor.git
cd crossbuild-assessor
docker build -t crossbuild-assessor .
```

## Test installation

```bash
docker run crossbuild-assessor python db_loader.py --help
```

You should see the help message.

## Next steps

1. Prepare your input data (see [user guide](user-guide.md))
2. Create a configuration file
3. Run the pipeline
