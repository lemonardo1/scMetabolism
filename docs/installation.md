# Installation Guide

## Quick Installation

### From PyPI (Recommended)

```bash
pip install scmetabolism
```

### From GitHub

```bash
# Latest stable release
pip install git+https://github.com/your-username/scMetabolism-python.git

# Development version
pip install git+https://github.com/your-username/scMetabolism-python.git@develop
```

## Installation Options

### Basic Installation

For basic functionality:

```bash
pip install scmetabolism
```

### With Scanpy Support

For single-cell analysis with Scanpy integration:

```bash
pip install "scmetabolism[scanpy]"
```

### Full Installation

For all optional features:

```bash
pip install "scmetabolism[all]"
```

### Development Installation

For contributing to the project:

```bash
git clone https://github.com/your-username/scMetabolism-python.git
cd scMetabolism-python
pip install -e .[dev]
```

## System Requirements

### Python Version
- Python 3.7 or higher
- Tested on Python 3.7, 3.8, 3.9, 3.10, 3.11

### Operating Systems
- Linux (Ubuntu 18.04+, CentOS 7+)
- macOS (10.14+)
- Windows (10+)

### Memory Requirements
- Minimum: 4GB RAM
- Recommended: 8GB+ RAM for large datasets
- Large datasets (>50k cells): 16GB+ RAM recommended

## Dependencies

### Required Dependencies

These are automatically installed:

- numpy >= 1.19.0
- pandas >= 1.2.0
- scipy >= 1.6.0
- scikit-learn >= 0.24.0
- matplotlib >= 3.3.0
- seaborn >= 0.11.0

### Optional Dependencies

Install as needed:

```bash
# For Scanpy integration
pip install scanpy>=1.8.0 anndata>=0.7.0

# For enhanced gene set analysis
pip install gseapy>=0.10.0

# For interactive plots
pip install plotly>=5.0.0
```

## Verification

Test your installation:

```python
import scmetabolism
print(scmetabolism.__version__)

# Run a quick test
from scmetabolism import ScMetabolism
sc_metab = ScMetabolism()
print("Installation successful!")
```

## Troubleshooting

### Common Issues

#### 1. ImportError: No module named 'scmetabolism'

**Solution**: Make sure you installed the package correctly:
```bash
pip list | grep scmetabolism
```

If not found, reinstall:
```bash
pip install --upgrade scmetabolism
```

#### 2. Version conflicts

**Solution**: Create a fresh virtual environment:
```bash
python -m venv scmetab_env
source scmetab_env/bin/activate  # Linux/Mac
# or
scmetab_env\Scripts\activate  # Windows

pip install scmetabolism
```

#### 3. Memory errors with large datasets

**Solutions**:
- Increase system memory
- Use data subsampling
- Process data in chunks
- Use `n_cores=1` to reduce memory usage

#### 4. Slow performance

**Solutions**:
- Increase `n_cores` parameter
- Use SSD storage for data
- Consider data preprocessing to reduce size

### Platform-Specific Issues

#### macOS

If you encounter compilation errors:
```bash
# Install Xcode command line tools
xcode-select --install

# Update pip and setuptools
pip install --upgrade pip setuptools
```

#### Windows

For Windows users, consider using Anaconda:
```bash
conda install -c conda-forge scmetabolism
```

#### Linux

On older Linux distributions, you might need:
```bash
# Update system packages
sudo apt-get update
sudo apt-get install python3-dev build-essential

# Or for CentOS/RHEL
sudo yum install python3-devel gcc gcc-c++
```

## Virtual Environments

### Using venv

```bash
# Create virtual environment
python -m venv scmetab_env

# Activate
source scmetab_env/bin/activate  # Linux/Mac
scmetab_env\Scripts\activate     # Windows

# Install
pip install scmetabolism

# Deactivate when done
deactivate
```

### Using conda

```bash
# Create environment
conda create -n scmetab python=3.9

# Activate
conda activate scmetab

# Install
pip install scmetabolism

# Deactivate when done
conda deactivate
```

## Docker Installation

For reproducible environments:

```dockerfile
FROM python:3.9-slim

RUN pip install scmetabolism[all]

# Your analysis code here
COPY . /app
WORKDIR /app

CMD ["python", "your_analysis.py"]
```

Build and run:
```bash
docker build -t scmetab-analysis .
docker run -v $(pwd):/data scmetab-analysis
```

## Jupyter Notebook Setup

For interactive analysis:

```bash
# Install Jupyter
pip install jupyter

# Install scmetabolism
pip install scmetabolism[all]

# Start Jupyter
jupyter notebook
```

In your notebook:
```python
import scmetabolism
import pandas as pd
import matplotlib.pyplot as plt

# Your analysis here
```

## Performance Optimization

### For Large Datasets

1. **Use appropriate data types**:
```python
# Use sparse matrices when possible
import scipy.sparse as sp
count_matrix = sp.csr_matrix(count_matrix)
```

2. **Optimize memory usage**:
```python
# Process in chunks
chunk_size = 1000
for i in range(0, n_cells, chunk_size):
    chunk = count_matrix[:, i:i+chunk_size]
    # Process chunk
```

3. **Use multiple cores**:
```python
sc_metab.compute_metabolism(
    count_matrix,
    n_cores=4  # Adjust based on your system
)
```

## Getting Help

If you encounter issues:

1. Check the [FAQ](https://github.com/your-username/scMetabolism-python/wiki/FAQ)
2. Search [existing issues](https://github.com/your-username/scMetabolism-python/issues)
3. Create a [new issue](https://github.com/your-username/scMetabolism-python/issues/new) with:
   - Your system information
   - Complete error message
   - Minimal reproducible example