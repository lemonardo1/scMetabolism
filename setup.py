from setuptools import setup, find_packages

with open("README.md", "r", encoding="utf-8") as fh:
    long_description = fh.read()

setup(
    name="scmetabolism",
    version="0.2.1",
    author="Qiang Gao, Yingcheng Wu",
    author_email="gaoqiang@fudan.edu.cn",
    description="Quantifying the single-cell metabolism activity",
    long_description=long_description,
    long_description_content_type="text/markdown",
    url="https://github.com/wu-yc/scMetabolism",
    packages=find_packages(),
    classifiers=[
        "Development Status :: 4 - Beta",
        "Intended Audience :: Science/Research",
        "License :: OSI Approved :: GNU General Public License v3 (GPLv3)",
        "Operating System :: OS Independent",
        "Programming Language :: Python :: 3",
        "Programming Language :: Python :: 3.7",
        "Programming Language :: Python :: 3.8",
        "Programming Language :: Python :: 3.9",
        "Programming Language :: Python :: 3.10",
        "Topic :: Scientific/Engineering :: Bio-Informatics",
    ],
    python_requires=">=3.7",
    install_requires=[
        "numpy>=1.19.0",
        "pandas>=1.2.0",
        "scipy>=1.6.0",
        "scikit-learn>=0.24.0",
        "matplotlib>=3.3.0",
        "seaborn>=0.11.0",
        "scanpy>=1.8.0",
        "anndata>=0.7.0",
        "gseapy>=0.10.0",
        "plotly>=5.0.0",
    ],
    extras_require={
        "dev": [
            "pytest>=6.0",
            "pytest-cov>=2.0",
            "black>=21.0",
            "flake8>=3.8",
        ],
    },
    include_package_data=True,
    package_data={
        "scmetabolism": ["data/*.gmt"],
    },
)