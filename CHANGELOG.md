# Changelog

All notable changes to this project will be documented in this file.

The format is based on [Keep a Changelog](https://keepachangelog.com/en/1.0.0/),
and this project adheres to [Semantic Versioning](https://semver.org/spec/v2.0.0.html).

## [Unreleased]

### Added
- Initial Python implementation of scMetabolism
- Support for AUCell, ssGSEA, and GSVA methods
- ALRA imputation functionality
- Comprehensive visualization suite
- Scanpy/AnnData integration
- KEGG and REACTOME gene sets
- Complete test suite
- CI/CD pipeline

## [0.2.1] - 2024-01-XX

### Added
- **Core functionality**
  - `ScMetabolism` class for metabolism analysis
  - Support for multiple scoring methods (AUCell, ssGSEA, GSVA)
  - ALRA imputation for single-cell data
  - KEGG and REACTOME metabolism gene sets (85 and 82 pathways respectively)

- **Visualization**
  - `MetabolismPlotter` class with multiple plot types
  - Dimension reduction plots (UMAP/t-SNE)
  - Dot plots for pathway comparison
  - Box plots for distribution analysis
  - Heatmaps for comprehensive overview

- **Integration**
  - Scanpy/AnnData compatibility
  - Pandas DataFrame support
  - Multi-core processing support

- **Utilities**
  - Data preprocessing functions
  - Gene set loading utilities
  - ALRA imputation implementation

- **Documentation**
  - Comprehensive README with examples
  - API documentation
  - Usage examples in multiple formats

- **Testing**
  - Unit tests for core functionality
  - Integration tests
  - Continuous integration setup

### Technical Details
- **Dependencies**: numpy, pandas, scipy, scikit-learn, matplotlib, seaborn
- **Optional dependencies**: scanpy, anndata for single-cell integration
- **Python support**: 3.7, 3.8, 3.9, 3.10, 3.11
- **License**: GPL-3.0

### Migration from R
- Complete port of original R package functionality
- Maintained API compatibility where possible
- Enhanced Python-specific features
- Improved performance with vectorized operations

## [0.1.0] - 2024-01-XX

### Added
- Initial project setup
- Basic package structure
- Core algorithm implementations

---

## Notes

### Versioning Strategy
- **Major version** (X.0.0): Breaking changes to public API
- **Minor version** (0.X.0): New features, backward compatible
- **Patch version** (0.0.X): Bug fixes, backward compatible

### Release Process
1. Update version numbers in relevant files
2. Update this CHANGELOG.md
3. Create GitHub release with tag
4. Automatic PyPI deployment via GitHub Actions

### Contributing
See [CONTRIBUTING.md](CONTRIBUTING.md) for guidelines on contributing to this project.