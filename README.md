# scMetabolism

`scMetabolism`은 단일 세포 해상도에서 대사 활동을 정량화하는 Python 패키지입니다.

![Screenshot](logo.jpg)

## 🚀 설치

### PyPI에서 설치 (권장)

```bash
pip install scmetabolism
```

### GitHub에서 개발 버전 설치

```bash
pip install git+https://github.com/lemonardo1/scMetabolism.git
```

### 로컬 개발 설치

```bash
git clone https://github.com/lemonardo1/scMetabolism.git
cd scMetabolism
pip install -e .
```

## 📋 필요 패키지

### 필수 패키지
- numpy >= 1.19.0
- pandas >= 1.2.0
- scipy >= 1.6.0
- scikit-learn >= 0.24.0
- matplotlib >= 3.3.0
- seaborn >= 0.11.0
- numba >= 0.56.0 (성능 최적화)
- networkx >= 2.6.0 (네트워크 분석)
- pyyaml >= 5.4.0 (설정 관리)

### 선택적 패키지
```bash
# Scanpy 지원 (AnnData 객체 사용)
pip install "scmetabolism[scanpy]"

# 인터랙티브 시각화 (Plotly)
pip install "scmetabolism[plotly]"

# 모든 기능 (권장)
pip install "scmetabolism[all]"

# 개발자 도구
pip install "scmetabolism[dev]"
```

### 시스템 요구사항
- Python 3.7 이상
- 메모리: 최소 4GB (대용량 데이터의 경우 8GB 이상 권장)
- 디스크 공간: 1GB (GO 데이터베이스 캐시 포함)

## 🔬 빠른 시작

### 1. 기본 사용법

```python
import pandas as pd
import numpy as np
from scmetabolism import ScMetabolism, MetabolismPlotter

# 데이터 로드 (genes x cells 형태의 count matrix)
count_matrix = pd.read_csv("your_count_matrix.csv", index_col=0)

# ScMetabolism 객체 생성
sc_metab = ScMetabolism()

# 대사 점수 계산
metabolism_scores = sc_metab.compute_metabolism(
    count_matrix=count_matrix,
    method="aucell",  # "aucell", "ssgsea", "gsva" 중 선택
    imputation=False,  # ALRA imputation 사용 여부
    n_cores=2,
    metabolism_type="KEGG"  # "KEGG", "REACTOME", "GO_metabolism" 등
)

print(f"계산된 대사 경로 수: {metabolism_scores.shape[0]}")
print(f"분석된 세포 수: {metabolism_scores.shape[1]}")
```

### 2. Scanpy/AnnData와 함께 사용

```python
import scanpy as sc
import anndata as ad

# AnnData 객체로 작업
adata = sc.read_h5ad("your_data.h5ad")

# 대사 점수 계산 및 AnnData에 추가
adata = sc_metab.compute_metabolism_scanpy(
    adata=adata,
    method="aucell",
    metabolism_type="KEGG"
)

# 결과는 adata.obsm['metabolism']에 저장됩니다
print("대사 점수가 adata.obsm['metabolism']에 저장되었습니다")
```

### 3. 시각화

```python
# 시각화 객체 생성
plotter = MetabolismPlotter()

# 차원 축소 플롯
fig = plotter.dim_plot(
    embedding=umap_coords,  # 2D 좌표 (cells x 2)
    metabolism_scores=metabolism_scores,
    pathway="Glycolysis / Gluconeogenesis",
    embedding_type="umap"
)

# 도트 플롯
pathways_of_interest = [
    "Glycolysis / Gluconeogenesis",
    "Oxidative phosphorylation", 
    "Citrate cycle (TCA cycle)"
]

fig = plotter.dot_plot(
    metabolism_scores=metabolism_scores,
    metadata=cell_metadata,
    pathways=pathways_of_interest,
    group_by="cell_type"
)

# 박스 플롯
fig = plotter.box_plot(
    metabolism_scores=metabolism_scores,
    metadata=cell_metadata,
    pathways=pathways_of_interest,
    group_by="cell_type"
)
```

## 📊 지원하는 방법론

### 점수 계산 방법
- **AUCell**: Area Under the Curve 기반 방법 (권장)
- **ssGSEA**: Single-sample Gene Set Enrichment Analysis  
- **GSVA**: Gene Set Variation Analysis

### 유전자 세트
- **KEGG**: 85개 대사 경로
- **REACTOME**: 82개 대사 경로
- **Gene Ontology (GO)**: 
  - GO_metabolism: 대사 관련 GO 용어
  - GO_BP: 생물학적 과정 (Biological Process)
  - GO_MF: 분자 기능 (Molecular Function)
  - GO_CC: 세포 구성요소 (Cellular Component)
  - GO_all: 모든 GO 용어

## 🔧 고급 사용법

### 데이터 품질 관리

```python
from scmetabolism import DataValidator, QualityMetrics

# 데이터 검증
validator = DataValidator()
validation_results = validator.validate_count_matrix(
    count_matrix,
    min_genes_per_cell=200,
    min_cells_per_gene=3,
    max_mito_percent=20.0
)

# 검증 결과 시각화
validator.plot_validation_summary()

# 품질 지표 계산
quality_metrics = QualityMetrics.calculate_score_quality(metabolism_scores)
```

### 성능 최적화

```python
from scmetabolism import SparseMatrixHandler, MemoryOptimizer, BatchProcessor

# 희소 행렬 처리
sparse_handler = SparseMatrixHandler()
sparse_matrix = sparse_handler.to_sparse(count_matrix, threshold=0.7)

# 메모리 최적화
optimizer = MemoryOptimizer()
optimized_df = optimizer.optimize_dtypes(count_matrix)

# 배치 처리 (대용량 데이터)
batch_processor = BatchProcessor(batch_size=5000)
results = batch_processor.process_cells_in_batches(count_matrix, processing_func)
```

### 설정 관리

```python
from scmetabolism import get_config, set_config_value

# 기본 설정 확인
config = get_config()
print(f"기본 방법: {config.get('analysis.default_method')}")

# 설정 변경
set_config_value('analysis.default_method', 'ssgsea')
set_config_value('visualization.default_colormap', 'plasma')

# 설정 저장
config.save_config()
```

### 데이터 전처리

```python
from scmetabolism.utils import preprocess_data

processed_matrix = preprocess_data(
    count_matrix=raw_count_matrix,
    min_genes=200,  # 세포당 최소 유전자 수
    min_cells=3,    # 유전자당 최소 세포 수
    normalize=True, # CPM 정규화
    log_transform=True  # 로그 변환
)
```

### ALRA Imputation

```python
from scmetabolism.utils import alra_imputation

imputed_matrix = alra_imputation(count_matrix, k=50)
```

### 인터랙티브 시각화

```python
# Plotly 기반 인터랙티브 플롯
interactive_fig = plotter.interactive_dim_plot(
    embedding=umap_coords,
    metabolism_scores=metabolism_scores,
    pathway="Glycolysis / Gluconeogenesis",
    metadata=cell_metadata,
    color_by="cell_type"
)
interactive_fig.show()

# 네트워크 플롯
network_fig = plotter.pathway_network_plot(
    metabolism_scores=metabolism_scores,
    correlation_threshold=0.5
)
network_fig.show()

# 인터랙티브 히트맵
heatmap_fig = plotter.interactive_heatmap(
    metabolism_scores=metabolism_scores,
    pathways=pathways_of_interest
)
heatmap_fig.show()
```

### Gene Ontology (GO) 분석

```python
from scmetabolism import GOAnalysis

# GO 분석 객체 생성
go_analysis = GOAnalysis(organism="human")

# 대사 관련 GO 용어만 사용
metabolism_scores = sc_metab.compute_metabolism(
    count_matrix=count_matrix,
    method="aucell",
    metabolism_type="GO_metabolism"  # GO 대사 용어
)

# 생물학적 과정 (Biological Process) 분석
bp_scores = sc_metab.compute_metabolism(
    count_matrix=count_matrix,
    method="aucell", 
    metabolism_type="GO_BP"
)

# 직접 GO 유전자 세트 생성
go_gene_sets = go_analysis.create_go_gene_sets(
    aspects=["biological_process"],
    min_genes=10,
    max_genes=200
)

# 사용자 정의 GO 분석
sc_metab.gene_sets = go_gene_sets
custom_scores = sc_metab._compute_aucell(count_matrix, n_cores=2)
```

### 사용자 정의 유전자 세트

```python
custom_gene_sets = {
    "Custom_Pathway_1": ["GENE1", "GENE2", "GENE3"],
    "Custom_Pathway_2": ["GENE4", "GENE5", "GENE6"]
}

sc_metab.gene_sets = custom_gene_sets
scores = sc_metab._compute_aucell(count_matrix, n_cores=2)
```

## 🖥️ 명령줄 인터페이스 (CLI)

scMetabolism은 명령줄에서도 사용할 수 있습니다:

```bash
# 기본 분석
scmetabolism analyze --input data.csv --output results/ --method aucell

# GO 분석
scmetabolism analyze --input data.csv --gene-sets GO_metabolism --output results/

# 품질 관리
scmetabolism qc --input data.csv --output qc_report.html

# 설정 관리
scmetabolism config --show
scmetabolism config --set analysis.default_method=ssgsea
```

## 📖 예제

완전한 사용 예제는 `examples/basic_usage.py`를 참조하세요:

```bash
python examples/basic_usage.py
```

## 🧪 테스트

```bash
# 테스트 실행
pytest tests/

# 커버리지 포함
pytest tests/ --cov=scmetabolism
```

## 📚 인용

이 패키지를 사용하시면 다음 논문을 인용해 주세요:

**scMetabolism**
```
Yingcheng Wu, Shuaixi Yang, Jiaqiang Ma, et al. 
Spatiotemporal Immune Landscape of Colorectal Cancer Liver Metastasis at Single-Cell Level. 
Cancer Discovery. 2021.
```

**알고리즘 및 유전자 세트**
1. Aibar S, et al. AUCell: predicting transcription factor targets from single-cell RNA-seq data. Nat Methods. 2017.
2. Hänzelmann S, et al. GSVA: gene set variation analysis for microarray and RNA-seq data. BMC Bioinformatics. 2013.
3. Linderman GC, et al. Zero-preserving imputation of single-cell RNA-seq data. Nat Commun. 2022.

## 🤝 기여하기

1. Fork the repository
2. Create your feature branch (`git checkout -b feature/AmazingFeature`)
3. Commit your changes (`git commit -m 'Add some AmazingFeature'`)
4. Push to the branch (`git push origin feature/AmazingFeature`)
5. Open a Pull Request

## 📄 라이선스

이 프로젝트는 GPL-3.0 라이선스 하에 배포됩니다. 자세한 내용은 [LICENSE](LICENSE) 파일을 참조하세요.

## 📞 문의

- **기술적 문의**: [GitHub Issues](https://github.com/lemonardo1/scMetabolism/issues)
- **이메일**: gaoqiang@fudan.edu.cn

## 🙏 감사의 말

Original R package developers:
- Qiang Gao (gaoqiang@fudan.edu.cn)
- Yingcheng Wu (wuyc@mail.com)

Copyright (C) 2020-2024 Gao Lab @ Fudan University.