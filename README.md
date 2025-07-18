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
pip install git+https://github.com/your-username/scMetabolism-python.git
```

### 로컬 개발 설치

```bash
git clone https://github.com/your-username/scMetabolism-python.git
cd scMetabolism-python
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

### 선택적 패키지
```bash
# Scanpy 지원
pip install "scmetabolism[scanpy]"

# 모든 기능
pip install "scmetabolism[all]"
```

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
    metabolism_type="KEGG"  # "KEGG" 또는 "REACTOME"
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

## 🔧 고급 사용법

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

### 사용자 정의 유전자 세트

```python
custom_gene_sets = {
    "Custom_Pathway_1": ["GENE1", "GENE2", "GENE3"],
    "Custom_Pathway_2": ["GENE4", "GENE5", "GENE6"]
}

sc_metab.gene_sets = custom_gene_sets
scores = sc_metab._compute_aucell(count_matrix, n_cores=2)
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

- **기술적 문의**: [GitHub Issues](https://github.com/your-username/scMetabolism-python/issues)
- **이메일**: gaoqiang@fudan.edu.cn

## 🙏 감사의 말

Original R package developers:
- Qiang Gao (gaoqiang@fudan.edu.cn)
- Yingcheng Wu (wuyc@mail.com)

Copyright (C) 2020-2024 Gao Lab @ Fudan University.