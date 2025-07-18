# Contributing to scMetabolism

scMetabolism에 기여해 주셔서 감사합니다! 이 문서는 프로젝트에 기여하는 방법을 안내합니다.

## 🚀 시작하기

### 개발 환경 설정

1. Repository를 fork하고 clone합니다:
```bash
git clone https://github.com/your-username/scMetabolism-python.git
cd scMetabolism-python
```

2. 개발 환경을 설정합니다:
```bash
# 가상환경 생성 (권장)
python -m venv venv
source venv/bin/activate  # Linux/Mac
# 또는
venv\Scripts\activate  # Windows

# 개발 의존성과 함께 설치
pip install -e .[dev]
```

3. pre-commit hooks 설정 (선택사항):
```bash
pip install pre-commit
pre-commit install
```

## 📝 기여 방법

### 1. 이슈 생성

- 버그 리포트, 기능 요청, 질문 등은 [GitHub Issues](https://github.com/your-username/scMetabolism-python/issues)에 등록해 주세요.
- 이슈 템플릿을 사용하여 명확하고 상세한 정보를 제공해 주세요.

### 2. Pull Request 과정

1. **브랜치 생성**:
```bash
git checkout -b feature/your-feature-name
# 또는
git checkout -b bugfix/issue-number
```

2. **코드 작성**:
   - 코딩 스타일 가이드를 따라주세요
   - 테스트를 작성해 주세요
   - 문서를 업데이트해 주세요

3. **테스트 실행**:
```bash
# 모든 테스트 실행
pytest tests/

# 커버리지 확인
pytest tests/ --cov=scmetabolism --cov-report=html

# 코드 스타일 확인
flake8 scmetabolism
black --check scmetabolism
```

4. **커밋**:
```bash
git add .
git commit -m "Add: 새로운 기능 설명"
```

5. **Push 및 PR 생성**:
```bash
git push origin feature/your-feature-name
```

GitHub에서 Pull Request를 생성하고 템플릿을 따라 작성해 주세요.

## 🎯 코딩 가이드라인

### 코드 스타일

- **Python 스타일**: PEP 8을 따릅니다
- **포매터**: Black (line length: 88)
- **Linter**: flake8
- **타입 힌트**: 가능한 한 타입 힌트를 사용합니다

### 네이밍 컨벤션

- **함수/변수**: snake_case
- **클래스**: PascalCase
- **상수**: UPPER_SNAKE_CASE
- **Private 멤버**: _leading_underscore

### 문서화

- **Docstring**: Google 스타일을 사용합니다
- **주석**: 복잡한 로직에는 명확한 주석을 추가합니다

예시:
```python
def compute_metabolism(
    self,
    count_matrix: pd.DataFrame,
    method: str = "aucell",
    imputation: bool = False
) -> pd.DataFrame:
    """
    Compute metabolism scores for single cells.
    
    Args:
        count_matrix: Gene expression count matrix (genes x cells)
        method: Method for scoring ("aucell", "ssgsea", "gsva")
        imputation: Whether to perform ALRA imputation
        
    Returns:
        Metabolism scores (pathways x cells)
        
    Raises:
        ValueError: If method is not supported
    """
```

## 🧪 테스트

### 테스트 작성

- 새로운 기능에는 반드시 테스트를 작성해 주세요
- `tests/` 디렉토리에 `test_*.py` 형식으로 파일을 생성합니다
- pytest를 사용합니다

### 테스트 실행

```bash
# 전체 테스트
pytest

# 특정 파일
pytest tests/test_core.py

# 특정 함수
pytest tests/test_core.py::test_function_name

# 커버리지 포함
pytest --cov=scmetabolism
```

## 📚 문서화

### README 업데이트

- 새로운 기능을 추가했다면 README.md를 업데이트해 주세요
- 예제 코드를 포함해 주세요

### 예제 추가

- `examples/` 디렉토리에 실용적인 예제를 추가해 주세요
- Jupyter notebook 형태도 환영합니다

## 🐛 버그 리포트

버그를 발견하셨다면 다음 정보를 포함해서 이슈를 등록해 주세요:

- **환경 정보**: Python 버전, OS, 패키지 버전
- **재현 방법**: 최소한의 재현 가능한 예제
- **예상 결과**: 어떤 결과를 기대했는지
- **실제 결과**: 실제로 어떤 일이 일어났는지
- **에러 메시지**: 전체 traceback 포함

## 💡 기능 요청

새로운 기능을 제안하실 때는:

- **사용 사례**: 왜 이 기능이 필요한지
- **제안 구현**: 어떻게 구현할지에 대한 아이디어
- **대안**: 다른 해결 방법이 있는지

## 📋 체크리스트

Pull Request를 제출하기 전에 확인해 주세요:

- [ ] 코드가 PEP 8을 따르는가?
- [ ] 모든 테스트가 통과하는가?
- [ ] 새로운 기능에 테스트를 추가했는가?
- [ ] 문서를 업데이트했는가?
- [ ] 커밋 메시지가 명확한가?
- [ ] 변경사항이 기존 API를 깨뜨리지 않는가?

## 🏷️ 릴리스 프로세스

1. 버전 번호 업데이트 (`setup.py`, `pyproject.toml`, `__init__.py`)
2. CHANGELOG.md 업데이트
3. 태그 생성 및 GitHub Release
4. 자동으로 PyPI에 배포됨

## 📞 도움이 필요하신가요?

- [GitHub Issues](https://github.com/your-username/scMetabolism-python/issues)에 질문을 올려주세요
- [GitHub Discussions](https://github.com/your-username/scMetabolism-python/discussions)에서 토론해 주세요

## 🙏 감사합니다

모든 기여자분들께 감사드립니다! 여러분의 기여가 scMetabolism을 더 좋은 도구로 만들어 갑니다.