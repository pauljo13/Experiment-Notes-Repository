---

# TCGA

---

# TCGA

---

TCGA(The Cancer Genome Atlas) 암 유전체 지도 사업부로 2006년 미국 국립암연구소(NCI)와 국립인간게놈연구소(NHGRI)가 주도해서 시작한 대규모 프로젝트

### 목적

- 다양한 암종의 유전체 정보를 체계적으로 수집하고 분석
- 암의 발생과 진행, 치료 반응 등을 유전체 수준에서 이해하려는 목적
- 개인 맞춤형 치료(Precision medicine)의 기반을 마련하는 데 기여

### 제공하는 주요 데이터 유형

|데이터 종류|설명|
|---|---|
|DNA 시퀀싱|유전자 돌연변이, 복제 수 변이 (CNV) 등|
|RNA 시퀀싱|유전자 발현량 (expression levels)|
|miRNA 시퀀싱|마이크로RNA 발현|
|DNA 메틸화|후성유전적 조절 메커니즘 정보|
|단백질 발현 (RPPA)|역상 단백질 배열 분석 (protein expression data)|
|임상 데이터|환자 생존율, 치료 반응, 병기 등|
|병리 이미지|H&E 슬라이드 이미지 등|

→ 대부분 Bulk-RNA sep

### TCGA 데이터의 특징

- 공개 데이터
- 표준화된 분석 파이프라인으로 품질이 우수
- 암종 간 비교 분석(pan-cancer analysis)이 가능

# TCGA data analysis

---

### TCGAbiolinks

R에서 TCGA 데이터를 검색, 다운로드, 전처리, 분석할 수 있게 도와주는 매우 강력한 패키지

### 설치

```r
if (!requireNamespace("BiocManager", quietly = TRUE))
    install.packages("BiocManager")

BiocManager::install("TCGAbiolinks")

library(TCGAbiolinks)
```

### 데이터 검색

```r
query <- GDCquery(
    project = "TCGA-BRCA",              # 분석할 암종
    data.category = "Transcriptome Profiling",  
    data.type = "Gene Expression Quantification",  
    workflow.type = "HTSeq - Counts"
)
```

- `project` : TCGA 암종 어떤 것을 분석할지, 프로젝트 이름
    
- `data.category` : 데이터를 크게 분류하는 상위 카테고리
    
    |`data.category`|설명|대표 사용 예시|
    |---|---|---|
    |**Transcriptome Profiling**|RNA-seq 기반 유전자 발현 데이터|유전자 발현 분석, DEG, 클러스터링|
    |**Genomic Alteration** (구 `Simple Nucleotide Variation`)|돌연변이(SNV), Indel 정보|변이 탐색, 바이오마커 분석|
    |**Copy Number Variation**|유전자 복제수 변이(CNV) 데이터|CNV 기반 암 유전자 탐색|
    |**DNA Methylation**|DNA 메틸화 패턴 정보|후성유전학적 조절 연구|
    |**Clinical**|환자 정보 (진단, 병기, 생존 등)|생존 분석, 임상적 특성 비교|
    |**Biospecimen**|샘플 준비 정보 (조직 수집, 처리 등)|전처리 QC|
    |**Proteome Profiling**|단백질 발현 데이터 (RPPA 등)|단백질 기반 바이오마커 연구|
    |**Imaging**|병리 이미지 (H&E 슬라이드)|조직학 분석, AI 기반 이미지 분석|
    |**Raw sequencing data**|FASTQ 등의 원시 시퀀싱 데이터|재분석, 맞춤형 파이프라인 적용|
    |**Structural Variation**|게놈 구조적 변이 (SV) 정보|유전체 Rearrangement 연구|
    
- `data.type` : 실험 형태
    
    - **Transcriptome Profiling**
        
        |`data.type`|설명|
        |---|---|
        |**Gene Expression Quantification**|유전자 단위의 발현량 데이터|
        |Isoform Expression Quantification|전사체(isoform) 단위의 발현량|
        |miRNA Expression Quantification|마이크로 RNA 발현량|
        
    - **Genomic Alteration**
        
        |`data.type`|설명|
        |---|---|
        |**Masked Somatic Mutation**|암세포 특이적 돌연변이 데이터 (보안 처리된)|
        |Aggregated Somatic Mutation|TCGA 전체에서 합쳐진 돌연변이|
        |Annotated Somatic Mutation|주석이 붙은 SNV/INDEL 목록|
        
    - **Copy Number Variation**
        
        |`data.type`|설명|
        |---|---|
        |**Copy Number Segment**|게놈 상의 복제수 변화 영역 정보|
        |GISTIC2 Copy Number|CNV 분석 툴 GISTIC2 결과 (이득/손실 등)|
        
    - **DNA Methylation**
        
        |`data.type`|설명|
        |---|---|
        |**Methylation Beta Value**|프로브별 메틸화 수치 (0~1 사이)|
        |Methylation Raw Value|비정규화된 메틸화 수치|
        
    - **Clinical**
        
        |`data.type`|설명|
        |---|---|
        |**Clinical Supplement**|임상 부가정보 (치료, 병리 등)|
        |Clinical Biospecimen Supplement|샘플 관련 임상 정보|
        |Clinical Drug|약물 투여 관련 정보|
        
- `workflow.type` : 데이터 처리 방법
    
    - `"HTSeq - Counts”` : HTSeq라는 툴로 read를 유전자에 카운트
    - `"HTSeq - FPKM"` : 길이 정규하된 발현량 (유전자 길이와 전체 read 수)
    - `"HTSeq - TPM"` : FPKM 보다 상대 비교에 적합함

### 데이터 다운로드

```r
GDCdownload(query)
```

### 데이터 준비

```r
data <- GDCprepare(query)
```

## SummarizedExperiment

```r
data 구조
# class : SummarizedExperiment
# dim : (gene 수, sample 수)
```

### assay(data)

- 유전자 발현 행렬 (또는 메틸화 수치, 카운트 등)
- rownames : gene ID
- colnames : sample ID

```r
assay(data)
df <- as.data.frame(assay(data))
```

### colData(data)

- 샘플 메타데이터(임상 정보, 조직 정보 등)
- 각 열이 샘플 하나

```r
colData(data)
phnoe_df <- as.data.frame(colData(data))
```

### rowData(data)

- 유전자 메타데이터(위치, 심볼 등)
- gene_id, gene_name, chromosome

```r
rowData(data)
```

## EDA

---

### 1. Annotation

gene expression data의 유전자는 유전자 이름이 아닌 데이터를 만드는 장치에서 사용하는 id 코드(Ensembl Gene ID)로 되어 있다. 그렇기에 이를 데이터 분석에 편리하게 gene symbol로 바꿔줘야 한다.

```r
# BiocManager::install("biomaRt") 먼저 설치 필요
library(biomaRt)

# Ensembl에서 사람 유전자 데이터 가져오기
mart <- useMart("ensembl", dataset = "hsapiens_gene_ensembl")

# 변환 실행
genes <- c("ENSG00000141510", "ENSG00000139618")
getBM(
  attributes = c("ensembl_gene_id", "hgnc_symbol"),
  filters = "ensembl_gene_id",
  values = genes,
  mart = mart
)
```

이와 같은 방법으로 바꿀 수도 있지만,

```r
genes <- as.data.frame(rowData(data)) # 유전자 메타데이터

gene_names <- genes$gene_name[match(rownames(df), genes$gene_id)]
gene_names[is.na(gene_names)] <- "NA"
rownames(df) <- make.unique(gene_names)
```

이런식으로 data에 있는 데이터를 사용하여 바꿀 수도 있다.

### 2. Batch

**같은 조건의 샘플이어도**, 실험을 수행한 **시간, 장소, 장비, 실험자** 등의 **비생물학적 요인** 때문에 생기는 **데이터 간 차이**를 말한다.

이 효과는 우리가 실제로 알고자 하는 생물학적 차이(예: 암 vs 정상)와 **무관한 "잡음"**이기 때문에,

제거하지 않으면 분석 결과가 왜곡될 수 있어.

- Batch 효과의 영향
    - 가짜 차이(false positive)를 유발하거나,
    - 진짜 차이를 **가려버릴 수도** 있음
    - 클러스터링, PCA, DEG 결과가 왜곡됨

```r
batch <- pheno_df$preservation_method
mod <- model.matrix(~1, data = pheno_df)
log_df <- log2(df + 1)
df_combat <- ComBat(dat = as.matrix(log_df), batch = batch, mod = mod)
```

- ComBat
    
    ```r
    library(sva)
    
    # 데이터: exprs (gene x sample), batch 정보, group 정보
    combat_exprs <- ComBat(dat = exprs, batch = batch_vector, mod = model.matrix(~group))
    
    ```
    
- limma
    
    ```r
    library(limma)
    
    clean_exprs <- removeBatchEffect(exprs, batch = batch_vector, design = model.matrix(~group))
    
    ```
    
- PCA
    
    ```r
    library(uwot)
    umap_result <- umap(t(clean_exprs))
    plot(umap_result, col = batch_vector)
    ```
    

|방법|배치를 "진짜로 제거"하나?|설명|
|---|---|---|
|**PCA**|❌ No|단지 시각화/진단용|
|**PCA 후 PC 제거**|⚠️ 제한적|PC1 등 제거 → 정보 손실 가능성 있음|
|**ComBat / removeBatchEffect**|✅ Yes|배치 정보를 모델링해서 제거함|

### 3. DESeq2

**DESeq2**는 RNA-seq 등의 count 기반 데이터에서 **차등 발현 유전자(DEG)** 를 찾기 위한 **Bioconductor 패키지**

**두 그룹(또는 조건) 간 유전자 발현 차이(log2 Fold Change)** 를 **통계적으로 신뢰성 있게 검정**하는 것

**순서**

- **Raw count 입력**
    - 예: gene × sample 형태의 매트릭스
- **정규화 (Normalization)**
    - 각 샘플 간 **시퀀싱 깊이(library size)** 차이를 보정
    - `estimateSizeFactors()`
- **분산 추정 (Dispersion estimation)**
    - 유전자별로 발현의 **변동성**(생물학적 + 기술적)을 추정
- **통계 모델 적합**
    - **음이항 분포 모델**을 기반으로 조건 간 차이를 평가
- **DEG 결과 추출**
    - log2 Fold Change
    - p-value 및 FDR (adjusted p-value)

시각화 중심 or 다양한 배치 요인 있을 때:

```r
DESeq2 정규화
→ log2 변환 or vst()
→ ComBat 등으로 배치 보정
→ PCA / UMAP / clustering
```

DEG 중심이면:

```r
DESeq2 with design = ~ batch + condition
→ 바로 DEA 수행
```

**주요 함수**

|함수|설명|
|---|---|
|`DESeqDataSetFromMatrix()`|분석용 객체 생성|
|`DESeq()`|전체 분석 수행 (정규화 + 모델 적합)|
|`results()`|DEG 결과 추출|
|`plotMA()`|DEG 시각화|
|`lfcShrink()`|log2FC 안정화 (작은 샘플일수록 추천)|

**DESeqDataSetFromMatrix**

|인자 이름|설명|예시|
|---|---|---|
|**`countData`**|유전자 발현 **raw count matrix**(행 = 유전자, 열 = 샘플)|`df`|
|**`colData`**|각 샘플에 대한 **메타데이터** (DataFrame, 열 수 = countData 열 수)|`pheno_df`|
|**`design`**|분석할 **모델 공식 (formula)**조건, 배치 등 지정|`~ condition`, `~ batch + condition`|
|**`tidy`** _(기본: FALSE)_|`countData`가 tidy 형태일 경우 TRUE로 설정(보통 사용 X)|사용하지 않음|

### 4. VST

**VST(Variance Stabilizing Transformation는 RNA-seq count 데이터에서 유전자 간 발현의 분산을 안정화(stabilize) 시켜주는 변환 방법 → DESeq2 패키지에서 제공**

**RNA-seq 데이터의 원래 문제:**

- count 데이터는 **low count일수록 분산이 작고**, **high count일수록 분산이 큼**
    
- 즉, **분산이 발현값에 따라 종속됨 (heteroscedasticity)**
    
    → 이 상태로는 PCA나 UMAP 같은 분석에서 **과하게 표현되거나 묻히는 유전자**가 생김
    

**VST**

- log2 변환처럼 보이지만, **발현량에 따라 달라지는 분산을 수학적으로 평탄화**해줌
- 결과적으로 각 유전자의 값이 **비슷한 스케일에서 비교 가능**해짐

**원리**

**음이항 분포 기반으로 추정된 분산 함수를 이용해 변환함**

→ 단순히

```
log2(count + 1)
```

보다 더 정교함

|항목|설명|
|---|---|
|🎯 목적|분산 안정화 (PCA/UMAP에 적합)|
|💡 특징|count 값이 크든 작든 **비슷한 수준의 변동성 유지**|
|🔄 결과|**정규화 + log + 분산 보정**된 continuous 데이터|
|📊 사용처|PCA, clustering, UMAP, heatmap 등 시각화|
|🧪 DEG 분석|❌ 사용 X → raw count + DESeq2 사용해야 함|

|방법|차이|
|---|---|
|`log2(normalized + 1)`|단순 log 변환 → 분산 불안정|
|`vst(dds)`|정규화 + 분산 안정화 → 더 부드러운 분포|

### 유전자 추출 분석

무조건 _DESeq2 + VST + ComBat 이후_ 에 유전자 필터링을 해야 한다.

```r
Raw Count
↓
DESeq2 정규화 (size factor)
↓
VST 또는 log2(normalized + 1)
↓
ComBat (batch 정정)
↓
→ 특정 유전자 추출 (e.g., DEGs, top variable genes)
↓
→ PCA / UMAP / clustering 등 시각화
```

정규화(VST 포함)와 batch 정정은 전체 유전자에서 수행해야 하고,

그 이후에 필요한 유전자만 선택해서 PCA하는 것이 가장 정확하다.

# Analysis

---

## PCA analysis

<aside> 💡

고차원 데이터를 주요한 방향(성분)으로 압축해서,

데이터의 구조나 패턴을 시각적으로 파악할 수 있도록 해주는 차원 축소 기법이다.

</aside>

- 왜 쓰는가?
    
    1. **데이터 시각화 (2D/3D)**
        - 수천 개 유전자를 몇 개의 주성분으로 축약 → 2D 평면에서 군집 확인 가능
    2. **샘플 간 패턴 파악**
        - 암 vs 정상, batch 간 차이, outlier 탐색
    3. **노이즈 제거 / 특성 축소**
        - 분석 효율을 높이고 과적합 방지
- 동작법
    
    PCA는 **데이터에서 분산이 가장 큰 방향**(즉, 가장 정보가 많은 방향)을 찾아서 그걸 **PC1 (1st principal component)**라고 하고, 그 다음 큰 방향을 **PC2**, 그다음을 **PC3** … 이런 식으로 나눠.
    
    - PC1, PC2는 전체 데이터의 분산(정보)을 최대한 잘 보존하면서도 서로 직교(orthogonal)**함
    - 보통 PC1 + PC2만 시각화해도 데이터의 **주된 구조**를 이해할 수 있어

```r
variances <- apply(df_combat, 1, var)
zero_var_genes <- names(variances[variances == 0 | is.na(variances)])
constant_cols <- apply(df_combat, 1, function(x) length(unique(x)) == 1)
df_combat_filtered <- df_combat[variances > 0, ]

df_t <- t(df_combat_filtered)
df_scaled <- scale(df_t)

pca_res <- prcomp(df_scaled, scale. = F, center = T)
```

### 시각화

Scree plot (주성분 분산 설명 비율 시각화)

```r
plot(pca_res, type = 'l', main = 'Scree Plot')
```

![image.png](attachment:e447120e-f167-4e78-abd4-9bf334319438:image.png)

Plot

```r
pca_df <- as.data.frame(pca_res$x)
pca_df$group <- pheno_df$shortLetterCode[match(rownames(pca_df), pheno_df$barcode)]

plot(pca_df$PC1, pca_df$PC2,
     col = as.factor(pca_df$group),
     pch = 19,
     xlab = "PC1", ylab = "PC2",
     main = "PCA: PC1 vs PC2")
legend("topright", legend = unique(pca_df$group), col = 1:length(unique(pca_df$group)), pch = 19)

```

![image.png](attachment:504c67e8-0e0f-42b2-b174-9e3659a3538a:image.png)

ggplot2 활용

```r
library(ggplot2)

ggplot(pca_df, aes(x = PC1, y = PC2, color = group)) +
  geom_point(size = 2, alpha = 0.8) +
  theme_minimal() +
  labs(title = "PCA: PC1 vs PC2",
       x = paste0("PC1 (", round(summary(pca_res)$importance[2, 1]*100, 1), "%)"),
       y = paste0("PC2 (", round(summary(pca_res)$importance[2, 2]*100, 1), "%)")) +
  theme(plot.title = element_text(hjust = 0.5))
```

![image.png](attachment:2d8b1c68-681d-4d81-9b26-5a71938c26b7:image.png)

## UMAP

---

UMAP(Uniform Manifold Approximation and Projection)은 고차원 데이터를 저차원으로 효과적으로 **시각화**하거나 **차원 축소**하는 데 사용되는 기법이야. 비슷한 방식으로 쓰이는 t-SNE와 비교되기도 하는데, UMAP은 **더 빠르고**, **더 잘 일반화**된다는 장점이 있다.

_핵심 개념_

- Manifold assumption : 데이터는 고차원 공간에 존재하지만, 그 본질은 더 낮은 차원의 다양체(manifold)에 있다는 가정.
- Topology preservation : 데이터 간의 국소적 구조(local structure)를 최대한 본존하면서 저 차원으로 투영함

**원리**

1. 고차원 공간에서 neighbors 관계 그래프 구성
    
    → 유사한 샘플 간 연결성을 확률로 표현
    
2. 저차원 공간에서 비슷한 연결성 가지도록 최적화
    
    → 고차원 이웃 관계가 저차원에서도 잘 보존되도록 조정
    

**장점**

- 속도 빠름
- 전역 구조(Global structure)도 어느 정도 보존 가능
- Supervised UMAP: 라벨 정보도 활용 가능
- scikit-learn과 잘 통합됨(umap-learn)

|항목|UMAP|t-SNE|
|---|---|---|
|📌 기본 철학|Topology 보존 (거리 + 연결 구조)|확률 분포 보존 (국소 거리)|
|⚙️ 수학 기반|리만 다양체 이론 + 범주론|조건부 확률 분포 (KL divergence 최소화)|
|⏱ 속도|빠름 (대용량에 적합)|느림 (특히 샘플 수 많을 때)|
|🔁 재현성|동일 파라미터로 반복 가능|초기값 의존 → 결과 변동 있음|
|🌐 전역 구조 보존|**상대적으로 우수함**|전역 구조는 잘 안 보존됨|
|🔄 역변환 가능성|이론적으로 가능 (semi-invertible)|불가능|
|🎯 지도 학습 가능|가능 (Supervised UMAP)|불가능|
|🧪 하이퍼파라미터|`n_neighbors`, `min_dist` 등 (의미 해석 쉬움)|`perplexity` (직관적이지 않음)|

**설치**

```r
install.packages("uwot")
```

**기본 사용**

```r
library(uwot)

umap_result <- umap(expression_matrix)
plot(umap_result, col = sample_group, pch = 19, main = "UMAP")
```

- `n_neighbors` : 국소 구조를 얼마나 고려할지 (기본 15)
    
    |설정|설명|
    |---|---|
    |`5 ~ 15`|**작은 국소 군집** 강조 (작은 클러스터 잘 드러남)|
    |`15 ~ 30`|**일반적으로 많이 사용**, 균형 좋음|
    |`30 ~ 100`|전체적인 구조 중시, 클러스터 경계가 덜 명확해질 수 있음|
    
- `min_dist` : 군집 간 거리 조절 (작을수록 조밀한 클러스터)
    
    |설정|설명|
    |---|---|
    |`0.001 ~ 0.1`|매우 조밀한 군집 (빽빽하게 표현됨)|
    |`0.1 ~ 0.5`|기본 추천 범위 (군집 구분 잘됨)|
    |`0.5 ~ 0.9`|군집 사이 간격이 넓어짐 (덜 조밀, 더 퍼짐)|
    
- `metric` : 거리 측정 방식 (`"euclidean"`, `"cosine"` 등)
    

|파라미터|설명|추천 범위 / 기본값|
|---|---|---|
|`n_neighbors`|이웃의 개수 (local vs global 구조 결정)|`5 ~ 100` (기본: `15`)|
|`min_dist`|점들 간 최소 거리 (군집의 조밀도 조절)|`0.001 ~ 0.9` (기본: `0.1`)|
|`metric`|거리 계산 방식|`"euclidean"` (기본), `"cosine"`, `"manhattan"`, `"correlation"` 등|
|`n_components`|축소할 차원 수|`2` (기본), 3도 가능|
|`n_epochs`|학습 반복 횟수|기본: 500 (`uwot`), 200 (`Python`)|
|`init`|초기 임베딩 방법|`"spectral"` (기본), `"random"`|
|`scale`|입력 데이터 정규화 여부|`FALSE` (기본)|
|`spread`|전체 클러스터 분산도 조절|기본: `1.0`|
|`set_op_mix_ratio`|local vs global 연결 균형|기본: `1.0` (local 우선)|
|`random_state` / `set.seed()`|랜덤성 제어 (결과 재현용)|예: `set.seed(42)`|

**PCA & UMAP**

|방식|요약|추천 상황|
|---|---|---|
|PCA 후 UMAP|노이즈 제거 + 속도 향상|유전자 수 많을 때|
|UMAP만|모든 정보 반영|소규모 데이터 / 변수 이미 줄였을 때|

**DESeq2 & PCA**

|항목|**DESeq2 정규화**|**PCA (Principal Component Analysis)**|
|---|---|---|
|🎯 목적|**라이브러리 크기/샘플 간 발현량 차이 조정**|**고차원 데이터를 저차원으로 요약**|
|⚙️ 처리 대상|**raw count** → normalized count|**normalized data**|
|🧪 수학적 기법|scaling by size factor|고차원 → 주성분 (선형 조합)|
|🧼 역할|데이터 전처리 (정규화 단계)|데이터 요약 (시각화/차원 축소용)|
|📦 사용 시기|분석 초반부 (DEA 전 필수)|분석 후반부 (시각화 등)|

- **DESeq2 정규화**:
    
    → 각 샘플마다 시퀀싱 깊이, 기술적 편차 등을 **보정**해서 비교 가능한 발현값을 만드는 과정.
    
    👉 _데이터를 공정하게 만드는 전처리_
    
- **PCA**:
    
    → 정규화된 데이터를 바탕으로 가장 큰 변동성 축을 찾아 **요약/시각화**하는 분석.
    
    👉 _데이터를 2~3차원으로 줄여서 구조를 파악_
    

**DESeq2 & log2**

|작업|기능|정규화 포함 여부|
|---|---|---|
|`log2(counts + 1)`|분포 안정화만|❌|
|`DESeq2 정규화`|시퀀싱 깊이 등 보정|✅|
|`vst()` / `rlog()`|정규화 + 분포 안정화|✅✅|

## DEA (Differential Expression Analysis)

---

## GSEA (Gene Set Enrichment Analysis)

---

## Co-expression Network Analysis

---

## Immune Cell Infiltration Estimation

---

## Cell Type Deconvolution

---

## Transcription Factor Activity Inference