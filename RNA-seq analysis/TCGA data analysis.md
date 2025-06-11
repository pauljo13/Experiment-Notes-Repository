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
    
    - PC1, PC2는 전체 데이터의 분산(정보)을 최대한 잘 보존하면서도 **서로 직교(orthogonal)** 함
    - 보통 PC1 + PC2만 시각화해도 데이터의 **주된 구조**를 이해할 수 있어
