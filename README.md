## Cohort Summary

The dataset consisted of human pancreatic islet RNA-seq samples derived from donors across multiple glycemic states. For this analysis, a subset of 57 samples was selected, including 39 Type 2 Diabetes (T2D) and 18 non-diabetic (ND) individuals.

### Limitations
- Observational dataset with no controlled intervention
- Tissue heterogeneity due to mixed islet cell populations
- Lack of longitudinal follow-up and survival outcomes
- Potential confounding from donor variability (age, clinical history)

## Biological Findings

### Differential Expression
Significant differentially expressed genes (DEGs) were identified between T2D and non-diabetic samples, indicating strong transcriptional differences associated with disease status.

### Biomarker Separation
A subset of top-ranked genes demonstrated strong discriminatory power between T2D and ND groups:
- PCA analysis showed clear clustering by disease group
- Logistic regression achieved high ROC-AUC, indicating predictive capability

### Enriched Pathways
Gene Ontology enrichment analysis revealed key biological processes:

- Glucose metabolism and insulin signaling pathways
- Beta-cell dysfunction and impaired insulin secretion
- Cellular stress response pathways
- Inflammatory signaling and immune activation

### Biological Interpretation
The findings are consistent with known T2D mechanisms:
- Reduced beta-cell function leading to impaired insulin production
- Increased metabolic stress and oxidative damage
- Activation of inflammatory pathways contributing to disease progression

## Dataset
- Source: GEO GSE164416
- Tissue: Human pancreatic islets
- Comparison: Type 2 Diabetes vs Non-Diabetic donors
- Analysis subset:
  - T2D: 39 samples
  - ND: 18 samples

## Methods
- Metadata cleaning using GEOquery
- Count filtering and sample matching
- Differential expression using DESeq2
- Gene annotation using org.Hs.eg.db
- Significant DEG filtering using adjusted p-value and log2 fold-change
- Biomarker model using logistic regression
- ROC-AUC model evaluation
- GO Biological Process enrichment using clusterProfiler

## Results Generated
- PCA plot
- Volcano plot
- Top DEG heatmap
- Significant and non-significant gene tables
- Top biomarker table
- ROC-AUC curve
- Biomarker PCA plot
- GO enrichment dotplot

  
This analysis demonstrates how transcriptomic data from observational datasets can be used to identify clinically relevant biomarkers. While the gene panel shows strong predictive performance, limitations such as lack of survival endpoints and confounding variables must be considered.

The pipeline reflects a translational approach where molecular signatures are integrated with statistical modeling to support disease classification and potential biomarker discovery.
