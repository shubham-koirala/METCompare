
# Deciphering Liver Metastases with Pan-Cancer Transcriptomics

We conducted a pan-cancer analysis of liver metastases, focusing on the transcriptional changes in cancer cells between primary and metastatic sites.

We first performed a systematic comparison of transcriptomic features across bulk and single-cell transcriptomics, highlighting the necessity of incorporating both types of transcriptomics. We then developed a novel tool called "DEBoost" for computing differentially expressed genes in cancer cells from bulk RNA-seq data. The comparison of nine cancer types revealed that cancers exhibit distinct transcriptional changes in liver metastasis while they share some patterns. Notably, cancer cells could partially mimic the secretome of hepatocytes by selectively expressing liver-specific genes in both primary and metastatic sites.

We also created a liver metastases single-cell atlas consisting of 750,000 cells. By comparing the single-cell transcriptome between primary and metastatic sites, we predicted drugs that could induce transformation from the primary state to the liver metastasis state.

Among these drugs, we discovered that the anti-diabetic drug sitagliptin might promote pancreatic cancer liver metastasis. Further analysis of drug adverse effect reports confirmed that sitagliptin has a significantly higher reporting odds ratio of pancreatic cancer and liver metastasis than other anti-diabetic drugs, including metformin.

## Repository Structure

There are three subfolders in the `Analysis` folder:

1. **BULK**
   - **SAMPLE_SELECTION**: Code used for selecting transcriptomic data of primary and metastatic sites for various cancers. The main datasets used are:
     - [TCGA](https://github.com/BioinformaticsFMRP/TCGAbiolinks/)
     - [MET500](https://xenabrowser.net/datapages/?cohort=MET500%20(expression%20centric))
     - [CCLE](https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE36133)
     - [GEPNET](https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE98894)
     - [Colorectal Cancer](https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE50760)
     - [Prostate Cancer](https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE147250)
       
   Some of these datasets are available in our OCTAD database ([link](https://github.com/Bin-Chen-Lab/octad)).

   - **DE (Differential Expression Analysis)**: Code used to perform DE analysis using the DEBoost method.
   - **DOWNSTREAM**: Pathway analysis done using differentially expressed genes.

2. **SINGLE_CELL**
   - Code files:
     - For building a liver metastasis single-cell atlas by integrating various single-cell datasets (scArches). The liver metastasis single-cell atlas can be downloaded from [here](https://chenlab-data-public.s3.us-west2.amazonaws.com/LIVER_METASTASIS_ATLAS/Chen_LiverMetastasis_new.RData).
     - Sample code for inferring malignant cells from the single-cell atlas (InferCNV).
     - Differential gene expression code used to generate the "Primary vs Metastatic signature" by comparing only malignant cells (predicted by InferCNV using the liver metastasis single-cell atlas).

3. **DRUG_STUDY**
   - Data and code files for:
     - Predicting drugs that could induce primary to metastasis transformation (utilizes signatures created using the liver metastasis single-cell atlas).
     - Investigating the possible role of the antidiabetic drug Sitagliptin in pancreatic cancer metastasis, which utilizes FDA Adverse Event Reporting System (FAERS) data.

All R scripts were tested on R (versions 4.1.1 and 3.5.1). Even R (version 3.5.1) is suitable to execute the codes.

Python (versions 3.9.1 and 3.9.2) was used for drug prediction and scanvi integration, respectively.




