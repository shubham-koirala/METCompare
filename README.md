Deciphering liver metastases with pan-cancer transcriptomics

We conducted a pan-cancer analysis of liver metastases, focusing on the transcriptional changes in cancer cells between primary and metastatic sites. 

We first performed a systematic comparison of transcriptomic features
across bulk and single-cell transcriptomics, highlighting the necessity of incorporating both
types of transcriptomics. We then developed a novel tool “DEBoost” for computing differentially expressed genes in cancer cells from bulk RNA-seq data.  

The comparison of nine cancer types revealed that cancers exhibit distinct transcriptional changes in liver metastasis while they   share   some   patterns.  Notably, cancer   cells   could   partially   mimic   the  secretome   of hepatocytes by selectively expressing liver-specific genes in both primary and metastatic sites. Among these, CPN1, a liver-specific gene, was found to be highly expressed in colorectal cancer tissues and cell lines. In vivo experiments demonstrated that overexpression of CPN1 promoted colorectal cancer liver metastasis. 

We also created a liver metastases single-cell atlas consisting of 750,000 cells. By comparing the single cell transcriptome between primary and metastatic site, we predicted drugs that could induce transformation from primary state to liver metastasis state. Among these drugs, we   discovered   that   the   anti-diabetic   drug   sitagliptin   might   promote pancreatic cancer liver metastasis. 

Further analysis of drug adverse effect reports confirmed that sitagliptin has a significantly higher reporting odds ratio of pancreatic cancer and liver metastasis than other anti-diabetic drugs including metformin. 

Repository structure
There are three subfolders in Anlaysis folder.
1.	BULK

•	SAMPLE_SELECTION

It includes code used for selecting transcriptomic data of primary and metastatic site for 
various cancers. The main datasets used are: 
          
a.	TCGA (https://github.com/BioinformaticsFMRP/TCGAbiolinks)
b.	MET500 (https://xenabrowser.net/datapages/?cohort=MET500%20(expression%20centric)

c.	CCLE(https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE36133)
d.	GEPNET (https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE98894)
e.	Colorectal Cancer (https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE50760)
f.	Prostate Cancer (https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE147250)


Some of these datasets are available in our OCTAD database. https://github.com/Bin-Chen-Lab/octad.

•	DE (Differential Expression Analysis)

It includes code used to perform DE analysis using DEBoost method.

•	DOWNSTREAM (ANALYSIS)

It includes pathway analysis done using differentially expressed genes.

2.	SINGLE_CELL

It includes code files:

a.	for building liver metastasis single cell atlas by integrating various single cell datasets. (scArches)

b.	sample code for inferring malignant cells from single cell atals (InferCNV)
c.	differential gene expression code used to generate ‘Primary vs Metastatic singnature’ by comparing only malignant cells (predicted by InferCNV using liver metastasis single cell atlas)

     The liver metastasis single cell atlas can be downloaded from  
https://chenlab-data-public.s3.us-west2.amazonaws.com/LIVER_METASTASIS_ATLAS/Chen_LiverMetastasis_new.RData.  


3.	DRUG_STUDY

It includes data and code files:

a.	For predicting drugs which could induce primary to metastasis transformation. 
(It utilizes signatures created using liver metastasis single cell atlas)

b.	Possible role of antidiabetic drug Sitagliptin in Pancreatic cancer metastasis which utilizes FDA Adverse Event Reporting System (FAERS) data.


All R scripts were tested on R (version 4.1.1). Even R (version 3.5.1) is suitable to execute the codes. 

Python (version 3.9.1 and 3.9.2 ) were used for drug prediction and scanvi integration respectively. 



