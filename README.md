# Aberrant splicing of SGK1 in DLBCL

Analysis of RNA-Seq data from the GOYA trial to identify splice-mutant cases and 
look for the presence of aberrantly spliced SGK1 neoisoforms.  

**Paper:**  
Jie Gao, Eirini Sidiropoulou, Ieuan Walker, Joanna Krupka, Karol Mizielinski, 
Shamith Samarajiwa & Daniel J Hodson
[SGK1 mutations in DLBCL generate hyperstable protein neoisoforms that promote AKT independence]()  

### Background     

Serum and Glucocorticoid-regulated Kinase-1 (SGK1) is one of the most frequently 
mutated genes in Diffuse Large B Cell Lymphoma (DLBCL). However, little is known 
about its function or the consequence of its mutation. The frequent finding of 
truncating mutations has led to the widespread assumption that these represent 
loss-of-function variants and accordingly, that SGK1 must act as a tumour suppressor. 
Here we show that instead, the most common SGK1 mutations lead to production of 
**aberrantly spliced mRNA neoisoforms** in which translation is initiated from 
downstream methionines. 

Here, we are reanalysing RNA-Seq dataset of 553 DLBCL patients in order to 
find cases with SGK1 mutations and investigate the scale of aberrant splicing in 
these samples. Raw `FASTQ` files are avaliable from GEO with accession number 
[GSE125966](https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE125966) 
(date accessed 28 October 2019).     

### Workflow  

We utilized the open source GOYA DLBCL clinical trial – paired-end RNA-Seq dataset 
of 553 patient samples, which is publicly available from