# psiMNB

## Multinomial Naive Bayes for finding PUS-dependent pseudouridylation

![Alt text](psiMNB.png)

## Contents
- [pre-installation](#pre-installation)
- [test data](#test-data)
- [Usage](#Usage)
  - [generate training dataset](#generate-training-dataset)

### pre-installation
psiMNB is a training workflow for building Multinomial Naive Bayes (MNB) models of PUS label/probability prediction. psiMNB requires homer/motifFinding/motifStack pre-installation and predominantly used in unix-based operating systems. Therefore, for the usability of psiMNB, we recommend running homer/motifFinding/motifStack and the scripts (build_psi_MNB_overall_test_kmer.py/build_psi_MNB_overall.py/build_psi_MNB_TRUB1.py/build_psi_MNB_PUS3.py/build_psi_MNB_PUS1.py) in WSL2 (WSL2 installation guide: https://pureinfotech.com/install-windows-subsystem-linux-2-windows-10/) or unix-based system with R and python.

```R
#install motifStack in R
if (!require("BiocManager", quietly = TRUE))
    install.packages("BiocManager")

BiocManager::install("motifStack")
```
### test data
Test data: human_PUS_MNB_input_k-mer_overall.txt, human_PUS_MNB_input_k-mer_TRUB1.txt, human_PUS_MNB_input_k-mer_PUS3.txt, human_PUS_MNB_input_k-mer_PUS1.txt. (i.e. training dataset, result generated by homer/motifFinding/motifStack)

### Usage

#### generate training dataset
```shell
#run motifFinding: full_gene_contigs.bed6 is the bed6 file containing all the input Ψ-sites
./motifFinding ePSI_seq_total_polyA_Day0_mix_motifs/full_gene_HomerGenomeResults/full_gene_contigs.bed6 hg38 ePSI_seq_total_polyA_Day0_mix_motifs/full_gene_HomerGenomeResults -norevopp -noknown -rna -len 4,5,6,7,8,9 -p 20 -size given -dumpFasta
```

#### determine k-mer
```python
python build_psi_MNB_overall_test_kmer.py
```

#### build MNB model
```python
python build_psi_MNB_overall.py
python build_psi_MNB_TRUB1.py
python build_psi_MNB_PUS3.py
python build_psi_MNB_PUS1.py
```