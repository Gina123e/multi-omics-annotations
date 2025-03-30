# Multi-omics annotations
We usually adopt different types of data from different sources/labs/technologies, such as ATAC-seq, m6A-seq, eCLIP-seq, icSHAPE-seq etc. Inegrating multi-omics data is paramount for people who conduct# research in bioinformatics. And the dataset we fetch is the *bed file format. It contains the absolute coordinates of each modification sites. But before training the deep learning model, we usually need to align thess datasets onto whole transcriptome, and get the relative coordinates onto the transcript. Here, I shared a pipeline to teach beginners in bioinformatics and compbio how to integrate different types of omics.

## Illustration
Fig1. 
![Annoate m6A sites on genome](../images/m6A genome.drawio.png)
## preparation
Firstly, you should download human annotation file-gencode.v38.annotation.gtf from ensembl and hg38.fa from NCBI. 
Then, you should use **star** to fetch transcript sequence. Notice! The RNA sequence of gene on minus strand must be reversed.

## modles
- eCLIP-seq + m6A integration.
```

```
- ATAC-seq + eCLIP overlap

## Example
See jupyter notebook

## License

This project is licensed under the MIT License - see the [LICENSE](LICENSE) file for details.
