# Single-Cell Multiomics Tools

Work-in-progress list of tools for scRNA-seq and scATAC-seq data processing.

## Other lists of single-cell tools 
- [Awesome Single Cell](https://github.com/seandavi/awesome-single-cell)

## My favorites

### Preprocessing & QC
- [Cell Ranger](https://support.10xgenomics.com/) – 10x alignment & quantification
- [STARsolo](https://github.com/alexdobin/STAR) – fast droplet-based alignment
- [alevin-fry](https://github.com/COMBINE-lab/alevin-fry) – lightweight pseudoalignment

### Normalization & Integration
- [Seurat](https://satijalab.org/seurat/) – R ecosystem, CCA/RPCA integration
- [Scanpy](https://scanpy.readthedocs.io/) – Python ecosystem
- [scVI](https://scvi-tools.org/) – deep generative model for integration

### Trajectory / RNA Velocity
- [scVelo](https://scvelo.readthedocs.io/)
- [Monocle3](https://cole-trapnell-lab.github.io/monocle3/)

### Cell type proportions
- scCoda: identifies compositional changes (i.e. cell type proportions changes) in scRNA-seq datasets.
  - [Manuscript](https://www.nature.com/articles/s41467-021-27150-6)
  - [GitHub](https://github.com/theislab/scCODA)
  - [Documentation](https://pertpy.readthedocs.io/en/stable/tutorials/notebooks/sccoda_extended.html)
- tascCoda: extension of scCoda accounting a tree structure, such a cell lineage hierarchy or a taxonomic tree.
  - [Manuscript](https://www.frontiersin.org/journals/genetics/articles/10.3389/fgene.2021.766405/full)
  - [GitHub](https://github.com/bio-datascience/tascCODA)
  - [Tutorial](https://pertpy.readthedocs.io/en/stable/tutorials/notebooks/tasccoda.html)
- Milo: differential abundance analysis using KNN.
  - [Manuscript](https://www.nature.com/articles/s41587-021-01033-z)
  - [GitHub](https://github.com/MarioniLab/miloR) 
  - [Tutorial](https://pertpy.readthedocs.io/en/stable/tutorials/notebooks/milo.html)

### scATAC-seq Preprocessing
- [Cell Ranger ATAC](https://support.10xgenomics.com/)

### scTAC-seq Peak Calling & Annotation
- [SnapATAC2](https://github.com/kaizhang/SnapATAC2)
- [ArchR](https://www.archrproject.com/)
- [Signac](https://stuartlab.org/signac/) – Seurat companion for ATAC

## Multimodal Integration
- WNN (Seurat v4/v5)
- MOFA+
  - [Manuscript](https://link.springer.com/article/10.1186/s13059-020-02015-1)
  - [GitHub](https://github.com/bioFAM/MOFA2)
  - [Documentation](https://biofam.github.io/MOFA2/)
- [MultiVI](https://scvi-tools.org/)
  - [Manuscript](https://www.nature.com/articles/s41592-023-01909-9)
- scGlue
  - [Manuscript](https://www.nature.com/articles/s41587-022-01284-4)
  - [GitHub](https://github.com/gao-lab/GLUE)
  - [Documentation](https://scglue.readthedocs.io/en/latest/)

## Multimodal analysis
- SCENIC+: inference of enhancers, TFs and target to reconstruct gene regulatory networks (GRN)
  - [Manuscript](https://www.nature.com/articles/s41592-023-01938-4)
  - [GitHub](https://github.com/aertslab/scenicplus)
  - [Documentation](https://scenicplus.readthedocs.io/en/latest/)
- monaLisa: enrichment of TF motifs
  - [Manuscrip](https://academic.oup.com/bioinformatics/article/38/9/2624/6535228)
  - [GitHub](https://github.com/fmicompbio/monaLisa)
  - [Documentation](https://fmicompbio.github.io/monaLisa/articles/monaLisa.html) 

## Visualization
- [cellxgene](https://chanzuckerberg.github.io/cellxgene/)
- [UCSC Cell Browser](https://cells.ucsc.edu/)
