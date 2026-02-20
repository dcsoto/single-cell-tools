# Single-Cell Multiomics Tools

Work-in-progress list of my favorite tools for scRNA-seq and scATAC-seq data processing.

## Table of Contents
- [Tools](#tools)
  - [scRNA-seq preprocessing](#scrna-seq-preprocessing)
  - [scRNA-seq processing](#scrna-seq-processing)
  - [scRNA-seq batch correction](#scrna-seq-batch-correction)
  - [scRNA-seq doublet detection](#scrna-seq-doublet-detection)
  - [Trajectory / RNA Velocity](#trajectory--rna-velocity)
  - [Cell type proportions](#cell-type-proportions)
  - [scATAC-seq preprocessing](#scatac-seq-preprocessing)
  - [scATAC-seq processing](#scatac-seq-processing)
  - [Multimodal analysis](#multimodal-analysis)
  - [Visualization](#visualization)
  - [Multimodal integration](#multimodal-integration)
- [Tool ecosystems](#tool-ecosystems)
- [Other lists of single-cell tools](#other-lists-of-single-cell-tools)

## Tools

### scRNA-seq preprocessing
- Cell Ranger: 10x alignment & quantification
  - [Website](https://support.10xgenomics.com/)
- STARsolo: fast droplet-based alignment
  - [GitHub](https://github.com/alexdobin/STAR)

### scRNA-seq processing 
- Seurat: R ecosystem
  - [Website](https://satijalab.org/seurat/)
- Scanpy: Python ecosystem
  - [Documentation](https://scanpy.readthedocs.io/)

### scRNA-seq batch correction
- Harmony
- scVI

### scRNA-seq doublet detection
- DoubletFinder
  - [GitHub](https://github.com/chris-mcginnis-ucsf/DoubletFinder)  

### Trajectory / RNA Velocity
- scVelo
  - [Documentation](https://scvelo.readthedocs.io/)
- Monocle3
  - [GitHub](https://cole-trapnell-lab.github.io/monocle3/)

### Cell type proportions
- Speckle/Propeller
  - [Manuscript](https://academic.oup.com/bioinformatics/article/38/20/4720/6675456)
  - [GitHub](https://github.com/phipsonlab/speckle)
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

### scATAC-seq preprocessing
- [Cell Ranger ATAC](https://support.10xgenomics.com/)

### scATAC-seq processing
- Signac: Seurat companion for ATAC
  - [Website](https://stuartlab.org/signac/) 
- SnapATAC2
  - [GitHub](https://github.com/kaizhang/SnapATAC2)
- ArchR
  - [Website](https://www.archrproject.com/)

### Multimodal integration
- WNN (Seurat v4/v5)
- MOFA+
  - [Manuscript](https://link.springer.com/article/10.1186/s13059-020-02015-1)
  - [GitHub](https://github.com/bioFAM/MOFA2)
  - [Documentation](https://biofam.github.io/MOFA2/)
- MultiVI
  - [Manuscript](https://www.nature.com/articles/s41592-023-01909-9)
- scGlue
  - [Manuscript](https://www.nature.com/articles/s41587-022-01284-4)
  - [GitHub](https://github.com/gao-lab/GLUE)
  - [Documentation](https://scglue.readthedocs.io/en/latest/)
- Cobolt
  - [Manuscript](https://link.springer.com/article/10.1186/s13059-021-02556-z)
  - [GitHub](https://github.com/epurdom/cobolt)
  - [Tutorial](https://github.com/epurdom/cobolt/blob/master/docs/tutorial.ipynb)  

### Multimodal analysis
- SCENIC+: inference of enhancers, TFs and target to reconstruct gene regulatory networks (GRN)
  - [Manuscript](https://www.nature.com/articles/s41592-023-01938-4)
  - [GitHub](https://github.com/aertslab/scenicplus)
  - [Documentation](https://scenicplus.readthedocs.io/en/latest/)
- monaLisa: enrichment of TF motifs
  - [Manuscript](https://academic.oup.com/bioinformatics/article/38/9/2624/6535228)
  - [GitHub](https://github.com/fmicompbio/monaLisa)
  - [Documentation](https://fmicompbio.github.io/monaLisa/articles/monaLisa.html) 

### Visualization
- [cellxgene](https://chanzuckerberg.github.io/cellxgene/)
- [UCSC Cell Browser](https://cells.ucsc.edu/)

## Tool ecosystems 
- scverse: includes anndata, scanpy, SnapATAC2, and scvi-tools among others.
  - [Manuscript](https://www.nature.com/articles/s41587-023-01733-8)
  - [GitHub](https://github.com/scverse)
  - [Website](https://scverse.org/)
- scvi-tools: Python library for probabilistic modeling of single-cell data. Includes scVI and MultiVI.
  - [Manuscript](https://www.nature.com/articles/s41587-021-01206-w)
  - [GitHub](https://github.com/scverse/scvi-tools)
  - [Documentation](https://docs.scvi-tools.org/en/stable/index.html)
  - [Website](https://scvi-tools.org/)

## Other lists of single-cell tools 
- [Awesome Single Cell](https://github.com/seandavi/awesome-single-cell)
