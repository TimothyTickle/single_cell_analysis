

Single-cell RNA-Seq Analysis
========================================================
author: Timothy Tickle and Brian Haas
css: single_cell_analysis.css
date: October 10, 2016

We Know Tissues are Complex
===

![complex](images/complex_tissues.png)

We can Separate Populations of Cells
===

<div align="center">
<img src="images/dropseq_paper.png" height=500>
</div>

We can Separate Populations of Cells
===

<div align="center">
<img src="images/dropseq_seperating_cells.png" height=500>
</div>

Heterogeneity can Help Describe Cell Type
===

<div align="center">
<img src="images/karthik_retina_paper.png" height=500>
</div>

Heterogeneity can Help Describe Cell Type
===

<div align="center">
<img src="images/shekhar_cell_type.png" height=500>
</div>

Heterogeneity can Describe Phenotypes
===

<div align="center">
<img src="images/itay_paper.png" height=500>
</div>

Heterogeneity can Describe Phenotypes
===

<div align="center">
<img src="images/single_cell_clinical.png" height=500>
</div>

Before We Get Started
===

- Single-cell RNA-Seq analysis is developing.
  - Give you a feel for the data.
  - Walk different analysis in a couple data sets.
- There is a vivid diversity of methodology
  - These techniques will grow as the field does.
  - Why these tools were chosen.
- This is a guided conversation through scRNA-Seq Analysis
  - Breadth and targeted depth
  - If you have expertise to add, let's have that conversation!

Before We Get Started
===

- This is VERY hands-on.
  - Much can be applied to other analyses.
  - Strengthen those R ninja skills!
  - If you need, cut and pasting is available (cut_and_paste.txt).
  - There will many cute corgi pictures.

---

![ninja corgis](images/corgis/ninja_corgi.jpeg)

What We Will Attempt to Cover
===

### Section 1  
__Part 1:__ Overview of laboratory prep and sequence analysis.  
__Part 2:__ Characteristics of expression data and QC.  

### Section 2  
__Part 1:__ Plotting Single Cell RNA-Seq data.  
__Part 2:__ Evaluating and defining cell populations.  

Welcome to Section 1
===

- Overview of assays (SmartSeq2, Dropseq, 10X).
- Overview of sequence data (SmartSeq2, Dropseq, 10X).
- Overview of Dropseq sequence analysis pipeline.
- Overview of 10X sequence analysis pipeline.
- Sequence level quality control.
- What does single cell expression data look like and why?
- Introduction to RStudio.
- Initial data exploration.
- Quality control for expression matrices.
- Filtering genes and cells.
- Data analysis using UMIs.

Overview of assays
===

<div align="center">
<img src="images/corgis/scientist.jpg" height=500>
</div>

SmartSeq2
===
  
![smartseq_paper](images/smartseq_paper.png)

---

![smartseq_protocols](images/smartseq_protocols.png)

SmartSeq2: Description
===

- Developed for single cell but can performed using total RNA.
- Selects for poly-A tail.
- Full transcript assay.
  - Uses template switching for 5' end capture.
- Standard illumina sequencing.
  - Off-the-shelf products.
- Hundreds of samples.
- Often do not see UMI used.

SmartSeq2: Assay Overview
===

<div> 
<img src="images/smart_seq2.png" height=400>
</div>

---

- Poly-A capture with 30nt polyT and 25nt 5' anchor sequence.
- RT adding untemplated C
- Template switching
  - Locked Nucleic Acid binds to untemplated C
  - RT switches template
- Preamplification / cleanup
- DNA fragmentation and adapter ligation together.
- Gap Repair, enrich, purify

SmartSeq2: Equipment
===

![sequencer](images/nextseq-large.jpg)

---

![plate](images/plate.jpg)

Dropseq
===

<div>
<img src="images/dropseq_paper.png" height=500>
</div>

Dropseq: Description
===

- Moved throughput from hundreds to thousands.
- Droplet-based processing using microfluidics
  - Nanoliter scale aqueous drops in oil.
  - 3' End
- Bead based (STAMPs).
  - Single-cell transcriptomes attached to microparticles.
  - Cell barcodes use split-pool synthesis.
- Uses UMI (Unique Molecular Identifier).
  - RMT (Random Molecular Tag).
  - Degenerate synthesis.

Dropseq: Assay Overview
===

[Click Here for Dropseq Video Abstract](http://www.sciencedirect.com/science/article/pii/S0092867415005498)

Dropseq: Assay and UMI Construction
===

![drop_umi](images/dropseq_assay_umi.png)   
Macosko et al. 2015

Dropseq: Sequencing Overview
===

<div>
<img src="images/dropseq_assay_1.png" height=500>
</div>

Dropseq: Microfluidics
===

![drop_beads](images/dropseq_drops.gif)

Dropseq: Equipment
===

![dropseq_eq](images/dropseq_equipment.png)

Dropseq: Pointers
===

- Droplet-based assays can have leaky RNA.
- Before library generation wash off any medium (inhibits library generation).
- Adding PBS and BSA (0.05-0.01%) can protect the cell.
  - Too much produces a residue making harvesting the beads difficult.
- Filter all reagent with a 80 micron strainer before microfluidics.
- Some purchased devices add a hydrophobic coating.
  - Can deteriorate (2 months at best).
  - Recoating does work (in-house).

10X: Massively Parallel Sequencing
===

<div>
<img src="images/10X_paper.png" height=500>
</div>

10X: Description
===

- Droplet-based, 3' mRNA.
  - GEM (Gel Bead in Emulsion)
- Standardized instrumentation and reagents.
- More high-throughput scaling to tens of thousands.
- Less processing time.
- Up to 50% capture rate.
- Cell Ranger software is available for install.

10X: Assay Overview
===

![10x_assay](images/10X_wet.png)   
Zheng et al.

10X: Sequences
===

![10x_sequence](images/10X_wet_mol.png)

- 14nt barcode, 10nt UMI, 30nt polyT + sequencing adaptor and primers.

10X: Equipment
===

![10x_gem](images/10X_gem.png)

---

![10x_chip](images/dropseq_fluidics_big.png)

Comparing Assays
===

<div>
<img src="images/compare_power_paper.png" height=500>
</div>

Comparing Assays
===

Based on ERCC spike-ins.   

**Accuracy:** How well the abundance levels correlated with known spiked-in amounts.   

**Sensitivity:** Minimum number of input RNA molecules required to detect a spike-in.   

Comparing Assays
===

Svensson et al. 2016
<div>
<img src="images/compare_assays.png" height=500>
</div>

Comparing Assays
===

- Different assays have different throughput.
  - SmartSeq2 < Dropseq < 10X
- SmartSeq2 is full transcript.
- Plate-based methods get lysed in wells and so do not leak.
  - Droplet-based can have leaky RNA.
- In Dropseq assays RT happens outside the droplets
  - Can use harsher lysis buffers.
  - 10X needs lysis buffers compatible with the RT enzyme.
- 10X is more standardized and comes with a pipeline.
  - Dropseq is more customizable but more hands-on.
- Cost per library varies greatly.

A Word on Sorting
===

- After disassociating cells cells can be performed.
- Know your cells, are they sticky, are they big?
  - Select an appropriate sized nozzle.
- Don't sort too quickly (1-2k cells per second or lower)
  - the slower the more time cells sit in lysis after sorting
    - 10 minutes max in lysis (some say 30 minutes)
- Calibrate speed of instrument with beads
  - Check alignment every 5-6 plates
- Afterwards spin down to make sure cells are in lysis buffer
  - Flash freeze
- Chloe Villani on sorting https://youtu.be/HiKY0eHDhIs

Sequence Pipeline Overviews
===

![corgi_can](images/corgis/corgi_can.jpg)

Starting Sequences
===

<div align="center">
<img src="images/fastq.png" height=500>
</div>

Assays have Different Formats
===

![fastq_seq](images/sc_sequences.png)

Smartseq2: Pipeline Overview
===

![smart_seq_seq](images/population_rnaseq.png)

Dropseq: Pipeline Overview
===

![drop_seq_seq](images/dropseq_rnaseq.png)

---

- Common functionality: trimming, alignment, generating count matrix.
- Adds book keeping for cell barcodes and UMIs, bead error detection, cell barcode collapsing, UMI collapsing.

Dropseq: Further Help
===

![drop_seq_web](images/drop_seq_web.png)

10X: Pipeline Overview
===

![10_seq](images/10x_rnaseq.png)

---

- Steps very similar to the Dropseq pipeline.

10X: Further Help
===

![10_help](images/10x_web.png)

Sequence Level Quality Control
===

![seq_qc](images/population_rnaseq_qc.png)

---

- Much of the QC that is performed is using traditional tools.

Section Summary
===

- Single-cell RNA-Seq is a diverse ecosystem of assays.
  - Each assay has pros and cons.
- Sequences derived from these assays are complex and vary.
- Different pipelines are needed to address different sequence formats.
  - Common steps like aligning, QC, and reads counting are in common.

Let's Describe scRNA-Seq Count Data
===

![driving](images/corgis/driving_corgi.jpg)

We Start with a Matrix of Expression
===

![starting_matrix](images/starting_matrix.png)

Genes Have Different Distributions
===
class:small-code

![genes_log_normal](images/scdistributions_log_normal.png)

Genes Have Different Distributions
===
class:small-code

![genes_bimodal](images/scdistributions_bimodal.png)

Genes Have Different Distributions
===
class:small-code

![genes_zeroinfl](images/scdistributions_zeroinfl.png)

Genes Have Different Distributions
===
class:small-code

![genes_noise](images/scdistributions_noise.png)

Genes Have Different Distributions
===
class:small-code

![professor corgi](images/scdistributions_all.png)

Underlying Biology
===

Zero inflation.
- Drop-out event during reverse-transcription.
- Genes with more expression have less zeros.
- Complexity varies.

Transcription stochasticity.
- Transcription bursting.
- Coordinated transcription of multigene networks.
- Overdispersed counts.

Higher Resolution.
- More sources of signal
  - Cell cycle

---

![scde_bio](images/scde_bio.png)  
Karchenko et al.

Data analysis when using UMIs
===

**Read Counts**
 
![reads](images/reads.png)

---

 **Counts by UMIs**
 
![umi](images/umi.png)

Summary of the Data
===

- We are still understanding scData and how to apply it.
  - Data can be NOT normal.
  - Data can be Zero-inflated.
  - Data can be very noisy.
  - Cells vary in library complexity.
- Keeping these characteristics in analysis assumptions.
- Trend to filter more conservatively with UMIs.

Introduction to RStudio
===

![bed](images/corgis/dogbutt.jpg)

Motivation: Why Am I Using R?
===

- A lot of method development is happening in R.
- Free / open source / open science.
- Many supplemental computational biology packages.
- Data science is an art.
  - Data often requires one to create and manipulate analysis.
- This will allow you to experience key concepts in analysis.

RStudio (IDE)
===

![rstudio_overview](images/Rstudio_overview.png)

Initial data exploration
===

![professor corgi](images/corgis/professor_corgi.jpg)

Seurat
===

![suerat](images/seurat_pubs.png)

Seurat: Today's Data
===

- Using data provided by 10X (2700 cells).
- Peripheral blood mononuclear cells (PBMCs).
  - Any peripheral blood cell with a round nucleus.
  - Lymphocytes (T cells, B cells, NK cells)
  - Monocytes
- Complex mix cell populations used in immunity.

Logistics
===
class:small-code

- R allows methodology written by others to be imported.
  - Leverage other code.
  - Make your code available to others.  


```r
# Load libraries
library(dplyr) # Dataframe manipulation
library(Matrix) # Sparse matrices
library(useful) # Corner function
library(vioplot) # Violin pots
library(scater) # Single Cell QC
library(scde) # Differential Expression
library(org.Hs.eg.db) # Gene name manipulation
library(Seurat) # Single cell General Analysis
```

Representing Sparse Matrices
===


```r
# Load 10X data
pbmc.10X <- Read10X("./data/filtered_gene_bc_matrices/hg19")

# Memory use as a sparse matrix
object.size(pbmc.10X)
```

```
38715120 bytes
```

```r
# Memory use as a dense matrix
# 18 X more
object.size(as.matrix(pbmc.10X))
```

```
709264728 bytes
```

How are Sparse Representations Made?
===

- Sparse Matrix: A matrix where most of the elements are 0.
- Dense Matrix: A matrix where most elements are not 0.
- In this case, underlying data structure is a coordinate list.

2D Arrays vs Coordinate List
===

![sparse](images/sparseMatrices.png)

Create a Seurat Object
===


```r
# Expected raw counts (non-normalized data)
# Can give log transformed data but do not transform in setup method
pbmc.seurat <- new("seurat", raw.data=pbmc.10X)
```

What is in a Seurat Object?
===


```r
# Display the internal pieces of the Seurat Object
slotNames(pbmc.seurat)
```

```
 [1] "raw.data"            "data"                "scale.data"         
 [4] "var.genes"           "is.expr"             "ident"              
 [7] "pca.x"               "pca.rot"             "emp.pval"           
[10] "kmeans.obj"          "pca.obj"             "gene.scores"        
[13] "drop.coefs"          "wt.matrix"           "drop.wt.matrix"     
[16] "trusted.genes"       "drop.expr"           "data.info"          
[19] "project.name"        "kmeans.gene"         "kmeans.cell"        
[22] "jackStraw.empP"      "jackStraw.fakePC"    "jackStraw.empP.full"
[25] "pca.x.full"          "kmeans.col"          "mean.var"           
[28] "imputed"             "mix.probs"           "mix.param"          
[31] "final.prob"          "insitu.matrix"       "tsne.rot"           
[34] "ica.rot"             "ica.x"               "ica.obj"            
[37] "cell.names"          "cluster.tree"        "snn.sparse"         
[40] "snn.dense"           "snn.k"              
```

What is in a Seurat Object?
===


```
 [1] "raw.data"            "data"                "scale.data"         
 [4] "var.genes"           "is.expr"             "ident"              
 [7] "pca.x"               "pca.rot"             "emp.pval"           
[10] "kmeans.obj"          "pca.obj"             "gene.scores"        
[13] "drop.coefs"          "wt.matrix"           "drop.wt.matrix"     
[16] "trusted.genes"       "drop.expr"           "data.info"          
[19] "project.name"        "kmeans.gene"         "kmeans.cell"        
[22] "jackStraw.empP"      "jackStraw.fakePC"    "jackStraw.empP.full"
[25] "pca.x.full"          "kmeans.col"          "mean.var"           
[28] "imputed"             "mix.probs"           "mix.param"          
[31] "final.prob"          "insitu.matrix"       "tsne.rot"           
[34] "ica.rot"             "ica.x"               "ica.obj"            
[37] "cell.names"          "cluster.tree"        "snn.sparse"         
[40] "snn.dense"           "snn.k"              
```

What is in a Seurat Object?
===

Raw sparse matrix


```r
head(pbmc.seurat@raw.data)
```

```
6 x 2700 sparse Matrix of class "dgTMatrix"
                                                                          
MIR1302-10   . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . .
FAM138A      . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . .
OR4F5        . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . .
RP11-34P13.7 . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . .
RP11-34P13.8 . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . .
AL627309.1   . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . .
                                                                          
MIR1302-10   . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . .
FAM138A      . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . .
OR4F5        . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . .
RP11-34P13.7 . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . .
RP11-34P13.8 . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . .
AL627309.1   . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . .
                                                                          
MIR1302-10   . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . .
FAM138A      . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . .
OR4F5        . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . .
RP11-34P13.7 . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . .
RP11-34P13.8 . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . .
AL627309.1   . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . .
                                                                          
MIR1302-10   . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . .
FAM138A      . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . .
OR4F5        . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . .
RP11-34P13.7 . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . .
RP11-34P13.8 . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . .
AL627309.1   . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . .
                                                                          
MIR1302-10   . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . .
FAM138A      . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . .
OR4F5        . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . .
RP11-34P13.7 . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . .
RP11-34P13.8 . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . .
AL627309.1   . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . .
                                                                          
MIR1302-10   . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . .
FAM138A      . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . .
OR4F5        . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . .
RP11-34P13.7 . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . .
RP11-34P13.8 . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . .
AL627309.1   . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . .
                                                                          
MIR1302-10   . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . .
FAM138A      . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . .
OR4F5        . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . .
RP11-34P13.7 . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . .
RP11-34P13.8 . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . .
AL627309.1   . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . .
                                                                          
MIR1302-10   . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . .
FAM138A      . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . .
OR4F5        . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . .
RP11-34P13.7 . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . .
RP11-34P13.8 . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . .
AL627309.1   . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . .
                                                                          
MIR1302-10   . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . .
FAM138A      . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . .
OR4F5        . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . .
RP11-34P13.7 . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . .
RP11-34P13.8 . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . .
AL627309.1   . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . .
                                                                          
MIR1302-10   . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . .
FAM138A      . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . .
OR4F5        . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . .
RP11-34P13.7 . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . .
RP11-34P13.8 . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . .
AL627309.1   . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . .
                                                                          
MIR1302-10   . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . .
FAM138A      . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . .
OR4F5        . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . .
RP11-34P13.7 . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . .
RP11-34P13.8 . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . .
AL627309.1   . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . .
                                                                          
MIR1302-10   . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . .
FAM138A      . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . .
OR4F5        . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . .
RP11-34P13.7 . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . .
RP11-34P13.8 . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . .
AL627309.1   . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . .
                                                                          
MIR1302-10   . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . .
FAM138A      . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . .
OR4F5        . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . .
RP11-34P13.7 . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . .
RP11-34P13.8 . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . .
AL627309.1   . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . .
                                                                          
MIR1302-10   . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . .
FAM138A      . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . .
OR4F5        . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . .
RP11-34P13.7 . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . .
RP11-34P13.8 . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . .
AL627309.1   . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . .
                                                                          
MIR1302-10   . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . .
FAM138A      . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . .
OR4F5        . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . .
RP11-34P13.7 . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . .
RP11-34P13.8 . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . .
AL627309.1   . . . . . . . . . . 1 . . . . . . . . . . . . . . . . . . . .
                                                                          
MIR1302-10   . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . .
FAM138A      . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . .
OR4F5        . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . .
RP11-34P13.7 . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . .
RP11-34P13.8 . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . .
AL627309.1   . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . .
                                                                          
MIR1302-10   . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . .
FAM138A      . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . .
OR4F5        . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . .
RP11-34P13.7 . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . .
RP11-34P13.8 . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . .
AL627309.1   . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . .
                                                                          
MIR1302-10   . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . .
FAM138A      . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . .
OR4F5        . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . .
RP11-34P13.7 . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . .
RP11-34P13.8 . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . .
AL627309.1   . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . .
                                                                          
MIR1302-10   . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . .
FAM138A      . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . .
OR4F5        . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . .
RP11-34P13.7 . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . .
RP11-34P13.8 . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . .
AL627309.1   . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . .
                                                                          
MIR1302-10   . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . .
FAM138A      . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . .
OR4F5        . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . .
RP11-34P13.7 . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . .
RP11-34P13.8 . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . .
AL627309.1   . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . .
                                                                          
MIR1302-10   . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . .
FAM138A      . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . .
OR4F5        . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . .
RP11-34P13.7 . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . .
RP11-34P13.8 . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . .
AL627309.1   . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . .
                                                                          
MIR1302-10   . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . .
FAM138A      . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . .
OR4F5        . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . .
RP11-34P13.7 . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . .
RP11-34P13.8 . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . .
AL627309.1   . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . .
                                                                          
MIR1302-10   . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . .
FAM138A      . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . .
OR4F5        . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . .
RP11-34P13.7 . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . .
RP11-34P13.8 . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . .
AL627309.1   . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . .
                                                                          
MIR1302-10   . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . .
FAM138A      . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . .
OR4F5        . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . .
RP11-34P13.7 . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . .
RP11-34P13.8 . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . .
AL627309.1   . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . .
                                                                          
MIR1302-10   . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . .
FAM138A      . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . .
OR4F5        . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . .
RP11-34P13.7 . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . .
RP11-34P13.8 . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . .
AL627309.1   . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . .
                                                                          
MIR1302-10   . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . .
FAM138A      . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . .
OR4F5        . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . .
RP11-34P13.7 . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . .
RP11-34P13.8 . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . .
AL627309.1   . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . .
                                                                          
MIR1302-10   . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . .
FAM138A      . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . .
OR4F5        . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . .
RP11-34P13.7 . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . .
RP11-34P13.8 . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . .
AL627309.1   . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . .
                                                                          
MIR1302-10   . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . .
FAM138A      . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . .
OR4F5        . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . .
RP11-34P13.7 . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . .
RP11-34P13.8 . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . .
AL627309.1   . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . .
                                                                          
MIR1302-10   . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . .
FAM138A      . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . .
OR4F5        . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . .
RP11-34P13.7 . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . .
RP11-34P13.8 . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . .
AL627309.1   . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . .
                                                                          
MIR1302-10   . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . .
FAM138A      . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . .
OR4F5        . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . .
RP11-34P13.7 . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . .
RP11-34P13.8 . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . .
AL627309.1   . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . .
                                                                          
MIR1302-10   . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . .
FAM138A      . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . .
OR4F5        . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . .
RP11-34P13.7 . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . .
RP11-34P13.8 . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . .
AL627309.1   . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . .
                                                                          
MIR1302-10   . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . .
FAM138A      . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . .
OR4F5        . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . .
RP11-34P13.7 . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . .
RP11-34P13.8 . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . .
AL627309.1   . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . .
                                                                          
MIR1302-10   . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . .
FAM138A      . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . .
OR4F5        . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . .
RP11-34P13.7 . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . .
RP11-34P13.8 . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . .
AL627309.1   . . . . . . . . . . . . . . . . . . . . . . . . . 1 . . . . .
                                                                          
MIR1302-10   . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . .
FAM138A      . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . .
OR4F5        . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . .
RP11-34P13.7 . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . .
RP11-34P13.8 . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . .
AL627309.1   . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . .
                                                                          
MIR1302-10   . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . .
FAM138A      . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . .
OR4F5        . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . .
RP11-34P13.7 . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . .
RP11-34P13.8 . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . .
AL627309.1   . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . .
                                                                          
MIR1302-10   . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . .
FAM138A      . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . .
OR4F5        . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . .
RP11-34P13.7 . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . .
RP11-34P13.8 . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . .
AL627309.1   . . . . . . . . . . . . . . . . . . . . . 1 . . . . . . . . .
                                                                          
MIR1302-10   . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . .
FAM138A      . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . .
OR4F5        . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . .
RP11-34P13.7 . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . .
RP11-34P13.8 . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . .
AL627309.1   . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . .
                                                                          
MIR1302-10   . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . .
FAM138A      . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . .
OR4F5        . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . .
RP11-34P13.7 . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . .
RP11-34P13.8 . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . .
AL627309.1   . . . . . . . . . . . . . . 1 . . . . . . . . . . . . . . . .
                                                                          
MIR1302-10   . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . .
FAM138A      . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . .
OR4F5        . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . .
RP11-34P13.7 . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . .
RP11-34P13.8 . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . .
AL627309.1   . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . .
                                                                          
MIR1302-10   . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . .
FAM138A      . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . .
OR4F5        . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . .
RP11-34P13.7 . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . .
RP11-34P13.8 . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . .
AL627309.1   . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . .
                                                                          
MIR1302-10   . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . .
FAM138A      . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . .
OR4F5        . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . .
RP11-34P13.7 . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . .
RP11-34P13.8 . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . .
AL627309.1   . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . .
                                                                          
MIR1302-10   . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . .
FAM138A      . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . .
OR4F5        . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . .
RP11-34P13.7 . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . .
RP11-34P13.8 . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . .
AL627309.1   . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . .
                                                                          
MIR1302-10   . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . .
FAM138A      . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . .
OR4F5        . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . .
RP11-34P13.7 . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . .
RP11-34P13.8 . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . .
AL627309.1   . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . .
                                                                          
MIR1302-10   . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . .
FAM138A      . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . .
OR4F5        . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . .
RP11-34P13.7 . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . .
RP11-34P13.8 . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . .
AL627309.1   . . . . . . . . . . . . . . . . . . . . . . . . . . 1 . . . .
                                                                          
MIR1302-10   . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . .
FAM138A      . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . .
OR4F5        . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . .
RP11-34P13.7 . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . .
RP11-34P13.8 . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . .
AL627309.1   . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . .
                                                                          
MIR1302-10   . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . .
FAM138A      . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . .
OR4F5        . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . .
RP11-34P13.7 . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . .
RP11-34P13.8 . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . .
AL627309.1   . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . .
                                                                          
MIR1302-10   . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . .
FAM138A      . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . .
OR4F5        . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . .
RP11-34P13.7 . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . .
RP11-34P13.8 . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . .
AL627309.1   . . . . . . . . 1 . . . . . . . . . . . . . . . . . . . . . .
                                                                          
MIR1302-10   . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . .
FAM138A      . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . .
OR4F5        . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . .
RP11-34P13.7 . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . .
RP11-34P13.8 . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . .
AL627309.1   . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . .
                                                                          
MIR1302-10   . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . .
FAM138A      . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . .
OR4F5        . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . .
RP11-34P13.7 . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . .
RP11-34P13.8 . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . .
AL627309.1   . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . .
                                                                          
MIR1302-10   . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . .
FAM138A      . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . .
OR4F5        . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . .
RP11-34P13.7 . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . .
RP11-34P13.8 . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . .
AL627309.1   . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . .
                                                                          
MIR1302-10   . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . .
FAM138A      . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . .
OR4F5        . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . .
RP11-34P13.7 . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . .
RP11-34P13.8 . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . .
AL627309.1   . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . .
                                                                          
MIR1302-10   . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . .
FAM138A      . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . .
OR4F5        . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . .
RP11-34P13.7 . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . .
RP11-34P13.8 . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . .
AL627309.1   . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . .
                                                                          
MIR1302-10   . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . .
FAM138A      . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . .
OR4F5        . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . .
RP11-34P13.7 . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . .
RP11-34P13.8 . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . .
AL627309.1   . . . . . . . . . . . . . . . . 1 . . . . . . . . . . . . . .
                                                                          
MIR1302-10   . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . .
FAM138A      . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . .
OR4F5        . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . .
RP11-34P13.7 . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . .
RP11-34P13.8 . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . .
AL627309.1   . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . .
                                                                          
MIR1302-10   . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . .
FAM138A      . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . .
OR4F5        . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . .
RP11-34P13.7 . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . .
RP11-34P13.8 . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . .
AL627309.1   . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . .
                                                                          
MIR1302-10   . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . .
FAM138A      . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . .
OR4F5        . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . .
RP11-34P13.7 . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . .
RP11-34P13.8 . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . .
AL627309.1   . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . .
                                                                          
MIR1302-10   . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . .
FAM138A      . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . .
OR4F5        . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . .
RP11-34P13.7 . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . .
RP11-34P13.8 . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . .
AL627309.1   . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . .
                                                                          
MIR1302-10   . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . .
FAM138A      . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . .
OR4F5        . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . .
RP11-34P13.7 . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . .
RP11-34P13.8 . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . .
AL627309.1   . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . .
                                                                          
MIR1302-10   . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . .
FAM138A      . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . .
OR4F5        . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . .
RP11-34P13.7 . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . .
RP11-34P13.8 . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . .
AL627309.1   1 . . . . . . . . . . . . . 1 . . . . . . . . . . . . . . . .
                                                                          
MIR1302-10   . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . .
FAM138A      . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . .
OR4F5        . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . .
RP11-34P13.7 . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . .
RP11-34P13.8 . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . .
AL627309.1   . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . .
                                                                          
MIR1302-10   . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . .
FAM138A      . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . .
OR4F5        . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . .
RP11-34P13.7 . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . .
RP11-34P13.8 . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . .
AL627309.1   . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . .
                                                                          
MIR1302-10   . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . .
FAM138A      . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . .
OR4F5        . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . .
RP11-34P13.7 . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . .
RP11-34P13.8 . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . .
AL627309.1   . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . .
                                                                          
MIR1302-10   . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . .
FAM138A      . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . .
OR4F5        . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . .
RP11-34P13.7 . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . .
RP11-34P13.8 . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . .
AL627309.1   . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . .
                                                                          
MIR1302-10   . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . .
FAM138A      . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . .
OR4F5        . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . .
RP11-34P13.7 . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . .
RP11-34P13.8 . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . .
AL627309.1   . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . .
                                                                          
MIR1302-10   . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . .
FAM138A      . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . .
OR4F5        . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . .
RP11-34P13.7 . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . .
RP11-34P13.8 . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . .
AL627309.1   . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . .
                                                                          
MIR1302-10   . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . .
FAM138A      . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . .
OR4F5        . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . .
RP11-34P13.7 . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . .
RP11-34P13.8 . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . .
AL627309.1   . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . .
                                                                          
MIR1302-10   . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . .
FAM138A      . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . .
OR4F5        . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . .
RP11-34P13.7 . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . .
RP11-34P13.8 . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . .
AL627309.1   . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . .
                                                                          
MIR1302-10   . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . .
FAM138A      . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . .
OR4F5        . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . .
RP11-34P13.7 . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . .
RP11-34P13.8 . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . .
AL627309.1   . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . .
                                                                          
MIR1302-10   . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . .
FAM138A      . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . .
OR4F5        . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . .
RP11-34P13.7 . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . .
RP11-34P13.8 . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . .
AL627309.1   . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . .
                                                                          
MIR1302-10   . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . .
FAM138A      . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . .
OR4F5        . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . .
RP11-34P13.7 . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . .
RP11-34P13.8 . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . .
AL627309.1   . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . .
                                                                          
MIR1302-10   . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . .
FAM138A      . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . .
OR4F5        . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . .
RP11-34P13.7 . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . .
RP11-34P13.8 . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . .
AL627309.1   . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . .
                                                                          
MIR1302-10   . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . .
FAM138A      . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . .
OR4F5        . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . .
RP11-34P13.7 . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . .
RP11-34P13.8 . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . .
AL627309.1   . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . .
                                                                          
MIR1302-10   . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . .
FAM138A      . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . .
OR4F5        . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . .
RP11-34P13.7 . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . .
RP11-34P13.8 . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . .
AL627309.1   . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . .
                                                                          
MIR1302-10   . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . .
FAM138A      . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . .
OR4F5        . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . .
RP11-34P13.7 . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . .
RP11-34P13.8 . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . .
AL627309.1   . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . .
                                                                          
MIR1302-10   . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . .
FAM138A      . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . .
OR4F5        . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . .
RP11-34P13.7 . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . .
RP11-34P13.8 . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . .
AL627309.1   . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . .
                                                                          
MIR1302-10   . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . .
FAM138A      . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . .
OR4F5        . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . .
RP11-34P13.7 . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . .
RP11-34P13.8 . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . .
AL627309.1   . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . .
                                                                          
MIR1302-10   . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . .
FAM138A      . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . .
OR4F5        . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . .
RP11-34P13.7 . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . .
RP11-34P13.8 . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . .
AL627309.1   . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . .
                                                                          
MIR1302-10   . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . .
FAM138A      . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . .
OR4F5        . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . .
RP11-34P13.7 . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . .
RP11-34P13.8 . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . .
AL627309.1   . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . .
                                                                          
MIR1302-10   . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . .
FAM138A      . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . .
OR4F5        . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . .
RP11-34P13.7 . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . .
RP11-34P13.8 . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . .
AL627309.1   . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . .
                                                                          
MIR1302-10   . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . .
FAM138A      . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . .
OR4F5        . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . .
RP11-34P13.7 . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . .
RP11-34P13.8 . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . .
AL627309.1   . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . .
                                                                          
MIR1302-10   . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . .
FAM138A      . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . .
OR4F5        . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . .
RP11-34P13.7 . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . .
RP11-34P13.8 . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . .
AL627309.1   . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . .
                                                                          
MIR1302-10   . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . .
FAM138A      . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . .
OR4F5        . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . .
RP11-34P13.7 . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . .
RP11-34P13.8 . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . .
AL627309.1   . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . .
                                                                          
MIR1302-10   . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . .
FAM138A      . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . .
OR4F5        . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . .
RP11-34P13.7 . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . .
RP11-34P13.8 . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . .
AL627309.1   . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . .
                                                                          
MIR1302-10   . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . .
FAM138A      . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . .
OR4F5        . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . .
RP11-34P13.7 . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . .
RP11-34P13.8 . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . .
AL627309.1   . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . .
                                                                          
MIR1302-10   . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . .
FAM138A      . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . .
OR4F5        . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . .
RP11-34P13.7 . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . .
RP11-34P13.8 . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . .
AL627309.1   . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . .
                                                                          
MIR1302-10   . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . .
FAM138A      . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . .
OR4F5        . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . .
RP11-34P13.7 . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . .
RP11-34P13.8 . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . .
AL627309.1   . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . .
                                                                          
MIR1302-10   . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . .
FAM138A      . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . .
OR4F5        . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . .
RP11-34P13.7 . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . .
RP11-34P13.8 . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . .
AL627309.1   . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . .
                  
MIR1302-10   . . .
FAM138A      . . .
OR4F5        . . .
RP11-34P13.7 . . .
RP11-34P13.8 . . .
AL627309.1   . . .
```

What is in a Seurat Object?
===

- __var.genes__ Variable genes across cells
- __data.info__ Misc info including complexity (nGene) 
- __cell.names__ Column (cell) names 
- __gene.names__ Row (gene) names


```r
?seurat
```

Always Look at Your Data
===
class:midcenter

_**Hiding within those mounds of data is knowledge that could change the life of a patient, or change the world.**_  -– Atul Butte

What are Our Genes?
===
class:small-code


```r
# Gene names (row names)
head(row.names(pbmc.seurat@raw.data))
```

```
[1] "MIR1302-10"   "FAM138A"      "OR4F5"        "RP11-34P13.7"
[5] "RP11-34P13.8" "AL627309.1"  
```

```r
length(row.names(pbmc.seurat@raw.data))
```

```
[1] 32738
```

What are Our Cells?
===
class:small-code


```r
# Column names
# Sample / Cell names / Barcodes
head(colnames(pbmc.seurat@raw.data))
```

```
[1] "AAACATACAACCAC" "AAACATTGAGCTAC" "AAACATTGATCAGC" "AAACCGTGCTTCCG"
[5] "AAACCGTGTATGCG" "AAACGCACTGGTAC"
```

```r
length(colnames(pbmc.seurat@raw.data))
```

```
[1] 2700
```

How to Show Counts?
===
class:small-code


```r
# Only the corner
# The full data will be too large to see
corner(as.matrix(pbmc.seurat@raw.data))
```

How to Show Counts?
===
class:small-code


```
             AAACATACAACCAC AAACATTGAGCTAC AAACATTGATCAGC AAACCGTGCTTCCG
MIR1302-10                0              0              0              0
FAM138A                   0              0              0              0
OR4F5                     0              0              0              0
RP11-34P13.7              0              0              0              0
RP11-34P13.8              0              0              0              0
             AAACCGTGTATGCG
MIR1302-10                0
FAM138A                   0
OR4F5                     0
RP11-34P13.7              0
RP11-34P13.8              0
```

Complexity of Libraries
===

- Very important metadata.
- How many genes had a count of more than 1 in a cell.
  - Simplest, most common definition.

How Many Expressed Genes?
===
class:small-code


```r
# Plot genes per cell How many genes expressed per cells
complexity.per.cell <- apply(pbmc.seurat@raw.data, 2, function(x) sum(x > 0))
# Mean count per cell.
mean.count.per.cell <- apply(pbmc.seurat@raw.data, 2, function(x) mean(x))
# Gene prevalence
gene.prevalence <- apply(pbmc.seurat@raw.data, 1, function(x) sum(x > 0))
```

How Many Expressed Genes?
===
class:small-code


```r
# Plot genes per cell How many genes expressed per cell
vioplot(complexity.per.cell)
stripchart(complexity.per.cell, add = TRUE, vertical = TRUE, method = "jitter", 
    jitter = 0.3, pch = 19)
abline(h = 200, col = "red")
abline(h = 2500, col = "blue")
title("Study Complexity")
axis(side = 1, at = 1, labels = c("Study"))
```

How Many Expressed Genes?
===
class:small-code

![plot of chunk unnamed-chunk-14](single_cell_analysis-figure/unnamed-chunk-14-1.png)

---

- Shown complexity of all cells in a study. 
- Wide range of complexity.
- Prefiltered data at red line.
- Will define points about the blue line as outlier.

Box Plots vs Violin Plots
===
class:small-code

**Boxplot:**  
Robust representation of a distribution using quantiles

![plot of chunk unnamed-chunk-15](single_cell_analysis-figure/unnamed-chunk-15-1.png)

---

**Violin Plot:**  
Box plot with kernel density plot mirrored on sides.

![plot of chunk unnamed-chunk-16](single_cell_analysis-figure/unnamed-chunk-16-1.png)

Identifying Outliers?
===
class:small-code

- Cells that are unusually simple (or no counts)
- Cells that are unusually complex (doublets?)


```r
plot(complexity.per.cell, mean.count.per.cell + 1)
abline(v = 200, col = "red")
abline(h = log2(4))
```

Identifying Outliers?
===
class:small-code

- We will filter in a standard way, we are just describing the data.
- Red line shows a point of prefiltering.

![plot of chunk unnamed-chunk-18](single_cell_analysis-figure/unnamed-chunk-18-1.png)

Identifying Noise?
===
class:small-code


```r
# hist(gene.prevalence)
hist(log2(gene.prevalence))
abline(v = 3, col = "red")
```

Identifying Noise?
===
class:small-code

- People tend to filter very conservatively.
- 3 is where we will filter (red line).

![plot of chunk unnamed-chunk-20](single_cell_analysis-figure/unnamed-chunk-20-1.png)

Filter Cells: Removing the Outlier Cells
===
class:small-code

- Genes must be in 3 cells with alteast 200 genes.
- Scaled by 1000 (Total Sum Scaling)

```r
pbmc.seurat <- Setup(pbmc.seurat, min.cells = 3, min.genes = 200, do.logNormalize = TRUE, 
    total.expr = 10000, project = "Tutorial")
```

```
[1] "Performing log-normalization"
  |                                                                         |                                                                 |   0%  |                                                                         |======================                                           |  33%  |                                                                         |===========================================                      |  67%  |                                                                         |=================================================================| 100%
[1] "Scaling data matrix"
  |                                                                         |                                                                 |   0%  |                                                                         |=====                                                            |   7%  |                                                                         |=========                                                        |  14%  |                                                                         |==============                                                   |  21%  |                                                                         |===================                                              |  29%  |                                                                         |=======================                                          |  36%  |                                                                         |============================                                     |  43%  |                                                                         |================================                                 |  50%  |                                                                         |=====================================                            |  57%  |                                                                         |==========================================                       |  64%  |                                                                         |==============================================                   |  71%  |                                                                         |===================================================              |  79%  |                                                                         |========================================================         |  86%  |                                                                         |============================================================     |  93%  |                                                                         |=================================================================| 100%
```

Expression: Definition
===

**Population based RNA-Seq**
- Different ways to move from counts to expression.
- EM algorithms (RSEM)
- TPM with TMM normalization.
  - Control for transcript length and force the depth to be 10^6.
  - Trimmed means weighted focusing on invariant genes through samples.

**Today we are using**
- log2( ( cell_counts / cell_sum ) * 1000 + 1 )
- Account for sequencing depth.

Filter Genes: Averages are Less Useful
===
class:small-code

_**Say you were standing with one foot in the oven and one foot in an ice bucket. According to the percentage people, you should be perfectly comfortable.**_ --Bobby Bragan

Filter Genes: Averages are Less Useful
===
class:small-code

![average](images/prob_ave_expression.png)

Filter Genes: Using Prevalence
===

**Filtering with gene prevalence:**   
How many times a gene's count is greater than or equal to an expression threshold throughout cells.

![prevalence](images/prevalence.png)

Filter Genes: Using Prevalence
===
class:small-code

![occurence](images/prob_ave_expression2.png)

Filtering Using Metadata
===

![puppies](images/corgis/puppies.jpg)

What is Metadata?
===

Other information that describes your measurements.
- Patient information.
  - Life style (smoking), Patient Biology (age), Comorbidity
- Study information.
  - Treatment, Cage, Sequencing Site, Sequencing Date
- Sequence QC on cells.
  - Useful in filtering.

Filter Cells: Removing the Outlier Cells
===
class:small-code

- Outlier cells are not just measured by complexity
  - **Useful Metrics**
     - Percent Reads Mapping
     - Percent Mitochondrial Reads
     - Presence of marker genes
     - Intergenic/ exonic rate
     - 5' or 3' bias
     - other metadata ...
  - **Useful Tools**
     - Picard Tools and RNASeQC

Seurat: Filtering on Metadata
===
class:small-code

- Filtering on mitochondrial reads in Seurat.  

```r
# Get gene names
mito.gene.names <- grep("^MT-", rownames(pbmc.seurat@data), value=TRUE)

# Get TSS normalized mitochodrial counts
col.total.counts <- Matrix::colSums(expm1(pbmc.seurat@data))
mito.percent.counts <- Matrix::colSums(expm1(pbmc.seurat@data[mito.gene.names, ]))/col.total.counts

# Add to seurat object as a metadata
pbmc.seurat <- AddMetaData(pbmc.seurat, mito.percent.counts, "percent.mitochodrial")
```

Seurat: Filtering on Metadata
===
class:small-code

- Plot current metadata in Seurat Object.
 - Number gene.
 - Number UMI.
 - Percent mitochondrial counts.
 


Seurat: Filtering on Metadata
===
class:small-code

![plot of chunk unnamed-chunk-24](single_cell_analysis-figure/unnamed-chunk-24-1.png)

Seurat: Filtering on Metadata
===
class:small-code


```r
GenePlot(pbmc.seurat, "nUMI", "percent.mitochodrial")
```

Seurat: Filtering on Metadata
===
class:small-code

- Outlier percent mitochondria are very low expression.
- Very high expression low percent mitochondrial reads.

![plot of chunk unnamed-chunk-26](single_cell_analysis-figure/unnamed-chunk-26-1.png)

Seurat: Filtering on Metadata
===
class:small-code

![plot of chunk unnamed-chunk-27](single_cell_analysis-figure/unnamed-chunk-27-1.png)

Seurat: Filtering on Metadata
===
class:small-code


```r
dim(pbmc.seurat@data)
```

```
[1] 13714  2700
```

```r
pbmc.seurat <- SubsetData(pbmc.seurat, subset.name = "nGene", accept.high = 2500)
pbmc.seurat <- SubsetData(pbmc.seurat, subset.name = "percent.mitochodrial", accept.high = 0.05)
dim(pbmc.seurat@data)
```

```
[1] 13714  2638
```

Saving as an R Object
===
class:small-code

Saving the Seurat object
- Contains all manipulation so far.
- Can be loaded or shared.
- Does not contain environment.


```r
# How to save the intact object.
save(pbmc.seurat, file = "seurat_tutorial.Robj")
```

```r
# How to retrieve the intact object.
load("seurat_tutorial.Robj")
```

Saving as Text Files
===

You may need to export data to import into other applications.


```r
# Log-scale expression matrix
write.table(as.matrix(pbmc.seurat@data), file = "seurat_data.txt")

# Study metadata
write.table(pbmc.seurat@data.info, file = "seurat_metadata.txt")

# What is the metadata so far
head(pbmc.seurat@data.info)
```

Scater: Complex QC Plotting
===

![scater_paper](images/scater_paper.png)

Scater: Motivation
===

- Clinical studies are often rich with metadata.
- Most differential expression methodology is T-test.
  - Much like using covariates in linear regression to control here we stratify by metadata to describe the impact of the metadatum.
- More than what we are covering here.
- Online tutorial [Click Here](//bioconductor.org/packages/devel/bioc/vignettes/scater/inst/doc/vignette.html)

Scater: Skipping Ahead
===
class:small-code


```r
# Load Data
data("sc_example_counts")
data("sc_example_cell_info")
pd <- new("AnnotatedDataFrame", data=sc_example_cell_info)
rownames(pd) <- pd$Cell
example_sceset <- newSCESet(countData=sc_example_counts, phenoData=pd)
keep_feature <- rowSums(exprs(example_sceset)) > 0
example_sceset <- example_sceset[keep_feature,]
example_sceset <- calculateQCMetrics(example_sceset, feature_controls = 1:40)
```

Scater: Describe the Data
===

Small data set of 40 cells.
- Constructed metadata
- Metadata including
  - Mutation status (positive, negative)
  - Cell Cycle (S, G0, G1, G2M)
  - Treament (treat1, treat2)

Scater: What is the Data?
===
class:small-code


```r
corner(sc_example_counts)
```

```
          Cell_001 Cell_002 Cell_003 Cell_004 Cell_005
Gene_0001        0      123        2        0        0
Gene_0002      575       65        3     1561     2311
Gene_0003        0        0        0        0     1213
Gene_0004        0        1        0        0        0
Gene_0005        0        0       11        0        0
```

Scater: What is the info?
===
class:small-code


```r
corner(sc_example_cell_info)
```

```
             Cell Mutation_Status Cell_Cycle Treatment
Cell_001 Cell_001        positive          S    treat1
Cell_002 Cell_002        positive         G0    treat1
Cell_003 Cell_003        negative         G1    treat1
Cell_004 Cell_004        negative          S    treat1
Cell_005 Cell_005        negative         G1    treat2
```

Scater: Stratifying for Metadata
===
class:small-code


```r
plot(example_sceset, block1 = "Mutation_Status", block2 = "Treatment",
     colour_by = "Cell_Cycle", nfeatures = 300, exprs_values = "counts")
```

Scater: Stratifying for Metadata
===
class:small-code

![plot of chunk unnamed-chunk-36](single_cell_analysis-figure/unnamed-chunk-36-1.png)

Scater: Stratifying for Metadata
===
class:small-code


```r
plotExpression(example_sceset, rownames(example_sceset)[1:6],
               x = "Mutation_Status", exprs_values = "exprs", colour = "Treatment")
```

Scater: Stratifying for Metadata
===
class:small-code

![plot of chunk unnamed-chunk-38](single_cell_analysis-figure/unnamed-chunk-38-1.png)

Scater: Interactive Exploration
===
class:small-code


```r
scater_gui(example_sceset)
```

Section Summary
===

- scRNA-Seq is complex and can NOT be treated as population based measurements.
- Cell complexity is an important metadatum for QC.
- Conservative gene filtering based on gene count and prevalence is often used.
- Cells can be filtered by complexity and also other metadatum.

Where are we at?
===

We have now
- Loaded in data into R.
- Explored and plotted the data for QC.
- Filtered on counts and metadata.

Next
- Focus on plotting Genes and Cells.
- View cell groups to later use in DE.

Let's Get Started!
===

<div align="center">
<img src="images/corgis/pillows.jpg" height=500>
</div>

Welcome to section 2!
===

- Using Seurat to plot genes
- Dimensionality reduction and plotting (PCA, t-SNE).
- Plotting cells.
- Batch effects and study confounding.
- Moving from clusters to populations of cells.
- Differential Expression (SCDE).
- Pathway analysis.
- Overview of available methodology.

Using Seurat to plot genes
===

![tea](images/corgis/tea.jpg)

Seurat: Viewing Specific Genes
===
class:small-code

- Check the identity of the cells!!!
- Notice many zeros.


```r
VlnPlot(pbmc.seurat, c("GAPDH"))
```

Seurat: Viewing Specific Genes
===
class:small-code

- Check the identity of the cells!!!
- Notice many zeros.


```r
VlnPlot(pbmc.seurat, c("GAPDH"))
```

Seurat: Plotting Genes vs Genes
===
class:small-code


```r
# Plot a gene vs a gene
GenePlot(pbmc.seurat, "CD79A", "CD79B", cex.use = 1)
```

Seurat: Plotting Genes vs Genes
===
class:small-code


```r
GenePlot(pbmc.seurat, "CD79A", "CD79B", cex.use = 1)
```

<img src="single_cell_analysis-figure/unnamed-chunk-43-1.png" title="plot of chunk unnamed-chunk-43" alt="plot of chunk unnamed-chunk-43" style="display: block; margin: auto;" />

Dimensionality reduction and plotting
===

![outlier_corgi](images/corgis/outlier_corgi.jpg)

Dimensionality Reduction and Ordination
===

- Start with many measurements (high dimensional).
  - Want to reduce to a few features (lower-dimensional space).
- One way is to extract features based on capturing groups of variance.
- Another could be to preferentially select some of the current features.
  - We have already done this.
- We need this to plot the cells in 2D (or ordinate them)
- In scRNA-Seq PC1 may be complexity.

PCA: in Quick Theory
===

- Eigenvectors of covariance matrix.
- Find orthogonal groups of variance.
- Given from most to least variance.
  - Components of variation.
  - Linear combinations explaining the variance.
  
---

![pca_describe](images/PCA_plot_02.gif)

PCA: in Quick Theory
===

http://setosa.io/ev/principal-component-analysis

PCA: in Practice
===

Things to be aware of.
- Data with different magnitudes will dominate.
  - Zero center and divided by SD.
  - (Standardized).
- Can be affected by outliers.

t-SNE: Nonlinear Dimensional Reduction
===

![tsne](images/cell_cycle_1.png)

t-SNE: Collapsing the Visualization to 2D
===

![tsne_collapsed](images/cell_cycle_collapse.png)

t-SNE: What and Why?
===

![tsne](images/tsne.png)

PCA and t-SNE Together
===

Often t-SNE is performed on PCA components
- Liberal number of components.
- Removes mild signal (assumption of noise).
- Faster, on less data but, hopefully the same signal.

Plotting Cells
===

![fancy](images/corgis/fancy.jpg)

Seurat: Performing PCA
===
class:small-code

- Calculating PCA with the highly variable genes.   


```r
# Select highly variable genese
pbmc.seurat <- MeanVarPlot(pbmc.seurat,fxn.x=expMean,fxn.y=logVarDivMean,
                           x.low.cutoff=0.0125,x.high.cutoff=3,
                           y.cutoff=0.5,do.contour=FALSE,do.plot=FALSE)
```

```
[1] "Calculating gene dispersion"
  |                                                                         |                                                                 |   0%  |                                                                         |=====                                                            |   7%  |                                                                         |=========                                                        |  14%  |                                                                         |==============                                                   |  21%  |                                                                         |===================                                              |  29%  |                                                                         |=======================                                          |  36%  |                                                                         |============================                                     |  43%  |                                                                         |================================                                 |  50%  |                                                                         |=====================================                            |  57%  |                                                                         |==========================================                       |  64%  |                                                                         |==============================================                   |  71%  |                                                                         |===================================================              |  79%  |                                                                         |========================================================         |  86%  |                                                                         |============================================================     |  93%  |                                                                         |=================================================================| 100%
```

```r
pbmc.seurat <- PCA(pbmc.seurat,pc.genes=pbmc.seurat@var.genes,do.print=FALSE)
```

Seurat: Performing PCA
===
class:small-code


```r
# Calculate PCA projection
pbmc.seurat <- ProjectPCA(pbmc.seurat)
```

```
  |                                                                         |                                                                 |   0%  |                                                                         |=====                                                            |   7%  |                                                                         |=========                                                        |  14%  |                                                                         |==============                                                   |  21%  |                                                                         |===================                                              |  29%  |                                                                         |=======================                                          |  36%  |                                                                         |============================                                     |  43%  |                                                                         |================================                                 |  50%  |                                                                         |=====================================                            |  57%  |                                                                         |==========================================                       |  64%  |                                                                         |==============================================                   |  71%  |                                                                         |===================================================              |  79%  |                                                                         |========================================================         |  86%  |                                                                         |============================================================     |  93%  |                                                                         |=================================================================| 100%
[1] "PC1"
 [1] "MALAT1"  "RPS27A"  "PTPRCAP" "RPS27"   "RPL23A"  "RPS3A"   "RPL3"   
 [8] "IL32"    "LTB"     "CD3D"    "RPL21"   "RPSA"    "RPL9"    "RPL13A" 
[15] "RPS3"    "RPS6"    "RPL31"   "CD3E"    "RPL30"   "RPS15A"  "RPS25"  
[22] "RPS12"   "RPS18"   "LDHB"    "RPS23"   "RPS29"   "RPL27A"  "RPL13"  
[29] "RPLP2"   "CXCR4"  
[1] ""
 [1] "CST3"     "TYROBP"   "LST1"     "AIF1"     "FCN1"     "LYZ"     
 [7] "S100A9"   "FTL"      "FTH1"     "TYMP"     "FCER1G"   "CFD"     
[13] "LGALS1"   "LGALS2"   "S100A8"   "CD68"     "CTSS"     "SERPINA1"
[19] "IFITM3"   "SPI1"     "SAT1"     "IFI30"    "COTL1"    "PSAP"    
[25] "CFP"      "NPC2"     "GRN"      "S100A11"  "LGALS3"   "AP1S2"   
[1] ""
[1] ""
[1] "PC2"
 [1] "CD79A"     "MS4A1"     "TCL1A"     "HLA-DQA1"  "HLA-DQB1" 
 [6] "LINC00926" "RPL18A"    "CD79B"     "HLA-DRA"   "VPREB3"   
[11] "RPL32"     "LTB"       "RPL13A"    "RPL13"     "FCER2"    
[16] "HLA-DQA2"  "RPL11"     "CD74"      "HLA-DRB1"  "BANK1"    
[21] "CD37"      "RPS23"     "RPS27"     "RPL8"      "HLA-DPB1" 
[26] "HLA-DMA"   "RPL12"     "RPS18"     "TSPAN13"   "FCRLA"    
[1] ""
 [1] "NKG7"   "GZMB"   "PRF1"   "CST7"   "GZMA"   "FGFBP2" "GNLY"  
 [8] "CTSW"   "SPON2"  "GZMH"   "CCL4"   "B2M"    "FCGR3A" "CCL5"  
[15] "XCL2"   "KLRD1"  "CD247"  "CLIC3"  "GZMM"   "AKR1C3" "SRGN"  
[22] "TTC38"  "HLA-C"  "HCST"   "PRSS23" "HOPX"   "IGFBP7" "S1PR5" 
[29] "ITGB2"  "GPR56" 
[1] ""
[1] ""
[1] "PC3"
 [1] "PPBP"       "PF4"        "SDPR"       "GNG11"      "SPARC"     
 [6] "HIST1H2AC"  "NRGN"       "GP9"        "RGS18"      "TUBB1"     
[11] "CLU"        "AP001189.4" "CD9"        "ITGA2B"     "PTCRA"     
[16] "TMEM40"     "CA2"        "ACRBP"      "MMD"        "TREML1"    
[21] "F13A1"      "PGRMC1"     "SEPT5"      "MYL9"       "TSC22D1"   
[26] "MPP1"       "CMTM5"      "PTGS1"      "SNCA"       "RUFY1"     
[1] ""
 [1] "RPL10"  "RPS2"   "RPL11"  "RPL18A" "RPL32"  "RPL28"  "RPS19" 
 [8] "RPL12"  "RPL19"  "RPL13"  "RPS6"   "RPS14"  "RPS15"  "TMSB10"
[15] "RPLP1"  "RPL29"  "RPL6"   "RPL26"  "RPS4X"  "RPS16"  "RPS3"  
[22] "EEF1A1" "RPL15"  "RPL13A" "RPS7"   "RPL8"   "RPS12"  "RPL23A"
[29] "RPS18"  "RPLP2" 
[1] ""
[1] ""
[1] "PC4"
 [1] "CD79A"     "HLA-DQA1"  "CD79B"     "MS4A1"     "HLA-DQB1" 
 [6] "CD74"      "HLA-DPB1"  "HLA-DPA1"  "HLA-DRB1"  "TCL1A"    
[11] "HLA-DRA"   "LINC00926" "HLA-DQA2"  "HLA-DRB5"  "VPREB3"   
[16] "HLA-DMA"   "HLA-DMB"   "FCER2"     "BANK1"     "HVCN1"    
[21] "GZMB"      "HLA-DOB"   "PDLIM1"    "FCRLA"     "TSPAN13"  
[26] "FGFBP2"    "CD72"      "EAF2"      "PKIG"      "SPIB"     
[1] ""
 [1] "CD3D"    "LDHB"    "RPS14"   "IL7R"    "CD3E"    "RPL32"   "VIM"    
 [8] "IL32"    "RPS12"   "NOSIP"   "RPL28"   "GIMAP7"  "RPL11"   "RPL13"  
[15] "FYB"     "ZFP36L2" "RPL10"   "AQP3"    "JUNB"    "RPS25"   "RPLP1"  
[22] "MAL"     "LEF1"    "RGCC"    "S100A6"  "FOS"     "RPS3"    "CD2"    
[29] "RPL35A"  "RPL36"  
[1] ""
[1] ""
[1] "PC5"
 [1] "LTB"      "TMEM66"   "LDHB"     "HSPA8"    "PABPC1"   "NPM1"    
 [7] "RNASET2"  "JUNB"     "SOD1"     "CD52"     "NAP1L1"   "CALM2"   
[13] "EEF1A1"   "RPSA"     "HINT1"    "IL7R"     "VIM"      "FXYD5"   
[19] "AQP3"     "ARHGDIB"  "LDHA"     "HNRNPA1"  "NACA"     "HSP90AA1"
[25] "GSTK1"    "ITM2B"    "RPLP0"    "RPS19"    "ACTG1"    "RPS10"   
[1] ""
 [1] "S100A8"  "FGFBP2"  "GZMB"    "S100A9"  "NKG7"    "CCL4"    "S100A12"
 [8] "LGALS2"  "RBP7"    "GNLY"    "CST7"    "SPON2"   "GZMA"    "PRF1"   
[15] "CCL3"    "FOLR3"   "MS4A6A"  "GZMH"    "CD14"    "PRSS23"  "S1PR5"  
[22] "KLRD1"   "TYROBP"  "XCL2"    "CTSW"    "CCL5"    "CLIC3"   "CEBPD"  
[29] "TTC38"   "GSTP1"  
[1] ""
[1] ""
```

```r
# Can plot top genes for top components
PrintPCA(pbmc.seurat, pcs.print = 1:2, genes.print = 5, use.full = TRUE)
```

```
[1] "PC1"
[1] "MALAT1"  "RPS27A"  "PTPRCAP" "RPS27"   "RPL23A" 
[1] ""
[1] "CST3"   "TYROBP" "LST1"   "AIF1"   "FCN1"  
[1] ""
[1] ""
[1] "PC2"
[1] "CD79A"    "MS4A1"    "TCL1A"    "HLA-DQA1" "HLA-DQB1"
[1] ""
[1] "NKG7" "GZMB" "PRF1" "CST7" "GZMA"
[1] ""
[1] ""
```

Seurat: Performing PCA
===
class:small-code


```r
# Calculate PCA projection
pbmc.seurat <- ProjectPCA(pbmc.seurat)
```

```
  |                                                                         |                                                                 |   0%  |                                                                         |=====                                                            |   7%  |                                                                         |=========                                                        |  14%  |                                                                         |==============                                                   |  21%  |                                                                         |===================                                              |  29%  |                                                                         |=======================                                          |  36%  |                                                                         |============================                                     |  43%  |                                                                         |================================                                 |  50%  |                                                                         |=====================================                            |  57%  |                                                                         |==========================================                       |  64%  |                                                                         |==============================================                   |  71%  |                                                                         |===================================================              |  79%  |                                                                         |========================================================         |  86%  |                                                                         |============================================================     |  93%  |                                                                         |=================================================================| 100%
[1] "PC1"
 [1] "MALAT1"  "RPS27A"  "PTPRCAP" "RPS27"   "RPL23A"  "RPS3A"   "RPL3"   
 [8] "IL32"    "LTB"     "CD3D"    "RPL21"   "RPSA"    "RPL9"    "RPL13A" 
[15] "RPS3"    "RPS6"    "RPL31"   "CD3E"    "RPL30"   "RPS15A"  "RPS25"  
[22] "RPS12"   "RPS18"   "LDHB"    "RPS23"   "RPS29"   "RPL27A"  "RPL13"  
[29] "RPLP2"   "CXCR4"  
[1] ""
 [1] "CST3"     "TYROBP"   "LST1"     "AIF1"     "FCN1"     "LYZ"     
 [7] "S100A9"   "FTL"      "FTH1"     "TYMP"     "FCER1G"   "CFD"     
[13] "LGALS1"   "LGALS2"   "S100A8"   "CD68"     "CTSS"     "SERPINA1"
[19] "IFITM3"   "SPI1"     "SAT1"     "IFI30"    "COTL1"    "PSAP"    
[25] "CFP"      "NPC2"     "GRN"      "S100A11"  "LGALS3"   "AP1S2"   
[1] ""
[1] ""
[1] "PC2"
 [1] "CD79A"     "MS4A1"     "TCL1A"     "HLA-DQA1"  "HLA-DQB1" 
 [6] "LINC00926" "RPL18A"    "CD79B"     "HLA-DRA"   "VPREB3"   
[11] "RPL32"     "LTB"       "RPL13A"    "RPL13"     "FCER2"    
[16] "HLA-DQA2"  "RPL11"     "CD74"      "HLA-DRB1"  "BANK1"    
[21] "CD37"      "RPS23"     "RPS27"     "RPL8"      "HLA-DPB1" 
[26] "HLA-DMA"   "RPL12"     "RPS18"     "TSPAN13"   "FCRLA"    
[1] ""
 [1] "NKG7"   "GZMB"   "PRF1"   "CST7"   "GZMA"   "FGFBP2" "GNLY"  
 [8] "CTSW"   "SPON2"  "GZMH"   "CCL4"   "B2M"    "FCGR3A" "CCL5"  
[15] "XCL2"   "KLRD1"  "CD247"  "CLIC3"  "GZMM"   "AKR1C3" "SRGN"  
[22] "TTC38"  "HLA-C"  "HCST"   "PRSS23" "HOPX"   "IGFBP7" "S1PR5" 
[29] "ITGB2"  "GPR56" 
[1] ""
[1] ""
[1] "PC3"
 [1] "PPBP"       "PF4"        "SDPR"       "GNG11"      "SPARC"     
 [6] "HIST1H2AC"  "NRGN"       "GP9"        "RGS18"      "TUBB1"     
[11] "CLU"        "AP001189.4" "CD9"        "ITGA2B"     "PTCRA"     
[16] "TMEM40"     "CA2"        "ACRBP"      "MMD"        "TREML1"    
[21] "F13A1"      "PGRMC1"     "SEPT5"      "MYL9"       "TSC22D1"   
[26] "MPP1"       "CMTM5"      "PTGS1"      "SNCA"       "RUFY1"     
[1] ""
 [1] "RPL10"  "RPS2"   "RPL11"  "RPL18A" "RPL32"  "RPL28"  "RPS19" 
 [8] "RPL12"  "RPL19"  "RPL13"  "RPS6"   "RPS14"  "RPS15"  "TMSB10"
[15] "RPLP1"  "RPL29"  "RPL6"   "RPL26"  "RPS4X"  "RPS16"  "RPS3"  
[22] "EEF1A1" "RPL15"  "RPL13A" "RPS7"   "RPL8"   "RPS12"  "RPL23A"
[29] "RPS18"  "RPLP2" 
[1] ""
[1] ""
[1] "PC4"
 [1] "CD79A"     "HLA-DQA1"  "CD79B"     "MS4A1"     "HLA-DQB1" 
 [6] "CD74"      "HLA-DPB1"  "HLA-DPA1"  "HLA-DRB1"  "TCL1A"    
[11] "HLA-DRA"   "LINC00926" "HLA-DQA2"  "HLA-DRB5"  "VPREB3"   
[16] "HLA-DMA"   "HLA-DMB"   "FCER2"     "BANK1"     "HVCN1"    
[21] "GZMB"      "HLA-DOB"   "PDLIM1"    "FCRLA"     "TSPAN13"  
[26] "FGFBP2"    "CD72"      "EAF2"      "PKIG"      "SPIB"     
[1] ""
 [1] "CD3D"    "LDHB"    "RPS14"   "IL7R"    "CD3E"    "RPL32"   "VIM"    
 [8] "IL32"    "RPS12"   "NOSIP"   "RPL28"   "GIMAP7"  "RPL11"   "RPL13"  
[15] "FYB"     "ZFP36L2" "RPL10"   "AQP3"    "JUNB"    "RPS25"   "RPLP1"  
[22] "MAL"     "LEF1"    "RGCC"    "S100A6"  "FOS"     "RPS3"    "CD2"    
[29] "RPL35A"  "RPL36"  
[1] ""
[1] ""
[1] "PC5"
 [1] "LTB"      "TMEM66"   "LDHB"     "HSPA8"    "PABPC1"   "NPM1"    
 [7] "RNASET2"  "JUNB"     "SOD1"     "CD52"     "NAP1L1"   "CALM2"   
[13] "EEF1A1"   "RPSA"     "HINT1"    "IL7R"     "VIM"      "FXYD5"   
[19] "AQP3"     "ARHGDIB"  "LDHA"     "HNRNPA1"  "NACA"     "HSP90AA1"
[25] "GSTK1"    "ITM2B"    "RPLP0"    "RPS19"    "ACTG1"    "RPS10"   
[1] ""
 [1] "S100A8"  "FGFBP2"  "GZMB"    "S100A9"  "NKG7"    "CCL4"    "S100A12"
 [8] "LGALS2"  "RBP7"    "GNLY"    "CST7"    "SPON2"   "GZMA"    "PRF1"   
[15] "CCL3"    "FOLR3"   "MS4A6A"  "GZMH"    "CD14"    "PRSS23"  "S1PR5"  
[22] "KLRD1"   "TYROBP"  "XCL2"    "CTSW"    "CCL5"    "CLIC3"   "CEBPD"  
[29] "TTC38"   "GSTP1"  
[1] ""
[1] ""
```

Seurat: Performing PCA
===
class:small-code

- Top 5 genes in the loading for PC1 and PC2.


```r
# Can plot top genes for top components
PrintPCA(pbmc.seurat, pcs.print = 1:2, genes.print = 5, use.full = TRUE)
```

```
[1] "PC1"
[1] "MALAT1"  "RPS27A"  "PTPRCAP" "RPS27"   "RPL23A" 
[1] ""
[1] "CST3"   "TYROBP" "LST1"   "AIF1"   "FCN1"  
[1] ""
[1] ""
[1] "PC2"
[1] "CD79A"    "MS4A1"    "TCL1A"    "HLA-DQA1" "HLA-DQB1"
[1] ""
[1] "NKG7" "GZMB" "PRF1" "CST7" "GZMA"
[1] ""
[1] ""
```

Seurat: PCA Visualizations
===
class:small-code

- Top 30 genes associated with the first two components.


```r
VizPCA(pbmc.seurat, pcs.use=1:2)
```

Seurat: PCA Visualizations
===
class:small-code

![plot of chunk unnamed-chunk-49](single_cell_analysis-figure/unnamed-chunk-49-1.png)

Seurat: PCA Visualizations
===
class:small-code


```r
PCAPlot(pbmc.seurat, 1, 2)
```

Seurat: PCA Visualizations
===
class:small-code

![plot of chunk unnamed-chunk-51](single_cell_analysis-figure/unnamed-chunk-51-1.png)

Seurat: PCA Visualizations
===
class:small-code

- Plot top 30 genes in top 100 cells for PC1.


```r
PCHeatmap(pbmc.seurat, pc.use=1, cells.use=100, do.balanced=TRUE)
```

Seurat: PCA Visualizations
===
class:small-code

![plot of chunk unnamed-chunk-53](single_cell_analysis-figure/unnamed-chunk-53-1.png)

Seurat: Choosing Components
===
class:small-code


```r
# Time Intensive
# Jackstraw
# pbmc.seurat <- JackStraw(pbmc.seurat, num.replicate = 100, do.print = FALSE)
```

Seurat: Choosing Components
===
class:small-code

How do we choose how many components to use?   
 - When there is diminishing returns to include it.
 - Selection is performed more liberally in our setting.


```r
# Scree (elbow) plot
PCElbowPlot(pbmc.seurat)
```

Seurat: Choosing Components
===
class:small-code

![plot of chunk unnamed-chunk-56](single_cell_analysis-figure/unnamed-chunk-56-1.png)

Seurat: Store Clusters
===
class:small-code

- Determine subclusters for the plot.
- Separate graph based approach.
   - is not aware of the t-SNE projection.
- Using PC 1-10


```r
# 1 minute
pbmc.seurat <- FindClusters(pbmc.seurat, pc.use = 1:10, resolution = 0.6, print.output = 0, save.SNN = TRUE)
```

Seurat: Run t-SNE
===
class:small-code

- Calculate and plot t-SNE on PC 1- 10.   
- Uses Barnes-hut implementation.


```r
# Calculate t-SNE Ordination
pbmc.seurat <- RunTSNE(pbmc.seurat, dims.use = 1:10, do.fast = TRUE)
# Plot
TSNEPlot(pbmc.seurat)
```

Seurat: Run t-SNE
===
class:small-code

![plot of chunk unnamed-chunk-59](single_cell_analysis-figure/unnamed-chunk-59-1.png)

---

This is not PCA
- Distance is best understood in close neighbors.
- The measurement of distance is difficult to understand.

Seurat: Side by side
===
class:small-code

**PCA**  
![plot of chunk unnamed-chunk-60](single_cell_analysis-figure/unnamed-chunk-60-1.png)

---

**t-SNE**  
![plot of chunk unnamed-chunk-61](single_cell_analysis-figure/unnamed-chunk-61-1.png)

Seurat: Plotting Genes Through Clusters
===
class:small-code

Now that we have subclusters of cell populations plotting genes through subclusters is identical to before.
- Seurat stores and groups by subclusters automatically.


```r
VlnPlot(pbmc.seurat, c("MS4A1","CD79A"))
```

Seurat: Plotting Genes Through Clusters
===
class:small-code

![plot of chunk unnamed-chunk-63](single_cell_analysis-figure/unnamed-chunk-63-1.png)

Seurat: Plotting Genes on Clusters
===
class:small-code

You can also gene expression through out the cell ordination.
- Marker genes can help identify cell groups.
- Metadata can help visualize batch affects.


```r
FeaturePlot(pbmc.seurat, c("MS4A1","CD3E", "GNLY", "FCER1A"), cols.use = c("grey","blue"))
```

Seurat: Plotting Genes on Clusters
===
class:small-code

![plot of chunk unnamed-chunk-65](single_cell_analysis-figure/unnamed-chunk-65-1.png)

Confirm Cell Identity
===
class:small-code

It is important to know cells are expressing expected genes.
- Many cell populations have known expressing genes.
- Many assays use sorting.
- Could have contaminants or unexpected diversity.

QC the Clusters!
===
class:small-code


```r
FeaturePlot(pbmc.seurat, c("nGene"), cols.use = c("grey","blue"))
```

QC the Clusters!
===
class:small-code

![plot of chunk unnamed-chunk-67](single_cell_analysis-figure/unnamed-chunk-67-1.png)

QC the Clusters!
===
class:small-code

We are going to make a fake batch affect (site) and plot this as an example of how one can visualize unwanted signal.


```r
# Making Fake Data
fake.sites <- as.integer(pbmc.seurat@ident %in% c(5,2,8,7))
names(fake.sites) <- colnames(pbmc.seurat@data)
# Add metadata
pbmc.seurat <- AddMetaData(pbmc.seurat, fake.sites, "site")
# Plot feature
FeaturePlot(pbmc.seurat, c("site"), cols.use = c("green","orange"))
```

QC the Clusters!
===
class:small-code

![plot of chunk unnamed-chunk-69](single_cell_analysis-figure/unnamed-chunk-69-1.png)

Seurat: Getting your labels
===


```r
cell.labels <- pbmc.seurat@ident
corner(cell.labels)
```

Seurat: Getting your labels
===


```
AAACATACAACCAC AAACATTGAGCTAC AAACATTGATCAGC AAACCGTGCTTCCG AAACCGTGTATGCG 
             1              3              0              5              6 
Levels: 0 1 2 3 4 5 6 7 8
```

Section Summary
===

- Dimensionality reduction help reduce data while *hopefully* keeping important signal.
  - t-SNE on PCA is often used in analysis
- Created several types of plot often seen in publications.
  - Plotting genes (through subgroups).
  - Ordinating cells in t-SNE space.
  - Heat maps of genes associated with PC components.
  - Plotting metadata on projects of data is an important QC tool.

Where are We at Now?
===

- Can plot genes vs genes, and throughout cells.
- Can plot cells.
- Can break plotted cells into subgroups.
  - These can become labels for later DE.
- Can plot across subgroups (maker genes).
- Can use plots of cells for batch affect visualization?

Single Cell RNA-Seq and Batch Affects
===

![fall](images/corgis/fall.jpg)

Single Cell RNA-Seq and Batch Affects
===

![confounding_paper](images/confounding_paper.png)

What is Confounding?
===

<div align="center">
<img src="images/confounding.png" height=500>
</div>

Confounding by Design
===

![batch_effects](images/batch_effects_4.png)

Normalization and Batch Affect Correction
===

- The nature of scRNA-Seq assays can make them prone to confounding with batch affects.
  - Normalization and batch affect correction can help.
- Some are moving away from relying on a specific method.
  - Exploring the idea of combining or selecting from a collection of normalization or correction methods best for a specific study.
- Some believe UMI based analysis need not be normalized between samples given the absolute count of the molecules are being reported.
- For advanced analysis try [SCONE](//niryosef.wordpress.com/tools/scone) or another tool from the [awesome list](https://github.com/seandavi/awesome-single-cell).
  - Be careful not to remove biological signal with good experimental design (avoiding confounding by design).

Differential Expression (SCDE)
===

![clump](images/corgis/corgi_clump.jpg)

SCDE: What is the Data?
===

- Data is from [Islam et al. 2011](http://genome.cshlp.org/content/21/7/1160.long)
- Approximately 100 samples.
  - Caution we use less...
- Mouse cells
  - Embryonic fibroblasts (MEF; control)
  - Embryonic stem cells (ES)

SCDE: Single Cell Differential Expression
===

![scde_intro](images/scde_intro.png)

SCDE: in Quick Theory
===

For each group (ES or MEF).
- Genes are modeled to have two groups of counts.
  - Noisy area highly prone to dropout (Poisson distribution).
  - "Amplified" signal (Negative Binomial distribution).
- This makes the error model to understand how much one can trust counts.
- Pairwise within groups.

Differential Expression.
- Expected value * the probability of dropout in that cell for that expression level.

SCDE: Load Data
===
class:small-code

- Pull in the data for the tutorial.
- In the future you would read in the data.
  - Read.table or load Seurat object.


```r
data(es.mef.small)
```

SCDE: Filter
===
class:small-code

- Filter the noise.
  - Almost every package has a function for this.
  - Filtering on gene prevalence.


```r
dim(es.mef.small)
```

```
[1] 14897    40
```

```r
cd <- clean.counts(es.mef.small, min.lib.size=1000, min.reads = 1, min.detected = 1)
dim(cd)
```

```
[1] 12142    40
```

SCDE: Create Labels
===
class:small-code

- Let's labels for the two groups.


```r
## Setting up cells groups
data.groups <- rep(NA, ncol(es.mef.small))
data.groups[ grep("MEF", names(es.mef.small)) ] <- "MEF"
data.groups[ grep("ESC", names(es.mef.small)) ] <- "ESC"
data.groups <- factor(data.groups, levels = c("ESC","MEF"))
names(data.groups) <- colnames(es.mef.small)
table(data.groups)
```

```
data.groups
ESC MEF 
 20  20 
```

SCDE: Calculate Error Models
===
class:small-code

- Some steps in SCDE can be time intensive.   
  - For the tutorial we will read in the result of this step.


```r
## Calculate error models
## Time Intensive step
# o.ifm <- scde.error.models(counts=cd, groups= data.groups, n.cores=4,
#                           threshold.segmentation=TRUE, save.crossfit.plots=FALSE,
#                           save.model.plots=FALSE, verbose=1)
## Precomputed
data(o.ifm)
```

SCDE: Calculate Error Models
===
class:small-code

- How do we speed this up?


```r
## Calculate error models
## Time Intensive step
# o.ifm <- scde.error.models(counts=cd, groups= data.groups, n.cores=4,
#                           threshold.segmentation=TRUE, save.crossfit.plots=FALSE,
#                           save.model.plots=FALSE, verbose=1)
# Check number of cores
detectCores()
```

```
[1] 8
```

SCDE: Model Details
===
class:small-code


```r
# Error model coefficients (cells = rows)
# corr.a = slope of the correlated component fit
# Negative corr.a could be bad cells
# corr.b intercept of the correlated component fit
# corr.theta is NB over-dispersion
# fail.r background poisson rate
head(o.ifm)
```

SCDE: Model Details
===
class:small-code


```
          conc.b    conc.a    fail.r    corr.b    corr.a corr.theta
ESC_10 -1.449443 0.5639140 -2.302585 0.7148157 0.6496142  0.7732069
ESC_11 -3.244421 0.7327046 -2.302585 1.5918205 0.5351960  0.7070433
ESC_12 -4.472559 0.8073935 -2.302585 1.5203470 0.4909147  0.7372590
ESC_13 -5.208909 0.8804523 -2.302585 1.2539230 0.5242493  0.8215473
ESC_14 -4.124369 0.7794612 -2.302585 1.1127353 0.5620266  0.7456712
ESC_15 -5.410838 0.9324758 -2.302585 1.1571732 0.5482784  0.7712750
```

SCDE: Filter Out Poor Fits
===
class:small-code


```r
dim(o.ifm)
```

```
[1] 40  6
```

```r
valid.cells <- o.ifm$corr.a > 0
table(valid.cells)
```

```
valid.cells
TRUE 
  40 
```

```r
o.ifm <- o.ifm[valid.cells, ]
dim(o.ifm)
```

```
[1] 40  6
```

SCDE: Estimate a Prior to Start
===
class:small-code

- Bayesian technique
- Need a starting point or first good guess to start.   


```r
## Calculate the Prior (starting value)
o.prior <- scde.expression.prior(models=o.ifm, counts=cd, length.out=400, show.plot=FALSE)
```

SCDE: Perform Differential Tests
===

- Why are we doing this again?


```r
## Setting up cells groups
data.groups <- rep(NA, nrow(o.ifm))
data.groups[ grep("MEF", rownames(o.ifm)) ] <- "MEF"
data.groups[ grep("ESC", rownames(o.ifm)) ] <- "ESC"
data.groups <- factor(data.groups, levels = c("ESC","MEF"))
names(data.groups) <- row.names(o.ifm)

## Perform T-test like analysis
# 2 minutes on standard computer, 1 core
#ediff <- scde.expression.difference(o.ifm, cd, o.prior, groups=data.groups, n.randomizations=100, n.cores=2, verbose=1)
load("data/ediff.Robj")
```

Visualize Differentially Expressed Genes
===

- mle = log2 fold change (estimate).
- ub and lb = upper and lower bound on mle.
- ce = log2 fold change (conservative estimate).
- Z = Z-score .
- cZ = Z-score corrected for multiple hypothesis testing.

---

![scde_output](images/scde_output_2.png)

SCDE: Top Upregulated Genes
===

- What are the top upregulated genes?


```r
head(ediff[order(ediff$Z, decreasing  =  TRUE), ])
```

```
                    lb      mle        ub       ce        Z       cZ
Dppa5a        8.075220 9.984631 11.575807 8.075220 7.160813 5.989598
Pou5f1        5.370220 7.200073  9.189043 5.370220 7.160328 5.989598
Gm13242       5.688455 7.677425  9.785734 5.688455 7.159979 5.989598
Tdh           5.807793 8.075220 10.302866 5.807793 7.159589 5.989598
Ift46         5.449779 7.359190  9.228822 5.449779 7.150242 5.989598
4930509G22Rik 5.409999 7.478528  9.785734 5.409999 7.115605 5.978296
```

SCDE: Write Data to Text File
===

- Save your data.   


```r
write.table(ediff[order(abs(ediff$Z), decreasing = TRUE), ],
            file = "data/scde_results.txt", row.names = TRUE, col.names = TRUE, sep = "\t", quote = FALSE)
```

SCDE: Plot a Gene
===

- Let's look at our data.   


```r
scde.test.gene.expression.difference("Tdh", models = o.ifm, counts = cd, prior = o.prior)
```

SCDE: Plot a Gene
===

![plot of chunk unnamed-chunk-85](single_cell_analysis-figure/unnamed-chunk-85-1.png)

```
          lb     mle       ub       ce        Z       cZ
Tdh 5.728235 8.03544 10.30287 5.728235 7.151425 7.151425
```

SCDE: Interactive Exploration
===

- Interactive exploration.


```r
# scde.browse.diffexp(ediff, o.ifm, cd, o.prior, groups = groups, name = "diffexp1", port = 1299)
```

Pathway analysis
===

![saturation](images/corgis/corgi_saturation.jpg)

Pagoda
===

![pagoda_paper](images/pagoda_paper.png)

Pagoda: Overview Methodology
===

<div align="center">
<img src="images/pagoda_methodology.png" width=900>
</div>
Fan et al.

Pagoda: What is the Data?
===

- Data being used is [Pollen et al. 2014](//www.ncbi.nlm.nih.gov/pubmed/25086649)
- Approximately 300 samples.
  - We focus on 64 samples (human, cortex)
  - NPC - Neuronal Progenitor Cells
  - GW16 - Gestational week 16
  - GW21 - Gestational week 21
  - GW21+3 - Gestational week 21 + 3 weeks.

Pagoda: Load Data
===
class:small-code


```r
data(pollen)
# Original genes and cells (count matrix)
dim(pollen)
```

```
[1] 23710    64
```

```r
# Filter poor cells
pollen.clean <- clean.counts(pollen)
# Cleaned matrix dimensions
dim(pollen.clean)
```

```
[1] 11310    64
```

Pagoda: Create a Color Legend
===
class:small-code


```r
name.keys <- gsub("^Hi_(.*)_.*", "\\1", colnames(pollen.clean))
name.keys
```

```
 [1] "NPC"    "NPC"    "NPC"    "NPC"    "NPC"    "NPC"    "NPC"   
 [8] "NPC"    "NPC"    "NPC"    "NPC"    "NPC"    "NPC"    "NPC"   
[15] "NPC"    "GW16"   "GW16"   "GW21"   "GW21+3" "GW21+3" "GW16"  
[22] "GW21+3" "GW21+3" "GW16"   "GW16"   "GW16"   "GW16"   "GW16"  
[29] "GW16"   "GW16"   "GW16"   "GW16"   "GW16"   "GW21"   "GW21"  
[36] "GW16"   "GW16"   "GW21"   "GW16"   "GW16"   "GW16"   "GW16"  
[43] "GW16"   "GW21"   "GW21"   "GW21"   "GW21+3" "GW16"   "GW16"  
[50] "GW16"   "GW16"   "GW16"   "GW21"   "GW21+3" "GW21+3" "GW21+3"
[57] "GW21+3" "GW21+3" "GW21+3" "GW21+3" "GW21+3" "GW21+3" "GW21+3"
[64] "GW21+3"
```

Pagoda: Create a Color Legend
===
class:small-code


```r
l2cols <- c("coral4", "olivedrab3", "skyblue2", "slateblue3")[as.integer(factor(name.keys, 
    levels = c("NPC", "GW16", "GW21", "GW21+3")))]
l2cols
```

```
 [1] "coral4"     "coral4"     "coral4"     "coral4"     "coral4"    
 [6] "coral4"     "coral4"     "coral4"     "coral4"     "coral4"    
[11] "coral4"     "coral4"     "coral4"     "coral4"     "coral4"    
[16] "olivedrab3" "olivedrab3" "skyblue2"   "slateblue3" "slateblue3"
[21] "olivedrab3" "slateblue3" "slateblue3" "olivedrab3" "olivedrab3"
[26] "olivedrab3" "olivedrab3" "olivedrab3" "olivedrab3" "olivedrab3"
[31] "olivedrab3" "olivedrab3" "olivedrab3" "skyblue2"   "skyblue2"  
[36] "olivedrab3" "olivedrab3" "skyblue2"   "olivedrab3" "olivedrab3"
[41] "olivedrab3" "olivedrab3" "olivedrab3" "skyblue2"   "skyblue2"  
[46] "skyblue2"   "slateblue3" "olivedrab3" "olivedrab3" "olivedrab3"
[51] "olivedrab3" "olivedrab3" "skyblue2"   "slateblue3" "slateblue3"
[56] "slateblue3" "slateblue3" "slateblue3" "slateblue3" "slateblue3"
[61] "slateblue3" "slateblue3" "slateblue3" "slateblue3"
```

Pagoda: Create Error Models
===
class:small-code

- This is a time intensive step.
- For now we are going to read in the data.


```r
# knn <- knn.error.models(pollen.clean, k=ncol(pollen.clean)/4, n.cores=2,
# min.count.threshold=2, min.nonfailed=5, max.model.plots=10) Precomputed
# data
data(knn)
```

Pagoda: Normalize the Variance
===
class:small-code


```r
# varinfo <- pagoda.varnorm(knn, counts=pollen.clean,
# trim=3/ncol(pollen.clean), max.adj.var=5, n.cores=1, plot=TRUE)
load("data/varinfo.Robj")

# list top overdispersed genes
sort(varinfo$arv, decreasing = TRUE)[1:10]
# Control for complexity
varinfo <- pagoda.subtract.aspect(varinfo, colSums(pollen.clean[, rownames(knn)] > 
    0))
```

Pagoda: Normalize the Variance
===
class:small-code


```r
# varinfo <- pagoda.varnorm(knn, counts=pollen.clean,
# trim=3/ncol(pollen.clean), max.adj.var=5, n.cores=2, plot=TRUE)
load("data/varinfo.Robj")
```

Pagoda: Normalize the Variance
===
class:small-code

<div align="center">
<img src="images/varinfo.png" width=800>
</div>

Pagoda: Normalize the Variance
===
class:small-code


```r
# list top overdispersed genes
sort(varinfo$arv, decreasing = TRUE)[1:10]
```

```
     DCX     EGR1      FOS  IGFBPL1   MALAT1    MEF2C    STMN2    TOP2A 
5.000000 5.000000 5.000000 5.000000 5.000000 5.000000 5.000000 5.000000 
  BCL11A     SOX4 
4.739497 4.489101 
```

```r
# Control for complexity
varinfo <- pagoda.subtract.aspect(varinfo, colSums(pollen.clean[, rownames(knn)] > 
    0))
```

Pagoda: Load Precomputed Gene Groups
===
class:small-code

- This analysis path relies on known gene groups.   


```r
library(org.Hs.eg.db)
# translate gene names to ids
ids <- unlist(lapply(mget(rownames(pollen.clean), org.Hs.egALIAS2EG, ifnotfound = NA), 
    function(x) x[1]))
rids <- names(ids)
names(rids) <- ids
```

Pagoda: Load Precomputed Gene Groups
===
class:small-code


```r
# convert GO lists from ids to gene names
gos.interest <- unique(c(ls(org.Hs.egGO2ALLEGS)[1:100], "GO:0022008", "GO:0048699", 
    "GO:0000280", "GO:0007067"))
go.env <- lapply(mget(gos.interest, org.Hs.egGO2ALLEGS), function(x) as.character(na.omit(rids[x])))
go.env <- clean.gos(go.env)  # remove GOs with too few or too many genes
go.env <- list2env(go.env)  # convert to an environment
```

Pagoda: Known Gene Groups Overdispersion
===
class:small-code


```r
# pwpca <- pagoda.pathway.wPCA(varinfo, go.env, n.components = 1, n.cores =
# 1)
load("data/pwpca.Robj")
df <- pagoda.top.aspects(pwpca, return.table = TRUE, plot = TRUE, z.score = 1.96)
```

Pagoda: Known Gene Groups Overdispersion
===
class:small-code

![plot of chunk unnamed-chunk-97](single_cell_analysis-figure/unnamed-chunk-97-1.png)

Pagoda: Known Gene Groups Overdispersion
===
class:small-code


```r
head(df)
```

Pagoda: Novel Gene Group Overdispersion
===
class:small-code

- Unsupervised path.   


```r
# clpca <- pagoda.gene.clusters(varinfo, trim = 7.1/ncol(varinfo$mat),
# n.clusters = 50, n.cores = 1, plot = TRUE)
load("data/clpca.Robj")
```

Pagoda: Novel Gene Group Overdispersion
===
class:small-code

<div align="center">
<img src="images/clpca.png" width=500>
</div>

Pagoda: Novel Gene Group Overdispersion
===
class:small-code


```r
df <- pagoda.top.aspects(pwpca, clpca, return.table = TRUE, plot = TRUE, z.score = 1.96)
```

Pagoda: Novel Gene Group Overdispersion
===
class:small-code

![plot of chunk unnamed-chunk-101](single_cell_analysis-figure/unnamed-chunk-101-1.png)

Pagoda: Novel Gene Group Overdispersion
===
class:small-code


```r
head(df)
```

Pagoda: Novel Gene Group Overdispersion
===
class:small-code


```
             name npc   n    score         z     adj.z sh.z adj.sh.z
56  geneCluster.6   1 391 3.289412 13.814757 13.530139   NA       NA
49     GO:0000280   1 428 1.587359 11.592060 11.252131   NA       NA
50     GO:0007067   1 362 1.585331 10.873998 10.576419   NA       NA
79 geneCluster.29   1 171 1.583733  5.162516  4.523195   NA       NA
14     GO:0000070   1 116 1.446979  5.625549  5.247594   NA       NA
89 geneCluster.39   1 178 1.392980  4.137589  3.438496   NA       NA
```

Pagoda: Cell Clustering
===
class:small-code


```r
# Get full info on the top aspects tam <- pagoda.top.aspects(pwpca, clpca,
# n.cells = NULL, z.score = qnorm(0.01/2, lower.tail = FALSE))
load("data/tam.Robj")
# Determine overall cell clustering
hc <- pagoda.cluster.cells(tam, varinfo)
# tamr <- pagoda.reduce.loading.redundancy(tam, pwpca, clpca)
load("data/tamr.Robj")
# tamr2 <- pagoda.reduce.redundancy(tamr, distance.threshold = 0.9, plot =
# TRUE, cell.clustering = hc, labRow = NA, labCol = NA, box = TRUE, margins
# = c(0.5, 0.5), trim = 0)
load("data/tamr2.Robj")
col.cols <- rbind(groups = cutree(hc, 3))
```

Pagoda: Visualize Clustering
===
class:small-code


```r
pagoda.view.aspects(tamr2, cell.clustering = hc, box = TRUE, labCol = NA, margins = c(0.5, 
    20), col.cols = rbind(l2cols))
```

Pagoda: Visualize Clustering
===
class:small-code

![plot of chunk unnamed-chunk-106](single_cell_analysis-figure/unnamed-chunk-106-1.png)

Pagoda: Interactive Clustering
===
class:small-code


```r
## compile a browsable app, showing top three clusters with the top color bar
## app <- make.pagoda.app(tamr2, tam, varinfo, go.env, pwpca, clpca, col.cols
## = col.cols, cell.clustering = hc, title = 'NPCs') show app in the browser
## (port 1468) show.app(app, 'pollen', browse = TRUE, port = 1468)
```

Section Summary
===

- Generated groups of cells (labels) from subclusters.
- DE can be performed on cell subclusters.
  - Addresses the specific characteristics of the underlying data.
  - Is a t-test not an ANOVA
  - Can be time intensive.
- Pathway enrichment analysis can be performed using de novo and established gene groupings.

What Did We Miss? (So Much)
===

- Seurat
  - Data imputation, spatial inference, batch effect correction, ...
- SCDE
  - Batch Effect Correction, ...
- Pseudotemporal Analysis (time-series)
  - Monocle2, ...
- So many more tools are emerging, ...
-- https://github.com/seandavi/awesome-single-cell

What did we accomplish?
===

- Looked at the data and performed filtering and QC
- Plot genes and cells.
- Discovered subgroups of cells.
- Groups of cells were used as labels in DE.
- Groups of genes were analyzed for enrichment.

Thank You
===

- Aviv Regev
- Brian Haas
- Anindita Basu
- Asma Bankapur
- Chloe Villani
- Karthik Shekhar

---

- Matan Hofree
- Michel Cole
- Monika Kowalczyk
- Nir Yosef
- Regev Single Cell Working Group
- Today's Attendees

Questions?
===

![gradute corgi](images/corgis/graduate_corgi.jpg)

Notes: to Make a PDF
===
class:small-code

- Create a pdf file before you plot ( can plot multiple plots ).
- Close the plotting.


```r
# pdf( 'data/my_file.pdf', useDingbats = FALSE ) # Start pdf plot( 1:10,
# log(1:10 ) ) # plot in to the pdf file plot( seq(0,.9,.1), sin(0:9) ) #
# another plot for the pdf file dev.off() # Close pdf file ( very important
# )
```
