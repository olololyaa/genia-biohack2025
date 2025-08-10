# GAI - Gene Analysis & Insights

GAI - Gene Analysis & Insights is a new web-based data base. Here you can find anything you want to know about genes via intuitive GUI, construct a dataset using flexible filters, and download it for data analysis.

## Integrated data bases
GAI integrated the most valueable knowledge from a number of comprehensive data bases:

### Ensembl

Ensembl is a genome browser for vertebrate genomes that supports research in comparative genomics, evolution, sequence variation and transcriptional regulation.

### gnomAD

The Genome Aggregation Database (gnomAD™), originally launched in 2014 as the Exome Aggregation Consortium (ExAC), is the result of a coalition of investigators willing to share aggregate human exome and genome sequencing data from a variety of large-scale sequencing projects, and make summary data available for the wider scientific community.

### DGIdb

The Drug-Gene Interaction Database (DGIdb) streamlines the search for druggable therapeutic targets through the aggregation, categorization, and curation of drug and gene data from publications and expert resources. Containing over 10,000 genes and 20,000 drugs involved in over 70,000 drug-gene interactions, DGIdb facilitates research and clinical decision-making and acts as a comprehensive tool for exploring the druggable genome.

### DrugCentral

DrugCentral is online drug information resource created and maintained by Division of Translational Informatics at University of New Mexico in collaboration with the IDG. DrugCentral provides information on active ingredients chemical entities, pharmaceutical products, drug mode of action, indications, pharmacologic action.

### HOCOMOCO

HOmo sapiens COmprehensive MOdel COllection (HOCOMOCO) v11 provides transcription factor (TF) binding models for 680 human and 453 mouse TFs.

### dGTEx

The Developmental Genotype-Tissue Expression (dGTEx) Project is a new effort to study development-specific genetic effects on gene expression from donors at each of 4 developmental windows: infant, early childhood, pre-pubertal, post-pubertal.

### CORUM

Protein complexes are key molecular entities that integrate multiple gene products to perform cellular functions. The CORUM database is a collection of experimentally verified mammalian protein complexes.

### The Human Protein Atlas

The Human Protein Atlas is a Sweden-based program initiated in 2003 with the aim to map all the human proteins in cells, tissues, and organs using an integration of various omics technologies, including antibody-based imaging, mass spectrometry-based proteomics, transcriptomics, and systems biology. All the data in the knowledge resource is open access to allow scientists both in academia and industry to freely access the data for exploration of the human proteome.

### HPO

The Human Phenotype Ontology (HPO) project provides an ontology of medically relevant phenotypes, disease-phenotype annotations, and the algorithms that operate on these. The HPO can be used to support differential diagnostics, translational research, and a number of applications in computational biology by providing the means to compute over the clinical phenotype. The HPO is being used for computational deep phenotyping and precision medicine as well as integration of clinical data into translational research.

### Predicted Transcription Factors

Lambert et al., 2018, https://doi.org/10.1016/j.cell.2018.01.029

### Copy-number variations related to disease

Rayan L. et al., 2022, https://doi.org/10.1016/j.cell.2022.06.036

## Development Setup

### 1. Create a Virtual Environment

#### Windows (PowerShell)
```powershell
python -m venv venv
./venv/Scripts/Activate.ps1
```
If activation is blocked:
```powershell
# Run PowerShell as Administrator
Set-ExecutionPolicy -Scope CurrentUser -ExecutionPolicy RemoteSigned
```

#### macOS / Linux
```bash
python -m venv venv
source venv/bin/activate
```

### 2. Install Dependencies
```bash
pip install --upgrade pip
pip install -r requirements.txt
```

### 5. Run the Streamlit App
```bash
streamlit run app.py
```
The app will auto-reload on file save.

### 6. Debugging
Inline breakpoint:
```python
import pdb; pdb.set_trace()
```
