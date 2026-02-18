# Callus-co-cultivation
Code and datasets for the manuscript:“Augmenting the culturability of rice root bacteria through callus co-cultivation in high throughput culturomics”

## 01_Scripts
R scripts and QIIME 2 pipeline used for analysis.

- `qiime2.sh`: QIIME 2 processing pipeline.
- `Figure1B_Metabolite class 2 distribution.R`
- `Figure1C_Top compounds.R`
- `Figure1D_Pathway enrichment analysis.R`
- `Figure2B_Bacterial_Growth.R`
- `Figure2CD_Diversity.R`

## 02_Data
Raw and processed data tables.

- `manifest202406.csv`
- `metadata.tsv`
- `Figure2A_RawData（Number of bacterial growing wells per plate）.csv`
- `Figure2A_RawData_Analysis.csv`
- `Figure2B_Bacterial_Growth_RawData.csv`
- `Figure2B_RawData_Analysis.csv`
- `Figure2CD_Diversity_analysis.csv`
- `Figure2E_RawData(Family_Frequency of bacterial occurrence).csv`
- `Figure2E_RawData(Genus_Frequency of bacterial occurrence).csv`
- `Figure2F_RawData(R10723_Bacterial growth curve) .csv`
- `Figure2F_RawData(R10727_Bacterial growth curve) .csv`
- `Figure2H_(Laboratory system_plant height).csv`
- `Figure2I(Laboratory system_root length).csv`
- `Figure2K_RawData(Pot_Plant Height).csv`
- `Figure2L_RawData(Pot_Root length).csv`
- `Figure2M_RawData(Pot_Fresh Weight).csv`

## 03_Figures
Generated figures and original GraphPad Prism files.

- **PDF Output**:
  - `Figure1B_ClassII_PieChart.pdf`
  - `Figure1C_Top Compounds.pdf`
  - `Figure1D_Pathway enrichment analysis.pdf`
  - `Figure2B_(Number of bacterial growing wells).pdf`
  - `Figure2C_Alpha_Diversity.pdf`
  - `Figure2D_Beta_Diversity.pdf`

- **GraphPad Prism 8.0.2 Source Files** (Located in subfolder):
  - `Figure2A(Number of bacterial growing wells per plate).pzfx`
  - `Figure2E(Frequency of bacterial occurrence).pzfx`
  - `Figure2F(Bacterial growth curve).pzfx`
  - `Figure2H-I(Laboratory inoculation phenotype).pzfx`
  - `Figure2K-M(Potted plant experimental phenotypes).pzfx`

## Software Versions

**R Environment**
- Version: R 4.4.1
- Platform: x86_64-w64-mingw32

**Microbiome Analysis**
- QIIME 2: q2cli version 2024.5.0

**Statistical Plotting**
- GraphPad Prism: 8.0.2
