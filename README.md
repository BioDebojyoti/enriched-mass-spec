# enriched-mass-spec

# EnrichedMassSpec

A **web-based enrichment application** built to process Spectranaut mass-spectrometry outputs and perform downstream differential expression and enrichment analyses with interactive visualizations. The application is freely available and is hosted at SciLifeLab Serve: [https://enrichedmassspec.serve.scilifelab.se/](https://enrichedmassspec.serve.scilifelab.se/)

---

##  Features

- **Data upload**: Accepts Spectranaut output files (tab-separated).
- **Column selection**: Dynamic dropdowns to map key data columns.
- **Volcano plot visualization**: Interactive thresholding with user-defined cutoffs.
- **Downstream enrichment**: Functional analysis based on selected thresholds.

---

##  Input Data Format

- **Accepted file type**: Tab-separated values (e.g., `.tsv` or `.txt`), as output from Spectranaut.
- **Required columns** (identified via dropdowns):

  1. **Log₂ fold-change** — numeric values for differential expression (log2FC).
  2. **P-value** — statistical significance values (numeric).
  3. **Protein** — identifier, e.g., protein name or accession.
  4. **Gene** — gene symbol/name.
  5. **UniProt ID** — UniProt accession string.
  6. **Protein description** — textual description for context.

> The app dynamically inspects the uploaded file to populate dropdown menus, allowing users to assign their columns to these required data types.

---

##  Workflow Overview

### 1. Data Upload & Column Mapping
- Upload your Spectranaut `.tsv`.
- Use dropdowns to specify:
  - Log₂FC column  
  - P-value column  
  - Protein identifier column  
  - Gene symbol column  
  - UniProt ID column  
  - Protein description column  

### 2. Volcano Plot Parameters
- Define numerical thresholds via input fields:
  - **Log₂FC cutoff** — e.g., `±1`
  - **P-value cutoff** — e.g., `0.05`

### 3. Volcano Plot Generation
- Based on thresholds, the app generates an interactive volcano plot:
  - **X-axis**: log₂FC  
  - **Y-axis**: –log₁₀(p-value)  
  - Visual highlighting of significant up- and down-regulated proteins.

### 4. Downstream Enrichment Analysis
- The set of significant proteins (over cutoff thresholds) is subjected to enrichment analysis:
  - Selected functional databases (e.g., GO terms, pathways, etc.; adapt as available in the app).
  - Results typically include enriched categories with p-values, fold-enrichment, and visualizations (e.g., bar charts, tables).

---

##  Tab Descriptions *(Based on Application UI)*

> _As the webpage wasn’t accessible directly, these tab descriptions are inferred from typical Shiny/Dash layout and your description. If there are additional tabs (e.g., filtering, plots, export), you can update accordingly._

| Tab Name             | Purpose |
|----------------------|---------|
| **Upload / Input**   | Upload Spectranaut file and map columns via dropdowns. |
| **Volcano Plot**     | Visualize differential expression; apply log₂FC and p-value cutoffs and interactively explore significant hits. |
| **Enrichment**       | Run enrichment analysis on proteins passing thresholds; displays enriched terms, statistics, and visual summaries. |
| **Settings / Parameters** | (If present) Fine-tune analysis choices, e.g., select which databases to use for enrichment, adjust multiple test correction method, etc. |
| **Results / Export** | (If present) View and download tables or plots—e.g., volcano plot PNG, enrichment results CSV or JSON. |

---

##  Example Usage

```bash
# 1. Prepare your Spectranaut tsv file
# 2. Navigate to the app and upload the file
# 3. Map columns:
#    - Log₂FC → choose column A
#    - P-value → choose column B
#    - Protein → choose column C
#    - Gene → choose column D
#    - UniProt ID → choose column E
#    - Protein description → choose column F
# 4. Set thresholds:
#    - Log₂FC cutoff = 1
#    - P-value cutoff = 0.05
# 5. Generate volcano plot; inspect highlighted proteins
# 6. Run enrichment analysis on selected proteins
# 7. Explore results and export plots or tables
