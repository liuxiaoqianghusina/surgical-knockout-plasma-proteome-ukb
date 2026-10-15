# Surgical organ removal and the tissue origins of plasma proteins in UK Biobank

Analysis code for the manuscript *"Surgical organ removal reveals the tissue origins of plasma proteins and challenges proteomic organ clocks"*.

The study treats 24 types of organ removal surgery as natural human knockouts and relates them to 2,923 plasma proteins (Olink Explore 3072) in 52,054 UK Biobank participants. Organ loss is separated from the indication for surgery by comparing people operated on before blood sampling with people who underwent the same operation only after sampling. Further analyses compare resections of different extent within the same organ, examine within-person change in participants operated on between two Olink blood draws, compare the changes with tissue RNA specificity from the Human Protein Atlas, and test published proteomic organ-ageing clocks.

## Data

This research was conducted using the UK Biobank Resource under application number 1099827. No individual-level data, results or manuscript files are included in this repository. UK Biobank data are available to approved researchers through https://www.ukbiobank.ac.uk.

Public resources used:

| Resource | Use | Source |
|---|---|---|
| Human Protein Atlas, RNA tissue and single-cell type specificity | Tissue-origin analysis (`12_hpa_tissue_specificity.py`) | https://www.proteinatlas.org (download via the search API with columns `g,eg,up,rnats,rnatd,rnatsm,rnascs,rnascsm,rnabcs,rnabcsm`) |
| Oh et al., *Nat Med* 2025, Supplementary Tables (ST2 organ-enriched Olink proteins; ST3 organ clock weights) | Organ clocks | doi:10.1038/s41591-025-03798-1 |
| Wang et al., *Nat Aging* 2025, Supplementary Tables (ST4 clock proteins; ST5 reported correlations with age) | Second organ clock set (`21_wang_clocks.py`) and organ-enriched protein lists | doi:10.1038/s43587-025-01016-8 |

## Setup

Paths are written as placeholders. Before running, replace `/path/to/ukb_data` with the folder holding the UK Biobank extract and `/path/to/project/cache` with the project's `cache/` folder. Run the scripts from `scripts/`; intermediate files are written to `../cache/`, results to `../results/` and figures to `../figures/`. Scripts that read downloaded public files take their paths as command-line arguments, and should be run with `python3 -P` so that no module is imported from the data folders.

The scripts expect comma-separated tables keyed by `eid`, with columns named by field (`p<field>_i<instance>_a<array>`), as exported from the UK Biobank Research Analysis Platform. The file names below are generic; rename your own exports or edit the paths.

| File | Contents |
|---|---|
| `ukb_assessment_center.csv` | Baseline assessment fields: self-reported operations (20004, 20011), medications (20003), anthropometry, smoking, fasting time, imaging-visit dates |
| `ukb_population_characteristics.csv` | Sex, year and month of birth, Townsend deprivation index |
| `ukb_biological_samples.csv` | Serum creatinine and cystatin C; Olink plate and consortium-selection fields |
| `ukb_health_outcomes.csv` | Hospital diagnoses (ICD-10) and procedures (OPCS-4) with dates, cancer registry, death |
| `olink_instance0.csv`, `olink_instance2.csv`, `olink_instance3.csv` | Olink NPX values at baseline and the two imaging visits (`Participant.ID` plus one column per protein) |
| `ukb_codings.csv`, `ukb_data_dictionary.csv` | Coding tables and field dictionary |

Analyses were run with Python 3.14 (pandas 2.3.3, NumPy 2.4.6, SciPy 1.17.1, statsmodels 0.14.6, lifelines 0.30.3, matplotlib 3.11.0, openpyxl, Pillow).

## Pipeline

| Script | Purpose |
|---|---|
| `01_extract_medications.py` – `07_extract_repeat_visit_dates.py` | Extract medications, self-reported operations and covariates, extra covariates and biochemistry, OPCS-4 procedure dates (including resection extent), hospital outcomes and repeat-visit dates |
| `b_common.py` | Operation definitions, cohort construction (removed before blood draw, removed after blood draw, never operated), covariate design and exact complete-case regression grouped by missingness pattern |
| `08_atlas.py` | Operation × protein atlas with global FDR control, including within-organ contrasts |
| `09_indication_comparator_sensitivity.py` | Organ loss versus indication (operated-after comparator), sensitivity analyses and time since surgery |
| `10_repeat_olink_replication.py` | Replication in repeat Olink measurements |
| `11_within_person.py` | Within-person change in participants operated on between two blood draws |
| `12_hpa_tissue_specificity.py` | Specificity and sensitivity of tissue RNA enrichment |
| `13_organ_clocks.py` | Organ ages from published clock weights and their response to organ removal |
| `14_bowel_extent_kidney_function.py` | Extent of bowel resection; kidney function (filtration signature, eGFR adjustment) |
| `15_tonsillectomy_proteome.py`, `16_tonsillectomy_outcomes.py` | Childhood tonsillectomy: plasma proteome and incident outcomes |
| `17_clock_and_tonsil_supplementary.py` | eGFR-adjusted clocks, clock decomposition, organ-enriched protein fractions, tonsillectomy sensitivity analyses |
| `18_clock_specificity_index.py` | Organ-clock specificity index with confidence intervals |
| `19_table1.py` | Table 1 |
| `21_wang_clocks.py` | Second organ clock set: weights re-estimated for the clock proteins of Wang et al. (tenfold cross-fitting), organ age gaps, removal effects and specificity index; run before `20_figures.py` |
| `20_figures.py` | Main and supplementary figures (PDF and 600-dpi TIFF) |

## License

MIT License (see `LICENSE`).
