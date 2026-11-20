# NACC data guide for Track C

## Source and receipt

The National Alzheimer's Coordinating Center (NACC) is the official data source. The statistics course teaching assistants provided the three CSVs below as the Track C teaching subset. The extract version, data collection period, and preparation methods, including how UDS values were imputed, are unknown.

[NACC Researcher's Guide](https://www.naccdata.org/the-nacc-researchers-guide/) · [Course guide, p. 4](../docs/amarpr.pdf#page=4)

The files were received on **October 10**. The year and the NACC subset release date are undocumented. Participants were reported to have been selected at random. The sampling frame and eligibility criteria for that selection are unavailable.

## Supplied files

Place these files in `data/` with the names shown below.

| File | Contents and unit of observation |
|---|---|
| `udslong_imputed.csv` | Clinical, cognitive, and demographic UDS data with 1,237 participant rows and 180 columns. |
| `mriscanclean.csv` | MRI measurements for the same 1,237 participant IDs, with 126 columns. The primary thickness measure is in millimeters. |
| `combined_adjacency_matrix.csv` | A symmetric 62 × 62 binary adjacency matrix with an additional row-label column. The audit uses it to match the names of brain regions. |

The UDS data dictionary defines the `ANIMALS` score and its special codes. The MRI data dictionary defines the left cortical thickness measure and its unit. [UDS data dictionary](https://r2d2.kumc.edu/ADC/Protocols/UDS2_protocol/ResearchDataDict_2015Mar.pdf#page=205) · [MRI data dictionary](https://files.alz.washington.edu/documentation/SCAN-MRI-Imaging-RDD.pdf#page=17)

## Reconstruct the analysis dataset

Run these commands from the project root with Python 3.10 or later.

```shell
python3 -m src.data_audit > reports/phase1/data_audit.json
python3 -m src.build_dataset
```

The first command checks special codes, value ranges, unique IDs, and file structure. The second joins the files one-to-one on `NACCID` and creates `data/processed/participants.csv` with 1,237 rows. The derived table retains the original values, the usable values, reasons for unusable entries, and sample inclusion flags. The primary modeling sample includes **842 participants** with normal cognition. The files provide neither the scan date nor the cognitive visit date or number, so the time between measurements cannot be calculated. [Audit output](../reports/phase1/data_audit.json) · [Sample flow](../reports/phase1/sample_flow.json) · [Phase 1 build instructions](../README.md)
