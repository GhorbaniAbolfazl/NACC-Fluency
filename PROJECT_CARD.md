# Project card | Track C, NACC

## 1. Title and primary question

**Title.** Left inferior parietal cortical thickness and semantic verbal fluency in a NACC teaching subset.

**Research question.** Among participants with normal cognition, how many more or fewer animals are named correctly in 60 seconds per 0.1 mm greater left inferior parietal cortical thickness, after accounting for age, sex, and years of education?

This is a secondary observational analysis at the participant level. The working hypothesis is that greater cortical thickness is associated with a higher animal-naming score. The modeling phase will report the estimated association and its 95% confidence interval. Estimates for the full sample and comparisons across cognitive status groups are supplementary. The question and planned analyses were documented before fitting models for this project.

## 2. Data source and unit of observation

`udslong_imputed.csv` and `mriscanclean.csv` come from the teaching subset supplied for the course. The UDS file contains clinical assessments and cognitive test results. The MRI file contains measurements derived from imaging. `combined_adjacency_matrix.csv` describes 62 brain regions. It is documented in the data inventory and is used to check region labels. The primary model uses left cortical thickness from the MRI file.

| Data source | Data rows | Columns | Unique participant IDs |
|---|---:|---:|---:|
| UDS | 1,237 | 180 | 1,237 |
| MRI | 1,237 | 126 | 1,237 |

Both CSVs contain the same 1,237 participant IDs. These counts are taken directly from the input files before eligibility rules are applied. Each supplied file has one record per participant, so the join on `NACCID` produces one matched row per person. The reviewed columns contain no cognitive test date, scan date, or index visit number. The number of visits in the source database and the interval between the MRI and cognitive assessment are unknown.

The student reports that the participants in the teaching subset were selected at random. The files do not specify the sampling frame or the eligibility criteria for that selection. Interpretation is limited to the supplied subset and the defined analysis samples. Every person in the UDS file has a corresponding MRI record. The sample flow tracks MRI record coverage and whether the thickness measurement is usable as two separate quantities.

## 3. Variables and their roles

| Variable | Source | Definition and scale | Role |
|---|---|---|---|
| `NACCID` | UDS and MRI | Nominal participant identifier | Record linkage |
| `ANIMALS` | UDS | Number of correct responses in 60 seconds, discrete | Primary outcome |
| `LH_INFERIORPARIETAL_AVGTH` | MRI | Mean thickness of the left inferior parietal cortex, continuous, in mm | Primary explanatory variable |
| `NACCAGE` | UDS | Age at the visit, quantitative, in years | Covariate |
| `SEX` | UDS | Recorded sex, categorical | Covariate |
| `EDUC` | UDS | Years of education, discrete | Covariate |
| `NACCUDSD` | UDS | Recorded cognitive status, categorical | Defines the normal-cognition sample and comparison groups |

For the proposed sensitivity analyses, `VEG` provides an alternative measure of semantic fluency, and `RH_INFERIORPARIETAL_AVGTH` can be averaged with the left-side measurement. The provenance and preparation of individual values in the supplied UDS file are undocumented.

## 4. Audit rules and findings

The [audit script](src/data_audit.py) produces the [audit results and column inventory](reports/phase1/data_audit.json). The inventory covers all 180 UDS columns. For the MRI file it lists gray matter volume in cubic millimeters and cortical thickness in millimeters for each of the 62 regions across the two hemispheres. The region labels agree with the supplied matrix. Value and range checks focus on the primary variables and the measures designated for sensitivity analyses.

| Variable | Rule for usable values | Findings in the 1,237 records per file |
|---|---|---|
| `ANIMALS` | Whole numbers from 0 to 77. Codes 95–98 and −4 are recorded separately. | 1,225 usable. Code 96 occurred twice, code 97 eight times, and code 98 once. One value was 21.6. |
| `VEG` | Same count range and special codes as `ANIMALS`. | 1,222 usable, with 14 special codes and one value of 16.6. |
| `NACCAGE` | Age at the visit from 18 to 120 years. | All usable. Observed range 41–91 years. |
| `SEX` | Code 1 for male and 2 for female. | All usable. |
| `EDUC` | Whole years from 0 to 36. Code 99 means unknown. | 1,235 usable and two records coded 99. |
| `NACCUDSD` | 1 for normal cognition, 2 for impairment without MCI, 3 for MCI, and 4 for dementia. | All in range. The UDS file contains 847 participants classified as cognitively normal before eligibility filters. |
| Left and right thickness | Finite positive measurements. The rule uses no upper cutoff derived from the sample. | All passed the numerical checks. Left thickness ranged from 1.349 to 2.770 mm. |

For the fluency tests, codes 96, 97, and 98 indicate a cognitive or behavioral problem, another problem, and refusal, respectively. Fractional scores are flagged separately from missing codes because they do not fit the count scale. They are kept as supplied, with no rounding, and excluded from the primary analysis. Their origin is unknown. The thickness checks establish numerical usability. The teaching subset has no documented scan-level quality control information.

Both files had 1,237 unique, nonblank IDs. The audit table reports **variable-specific denominators before the analysis join**. The sample flow gives the denominators used for modeling.

## 5. Record linkage and analysis samples

The [dataset build](src/build_dataset.py) joins the files one-to-one on `NACCID`. All 1,237 UDS participants have a corresponding MRI record, and the two files agree on cognitive status for every matched participant. Group labels come from UDS. The build stops before writing outputs if it finds a blank or duplicate ID, or if two usable cognitive status codes disagree for a matched participant. If only the MRI status is unusable, the build keeps the UDS group label and flags the cross-file comparison as unavailable.

| Sample | Eligibility | Participants |
|---|---|---:|
| Matched participant table | ID present in both files | 1,237 |
| Cognitive status and MRI | Usable cognitive status and left thickness | 1,237 |
| Paired sample, all groups | Usable `ANIMALS` score and left thickness | 1,225 |
| Paired sample, normal cognition | Both measures available and cognitive status classified as normal | 843 |
| Supplementary modeling sample, all groups | Both measures plus usable age, sex, and education | 1,223 |
| Primary modeling sample, normal cognition | Both measures and all covariates available, with cognitive status classified as normal | 842 |

Of the 847 participants with normal cognition, four lacked a usable animal-naming score and one additional participant had unknown education. Across the full sample, 12 participants lacked a usable score and two additional participants had unknown education. All three Phase 4 models will use the same 842 participants in the primary comparison, so changes in fit and estimates reflect the model specifications on a common sample.

The [sample flow](reports/phase1/sample_flow.json) gives sample sizes by cognitive status and the first reason for exclusion from each sample. The first-reason counts count each excluded person once. The report also tallies unusable values by column so overlapping reasons remain visible. Every UDS participant has an MRI record and a usable left thickness measurement. The supplied files therefore provide no comparison group of UDS participants without MRI.

## 6. Descriptive statistics and group comparisons

The [descriptive results](reports/phase1/descriptive_statistics.json) report means, sample standard deviations, minima, quartiles, medians, maxima, and interquartile ranges for animal-naming scores, left cortical thickness, age, and years of education. Sex and cognitive status are summarized using counts and percentages. Percentages use denominators for both the full group and participants with usable values. The output also includes observed and cumulative frequencies for animal-naming scores and education in the full sample and the normal-cognition group.

Means and standard deviations summarize center and spread in the original units, which makes the groups easier to compare. The median and interquartile range describe the center and the middle 50% with less sensitivity to extreme observations. The mean and median of the animal-naming scores are close. The thickness distribution has a tail toward lower values, so the median and interquartile range give a useful description of its central portion.

The counts beside the summaries below are the denominators for those variables. The total in each row counts all matched participants with that cognitive status.

| Cognitive status | Group total | Usable scores | Animal-naming score, mean ± SD | Usable thickness values | Left thickness, mean ± SD (mm) |
|---|---:|---:|---:|---:|---:|
| Normal cognition | 847 | 843 | 22.32 ± 5.43 | 847 | 2.336 ± 0.121 |
| Impaired without MCI | 7 | 7 | 22.00 ± 4.16 | 7 | 2.314 ± 0.105 |
| MCI | 252 | 249 | 17.83 ± 5.56 | 252 | 2.258 ± 0.173 |
| Dementia | 131 | 126 | 12.27 ± 5.26 | 131 | 2.154 ± 0.214 |

Among participants with normal cognition, the middle half of the scores ranged from 18 to 26 correct responses, and the middle half of left thickness measurements ranged from 2.269 to 2.417 mm. The mean age was 67.27 years in the normal-cognition group and 70.83 years in the MCI group. The group summaries are unadjusted. The impairment-without-MCI group contains seven people, so a few observations can noticeably affect its summaries.

Quartiles use linear interpolation, and the sample standard deviation uses `n−1`. Values more than 1.5 interquartile ranges beyond the first or third quartile are flagged for review and kept in the summaries. In the normal-cognition group, this rule flagged four scores and fifteen thickness measurements. The output gives separate denominators for the paired samples and the modeling samples.

## 7. Phase 1 figures

| Figure | Sample and display | Descriptive reading |
|---|---|---|
| [Animal-naming score distribution](figures/phase1/01_animal_naming.pdf) | 843 usable scores in the normal-cognition group. Each bar represents one whole-number score. | Scores are concentrated around 22 correct responses. |
| [Left cortical thickness distribution](figures/phase1/02_left_parietal_thickness.pdf) | 847 usable measurements in the normal-cognition group, in bins 0.05 mm wide. | Values cluster around 2.3 mm, with a tail toward lower thicknesses. |
| [Cognitive status frequencies](figures/phase1/03_cognitive_status.pdf) | Counts and percentages for all 1,237 matched participants. | Normal cognition is the largest group. Seven participants are classified as impaired without MCI. |
| [Left cortical thickness by cognitive status](figures/phase1/04_thickness_by_status.pdf) | Box plots with the number of usable values below each group. | Median thickness is higher in the normal-cognition group than in the MCI and dementia groups. The distributions overlap. |

The line inside each box is the median, and the box edges are the first and third quartiles. The whiskers extend to the most extreme observed values within 1.5 interquartile ranges of the box. Circles mark values outside those limits, which are also included in the summaries. The comparisons describe the recorded groups.

## 8. Phase 1 work plan

| Step | Work | Completion criterion |
|---|---|---|
| 1. Set up the project and document the data | Establish a separate project folder and record the question, provenance, unit of observation, and limitations. | The four inputs match the supplied files, and the project card identifies the variables and their roles. |
| 2. Map and audit the data | Classify UDS and MRI columns and check special codes, ranges, fractional counts, missingness, and inconsistencies. | The checks can be reproduced and their rationale is documented. |
| 3. Link records and define the samples | Validate the ID-based join and define the primary normal-cognition sample, the full matched sample, and the denominator for each output. | Input counts, matches, exclusions, and their reasons can be regenerated. |
| 4. Summarize the data | Report frequencies, percentages, measures of center and spread, and descriptive group comparisons. | Tables state their denominators and explain the data quality and observed differences. |
| 5. Visualize the data | Plot the score and thickness distributions, cognitive status frequencies, and thickness by group. | At least four figures have readable titles, units, and captions. |
| 6. Assemble and review the outputs | Complete the notebook, summarize the research question and limitations, and build and review the LaTeX report. | Phase 1 can be reproduced from the supplied CSVs, and `reports/phase1.pdf` is ready. |

The final outputs for this phase are the [English report](reports/phase1.pdf) and [executed English notebook](notebooks/phase1.ipynb). The tables and figures were generated from this project's analysis.

Phase 2 examines distributions, conditional probabilities, and the unadjusted association. Phase 3 covers estimation, statistical tests, and participant-level bootstrap resampling. Phase 4 compares UDS-only, MRI-only, and combined models on a common sample. Phase 5 examines robustness and the limits of interpretation.

## References

- [Official statistics course project guide, pp. 4–8 and 18–19](docs/amarpr.pdf)
- [NACC Researcher's Guide](https://www.naccdata.org/the-nacc-researchers-guide/)
- [UDS data dictionary for ANIMALS and VEG](https://r2d2.kumc.edu/ADC/Protocols/UDS2_protocol/ResearchDataDict_2015Mar.pdf#page=205)
- UDS data dictionary for [sex](https://r2d2.kumc.edu/ADC/Protocols/UDS2_protocol/ResearchDataDict_2015Mar.pdf#page=41), [education](https://r2d2.kumc.edu/ADC/Protocols/UDS2_protocol/ResearchDataDict_2015Mar.pdf#page=45), [age at the visit](https://r2d2.kumc.edu/ADC/Protocols/UDS2_protocol/ResearchDataDict_2015Mar.pdf#page=290), and [cognitive status](https://r2d2.kumc.edu/ADC/Protocols/UDS2_protocol/ResearchDataDict_2015Mar.pdf#page=221)
- [MRI data dictionary for inferior parietal cortical thickness](https://files.alz.washington.edu/documentation/SCAN-MRI-Imaging-RDD.pdf#page=17)
- [Eastman et al., cortical thickness and semantic verbal fluency](https://pmc.ncbi.nlm.nih.gov/articles/PMC4207059/)
