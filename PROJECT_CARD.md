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

## 9. Phase 3 inference plan

The [analysis plan](reports/phase3/analysis_plan.json) is generated by [the Phase 3 planning script](src/phase3_plan.py). We set the plan after reviewing the Phase 1 and Phase 2 descriptions and before calculating the Phase 3 tests or confidence intervals. The analyses of this teaching extract are therefore exploratory.

We assume that participants are independent and target the source population represented by this extract, subject to each sample's eligibility rules. The undocumented sampling frame limits how far we can generalize, and interval coverage depends on the sampling assumptions. Each person contributes their single supplied UDS and MRI records. We retain every measurement flagged by the IQR review.

| Analysis | Sample and parameter | Planned method |
|---|---|---|
| Mean animal-naming score | Population mean ANIMALS score under the paired-normal eligibility rules, estimated from 843 participants | Arithmetic mean and an individual 95% Student t interval |
| Mean left cortical thickness | Population mean thickness under the same rules, estimated from the same 843 participants | Arithmetic mean and an individual 95% Student t interval |
| Score difference across MRI groups | Mean score at or above the fixed cutoff minus mean score below it, within normal cognition | Welch t test and a 95% Welch interval for the raw mean difference |
| Thickness across cognitive groups | The four mean thicknesses defined by UDS status, using 1,237 participants | Welch one-way ANOVA, group means with individual 95% t intervals, and descriptive sample eta-squared |

The MRI cutoff is fixed at **2.319 mm**, the pooled median used in Phase 2. Among the 843 paired participants with normal cognition, **486** are at or above the cutoff and **357** are below it. Values equal to the cutoff enter the higher group. The cutoff stays fixed during bootstrap resampling.

The first test has the null hypothesis that the score difference is zero and the two-sided alternative that it is nonzero. The second has the null hypothesis that all four mean thicknesses are equal and the alternative that at least one differs. The two planned tests form one family, with Holm-adjusted p-values evaluated at **0.05**. We will report both raw and adjusted p-values. The 95% confidence intervals are individual intervals.

In code order 1 through 4, the cognitive groups contain 847, 7, 252, and 131 participants. Alongside ANOVA, we will examine their sizes, distributions, spreads, and extreme values. We retain the seven-person group classified as impaired without MCI and make its limited precision explicit.

Both bootstrap analyses use **10,000 replicates** and percentile intervals with linear interpolation at the 2.5th and 97.5th percentiles. Whole-sample resampling uses seed **314159**, draws 843 complete participant rows with replacement, and computes both means from each draw. It illustrates sampling distributions and standard errors. The thickness mean's bootstrap interval will be compared with its t interval because Phase 2 showed a lower tail in the thickness measurements. Resampling for the score contrast uses seed **271828**, draws complete rows separately within each MRI group, and preserves the group sizes. Its interval will be compared with the Welch interval.

We will interpret interval widths and raw effects in correct responses and millimeters. Under the stated sampling assumptions, the nominal 95% interval procedure covers the population parameter in about 95% of repeated samples. The p-value describes how unusual the observed test statistic would be under the null and the test assumptions. Phase 4 will estimate the adjusted score difference per 0.1 mm of thickness using age, sex, and education.

### Phase 3, Step 2 estimates

The [mean estimates and confidence intervals](reports/phase3/mean_intervals.json) use the same 843 participants with normal cognition and usable values for both measures. We calculate the sample SD with denominator `n−1` and divide it by the square root of 843 to get the SE. For each individual 95% interval, we add and subtract the SE multiplied by 1.962785394, the Student t critical value with 842 degrees of freedom.

| Measure | Participants | Mean | Sample SD | Standard error | Individual 95% confidence interval |
|---|---:|---:|---:|---:|---|
| Animal-naming score, correct responses in 60 seconds | 843 | 22.3215 | 5.4350 | 0.18719 | 21.9541 to 22.6889 |
| Left inferior parietal cortical thickness, mm | 843 | 2.336446 | 0.121509 | 0.004185 | 2.328232 to 2.344660 |

These intervals describe our uncertainty about the two population means, subject to the eligibility rules and sampling assumptions above. The t procedure has exact coverage for independent normal observations; here we use it as a large-sample approximation for the means. The [completed bootstrap comparison in Step 3](reports/phase3/sampling_distributions.json) compares the thickness interval with a resampling-based interval given the observed lower tail.

Reproduce the results with `python3 -m src.phase3_mean_intervals` after installing `requirements.txt`. The [analysis script](src/phase3_mean_intervals.py) checks the saved plan and input data before calculating the estimates. It records the participant-table and plan hashes with the full-precision results. The Student t quantile comes from [SciPy](https://docs.scipy.org/doc/scipy/reference/generated/scipy.stats.t.html).

### Phase 3, Step 3 sampling distributions

The [bootstrap results](reports/phase3/sampling_distributions.json) contain 10,000 resamples of the same 843 paired participants with normal cognition, using seed 314159. In each draw, we select 843 complete participant rows with replacement and calculate both means from those rows. This keeps each person's UDS and MRI measurements together. We use the row order in the processed participant table and record its hash in the output.

| Measure | Participant sample SD | Analytic SE of the mean | Bootstrap SE of the mean |
|---|---:|---:|---:|
| Animal-naming score, correct responses | 5.43498 | 0.18719 | 0.18587 |
| Left inferior parietal thickness, mm | 0.121509 | 0.004185 | 0.004185 |

The participant SD shows how much measurements differ across people. The SE describes how the mean would vary if we repeatedly sampled 843 people. Here the bootstrap SE is the sample standard deviation of the 10,000 resampled means, using denominator `B−1`. It is 0.71% lower than the analytic SE for the score and 0.004% higher for thickness. The empirical distribution's exact with-replacement SE for a mean is the analytic SE multiplied by `sqrt((n−1)/n)`. Estimates from 10,000 replicates also have Monte Carlo variation. The mean of the resampled means estimates the observed mean; its small offset in this run reflects finite-replicate variation.

| Measure | Individual 95% Student t interval | Individual 95% bootstrap percentile interval |
|---|---|---|
| Animal-naming score, correct responses | 21.9541 to 22.6889 | 21.9561 to 22.6821 |
| Left inferior parietal thickness, mm | 2.328232 to 2.344660 | 2.328230 to 2.344604 |

Percentile limits use linear interpolation at the 2.5th and 97.5th percentiles. The thickness percentile interval is 0.33% narrower than its t interval, with closely aligned endpoints. Although individual thickness measurements have a lower tail, their resampled means are approximately bell-shaped at this sample size. Interval coverage remains approximate and depends on the sampling assumptions and eligibility rules above. Comparing the methods cannot resolve the missing information about the sampling frame.

The [score sampling distribution](figures/phase3/01_mean_score_sampling.pdf) and [thickness sampling distribution](figures/phase3/02_mean_thickness_sampling.pdf) are grayscale vector figures. Each histogram shows the density of the resampled means. The dashed curve is a normal reference centered at the observed mean, with SD equal to the analytic SE. The solid vertical line marks the observed mean, and the dotted lines mark the bootstrap percentile limits.

Reproduce this step with `python3 -m src.phase3_sampling`. The [script](src/phase3_sampling.py) validates the saved plan and mean intervals before resampling. The output stores all 10,000 pairs of means in aligned arrays, along with the input hashes, seed, generator and Python version. Resampling uses [Python's `random.choices`](https://docs.python.org/3/library/random.html#random.choices) with equal weights. The SE and percentile interval follow the [bootstrap definitions in SciPy's documentation](https://docs.scipy.org/doc/scipy/reference/generated/scipy.stats.bootstrap.html).

### Phase 3, Step 4 score contrast across MRI groups

The [Welch contrast](reports/phase3/score_difference_by_mri.json) compares animal-naming scores among the same 843 eligible participants with normal cognition. The fixed 2.319 mm thickness cutoff places 486 people in the at-or-above group and 357 in the below-cutoff group. Three participants are exactly at the cutoff and enter the higher group. Each person belongs to one group, so this is an independent-groups comparison.

| Left thickness group | Participants | Mean score | Sample SD |
|---|---:|---:|---:|
| At or above 2.319 mm | 486 | 22.4588 | 5.5432 |
| Below 2.319 mm | 357 | 22.1345 | 5.2861 |

The contrast is the higher group's mean minus the lower group's mean, in correct animal names in 60 seconds. The null hypothesis is a zero population mean difference; the alternative is a nonzero difference. The observed difference is **0.3244 correct responses**, with SE **0.37616** and an individual 95% Welch interval from **-0.4140 to 1.0628**. The two-sided test gives **t = 0.8624**, **df = 786.6565** and **raw p = 0.38874**.

If the population mean difference were zero and the Welch assumptions held, the approximate probability of an absolute t statistic at least as large as 0.8624 would be 0.38874. At the unadjusted 0.05 reference level, we therefore have insufficient evidence of a mean score difference. The interval includes effects ranging from about 0.414 fewer to 1.063 more correct responses in the higher group. Both the size and direction of the difference remain uncertain. To assess practical equivalence, we would need a prespecified margin, which this plan does not include.

Both groups have a median score of 22. The middle halves span 18.25 to 26.00 in the higher group and 18.00 to 26.00 in the lower group. Their sample SDs are 5.5432 and 5.2861. Separate group IQR fences flag four scores in the higher group and one in the lower group; all five remain in the analysis. Group-specific fences can flag a different number of values than the pooled review in Phase 1. The registered Welch method allows unequal group sizes and variances, with a large-sample approximation for the means.

This is an unadjusted exploratory contrast conditional on the fixed Phase 2 cutoff. Age, sex and education will enter the Phase 4 models. Interpretation also depends on the participant-independence assumption and the stated sampling-frame limitations. Holm adjustment belongs to the registered two-test family. The Step 4 JSON records the earlier stage, when the second test was unavailable, with `null` Holm fields. Final adjusted values and decisions are now available in the [completed two-test table](reports/phase3/hypothesis_tests.json). The 95% interval above is an individual interval.

Reproduce this step with `python3 -m src.phase3_mri_contrast`. The [script](src/phase3_mri_contrast.py) validates the saved plan, sample and group counts, and records the plan and participant-table hashes. The test uses [SciPy's independent-samples t procedure](https://docs.scipy.org/doc/scipy/reference/generated/scipy.stats.ttest_ind.html) with `equal_var=False`; the standard error and degrees of freedom follow the [Welch-Satterthwaite formulas documented by NIST](https://www.itl.nist.gov/div898/handbook/eda/section3/eda353.htm).

### Phase 3, Step 5 bootstrap comparison of the score contrast

The [contrast bootstrap](reports/phase3/score_difference_bootstrap.json) uses 10,000 replicates and seed 271828. Each replicate draws 486 complete participant rows with replacement from the at-or-above-cutoff group and 357 from the below-cutoff group. We resample the groups independently while keeping each person's UDS and MRI values together. Both group sizes and the 2.319 mm cutoff stay fixed. Within each replicate, we draw the higher group first and the lower group second, using one local random stream.

The observed higher-minus-lower mean difference remains **0.3244 correct responses**. The average resampled difference is 0.3243. The bootstrap SE is the sample standard deviation of the 10,000 differences, using denominator `B−1`. Percentile limits use linear interpolation at the 2.5th and 97.5th percentiles.

| Method | Standard error | Individual 95% confidence interval | Interval width |
|---|---:|---|---:|
| Welch | 0.37616 | -0.4140 to 1.0628 | 1.4768 |
| Within-group bootstrap, percentile | 0.37306 | -0.4119 to 1.0489 | 1.4608 |

All values in the table are in correct animal names in 60 seconds. The bootstrap SE is 0.82% lower than the Welch SE, and its interval is 1.09% narrower. The two methods give similar uncertainty estimates for this sample. Both intervals include zero and leave uncertainty about the size and direction of the population contrast.

For these empirical group distributions, the exact SE of the difference under with-replacement resampling is 0.37569. This calculation uses within-group variances with denominator `n`, while the Welch SE uses sample variances with denominator `n−1`. The 10,000-draw estimate also has Monte Carlo variation. For a difference of arithmetic means, the empirical-bootstrap expectation is the observed difference, so the small offset in this run comes from using a finite number of draws.

The resampling distribution is centered near the observed contrast and describes its uncertainty. The hypothesis-test p-value remains the Welch result from Step 4. Holm-adjusted values for both planned tests are now available in the [completed two-test table](reports/phase3/hypothesis_tests.json). These are individual confidence intervals conditional on the fixed cutoff and group sizes. Inference continues to depend on independent participants, the eligibility rules and the sampling-frame limitations. Agreement between the intervals does not establish population coverage or practical equivalence. Age, sex and education will be addressed in Phase 4.

The [comparison figure](figures/phase3/03_score_difference_bootstrap.pdf) shows the density of the resampled differences above the two 95% intervals. Its dashed curve is a normal reference centered at the observed difference, with SD equal to the Welch SE. The dotted vertical line marks zero. The point on each interval marks the same observed contrast.

Reproduce this step with `python3 -m src.phase3_contrast_bootstrap`. The [script](src/phase3_contrast_bootstrap.py) checks the saved Welch result against the current data and plan before resampling. The output retains all 10,000 differences, input and Welch-result hashes, seed, draw order, generator and Python version. Sampling uses [Python's `random.choices`](https://docs.python.org/3/library/random.html#random.choices) with equal weights; the percentile interval and bootstrap SE follow the [definitions documented by SciPy](https://docs.scipy.org/doc/scipy/reference/generated/scipy.stats.bootstrap.html).

### Phase 3, Step 6 cognitive-group comparison and completed tests

For the [thickness comparison](reports/phase3/thickness_by_status.json), we use all 1,237 participants with usable UDS cognitive status and left cortical thickness, including 847 with normal cognition. This MRI comparison does not require a usable animal-naming score. The four people excluded from the earlier 843-person paired-normal score analyses therefore remain eligible here.

| Cognitive status | Participants | Mean thickness (mm) | Sample SD (mm) | Individual 95% t interval (mm) |
|---|---:|---:|---:|---|
| Normal cognition | 847 | 2.3364 | 0.1213 | 2.3282 to 2.3446 |
| Impaired without MCI | 7 | 2.3143 | 0.1051 | 2.2171 to 2.4115 |
| MCI | 252 | 2.2576 | 0.1728 | 2.2361 to 2.2790 |
| Dementia | 131 | 2.1541 | 0.2143 | 2.1170 to 2.1911 |

The registered Welch one-way ANOVA tests the null hypothesis that all four population mean thicknesses are equal, against the alternative that at least one differs. It gives **F = 40.4282**, with **3** and **28.7651** degrees of freedom and **raw p = 1.914e-10**. Under the null and the Welch model assumptions, this p-value is the approximate upper-tail probability of an F statistic at least as large as the observed value.

Descriptive sample **eta-squared is 0.1439**, meaning that approximately **14.39%** of the observed squared thickness variation lies between the group means. It is calculated as the between-group sum of squares divided by the total sum of squares. This descriptive measure depends on the sample composition and is separate from the precision weights used in the Welch test. We use the group means and their intervals to show the effects in millimeters.

Group sizes differ by a factor of 121, and the largest sample variance is 4.159 times the smallest. Group-specific Q-Q plots show departures from a normal reference, including lower tails. Empirical skewness in code order is -0.950, -1.373, -1.426 and -0.556. The IQR review flags 15, 1, 7 and 2 measurements, respectively, and all are retained. Welch allows unequal variances, but its normal-theory approximation still depends on these distributions. The seven-person group classified as impaired without MCI has a wide mean interval and a low-tail observation. Its small size limits both precision and robustness to non-normality. We therefore interpret the omnibus result as an exploratory comparison under the recorded model assumptions.

The [completed two-test table](reports/phase3/hypothesis_tests.json) is the final source for raw and Holm-adjusted p-values at familywise alpha 0.05.

| Planned test | Raw p | Holm-adjusted p | Decision at 0.05 after Holm |
|---|---:|---:|---|
| Score difference across MRI groups | 0.38874 | 0.38874 | Do not reject the zero-difference null |
| Thickness across cognitive groups | 1.914e-10 | 3.828e-10 | Reject equality of all four means |

Holm orders the two raw p-values, applies the remaining-test multipliers and cumulative maxima, and caps adjusted values at one. It doubles the smaller ANOVA p-value in this family; the larger score-contrast p-value remains unchanged. Holm can control familywise error with dependent tests when their raw p-values are valid. The 95% mean intervals remain individual intervals. The omnibus rejection supports a difference in at least one population mean under the method's assumptions; identifying a specific pair would require a separately defined comparison. This analysis is unadjusted, with age, sex and education reserved for Phase 4. The sampling-frame limitations continue to restrict generalization.

The [mean-and-interval figure](figures/phase3/04_thickness_group_means.pdf) displays each group size beside its mean and individual 95% t interval. The [group Q-Q figure](figures/phase3/05_thickness_group_qq.pdf) shows all observations against a fitted normal reference within each group, including all seven points in the smallest group. Both figures are grayscale vector PDFs.

Reproduce the comparison results and figures with `python3 -m src.phase3_group_comparison`. The [script](src/phase3_group_comparison.py) validates the samples and the saved first test before computing the completed family. It uses the [Welch ANOVA formulation documented by SciPy](https://docs.scipy.org/doc/scipy/reference/generated/scipy.stats.f_oneway.html), with explicit degrees of freedom and the F survival probability. The [official R documentation](https://stat.ethz.ch/R-manual/R-devel/library/stats/html/p.adjust.html) describes Holm adjustment and its familywise error control. Input hashes identify the source data, plan and earlier contrast. The earlier Step 4 and Step 5 JSON outputs retain their stage-specific pending fields; the completed table above records the final adjustment.

### Phase 3, Step 7 integrated interpretation and review

The [inference review](reports/phase3/inference_review.json) passed 11 numerical and structural consistency checks. The participant table and sample flow were rebuilt from the supplied CSVs in memory. The saved plan and sample rules agree with the data. All numerical results, retained bootstrap draws and completed test-family results were checked, and all five figures reproduced byte for byte. These checks establish that the calculations are consistent. They cannot validate the sampling and measurement assumptions discussed below.

Within normal cognition, the higher-thickness group's observed mean score is 0.3244 correct responses higher than the lower group's. The Welch interval runs from -0.4140 to 1.0628, the bootstrap interval from -0.4119 to 1.0489, and the Holm-adjusted p-value is 0.38874. This grouped comparison gives insufficient evidence of a mean score difference at familywise alpha 0.05. Both intervals include zero and accommodate effects in either direction. A practical equivalence margin has not been specified, so the result leaves both the size and direction uncertain.

Across the four cognitive groups, Welch ANOVA gives F = 40.4282 with 3 and 28.7651 degrees of freedom and Holm-adjusted p = 3.828e-10. Under the method's assumptions, this supports a difference in at least one mean thickness. The descriptive sample eta-squared of 0.1439 partitions about 14.39% of the observed squared variation between group means. Individual group intervals describe the precision of each mean; identifying a specific pair requires a separately defined comparison. This supplementary analysis uses 1,237 people across cognitive statuses, while the score contrast uses 843 people with normal cognition. The analyses answer different questions and use different eligibility samples, so their findings need to be read in that context.

For the paired-normal thickness mean, the percentile interval is 0.33% narrower than the t interval. For the MRI-group score contrast, it is 1.09% narrower than the Welch interval. The procedures give similar uncertainty estimates conditional on the observed extract. Both bootstrap analyses reuse the supplied data. Their agreement is therefore an internal comparison of methods, with Monte Carlo variation from the finite 10,000 replicates. The two mean estimands have different units, correct responses and millimeters, and their SEs are interpreted on their own scales. All reported 95% intervals are individual intervals with nominal repeated-sampling coverage under the assumptions; Holm applies to the two planned test p-values.

The group of seven people classified as impaired without MCI has limited precision and a lower-tail observation. The group Q-Q plots show non-normal tails, so the Welch and t calibrations remain approximate, especially for that small group. All IQR-flagged measurements stay in the analysis. Using one supplied record per ID keeps the analysis at the participant level, while independence across people remains a model assumption. The student reports random selection, but its sampling frame is undocumented. Source preparation and imputation details and value-level masks are unavailable, so inference conditions on the supplied values and cannot quantify additional preparation uncertainty. Clinical/MRI timing and scan-level QC are also undocumented. These limitations affect the accuracy of inference and how far the findings can be generalized, so they remain part of our interpretation.

Having reviewed the extract descriptively, we use Phase 3 to estimate means and compare groups without covariate adjustment. These analyses are exploratory. The primary research parameter is the adjusted score difference per 0.1 mm greater thickness. Phase 4 will estimate that continuous association with age, sex and education on the common eligible sample of 842 participants with normal cognition. The current two-group comparison estimates a different parameter from that continuous association. We interpret p-values alongside effect sizes, units and interval widths, consistent with the [ASA principles on significance and p-values](https://www.amstat.org/asa/files/pdfs/P-ValueStatement.pdf). The p-values describe compatibility with the specified null models; their interpretation includes the modeling assumptions and the practical context.

Reproduce the review with `python3 -m src.phase3_review`. The [review script](src/phase3_review.py) records source-file hashes and environment versions and writes only its review output. The final p-values and decisions remain in the [completed two-test table](reports/phase3/hypothesis_tests.json). Step 8 assembles the shared tables, LaTeX report and notebook and documents the final requirement checks below.

### Phase 3, Step 8 assembled deliverables

The [report](reports/phase3.pdf) uses the [generated numeric source](reports/phase3/tables.tex) and the existing LaTeX template. It is grayscale. The [notebook](notebooks/phase3.ipynb) executes the calculations, shows the five shared figures and validates its results against the saved outputs. Section 10 of the report maps the nine general and five Track C requirements to the corresponding evidence.

The [table generator](src/phase3_tables.py) checks that the latest review passed, then formats each number. The integrated review was rerun after adding it, updating the source-code provenance while preserving the earlier numeric estimates and figure files. The [README](README.md#reproduce-phase-3) gives the full ordered analysis, LaTeX and notebook commands. The earlier step outputs remain reproducible snapshots, and the [completed two-test table](reports/phase3/hypothesis_tests.json) supplies final Holm values and decisions.

The next analysis phase estimates the continuous thickness association and compares the planned models on the common eligible sample.

## 10. Phase 4 model specification, Step 1

The [pre-fit analysis plan](reports/phase4/analysis_plan.json) is generated by [the planning script](src/phase4_plan.py). It records the question, eligibility rules, coding, three models and interpretation policy before any Phase 4 regression is fitted. Earlier phases have already described and compared these data, so the work remains exploratory. Recording this specification does not make it an external preregistration or an independent confirmation.

All three models use the same **842 participants with normal cognition** and usable animal-naming score, left thickness, age, recorded sex and education. The plan rebuilds the participant table from the raw CSVs in memory and checks it against the existing table and sample flow. It verifies the completed Phase 3 review and records input hashes, including the eligible-ID membership hash without exporting the IDs. Of 847 participants with normal cognition, four have unusable scores and one additional person has unknown education. The common sample contains 298 recorded male and 544 recorded female participants. Valid IQR-flagged measurements remain included.

The outcome is `ANIMALS`, the number of correct animal names in 60 seconds. Usable counts range from 0 to 77, and larger scores indicate better task performance. Special codes and fractional counts keep the Phase 1 exclusions. OLS will describe a linear conditional-mean association; it does not turn the count into a diagnosis or a clinically calibrated change measure.

To make intercepts meaningful, continuous predictors are centered at a jointly observed, central profile selected using predictors only: **recorded male, age 68 years, education 16 years and left thickness 2.344 mm**. Among male participants in the common sample, the script selects the profile closest to the sample's age, education and thickness medians, using the sum of squared deviations divided by the respective sample variances. Ties are resolved by ascending age, education and thickness. The selected person's score and ID are not recorded in the plan. This fixes an observed-support prediction point before fitting, without choosing it from favorable model results.

The coding is `age_c = NACCAGE - 68`, `education_c = EDUC - 16`, `female = 0` for recorded male (`SEX=1`) and `1` for recorded female (`SEX=2`), and `thickness_0_1mm = (LH_INFERIORPARIETAL_AVGTH - 2.344) / 0.1`. A thickness coefficient is therefore in correct responses per 0.1 mm, while age and education coefficients are per year. The sex coefficient compares recorded female with recorded male at the same values of other included predictors. We keep thickness continuous and do not reuse the Phase 3 median split.

| Model | Prespecified explanatory terms, with an intercept | Parameters; residual df if full rank |
|---|---|---|
| UDS-only | `age_c + female + education_c` | 4; 838 |
| MRI-only | `thickness_0_1mm` | 2; 840 |
| UDS + MRI | `thickness_0_1mm + age_c + female + education_c` | 5; 837 |

The primary estimand is the thickness coefficient in the combined model: the adjusted mean score difference per 0.1 mm greater thickness at the same age, recorded sex and education. Its prespecified test is two-sided, `H0: beta_thickness = 0` versus `H1: beta_thickness != 0`, at **alpha 0.05**, with an individual **95%** interval. Model-based OLS t inference is primary and depends on its stated error and mean-form assumptions. The simple-regression slope test and other coefficient p-values are supporting, nominal results. We do not claim familywise control across coefficients or phases; the completed Phase 3 Holm family is unchanged. Effects and interval widths take priority over significance labels, and a nonsignificant result will not be called equivalent without a justified margin.

The plan also fixes an **HC3 companion** for coefficient uncertainty in all three models, reported alongside the model-based results regardless of diagnostic significance. It uses the same OLS estimates and a heteroskedasticity-robust covariance, following the [official HC3 definition](https://www.statsmodels.org/stable/generated/statsmodels.regression.linear_model.OLSResults.HC3_se.html). Its t reference with residual degrees of freedom is an explicit approximation. We will not choose the method with the smaller p-value, count agreement as independent replication, or use HC3 to justify an individual prediction interval under changing residual variance.

Model comparisons will show the sample size, coefficients and uncertainty, R-squared, adjusted R-squared, SSE, residual SD and in-sample RMSE. Residual SD uses `sqrt(SSE/(n-p))`; RMSE uses `sqrt(SSE/n)`, with `p` including the intercept. These are in-sample association and fit summaries, not externally validated prediction performance. Pearson correlation and its comparison with simple regression will use these same 842 participants; the earlier 843-person correlation is not a substitute.

The planned checks include residual-versus-fitted, normal Q-Q and scale-location plots for all models. The combined model also receives predictor-specific residual checks, internally studentized residuals, leverage and Cook's distance review. The thresholds `|studentized residual| > 3`, `leverage > 2p/n` and `Cook's distance > 4/n` are review heuristics, not deletion rules. VIF above 5 prompts a documented collinearity review in the multiple regressions. All eligible participants stay in the primary fit; any leave-one-out influence checks are labeled diagnostic and summarized without participant IDs. A documented curvature pattern may justify one exploratory combined-model extension with squared centered thickness. Any other transformation needs an explicit diagnostic and substantive rationale before fitting. The original models remain visible.

At the fixed reference profile, we will report both the **95% interval for the mean score** and the **95% prediction interval for a new independent eligible participant**. The latter adds individual residual variation. Classical intervals require an adequate mean form and constant conditional error variance, with normal errors for exact t calibration; count-scale coverage is an approximation. Their limits will be reported without rounding to integer scores or clipping to 0-77, with any scale conflict noted. This is an educational model-based illustration, not a clinical prediction tool or an out-of-sample validation observation.

The single supplied record per person is verified, but the source visit count, clinical visit date and number, MRI date and interval between measurements remain unavailable. No first visit, baseline visit, nearest scan or temporal alignment is asserted. **Track C5 is only partly documentable:** participant-level structure and the supplied-record rule are documented, while temporal proximity would require authoritative metadata or clarification from the course team. Independence across people, the sampling frame, imputation/preparation uncertainty and scan quality remain limitations.

The plan maps all **11 general and 6 Track C requirements** to later work. These are planned coverage items, not completed checks. Step 2 will verify the common-sample design matrices; Step 3 covers simple regression and correlation; Step 4 fits the demographic and combined models; Step 5 compares models and collinearity; Step 6 examines residuals and influence; Step 7 documents justified extensions and profile intervals; Step 8 checks numerical consistency, interpretation and deliverables. No Phase 4 regression, coefficient test, prediction, report, notebook or figure has been produced in Step 1.

Run `python3 -m src.phase4_plan` after reproducing Phases 1-3. The command creates the plan if absent, leaves an identical existing plan untouched and refuses to overwrite a different one. Later analysis will use `load_registered_plan()` to check the frozen specification against its inputs. An intentional change requires an explicit reviewed plan revision before fitting. No dependency is added at this planning stage.

### Phase 4, Step 2 common-sample design matrices

The [design preparation report](reports/phase4/design_matrices.json) is generated by [src.phase4_design](src/phase4_design.py). It validates the frozen Step 1 plan, rebuilds eligibility from the raw CSVs in memory and prepares the three model matrices with the same **842 unique participants**, the same response vector and the same ascending-ID row order. The sample contains 298 recorded male and 544 recorded female participants. All response values are usable integer counts; their observed range in this sample is 7 to 43.

| Model | Rows × columns, including intercept | Full column rank | Residual df | Standardized-design condition number |
|---|---:|---:|---:|---:|
| UDS-only | 842 × 4 | 4 | 838 | 1.1289 |
| MRI-only | 842 × 2 | 2 | 840 | 1.0006 |
| UDS + MRI | 842 × 5 | 5 | 837 | 1.3508 |

The matrices preserve the registered column order, male-reference sex indicator, age and education centering, and continuous thickness per 0.1 mm. Decimal subtraction and scaling are applied before conversion to floating-point values. A reverse check recovers each supplied continuous predictor within an explicit floating-point tolerance and the sex code exactly. Shared columns are identical across models. The observed reference profile maps to the combined-model row `[1, 0, 0, 0, 0]` in its registered column order.

SVD and pivoted QR agree that all three matrices have full column rank. The tolerance and singular values are recorded, so a rank-deficient or numerically ambiguous design cannot silently lose a column. Conditioning is checked both in the registered units and after centering and dividing non-intercept predictors by their sample SDs. This standardization is only a numerical diagnostic: it does not replace the model matrices or change coefficient units. Full rank and moderate condition numbers establish numerical readiness, not valid regression assumptions or a substantive MRI association. VIF, residuals and influence remain later tasks.

The report contains aggregate checks and fingerprints only. Participant IDs, design rows and individual responses remain in memory and are not written to a new CSV or other participant-level output. After `python3 -m src.phase4_design`, later fitting code will call `load_prepared_designs()` to reconstruct the matrices and require the saved preparation report to match the current plan, source code, inputs and numerical environment. The command is idempotent for unchanged inputs. It uses the existing SciPy dependency; no requirement, frozen plan or previous-phase output was changed.

Step 2 prepares and validates inputs only. No regression coefficients, score-thickness correlation, p-values, residual diagnostics or prediction intervals have been calculated. The next step is the prespecified MRI-only simple regression and Pearson correlation on these same 842 participants.

### Phase 4, Step 3 simple regression and correlation

The [MRI-only results](reports/phase4/simple_regression.json) use the same **842 participants** and the validated two-column matrix from Step 2. The [shared OLS calculation](src/phase4_regression.py) solves a pivoted QR system and returns coefficients to the registered order. The [Step 3 analysis](src/phase4_simple.py) fits only MRI-only and calculates Pearson correlation on these same people. UDS-only and the integrated model are estimated in Step 4 below.

For centered thickness `T = (LH_INFERIORPARIETAL_AVGTH - 2.344) / 0.1`, the estimated line is `fitted_ANIMALS = 22.3376 + 0.0810 * T`. The intercept is the fitted mean score at **2.344 mm**, within the observed thickness range. Its individual model-based 95% interval is **21.9694 to 22.7059**. Because this model contains thickness alone, the intercept and slope do not condition on age, recorded sex or education.

| Uncertainty method | Slope per 0.1 mm | SE | Individual 95% slope interval | Two-sided nominal p |
|---|---:|---:|---|---:|
| Model-based OLS | 0.08098 | 0.15423 | -0.22173 to 0.38370 | 0.59965 |
| HC3 companion | 0.08098 | 0.15425 | -0.22178 to 0.38375 | 0.59971 |

The classical slope test gives **t = 0.5251** with **840 residual degrees of freedom**. Under zero population slope and the classical model assumptions, p = 0.59965 is the approximate probability of a t statistic at least as far from zero as observed. At nominal alpha 0.05, this secondary simple-model test does not reject zero slope. The intervals allow differences in either direction. We have not specified a practical-equivalence margin, so this result cannot establish equivalence or the absence of an association. Model-based uncertainty relies on the stated mean/error assumptions; the HC3 t reference is an explicit approximation and does not validate the functional form or sampling design.

Pearson **r = 0.01811** and **R-squared = 0.00032814**. About **0.0328%** of the observed squared score variation is accounted for by this unadjusted line. Adjusted R-squared is **-0.00086194**: its parameter penalty exceeds this small fit improvement, so the negative value is retained rather than clipped. Residual SD is **5.43297** correct responses, using `SSE/(n-p)`; in-sample RMSE is **5.42651**, using `SSE/n`. None of these measures establishes predictive performance on new data or a causal contribution of thickness.

The calculation checks `R_squared = r*r`, the simple-slope scaling identity and equality of the classical slope and Pearson nominal p-values. Those are linked identities for this one-predictor model, not independent scientific confirmations. The earlier Phase 2 correlation used 843 participants and is not substituted for this common-sample result. QR residual orthogonality, the sum-of-squares decomposition and the hat-matrix trace are also checked. HC3 is reported alongside model-based uncertainty as planned; the method with the smaller p-value is not selected.

The [grayscale scatterplot and fitted line](figures/phase4/01_score_thickness_simple_regression.pdf) show all 842 observations without jitter or exclusions. Duplicate coordinates can overlap. Thickness is displayed in millimeters, scores in correct responses, and the line is drawn over the observed thickness range. No mean-response or individual prediction ribbon is shown. The PDF is vector output, with participant IDs omitted.

Run `python3 -m src.phase4_simple` after the plan and design preparation. The result stores aggregate coefficients, covariance matrices, fit summaries and source hashes; fitted values, residuals, leverages and participant IDs stay in memory. Extremely small computed coefficient p-values carry an underflow flag rather than being interpreted as literal zero probabilities. `load_simple_results()` checks that the saved result still matches the current inputs and code.

This completes the unadjusted simple-regression and correlation calculation. The small observed point estimate and wide uncertainty relative to that estimate do not answer the adjusted research question. Step 4 below now adds UDS-only and the combined model. Residual/influence diagnostics, justified model extensions, profile prediction intervals and the assembled Phase 4 report remain later work. The frozen plan, preparation report, earlier phases and requirements are unchanged.

### Phase 4, Step 4 multiple regression and adjusted thickness association

The [multiple-regression analysis](src/phase4_multiple.py) fits the registered UDS-only and integrated models through the unchanged shared OLS engine. It first requires the frozen plan, current design report and reproducible MRI-only result. Both models retain all **842 participants**, identical response values and row order. No predictor is added, removed or selected by its p-value. The [aggregate results](reports/phase4/multiple_regression.json) contain every coefficient, both covariance matrices, SEs, t statistics, residual degrees of freedom, two-sided p-values, individual 95% intervals and fit summaries.

Coding remains `A = NACCAGE - 68`, `E = EDUC - 16`, `F = 0` for recorded male and `F = 1` for recorded female, and `T = (LH_INFERIORPARIETAL_AVGTH - 2.344) / 0.1`. The fitted mean-score equations, rounded here for readability, are:

- UDS-only: `fitted_ANIMALS = 22.11592 - 0.13790*A - 0.19393*F + 0.60271*E`.
- Integrated: `fitted_ANIMALS = 22.06945 - 0.12570*T - 0.14229*A - 0.14241*F + 0.60324*E`.

**UDS-only coefficient estimates** use **838 residual degrees of freedom**. All p-values in this table are nominal supporting descriptions.

| Term | Estimate | Model-based SE | Model-based 95% interval | Two-sided p | HC3 SE | HC3 95% interval | Nominal HC3 p |
|---|---:|---:|---|---:|---:|---|---:|
| `intercept` | 22.11592 | 0.29912 | 21.52880 to 22.70304 | Numerical underflow | 0.29057 | 21.54559 to 22.68626 | Numerical underflow |
| `age_c` | -0.13790 | 0.02216 | -0.18139 to -0.09441 | 7.66719e-10 | 0.02377 | -0.18455 to -0.09125 | 9.28627e-09 |
| `female` | -0.19393 | 0.36930 | -0.91879 to 0.53092 | 0.59963 | 0.36711 | -0.91450 to 0.52663 | 0.59745 |
| `education_c` | 0.60271 | 0.06829 | 0.46868 to 0.73675 | 6.23161e-18 | 0.06375 | 0.47759 to 0.72784 | 3.15708e-20 |

**Integrated coefficient estimates** use **837 residual degrees of freedom**. The model-based thickness test is the single registered primary test; other coefficient tests and all HC3 companion tests are nominal supporting descriptions. All intervals in both tables are individual rather than simultaneous.

| Term | Estimate | Model-based SE | Model-based 95% interval | Two-sided p | HC3 SE | HC3 95% interval | Nominal HC3 p |
|---|---:|---:|---|---:|---:|---|---:|
| `intercept` | 22.06945 | 0.30433 | 21.47212 to 22.66679 | Numerical underflow | 0.29532 | 21.48980 to 22.64911 | Numerical underflow |
| `thickness_0_1mm` | -0.12570 | 0.15078 | -0.42165 to 0.17024 | 0.40468 | 0.14523 | -0.41077 to 0.15936 | 0.38699 |
| `age_c` | -0.14229 | 0.02278 | -0.18699 to -0.09758 | 6.65644e-10 | 0.02480 | -0.19097 to -0.09361 | 1.34526e-08 |
| `female` | -0.14241 | 0.37450 | -0.87748 to 0.59265 | 0.70384 | 0.36898 | -0.86664 to 0.58182 | 0.69962 |
| `education_c` | 0.60324 | 0.06830 | 0.46918 to 0.73731 | 5.95577e-18 | 0.06386 | 0.47790 to 0.72859 | 3.40506e-20 |

For the primary adjusted association, the estimated score difference is **-0.12570 correct responses per 0.1 mm greater thickness**, at the same age, recorded sex and education. The model-based test gives **t = -0.83372**, **p = 0.40468** and a **95% interval from -0.42165 to 0.17024**. Its null is a zero integrated-model thickness coefficient; the alternative is two-sided. At the registered alpha of 0.05, the primary null is not rejected. The HC3 companion gives **p = 0.38699** and an interval from **-0.41077 to 0.15936**, using the prespecified approximate t reference with 837 degrees of freedom. Both intervals include zero and allow differences in either direction. No practical-equivalence margin has been supplied, so non-rejection does not establish no association or equivalence. These estimates are conditional associations, not within-person changes or causal effects. Residual and influence diagnostics are still pending.

The intercept describes a fitted mean for recorded males aged 68 with 16 years of education; the integrated model also fixes thickness at 2.344 mm. In the integrated model, one year older age corresponds to a fitted mean difference of **-0.14229** at the same thickness, recorded sex and education. One additional year of education corresponds to **0.60324** at the same thickness, age and recorded sex. The `female` coefficient, **-0.14241**, is recorded female minus recorded male at the same thickness, age and education. UDS-only coefficients have the same units and reference coding, with adjustment for the other UDS predictors and no thickness term. These are between-person comparisons over the observed support; the age coefficient is not a longitudinal aging rate, and the education coefficient is not a treatment effect. We report their uncertainty without making a collection of unadjusted significance claims.

UDS-only has **SSE = 21686.35393**, **R-squared = 0.12564023**, **adjusted R-squared = 0.12251006**, **residual SD = 5.08711** and **in-sample RMSE = 5.07501**. The integrated model has **SSE = 21668.35941**, **R-squared = 0.12636574**, **adjusted R-squared = 0.12219067**, **residual SD = 5.08804** and **in-sample RMSE = 5.07291**. Residual SD uses `sqrt(SSE/(n-p))`; RMSE uses `sqrt(SSE/n)`, with the intercept included in `p`. These are in-sample summaries. Step 5 below now assembles the three-model comparison and reviews coefficient changes and collinearity; no best-model or external prediction claim is made here.

The implementation checks the QR estimates and SSE against an independent SVD-based LAPACK `gelsd` solve for each multiple model. Earlier source fingerprints remain required. All fit arrays, participant IDs and responses stay in memory; the public result contains aggregates only. A computed zero p-value for the intercept carries a numerical-underflow flag and is not a literal zero probability. `load_multiple_results()` rebuilds and checks the saved result before later steps can use it. To reproduce this step after the preceding commands, run:

```shell
python3 -m src.phase4_multiple
```

The registered model-based primary test is retained regardless of the HC3 p-value. Classical t inference depends on the stated mean/error assumptions; the HC3 t convention is approximate. Neither method resolves unknown source preparation or imputation uncertainty, participant dependence, selection, unverified UDS/MRI timing or an inadequate mean specification. This step completes fitting and coefficient interpretation. Step 5 below now completes the collinearity review. Residual/influence review, any justified extension, profile intervals, the executed notebook and the assembled LaTeX report remain later work. The frozen plan, earlier analytical artifacts and requirements are unchanged.

### Phase 4, Step 5 model comparison and collinearity

The [comparison code](src/phase4_comparison.py) requires reproducible results from Steps 1-4 before producing the [three-model comparison](reports/phase4/model_comparison.json). The comparison uses the same **842 participants**, response values, shared predictor values and row order. It reuses the validated coefficient/fit results and checks that both baseline designs are nested in the integrated design. No participant or covariate is removed.

| Model | n | p | Residual df | SSE | R-squared | Adjusted R-squared | Residual SD | In-sample RMSE |
|---|---:|---:|---:|---:|---:|---:|---:|---:|
| UDS-only | 842 | 4 | 838 | 21686.35393 | 0.12564023 | 0.12251006 | 5.08711 | 5.07501 |
| MRI-only | 842 | 2 | 840 | 24794.41346 | 0.00032814 | -0.00086194 | 5.43297 | 5.42651 |
| UDS + MRI | 842 | 5 | 837 | 21668.35941 | 0.12636574 | 0.12219067 | 5.08804 | 5.07291 |

Here `p` includes the intercept. Residual SD uses `sqrt(SSE/(n-p))`, while in-sample RMSE uses `sqrt(SSE/n)`. All models have the same `SST = 24802.55226`. Adding thickness to UDS-only reduces SSE by **17.99453** and raises R-squared by **0.00072551**, about **0.07255 percentage points** of observed squared score variation. Adjusted R-squared decreases by **0.00031940**, and residual SD increases by **0.00092574**; RMSE decreases by **0.00210597**. This is numerically coherent: adding a column cannot increase nested-model SSE, but the residual-variance estimate uses one fewer residual degree of freedom. Compared with MRI-only, adding the demographic terms reduces SSE by **3126.05405** and raises R-squared by **0.12603760**. These are descriptive in-sample changes. They do not identify a best predictive model or establish external prediction performance.

| Term | Baseline | Baseline estimate | Integrated estimate | Integrated minus baseline | Model-based SE ratio | HC3 SE ratio |
|---|---|---:|---:|---:|---:|---:|
| `thickness_0_1mm` | MRI-only | 0.08098 | -0.12570 | -0.20669 | 0.97762 | 0.94153 |
| `age_c` | UDS-only | -0.13790 | -0.14229 | -0.00438 | 1.02795 | 1.04346 |
| `female` | UDS-only | -0.19393 | -0.14241 | 0.05152 | 1.01408 | 1.00509 |
| `education_c` | UDS-only | 0.60271 | 0.60324 | 0.00053 | 1.00023 | 1.00178 |

SE ratios above are integrated divided by baseline. The [grayscale coefficient figure](figures/phase4/02_coefficient_comparison.pdf) displays both individual 95% intervals for each estimate; full interval endpoints, widths and width ratios are also saved in the JSON. Thickness is per 0.1 mm, age and education per year, and recorded sex is female minus male. Panels use their own x-axis ranges. Intercepts are not treated as a common slope comparison because each model conditions on a different predictor set. The coefficients use the same participants and are dependent. A change in sign, interval overlap or a changed p-value is not an independent test that the coefficients differ. We do not express the change as a percentage of a near-zero unadjusted slope.

The thickness sign change can be checked directly. The marginal slope equals the integrated thickness coefficient plus the sum of each integrated covariate coefficient times that covariate's slope on thickness alone. On the common sample this gives `0.08098479 = -0.12570421 + 0.22106636 - 0.00978953 - 0.00458784`, up to displayed rounding. The contributions shown are age, recorded sex and education, respectively. Thickness and age have correlation **-0.23787**; the negative integrated age coefficient contributes positively to the marginal thickness slope in this algebraic identity. This explains the observed numerical change without attributing it to a calculation error. It is not a causal decomposition or proof of confounding, and the two thickness intervals include zero.

**Predictor dependence and VIF.** Pearson correlations below include the recorded-sex indicator and exclude the response and intercept. Both predictor-correlation matrices are saved in the JSON; the UDS-only matrix is the corresponding demographic submatrix on these same observations.

| Predictor | `thickness_0_1mm` | `age_c` | `female` | `education_c` |
|---|---:|---:|---:|---:|
| `thickness_0_1mm` | 1.00000 | -0.23787 | 0.17452 | -0.00358 |
| `age_c` | -0.23787 | 1.00000 | -0.06146 | -0.01728 |
| `female` | 0.17452 | -0.06146 | 1.00000 | -0.10212 |
| `education_c` | -0.00358 | -0.01728 | -0.10212 | 1.00000 |

| Predictor | UDS-only VIF | Integrated VIF |
|---|---:|---:|
| `thickness_0_1mm` | Not included | 1.08971 |
| `age_c` | 1.00436 | 1.06090 |
| `female` | 1.01464 | 1.04304 |
| `education_c` | 1.01111 | 1.01120 |

For each non-intercept predictor, an auxiliary SVD-based OLS regression includes an intercept and all the other predictors in that model. `VIF = SST_j/SSE_j = 1/(1-R_j²)`. Each result agrees with `SST_j * diagonal_j((X'X)^(-1))` from the earlier QR solution. The largest integrated VIF is **1.08971**, for thickness; the largest UDS-only VIF is **1.01464**, for recorded sex. No term exceeds the frozen **VIF > 5** review heuristic. The matrices retain full column rank, **4** for UDS-only and **5** for the integrated model. Their standardized-design condition numbers are **1.12887** and **1.35077**, respectively. Raw-unit condition numbers are recorded separately; diagnostic standardization does not change the reported coefficient units. MRI-only has a single explanatory term, so dependence among multiple MRI predictors is inapplicable.

These diagnostics show no indication of severe linear collinearity under the registered review rules in this sample. VIF measures the classical coefficient-variance factor associated with predictor dependence at fixed residual variance and predictor spread; its square root is the corresponding SE factor. It is not the HC3/model-based SE ratio, and low VIF does not validate residual assumptions, guarantee a narrow thickness interval or establish an absent association. The registered covariates and the primary integrated model are retained. No additional hypothesis test, automatic selection or alternative primary test is introduced.

After the preceding commands, run `python3 -m src.phase4_comparison`. `load_comparison_results()` rebuilds the comparison and rejects stale or edited results. The output stores aggregate summaries, small correlation matrices and source hashes; individual design, response, residual and influence arrays are not exported. The primary test and both uncertainty methods remain as reported in Step 4. Residual and influence review is now completed in Step 6 below; model retention or a justified extension, profile intervals and final report/notebook assembly remain later work. The frozen plan, earlier analytical artifacts and dependencies are unchanged.

### Phase 4, Step 6 residual and influence review

The [diagnostic calculations](src/phase4_diagnostics.py) rebuild the full-sample fits and require exact agreement with the saved coefficients, both uncertainty methods and fit summaries. [Figure code](src/phase4_diagnostic_figures.py) uses the resulting private arrays; the [aggregate diagnostic report](reports/phase4/residual_diagnostics.json) stores counts, distribution summaries, descriptive bins and sensitivity ranges. The registered primary models still contain all **842 participants**.

Three grayscale vector figures cover the planned review:

- [Residual, normal Q-Q and scale-location panels for all three models](figures/phase4/03_residual_diagnostics.pdf).
- [Integrated residuals versus thickness, age and education, with recorded-sex residual boxplots](figures/phase4/04_integrated_predictor_residuals.pdf).
- [Integrated residual/leverage and Cook/leverage panels, plus the flagged single-omission coefficient histogram](figures/phase4/05_integrated_influence.pdf).

Residuals are `e_i = y_i - fitted_i`. Internally studentized residuals are `r_i = e_i/(s*sqrt(1-h_ii))`, where `s = sqrt(SSE/(n-p))` comes from the full model. They are not externally studentized residuals. Cook's distance is `D_i = (r_i²/p)*h_ii/(1-h_ii)`, with `p` including the intercept. The registered strict review flags are `abs(r_i) > 3`, `h_ii > 2p/n` and `D_i > 4/n`. Cook's distance concerns the full coefficient vector, not just thickness. These are review heuristics, not tests or evidence that a record is erroneous.

| Model | Residual flags | Leverage flags | Cook flags | Union | Max abs(r) | Max leverage | Max Cook |
|---|---:|---:|---:|---:|---:|---:|---:|
| UDS-only | 5 | 39 | 34 | 65 | 4.27812 | 0.04362 | 0.03191 |
| MRI-only | 4 | 57 | 43 | 81 | 3.79262 | 0.06110 | 0.03977 |
| UDS + MRI | 5 | 45 | 37 | 70 | 4.29831 | 0.06470 | 0.02897 |

The integrated leverage threshold is **0.01187648**, and the Cook threshold is **0.00475059**. The three integrated counts overlap: their union is **70**, not their sum. The JSON records all eight mutually exclusive flag combinations, including the unflagged group. The influence scatterplots show all observations without IDs or source-row labels.

**Graphical evidence.** Normal Q-Q panels use standard-normal quantiles at `(i-0.5)/n` and a reference through the observed 25th/75th percentiles. Central points are near that line, while upper positive tails rise above it and lower tails flatten. Positive skewness and excess kurtosis support qualifying exact normal-error inference. The residual bands in the scatterplots partly reflect the integer-valued task score; they are not a reason to jitter, transform or discard scores.

| Model | Fitted-score range | Residual skewness | Excess kurtosis | Within-bin residual SD range |
|---|---|---:|---:|---|
| UDS-only | 13.29501 to 27.95405 | 0.37679 | 0.41938 | 4.64338 to 5.55121 |
| MRI-only | 21.63305 to 22.59839 | 0.36046 | 0.29797 | 4.99373 to 5.79958 |
| UDS + MRI | 13.28766 to 27.99339 | 0.38285 | 0.43383 | 4.82353 to 5.60290 |

Residual-versus-fitted and continuous-predictor panels retain every point and connect means from eight requested empirical-quantile bins. Duplicate edges are removed and tied predictor values stay together, so a panel can have fewer bins. These lines are descriptive, without confidence ribbons. In the integrated model, fitted-value bin residual means range from **-0.51646 to 0.56164**, small relative to residual SD **5.08804**. The central predictor support shows no obvious sustained broad curvature; sparse tails remain less informative. Scale-location bin means and within-bin SDs vary, with somewhat greater spread at the upper fitted end, but no dramatic funnel dominates. This does not establish linearity or constant variance. All fitted scores are within the usable 0-77 range; the fixed-profile prediction intervals are still reserved for Step 7. The planned HC3 companion remains visible alongside model-based uncertainty, without choosing it by a diagnostic p-value. No new normality or variance test is introduced.

**Diagnostic influence follow-up.** Each of the 70 flagged observations is omitted in turn, giving 70 separate fits with **841 participants** and **836 residual degrees of freedom**. Direct QR estimates are checked against `beta_(-i) = beta - (X'X)^(-1) x_i e_i/(1-h_ii)`. The thickness coefficients range from **-0.15084 to -0.08599**; the largest absolute change from the full-sample coefficient is **0.03971** correct responses per 0.1 mm. All 70 estimates remain negative. Both methods' recomputed intervals include zero in all 70 fits. The JSON also gives SE/endpoint ranges and summaries by each original flag category, without per-person rows.

| Fit | n | Residual df | Thickness coefficient | Model-based SE | Model-based interval | HC3 SE | HC3 interval |
|---|---:|---:|---:|---:|---|---:|---|
| Primary, full sample | 842 | 837 | -0.12570 | 0.15078 | -0.42165 to 0.17024 | 0.14523 | -0.41077 to 0.15936 |
| Diagnostic joint omission | 772 | 767 | -0.05327 | 0.16585 | -0.37884 to 0.27231 | 0.16373 | -0.37468 to 0.26815 |

The joint omission temporarily removes the **70-observation union**, leaving **772** participants and **767** residual degrees of freedom. Its coefficient shifts toward zero, to **-0.05327**, and its uncertainty remains substantial. The full 842-person primary coefficient stays **-0.12570**. The flagged set is fixed from the full model; it is not repeatedly recalculated to trim the data. Both methods use the usual individual 95% formulas in these diagnostic fits, but the flags depend on fitted outcomes. Their post-selection coverage is not justified, and ranges of interval endpoints are not simultaneous confidence bands. We report no new confirmatory p-values and do not adopt an omitted-observation model as the main analysis. These checks describe sensitivity, not equivalence, clinical importance or proof of an absent association.

Run `python3 -m src.phase4_diagnostics` after the preceding steps. `load_diagnostic_results()` rebuilds and checks the aggregate result before later work. Numerical validation covers the original fits, diagnostic formulas, flag partitions, bin membership, analytic deletion identity and source fingerprints. The figures were visually reviewed for this frozen sample on **2026-11-26**. The graphical assessment qualifies the classical error assumptions; it cannot verify the sampling frame, participant independence, source preparation/imputation, MRI quality or UDS/MRI timing. Step 7 will record the justified model-retention/extension decision before profile intervals or any optional extension. No quadratic term, response transformation or replacement primary model is fitted in this step. The frozen plan, earlier analytical artifacts and dependencies remain unchanged.

### Phase 4, Step 7 model decision and profile intervals

The residual plots support keeping the registered additive linear models for this analysis. In the integrated fit, mean residuals in the fitted-score bins range from -0.516 to 0.562 correct responses, against a residual SD of 5.088. The continuous-predictor plots show no clear, sustained curvature in the central range that would justify a quadratic thickness term. Sparse tails limit that assessment. Residual spread varies across bins and the normal Q-Q plot shows an asymmetric upper tail. This evidence qualifies the classical intervals and gives little guidance for choosing a specific score transformation. We retain the original score scale and report the prespecified HC3 companion for coefficient uncertainty. The full 842-person fit supplies the profile calculation, with the omission fits documenting sensitivity to influential observations. [Step 6 diagnostic results](reports/phase4/residual_diagnostics.json)

The profile fixed before fitting is a recorded male aged **68 years**, with **16 years of education** and **2.344 mm** of left inferior parietal thickness. Its values occur together in the modeling sample. The selection used predictor values, and the published profile contains those values alone. In registered column order, the profile row is `[1, 0, 0, 0, 0]`. The fitted conditional mean therefore equals the integrated intercept, **22.06945 correct animal names in 60 seconds**. [Frozen analysis plan](reports/phase4/analysis_plan.json)

| Nominal 95% interval | Standard error or prediction SD | Limits | Width |
|---|---:|---|---:|
| Conditional mean | 0.30433 | 21.47212 to 22.66679 | 1.19466 |
| One new independent observation | 5.09713 | 12.06480 to 32.07411 | 20.00932 |

With `g = x0' (X'X)^(-1) x0` and `s² = SSE/(n-p)`, the mean variance estimate is `s² g`. The prediction error variance estimate is `s² (1+g)` because an individual score adds residual variation to the uncertainty in the fitted mean. Both calculations use 837 residual degrees of freedom and a critical t value of 1.96280. The mean interval agrees with the classical intercept interval, and the individual interval is about 16.75 times as wide. Both limits lie in the usable score range of 0 to 77. The output keeps the continuous-scale limits at full numerical precision. [Prediction interval method](https://www.itl.nist.gov/div898/handbook/pmd/section5/pmd512.htm)

Exact classical calibration depends on the linear mean and independent normal errors with constant conditional variance. The prediction target also assumes a new eligible observation independent of the fitting sample and from the same conditional distribution. Given the observed residual tails and spread, interpret the numerical limits as model-based approximations, particularly for the individual score. Predictive coverage in new data remains unverified. HC3 addresses coefficient uncertainty and the individual interval uses the registered pooled residual variance.

Run `python3 -m src.phase4_prediction` from the project root after Step 6. The result records the evidence for retaining the model form (general requirement 8) and gives both intervals for the same profile (general requirement 9). Step 8 assembles and verifies the LaTeX report and notebook. [Code](src/phase4_prediction.py) · [Profile results](reports/phase4/profile_prediction.json)

### Phase 4, Step 8 report, notebook and verification

The final Phase 4 deliverables present the same 842-person analysis. The numerical review reproduces Steps 1–7, checks the three OLS fits and both covariance methods against an independent SVD solution, and reproduces the five figure PDFs. The executed notebook displays aggregate outputs and the report uses the generated numeric tables. [Final review](reports/phase4/regression_review.json) · [Report](reports/phase4.pdf) · [Notebook](notebooks/phase4.ipynb)

Section 7 of the report maps the 11 general and six Track C requirements to the relevant sections. The one-record-pair part of C5 is documented. Its temporal-proximity part requires dates or authoritative visit and scan metadata beyond the supplied files. The first seven steps retain their original analytical records.

## References

- [Official statistics course project guide, pp. 4–8 and 18–19](docs/amarpr.pdf)
- [NACC Researcher's Guide](https://www.naccdata.org/the-nacc-researchers-guide/)
- [UDS data dictionary for ANIMALS and VEG](https://r2d2.kumc.edu/ADC/Protocols/UDS2_protocol/ResearchDataDict_2015Mar.pdf#page=205)
- UDS data dictionary for [sex](https://r2d2.kumc.edu/ADC/Protocols/UDS2_protocol/ResearchDataDict_2015Mar.pdf#page=41), [education](https://r2d2.kumc.edu/ADC/Protocols/UDS2_protocol/ResearchDataDict_2015Mar.pdf#page=45), [age at the visit](https://r2d2.kumc.edu/ADC/Protocols/UDS2_protocol/ResearchDataDict_2015Mar.pdf#page=290), and [cognitive status](https://r2d2.kumc.edu/ADC/Protocols/UDS2_protocol/ResearchDataDict_2015Mar.pdf#page=221)
- [MRI data dictionary for inferior parietal cortical thickness](https://files.alz.washington.edu/documentation/SCAN-MRI-Imaging-RDD.pdf#page=17)
- [Eastman et al., cortical thickness and semantic verbal fluency](https://pmc.ncbi.nlm.nih.gov/articles/PMC4207059/)
