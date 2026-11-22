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

## References

- [Official statistics course project guide, pp. 4–8 and 18–19](docs/amarpr.pdf)
- [NACC Researcher's Guide](https://www.naccdata.org/the-nacc-researchers-guide/)
- [UDS data dictionary for ANIMALS and VEG](https://r2d2.kumc.edu/ADC/Protocols/UDS2_protocol/ResearchDataDict_2015Mar.pdf#page=205)
- UDS data dictionary for [sex](https://r2d2.kumc.edu/ADC/Protocols/UDS2_protocol/ResearchDataDict_2015Mar.pdf#page=41), [education](https://r2d2.kumc.edu/ADC/Protocols/UDS2_protocol/ResearchDataDict_2015Mar.pdf#page=45), [age at the visit](https://r2d2.kumc.edu/ADC/Protocols/UDS2_protocol/ResearchDataDict_2015Mar.pdf#page=290), and [cognitive status](https://r2d2.kumc.edu/ADC/Protocols/UDS2_protocol/ResearchDataDict_2015Mar.pdf#page=221)
- [MRI data dictionary for inferior parietal cortical thickness](https://files.alz.washington.edu/documentation/SCAN-MRI-Imaging-RDD.pdf#page=17)
- [Eastman et al., cortical thickness and semantic verbal fluency](https://pmc.ncbi.nlm.nih.gov/articles/PMC4207059/)
