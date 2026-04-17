# Remaining Work

The manuscript can now be supported by the figures and tables generated from the existing `Analyses/full_run` outputs, but several project-level items still need additional analysis or publication metadata before the manuscript is fully final.

## Additional analyses that would strengthen the manuscript

1. Host-response analysis:
Differential-expression, pathway-enrichment, and interferon/innate-immunity summaries were not present in the finalized project outputs. These can likely be built from the `Analyses/full_run/tables/04_host_counts_*` matrices, ideally with an explicit design formula that accounts for protocol and inoculum group.

2. Viral mutation and evolution analysis:
Variant-calling, consensus, and phylogenetic outputs suitable for manuscript reporting were not present in the repository snapshot reviewed here. Those analyses would require defining the viral BAM inputs to use, coverage thresholds, allele-frequency filters, and the expected set of viruses/samples for follow-up.

3. Evidence-rich coverage examples:
If publication space allows, it would be valuable to add one or two direct genome-coverage plots from representative BAM files to complement the summary-style Figure 4. The current manuscript uses ViraQuant breadth metrics as a proxy because ready-made coverage panels were not available in the finalized outputs.

## Publication metadata still to finalize

1. Add the public GitHub or internal repository URL, license, and release tag in the Code Availability section.

2. Replace PMID placeholders with a formatted reference list in the `References` section.

3. Confirm whether the targeted reference bundle should be described in the text as a legacy `7viruses` label containing eight virus groups, or whether the bundle name itself should be updated before submission.
