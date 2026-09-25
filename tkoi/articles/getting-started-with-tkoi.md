# Transcriptomic Knowledge Graph Integration with tKOI

## Overview

`tKOI` (Transcriptomic Knowledge-graph Omics Integration) is an R
package for integrating transcriptomic data with a biological knowledge
graph to identify enriched biological concepts and pathways. It combines
PageRank propagation with permutation testing and ontology-based
annotations for high-resolution network interpretation.

Since version 1.1.0, the core of
[`run_tkoi()`](../reference/run_tkoi.md) is written in C++ and runs on
multiple threads. Personalized PageRank for the observed data and every
permutation is solved in batches, so a full analysis with 100
permutations on the complete knowledge graph (`tkoi_net`) took 25.5 to
38.1 s on an 8-core Apple M2, compared with 244.3 s for `tkoi` 1.0.0.
The package therefore needs a C++ compiler toolchain to install from
source (Xcode Command Line Tools on macOS, Rtools on Windows):

``` r

devtools::install_github("BaranziniLab/tkoi")
```

## Step-by-Step Example

This vignette demonstrates a complete workflow using the `tkoi` package.

### Load Example Expression Data

``` r

library(tkoi)

file_path = system.file("extdata", "example_data.csv", package = "tkoi")
expression_data = data.table::fread(file_path)
head(expression_data)
```

The table needs the columns `gene_name` (Ensembl gene IDs), `logfc`, and
`pvalue`. Rows with a missing or blank `gene_name` are ignored, and only
the first row of a duplicated gene is used. The downstream functions
([`export_gene_exploration_data()`](../reference/export_gene_exploration_data.md),
[`make_gene_exploration_plot()`](../reference/make_gene_exploration_plot.md),
and [`plot_network()`](../reference/plot_network.md)) read the table the
same way.

### Run tKOI Analysis

``` r

set.seed(1)

tkoi_result = run_tkoi(
  expression_data = expression_data,
  subnetwork = tkoi::tkoi_net,
  pvalue_threshold = 0.05,
  logfc_threshold = 0.25,
  indirect_link_threshold = 3,
  topology_similarity = 0.9,
  n_permutation = 100,
  damping_factor = 0.85,
  maximum_iteration = 500,
  n_cores = NULL,
  keep_permutations = TRUE
)
```

Genes are kept as seeds when `pvalue <= pvalue_threshold` and
`abs(logfc) >= logfc_threshold`, and each seed is weighted by
`abs(logfc)`. `indirect_link_threshold` does not filter anything: it
only orders each result table, listing first the nodes within two hops
of at least that many seed genes. No node is excluded.

The arguments added in version 1.1.0 control how the analysis runs and
what it stores; the statistics themselves do not depend on them beyond
numerical precision:

- `n_cores`: number of CPU threads. `NULL` (the default) uses every
  available core, or `getOption("tkoi.n_cores")` when that option is
  set. Requests above the number of available cores (including Slurm and
  Linux container CPU limits) are capped.
- `keep_permutations`: if `TRUE` (the default), the `pagerank_data` slot
  holds every permutation’s PageRank vector (`perm.1`, `perm.2`, …). If
  `FALSE`, it holds only the null mean and standard deviation
  (`null_mean`, `null_sd`), which saves memory for large
  `n_permutation`.
- `tolerance`: relative residual at which the PageRank solver stops
  (default `1e-14`). The solver needs about 55 iterations per PageRank
  vector on `tkoi_net`; `maximum_iteration` caps the number of
  iterations, and a warning is raised if a PageRank vector has not
  converged by then.
- `verbose`: set to `FALSE` to silence the progress messages.

### Performance and Reproducibility

The permutation null is drawn from R’s random number generator, so
[`set.seed()`](https://rdrr.io/r/base/Random.html) fixes the result, and
the same seed gives bit-identical statistics for any `n_cores`. On
`tkoi_net`, the same seed draws the same null gene sets as a sequential
`tkoi` 1.0.0 run, and the statistics agree with `tkoi` 1.0.0: `beta` to
within 1e-8 and PageRank vectors to within 6e-15. The result tables are
more complete than in 1.0.0: every node is reported (unannotated nodes
appear with missing annotation columns), the Compound FDR is adjusted
over the reported human metabolites only, nodes whose null has no
variation are untestable (`NaN`, listed last) instead of infinitely
significant, and the `human_metabolites` list was repaired. Rows can
otherwise be ordered differently only where FDR and `beta` agree to 10
significant digits, because version 1.1.0 orders such ties
deterministically.

The PageRank solver stops when both the residual of the symmetric system
it solves and the 1-norm residual of the original PageRank system are at
most `tolerance`. At the default `tolerance = 1e-14`, compared with a
reference solution computed at `1e-16` on `tkoi_net`, the 1-norm error
of a PageRank vector was 2.3e-15, and the 99th and 99.9th percentiles of
the per-node relative error were 1.5e-11 and 1.5e-10 (largest 7.5e-8).
[`igraph::page_rank()`](https://r.igraph.org/reference/page_rank.html)
gave 4.1e-12, 1.2e-9, 4.4e-9, and 2.7e-7 on the same measures, so the
default is at least as accurate both per node and overall.

On an Apple M2 (8 cores), a 21,414-gene dataset with 1,201 seed genes on
the full `tkoi_net` took 5.3 s with 10 permutations (48.5 s with `tkoi`
1.0.0) and 25.5 to 38.1 s with 100 permutations (244.3 s with `tkoi`
1.0.0), about 6 to 10 times faster; the range reflects run-to-run
variation from thermal throttling. The first run in a session also
builds the network matrix (8.8 s for 10 permutations). With
`n_cores = 1`, 10 permutations took 10.9 s, still faster than `tkoi`
1.0.0.

Loading [`tkoi::tkoi_net`](../reference/tkoi_net.md) takes a few seconds
and about 0.5 GB of memory the first time it is used in a session.
Keeping all 100 permutations adds about 0.75 GB. Before any PageRank is
computed, [`run_tkoi()`](../reference/run_tkoi.md) checks the memory the
run needs against the memory available (physical memory, or a Linux
container limit) and stops with advice, such as
`keep_permutations = FALSE`, if the run would not fit.

### Perform Gene Ontology Enrichment

``` r

tkoi_result = run_gene_enrichment(tkoi_result)
```

### Visualize GO vs Graph Enrichment

``` r

tkoi_result@gene_enrichment_comparison$comparison_scatter1
```

``` r

tkoi_result@gene_enrichment_comparison$comparison_scatter2
```

### Gene-Level Network Visualization

``` r

plt1 = make_gene_exploration_plot(
  tkoi_list = tkoi_result,
  sig_color = "#F39B7FB2",
  non_sig_color = "gray"
)
plt1
```

### Export Enrichment Summary Table

``` r

gene_data = export_gene_exploration_data(tkoi_result)
head(gene_data)
```

The node-level tables for every node type can be written to an Excel
workbook, one sheet per node type:

``` r

export_network_summary_statistics(tkoi_result, filename = "tkoi_network_statistics.xlsx")
```

### Visualize Top Enriched Genes

``` r

plt2 = visualize_topn(
  tkoi_list = tkoi_result,
  category = "Gene",
  top_n = 25,
  high_color = "#FF5733",
  low_color = "#154360"
)
plt2
```

### Plot the Network Around an Enriched Node

[`plot_network()`](../reference/plot_network.md) draws the part of the
knowledge graph that links a target node to the significant genes within
`degree_expansion` hops. Genes are colored by log fold change (blue
down, red up), the target is orange, and node size follows `beta`. Pass
the same `subnetwork` that was used in
[`run_tkoi()`](../reference/run_tkoi.md).

``` r

top_term = tkoi_result@network_summary_statistics$BiologicalProcess$node_id[1]

plot_network(
  tkoi_result = tkoi_result,
  target_node_id = top_term,
  degree_expansion = 2,
  network_layout_type = "kk",
  subnetwork = tkoi::tkoi_net
)

# Every node within one hop of the term, including the term itself
neighbors = get_neighboring_nodes(top_term, degree_expansion = 1)
```

### Save Analysis Result (Optional)

``` r

save(tkoi_result, file = "tkoi_result.rda")
```

## Session Info

``` r

sessionInfo()
#> R version 4.6.1 (2026-06-24)
#> Platform: x86_64-pc-linux-gnu
#> Running under: Ubuntu 24.04.5 LTS
#> 
#> Matrix products: default
#> BLAS:   /usr/lib/x86_64-linux-gnu/openblas-pthread/libblas.so.3 
#> LAPACK: /usr/lib/x86_64-linux-gnu/openblas-pthread/libopenblasp-r0.3.26.so;  LAPACK version 3.12.0
#> 
#> locale:
#>  [1] LC_CTYPE=C.UTF-8       LC_NUMERIC=C           LC_TIME=C.UTF-8       
#>  [4] LC_COLLATE=C.UTF-8     LC_MONETARY=C.UTF-8    LC_MESSAGES=C.UTF-8   
#>  [7] LC_PAPER=C.UTF-8       LC_NAME=C              LC_ADDRESS=C          
#> [10] LC_TELEPHONE=C         LC_MEASUREMENT=C.UTF-8 LC_IDENTIFICATION=C   
#> 
#> time zone: UTC
#> tzcode source: system (glibc)
#> 
#> attached base packages:
#> [1] stats     graphics  grDevices utils     datasets  methods   base     
#> 
#> loaded via a namespace (and not attached):
#>  [1] digest_0.6.39     desc_1.4.3        R6_2.6.1          fastmap_1.2.0    
#>  [5] xfun_0.61         cachem_1.1.0      knitr_1.52        htmltools_0.5.9  
#>  [9] rmarkdown_2.32    lifecycle_1.0.5   cli_3.6.6         sass_0.4.10      
#> [13] pkgdown_2.2.1     textshaping_1.0.5 jquerylib_0.1.4   systemfonts_1.3.2
#> [17] compiler_4.6.1    tools_4.6.1       ragg_1.5.2        bslib_0.12.0     
#> [21] evaluate_1.0.5    yaml_2.3.12       otel_0.2.0        jsonlite_2.0.0   
#> [25] rlang_1.3.0       fs_2.1.0          htmlwidgets_1.6.4
```
