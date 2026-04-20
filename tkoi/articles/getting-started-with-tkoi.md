# Transcriptomic Knowledge Graph Integration with tKOI

## Overview

`tKOI` (Transcriptomic Knowledge-graph Omics Integration) is an R
package for integrating transcriptomic data with a biological knowledge
graph to identify enriched biological concepts and pathways. It combines
PageRank propagation with permutation testing and ontology-based
annotations for high-resolution network interpretation.

## Step-by-Step Example

This vignette demonstrates a complete workflow using the `tkoi` package.

### Load Example Expression Data

``` r
file_path = system.file("extdata", "example_data.csv", package = "tkoi")
expression_data = data.table::fread(file_path)
head(expression_data)
```

### Run tKOI Analysis

``` r
tkoi_result = run_tkoi(
  expression_data = expression_data,
  subnetwork = tkoi::tkoi_net,
  pvalue_threshold = 0.05,
  logfc_threshold = 0.25,
  indirect_link_threshold = 3,
  topology_similarity = 0.9,
  n_permutation = 100,
  damping_factor = 0.85,
  maximum_iteration = 500
)
```

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

### Save Analysis Result (Optional)

``` r
save(tkoi_result, file = "tkoi_result.rda")
```

## Session Info

``` r
sessionInfo()
#> R version 4.5.3 (2026-03-11)
#> Platform: x86_64-pc-linux-gnu
#> Running under: Ubuntu 24.04.4 LTS
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
#>  [5] xfun_0.57         cachem_1.1.0      knitr_1.51        htmltools_0.5.9  
#>  [9] rmarkdown_2.31    lifecycle_1.0.5   cli_3.6.6         sass_0.4.10      
#> [13] pkgdown_2.2.0     textshaping_1.0.5 jquerylib_0.1.4   systemfonts_1.3.2
#> [17] compiler_4.5.3    tools_4.5.3       ragg_1.5.2        bslib_0.10.0     
#> [21] evaluate_1.0.5    yaml_2.3.12       otel_0.2.0        jsonlite_2.0.0   
#> [25] rlang_1.2.0       fs_2.1.0          htmlwidgets_1.6.4
```
