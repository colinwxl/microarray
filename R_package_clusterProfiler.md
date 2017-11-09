## [clusterProfiler](https://bioconductor.org/packages/release/bioc/html/clusterProfiler.html) [Bioconductor Package]
2017-6-19
[Home Page](https://guangchuangyu.github.io/clusterProfiler/)
### Abstract
clusterProfiler implements methods to analyze and visualize functional profiles of genomic coordinates (supported by ChIPseeker), gene and gene clusters.

### Supported Analysis
- Over-Representation Analysis
- Gene Set Enrichment Analysis
- Biological theme comparison

### Supported ontologies/pathways
- Disease Ontology (via DOSE)
- Network of Cancer Gene (via DOSE)
- DisGeNET (via DOSE)
- Gene Ontology (supports many species with GO annotation query online via AnnotationHub)
- KEGG Pathway and Module with latest online data (supports more than 4000 species listed in http://www.genome.jp/kegg/catalog/org_list.html)
- Reactome Pathway (via ReactomePA)
- DAVID (via RDAVIDWebService)
- Molecular Signatures Database
	- hallmark gene sets
	- positional gene sets
	- curated gene sets
	- motif gene sets
	- computational gene sets
	- GO gene sets
	- oncogenic signatures
	- immunologic signatures
- Other Annotations
	- from other sources (e.g. DisGeNET as an example)
	- user’s annotation
	- customized ontology
	- and many others

### Visulization
- barplot
- cnetplot
- dotplot
- enrichMap
- gseaplot
- plotGOgraph (via topGO package)
- upsetplot

### Citation
G Yu, LG Wang, Y Han, QY He. clusterProfiler: an R package for comparing biological themes among gene clusters. OMICS: A Journal of Integrative Biology 2012, 16(5):284-287. doi:[10.1089/omi.2011.0118](http://dx.doi.org/10.1089/omi.2011.0118)

### Introduction
In recently years, high-throughput experimental techniques such as microarray, RNA-Seq and mass spectrometry can detect cellular molecules at systems-level. These kinds of analyses generate huge quantitaties of data, which need to be given a biological interpretation. A commonly used approach is via clustering in the gene dimension for grouping different genes based on their similarities1.

To search for shared functions among genes, a common way is to incorporate the biological knowledge, such as Gene Ontology (GO) and Kyoto Encyclopedia of Genes and Genomes (KEGG), for identifying predominant biological themes of a collection of genes.

After clustering analysis, researchers not only want to determine whether there is a common theme of a particular gene cluster, but also to compare the biological themes among gene clusters. The manual step to choose interesting clusters followed by enrichment analysis on each selected cluster is slow and tedious. To bridge this gap, we designed clusterProfiler2, for comparing and visualizing functional profiles among gene clusters.

### Installation
```
BiocInstaller::biocLite("AnnotationDbi")
BiocInstaller::biocLite("IRanges")
BiocInstaller::biocLite("tibble")
BiocInstaller::biocLite("DO.db")
BiocInstaller::biocLite("clusterProfiler")
```

### bir: Biological Id Translator
```
bitr(geneID, fromType, toType, OrgDb="org.Hs.eg.db", drop = TRUE)
library(org.Hs.eg.db)
keytypes(org.Hs.eg.db)
## Examples:
x <- c("GPX3",  "GLRX",   "LBP",   "CRYAB", "DEFB1", "HCLS1",   "SOD2",   "HSPA2")
eg = bitr(x, fromType="SYMBOL", toType="ENTREZID", OrgDb="org.Hs.eg.db")
head(eg)
ids <- bitr(x, fromType="SYMBOL", toType=c("UNIPROT", "ENSEMBL"), OrgDb="org.Hs.eg.db")
head(ids)
```
> For GO analysis, user don’t need to convert ID, all ID type provided by `OrgDb` can be used in `groupGO`,  `enrichGO` and `gseGO` by specifying keytype parameter.

### bitr_kegg: converting biological IDs using KEGG API
```
data(gcSample)
hg <- gcSample[[1]]
head(hg)
eg2np <- bitr_kegg(hg, fromType='kegg', toType='ncbi-proteinid', organism='hsa')
head(eg2np)
```
> The ID type (both fromType & toType) should be one of ‘kegg’, ‘ncbi-geneid’, ‘ncbi-proteinid’ or ‘uniprot’. The ‘kegg’ is the primary ID used in KEGG database. The data source of KEGG was from NCBI. A rule of thumb for the ‘kegg’ ID is entrezgene ID for eukaryote species and Locus ID for prokaryotes.

### GO Analysis
**Any gene ID type that supported in OrgDb can be directly used in GO analyses.**
#### Supported organisms
GO analyses (`groupGO()`, `enrichGO()` and `gseGO()`) support organisms that have an `OrgDb` object available.
If user have GO annotation data (in `data.frame` format with first column of gene ID and second column of GO ID), they can use `enricher()` and `gseGO()` functions to perform over-representation test and gene set enrichment analysis.
If genes are annotated by direction annotation, it should also annotated by its ancestor GO nodes (indirect annation). If user only has direct annotation, they can pass their annotation to `buildGOmap` function, which will infer indirection annotation and generate a `data.frame` that suitable for both `enricher()` and  `gseGO()`.

#### GO classification
In clusterProfiler, `groupGO` is designed for gene classification based on GO distribution at a specific level. 
```
data(geneList, package = "DOSE")
gene = names(geneList[abs(geneList) > 2])
ggo <- groupGO(gene     = gene, 
               OrgDb    = org.Hs.eg.db,
               ont      = "CC",
               level    = 3,
               readable = TRUE)
```
> The input parameters of gene is a vector of gene IDs (can be any ID type that supported by corresponding  OrgDb).If readable is setting to TRUE, the input gene IDs will be converted to gene symbols.

#### GO over-representation test
```
ego <- enrichGO(gene          = gene,
                universe      = names(geneList),
                OrgDb         = org.Hs.eg.db,
                keytype		  = 'ENTREZID',
                ont           = "CC",
                pAdjustMethod = "BH",
                pvalueCutoff  = 0.01,
                qvalueCutoff  = 0.05,
        readable      = TRUE)
head(ego)
```
- drop specific GO terms or level
`enrichGO` test the whole GO corpus and enriched result may contains very general terms. With `dropGO` function, user can remove specific GO terms or GO level from results obtained from both `enrichGO` and  `compareCluster`.

- test GO at sepcific level
`enrichGO` doesn’t contain parameter to restrict the test at specific GO level. Instead, we provide a function  `gofilter` to restrict the result at specific GO level. It works with results obtained from both `enrichGO` and  `compareCluster`.

- reduce redundancy of enriched GO terms
According to [issue #28](https://github.com/GuangchuangYu/clusterProfiler/issues/28), I implement a `simplify` method to remove redundant GO terms obtained from  enrichGO. An example can be found in [the blog post](https://guangchuangyu.github.io/2015/10/use-simplify-to-remove-redundancy-of-enriched-go-terms/). It internally call `GOSemSim` to calculate similarities among GO terms and remove those highly similar terms by keeping one representative term. The `simplify` method also works with both outputs from `enrichGO` and `compareCluster`.

#### GO Gene Set Enrichment Analysis
All genes can be used in `GSEA`; `GSEA` aggregates the per gene statistics across genes within a gene set, therefore making it possible to detect situations where all genes in a predefined set change in a small but coordinated way. Since it is likely that many relevant phenotypic differences are manifested by small but consistent changes in a set of genes.
```
ego3 <- gseGO(geneList     = geneList,
              OrgDb        = org.Hs.eg.db,
              ont          = "CC",
              nPerm        = 1000,
              minGSSize    = 100,
              maxGSSize    = 500,
              pvalueCutoff = 0.05,
              verbose      = FALSE)
```
> GSEA use permutation test, user can set nPerm for number of permutations. Only gene Set size in  [minGSSize, maxGSSize] will be tested.

#### GO Semantic Similarity Analysis
GO semantic similarity can be calculated by `GOSemSim`. We can use it to cluster genes/proteins into different clusters based on their functional similarity and can also use it to measure the similarities among GO terms to reduce the redundancy of GO enrichment results.

### KEGG analysis
The annotation package, `KEGG.db`, is not updated since 2012. It’s now pretty old and in clusterProfiler,  `enrichKEGG` (for KEGG pathway) and `enrichMKEGG` (for KEGG module) supports downloading latest online version of KEGG data for enrichment analysis. Using `KEGG.db` is also supported by explicitly setting `use_internal_data` parameter to TRUE, but it’s not recommended.

With this new feature, organism is not restricted to those supported in previous release, it can be any species that have KEGG annotation data available in KEGG database. User should pass abbreviation of academic name to the organism parameter. The full list of KEGG supported organisms can be accessed via http://www.genome.jp/kegg/catalog/org_list.html.

clusterProfiler provides `search_kegg_organism()` function to help searching supported organisms.
```
search_kegg_organism('ece', by='kegg_code')
ecoli <- search_kegg_organism('Escherichia coli', by='scientific_name')
```

#### KEGG over-representation test
```
kk <- enrichKEGG(gene         = gene,
                 organism     = 'hsa',
                 pvalueCutoff = 0.05)
```
> **Input ID type can be kegg, ncbi-geneid, ncbi-proteinid or uniprot**, an example can be found in [the post](https://guangchuangyu.github.io/2016/05/convert-biological-id-with-kegg-api-using-clusterprofiler/).

#### KEGG Gene Set Enrichment Analysis
```
kk2 <- gseKEGG(geneList     = geneList,
               organism     = 'hsa',
               nPerm        = 1000,
               minGSSize    = 120,
               pvalueCutoff = 0.05,
               verbose      = FALSE)
```

#### KEGG Module over-representation test
KEGG Module is a collection of manually defined function units. In some situation, KEGG Modules have a more straightforward interpretation.
```
mkk <- enrichMKEGG(gene = gene,
                   organism = 'hsa')
```

#### KEGG Module Gene Set Enrichment Analysis
```
mkk2 <- gseMKEGG(geneList = geneList,
                 species = 'hsa')
```

### Diease analysis
DOSE supports Disease Ontology (DO) Semantic and Enrichment analysis. The `enrichDO` function is very useful for identifying disease association of interesting genes, and function `gseDO` function is designed for gene set enrichment analysis of DO.

In addition, DOSE also supports enrichment analysis of Network of Cancer Gene (NCG) and Disease Gene Network, please refer to the DOSE vignettes.

### Reactome pathway analysis
ReactomePA uses Reactome as a source of pathway data. The function call of `enrichPathway` and  `gsePathway` in ReactomePA is consistent with `enrichKEGG` and `gseKEGG`.

### DAVID functional analysis
To bridge the gap between DAVID and clusterProfiler, we implemented `enrichDAVID`. This function query enrichment analysis result from DAVID webserver via `RDAVIDWebService` and stored the result as an `enrichResult` instance, so that we can use all the visualization functions in clusterProfiler to visualize DAVID results. `enrichDAVID` is fully compatible with  `compareCluster` function and comparing enrichment results from different gene clusters is now available with DAVID.
```
david <- enrichDAVID(gene = gene,
                     idType = "ENTREZ_GENE_ID",
                     listType = "Gene",
                     annotation = "KEGG_PATHWAY",
                     david.user = "clusterProfiler@hku.hk")
```
DAVID Web Service has the following limitations:

- A job with more than 3000 genes to generate gene or term cluster report will not be handled by DAVID due to resource limit.
- No more than 200 jobs in a day from one user or computer.
- DAVID Team reserves right to suspend any improper uses of the web service without notice.

For more details, please refer to http://david.abcc.ncifcrf.gov/content.jsp?file=WS.html.

As user has limited usage, please register and use your own user account to run enrichDAVID.

### Universal enrichment analysis
clusterProfiler supports both hypergeometric test and gene set enrichment analyses of many ontology/pathway, but it’s still not enough for users may want to analyze their data with unsupported organisms, slim version of GO, novel functional annotation (e.g. GO via BlastGO or KEGG via KAAS), unsupported ontologies/pathways or customized annotations.

clusterProfiler provides `enricher` function for hypergeometric test and `GSEA` function for gene set enrichment analysis that are designed to accept user defined annotation. They accept two additional parameters `TERM2GENE` and `TERM2NAME`. As indicated in the parameter names, `TERM2GENE` is a `data.frame` with first column of term ID and second column of corresponding mapped gene and `TERM2NAME` is a `data.frame` with first column of term ID and second column of corresponding term name. `TERM2NAME` is optional.

An example of using `enricher` and `GSEA` to analyze DisGeNet annotation is presented in the post, [use clusterProfiler as an universal enrichment analysis tool](http://guangchuangyu.github.io/2015/05/use-clusterprofiler-as-an-universal-enrichment-analysis-tool/).

### Using MSigDB gene set collections
The MSigDB is a collection of annotated gene sets, it include 8 major collections:
- H: hallmark gene sets
- C1: positional gene sets
- C2: curated gene sets
- C3: motif gene sets
- C4: computational gene sets
- C5: GO gene sets
- C6: oncogenic signatures
- C7: immunologic signatures

Users can use `enricher` and `GSEA` function to analyze gene set collections downloaded from Molecular Signatures Database (MSigDb). clusterProfiler provides a function, `read.gmt`, to parse the gmt file into a `TERM2GENE` data.frame that is ready for both `enricher` and `GSEA` functions.

```
gmtfile <- system.file("extdata", "c5.cc.v5.0.entrez.gmt", package="clusterProfiler")
c5 <- read.gmt(gmtfile)

egmt <- enricher(gene, TERM2GENE=c5)

egmt2 <- GSEA(geneList, TERM2GENE=c5, verbose=FALSE)
```

### Functional analysis of NGS data
Functional analysis using NGS data (eg, RNA-Seq and ChIP-Seq) can be performed by linking coding and non-coding regions to coding genes via ChIPseeker package, which can annotates genomic regions to their nearest genes, host genes, and flanking genes respectivly. In addtion, it provides a function, `seq2gene`, that simultaneously considering host genes, promoter region and flanking gene from intergenic region that may under control via cis-regulation. This function maps genomic regions to genes in a many-to-many manner and facilitate functional analysis. For more details, please refer to ChIPseeker.

### Visualization
The function calls of `groupGO`, `enrichGO`, `enrichKEGG`, `enrichDO`, `enrichPathway` and `enricher` are consistent and all the output can be visualized by bar plot, enrichment map and category-gene-network plot. It is very common to visualize the enrichment result in bar or pie chart. We believe the pie chart is misleading and only provide bar chart.

#### barplot
```
barplot(ggo, drop=TRUE, showCategory=12)
barplot(ego, showCategory=8)
```

#### dotplot
`dotplot(ego)`

#### enrichMap
Enrichment map can be viusalized by `enrichMap`, which also support results obtained from hypergeometric test and gene set enrichment analysis.
`enrichMap(ego)`

#### cnetplot
In order to consider the potentially biological complexities in which a gene may belong to multiple annotation categories and provide information of numeric changes if available, we developed `cnetplot` function to extract the complex association.
```
## categorySize can be scaled by 'pvalue' or 'geneNum'
cnetplot(ego, categorySize="pvalue", foldChange=geneList)
```

#### plotGOgraph
`plotGOgraph`, which is based on topGO, can accept output of enrichGO and visualized the enriched GO induced graph.
`plotGOgraph(ego)`

#### gseaplot
Running score of gene set enrichment analysis and its association of phenotype can be visualized by `gseaplot`.
`gseaplot(kk2, geneSetID = "hsa04145")`

#### browseKEGG
To view the KEGG pathway, user can use `browseKEGG` function, which will open web browser and highlight enriched genes.
`browseKEGG(kk, 'hsa04110')`

#### pathview from pathview package
clusterProfiler users can also use `pathview` from the pathview to visualize KEGG pathway.

The following example illustrate how to visualize “hsa04110” pathway, which was enriched in our previous analysis.
```
library("pathview")
hsa04110 <- pathview(gene.data  = geneList,
                     pathway.id = "hsa04110",
                     species    = "hsa",
                     limit      = list(gene=max(abs(geneList)), cpd=1))
```

### Biological theme comparison
clusterProfiler was developed for biological theme comparison, and it provides a function, `compareCluster`, to automatically calculate enriched functional categories of each gene clusters.
```
data(gcSample)
lapply(gcSample, head)
ck <- compareCluster(geneCluster = gcSample, fun = "enrichKEGG")
head(as.data.frame(ck))
# The input for geneCluster parameter should be a named list of gene IDs. To speed up the compilation of this document, we set use_internal_data = TRUE.
```

#### Formula interface of compareCluster
`compareCluster` also supports passing a formula (the code to support formula has been contributed by Giovanni Dall’Olio) of type `Entrez∼group` or `Entrez∼group+othergroup`.
```
mydf <- data.frame(Entrez=names(geneList), FC=geneList)
mydf <- mydf[abs(mydf$FC) > 1,]
mydf$group <- "upregulated"
mydf$group[mydf$FC < 0] <- "downregulated"
mydf$othergroup <- "A"
mydf$othergroup[abs(mydf$FC) > 2] <- "B"

formula_res <- compareCluster(Entrez~group+othergroup, data=mydf, fun="enrichKEGG")

head(as.data.frame(formula_res))
```

#### Visualization of profile comparison
`dotplot(ck)` or `dotplot(formula_res)`
`dotplot(formula_res, x=~group) + ggplot2::facet_grid(~othergroup)`
![](https://i.imgur.com/iw5Xm88.png)

By default, only top 5 (most significant) categories of each cluster was plotted. User can changes the parameter `showCategory` to specify how many categories of each cluster to be plotted, and if `showCategory` was set to **NULL**, the whole result will be plotted.

The plot function accepts a parameter `by` for setting the scale of dot sizes. The default parameter `by` is setting to “geneRatio”, which corresponding to the “GeneRatio” column of the output. If it was setting to count, the comparison will be based on gene counts, while if setting to rowPercentage, the dot sizes will be normalized by count/(sum of each row)

Color gradient ranging from red to blue correspond to in order of increasing p-values. That is, red indicate low p-values (high enrichment), and blue indicate high p-values (low enrichment). P-values and adjusted p-values were filtered out by the threshold giving by parameter `pvalueCutoff`, and FDR can be estimated by `qvalue`.