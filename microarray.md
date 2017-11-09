## Microarray Analysis
1. NCBI-GEO
* GDS# 数据集编号
* GPL# 检测平台类型编号 ["bioconductor中对应平台数据库"](http://www.bio-info-trainee.com/1399.html )
* GSE# 研究项目编号
* GSM# 样本编号

2. 数据urls
* 表达矩阵：（加粗处需修改）
ftp://ftp.ncbi.nlm.nih.gov/geo/series/**GSE1nnn**/**GSE1009**/matrix/**GSE1009**_series_matrix.txt.gz
* 芯片原始数据：
ftp://ftp.ncbi.nlm.nih.gov/geo/series/**GSE1nnn**/**GSE1009**/suppl/**GSE1009**_RAW.tar
* 样本分组信息：
https://www.ncbi.nlm.nih.gov/geo/browse/?view=samples&search=**GSE51068**&zsort=date&mode=csv&page=undefined&display=5000
* GPL探针注释文件：
1）ftp://ftp.ncbi.nlm.nih.gov/geo/platforms/**GPL13nnn**/**GPL13667**/soft/**GPL13667**_family.soft.gz
2）https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?mode=raw&is_datatable=true&acc=GPL13667&id=15572&db=GeoDb_blob92

3. ExpressionSet对象：表达矩阵加上样本信息的一个封装
组成：
* assayData：一个matrix类型或environment类型的数据。用于保存表达数据值。(不能缺失)
* 头文件：用于描述实验平台相关数据，其中包括phenoData，featureData，protocolData以及annotation等。
> phenoData是一个存放样品信息的data.frame或者AnnotatedDataFrame类型的数据。如果有行号的话，其行号必须与assayData的列号一致（也就是样品名）。如果没有行号，则其行数必须与assayData的列数一致。
featureData是一个存放features的data.frame或者AnnotatedDataFrame类型的数据。它的行数必须与assayData的行数一致。如果有行号的话，那么它的行号必须和assayData的行号一致。
annotation是用于存放芯片类型的字符串，比如hgu95av2之类。
protocolData用于存放设备相当的数据。它是AnnotatedDataFrame类型。它的维度必须与assayData的维度一致。
* experimentData: 一个MIAME类型的数据，它用于保存和实验设计相关的资料，比如实验室名，发表的文章，等等。那么什么是MIAME类呢？MIAME是Minimum Information About a Microarray Experiment的首字母缩写，它包括以下一些属性（slots）：
> name: 字符串，实验名称
lab: 字符串，实验室名称
contact: 字符串，联系方式
title: 字符串，一句话描述实验的内容
abstract: 字符串，实验摘要
url: 字符串，实验相关的网址
samples: list类，样品的信息
hybridizations: list类，杂交的信息
normControls: list类，对照信息，比如一些持家基因（house keeping genes）
preprocessing: list类，原始数据的预处理过程
pubMedIds: 字符串，pubMed索引号
others: list类，其它相关的信息
```
	library(CLL)
	data(sCLLex)
	exprMatrix=exprs(sCLLex)
	meta=pData(sCLLex)
	# 自己构造 ExpressionSet 对象
 	metadata <- data.frame(labelDescription=c('SampleID', 'Disease'),row.names=c('SampleID', 'Disease'))
	phenoData <- new("AnnotatedDataFrame",data=meta,varMetadata=metadata)
	myExpressionSet <- ExpressionSet(assayData=exprMatrix,
                                     phenoData=phenoData,
                                     annotation="hgu95av2")
```

4. R
* 对象函数查询
```
	class(object)
	[1] 'type'
	methods(class='type')
```
* [R语言里面的一个数据集ALL(Acute Lymphoblastic Leukemia)简介](http://www.bio-info-trainee.com/693.html)
> 这个数据集是对ALL这个病的研究数据，共涉及到了128个ALL病人，其中95个是B细胞的ALL，剩余33个是T细胞的ALL。
> ALL$BT：记录病人的分组信息，bcell = grep("^B", as.character(ALL$BT))通过这句话可以挑选出B细胞病人。
> ALL$mol.biol：记录病人的几种突变情况（molecular biology testing）
> types = c("NEG", "BCR/ABL")；moltyp = which(as.character(ALL$mol.biol) %in% types)
> 包括sex，age，cod，diagnosis，等等，这个'data.frame':共有128 obs. of  21 variables
> str(exprs(ALL)) 查看表达数据集

* [affy包读取affymetrix基因表达芯片数据](http://www.bio-info-trainee.com/1580.html)：
```
	library(affy)
	setwd("...")
	data.raw <- ReadAffy(celfile.path = getwd())
	# normalization:
	eset.mas5 <- mas5(data.raw)
	eset.rma <- rma(data.raw)
	# filtering:
	calls  <- mas5calls(data.raw)
	calls <- exprs(calls)
	absent <- rowSums(calls == 'A') # how may samples are each gene 'absent' in all samples
	absent <- which (absent == ncol(calls)) # which genes are 'absent' in all samples
	rmaFiltered <- eset.rma[-absent,] # filters out the genes 'absent' in all samples
	write.exprs(eset.rma,file="Exprs.txt")
```

* [oligo包读取affymetrix基因表达芯片数据](http://www.bio-info-trainee.com/1586.html):
> [HuGene-1_1-st] Affymetrix Human Gene 1.1 ST Array这个平台虽然也是affymetrix公司的，但是affy包就无法处理 了，这时候就需要oligo包了！
```
	geneCELs=list.celfiles('/path/GSE48452/cel_files/',listGzipped=T,full.name=T)
	#用全路径，一般cel文件也是压缩包形式，没必要解压
	affyGeneFS <- read.celfiles(geneCELs)  ##读取ｃｅｌ文件
	geneCore <- rma(affyGeneFS, target = "core")　 ##这一步是normalization，会比较耗时
	genePS <- rma(affyGeneFS, target = "probeset")
	#两种normlization的方法，##一般我们会选择transcript相关的
	## 这个芯片平台还需要自己把探针ID赋值给表达矩阵
	featureData(genePS) <- getNetAffx(genePS, "probeset")
	featureData(geneCore) <- getNetAffx(geneCore, "transcript")
```
![Gene 2.0 ST Arrays](http://i.imgur.com/smSDlKR.png)

* [lumi包处理illumina的bead系列表达芯片](http://www.bio-info-trainee.com/1944.html): 
> Paper:[lumi: a pipeline for processing Illumina microarray](https://oup.silverchair-cdn.com/oup/backfile/Content_public/Journal/bioinformatics/24/13/10.1093/bioinformatics/btn224/2/btn224.pdf?Expires=1498901890&Signature=SBbblMiTJ0QYMnKOg~keIzMs6zAbDnmkId7CCG9B6c61zWSQbe9-NofG0QXf~QfFuDMMjZDy69vqi4OPzfoCqHkOXUMxTS~gI~JNDv-aAX3QGpcPWAkDmfjK19NH6HpEF47~oTHlnicmKY2mkCGPBaw6icjAmSu3YPq5WFp3Ge985Mqu2BELzSqEOsFN5AcGGtRqtsXFAB9wiCQv4Cx2Lpt~UTZDxwfxJiIiXkE5AF1NbqIxgoBOJlaJeSsyha~WbbUtBE7KISqbzrgnu387ywv~xK7NDi23ODYuKP7VpM5CLX4nrJ5RQOKJu66u0ZXTlcGLoKiZrwLPaKWTEMSEjA__&Key-Pair-Id=APKAIUCZBIA4LVPAVW3Q)
> 对象：LumiBatch
> `lumiExpresso()` 可以直接处理LumiBatch对象，这个函数结合了,N,T,B,Q(normalization,transformation,backgroud correction,qulity control)四个步骤。
```
	x.lumi <- lumiR.batch() # 读取
	pData(phenoData(x.lumi))
	lumi.N.Q <- lumiExpresso(x.lumi)
	dataMatrix <- exprs(lumi.N.Q)
```

* [芯片探针注释基因ID或者gene symbol，并对每个基因挑选最大表达量探针](http://www.bio-info-trainee.com/1502.html)
1> 手动注释： 自己准备`anno.txt`(一列探针，一列gene ID or gene symbol)，R 用`match`函数：`anno$gene.symbol[match(rownames(exprSet),anno$probeID)]`即可实现。
注释结果一般含有 "---" or "///"
2> bioconductor注释：下载GPL相应的`xxx.db`: 如`BiocInstaller::biocLite("hgu95av2.db")`;
```
	library(hgu95av2.db)
	probeset = rownames(exprSet)
	symbol = as.character(as.list(hgu95av2SYMBOL[probeset]))
	# annotate 包提供 getSYMBOL(probeset,"hgu95av2")
	# 还可以用 lookUp(probeset, "hgu95av2", "SYMBOL")
	exprSet = cbind.data.frame(symbol, exprSet)
```
3> 对每个基因挑选最大表达量探针
```
	rmDupID <- function(dup_exprSet){
  
  print(paste('Total input records:',nrow(dup_exprSet),sep = " "))
  
  exprSet=dup_exprSet[,-1] ## first column is the ID needed to remove duplicate.probably is gene symbol or entrez ID.
  rowMeans=apply(exprSet,1,function(x) mean(as.numeric(x),na.rm=T))
  dup_exprSet=dup_exprSet[order(rowMeans,decreasing=T),]
  exprSet=dup_exprSet[!duplicated(dup_exprSet[,1]),]
  #exprSet=apply(exprSet,2,as.numeric)
  exprSet=exprSet[!is.na(exprSet[,1]),]
  exprSet=exprSet[exprSet[,1] !='NA',]
  #exprSet=exprSet[exprSet[,1] !='---',]
  #exprSet=exprSet[!str_detect(exprSet[,1],'///'),]
  rownames(exprSet)=exprSet[,1]
  exprSet=exprSet[,-1]
  #str(exprSet)
  rn=rownames(exprSet)
  exprSet=apply(exprSet,2,as.numeric)
  rownames(exprSet)=rn
  #exprSet[1:4,1:4]
  print(paste("Total output records:",nrow(exprSet),sep = " "))
  return(exprSet)
}
	exprSet = rmDupID(exprSet)
```

* Others
>\>90% of all studies using Affy chips will use RMA for normalization.

5. [GSEA, Gene Set Enrichment Analysis](http://software.broadinstitute.org/gsea/index.jsp)
> 下载地址：http://software.broadinstitute.org/gsea/downloads.jsp
> 好像只有32-bit的版本可以用，包括gsea2-2.2.4.jar都不能用-d64
> GSEA supported data files are simply tab delimited ASCII text files, which have special file extensions that identify them. For example, expression data usually has the extension *.gct, phenotypes *.cls, gene sets *.gmt, and chip annotations *.chip. Click the More on file formats help button to view detailed descriptions of all the data file formats.
并且提供了测试数据：http://software.broadinstitute.org/gsea/datasets.jsp
> 使用手册：[GSEA Useer Guide](http://software.broadinstitute.org/gsea/doc/GSEAUserGuideFrame.html?_GSEAPreranked_Page)
* [Exmaple1：Desktop Application](http://www.bio-info-trainee.com/1282.html)
> 注：需要将表格中的“NA”删除，不然读入会报错
> * GSEA Statistics：
GSEA computes four key statistics for the gene set enrichment analysis report:
Enrichment Score (ES)
Normalized Enrichment Score (NES)
False Discovery Rate (FDR)
Nominal P Value
> * GSEA Report
This section discusses the content of the report generated by the gene set enrichment analysis:
Enrichment in Phenotype
Dataset Details
Gene Set Details
Gene Markers
Global Statistics and Plots
Other
Detailed Enrichment Results
Gene Set Details Report

* [Example 2: command line or offline](http://www.bio-info-trainee.com/1334.html)
> 芯片注释文件：ftp://ftp.broadinstitute.org/pub/gsea/annotations
> gene sets 文件：http://software.broadinstitute.org/gsea/downloads.jsp#msigdb
```
java -cp gsea2-2.2.4.jar -Xmx1024m xtools.gsea.Gsea --help #查看帮助文档
java -cp gsea2-2.2.4.jar -Xmx1024m xtools.gsea.Gsea -gmx [gene sets annotation].gmt -res [exprs file].gct -cls [groups info].cls -chip [chip anno file].chip -out [out dir] -rpt_label [label name] # 如果gct文件的第一列是芯片探针，才需要使用-chip
```

* [原理](http://note.youdao.com/share/?token=DBDB0277A315444BBBAB2024190208AE&gid=23785548#/)：
> 首先对每个样本里面的基因的表达值在样本内部进行排序，本质是是根据该基因在两个group之间的差异来排序！但是差异如何量化，就有多种方法了，可以是Signal2Noise 值，或者是Ttest值，或者是fold change，logFC。默认的，GSEA会根据signal-to-noise metric 来对基因进行排序。`-metric <Metric>`修改，Hints: Signal2Noise,tTest,Cosine,Euclidean,Manhatten,Pearson,Ratio_of_Classes,Diff_of_Classes,log2_Ratio_of_Classes。如果是自己已经排序好了的基因，也可以直接拿来做GSEA分析了，见： [GSEAPreranked Page](http://software.broadinstitute.org/gsea/doc/GSEAUserGuideFrame.html?_GSEAPreranked_Page) in the GSEA User Guide.
> 如果是affymetrix的表达矩阵，不需要提前进行`Present/Marginal/Absent Calls`来过滤掉一些表达探针。
> 如果是gct and pcl 的表达矩阵，缺失值空着就好了。但是如果缺失值太多了，这样在计算signal-to-noise的时候，不同group的样本数就不一致了，mean和sd都会变好，最好是避免这样的情况，可以考虑进行插值，或者过滤掉这样的探针。
> 不需要提前过滤掉低表达量的探针或者低variance的探针。它们都会在我们算好的 ranked gene list 的中间部分，增强我们的统计效应。完全不用担心数据量计算时间的问题。
> 果要想计算Signal2Noise ，每个group必须要有3个及以上的samples。
除了两个group之间的比较可以做gsea之外，还可以针对连续性的phenotypes和time-course数据。

6. [差异分析](http://note.youdao.com/share/?token=8D29A5320C934AB29BFAF4052B890DBD&gid=23785548#/)
* [用excel表格做差异分析](http://www.bio-info-trainee.com/1205.html)
* [用samr包对芯片数据做差异分析](http://www.bio-info-trainee.com/1608.html) or http://blog.qiubio.com:8080/archives/2083/3
* [用limma包对芯片数据做差异分析](http://www.bio-info-trainee.com/1194.html)
![](http://i.imgur.com/53si2r5.png)
> 基因表达矩阵：
```
	exprSet=read.table("GSE63067_series_matrix.txt.gz",comment.char = "!",stringsAsFactors=F,header=T)
	rownames(exprSet)=exprSet[,1]
	exprSet=exprSet[,-1]
```
> 分组矩阵：
![](http://i.imgur.com/PjvYX2B.png)
> 差异比较矩阵：
![](http://i.imgur.com/lnOk6P2.png)

7. [富集分析](http://note.youdao.com/share/?token=5C72F0011F2E4CE3B11FA4981602E94D&gid=23785548#/)
> 超几何分布很简单，球分成黑白两色，总数为N，M个白球，那么你随机抽有n个球，应该抽多少白球的问题！
公式就是 exp\_count=n\*M\/N（理论值）然后你实际上抽了k白球，就可以计算一个概率值（OddsRatio=k\/(n*(M\/N))）！
> 换算成通路的富集概念就是，总共有多少基因，你的通路有多少基因，你的通路被抽中了多少基因（在差异基因里面属于你的通路的基因）
> 原始注释文件内为map geneID to GO ID
* [GOstats 方法](https://www.plob.org/article/7126.html)
```
	library("org.Hs.eg.db")
	library("GSEABase")
	library("GOstats")
	genes <- c("AREG", "FKBP5", "CXCL13", "KLF9", "ZC3H12A", "P4HA1", "TLE1", "CREB3L2", "TXNIP", "PBX1", "GJA1", "ITGB8", "CCL3", "CCND2", "KCNJ15", "CFLAR", "CXCL10", "CYSLTR1", "IGFBP7", "RHOB", "MAP3K5", "CAV2", "CAPN2", "AKAP13", "RND3", "IL6ST", "RGS1", "IRF4", "G3BP1", "SEL1L", "VEGFA", "SMAD1", "CCND1", "CLEC3B", "NEB", "AMD1", "PDCD4", "SCD", "TM2D3", "BACH2", "LDLR", "BMPR1B", "RFXAP", "ASPH", "PTK2B", "SLC1A5", "ENO2", "TRPM8", "SATB1", "MIER1", "SRSF1", "ATF3", "CCL5", "MCM6", "GCH1", "CAV1", "SLC20A1")
	goAnn <- get("org.Hs.egGO")
	universe <- Lkeys(goAnn)
	entrezIDs <- mget(genes, org.Hs.egSYMBOL2EG, ifnotfound=NA)
	entrezIDs <- as.character(entrezIDs)
	params <- new("GOHyperGParams",
		geneIds=entrezIDs,
		universeGeneIds=universe,
		annotation="org.Hs.eg.db",
		ontology="BP", # ontology = c("BP","CC","MF")
		pvalueCutoff=0.01,
		conditional=FALSE,
		testDirection="over")
	over <- hyperGTest(params)
	library(Category)
	glist <- geneIdByCategory(over)
	glist <- sapply(glist, function(.ids) {
		.sym <- mget(.ids, envir=org.Hs.egSYMBOL, ifnotfound=NA)
		.sym[is.na(.sym)] <- .ids[is.na(.sym)]
		paste(.sym, collapse=";")
		})
	bp <- summary(over)
	bp$Symbols <- glist[as.character(bp$GOBPID)]
	## KEGG enrichment analysis
	library(KEGG.db)
	keggAnn <- get("org.Hs.egPATH")
	universe <- Lkeys(keggAnn)
	params <- new("KEGGHyperGParams", 
		geneIds=entrezIDs, 
		universeGeneIds=universe, 
		annotation="org.Hs.eg.db", 
		categoryName="KEGG", 
		pvalueCutoff=0.01,
		testDirection="over")
	over <- hyperGTest(params)
	kegg <- summary(over)
	library(Category)
	glist <- geneIdByCategory(over)
	glist <- sapply(glist, function(.ids) {
		.sym <- mget(.ids, envir=org.Hs.egSYMBOL, ifnotfound=NA)
		.sym[is.na(.sym)] <- .ids[is.na(.sym)]
		paste(.sym, collapse=";")
		})
	kegg$Symbols <- glist[as.character(kegg$KEGGID)]
	library("pathview")
	gIds <- mget(genes, org.Hs.egSYMBOL2EG, ifnotfound=NA)
	gEns <- unlist(gIds)
	gene.data <- rep(1, length(gEns)
	names(gene.data) <- gEns
	for(i in 1:3){pv.out <- pathview(gene.data, pathway.id=as.character(kegg$KEGGID)[i], species="hsa", out.suffix="pathview", kegg.native=T)}
```
函数化：
```
enrich <- function(entrezIDs, orgDbName="org.Hs.eg.db", pvalueCutoff=.01){
	require(orgDbName, character.only=TRUE)
	require("GSEABase")
	require("GOstats")
	require("Category")
	require("KEGG.db")
	goAnn <- get(gsub(".db", "GO", orgDbName))
	universe <- Lkeys(goAnn)
	onto <- c("BP", "MF", "CC")
	res <- lapply(onto, function(.onto){
		param <- new('GOHyperGParams',
					 geneIds= entrezIDs,
					 universeGeneIds=universe,
					 annotation=orgDbName,
					 ontology=.onto,
					 pvalueCutoff=pvalueCutoff,
					 conditional=FALSE,
					 testDirection="over")
		over <- hyperGTest(param)
		glist <- geneIdsByCategory(over)
		glist <- sapply(glist, function(.ids) {
			.sym <- mget(.ids, envir=get(gsub(".db", "SYMBOL", orgDbName)), ifnotfound=NA)
			.sym[is.na(.sym)] <- .ids[is.na(.sym)]
			paste(.sym, collapse=";")
		})
		summary <- summary(over)
		if(nrow(summary)>1) summary$Symbols <- glist[as.character(summary[, 1])]
		summary
	})
	names(res) <- onto
	keggAnn <- get(gsub(".db", "PATH", orgDbName))
	universe <- Lkeys(keggAnn)
	param <- new("KEGGHyperGParams",
				 geneIds=entrezIDs,
				 universeGeneIds=universe,
				 annotation=orgDbName,
				 categoryName="KEGG",
				 pvalueCutoff=pvalueCutoff,
				 testDirection="over")
	over <- hyperGTest(param)
	kegg <- summary(over)
	glist <- geneIdsByCategory(over)
	glist <- sapply(glist, function(.ids) {
		.sym <- mget(.ids, envir=get(gsub(".db", "SYMBOL", orgDbName)), ifnotfound=NA)
		.sym[is.na(.sym)] <- .ids[is.na(.sym)]
		paste(.sym, collapse=";")
	})
	kegg$Symbols <- glist[as.character(kegg$KEGGID)]
	res[["kegg"]] <- kegg
	res
}
```
8. PPI Analysis
* [STINGdb](http://www.bioconductor.org/packages/release/bioc/html/STRINGdb.html)
> 它只需要一个3列的data.frame，分别是logFC,p.value,gene ID,就是标准的差异分析的结果。
> 然后用string_db$map函数给它加上一列是 string 数据库的蛋白ID，然后用string_db$add_diff_exp_color函数给它加上一列是color。用string_db$plot_network函数画网络图，只需要 string 数据库的蛋白ID，如果需要给蛋白标记不同的颜色，需要用string_db$post_payload来把color对应到每个蛋白，然后再画网络图。也可以直接用get_interactions函数得到所有的PPI数据，然后写入到本地，再导入到cytoscape进行画图。
* Basic Use：
```
	library(STRINGdb)
	string_db <- STRINGdb$new(version="10", species=9606,score_threshold=0,input_directory="")
	## Besides, if you specify a local directory to the parameter input-directory, all the database files will be downloaded into this directory and the package can then be used off-line.
	## retrieve/search the species table:
	get_STRING_species(version="10", species_name=NULL)

	STRINGdb$methods() # To list all the methods available.
	STRINGdb$help("get_graph") # To visualize their documentation.
	
	data(diff_exp_example1) ## from GSE9008, containing three columns: pvalue, logFC, gene
	
	## map the gene names to the STRING database identifiers using the "map" method.mapping function supports several common identifiers (e.g. HUGO, Entrez GeneID, ENSEMBL proteins, RefSeq transcripts ... etc.).
	example1_mapped <- string_db$map(diff_exp_example1, "gene", removeUnmappedRows=TRUE ) ## returns the input dataframe with the "STRING_id" additional column.
	hits <- example1_mapped$STRING_id[1:200]
	string_db$plot_network(hits)
	## filter by p-value and add a color column
	example1_mapped_pval05 <- string_db$add_diff_exp_color(subset(example1_mapped, pvalue < 0.05), logFcColStr="logFC")
	## post payload information to the STRING server
	payload_id <- string_db$post_payload(example1_mapped_pval05$STRING_id, colors=example1_mapped_pval05$color)
	## display a STRING network png with the "halo"
	string_db$plot_network(hits, payload_id=payload_id)
	
	## perform enrichment analysis
	## plot the enrichment for the best 1000 genes
	string_db$plot_ppi_enrichment( example1_mapped$STRING_id[1:1000], quiet=TRUE )
	## you should see more enrichment at the beginning of the list than at the end. you can also use the enrichment graph to help you to define a threshold on the number of proteins to consider.

	##  Gene Ontology, KEGG pathway and Interpro domains enrichment
	enrichmentGO <- string_db$get_enrichment( hits, category = "Process", methodMT = "fdr", iea = TRUE)
	enrichmentKEGG <- string_db$get_enrichment( hits, category = "KEGG", methodMT = "fdr", iea = TRUE )
	## If you have performed your experiment on a predefined set of proteins, it is important to run the enrichment statistics using that set as a background (otherwise you would get a wrong p-value !)
	backgroundV <- example1_mapped$STRING_id[1:2000]
	string_db$set_background(backgroundV)
	## You can also set the background when you instantiate the STRINGdb object:
	string_db <- STRINGdb$new( score_threshold=0, backgroundV = backgroundV )
	## compare the enrichment of two or more lists of genes
	eh <- string_db$enrichment_heatmap( list( hits[1:100], hits[101:200]), list("list1","list2"), title="My Lists" )

	## Clustering
	# get clusters
	clustersList <- string_db$get_clusters(example1_mapped$STRING_id[1:600])
	options(SweaveHooks=list(fig=function() par(mar=c(2.1, 0.1, 4.1, 2.1))))
	# plot first 4 clusters
	par(mfrow=c(2,2))
	for(i in seq(1:4)){string_db$plot_network(clustersList[[i]])}
	
	## additional protein information
	string_proteins <- string_db$get_proteins()
	# get the STRING identifier:
	tp53 = string_db$mp( "tp53" )
	atm = string_db$mp( "atm" )
	# see the proteins that interact with one or more of your proteins:
	string_db$get_neighbors( c(tp53, atm) )
	# retrieve the interactions that connect certain input proteins between each other
	string_db$get_interactions( c(tp53, atm) )
	# retrieve the pubmed identifiers of the articles that contain the name of both the proteins (if any)
	string_db$get_pubmed_interaction( tp53, atm )
	# get the reciprocal best hits of the following protein in all the STRING species
	string_db$get_homologs_besthits(tp53, symbets = TRUE)
	# get the homologs of the following two proteins in the mouse (i.e. species_id=10090)
	string_db$get_homologs(c(tp53, atm), target_species_id=10090, bitscore_threshold=60 )

	## benchmarkong protein-protein interaction
	# you need to provide as input a sorted interaction data frame (with the columns "proteinA", "proteinB", "score").
	data(interactions_example)
	interactions_benchmark = string_db$benchmark_ppi(interactions_example, pathwayType = "KEGG", max_homology_bitscore = 60, precision_window = 400, exclude_pathways = "blacklist")
	# plot the precision vs the number of sorted interactions
	plot(interactions_benchmark$precision, ylim=c(0,1), type="l", xlim=c(0,700), xlab="interactions", ylab="precision")
	interactions_pathway_view = string_db$benchmark_ppi_pathway_view(interactions_benchmark, precision_threshold=0.2, pathwayType = "KEGG")
```

9. [Cytoscape](http://www.cytoscape.org/index.html)
* `igraph` Package
> Vertices and edge are indexed from one in R igraph, since version 0.6. Indices are continual, from 1 to |V|.

V = {A, B, C, D, E}
E = ((A, B),(A, C),(B, C),(C, E)).
A = 1, B = 2, C = 3, D = 4, E = 5

`g <- graph( c(1,2, 1,3, 2,3, 3,5), n=5 )`
1> igraph objects
2> print(), summary(), is.igraph(), is.directed(), vcount(), ecount()
3> Visualization: `g <- graph.tree(40, 4); plot(g, layout=layout.circle)`
```
# Force directed layouts
	plot(g, layout=layout.fruchterman.reingold)
	plot(g, layout=layout.graphopt)
	plot(g, layout=layout.kamada.kawai)
```
4> Interactive: `tkplot(g, layout=layout.kamada.kawai); l <- layout.kamada.kawai(g)`
5> 3D: `require(rgl); rglplot(g, layout=l)`
6> Visual properties: `plot(g, layout=l, vertex.color="cyan")`
7> Naming vertices: `V(g)$name <- sample(letters, vcount(g))`
8> Graph Types: `g <- graph.ring(10)`[U]; `g <- graph.formula( Alice-Bob-Cecil-Alice, Daniel-Cecil-Eugene, Cecil-Gordon )`[U]; `g2 <- graph.formula( Alice-Bob:Cecil:Daniel, Cecil:Daniel-Eugene:Gordon )`[U]; `g3 <- graph.formula( Alice +-+ Bob --+ Cecil+-- Daniel, Eugene --+ Gordon:Helen )`[D]; number of '-' can be arbitrary.
9> Vertex/Edge sets, attributes: V(g); E(g)
```
	## Smart indexing
	V(g)[color=="white"]
	g <- erdos.renyi.game(100, 1/100)
	V(g)$color <- sample( c("red", "black"), vcount(g), rep=TRUE)
	E(g)$color <- "grey"
	red <- V(g)[ color == "red" ]
	bl <- V(g)[ color == "black" ]
	E(g)[ red %--% red ]$color <- "red"
	E(g)[ bl %--% bl ]$color <- "black"
	plot(g, vertex.size=5, layout=layout.fruchterman.reingold)
```
10> Centrality in networks: degree; closeness; betweenness; eigenvector centrality; page rank
11> Community structure in networks: How to define what is modular?
12> Cohensive blocks: A collectivity is structurally cohesive to the extent that the social relations of its members hold it together.A group is structurally cohesive to the extent that multiple independent relational paths among all pairs of members hold it together.

* [Network visualization with R](http://kateto.net/network-visualization)

* cytoscape
1> Input formats: http://wiki.cytoscape.org/Cytoscape_User_Manual/Network_Formats
2> Sample Data: C:\Program Files\Cytoscape_v3.5.1\sampleData
3> APP installation: 
   Target directory: C:\Users\colin\CytoscapeConfiguration\3\apps\installed
   Web: http://apps.cytoscape.org/

10. Subnetwork Model
* [BioNet Package to find maximal-scoring subgraph](http://www.bio-info-trainee.com/2071.html)
> Heuristically, 整合多种数据结果给一个网络打分。
> `subNetwork()`从一个大网络提提取a interested list的子网络，graph objects可以是graphNet或者是igraph。
> `aggrPvals()` aggregates several p-values into one p-value of p-values based on the order statistics of p-values. An overall p-value is given by the ith order statistic.
> `fitBumModel` The function fits a beta-uniform mixture model to a given p-value distribution. The BUM method was introduced by Stan Pounds and Steve Morris to model the p-value distribution as a signal-noise decompostion. The signal component is assumed to be B(a,1)-distributed, whereas the noise component is uniform-distributed under the null hypothesis.
> `datadir <- file.path(path.package("..."), "extdata")` # 可以获得extdata的文件夹路径