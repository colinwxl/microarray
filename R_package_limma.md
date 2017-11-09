## limma package( R, Bioconductor)
### focus on dealing affy Microarray
> **Limma** is a package for the analysis of gene expression data arising from microarray or RNA-Seq technologies. A core capability is the use of linear models to assess differential expression in the context of multifactor designed experiments.
> **Limma** provides a strong suite of functions for reading, exploring and pre-processing data from two-color microarrays. Alternation: `marray` package.
> **Limma** can read output data from a variety of image analysis software platforms, including GenePix, ImaGene etc. Either one-channel or two-channel formats can be processed.
> Functions for reading and pre-processing expression data from Illumina BeadChips were introduced in limma 3.0.0.
> From version 3.9.19, limma includes functions to analyse RNA-Seq experiments.
> `read.maimages()` function with arguement `wt.fun`: Imgaes-derived Spot Quality Weights.
> `help(package=limma)` to see LIMMA User's Guide (pdf).
> http://www.bioconductor.org/help/support/posting-guide/

1. Installation
```
	source("http://www.bioconductor.org/biocLite.R")
	biocLite("limma") # BiocInstaller::biocLite("limma")
	biocLite("statmod")
```
2. Quick Start
```
	library(gcrma)
	library(limma)
	targets <- readTargets("targets.txt") 
	# The targets frame normally contains a FileName column, giving the name of the image-analysis output file, 
	# a Cy3 column giving the RNA type labelled with Cy3 dye for that slide and 
	# a Cy5 column giving the RNA type labelled with Cy5 dye for that slide.
	ab <- ReadAffy(filenames = targets$FileName)
	eset <- gcrma(ab)
	design <- model.matrix()
	fit <- lmFit(eset, design)
	contrast.matrix <- makeContrasts(..., levels = design)
	fit2 <- contrasts.fit(fit, contrast.matrix)
	fit2 <- eBayes(fit2)
	topTable(fit, coef=) # show statistics for top 10 genes
	results <- decideTests(fit2)
	vennDiagram(results)
	# one-way ANOVA for each gene except that the residual mean squares have been moderated between genes
	topTable(fit2, number=30)
```
3. Data Objects

|Object   |  Description                                                                   |
|:-------:|:------------------------------------------------------------------------------:|
|EListRaw |Raw Expression list. A class used to store single-channel raw intensities prior to normalization. Intensities are unlogged. Objects of this class contain one row for each probe and one column for each array. The function read.ilmn() for example creates an object of this class.|
|EList|Expression list. Contains background corrected and normalized log-intensities. Usually created from an EListRaw objecting using normalizeBetweenArrays() or neqc().|
|RGList|Red-Green list. A class used to store raw two-color intensities as they are read in from an image analysis output file, usually by read.maimages().|
|MAList|Two-color intensities converted to M-values and A-values, i.e., to within-spot and wholespot contrasts on the log-scale. Usually created from an RGList using MA.RG() or normalizeWithinArrays(). Objects of this class contain one row for each spot. There may be more than one spot and therefore more than one row for each probe.|
|MArrayLM|MicroArray Linear Model. Store the result of fittting gene-wise linear models to the normalized intensities or log-ratios. Usually created by lmFit(). Objects of this class normally contain one row for each unique probe.|
|TestResults|Store the results of testing a set of contrasts equal to zero for each probe. Usually created by decideTests(). Objects of this class normally contain one row for each unique probe.|

4. Linear Models Overview
* `design matrix`: indicates in effect which RNA samples have been applied to each array.
* `contrast matrix`: specifies which comparisons you would like to make between the RNA samples.
* **philosophy**: 
> You have to start by fitting a linear model to your data which fully models the systematic part of your data. The model is specified by the design matrix. Each row of the design matrix corresponds to an array in your experiment and each column corresponds to a coefficient that is used to describe the RNA sources in your experiment. With Affymetrix or single-channel data, or with two-color with a common reference, you will need as many coefficients as you have distinct RNA sources, no more and no less. With direct-design twocolor data you will need one fewer coefficient than you have distinct RNA sources, unless you wish to estimate a dye-effect for each gene, in which case the number of RNA sources and the number of coefficients will be the same. Any set of independent coefficients will do, providing they describe all your treatments. The main purpose of this step is to estimate the variability in the data, hence the systematic part needs to be modelled so it can be distinguished from random variation.

5. Single-Channel Experimental Design
> refer to limma(User's Guide) chapter9

6. Statistics for Differential Expression
> moderated t-statistic: standard errors have been moderated across genes, i.e., squeezed towards a common value, using a simple Bayesian model.
> multiple testing ajustment: "BH"(Benjamini and Hochberg's method)
> refer to limma(User's Guide) chapter13

7. Array Quality Weights
`arrayWeights()` function
> Array weights are generally useful when there is some reason to expect variable array quality. <human in vivo data>
> refer to limma(User's Guide) chapter14

8. RNA-seq Data
> refer to limma(User's Guide) chapter15

9. Single-Channel Case Studies
* Effect of Estrogen on Breast Cancer Tumor Cells: A 2x2 Factorial Experiment with Affymetrix Arrays
> packages needed； limma, affay, estrogen, hgu95av2cdf, annotate, hgu95av2.db
> The data gives results from a 2x2 factorial experiment on MCF7 breast cancer cells using Affymetrix HGU95av2 arrays. The factors in this experiment were estrogen (present or absent) and length of exposure (10 or 48 hours). The aim of the study is the identify genes which respond to estrogen and to classify these into early and late responders. Genes which respond early are putative direct-target genes while those which respond late are probably downstream targets in the molecular pathway.

```
	library(limma)
	library(affy)
	library(hgu95av2cdf)
	datadir <- file.path(find.package("estrogen"),"extdata")
	dir(datadir)
	targets <- readTargets("phenoData.txt",path=datadir,sep="",row.names="filename")
	View(targets)
	ab <- ReadAffy(filenames=targets$filename, celfile.path=datadir)
	eset <- rma(ab)
	library(annotate)
	library(hgu95av2.db)
	ID <- featureNames(eset)
	Symbol <- getSYMBOL(ID,"hgu95av2.db")
	fData(eset) <- data.frame(Symbol=Symbol)
	# Given that we are interested in the early and late estrogen responders, we can choose a parametrization which includes these two contrasts.
	treatments <- factor(c(1,1,2,2,3,3,4,4),labels=c("e10","E10","e48","E48"))
	contrasts(treatments) <- cbind(Time=c(0,0,1,1),E10=c(0,1,0,0),E48=c(0,0,0,1))
	design <- model.matrix(~treatments)
	colnames(design) <- c("Intercept","Time","E10","E48")
	# The second coefficient picks up the effect of time in the absence of estrogen. The third and fourth coefficients estimate the log2-fold change for estrogen at 10 hours and 48 hours respectively.
	fit <- lmFit(eset,design)
	# We are only interested in the estrogen effects, so we choose a contrast matrix which picks these two coefficients out:
	cont.matrix <- cbind(E10=c(0,1,0,0),E48=c(0,0,0,1)) # not agree with the original
	fit2 <- contrasts.fit(fit, cont.matrix)
	fit2 <- eBayes(fit2)
	# One way to decide which changes are significant for each gene would be to use Benjamini and Hochberg's method to control the false discovery rate across all the genes and both tests:
	results <- decideTests(fit2, method = "global")
	# Another method would be to adjust the F-test p-values rather than the t-test p-values:
	results <- decideTests(fit2, method = "nestedF")
	results <- classifyTestsF(fit2, p.value=0.0001)
```
> refer to limma(User's Guide) chapter17 P100-P105

* Comparing Mammary Progenitor Cell Populations with Illumina BeadChips
> This case study examines the expression profiles of adult mammary stem cells and of progenitor and mature mammary lumina cells. The data was first published on <Aberrant luminal progenitors as the candidate target population for basal tumor development in BRCA1 mutation carriers>, which used the expression profiles to show that limina progenitor cells are the likely cell of origina for basal-like breast cancer.The data files used in this case study can be downloaded from http://bioinf.wehi.edu.au/marray/IlluminaCaseStudy.
> Breast tissue was obtained from three healthy human donors who were undergoing reduction mammoplasties. Epithelial cells were sorted into three subpopulations enriched for mammary stem cells (MS), luminal progenitor cells (LP) and mature luminal cells (ML) [19]. The MS, LP and ML cells representative a lineage of luminal cells use to construct the ducts used to transport milk in the breast:
![](http://i.imgur.com/B1ogLPV.png)
> Stromal cells were also profiles as a comparison group. There were therefore four cell populations from each person.
> The file `probe profile.txt` contains the expression profiles for regular probes, designed to interrogate the expression levels of genes. `control probe profile.txt` contains the profiles of control probes, including negative control probes.

```
	dir()
	[1]"control probe profile.txt"
	[2] "probe profile.txt"
	[3] "Targets.txt"

	library(limma)
	targets <- readTargets("Targets.txt")
	targets <- targets[,c(3,1,2)]
	colnames(targets) <- c("Donor", "Age", "CellType")
	# We read in the expression profiles for both regular and control probes, telling read.ilm that we wish to read the detection p-values as well as the expression values:
	x <- read.ilmn(files="probe profile.txt",ctrlfiles="control probe profile.txt", other.columns="Detection")
	# This reads a EListRaw object. There are about 750 negative probes and about 49,000 regular probes:
	table(x$genes$Status)
	options(digits=3)
	head(x$E)
	boxplot(log2(x$E),range=0,ylab="log2 intensity")
	# The detection values contain p-values for testing whether each probe is more intense than the negative control probes. Small values are evidence that the probe corresponds to a truly expressed gene:
	head(x$other$Detection)
	# We can go further than this and estimate the overall proportion of the regular probes that correspond	to expressed transcript, using the method in <Estimating the proportion of microarray probes expressed in an RNA sample>.
	pe <- propexpr(x)
	dim(pe) <- c(4,3)
	dimnames(pe) <- list(CellType=c("MS","Stroma","ML","LP"),Donor=c(1,2,3))
	# background correction and normalization:
	y <- neqc(x) 
	# The neqc functions performs normexp background correction using negative controls, then quantile normalizes and finally log2 transforms.
	# It also automatically removes the control probes, leaving only the regular probes in y: dim(y) <48803	12>
	
	# Filter out probes that are not expressed. We keep probes that are expressed in at least three arrays according to a detection p-values of 5%:
	expressed <- rowSums(y$other$Detection < 0.05) >= 3
	y <- y[expressed,] # dim(y) <24691	12>
	
	# A multi-dimensional scaling plot shows that the cell types are cell separated:
	plotMDS(y,labels=targets$CellType)

	# Arrays from the same donor are not independent, so we need to estimate the within-dinor correlation:
	ct <- factor(targets$CellType)
	design <- model.matrix(~0+ct)
	colnames(design) <- levels(ct)
	dupcor <- duplicateCorrelation(y,design,block=targets$Donor)
	dupcor$consensus.correlation
	
	# We make all possible pairwise comparisons between the epithelial cell types, allowing for the correlation within donors:
	fit <- lmFit(y,design,block=targets$Donor,correlation=dupcor$consensus.correlation)
	contrasts <- makeContrasts(mL-MS, pL-MS, mL-pL, levels=design)
	fit2 <- contrasts.fit(fit, contrasts)
	fit2 <- eBayes(fit2, trend=TRUE)
	summary(decideTests(fit2, method="global"))
	
	# Now we find genes uniquely expressed in LP cells, as compared to MS and ML. We refit the linear model, making LP the reference cell type:
	ct <- relevel(ct, ref="pL")
	design <- model.matrix(~ct)
	fit <- lmFit(y,design,block=targets$Donor,correlation=dupcor$consensus.correlation)
	fit2 <- fit[,c("ctMS","ctmL")]
	fit2 <- eBayes(fit2, trend=TRUE)
	# The we find all those genes that are up-regulated in LP vs both MS and ML, using a 2-fold-change and 5% FDR:
	results <- decideTests(fit2, lfc=1)
	vennDiagram(results, include=c("up","down"))

	LP.sig <- rowSums(results>0)==2
	topTable(fit2[LP.sig,])
```
![](http://i.imgur.com/Ns2veVJ.png)

* Time Course Effects of Corn Oil on Rat Thymus with Agilent 4x44K Arrays
> This case study analyses a time-course experiment using single-channel Agilent Whole Rat Genome Microarray 4x44K v3 arrays.
> The experiment concerns the effect of corn oil on gene expression in the thymus of rats. The data was submitted by Hong Weiguo to ArrayExpress as series E-GEOD-33005.
> All files were downloaded from http://www.ebi.ac.uk/arrayexpress/experiments/E-GEOD-33005.

```
	SDRF <- read.delim("E-GEOD-33005.sdrf.txt",check.names=FALSE,stringsAsFactors=FALSE)
	x <- read.maimages(SDRF[,"Array Data File"],source="agilent",green.only=TRUE)
	y <- backgroundCorrect(x,method="normexp")
	y <- normalizeBetweenArrays(y,method="quantile")
	neg95 <- apply(y$E[y$genes$ControlType==-1,],2,function(x) quantile(x,p=0.95))
	cutoff <- matrix(1.1*neg95,nrow(y),ncol(y),byrow=TRUE)
	isexpr <- rowSums(y$E > cutoff) >= 4
	table(isexpr)
	y0 <- y[y$genes$ControlType==0 & isexpr,]
	> Treatment <- SDRF[,"Characteristics[treatment]"]
	levels <- c("10 ml/kg saline","2 ml/kg corn oil","5 ml/kg corn oil","10 ml/kg corn oil")
	Treatment <- factor(Treatment,levels=levels)
	design <- model.matrix(~Treatment)
	fit <- lmFit(y0,design)
	fit <- eBayes(fit,trend=TRUE)
	plotSA(fit, main="Probe-level")
	summary(decideTests(fit[,-1]))
	yave <- avereps(y0,ID=y0$genes[,"SystematicName"])
	fit <- lmFit(yave,design)
	fit <- eBayes(fit,trend=TRUE)
	plotSA(fit, main="Gene-level")
	summary(decideTests(fit[,-1]))
```

10. RNA-Seq Case Studies
> refer to limma(User's Guide) chapter18

11. Output function
`topTable()` set parameter **number** to be **Inf** toget all results 
`write.fit()`