# Recurrence analysis

Setup your environment. 
```{r echo==FALSE}
source("outliers/bin/helper_redBlocks.r")
```

Load your data. You will need either multiple expression experiments or multiple differentially expressed gene (DEGs) lists. In this example, we will load a set of DEGs from Huntington's disease studies. 
```{r echo=FALSE}
load("hd.DE.Rdata")
```
This file has two large matrices, genes by study, one with downregulated genes, and one with the upregulated genes. We analyze them independently. 

These are the five studies and the number of genes that were significantly upregulated in each.
```{r}
colSums(genes.up)
subgenesets = (genes.up[fg,])*1
```

To find gene-level recurrence, we simply sum the rows, and permute to calculate signficance (or use a bionomial test).  
```{r}
recur = rowSums(subgenesets)
fdrs = calc_fdrs_recur( subgenesets )
fdrs.bin = calc_binom_recur(subgenesets)
hist( recur[recur>0]+0.01 , col=cols.rec[3], main="Recurrence", xlab="Count",xlim=c(1,4))
abline( v = fdrs$Pt)
```
![step1](https://github.com/sarbal/redBlocks/blob/master/outliers/imgs/step1.png "gene-level recurrence")

Next, we look for pathway enrichment, and pathway-level recurrence. 
```{r}
library(EGAD)
load("goslim.Rdata")
n = dim(subgenesets)[2]
annotations = make_annotations(GO.human[,c(1,3)], (unique(GO.human[,1])), go.slim[ff,1])
go.enrich = lapply(1:n, function(i) gene_set_enrichment( names(which(subgenesets[,i]==1)), annotations, go.slim[ff,1:2]))
paths = sapply(1:n, function(i) (go.enrich[[i]]$padj<0.05)*1 )
paths.padj = sapply(1:n, function(i) (go.enrich[[i]]$padj) )
rownames(paths) = go.enrich[[1]][,1]
rownames(paths.padj) = go.enrich[[1]][,1]


f = rowSums( paths) > 1
sigtemp= paths
sigtemp[sigtemp==1] = "*"
sigtemp[sigtemp==0] = ""

heatmap.3( -log10(paths.padj[f,]), Colv=F, Rowv=F, col=cols9,cexRow = 0.7, cellnote=sigtemp[f,], notecol="black", notecex=2 )

```
![step2](https://github.com/sarbal/redBlocks/blob/master/outliers/imgs/step2.png "pathway enrichment")

And recurrence of the pathways. 
```{r}
fdrs.paths = calc_fdrs_recur(paths)
recur.path = rowSums(paths)
hist( recur.path[recur.path>0]+0.01 , col=cols.rec[3], main="Recurrence", xlab="Count",xlim=c(1,4))
abline( v = fdrs.paths$Pt)
```
![step3](https://github.com/sarbal/redBlocks/blob/master/outliers/imgs/step3.png "pathway-level recurrence")

```{r}
genes.sub = subgenesets[recur>0,]
m = match(rownames(genes.sub), rownames(annotations))
f.r = !is.na(m)
f.a = m[f.r]
go.sub = annotations[f.a,]
pheatmap(go.sub)
pheatmap(genes.sub[f.r,] )
pheatmap(paths.padj)
```
![step4a](https://github.com/sarbal/redBlocks/blob/master/outliers/imgs/step4a.png "GO-gene matrix")
![step4b](https://github.com/sarbal/redBlocks/blob/master/outliers/imgs/step4b.png "Gene-GO matrix")
![step4c](https://github.com/sarbal/redBlocks/blob/master/outliers/imgs/step4c.png "Pathway matrix")


We can look at the enrichment of the recurrent genes. 
```{r}
go.enrich.recur = lapply(1:n, function(i) gene_set_enrichment( names(recur[recur>(i-1)]), annotations, go.slim[ff,1:2]))
pathsrec = sapply(1:n, function(i) (go.enrich.recur[[i]]$padj<0.05)*1 )
pathsrec.padj = sapply(1:n, function(i) go.enrich.recur[[i]]$padj )
rownames(pathsrec) = go.enrich.recur[[1]][,1]
rownames(pathsrec.padj) = go.enrich.recur[[1]][,1]

heatmap.3( -log10(pathsrec.padj), Colv=F, Rowv=F, col=cols9,cexRow = 0.7, notecol="black", notecex=2 )
```
![step5](https://github.com/sarbal/redBlocks/blob/master/outliers/imgs/step5.png "recurrent gene enrichment")

For our last assessment, we look at co-expression as a secondary take at systems-level analyses. 
```{r}
load("gene_annotations_v22.Rdata")
load("freq.75.net.Rdata")
network = make_net_mat(data,n,genes[,2], filt)
nettype = "tally75"
# clean up 
rm(data)

# set up gene sets to extract network
# since genes.up and genes.down are from the same data the filtering indices are the same
res.up = run_filtering( subgenesets, "up", network, outputflag=T )

recur = res.up[[1]]
recur.filt = res.up[[2]]
fdrs = res.up[[3]]
fdrs.filt =res.up[[4]]
recurs =cbind(recur, recur.filt)
pre.post.mat = recur_mat(recurs)
plot_2D_hist(pre.post.mat, fdrs$Pt, fdrs.filt$Pt, col=cols11)

```
![step6](https://github.com/sarbal/redBlocks/blob/master/outliers/imgs/step6.png "pre/post")


Finally, we can classify our genes of interest based on the gene-level and pathways-level properties. 
```{r}
recurs[recurs[,2]>= fdrs$Pt,]
```


