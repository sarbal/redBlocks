# Figures for manuscript
These are the scripts to draw the manuscript figures. 
Some additional analyses are interspersed throughout but most are pre-run and stored in their respective directories as indicated. 

## Set up environment 
source("outliers/bin/helper_redBlocks.r")

## Figure 3 
![fig3](https://github.com/sarbal/redBlocks/blob/master/outliers/imgs/figure3_genes.png "Figure 3")
Caption: Disease expression analysis with a family-based approach. (A) The expression fold change for each gene is calculated within each family (top 100 up and down regulated genes are shown). (B) Overlaps in DE gene sets between the individual families (numbers in boxes), and the significance of this overlap (colored corresponding to -log10 P-value of the hypergeometric test). Overlaps are mostly small. (C) The replicable genes are those that are recurrent across families. The recurrence distributions for both up- and downregulated genes across the 6 families are shown. Using the binomial test, we find that genes recurring 3 or more times are significant (FDR<0.05). These genes are listed, with 4 up- and 14 downregulated genes significantly recurrent. (D) Robustness assessment of the DE threshold. The plot shows the number of recurrent genes as a function of the number of differentially expressed genes and the significance of the recurrence in grey.

#### Panel A
```{r}
load("outliers/figures/fig3.panelA.Rdata")
pdf("outliers/figures/fig3.panelA.gene.recur.pdf")
heatmap.3(up.plot, Colv=F, col=rev(viridis(100)), Rowv=F, RowSideCol=rev(magma(6))[up.bar], main="Up")
heatmap.3(down.plot, Colv=F, col=rev(viridis(100)), Rowv=F, RowSideCol=rev(magma(6))[down.bar], main="Down")
leg2 = cbind(rev(magma(6)), 1:6)
plot(0,1)
legend( "bottom", leg=leg2[,2], col=leg2[,1], lwd=3 )
dev.off()

# To reproduce analysis:
load("outliers/taf1_pedigrees/summary.deg.Rdata")
fcc = c(7:9,12,11,10)
all = unique(array(sapply(fcc, function(i) head(order(genesets.fcranks[,i]), n=100 ))  ))
freq = count(array(sapply(fcc, function(i) head(order(genesets.fcranks[,i]), n=100 ))  ))
m = match( all, freq[,1])
om = order(-freq[m,2])
down.plot = genesets.fcranks[all,fcc][om,]/dim(genesets.fcranks)[1]
down.bar = freq[m,2][om]

freq = count(array(sapply(fcc, function(i) tail(order(genesets.fcranks[,i]), n=100 ))  ))
all = unique(array(sapply(fcc, function(i) tail(order(genesets.fcranks[,i]), n=100 ))  ))
m = match( all, freq[,1])
om = order(-freq[m,2])
up.plot = 1- (genesets.fcranks[all,fcc][om,]/dim(genesets.fcranks)[1])
up.bar = freq[m,2][om]

leg2 = cbind(rev(magma(6)), 1:6)
save( up.plot, down.plot, up.bar, down.bar, leg2, file="outliers/figures/fig3.panelA.Rdata")
```

#### Panel B
```{r}
load("outliers/figures/fig3.panelB.Rdata")
pdf("outliers/figures/fig3.panelB.common.overlaps.pdf")
heatmap.3(mat.comb, col=cols, Rowv=F, Colv=F, cellnote=sigtemp2, notecol="black", notecex=2 )
dev.off()

# To reproduce analysis:
load("outliers/taf1_pedigrees/summary.deg.Rdata")
load("outliers/taf1_pedigrees/taf1.DE.Rdata")
reorder = c(1:3,6,5,4)
down100overlaps = sapply(reorder, function(i) sapply(reorder, function(j) length(intersect(which(genes.down[,i]>0), which(genes.down[,j]>0) )) ) )
up100overlaps = sapply(reorder, function(i) sapply(reorder, function(j)  length(intersect(which(genes.up[,i]>0), which(genes.up[,j]>0) )) ) )

up100overlaps.pvals = up100overlaps * 0
down100overlaps.pvals = down100overlaps * 0

for(i in 1:6){
	for(j in 1:6){

		q = up100overlaps[i,j] - 1
		m = 100
		n = dim(genesets.fcranks)[1] - m
		k = 100
		up100overlaps.pvals[i,j] = phyper(q, m, n , k, lower.tail = FALSE, log.p = FALSE)

		q = down100overlaps[i,j] -1
		m = 100
		n = dim(genesets.fcranks)[1] - m
		k = 100
		down100overlaps.pvals[i,j] = phyper(q, m, n , k, lower.tail = FALSE, log.p = FALSE)
	}
}

bot = row(down100overlaps.pvals) > col( down100overlaps.pvals )
top  = row(down100overlaps.pvals) <  col( down100overlaps.pvals )
mat.down = -log10(down100overlaps.pvals)
mat.up = -log10(up100overlaps.pvals)
mat.comb = mat.down * 0
mat.comb[bot] =  -log10(p.adjust(down100overlaps.pvals[bot]))
mat.comb[top] =  -log10(p.adjust(up100overlaps.pvals[top]))
count.comb =  up100overlaps
count.comb[bot] = down100overlaps[bot]
diag(count.comb) = ""
sigtemp = (mat.comb > -log10(0.05) ) * 1
sigtemp[sigtemp==1] = "*"
sigtemp[sigtemp==0] = ""
sigtemp2 = matrix(paste(count.comb, sigtemp), ncol=6 )

save(mat.comb, sigtemp2, file="outliers/figures/fig3.panelB.Rdata")
```

#### Panel C
```{r} 
load("outliers/figures/fig3.panelC.Rdata")
pdf("outliers/figures/fig3.panelC.hist.pdf")
hist( recur.u[recur.u>0]+0.01 , col=cols.rec[2], main="Recurrence, up", xlab="Count",xlim=c(1,4))
abline( v = fdrs.up$Pt)
hist( recur.d[recur.d>0]+0.01 , col=cols.rec[3], main="Recurrence, down", xlab="Count",xlim=c(1,4))
abline( v = fdrs.down$Pt)
dev.off()

# To reproduce analysis:
load("outliers/analysis/fdr_calcs.genes.up.Rdata")
fdrs.up = fdrs
load("outliers/analysis/fdr_calcs.genes.down.Rdata")
fdrs.down= fdrs
load("taf1.DE.Rdata")
recur.d = rowSums(genes.down)
recur.u = rowSums(genes.up)
cols.rec  =c(magma(5)[2], "deeppink4", "darkcyan" )
save(cols.rec, recur.u ,recur.d, fdrs.down, fdrs.up, file="outliers/figures/fig3.panelC.Rdata")
```

#### Panel D
```{r} 
load("outliers/figures/fig3.panelD.Rdata")
pdf("outliers/figures/fig3.panelD.pdf")
plot_recurrence(xrange, n.sig.rec.up, p.tests.up, main="Up", xlab="Number of DE genes", ylab="Number of genes significantly recurrent")
plot_recurrence(xrange, n.sig.rec.down, p.tests.down,main="Down", xlab="Number of DE genes", ylab="Number of genes significantly recurrent")
dev.off()

# To reproduce analyses:
plot_recurrence <- function(xrange, n.sig, p.tests, ...){
      plot( xrange , n.sig , type="l", lwd=2, axes=F, ...)
      axis(1); axis(2);
      par(new=TRUE)
      plot( xrange , -log10(p.tests) , type="l",   axes=F, col="lightgrey", ylab="", xlab="")
      axis(4); mtext("-log10 adjusted P-value of significance threshold",4)
}

load("outliers/taf1_pedigrees/p.recur.pre.Rdata")
xrange = 10:1000
p.max = 0.05
n.sig.rec.up = sapply(xrange-9, function(i)  sum(p.recur.up[[i]][p.recur.up[[i]][,2]  <= p.max,3]) )
n.sig.rec.down = sapply(xrange-9, function(i)  sum(p.recur.down[[i]][p.recur.down[[i]][,2]  <= p.max,3]) )
p.tests.up = sapply(xrange-9, function(i)   p.recur.up[[i]][p.recur.up[[i]][,2]  <= p.max,2][1] )
p.tests.down = sapply(xrange-9, function(i) p.recur.down[[i]][p.recur.down[[i]][,2]  <= p.max,2][1])

save(p.tests.down, p.tests.up,n.sig.rec.down, n.sig.rec.up ,p.max ,xrange , plot_recurrence, file="outliers/figures/fig3.panelD.Rdata")
```

## Figure 4 
![fig4](https://github.com/sarbal/redBlocks/blob/master/outliers/imgs/figure4_pathways.png "Figure 4")
Caption: Gene set enrichment assessment of TAF1 cohort. (A) Top GO enrichment results for each family for up-, and downregulated genes. Significant terms (FDR < 0.05) are highlighted with an asterix. (B) The frequency and significance of recurrence of each GO term is plotted for the upregulated genes. (C) Gene-GO membership matrix for the upregulated genes. Each column is a gene, and each row a family. The colored bars below highlight the GO terms that these genes belong to. The signal associated with the recurrent GO terms is distributed across different genes, shown by low overlap across the families. (D) There are no significantly recurrent pathways with the downregulated genes. (E) However, the recurrent genes themselves are enriched for ribosomal pathways (p-adjusted<0.05), as shown in the gene-GO membership matrix. The three RP* genes seem to drive almost all the signal.
#### Panel A
```{r}
load("outliers/figures/fig4.panelA.Rdata")
pdf("outliers/figures/fig4.panelA.pdf")
plot_enrich_res(padj.up, col.size)
plot_enrich_res(padj.down, col.size)
dev.off()
```

#### Panel B
```{r}
load("outliers/figures/fig4.panelB.Rdata")
pdf("outliers/figures/fig4.panelB.pdf")
hist( recur.up.path[recur.up.path>0]+0.01 , col=cols.rec[2], main="Recurrence, up", xlab="Recurrence of GO group",xlim=c(1,4))
abline( v = Pt.up, lwd=2, col="lightgrey")
dev.off()
```

#### Panel C
```{r}
load("outliers/figures/fig4.panelC.Rdata")
pdf("outliers/figures/fig4.panelC.pdf")
pheatmap((gene.mat2[o,fo]), col=rev(magma(2)), cluster_rows=F, cluster_cols=F, annotation_row=row_df, border=NA )
dev.off()
```


#### Panel D
```{r}
load("outliers/figures/fig4.panelD.Rdata")
pdf("outliers/figures/fig4.panelD.pdf")
hist( recur.down.path[recur.down.path>0]+0.01 , col=cols.rec[3], main="Recurrence, down", xlab="Recurrence of GO group",xlim=c(1,4))
abline( v = Pt.down, lwd=2, col="lightgrey")
dev.off()
```

#### Panel E
```{r}
load("outliers/figures/fig4.panelE.Rdata")
pdf("outliers/figures/fig4.panelE.pdf")
pheatmap((gene.mat2[o1,][f1,fo]), col=rev(magma(2)), cluster_rows=F, cluster_cols=F, annotation_row=row_df[f1,], border=NA )
dev.off()
```

```{r}
# To reproduce analysis:
load("outliers/taf1_pedigrees/analysis/fdr_calcs.goslim.paths.all.Rdata")
load("outliers/taf1_pedigrees/analysis/summary.goslim.paths.Rdata")
load("outliers/taf1_pedigrees/analysis/summary.goslim.paths.recurs.Rdata")
load("outliers/taf1_pedigrees/analysis/goslim.summary2genes.Rdata")
load("outliers/taf1_pedigrees/analysis/goslim.enrich.all.Rdata")

fcc = c(1:3,6,5,4)
col.size = go.slim[ff,3]
padj.up = apply(pvals.up[ff,fcc],2, p.adjust)
rownames(padj.up) = go.slim[ff,1]
padj.down = apply(pvals.down[ff,fcc],2, p.adjust)
rownames(padj.down) = go.slim[ff,1]
save(padj.down, padj.up,plot_enrich_res,col.size,ff,go.slim, cols4b, file="outliers/figures/fig4.panelA.Rdata")

recur.up.path = recurs.paths[,2]
Pt.up = pt.p.all[[2]]
save(recur.up.path,  Pt.up , cols.rec, file="outliers/figures/fig4.panelB.Rdata")
recur.down.path = recurs.paths[,3]
Pt.down = pt.p.all[[3]]
save(recur.down.path,  Pt.down , cols.rec, file="outliers/figures/fig4.panelD.Rdata")

plot_enrich_res <- function(padj.res, col.size,  p=0.05  ){
   f = rowSums( 1* (padj.res < p )) > 0
   sigtemp = (padj.res < p) * 1
   sigtemp[sigtemp==1] = "*"
   sigtemp[sigtemp==0] = ""
   recur.p = rowSums(padj.res < p)
   recur.col = c(0,rev(magma(6)))[recur.p+1]
   cols4b = colorpanel(5000, "lightgrey", "blue", "darkblue")
   o = order(go.slim[ff,3][f] )
   heatmap.3( -log10(padj.res[f,][o,]), Colv=F, Rowv=F, RowSideColors = recur.col[f][o], col=cols9,cexRow = 0.7, cellnote=sigtemp[f,][o,], notecol="black", notecex=2 )
   heatmap.3( -log10(padj.res[f,][o,]), Colv=F, Rowv=F, RowSideColors = cols4b[col.size][f][o], col=cols9,cexRow = 0.7, cellnote=sigtemp[f,][o,], notecol="black", notecex=2 )
}


# Cleaning up results
fo = c(1:3,6,5,4)
load("temp2slim.Rdata")
library(pheatmap)
f11 = order(gene.mat2[,1] ,gene.mat2[,2], gene.mat2[,3], gene.mat2[,6], gene.mat2[,5], gene.mat2[,4] )
gene.mat2 = gene.mat2[f11,]
gene.mat = gene.mat[f11,]
m = match(rownames(gene.mat), rownames(annot.sub) )
o1 = order(annot.sub[m,1] )
o2 = order(annot.sub[m,2] )
o3 = order(annot.sub[m,3] )
o = o1
rec2 = rowSums(gene.mat2)
row_df = data.frame( cbind(annot.sub[m,][o,] + 1 , rec2[o]))
row_df[,1] = as.factor(row_df[,1])
row_df[,2] = as.factor(row_df[,2])
row_df[,3] = as.factor(row_df[,3])
row_df[,4] = as.factor(row_df[,4])

rownames(row_df) = rownames(gene.mat2[o,])
pheatmap((gene.mat2[o,fo]), col=rev(magma(2)), cluster_rows=F, cluster_cols=F, annotation_row=row_df, border=NA )
save( gene.mat2, o , fo, row_df, file="outliers/figures/fig4.panelC.Rdata")

load("recur4slim.Rdata")
m = match(rownames(gene.mat), rownames(annot.sub) )
o1 = order(annot.sub[m,1] )
o2 = order(annot.sub[m,2] )
rec2 = rowSums(gene.mat2[o1,])
f1 = rec2  >=3
row_df = data.frame(cbind(annot.sub[m,][o1,]  + 1 , rec2  ))
rownames(row_df) = rownames(gene.mat2[o1,])
row_df = apply(row_df, 2, as.factor)
row_df = as.data.frame(row_df)
pheatmap((gene.mat2[o1,][f1,]), col=rev(magma(2)), cluster_rows=F, cluster_cols=F, annotation_row=row_df[f1,], border=NA )
fo = c(1,2,3,6,5,4)
pheatmap((gene.mat2[o1,][f1,fo]), col=rev(magma(2)), cluster_rows=F, cluster_cols=F, annotation_row=row_df[f1,], border=NA )

save(gene.mat2, o1, f1, fo, row_df , file="outliers/figures/fig4.panelE.Rdata)
```



## Figure 5 
![fig5](https://github.com/sarbal/redBlocks/blob/master/outliers/imgs/figure5_coexpfilt.png "Figure 5")
Caption: Co-expression of differentially expressed genes generates enrichment.  (A) As an example from family 1, we show the co-expression frequency sub-network as a heatmap, where genes showing decreased expression show co-expression. Co-expression blocks define modules as determined by the clustering (see rows). The modules are enriched for particular genes, mainly ribonucleoproteins. Performing a gene set enrichment analysis on these genes (Fisher’s exact test on GO groups), genes (rows) that generate the enrichment (columns are enriched GO terms) almost exclusively overlap with the co-expression blocks. The prominent pathways are ribosome related. (B) The significantly recurrent genes can be divided into those present within co-expression modules (joint) and those not (disjoint). The genes in bold are the functional outliers and the venn diagrams summarizes the number of genes in each category. (C) If we look at the enrichment of these DE gene sets (pre-filtering dark line +/-SD shadow), we see that filtering off the modules removes all but a few significant terms (lighter line, +/-SD shadow).


#### Panel A
```{r}
# To plot panel:
load("outliers/figures/fig5.panelA2.Rdata")
load("outliers/figures/fig5.panelA.Rdata")
pdf("outliers/figures/fig5.panelA.heatmap.f1.pdf")
heatmap.2(n.coexp, col=viridis(100), density="none", trace="none", Rowv=consDend, Colv=consDend, RowSideColors=unmergedColors3, ColSideColors=unmergedColors3,cexRow=0.5, cexCol=0.5, main="Downregulated genes")
heatmap.3( cbind(go.figslim, go.msig), col=cols6, Rowv=consDend, ,Colv=F, RowSideColors=unmergedColors3, ColSideCol =cols4b[c(col.figslim, col.msig)], cexRow=0.5, cexCol=0.5, main="")
heatmap.3( cbind(go.figslim2, go.msig2), col=cols6, Rowv=consDend, ,Colv=F, RowSideColors=unmergedColors3, ColSideCol =cols4b[c(col.figslim, col.msig)], cexRow=0.5, cexCol=0.5, main="")
dev.off()

# To reproduce analysis:
load("outliers/taf1_pedigrees/1adj.tally_75.subnet.down1000.Rdata")
n.coexp = subnet[1:100,1:100]
medK = 0.28
filtMin=6
temp = n.coexp
temp[temp >medK] = 1
temp[temp <= medK] = 0
consTree = hclust(as.dist(temp), method = "average");
consDend = as.dendrogram(consTree)
unmergedLabels3 = cutreeDynamic(dendro = consTree, distM = temp, deepSplit = 2, cutHeight = 0.995, minClusterSize = 2, pamRespectsDendro = FALSE );
nsclust = as.numeric(unmergedLabels3)+1
unmergedColors3 = magma( max(nsclust))[nsclust]
n.coexp = 1-n.coexp
save(n.coexp, consDend, unmergedColors3 , file="outliers/figures/fig5.panelA.Rdata")
```


#### Panel B
```{r}
load("outliers/taf1_pedigrees/analysis/summary.gene.recurs.int.Rdata")
run_recur_compare(recurs,pt.all,  "outliers/figures/fig5.panelB.Rdata")
load("outliers/figures/fig5.panelB.Rdata")
pdf("outliers/figures/fig5.panelB.pdf")
plot_2D_hist (up.mat10, Pt.pre.up, Pt.post.up, xlab="Recurrence pre-filtering", ylab="Recurrence post-filtering", col=cols11)
plot_2D_hist (down.mat10, Pt.pre.down, Pt.post.down, xlab="Recurrence pre-filtering", ylab="Recurrence post-filtering", col=cols11)
dev.off()
``` 
#### Panel C
```{r}
load("outliers/figures/fig5.panelC.Rdata")
pdf("outliers/figures/fig5.panelC.pdf")
plot_de_range(counts_sig.pre.down, counts_sig.post.down, cols.rec[3], ylim=c(0,15), xlim=c(10,200) )
plot_de_range(counts_sig.pre.up, counts_sig.post.up, cols.rec[2], ylim=c(0,15), xlim=c(10,200) )
dev.off()
```

```{r}
cols.rec  =c(magma(5)[2], "deeppink4", "darkcyan" )
plot_de_range <- function( counts.pre, counts.post, col, ...){
        #  yrange = range( cbind(colMeans(counts.pre)+colSD(counts.pre), colMeans(counts.pre)-colSD(counts.pre),
        #                        colMeans(counts.post)+colSD(counts.post), colMeans(counts.post)-colSD(counts.post) ))
          X_c = 10:1000
          Y_c = colMeans(counts.pre)
          std_Y_c = colSD(counts.pre)

          plot(X_c, Y_c, lwd=2, col=col,   type="l", ... )
          polygon(c(X_c, rev(X_c)), c(Y_c - std_Y_c, rev(Y_c + std_Y_c)), col = "lightgrey", border = NA)
          lines(X_c, Y_c, lwd=2, col=col)

          Y_c = colMeans(counts.post)
          std_Y_c = colSD(counts.post)
          polygon(c(X_c, rev(X_c)), c(Y_c - std_Y_c, rev(Y_c + std_Y_c)), col = makeTransparent(col), border = NA)
          lines(X_c, Y_c, lwd=2, col=col)

          

}
```

## Figure 6 
![fig6](https://github.com/sarbal/redBlocks/blob/master/outliers/imgs/figure6_meta_mod.png "Figure 6")
Caption: Differential expression meta-analysis in three other disorders. (A) Recurrence of genes in Huntington’s disease (HD), (B) Parkinson’s disease (PD) and (C) schizophrenia (SCZ), and whether they occur in groups (joint) or not (disjoint). The venn diagrams summarize the number of recurrent genes and their joint or disjoint designation.

#### Panel A
```{r}
library(tidyverse)
library(plyr)
load("outliers/figures/fig6.panelA.Rdata")
pdf("outliers/figures/fig6.panelA.pdf")
plot_2D_hist (up.mat10, Pt.pre.up, Pt.post.up, xlab="Recurrence pre-filtering", ylab="Recurrence post-filtering", col=cols11)
plot_2D_hist (down.mat10, Pt.pre.down, Pt.post.down, xlab="Recurrence pre-filtering", ylab="Recurrence post-filtering", col=cols11)
dev.off()
```
#### Panel B
```{r}
load("outliers/figures/fig6.panelB.Rdata")
pdf("outliers/figures/fig6.panelB.pdf")
plot_2D_hist (up.mat10, Pt.pre.up, Pt.post.up, xlab="Recurrence pre-filtering", ylab="Recurrence post-filtering", col=cols11)
plot_2D_hist (down.mat10, Pt.pre.down, Pt.post.down, xlab="Recurrence pre-filtering", ylab="Recurrence post-filtering", col=cols11)
dev.off()
```

#### Panel C
```{r}
load("outliers/figures/fig6.panelC.Rdata")
pdf("outliers/figures/fig6.panelC.pdf")
plot_2D_hist (up.mat10, Pt.pre.up, Pt.post.up, xlab="Recurrence pre-filtering", ylab="Recurrence post-filtering", col=cols11)
plot_2D_hist (down.mat10, Pt.pre.down, Pt.post.down, xlab="Recurrence pre-filtering", ylab="Recurrence post-filtering", col=cols11)
dev.off()

#### Venn insets 
```

```{r}
# To reproduce analysis:
load("outliers/other_diseases/HD/summary.gene.recurs.int.Rdata")
run_recur_compare(recurs,pt.all, "outliers/figures/fig6.panelA.Rdata")

load("outliers/other_diseases/PD/summary.gene.recurs.int.Rdata")
run_recur_compare(recurs, pt.all, "outliers/figures/fig6.panelB.Rdata")

load("outliers/other_diseases/SCZ/summary.gene.recurs.int.Rdata")
run_recur_compare(recurs,pt.all,  "outliers/figures/fig6.panelC.Rdata")

run_recur_compare <- function(recurs, pt.all, filename){
      require(plyr)
      require(tidyverse)
      # Set up upregulated
      temp1 = recurs[,c(4,3)]
      temp = plyr::count( temp1)
      temp.mat = spread(temp, key=x.recur.u.filt, value=freq)
      rownames(temp.mat) = temp.mat[,1]
      temp.mat = temp.mat[,-1]
      temp.mat[is.na(temp.mat)] = 0
      temp.mat = as.matrix(temp.mat )
      
      up.mat10 =   log10(temp.mat )  + 1
      up.mat10[!is.finite(up.mat10)] = 0
      Pt.pre.up = pt.all[3]
      Pt.post.up = pt.all[4]
      
      
      # Set up downregulated
      temp1 = recurs[,c(6,5)]
      temp = plyr::count( temp1)
      temp.mat = spread(temp, key=x.recur.d.filt, value=freq)
      rownames(temp.mat) = temp.mat[,1]
      temp.mat = temp.mat[,-1]
      temp.mat[is.na(temp.mat)] = 0
      temp.mat = as.matrix(temp.mat )
      
      down.mat10 =   log10(temp.mat )  + 1
      down.mat10[!is.finite(down.mat10)] = 0
      Pt.pre.down = pt.all[5]
      Pt.post.down = pt.all[6]
      
      save(up.mat10, Pt.pre.up, Pt.post.up, down.mat10, Pt.pre.down ,Pt.post.down, plot_2D_hist, file=filename)

}

plot_2D_hist <- function(mat,pt.x, pt.y, ...){
    image( mat, axes=F, ... )
    n=dim(mat)[2] -1
    ni = diff((0:n/n))[1]
    abline( h= (0:n/n)[pt.y]+ni/2, lwd=3, col="grey")
    axis(2, at=0:n/n, lab=0:n)
    n=dim( mat)[1] -1
    ni = diff((0:n/n))[1]
    axis(1, at=0:n/n, lab=(0:n)  )
    abline( v= (0:n/n)[pt.x]+ni/2, lwd=3, col="grey")
}
```
