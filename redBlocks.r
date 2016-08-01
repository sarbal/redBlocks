# Helper functions 
make_net_mat <- function(data, n, genes, filt){
	net = diag(n) * 0
        rownames(net) = genes
	colnames(net) = genes
        net[filt] = data
	net = (net + t(net))
	return(net)
}

calc_cpm <-function(X){
	K  = colSums(X)
	X.cpm = sapply(1:length(K), function(k) 10^6*X[,k]/K[k] )
	return(X.cpm)
}

run_edgeR <- function(X, f.a, filt, group){
    y <- DGEList(counts=X[f.a,filt], group=group[filt])
    y <- estimateGLMCommonDisp(y)
    design <- model.matrix(~group[filt])
    fit <- glmFit(y, design)
    degs <- glmLRT(fit, coef=2)
    pvals = p.adjust(degs$table[,4])
    degs = list(degs, pvals)
    return(degs)
}

run_edgeR_filter <- function(X, f.a, filt, group){
    y <- DGEList(counts=X[,filt], group=group[filt])
    keep <- rowSums(cpm(y)>1) >= 2
    yk = y[keep, keep.lib.sizes=FALSE]
    yk <- estimateGLMCommonDisp(yk)
    design <- model.matrix(~group[filt])
    fit <- glmFit(yk, design)
    degs <- glmLRT(fit, coef=2)
    pvals = p.adjust(degs$table[,4])
    degs = list(degs, pvals)
    return(degs)
}


run_DESeq2 <- function(X, f.a, filt, group){
    conditions = group[filt] + 1
    samples = colnames(X)[filt]
    colData = as.data.frame(cbind( samples, conditions ) )
    colnames(colData) = c("samples", "conditions")
    dds = DESeqDataSetFromMatrix(countData=X[f.a,filt], colData=colData, design=~conditions)
    dds = DESeq(dds)
    res = results(dds, contrast=c("conditions", "2", "1"))
    return(res)
}



run_DESeq2_filter <- function(X, f.a, filt, group){
    conditions = group[filt] + 1
    samples = colnames(X)[filt]
    keep <- rowSums(cpm(X)>1) >= 2
    colData = as.data.frame(cbind( samples, conditions ) )
    colnames(colData) = c("samples", "conditions")
    dds = DESeqDataSetFromMatrix(countData=X[keep,filt], colData=colData, design=~conditions)
    dds = DESeq(dds)
    res = results(dds, contrast=c("conditions", "2", "1"))
    return(res)
}


calc_DE <- function(X, f.a, filt, group){
	X = X[f.a,filt]
	group = group[filt]
  	if( sum(group==1) < 2  ) {
           	m.X1 = (X[,group==1])
	} else {
		m.X1 = rowMeans(X[,group==1])
 	}
  	if( sum(group==2) < 2  ) {
           	m.X2 = (X[,group==2])
	} else {
		m.X2 = rowMeans(X[,group==2])
 	}
        m.X = rowMeans(X)
	fc = log2(m.X1/m.X2)
	X.ps = sapply(1:dim(X)[1], function(k) wilcox.test(X[k,group==1], X[k,group==2])$p.val)
	X.padj = p.adjust(X.ps , method = "BH")
	de = cbind(m.X, fc, X.ps, X.padj, m.X1, m.X2)
	return(de)
}

calc_DE_filter <- function(X, f.a, filt,group){
        keep <- rowSums(X>1) >= 2
	X = X[keep,filt]
	group = group[filt]
  	if( sum(group==1) < 2  ) {
           	m.X1 = (X[,group==1])
	} else {
		m.X1 = rowMeans(X[,group==1])
 	}
  	if( sum(group==2) < 2  ) {
           	m.X2 = (X[,group==2])
	} else {
		m.X2 = rowMeans(X[,group==2])
 	}
        m.X = rowMeans(X)
	fc = log2(m.X1/m.X2)
	X.ps = sapply(1:dim(X)[1], function(k) wilcox.test(X[k,group==1], X[k,group==2])$p.val)
	X.padj = p.adjust(X.ps , method = "BH")
	de = cbind(m.X, fc, X.ps, X.padj, m.X1, m.X2)
	return(de)
}



plot_deg_coexp <- function(deg, flag, runid, filtMin=5, freq.net){
	if(flag == "wilcox"){
		m =  match(rownames(deg), attr$ensemblID[f.a])
		f.m = !is.na(m)
		f.am = m[f.m]
		sub.Y = attr$chr[f.a][f.am] == "chrY"
		sub.X = attr$chr[f.a][f.am] == "chrX"
		xist =  which(rownames(deg)[f.m] == "ENSG00000229807")
		m.X  = log2(deg[,1])[f.m]
		fc   = deg[,2][f.m]
		padj = -log10(deg[,4])[f.m]


	} else if(flag == "deseq"){
		m =  match(rownames(deg), attr$ensemblID[f.a])
		f.m = !is.na(m)
		f.am = m[f.m]
		sub.Y = attr$chr[f.a][f.am] == "chrY"
		sub.X = attr$chr[f.a][f.am] == "chrX"
		xist =  which(rownames(deg)[f.m] == "ENSG00000229807")
		m.X  = log2(1+deg$baseMean)[f.m]
		fc   = deg$log2FoldChange[f.m]
		fc[is.na(fc) ]  = 0
		padj = -log10(deg$padj)[f.m]

	} else if(flag == "edgeR"){
		m =  match(rownames(deg[[1]]$table), attr$ensemblID[f.a])
		f.m = !is.na(m)
		f.am = m[f.m]
		sub.Y = attr$chr[f.a][f.am] == "chrY"
		sub.X = attr$chr[f.a][f.am] == "chrX"
		xist =  which(rownames(deg[[1]]$table)[f.m] == "ENSG00000229807")
		m.X  = deg[[1]]$table$logCPM[f.m]
		fc   = deg[[1]]$table$logFC[f.m]
		padj = -log10(deg[[2]])[f.m]

        } else {
		return()
	}
	n.coexp=list()
        n.filt = list()

        deg.top= list()
        fc.top=list()
	
	deg.top[[1]] = tail(order(fc), n=100)
	deg.top[[2]] = head(order(fc), n=100)
        
	print(deg.top[[1]])
        print(deg.top[[2]])
	
	fc.top[[1]] = cbind(m.X, fc, padj)[deg.top[[1]],]
        fc.top[[2]] = cbind(m.X, fc, padj)[deg.top[[2]],]

	temp.net = freq.net[f.m,f.m]

	n.coexp[[1]]= temp.net[deg.top[[1]], deg.top[[1]] ]
        n.filt[[1]] = plot_and_filter_coexp(n.coexp[[1]], paste("up", runid),filtMin)

        n.coexp[[2]]= temp.net[ deg.top[[2]] , deg.top[[2]] ]
        n.filt[[2]] = plot_and_filter_coexp(n.coexp[[2]], paste("down", runid),filtMin )

	return(list(deg.top, n.filt, fc.top))
}

plot_and_filter_coexp <- function(n.coexp, runid, filtMin=5, median=0.28, use_median=TRUE, plot_results=TRUE){


	# Use median
	temp = n.coexp
	temp[temp > median] = 1
	temp[temp <= median] = 0


	# Cluster
	if( use_median ){
		consTree = hclust(as.dist(temp), method = "average");
		consDend = as.dendrogram(consTree)
	        unmergedLabels3 = cutreeDynamic(dendro = consTree, distM = temp,
	                deepSplit = 2, cutHeight = 0.995,
	                minClusterSize = 2,
	                pamRespectsDendro = FALSE );
	 } else {
	 	consTree = hclust(as.dist(n.coexp), method = "average");
                consDend = as.dendrogram(consTree)
	        unmergedLabels3 = cutreeDynamic(dendro = consTree, distM = n.coexp,
	                deepSplit = 2, cutHeight = 0.995,
	                minClusterSize = 2,
	                pamRespectsDendro = FALSE );
	 }

	unmergedColors3 = cols2[ as.numeric(unmergedLabels3)+1]


	# Label and count modules
	i.prev = ""
	ki = 1
	ji = 0
	unmergedLabels.mod = as.numeric(unmergedLabels3[consTree$order]) * 0
	for( ii in as.numeric(unmergedLabels3[consTree$order]) ){
	   if( ii == i.prev){
	       unmergedLabels.mod[ki] = ji
	    } else {
	        i.prev = ii
	        ji = ji + 1
	        unmergedLabels.mod[ki] = ji

	    }
	    ki = ki + 1
	}

	# Filter modules
	f.freq = count(unmergedLabels.mod)[,2] < filtMin
	f.keep = count(unmergedLabels.mod)[f.freq,1]
	unmergedLabels.mod2 = unmergedLabels.mod[order(consTree$order)]
	m = match( unmergedLabels.mod2, f.keep)
	f.fm = !is.na(m)
	
	# Plots
	if(plot_results) {
		# Plot input
		heatmap.2(n.coexp, density="none", trace="none", main=runid )
		# Plot output of clustering and module detection
		heatmap.2(temp, density="none", trace="none", col=cols, Rowv=consDend, Colv=consDend, RowSideColors=unmergedColors3, ColSideColors=unmergedColors3,cexRow=0.5, cexCol=0.5, main=runid)
        	heatmap.2(n.coexp, density="none", trace="none", col=cols, Rowv=consDend, Colv=consDend, RowSideColors=unmergedColors3, ColSideColors=unmergedColors3,cexRow=0.5, cexCol=0.5, main=runid)
		# Plot filter
		heatmap.2(n.coexp[f.fm,f.fm], density="none", trace="none", cexRow=0.5, cexCol=0.5, main=runid)
	}
	# Return filter
	return(f.fm)
}



