# To load network
load("freq.75.net.Rdata")
freq.net = make_net_mat(data,n,genes[,2], filt)

# To run DE analysis
load("counts_data.Rdata")
X = counts
X.cpm = calc_cpm(X)

# Data filters
# row_filter - already
col_filter = c(prb, fam) # example, selecting all the probands and all the parents
ftr = famf
mtr = famm 

col_filter = labels[c(prb, fam),2]==1  # example, selecting  family 1 
ftr = labels[famf,2]==1
mtr = labels[famm,2]==1 

# Set up conditions here (1 - unaffected, 2 - affected)
group = as.numeric(labels[,6] == 1 ) + 1    

# Run DE 
deg.edgeR.filt = run_edgeR_filter(X, row_filter, col_filter ,group)
deg.deseq.filt = run_DESeq2_filter(X,row_filter, col_filter,group)
deg.wilcox.filt = calc_DE_filt(X.cpm+1,row_filter, col_filter,2-group, ftr, mtr, 100)

# Filter DE results for common co-expression
edgeR.results = plot_deg_coexp(deg.edgeR.filt, "edgeR", runid, filtMin=5, freq.net)
wilcox.results = plot_deg_coexp(deg.wilcox.filt, "wilcox", runid, filtMin=5, freq.net)
deseq.results = plot_deg_coexp(deg.deseq.filt, "deseq", runid, filtMin=5, freq.net)
