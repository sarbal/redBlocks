# To load network
load("freq.75.net.Rdata")
freq.net = make_net_mat(data,n,genes[,2], filt)

# To run DE analysis
load("counts_data.Rdata")
X = counts
X.cpm = calc_cpm(X)

# Data filters
f.a = row_filter
file = col_filter

# Set up conditions here (1 - unaffected, 2 - affected)
group = conditions

# Run DE 
deg.edgeR.filt = run_edgeR_filter(X, f.a, filt,group)
deg.deseq.filt = run_DESeq2_filter(X,f.a,filt,group)
deg.wilcox     = calc_DE(X.cpm+1,f.a,filt,2-group)

# Filter DE results for common co-expression
edgeR.results = plot_deg_coexp(deg.edgeR.filt, "edgeR", runid, filtMin=5, freq.net)
wilcox.results = plot_deg_coexp(deg.wilcox, "wilcox", runid, filtMin=5, freq.net)
deseq.results = plot_deg_coexp(deg.deseq.filt, "deseq", runid, filtMin=5, freq.net)
