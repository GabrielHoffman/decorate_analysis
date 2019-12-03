#! /usr/bin/env Rscript

library(getopt)

spec = matrix(c(
  'nbatches', 'n', 2, "integer",
  'batch'   , 'b', 0, "integer",
  'prefix'  , 'p', 1, "character"
), byrow=TRUE, ncol=4)
opt = getopt(spec)


suppressPackageStartupMessages(library(decorate))
suppressPackageStartupMessages(library(GenomicRanges))
suppressPackageStartupMessages(library(limma)     )  
suppressPackageStartupMessages(library(data.table)  )  
suppressPackageStartupMessages(library(PRROC)   )    
suppressPackageStartupMessages(library(mvtnorm) )
# library(EnsDb.Hsapiens.v86)    
suppressPackageStartupMessages(library(BiocParallel))
# register(SnowParam(4, progressbar=TRUE))
register(SerialParam(progressbar=TRUE))

set.seed(1)
# ensdb = EnsDb.Hsapiens.v86

# autocorrelation decaying with distance
autocorr.mat <- function(gr, rho = 0.5, s=10) {
    mat <- diag(length(gr))
    pos = start(gr)
    row_mat = row(mat)
    col_mat = col(mat)

    # convert to position
    row_mat[] = pos[row_mat] / s
    col_mat[] = pos[col_mat] / s	

    # correlation
    rho^abs(row_mat-col_mat)
}

# simulate features
simulate_features = function(gr, rho, s, info, beta_disease, beta_confounding, n_clusters, diffCorrScale=1 ){

	n_features = floor(length(gr) / n_clusters)

	# Generate dataset
	epiData = lapply( seq_len(n_clusters), function(i){

		# which position based on cluster
		idx = seq((i-1)*n_features+1, i*n_features, by=1)
		
		# simulate autocorrelation matrix
		C1 = autocorr.mat(gr[idx], rho, s)
		if( diffCorrScale == 1){
			C2 = C1
		}else{
			C2 = autocorr.mat(gr[idx], rho*diffCorrScale, s)
		}

		# simulate data based on correlation
		X = matrix(NA, nrow=nrow(info), ncol= n_features)
		X[info$Disease == 1,] = rmvnorm( sum(info$Disease == 1), sigma=C1)
		X[info$Disease != 1,] = rmvnorm( sum(info$Disease != 1), sigma=C2)
		
		confound = model.matrix(~ Disease+0, info)[,1] * beta_disease + model.matrix(~ Confound+0, info)[,1] * beta_confounding

		# return X with rows as features
		t(X + confound)
	})
	epiData = do.call("rbind", epiData)
	rownames(epiData) = names(gr)[1:nrow(epiData)]
	epiData
}

run_simulation = function( simLocation, sim_params, i, info, n_clusters){
	cat("job: ", i)

	info = data.frame(Disease = as.character(sample(2, sim_params$n_samples[i], replace=TRUE)), 
		Confound = as.character(sample(2, sim_params$n_samples[i], replace=TRUE)))

	# Simulate data
	epiData = simulate_features(
		gr 					= simLocation, 
		rho					= sim_params$rho[i],
		s 					= 1000, 
		info 				= info, 
		beta_disease		= sim_params$beta_disease[i], 
		beta_confounding 	= sim_params$beta_confounding[i],
		n_clusters 			= n_clusters,
		 diffCorrScale		= sim_params$diffCorrScale[i])

	# cor(t(epiData[,info$Disease==1]))[1:4, 1:4]
	# cor(t(epiData[,info$Disease==2]))[1:4, 1:4]

	gr = simLocation[rownames(epiData)]

	if( sim_params$useResid[i] ){
		# Compute residuals  
		design = model.matrix(~ Disease + Confound, info)
		fit = lmFit(epiData, design)
		residValues = residuals( fit, epiData)
	}else{
		residValues = epiData
	}

	# compute mean R^2 value of variables
	# summary(lm(epiData[1,] ~ Disease + Confound,info))$r.squared
	corrValues = sapply(1:nrow(residValues), function(k) 1-cor(residValues[k,], epiData[k,])^2)

	# Compute clustering
	treeList = runOrderedClusteringGenome( residValues, gr, method.corr="spearman", quiet=TRUE )

	# Plot Correlation vs Distance
	# dfDist = evaluateCorrDecay( treeList, gr, verbose=FALSE)
	# plotCorrDecay( dfDist, outlierQuantile=1e-5 )

	# Clustering
	treeListClusters = createClusters( treeList, method = "meanClusterSize", meanClusterSize=length(simLocation) / n_clusters )

	# Evaluate strength of correlation for each cluster
	clstScore = scoreClusters(treeList, treeListClusters )

	# # Filter to retain only strong clusters
	# # If lead eigen value fraction (LEF) > 30% then keep clusters
	# # LEF is the fraction of variance explained by the first eigen-value
	clustInclude = retainClusters( clstScore, "LEF", 0.0 )

	# # get retained clusters
	treeListClusters_filter = filterClusters( treeListClusters, clustInclude)

	# # collapse redundant clusters
	treeListClusters_collapse = collapseClusters( treeListClusters_filter, simLocation, jaccardCutoff=0.9)

	# Plot correlations and clusters in region defind by query
	# get single entry giving range of the region
	query = range(simLocation)
	 # plotDecorate( ensdb, treeList, treeListClusters_collapse, simLocation, query)  

	# Evaluate Differential Correlation between two subsets of data
	# Use Spearman correlation to reduce the effect of outliers
	resDiffCorr = evalDiffCorr( residValues, info$Disease, simLocation[rownames(residValues)], 
	  treeListClusters_collapse, method='Box',  method.corr="spearman")

	# Summarize results
	res = combineResults( resDiffCorr, clstScore, treeListClusters, simLocation)

	# combine results and simulation information
	with(res, data.frame(		i 					= i,
	n_samples 			= nrow(info),
	rho					= sim_params$rho[i],
	beta_disease		= sim_params$beta_disease[i], 
	beta_confounding 	= sim_params$beta_confounding[i],
	diffCorrScale		= sim_params$diffCorrScale[i], 
	mean_rsq 			= mean(corrValues),
	id, chrom, cluster, pValue, stat,N,LEF))
}



# Parameters for simulation
############################

pos = seq(1, 3e7, by=3000)[1:10000]

simLocation = GRanges("chr1", IRanges(pos, pos+1, names=paste0('peak_', 1:length(pos))))

n_clusters = 2000
n_features_per_cluster = length(simLocation) / n_clusters
n_features_per_cluster

sim_params = expand.grid( 	n_samples 	= c(100, 200),
							rho 			= c(.9), 
							beta_disease 	= c(0),
							beta_confounding= c(0, .05, .1, .4, .6), 
							diffCorrScale 	= seq(1, 1.06, length.out=5),
							useResid 		= c(TRUE, FALSE))

sim_params = unique(sim_params)
rownames(sim_params) = 1:nrow(sim_params)

if( opt$batch == 1){
	saveRDS(sim_params, file=paste0(opt$prefix, '_sim_params.RDS') )
}
# run
idx = floor( seq(0, nrow(sim_params), length.out=opt$nbatches+1) )

sim_results = lapply( (idx[opt$batch]+1):idx[opt$batch+1], function(i){

	cat("\r", i, ' / ', nrow(sim_params), '   ')

	run_simulation( simLocation, sim_params, i, info, n_clusters)
})

# combine results
resSim = data.table(do.call("rbind", sim_results))
resSim = resSim[!is.nan(stat)&is.finite(stat),]

saveRDS(resSim, file=paste0(opt$prefix, '_',opt$batch, ".RDS"))

