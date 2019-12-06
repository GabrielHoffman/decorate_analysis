
# q()
# R

cd /sc/orga/work/hoffmg01/decorate/sims/

# /hpc/users/hoffmg01/build2/decorate_analysis/simulations_full.R --nbatches 30 --batch 1 --prefix /hpc/users/hoffmg01/work/decorate/decorate

ml openssl udunits proj gdal geos
NBATCHES=60
export OMP_NUM_THREADS=1
rm -f decorate_*
seq 1 $NBATCHES | parallel -P60 /hpc/users/hoffmg01/build2/decorate_analysis/simulations_full.R --nbatches $NBATCHES --batch {} --prefix /sc/orga/work/hoffmg01/decorate/sims/decorate


# Evaluate PR curve
###################
 
suppressPackageStartupMessages(library(PRROC))   
suppressPackageStartupMessages(library(ggplot2))   
suppressPackageStartupMessages(library(data.table))    

path = '/Users/gabrielhoffman/link/decorate/'
# path = '/hpc/users/hoffmg01/work/decorate/'
files = dir(path, full.names=TRUE)
files = files[grep('sim_params', files, invert=TRUE)]

resSim = lapply( files, function(file){
	readRDS(file)
	})
resSim = data.table(do.call("rbind", resSim))
table(resSim$N)
# resSim = resSim[resSim$N==15,]

sim_params = readRDS(paste0(path, 'decorate_sim_params.RDS'))

prRes = sim_params
prRes$AUPR = rep(NA, nrow(prRes))

for( beta_disease in unique(sim_params$beta_disease) ){
for( beta_confounding in unique(sim_params$beta_confounding) ){
for( useResid in unique(sim_params$useResid) ){
for( rho in unique(sim_params$rho) ){
for( n_samples in unique(sim_params$n_samples) ){

	idxs = sim_params[(sim_params$beta_disease==beta_disease) & 
			(sim_params$beta_confounding==beta_confounding) &
			(sim_params$useResid==useResid) &
			(sim_params$rho==rho) &
			(sim_params$n_samples==n_samples),]
	idxs = as.numeric(rownames(idxs))

	idx_baseline = idxs[1]

	prList = lapply( idxs[-1], function(idx){

		res1 = resSim[i==idx_baseline,]
		res2 = resSim[i==idx,]

		pr.curve(abs(res2$stat), abs(res1$stat), curve=TRUE, rand.compute=TRUE)
		})

	# plot( prList[[2]], rand.plot=TRUE)

	# plot AUPR
	prRes$AUPR[idxs[-1]] = sapply( prList, function(pr) pr$auc.integral )
}
}
}
}
}

# ggplot(prRes[!is.na(prRes$AUPR),], aes(diffCorrScale, AUPR, fill=diffCorrScale)) + geom_bar(stat="identity") + theme_bw(17) + theme(aspect.ratio=1) + facet_wrap( ~ beta_disease + useResid)


ncol = length(unique(sim_params$beta_disease))

df = prRes[!is.na(prRes$AUPR),]

ggplot(df, aes(diffCorrScale, AUPR, color=useResid)) + geom_point() + geom_line() + theme_bw(12) + theme(aspect.ratio=1) + facet_wrap( ~ n_samples + beta_disease, ncol=ncol ) + xlab("Effect size") + geom_hline(yintercept=0.5, linetype="dashed") + ylim(0, 1) + ggtitle("AUPR vs diffCorrScale")



ncol = length(unique(sim_params$diffCorrScale)) - 1

# df = prRes[!is.na(prRes$AUPR)&(prRes$n_samples==n_samples),]
df = prRes[!is.na(prRes$AUPR),]

ggplot(df, aes(beta_disease, AUPR, color=useResid)) + geom_point() + geom_line() + theme_bw(12) + theme(aspect.ratio=1) + xlab("Effect size of confounding factor") + geom_hline(yintercept=0.5, linetype="dashed") + ylim(0, 1) + facet_wrap( ~ n_samples + diffCorrScale, ncol=ncol) + ggtitle("AUPR vs effect of confounder")




# simulate blocks of all diff corr
# 1) test for Disease
# 2) test for Disese + confounder

# Condition
# 1) diffCorrScale = 1
# 2) diffCorrScale= different

# also vary
# a) number of features
# b) rho


# N = 100
# x = rnorm(N, 100,12)
# beta = c(2,7)
# y1 = x * beta[1] + rnorm(N)
# y2 = x * beta[2] + rnorm(N)

# cov(y1, y2)

# cov(x * beta[1], x * beta[2] )

# X = cbind(1, x)
# Beta = t(cbind(1, c(1,2,3,4)))

# Sigma = autocorr.mat( simLocation[1:4], rho=.99) * 100

# Y = X %*% Beta + rmvnorm(100, rep(0,4), Sigma)

# cov(Y)[1,2]

# cov(X %*% Beta[,1], X %*% Beta[,2]) + Sigma[1,2]


