# Gabriel Hoffman
# December 4, 2019

##################
# Run Simulaions #
##################

cd /sc/orga/work/hoffmg01/decorate/sims/

# /hpc/users/hoffmg01/build2/decorate_analysis/simulations_full.R --nbatches 30 --batch 1 --prefix /sc/orga/work/hoffmg01/decorate/sims/decorate_1 --seed 1

ml openssl udunits proj gdal geos
NBATCHES=60
export OMP_NUM_THREADS=1
rm -f decorate_*
rm -f jobs.sh
for SEED in $(seq 1 50);
do
seq 1 $NBATCHES | parallel -P1 echo "/hpc/users/hoffmg01/build2/decorate_analysis/simulations_full.R --nbatches $NBATCHES --batch {} --prefix /sc/orga/work/hoffmg01/decorate/sims_v2/decorate_${SEED} --seed ${SEED}" >> jobs.sh
done

cat jobs.sh | parallel -P60 "echo {} | sh"




##########################
# Plots from simulations #
##########################
 
# suppressPackageStartupMessages(library(PRROC))   
suppressPackageStartupMessages(library(pROC))   
suppressPackageStartupMessages(library(ggplot2))  
suppressPackageStartupMessages(library(foreach))  
suppressPackageStartupMessages(library(doParallel))    
suppressPackageStartupMessages(library(data.table))    

# path = '/Users/gabrielhoffman/link/decorate/'
path = './'
files = dir(path, pattern="RDS$", full.names=TRUE)
files = files[grep('sim_params', files, invert=TRUE)]

resSim = mclapply( files, function(file){

	cat("\r", grep(file, files), ' / ', length(files), '      ')

	res = readRDS(file)
	res$seed = as.numeric(sapply(strsplit(basename(file), '_'), function(x) x[2]))

	# sort by LEF
	res = res[,.SD[order(LEF, decreasing=TRUE),], by=i]

	# keep top 1000 for each i
	res[,.SD[order(LEF, decreasing=TRUE)<=1000,], by=i]
	}, mc.cores=40)
resSim = data.table(do.call("rbind", resSim))
setkey(resSim,i)
table(resSim$id)
# resSim = resSim[resSim$N==15,]

sim_params = readRDS(paste0(path, 'decorate_1_sim_params.RDS'))

prRes = sim_params
prRes$i = rep(NA, nrow(prRes))
prRes$AUC = rep(NA, nrow(prRes))
prRes$AUC.low = rep(NA, nrow(prRes))
prRes$AUC.high = rep(NA, nrow(prRes))

count = 0
for( beta_disease in unique(sim_params$beta_disease) ){
for( beta_confounding in unique(sim_params$beta_confounding) ){
for( useResid in unique(sim_params$useResid) ){
for( rho in unique(sim_params$rho) ){
for( n_samples in unique(sim_params$n_samples) ){
for( n_features_per_cluster in unique(sim_params$n_features_per_cluster) ){

	idxs = sim_params[(sim_params$beta_disease==beta_disease) & 
			(sim_params$beta_confounding==beta_confounding) &
			(sim_params$useResid==useResid) &
			(sim_params$rho==rho) &
			(sim_params$n_samples==n_samples) &
			(sim_params$n_features_per_cluster==n_features_per_cluster),]

	idxs = as.numeric(rownames(idxs))

	idx_baseline = idxs[1]
	res1 = resSim[i==idx_baseline,]
	res2_pre = resSim[i%in%idxs[-1],]

	prList = lapply( idxs[-1], function(idx){
	# prList = foreach( idx = idxs[-1]) %dopar% {
		res2 = res2_pre[i==idx,]

		if( nrow(res2) == 0){
			res = NULL
		}else{
			# pr.curve(abs(res2$stat), abs(res1$stat), curve=TRUE, rand.compute=TRUE)
			pred = abs(c(res2$stat, res1$stat))
			response = c(rep(0, nrow(res1)), rep(1, nrow(res2)))
			curve = roc( response, pred, quiet=TRUE)
			res = t(data.frame(ci(curve)))
			colnames(res) = c("AUC.low", "AUC", "AUC.high")
		}
		res
	})
	res = data.frame(do.call("rbind", prList))

	count = count + 1
	cat("\r", count, '      ')
	# plot AUPR
	# prRes$AUPR[idxs[-1]] = sapply( prList, function(pr) pr$auc.integral )

	prRes$i[idxs] = idxs
	prRes$AUC[idxs[-1]] = res$AUC
	prRes$AUC.low[idxs[-1]] = res$AUC.low
	prRes$AUC.high[idxs[-1]] = res$AUC.high
}
}
}
}
}
}

# False positive rate
resFPR = resSim[,data.frame(FPR=sum(pValue<0.05)/length(pValue)),by=i]
resFPR = merge(resFPR, prRes, by='i', all.x=TRUE)


# save results
# save(list=c('prRes','sim_params', 'resFPR'), file="session.RDATA")

# load results
# load(file="decorate/session.RDATA")

suppressPackageStartupMessages(library(ggplot2))  

# False positive rate
######################
# ggplot(resFPR[resFPR$diffCorrScale==1,], aes(as.factor(n_samples), FPR, fill=as.factor(n_features_per_cluster))) + geom_bar(stat="identity", position='dodge') + facet_wrap(~beta_confounding + useResid) + theme_bw(8) + theme(aspect.ratio=1) + geom_hline(yintercept=0.05, color="grey60", linetype="dashed")

ncol = length(unique(sim_params$beta_confounding))

df = resFPR[resFPR$diffCorrScale==1,]

fig = ggplot(df, aes(as.factor(n_samples), FPR, fill=as.factor(useResid))) + geom_bar(stat="identity", position='dodge') + facet_wrap(~n_features_per_cluster+beta_confounding, ncol=ncol) + theme_bw(8) + theme(aspect.ratio=1, legend.position="bottom", strip.text.x = element_text(margin = margin(.05,0,.05,0, "cm"))) + geom_hline(yintercept=0.05, color="grey60", linetype="dashed") + scale_fill_manual("Use residuals", values=c("black", "dodgerblue")) + ylab("False Positive Rate") + xlab("Number of samples")
ggsave(fig, file="false_positives.pdf")



# AUC versus diffCorScale
#########################
ncol = length(unique(sim_params$beta_confounding))

df = prRes[!is.na(prRes$AUC),]
df = df[df$n_features_per_cluster == 20,]
# df = df[df$n_samples == 100,]
# df = df[df$beta_confounding == 0.9,]

fig = ggplot(df, aes(diffCorrScale, AUC, color=useResid)) + geom_point(size=.1) + geom_errorbar(aes(ymin=AUC.low, ymax=AUC.high, color=useResid), width=.0001) + geom_line() + theme_bw(8) + theme(aspect.ratio=1, legend.position="bottom", axis.text.x = element_text(size=5, angle=45, hjust=1), strip.text.x = element_text(margin = margin(.05,0,.05,0, "cm"))) + facet_wrap( ~ n_samples + beta_confounding, ncol=ncol ) + xlab("Effect size") + geom_hline(yintercept=0.5, linetype="dashed", color='grey60') + ylim(0, 1)  + scale_color_manual("Use residuals", values=c("black", "dodgerblue")) + ylim(0.5, 1) #+ ggtitle("AUC vs diffCorrScale")
ggsave(fig, file="fig_diffCorrScale.pdf")



# AUC vs confounding
####################
ncol = length(unique(sim_params$diffCorrScale)) - 1

# df = prRes[!is.na(prRes$AUC)&(prRes$n_samples==n_samples),]
df = prRes[!is.na(prRes$AUC),]
df = df[df$n_features_per_cluster == 10,]
df$diffCorrScale = round(df$diffCorrScale, 3)

fig = ggplot(df, aes(beta_confounding, AUC, color=useResid)) + geom_point(size=.1) + geom_errorbar(aes(ymin=AUC.low, ymax=AUC.high, color=useResid), width=.0001) + geom_line() + theme_bw(8) + theme(aspect.ratio=1, legend.position="bottom", strip.text.x = element_text(margin = margin(.05,0,.05,0, "cm"))) + xlab("Effect size of confounding factor") + geom_hline(yintercept=0.5, linetype="dashed", color='grey60') + ylim(0, 1) + facet_wrap( ~ n_samples + diffCorrScale, ncol=ncol, labeller=) + scale_color_manual("Use residuals", values=c("black", "dodgerblue")) + ylim(0.5, 1)# + ggtitle("AUC vs effect of confounder")
ggsave(fig, file="fig_confounding.pdf")




# AUC vs n_features_per_cluster
################################
ncol = length(unique(sim_params$beta_confounding)) 
df = prRes[!is.na(prRes$AUC),]
table(df$beta_confounding)
df = df[(df$useResid == TRUE),]
df$diffCorrScale = round(df$diffCorrScale, 3)

fig = ggplot(df, aes(n_features_per_cluster, AUC, color=as.factor(n_samples))) + geom_point(size=.1) + geom_errorbar(aes(ymin=AUC.low, ymax=AUC.high, color=as.factor(n_samples)), width=.0001) + geom_line() + theme_bw(8) + theme(aspect.ratio=1, legend.position="bottom", strip.text.x = element_text(margin = margin(.05,0,.05,0, "cm"))) + xlab("Number of features per cluster") + geom_hline(yintercept=0.5, linetype="dashed") + ylim(0, 1) + facet_wrap(~diffCorrScale+beta_confounding, ncol=ncol) + scale_color_discrete("# Samples") + ylim(0.5, 1) # + ggtitle("AUC vs number of features")
ggsave(fig, file="fig_features.pdf")



# # AUC vs confonding variable colored by n_features_per_cluster
# ncol = length(unique(sim_params$diffCorrScale)) - 1

# df = prRes[!is.na(prRes$AUC),]
# df = df[(df$useResid == FALSE),]

# ggplot(df, aes(beta_confounding, AUC, color=as.character(n_features_per_cluster))) + geom_point(size=.1) + geom_errorbar(aes(ymin=AUC.low, ymax=AUC.high, color=useResid), width=.0001) + geom_line() + theme_bw(12) + theme(aspect.ratio=1) + xlab("Effect size of confounding factor") + geom_hline(yintercept=0.5, linetype="dashed") + ylim(0, 1) + facet_wrap( ~ n_samples  + diffCorrScale, ncol=ncol) + ggtitle("AUC vs diffCorrScale")





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


