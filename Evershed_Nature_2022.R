#
#  ********************************
#  **  Written by Yoan Diekmann  **
#  **    y.diekmann@ucl.ac.uk    **
#  ********************************
#
#
#
library('parallel')
numCores <- detectCores()
library('maptools')
library('splancs')
library('stringr')
library('DEoptimR')
library('optimization')
library('learnPopGen')
library('latex2exp')
#
#
#
# http://stackoverflow.com/questions/2769510/numeric-comparison-difficulty-in-r
.TOL <- 1e-10
'%==%' <- function(x, y) return (abs(x-y) <= .TOL)
#
.add.alpha <- function(alpha, col=1) {
    if(missing(alpha)) stop("Please provide a vector of alpha values.")
    col_t = col2rgb(col)/255
    apply( as.matrix(alpha), 2, function(x) rgb(col_t[1], col_t[2], col_t[3], alpha=x) )
}
#
#
#
.convert_genotypes_to_binary <- function(time_points, genotypes) {
  lp = data.frame(allele=NA, mid=NA)
  for (i in 1:length(time_points)) {
    lp[(2*i - 1):(2*i), ]$mid = time_points[i]
    #
    lp[2*i - 1, ]$allele = as.numeric(substr(toString(genotypes[i]), 1, 1)=='A')
    lp[2*i, ]$allele = as.numeric(substr(toString(genotypes[i]), 2, 2)=='A')
  }
  return(lp)
}
#
fit_glm_sigmoidal <- function(time_points, genotypes, date_range=(12000:0)) {
    lp = .convert_genotypes_to_binary(time_points, genotypes)
    model <- glm(allele ~ mid, family='binomial', data=lp)
    pred <- predict(model, newdata=data.frame(mid=date_range), se.fit=T)
    fnc <- family(model)$linkinv
    pred$pred <- fnc(pred$fit)
    pred$LC <- fnc(pred$fit - 1.96*pred$se.fit)
    pred$UC <- fnc(pred$fit + 1.96*pred$se.fit)
    #    
    return(list(lp=lp, pred=pred))
}
#
#
#
plot_driver <- function(driver_name, driver, Britain_MAP=5, Britain_u95=4, Britain_l95=3,
                        NCentEur_MAP=8, NCentEur_u95=7, NCentEur_l95=6,
                        South_MAP=11, South_u95=10, South_l95=9,
                        Baltic_MAP=12, Baltic_u95=NA, Baltic_l95=NA) {
  #
  figure_name = paste('results/', driver_name, '_data.pdf', sep='')
  if (file.exists(figure_name)) { return() }
  #
  #
  subplot <- function(MAP, u95, l95, header) {
    col_indices = c(MAP, u95, l95)
    col_indices = col_indices[!is.na(col_indices)]
    plot(NA, xlim=c(8000, 2500), ylim=c(min(driver[, col_indices], na.rm = T),
                                        max(driver[, col_indices], na.rm = T)),
         xlab='', ylab='', main=header)
    #
    if (!(is.na(u95) | is.na(l95))) {
      polygon(c(driver$startBP[1],
                as.vector(rbind(driver$startBP, driver$endBP)), 
                driver$endBP[length(driver$endBP)],
                rev(as.vector(rbind(driver$startBP, driver$endBP)))),
              c(driver[1, l95],
                as.vector(rbind(driver[, l95], driver[, l95])),
                driver[length(driver[, u95]), u95],
                rev(as.vector(rbind(driver[, u95], driver[, u95])))),
              col=.add.alpha(0.2), border = NA)
    }
    #
    lines(as.vector(rbind(driver$startBP, driver$endBP)),
          as.vector(rbind(driver[, MAP], driver[, MAP])), lwd=2)
  }
  #
  #
  pdf(file = figure_name, height = 10)
  #
  par(mfrow=c(4,1))
  #
  #
  par(oma=c(3, 4, 0, 0))
  par(mar=c(3, 1, 4, 2))
  #
  subplot(Britain_MAP, Britain_u95, Britain_l95, 'Britain')
  subplot(NCentEur_MAP, NCentEur_u95, NCentEur_l95, 'North and Central Europe')
  subplot(South_MAP, South_u95, South_l95, 'South')
  subplot(Baltic_MAP, Baltic_u95, Baltic_l95, 'Baltic')
  #
  mtext(text = 'year [BP]', side = 1, line = 1, outer = T)
  mtext(text = driver_name, side = 2, line = 2, outer = T)
  #
  dev.off()
}
#
#
#
normalise_table <- function(table, first_col=3, range=c(0, 1)) {
    .max = max(table[, first_col:ncol(table)], na.rm = TRUE)
    .min = min(table[, first_col:ncol(table)], na.rm = TRUE)
    #
    table[, first_col:ncol(table)] = 
                        (((table[, first_col:ncol(table)] - .min) / (.max - .min)) * (range[2] - range[1])) + range[1]
  return(table)
}
#
#
#
getKMLnames <- function(file) {
  #
  names <- readLines(file)
  names <- names[grep('\t\t\t<name>', names)]
  names <- gsub('\t\t\t<name>', '', names)
  names <- gsub('</name>', '', names)
  #
  return(names)
}
#
#
#
# **************************
# ***       PART  I      ***
# **************************
#
dat = read.csv(file = 'data/AADRv44.subset.mpileups.tsv', sep = '\t', header=TRUE)
dat = dat[dat$genotype!='', ]
time_points = dat$mean_date
genotypes = dat$genotype
#
if (!file.exists('results/aDNA_data_8000_2500_BP_merged.pdf')) {
  date_range=(12000:0)
  res = fit_glm_sigmoidal(time_points, genotypes, date_range)
  pdf(file = 'results/aDNA_data_8000_2500_BP_merged.pdf', height = 5)
  plot(res$lp$mid, res$lp$allele, xlim=c(12000, 0), pch='+',
       xlab='year [BP]', ylab='LP allele frequency', main='aDNA data 8000-2500 BP')
  lines(date_range, res$pred$pred, col='red')
  lines(date_range, res$pred$LC, col='orange', lty='dashed')
  lines(date_range, res$pred$UC, col='orange', lty='dashed')
  dev.off()
}
#
polygons = getKMLcoordinates('data/polygons.model.6.kml', ignoreAltitude=T)
names(polygons) = getKMLnames('data/polygons.model.6.kml')
#
sel_Britain = as.numeric(rownames(pip(pts=data.frame(x=dat$long, y=dat$lat),
                                      poly=polygons$`British Isles`, bound=TRUE)))
time_points_Britain = time_points[sel_Britain]
genotypes_Britain = genotypes[sel_Britain]
#
sel_NCentEur = as.numeric(rownames(pip(pts=data.frame(x=dat$long, y=dat$lat),
                                       poly=polygons$`Rhine Danube axis`, bound=TRUE))) 
time_points_NCentEur = time_points[sel_NCentEur]
genotypes_NCentEur = genotypes[sel_NCentEur]
#
sel_South = as.numeric(rownames(pip(pts=data.frame(x=dat$long, y=dat$lat),
                                      poly=polygons$`Mediterranean Europe`, bound=TRUE))) 
time_points_South = time_points[sel_South]
genotypes_South = genotypes[sel_South]
#
sel_Baltic = as.numeric(rownames(pip(pts=data.frame(x=dat$long, y=dat$lat),
                                    poly=polygons$`Baltic region`, bound=TRUE))) 
time_points_Baltic = time_points[sel_Baltic]
genotypes_Baltic = genotypes[sel_Baltic]
#
if (!file.exists('results/aDNA_data_8000_2500_BP.pdf')) {
  date_range=(12000:0)
  res_Britain = fit_glm_sigmoidal(time_points_Britain, genotypes_Britain, date_range)
  res_NCentEur = fit_glm_sigmoidal(time_points_NCentEur, genotypes_NCentEur, date_range)
  res_South = fit_glm_sigmoidal(time_points_South, genotypes_South, date_range)
  res_Baltic = fit_glm_sigmoidal(time_points_Baltic, genotypes_Baltic, date_range)
  #
  pdf(file = 'results/aDNA_data_8000_2500_BP.pdf', height = 7, width = 9)
  #
  par(mfrow=c(2,2))
  #
  par(oma=c(3, 4, 0, 0))
  par(mar=c(3, 1, 4, 2))
  #
  plot(res_Britain$lp$mid, res_Britain$lp$allele, xlim=c(12000, 0), pch='+',
       xlab='', ylab='', main='Britain')
  lines(date_range, res_Britain$pred$pred, col='red')
  lines(date_range, res_Britain$pred$LC, col='orange', lty='dashed')
  lines(date_range, res_Britain$pred$UC, col='orange', lty='dashed')
  #
  plot(res_NCentEur$lp$mid, res_NCentEur$lp$allele, xlim=c(12000, 0), pch='+',
       xlab='', ylab='', main='North and Central Europe')
  lines(date_range, res_NCentEur$pred$pred, col='red')
  lines(date_range, res_NCentEur$pred$LC, col='orange', lty='dashed')
  lines(date_range, res_NCentEur$pred$UC, col='orange', lty='dashed')
  #
  plot(res_South$lp$mid, res_South$lp$allele, xlim=c(12000, 0), ylim = c(0, 1), pch='+',
       xlab='', ylab='', main='South')
  lines(date_range, res_South$pred$pred, col='red')
  lines(date_range, res_South$pred$LC, col='orange', lty='dashed')
  lines(date_range, res_South$pred$UC, col='orange', lty='dashed')
  #
  plot(res_Baltic$lp$mid, res_Baltic$lp$allele, xlim=c(12000, 0), ylim = c(0, 1), pch='+',
       xlab='', ylab='', main='Baltic')
  lines(date_range, res_Baltic$pred$pred, col='red')
  lines(date_range, res_Baltic$pred$LC, col='orange', lty='dashed')
  lines(date_range, res_Baltic$pred$UC, col='orange', lty='dashed')
  
  #
  mtext(text = 'year [BP]', side = 1, line = 1, outer = T)
  mtext(text = 'LP allele frequency', side = 2, line = 2, outer = T)
  #
  dev.off()
}
#
#
#
# NORMALISATION
#
cluster_normalised = read.csv(file = 'data/cluster stat.csv')
plot_driver('cluster_normalised', cluster_normalised, Britain_MAP=3, Britain_u95=NA, Britain_l95=NA,
                                                      NCentEur_MAP=5, NCentEur_u95=NA, NCentEur_l95=NA,
                                                      South_MAP=6, South_u95=NA, South_l95=NA,
                                                      Baltic_MAP=4, Baltic_u95=NA, Baltic_l95=NA)
cluster_normalised_inv = cbind(cluster_normalised[, 1:2], 1 - cluster_normalised[, 3:ncol(cluster_normalised)])
plot_driver('cluster_normalised_inv', cluster_normalised_inv, Britain_MAP=3, Britain_u95=NA, Britain_l95=NA,
                                                              NCentEur_MAP=5, NCentEur_u95=NA, NCentEur_l95=NA,
                                                              South_MAP=6, South_u95=NA, South_l95=NA,
                                                              Baltic_MAP=4, Baltic_u95=NA, Baltic_l95=NA)
#
domestics_normalised = read.csv(file = 'data/domestic animal proportion.csv')
plot_driver('domestics_normalised', domestics_normalised, Britain_MAP=3, Britain_u95=NA, Britain_l95=NA,
                                                          NCentEur_MAP=5, NCentEur_u95=NA, NCentEur_l95=NA,
                                                          South_MAP=6, South_u95=NA, South_l95=NA,
                                                          Baltic_MAP=4, Baltic_u95=NA, Baltic_l95=NA)
domestics_normalised_inv = cbind(domestics_normalised[, 1:2], 1 - domestics_normalised[, 3:ncol(domestics_normalised)])
plot_driver('domestics_normalised_inv', domestics_normalised_inv, Britain_MAP=3, Britain_u95=NA, Britain_l95=NA,
                                                                  NCentEur_MAP=5, NCentEur_u95=NA, NCentEur_l95=NA,
                                                                  South_MAP=6, South_u95=NA, South_l95=NA,
                                                                  Baltic_MAP=4, Baltic_u95=NA, Baltic_l95=NA)
#
milk_normalised = read.csv(file = 'data/milk proportion.csv')
plot_driver('milk_normalised', milk_normalised, Britain_MAP=3, Britain_u95=NA, Britain_l95=NA,
                                                NCentEur_MAP=5, NCentEur_u95=NA, NCentEur_l95=NA,
                                                South_MAP=6, South_u95=NA, South_l95=NA,
                                                Baltic_MAP=4, Baltic_u95=NA, Baltic_l95=NA)
#
pop_fluc_normalised = read.csv(file = 'data/pop fluctuations stat.csv')
plot_driver('pop_fluc_normalised', pop_fluc_normalised, Britain_MAP=3, Britain_u95=NA, Britain_l95=NA,
                                                        NCentEur_MAP=5, NCentEur_u95=NA, NCentEur_l95=NA,
                                                        South_MAP=6, South_u95=NA, South_l95=NA,
                                                        Baltic_MAP=4, Baltic_u95=NA, Baltic_l95=NA)
pop_fluc_normalised_inv = cbind(pop_fluc_normalised[, 1:2], 1 - pop_fluc_normalised[, 3:ncol(pop_fluc_normalised)])
plot_driver('pop_fluc_normalised_inv', pop_fluc_normalised_inv, Britain_MAP=3, Britain_u95=NA, Britain_l95=NA,
                                                                NCentEur_MAP=5, NCentEur_u95=NA, NCentEur_l95=NA,
                                                                South_MAP=6, South_u95=NA, South_l95=NA,
                                                                Baltic_MAP=4, Baltic_u95=NA, Baltic_l95=NA)
#
insolation = read.csv(file = 'data/midday insolation.csv')
insolation_normalised = data.frame(
                    'startBP' = c(8000, 5250),
                    'endBP' = c(5250, 2500),
                    'Britain' = c(insolation$British.Isles.MAP[1], insolation$British.Isles.MAP[1]),
                    'Baltic' = c(insolation$Baltic.region.MAP[1], insolation$Baltic.region.MAP[1]),
                    'NCentEur' = c(insolation$Rhine.Danube.axis.MAP[1], insolation$Rhine.Danube.axis.MAP[1]),
                    'South' = c(insolation$Mediterranean.Europe.MAP[1], insolation$Mediterranean.Europe.MAP[1]))
plot_driver('insolation_normalised', insolation_normalised, Britain_MAP=3, Britain_u95=NA, Britain_l95=NA,
                                                            NCentEur_MAP=5, NCentEur_u95=NA, NCentEur_l95=NA,
                                                            South_MAP=6, South_u95=NA, South_l95=NA,
                                                            Baltic_MAP=4, Baltic_u95=NA, Baltic_l95=NA)
insolation_normalised_inv = cbind(insolation_normalised[, 1:2],
                                  1 - insolation_normalised[, 3:ncol(insolation_normalised)])
plot_driver('insolation_normalised_inv', insolation_normalised_inv, Britain_MAP=3, Britain_u95=NA, Britain_l95=NA,
            NCentEur_MAP=5, NCentEur_u95=NA, NCentEur_l95=NA,
            South_MAP=6, South_u95=NA, South_l95=NA,
            Baltic_MAP=4, Baltic_u95=NA, Baltic_l95=NA)
#
input = list('cluster' = cluster_normalised, 'cluster_inv' = cluster_normalised_inv,
             'domestics' = domestics_normalised, 'domestics_inv' = domestics_normalised_inv,
             'milk' = milk_normalised,
             'pop_fluc_globally' = pop_fluc_normalised, 'pop_fluc_inv_globally' = pop_fluc_normalised_inv,
             'insolation_globally' = insolation_normalised, 'insolation_inv_globally' = insolation_normalised_inv)
#
#
#
# ****************************
# ***  FUNCTIONS, PART II  ***
# ****************************
#
.get_m <- function(ms, i) {
  stopifnot(ncol(ms)==3)
  # https://stackoverflow.com/questions/29388334/find-position-of-first-value-greater-than-x-in-a-vector/29388546
  ms[Position(function(x) x < i, ms$endBP), ]
}
#
#
#
logistic_fct_sel <- function(s, t, x, y0=0.00005, gen=28) {
  y0 / ( y0 + (1 - y0)*exp(-s*((t-x)/gen)) )
}
#
#
#
f <- function(p_t, s, t, x, y0, gen=28) {
  # So the time taken to change from p0 to pt in the dominant case is exactly
  # the time taken to change from (1 − pt) to (1 − p0) in the recessive case.
  (((1/s) * ( -(1/(1-y0)) + log((1-y0)/y0) + (1/(1-p_t)) - log((1-p_t)/p_t) - s*log(y0/p_t) )) - ((t-x)/gen))^2
}
#
#
#
dominant_sel <- Vectorize(function(s, t, x, y0=0.00005, gen=28) {
  if (s %==% 0) {
    return(y0)
  }
  # https://cooplab.github.io/popgen-notes/#one-locus-models-of-selection # FACTOR OF TWO, 1+2s -> 1+s
  optimize(f, c(0.000000001, 0.999999999), s=2*s, t=t, x=x, y0=y0, gen=gen, tol = .Machine$double.eps)$minimum
}, vectorize.args = c('x'))
#
#
#
ms_to_ss_breaks <- function(ms, s_1, a, start=8000, end=2500) {
  index_first = as.numeric(rownames(.get_m(ms, start)))
  index_last = as.numeric(rownames(.get_m(ms, end+1)))
  #
  breaks = c(start, ms$endBP[(index_first):(index_last-1)], end)
  ss = s_1 * ms[index_first:index_last, 3]^((1/a)-1)
  stopifnot(length(ss) + 1 == length(breaks))
  #
  return(list('breaks' = breaks, 'ss' = ss))
}
#
#
#
# for power analysis with drift
drift_additive_selection <- function(s, total_years, N_e, y0, gen=28) {
  stopifnot(y0 >= 1/(2*N_e) | y0 %==% (1/(2*N_e))) # y0 is bound by 1/(2*N_e)
  gens = ceiling(total_years / gen) - 1
  #
  p = list(c(0))
  while (! p[[1]][length(p[[1]])]>0) {
    pdf(file = NULL)
    p <- drift.selection(p0=y0, Ne=N_e, w=c(1, 1-s, 1-2*s), ngen=gens, nrep=1) # additive model
    dev.off()
  }
  #
  return(c(sapply(p[[1]][1:(length(p[[1]])-1)], function(x) rep(x, each=gen)),
           if (total_years%%gen==0) rep(p[[1]][length(p)], gen) else rep(p[[1]][length(p[[1]])], total_years%%gen)))
}
#
ss_breaks_to_der_al_freq <- function(ss, breaks, fct_selection, y0=0.00005, gen=28, N_e=0) {
  # https://stackoverflow.com/questions/30850400/using-identical-in-r-with-multiple-vectors
  stopifnot(length(ss) + 1 == length(breaks))
  #
  start_freq = y0
  res = numeric()
  for (i in 1:length(ss)) {
    if (N_e > 0) {
      freqs = drift_additive_selection(ss[i], breaks[i] - breaks[i+1] + 1, N_e, y0=start_freq, gen=gen)
    } else {
      freqs = fct_selection(ss[i], breaks[i], (breaks[i]-1):breaks[i+1], y0=start_freq, gen=gen)
    }
    res = c(res, start_freq, freqs[1:(length(freqs)-1)])
    start_freq = freqs[length(freqs)]
  }
  return(c(res, start_freq))
}
#
#
#
data_log_likelihood <- function(genotypes, time_points, derived_allele_freqs, start=8000, end=2500) {
  #
  indices = abs(time_points - start + 1)
  #
  lls = unlist(lapply(1:length(indices),
                      function(index) {
                        str_count(genotypes[index], 'A') * log(derived_allele_freqs[indices[index]]) +
                        str_count(genotypes[index], 'G') * log(1 - derived_allele_freqs[indices[index]]) } ))
  # 0 * -Inf = NaN (when ancestral allele freq. 1 cause of high s and homozygous derived )
  # therefore na.rm below! effectively treats the case as 0 (2*log(1))
  return(sum(lls, na.rm = TRUE))
}
#
#
#
full_likelihood <- function(initial_freq, s_1, a, fct_selection,
                            ms_Britain, genotypes_Britain, time_points_Britain,
                            ms_NCentEur, genotypes_NCentEur, time_points_NCentEur,
                            ms_South, genotypes_South, time_points_South,
                            ms_Baltic, genotypes_Baltic, time_points_Baltic) {
  #
  ss_breaks_Britain = ms_to_ss_breaks(ms_Britain, s_1, a)
  ss_breaks_NCentEur = ms_to_ss_breaks(ms_NCentEur, s_1, a)
  ss_breaks_South = ms_to_ss_breaks(ms_South, s_1, a)
  ss_breaks_Baltic = ms_to_ss_breaks(ms_Baltic, s_1, a)
  #
  return(data_log_likelihood(genotypes_Britain, time_points_Britain,
                            ss_breaks_to_der_al_freq(ss_breaks_Britain$ss, ss_breaks_Britain$breaks, fct_selection,
                                                     y0=initial_freq)) +
         data_log_likelihood(genotypes_NCentEur, time_points_NCentEur,
                            ss_breaks_to_der_al_freq(ss_breaks_NCentEur$ss, ss_breaks_NCentEur$breaks, fct_selection,
                                                     y0=initial_freq)) + 
         data_log_likelihood(genotypes_South, time_points_South,
                            ss_breaks_to_der_al_freq(ss_breaks_South$ss, ss_breaks_South$breaks, fct_selection,
                                                     y0=initial_freq)) +
         data_log_likelihood(genotypes_Baltic, time_points_Baltic,
                             ss_breaks_to_der_al_freq(ss_breaks_Baltic$ss, ss_breaks_Baltic$breaks, fct_selection,
                                                      y0=initial_freq))) 
}
#
#
#
optimise_JDEoptim <- function(stat_normalised, fct_selection, indices_regional_MAPs,
                              genotypes_Britain_=genotypes_Britain, time_points_Britain_=time_points_Britain,
                              genotypes_NCentEur_=genotypes_NCentEur, time_points_NCentEur_=time_points_NCentEur,
                              genotypes_South_=genotypes_South, time_points_South_=time_points_South,
                              genotypes_Baltic_=genotypes_Baltic, time_points_Baltic_=time_points_Baltic,
                              init_f_0=0.0001, init_s=0.03) {
  JDEoptim(lower=c(0, 0, 0.0000001), upper=c(0.05, 0.5, 1), add_to_init_pop=c(init_f_0, init_s, 1),
           fn=function(params)
              -full_likelihood(params[1], params[2], params[3], fct_selection,
                        stat_normalised[, c(1, 2, indices_regional_MAPs[1])], genotypes_Britain_, time_points_Britain_,
                        stat_normalised[, c(1, 2, indices_regional_MAPs[2])], genotypes_NCentEur_, time_points_NCentEur_,
                        stat_normalised[, c(1, 2, indices_regional_MAPs[3])], genotypes_South_, time_points_South_,
                        stat_normalised[, c(1, 2, indices_regional_MAPs[4])], genotypes_Baltic_, time_points_Baltic_),
           trace=TRUE, NP=200, fnscale=100, tol=5e-7, maxiter=400)
}
#
optimise_all <- function(robj_out_path, res_const, optimise_function, selection_function, init=TRUE) {
  #
  if (init==TRUE) {
    results <- mclapply(input, optimise_function, fct_selection=selection_function, indices_regional_MAPs=c(3, 5, 6, 4),
                        init_f_0=res_const$par[1], init_s=res_const$par[2],
                        mc.cores = numCores)
  } else {
    results <- mclapply(input, optimise_function, fct_selection=selection_function, indices_regional_MAPs=c(3, 5, 6, 4),
                        mc.cores = numCores)
  }
  #
  results[['const']] = res_const
  #
  save(results, file = robj_out_path)
  return(results)
}
#
#
#
# ****************************
# ***        PART II       ***
# ****************************
#
const = data.frame('startBP' = c(8000, 5250), 'endBP' = c(5250, 2500), 'tmp' = c(1, 1))
#
#
#
if (!file.exists('results/optim_JDE_additive_1.Robj')) {
  res_const = JDEoptim(lower=c(0, 0), upper=c(0.05, 0.2),
                       fn=function(params)
                         -full_likelihood(params[1], params[2], 1, logistic_fct_sel,
                                          const, genotypes_Britain, time_points_Britain,
                                          const, genotypes_NCentEur, time_points_NCentEur,
                                          const, genotypes_South, time_points_South,
                                          const, genotypes_Baltic, time_points_Baltic),
                       trace=TRUE, NP=200, fnscale=100, tol=5e-7, maxiter=400)
  results = optimise_all('results/optim_JDE_additive_1.Robj', res_const, optimise_JDEoptim, logistic_fct_sel)
} else {
  load(file = 'results/optim_JDE_additive_1.Robj')
}
res_JDE_additive_1 = results
#
if (!file.exists('results/optim_JDE_additive_2.Robj')) {
  res_const = JDEoptim(lower=c(0, 0), upper=c(0.05, 0.2),
                       fn=function(params)
                         -full_likelihood(params[1], params[2], 1, logistic_fct_sel,
                                          const, genotypes_Britain, time_points_Britain,
                                          const, genotypes_NCentEur, time_points_NCentEur,
                                          const, genotypes_South, time_points_South,
                                          const, genotypes_Baltic, time_points_Baltic),
                       trace=TRUE, NP=200, fnscale=100, tol=5e-7, maxiter=400)
  results = optimise_all('results/optim_JDE_additive_2.Robj', res_const, optimise_JDEoptim, logistic_fct_sel)
} else {
  load(file = 'results/optim_JDE_additive_2.Robj')
}
res_JDE_additive_2 = results
#
if (!file.exists('results/optim_JDE_additive_3.Robj')) {
  res_const = JDEoptim(lower=c(0, 0), upper=c(0.05, 0.2),
                       fn=function(params)
                         -full_likelihood(params[1], params[2], 1, logistic_fct_sel,
                                          const, genotypes_Britain, time_points_Britain,
                                          const, genotypes_NCentEur, time_points_NCentEur,
                                          const, genotypes_South, time_points_South,
                                          const, genotypes_Baltic, time_points_Baltic),
                       trace=TRUE, NP=200, fnscale=100, tol=5e-7, maxiter=400)
  results = optimise_all('results/optim_JDE_additive_3.Robj', res_const, optimise_JDEoptim, logistic_fct_sel)
} else {
  load(file = 'results/optim_JDE_additive_3.Robj')
}
res_JDE_additive_3 = results
#
#
if (!file.exists('results/optim_JDE_additive_noInit_1.Robj')) {
  res_const = JDEoptim(lower=c(0, 0), upper=c(0.05, 0.2),
                       fn=function(params)
                         -full_likelihood(params[1], params[2], 1, logistic_fct_sel,
                                          const, genotypes_Britain, time_points_Britain,
                                          const, genotypes_NCentEur, time_points_NCentEur,
                                          const, genotypes_South, time_points_South,
                                          const, genotypes_Baltic, time_points_Baltic),
                       trace=TRUE, NP=200, fnscale=100, tol=5e-7, maxiter=400)
  results = optimise_all('results/optim_JDE_additive_noInit_1.Robj', res_const, optimise_JDEoptim, logistic_fct_sel,
                         init=FALSE)
} else {
  load(file = 'results/optim_JDE_additive_noInit_1.Robj')
}
res_JDE_additive_noInit_1 = results
#
if (!file.exists('results/optim_JDE_additive_noInit_2.Robj')) {
  res_const = JDEoptim(lower=c(0, 0), upper=c(0.05, 0.2),
                       fn=function(params)
                         -full_likelihood(params[1], params[2], 1, logistic_fct_sel,
                                          const, genotypes_Britain, time_points_Britain,
                                          const, genotypes_NCentEur, time_points_NCentEur,
                                          const, genotypes_South, time_points_South,
                                          const, genotypes_Baltic, time_points_Baltic),
                       trace=TRUE, NP=200, fnscale=100, tol=5e-7, maxiter=400)
  results = optimise_all('results/optim_JDE_additive_noInit_2.Robj', res_const, optimise_JDEoptim, logistic_fct_sel,
                         init=FALSE)
} else {
  load(file = 'results/optim_JDE_additive_noInit_2.Robj')
}
res_JDE_additive_noInit_2 = results
#
if (!file.exists('results/optim_JDE_additive_noInit_3.Robj')) {
  res_const = JDEoptim(lower=c(0, 0), upper=c(0.05, 0.2),
                       fn=function(params)
                         -full_likelihood(params[1], params[2], 1, logistic_fct_sel,
                                          const, genotypes_Britain, time_points_Britain,
                                          const, genotypes_NCentEur, time_points_NCentEur,
                                          const, genotypes_South, time_points_South,
                                          const, genotypes_Baltic, time_points_Baltic),
                       trace=TRUE, NP=200, fnscale=100, tol=5e-7, maxiter=400)
  results = optimise_all('results/optim_JDE_additive_noInit_3.Robj', res_const, optimise_JDEoptim, logistic_fct_sel,
                         init=FALSE)
} else {
  load(file = 'results/optim_JDE_additive_noInit_3.Robj')
}
res_JDE_additive_noInit_3 = results
#
#
#
if (!file.exists('results/optim_JDE_dominant_1.Robj')) {
  res_const_dominant = JDEoptim(lower=c(0, 0), upper=c(0.05, 0.2),
                                fn=function(params)
                                  -full_likelihood(params[1], params[2], 1, dominant_sel,
                                                   const, genotypes_Britain, time_points_Britain,
                                                   const, genotypes_NCentEur, time_points_NCentEur,
                                                   const, genotypes_South, time_points_South,
                                                   const, genotypes_Baltic, time_points_Baltic),
                                trace=TRUE, NP=200, fnscale=100, tol=5e-7, maxiter=400)
  results = optimise_all('results/optim_JDE_dominant_1.Robj', res_const_dominant, optimise_JDEoptim, dominant_sel)
} else {
  load(file = 'results/optim_JDE_dominant_1.Robj')
}
res_JDE_dominant_1 = results
#
if (!file.exists('results/optim_JDE_dominant_2.Robj')) {
  res_const_dominant = JDEoptim(lower=c(0, 0), upper=c(0.05, 0.2),
                                fn=function(params)
                                  -full_likelihood(params[1], params[2], 1, dominant_sel,
                                                   const, genotypes_Britain, time_points_Britain,
                                                   const, genotypes_NCentEur, time_points_NCentEur,
                                                   const, genotypes_South, time_points_South,
                                                   const, genotypes_Baltic, time_points_Baltic),
                                trace=TRUE, NP=200, fnscale=100, tol=5e-7, maxiter=400)
  results = optimise_all('results/optim_JDE_dominant_2.Robj', res_const_dominant, optimise_JDEoptim, dominant_sel)
} else {
  load(file = 'results/optim_JDE_dominant_2.Robj')
}
res_JDE_dominant_2 = results
#
if (!file.exists('results/optim_JDE_dominant_3.Robj')) {
  res_const_dominant = JDEoptim(lower=c(0, 0), upper=c(0.05, 0.2),
                                fn=function(params)
                                  -full_likelihood(params[1], params[2], 1, dominant_sel,
                                                   const, genotypes_Britain, time_points_Britain,
                                                   const, genotypes_NCentEur, time_points_NCentEur,
                                                   const, genotypes_South, time_points_South,
                                                   const, genotypes_Baltic, time_points_Baltic),
                                trace=TRUE, NP=200, fnscale=100, tol=5e-7, maxiter=400)
  results = optimise_all('results/optim_JDE_dominant_3.Robj', res_const_dominant, optimise_JDEoptim, dominant_sel)
} else {
  load(file = 'results/optim_JDE_dominant_3.Robj')
}
res_JDE_dominant_3 = results
#
#
if (!file.exists('results/optim_JDE_dominant_noInit_1.Robj')) {
  res_const_dominant = JDEoptim(lower=c(0, 0), upper=c(0.05, 0.2),
                                fn=function(params)
                                  -full_likelihood(params[1], params[2], 1, dominant_sel,
                                                   const, genotypes_Britain, time_points_Britain,
                                                   const, genotypes_NCentEur, time_points_NCentEur,
                                                   const, genotypes_South, time_points_South,
                                                   const, genotypes_Baltic, time_points_Baltic),
                                trace=TRUE, NP=200, fnscale=100, tol=5e-7, maxiter=400)
  results = optimise_all('results/optim_JDE_dominant_noInit_1.Robj', res_const_dominant, optimise_JDEoptim,
                         dominant_sel, init=FALSE)
} else {
  load(file = 'results/optim_JDE_dominant_noInit_1.Robj')
}
res_JDE_dominant_noInit_1 = results
#
if (!file.exists('results/optim_JDE_dominant_noInit_2.Robj')) {
  res_const_dominant = JDEoptim(lower=c(0, 0), upper=c(0.05, 0.2),
                                fn=function(params)
                                  -full_likelihood(params[1], params[2], 1, dominant_sel,
                                                   const, genotypes_Britain, time_points_Britain,
                                                   const, genotypes_NCentEur, time_points_NCentEur,
                                                   const, genotypes_South, time_points_South,
                                                   const, genotypes_Baltic, time_points_Baltic),
                                trace=TRUE, NP=200, fnscale=100, tol=5e-7, maxiter=400)
  results = optimise_all('results/optim_JDE_dominant_noInit_2.Robj', res_const_dominant, optimise_JDEoptim,
                         dominant_sel, init=FALSE)
} else {
  load(file = 'results/optim_JDE_dominant_noInit_2.Robj')
}
res_JDE_dominant_noInit_2 = results
#
if (!file.exists('results/optim_JDE_dominant_noInit_3.Robj')) {
  res_const_dominant = JDEoptim(lower=c(0, 0), upper=c(0.05, 0.2),
                                fn=function(params)
                                  -full_likelihood(params[1], params[2], 1, dominant_sel,
                                                   const, genotypes_Britain, time_points_Britain,
                                                   const, genotypes_NCentEur, time_points_NCentEur,
                                                   const, genotypes_South, time_points_South,
                                                   const, genotypes_Baltic, time_points_Baltic),
                                trace=TRUE, NP=200, fnscale=100, tol=5e-7, maxiter=400)
  results = optimise_all('results/optim_JDE_dominant_noInit_3.Robj', res_const_dominant, optimise_JDEoptim,
                         dominant_sel, init=FALSE)
} else {
  load(file = 'results/optim_JDE_dominant_noInit_3.Robj')
}
res_JDE_dominant_noInit_3 = results
#
#
#
chooser <- function(list_of_res) {
  vector_of_lLs = as.vector(sapply(list_of_res, function(x) -x$value))
  max_val_index = which(vector_of_lLs==max(vector_of_lLs))
  list_of_res[[max_val_index[1]]] # always choose first element which achieves maximum
}
#
res_const = chooser(list('1' = res_JDE_additive_1$const, '2' = res_JDE_additive_2$const,
                         '3' = res_JDE_additive_3$const,
                         '4' = res_JDE_additive_noInit_1$const, '5' = res_JDE_additive_noInit_2$const,
                         '6' = res_JDE_additive_noInit_3$const))
res_cluster = chooser(list('1' = res_JDE_additive_1$cluster, '2' = res_JDE_additive_2$cluster,
                           '3' = res_JDE_additive_3$cluster,
                           '4' = res_JDE_additive_noInit_1$cluster, '5' = res_JDE_additive_noInit_2$cluster,
                           '6' = res_JDE_additive_noInit_3$cluster))
res_cluster_inv = chooser(list('1' = res_JDE_additive_1$cluster_inv, '2' = res_JDE_additive_2$cluster_inv,
                               '3' = res_JDE_additive_3$cluster_inv,
                               '4' = res_JDE_additive_noInit_1$cluster_inv, '5' = res_JDE_additive_noInit_2$cluster_inv,
                               '6' = res_JDE_additive_noInit_3$cluster_inv))
res_domestics = chooser(list('1' = res_JDE_additive_1$domestics, '2' = res_JDE_additive_2$domestics,
                             '3' = res_JDE_additive_3$domestics,
                             '4' = res_JDE_additive_noInit_1$domestics, '5' = res_JDE_additive_noInit_2$domestics,
                             '6' = res_JDE_additive_noInit_3$domestics))
res_domestics_inv = chooser(list('1' = res_JDE_additive_1$domestics_inv, '2' = res_JDE_additive_2$domestics_inv,
                                 '3' = res_JDE_additive_3$domestics_inv,
                                 '4' = res_JDE_additive_noInit_1$domestics_inv,
                                 '5' = res_JDE_additive_noInit_2$domestics_inv,
                                 '6' = res_JDE_additive_noInit_3$domestics_inv))
res_milk = chooser(list('1' = res_JDE_additive_1$milk, '2' = res_JDE_additive_2$milk,
                        '3' = res_JDE_additive_3$milk,
                        '4' = res_JDE_additive_noInit_1$milk, '5' = res_JDE_additive_noInit_2$milk,
                        '6' = res_JDE_additive_noInit_3$milk))
res_pop_fluc_globally = chooser(list('1' = res_JDE_additive_1$pop_fluc_globally,
                                     '2' = res_JDE_additive_2$pop_fluc_globally,
                                     '3' = res_JDE_additive_3$pop_fluc_globally,
                                     '4' = res_JDE_additive_noInit_1$pop_fluc_globally,
                                     '5' = res_JDE_additive_noInit_2$pop_fluc_globally,
                                     '6' = res_JDE_additive_noInit_3$pop_fluc_globally))
res_pop_fluc_inv_globally = chooser(list('1' = res_JDE_additive_1$pop_fluc_inv_globally,
                                         '2' = res_JDE_additive_2$pop_fluc_inv_globally,
                                         '3' = res_JDE_additive_3$pop_fluc_inv_globally,
                                         '4' = res_JDE_additive_noInit_1$pop_fluc_inv_globally,
                                         '5' = res_JDE_additive_noInit_2$pop_fluc_inv_globally,
                                         '6' = res_JDE_additive_noInit_3$pop_fluc_inv_globally))
res_insolation_globally = chooser(list('1' = res_JDE_additive_1$insolation_globally,
                                       '2' = res_JDE_additive_2$insolation_globally,
                                       '3' = res_JDE_additive_3$insolation_globally,
                                       '4' = res_JDE_additive_noInit_1$insolation_globally,
                                       '5' = res_JDE_additive_noInit_2$insolation_globally,
                                       '6' = res_JDE_additive_noInit_3$insolation_globally))
res_insolation_inv_globally = chooser(list('1' = res_JDE_additive_1$insolation_inv_globally,
                                           '2' = res_JDE_additive_2$insolation_inv_globally,
                                           '3' = res_JDE_additive_3$insolation_inv_globally,
                                           '1' = res_JDE_additive_noInit_1$insolation_inv_globally,
                                           '2' = res_JDE_additive_noInit_2$insolation_inv_globally,
                                           '3' = res_JDE_additive_noInit_3$insolation_inv_globally))
#
res_const_dominant = chooser(list('1' = res_JDE_dominant_1$const, '2' = res_JDE_dominant_2$const,
                                  '3' = res_JDE_dominant_3$const,
                                  '4' = res_JDE_dominant_noInit_1$const, '5' = res_JDE_dominant_noInit_2$const,
                                  '6' = res_JDE_dominant_noInit_3$const))
res_cluster_dominant = chooser(list('1' = res_JDE_dominant_1$cluster, '2' = res_JDE_dominant_2$cluster,
                                    '3' = res_JDE_dominant_3$cluster,
                                    '4' = res_JDE_dominant_noInit_1$cluster, '5' = res_JDE_dominant_noInit_2$cluster,
                                    '6' = res_JDE_dominant_noInit_3$cluster))
res_cluster_inv_dominant = chooser(list('1' = res_JDE_dominant_1$cluster_inv, '2' = res_JDE_dominant_2$cluster_inv,
                                        '3' = res_JDE_dominant_3$cluster_inv,
                                        '4' = res_JDE_dominant_noInit_1$cluster_inv,
                                        '5' = res_JDE_dominant_noInit_2$cluster_inv,
                                        '6' = res_JDE_dominant_noInit_3$cluster_inv))
res_domestics_dominant = chooser(list('1' = res_JDE_dominant_1$domestics, '2' = res_JDE_dominant_2$domestics,
                                      '3' = res_JDE_dominant_3$domestics,
                                      '4' = res_JDE_dominant_noInit_1$domestics,
                                      '5' = res_JDE_dominant_noInit_2$domestics,
                                      '6' = res_JDE_dominant_noInit_3$domestics))
res_domestics_inv_dominant = chooser(list('1' = res_JDE_dominant_1$domestics_inv,
                                          '2' = res_JDE_dominant_2$domestics_inv,
                                          '3' = res_JDE_dominant_3$domestics_inv,
                                          '4' = res_JDE_dominant_noInit_1$domestics_inv,
                                          '5' = res_JDE_dominant_noInit_2$domestics_inv,
                                          '6' = res_JDE_dominant_noInit_3$domestics_inv))
res_milk_dominant = chooser(list('1' = res_JDE_dominant_1$milk, '2' = res_JDE_dominant_2$milk,
                                 '3' = res_JDE_dominant_3$milk,
                                 '4' = res_JDE_dominant_noInit_1$milk, '5' = res_JDE_dominant_noInit_2$milk,
                                 '6' = res_JDE_dominant_noInit_3$milk))
res_pop_fluc_dominant_globally = chooser(list('1' = res_JDE_dominant_1$pop_fluc_globally,
                                              '2' = res_JDE_dominant_2$pop_fluc_globally,
                                              '3' = res_JDE_dominant_3$pop_fluc_globally,
                                              '4' = res_JDE_dominant_noInit_1$pop_fluc_globally,
                                              '5' = res_JDE_dominant_noInit_2$pop_fluc_globally,
                                              '6' = res_JDE_dominant_noInit_3$pop_fluc_globally))
res_pop_fluc_inv_dominant_globally = chooser(list('1' = res_JDE_dominant_1$pop_fluc_inv_globally,
                                                  '2' = res_JDE_dominant_2$pop_fluc_inv_globally,
                                                  '3' = res_JDE_dominant_3$pop_fluc_inv_globally,
                                                  '4' = res_JDE_dominant_noInit_1$pop_fluc_inv_globally,
                                                  '5' = res_JDE_dominant_noInit_2$pop_fluc_inv_globally,
                                                  '6' = res_JDE_dominant_noInit_3$pop_fluc_inv_globally))
res_insolation_dominant_globally = chooser(list('1' = res_JDE_dominant_1$insolation_globally,
                                                '2' = res_JDE_dominant_2$insolation_globally,
                                                '3' = res_JDE_dominant_3$insolation_globally,
                                                '4' = res_JDE_dominant_noInit_1$insolation_globally,
                                                '5' = res_JDE_dominant_noInit_2$insolation_globally,
                                                '6' = res_JDE_dominant_noInit_3$insolation_globally))
res_insolation_inv_dominant_globally = chooser(list('1' = res_JDE_dominant_1$insolation_inv_globally,
                                                    '2' = res_JDE_dominant_2$insolation_inv_globally,
                                                    '3' = res_JDE_dominant_3$insolation_inv_globally,
                                                    '4' = res_JDE_dominant_noInit_1$insolation_inv_globally,
                                                    '5' = res_JDE_dominant_noInit_2$insolation_inv_globally,
                                                    '6' = res_JDE_dominant_noInit_3$insolation_inv_globally))
#
#
#
# *****************************
# ***  FUNCTIONS, PART III  ***
# *****************************
#
plot_power_analysis_lLs <- function(out_infix, results) {
  figure_name = paste0('results/power_analysis_', out_infix, '_lLs.pdf')
  if (file.exists(figure_name)) { return() }
  #
  pdf(file = figure_name, width=5.6, height=5)
  par(mar=c(9, 4, 3, 3))
  #
  ys_additive = -c(results$const$value, results$milk$value,
                   results$insolation_globally$value, results$insolation_inv_globally$value,
                   results$pop_fluc_globally$value, results$pop_fluc_inv_globally$value,
                   results$cluster$value, results$cluster_inv$value,
                   results$domestics$value, results$domestics_inv$value)
  #
  plot(ys_additive, pch='-', col='black', cex=3, ylim = range(ys_additive), xaxt='n', xlab='', ylab='log likelihood')
  #
  lines(par('usr')[1:2], c(-results$const$value + qchisq(0.05, df=1, lower.tail = FALSE) / 2,
                           -results$const$value + qchisq(0.05, df=1, lower.tail = FALSE) / 2), col='red')
  #
  axis(side = 1, at = 1:10, labels = c('constant', 'milk', 'insolation', 'insolation inv.',
                                       'population fluc.', 'population fluc. inv.', 'cluster', 'cluster inv.', 
                                       'domesticates', 'wild'), las=2, cex=0.8)
  #
  par(xpd=TRUE)
  legend('topright', c('additive', 'significance'), 
         pch=c('-', NA), pt.cex = c(3, 1), lty=c(NA, 1),
         col=c('black', 'red'), inset = c(-0.10, -0.2), bg = 'white', cex = 0.8)
  #
  dev.off()
}
#
#
#
.plot_region <- function(region_name, time_points, genotypes, ms,
                         ss_breaks, der_al_freq, max_sel, max_freq, drifted_der_al_freq=NA, legend=FALSE) {
  par(mar=c(5, 4, 4, 4))
  #
  plot(as.vector(rbind(ms$startBP, ms$endBP)),
       as.vector(rbind(ss_breaks$ss, ss_breaks$ss)),
       xlim=c(8000, 2500), ylim=c(0, max_sel), 
       type='l', lwd=1, xlab='', ylab='',
       main=paste(region_name))
  #
  lp = .convert_genotypes_to_binary(time_points, genotypes)
  lp$cols = 'firebrick'
  der_al_indices = which(lp$allele==1)
  if (length(der_al_indices)>0) {
      lp[which(lp$allele==1), ]$cols = 'steelblue'
    lp[which(lp$allele==1), ]$allele = max_sel   
  }
  #
  points(lp$mid, lp$allele, pch=16, col=lp$cols, cex=0.8)
  #
  par(new = T)
  plot(der_al_freq, type='l', col = 'darkgreen', ylim=c(0, max_freq), axes=F, xlab=NA, ylab=NA)
  lines(drifted_der_al_freq, col = 'darkgreen', lwd=0.7) 
  axis(side = 4)
  #
  if (legend==TRUE) {
    legend('topleft', c('Lactase non-persistence allele', 'Lactase persistence allele'),
           pch=16, cex=0.8, col=c('firebrick', 'steelblue')) }
}
#
#
#
plot_result <- function(driver_name, selection_name, fct_selection, res,
                        time_points_Britain, genotypes_Britain, ms_Britain, 
                        time_points_NCentEur, genotypes_NCentEur, ms_NCentEur, 
                        time_points_South, genotypes_South, ms_South,
                        time_points_Baltic, genotypes_Baltic, ms_Baltic,
                        drifted_der_al_freq_Britain=NA,
                        drifted_der_al_freq_NCentEur=NA,
                        drifted_der_al_freq_South=NA, 
                        drifted_der_al_freq_Baltic=NA) {
  #
  figure_name = paste('results/', driver_name, '_', selection_name, '_res.pdf', sep='')
  if (file.exists(figure_name)) { return() }
  #
  #
  pdf(file = figure_name, width = 10)
  #
  par(mfrow=c(2,2))
  #
  par(oma=c(2, 4, 4, 4))
  #
  initial_freq = res$par[1]
  s_1 = res$par[2]
  a = res$par[3]
  #
  ss_breaks_Britain = ms_to_ss_breaks(ms_Britain, s_1, a)
  der_al_freq_Britain =
              ss_breaks_to_der_al_freq(ss_breaks_Britain$ss, ss_breaks_Britain$breaks, fct_selection, y0=initial_freq)
  ss_breaks_NCentEur = ms_to_ss_breaks(ms_NCentEur, s_1, a)
  der_al_freq_NCentEur =
              ss_breaks_to_der_al_freq(ss_breaks_NCentEur$ss, ss_breaks_NCentEur$breaks, fct_selection, y0=initial_freq)
  ss_breaks_South = ms_to_ss_breaks(ms_South, s_1, a)
  der_al_freq_South =
              ss_breaks_to_der_al_freq(ss_breaks_South$ss, ss_breaks_South$breaks, fct_selection, y0=initial_freq)
  ss_breaks_Baltic = ms_to_ss_breaks(ms_Baltic, s_1, a)
  der_al_freq_Baltic =
              ss_breaks_to_der_al_freq(ss_breaks_Baltic$ss, ss_breaks_Baltic$breaks, fct_selection, y0=initial_freq)
  #
  max_sel = max(c(ss_breaks_Britain$ss, ss_breaks_NCentEur$ss, ss_breaks_South$ss, ss_breaks_Baltic$ss))
  max_freq = max(c(der_al_freq_Britain, der_al_freq_NCentEur, der_al_freq_South, der_al_freq_Baltic,
                   drifted_der_al_freq_Britain, drifted_der_al_freq_NCentEur, drifted_der_al_freq_South,
                   drifted_der_al_freq_Baltic), na.rm = TRUE)
  #
  .plot_region('Britain', time_points_Britain, genotypes_Britain, ms_Britain, 
               ss_breaks_Britain, der_al_freq_Britain, max_sel, max_freq, drifted_der_al_freq_Britain)
  .plot_region('North and Central Europe', time_points_NCentEur, genotypes_NCentEur, ms_NCentEur, 
               ss_breaks_NCentEur, der_al_freq_NCentEur, max_sel, max_freq, drifted_der_al_freq_NCentEur)
  .plot_region('South', time_points_South, genotypes_South, ms_South, 
               ss_breaks_South, der_al_freq_South, max_sel, max_freq, drifted_der_al_freq_South, legend=TRUE)
  .plot_region('Baltic', time_points_Baltic, genotypes_Baltic, ms_Baltic, 
               ss_breaks_Baltic, der_al_freq_Baltic, max_sel, max_freq, drifted_der_al_freq_Baltic)
  #
  mtext(side = 1, line = 0, text = 'year [BP]', outer = T)
  mtext(side = 2, line = 1, text = 'selection coefficient', outer = T)
  mtext(side = 3, line = 1, text = 
                                TeX(sprintf(r'(\textbf{%s, %s ($f_0=%.4f, s=%.3f, a=%.3f$)})',
                                            driver_name, selection_name,
                                            round(initial_freq, digits = 4), round(s_1, 3), round(a, 3))),
        font=2, outer = T)
  mtext(side = 4, line = 1, 'predicted LP allele frequency', col='darkgreen', outer = T)
  #
  dev.off()
}
#
#
#
sample_at_time_points <- function(time_points, genotypes, derived_allele_freqs, start=8000, end=2500) {
  #
  stopifnot(length(time_points)==length(genotypes))
  #
  indices = abs(time_points - start + 1)
  #
  sampled_genotypes = vector(length = length(time_points))
  for (i in 1:length(time_points)) {
    sampled_genotypes[i] = paste0(sample(c('G', 'A'), nchar(genotypes[i]), replace=TRUE,
                                         prob=c(1 - derived_allele_freqs[indices[i]], derived_allele_freqs[indices[i]])),
                                  collapse='')
  }
  return(sampled_genotypes)
}
#
if (!file.exists('results/power_analysis_milk_noise.Robj')) {
  milk_normalised_10pc = cbind(milk_normalised[, 1:2], apply(X = milk_normalised[, 3:6], MARGIN = 2,
                                          FUN = function(x) x * runif(n = nrow(milk_normalised), min = 0.9, max = 1.1)))
  milk_normalised_20pc = cbind(milk_normalised[, 1:2], apply(X = milk_normalised[, 3:6], MARGIN = 2,
                                          FUN = function(x) x * runif(n = nrow(milk_normalised), min = 0.8, max = 1.2)))
  milk_normalised_30pc = cbind(milk_normalised[, 1:2], apply(X = milk_normalised[, 3:6], MARGIN = 2,
                                          FUN = function(x) x * runif(n = nrow(milk_normalised), min = 0.7, max = 1.3)))
  milk_normalised_50pc = cbind(milk_normalised[, 1:2], apply(X = milk_normalised[, 3:6], MARGIN = 2,
                                          FUN = function(x) x * runif(n = nrow(milk_normalised), min = 0.5, max = 1.5)))
  input_milk_noise = list('milk_10pc' = milk_normalised_10pc, 'milk_20pc' = milk_normalised_20pc, 
                          'milk_30pc' = milk_normalised_30pc, 'milk_50pc' = milk_normalised_50pc)
  save(input_milk_noise, file = 'results/power_analysis_milk_noise.Robj')
} else {
  load(file = 'results/power_analysis_milk_noise.Robj')
}
#
power_analysis_milk <- function(out_infix, fct_selection, initial_freq, s_1, a, N_e=0) {
  #
  sample_genotypes <- function(milk_normalised_region, time_points, genotypes) {
    ss_breaks = ms_to_ss_breaks(milk_normalised_region, s_1, a)
    der_al_freq = ss_breaks_to_der_al_freq(ss_breaks$ss, ss_breaks$breaks, fct_selection, y0=initial_freq)
    if (N_e > 0) {
      drifted_der_al_freq = ss_breaks_to_der_al_freq(ss_breaks$ss, ss_breaks$breaks, fct_selection,
                                                     y0=initial_freq, N_e=N_e)
      sampled_genotypes = sample_at_time_points(time_points, genotypes, drifted_der_al_freq)
    } else {
      drifted_der_al_freq = NA
      sampled_genotypes = sample_at_time_points(time_points, genotypes, der_al_freq)
    }
    return(list('drifted_der_al_freq' = drifted_der_al_freq, 'sampled_genotypes' = sampled_genotypes))
  }
  #
  if (!file.exists(file = paste0('results/power_analysis_', out_infix, '_genotypes_freqs.Robj'))) {
    tmp = sample_genotypes(milk_normalised[, c(1, 2, 3)], time_points_Britain, genotypes_Britain)
    drifted_der_al_freq_Britain = tmp$drifted_der_al_freq
    sampled_genotypes_Britain = tmp$sampled_genotypes
    #
    tmp = sample_genotypes(milk_normalised[, c(1, 2, 5)], time_points_NCentEur, genotypes_NCentEur)
    drifted_der_al_freq_NCentEur = tmp$drifted_der_al_freq
    sampled_genotypes_NCentEur = tmp$sampled_genotypes
    #
    tmp = sample_genotypes(milk_normalised[, c(1, 2, 6)], time_points_South, genotypes_South)
    drifted_der_al_freq_South = tmp$drifted_der_al_freq
    sampled_genotypes_South = tmp$sampled_genotypes
    #
    tmp = sample_genotypes(milk_normalised[, c(1, 2, 4)], time_points_Baltic, genotypes_Baltic)
    drifted_der_al_freq_Baltic = tmp$drifted_der_al_freq
    sampled_genotypes_Baltic = tmp$sampled_genotypes
    #
    save(drifted_der_al_freq_Britain, drifted_der_al_freq_NCentEur, drifted_der_al_freq_South,
         drifted_der_al_freq_Baltic,
         sampled_genotypes_Britain, sampled_genotypes_NCentEur, sampled_genotypes_South, sampled_genotypes_Baltic,
         file = paste0('results/power_analysis_', out_infix, '_genotypes_freqs.Robj'))
  } else {
    load(file = paste0('results/power_analysis_', out_infix, '_genotypes_freqs.Robj'))
  }
  #
  if (!file.exists(paste0('results/power_analysis_', out_infix, '.Robj'))) {
    res_const = JDEoptim(lower=c(0, 0), upper=c(0.05, 0.2),
                         fn=function(params)
                           -full_likelihood(params[1], params[2], 1, fct_selection,
                                            const, sampled_genotypes_Britain, time_points_Britain,
                                            const, sampled_genotypes_NCentEur, time_points_NCentEur,
                                            const, sampled_genotypes_South, time_points_South,
                                            const, sampled_genotypes_Baltic, time_points_Baltic),
                         trace=TRUE, NP=200, fnscale=100, tol=5e-7, maxiter=400)
    #
    results <- mclapply(input, optimise_JDEoptim,
                        fct_selection=fct_selection, indices_regional_MAPs=c(3, 5, 6, 4),
                        genotypes_Britain=sampled_genotypes_Britain, time_points_Britain=time_points_Britain,
                        genotypes_NCentEur=sampled_genotypes_NCentEur, time_points_NCentEur=time_points_NCentEur,
                        genotypes_South=sampled_genotypes_South, time_points_South=time_points_South,
                        genotypes_Baltic=sampled_genotypes_Baltic, time_points_Baltic=time_points_Baltic,
                        init_f_0=res_const$par[1], init_s=res_const$par[2],
                        mc.cores = numCores)
    #
    results[['const']] = res_const
    #
    save(results, file = paste0('results/power_analysis_', out_infix, '.Robj'))
  } else {
    load(file = paste0('results/power_analysis_', out_infix, '.Robj'))
  }
  #
  plot_result('power_analysis', out_infix, fct_selection, results$milk,
              time_points_Britain, sampled_genotypes_Britain, milk_normalised[, c(1, 2, 3)],
              time_points_NCentEur, sampled_genotypes_NCentEur, milk_normalised[, c(1, 2, 5)], 
              time_points_South, sampled_genotypes_South, milk_normalised[, c(1, 2, 6)],
              time_points_Baltic, sampled_genotypes_Baltic, milk_normalised[, c(1, 2, 4)], 
              drifted_der_al_freq_Britain=drifted_der_al_freq_Britain,
              drifted_der_al_freq_NCentEur=drifted_der_al_freq_NCentEur,
              drifted_der_al_freq_South=drifted_der_al_freq_South, 
              drifted_der_al_freq_Baltic=drifted_der_al_freq_Baltic)
}
#
#
#
# *****************************
# ***        PART III       ***
# *****************************
#
clean_version <- function(x) {
  # https://www.roelpeters.be/how-to-prevent-scientific-notation-in-r/
  str_remove(format(x, scientific=FALSE), "[.]")
}
#
power_analysis = data.frame()
#
counter = 1
for (a in c(0.25, 0.6, 0.8)) {
  for (s_1 in c(0.02, 0.035, 0.05)) {
    for (initial_freq_N_e in list(c(0.0005, 0), c(0.005, 0), c(0.01, 0),
                                  c(0.005, 10000), c(0.005, 5000), c(0.005, 1000), c(0.005, 500))) {
      N_e = initial_freq_N_e[2]
      #
      for (rep in 1:3) {
        out_infix = paste0('additive_f', clean_version(initial_freq_N_e[1]), '_s',
                           clean_version(s_1), '_a', clean_version(a), '_rep', rep,
                           if (N_e > 0) paste0('_Ne', clean_version(N_e)) else '')
        print(out_infix)
        #
        power_analysis[counter, c('a', 's_1', 'initial_freq', 'N_e', 'rep')] = c(a, s_1, initial_freq_N_e[1], N_e, rep)
        #
        power_analysis_milk(out_infix, fct_selection = logistic_fct_sel,
                            initial_freq = initial_freq_N_e[1], s_1 = s_1, a = a, N_e = N_e)
        #
        load(file = paste0('results/power_analysis_', out_infix, '_genotypes_freqs.Robj'))
        load(file = paste0('results/power_analysis_', out_infix, '.Robj'))
        plot_power_analysis_lLs(out_infix, results)
        #
        nb_der_allele = str_count(paste0(c(sampled_genotypes_Britain, sampled_genotypes_NCentEur,
                                           sampled_genotypes_South, sampled_genotypes_Baltic), collapse=''), 'A')
        power_analysis[counter, c('nb_der_allele', 'est_a', 'est_s_1', 'est_initial_freq', 'lL_const', 'lL_milk')] =
                                c(nb_der_allele, results$milk$par[3], results$milk$par[2], results$milk$par[1],
                                  results$const$value, results$milk$value)
        #
        counter = counter + 1
      }
    }
  }
}
#
#
#
pa = power_analysis
f = 5e-03
#
#
#
for (N_e in c(0, 10000, 5000, 1000, 500)) {
  #
  pdf(paste0('results/power_analysis_s_1_a_', if (N_e==0) 'no_drift' else paste0('Ne', N_e), '.pdf'), height=9)
  #
  par(mfrow=c(2, 1))
  #
  par(oma=c(0, 0, 4, 0))
  par(mar=c(3, 4, 1, 2))
  #
  plot(c(rep(0.7, 3), rep(1, 3), rep(1.3, 3), rep(1.7, 3), rep(2, 3), rep(2.3, 3), rep(2.7, 3), rep(3, 3), rep(3.3, 3)),
       c(pa[pa$a==0.25 & pa$s_1==0.02 & pa$initial_freq==f & pa$N_e==N_e, 'est_s_1'],
         pa[pa$a==0.25 & pa$s_1==0.035 & pa$initial_freq==f & pa$N_e==N_e, 'est_s_1'],
         pa[pa$a==0.25 & pa$s_1==0.05 & pa$initial_freq==f & pa$N_e==N_e, 'est_s_1'],
         pa[pa$a==0.6 & pa$s_1==0.02 & pa$initial_freq==f & pa$N_e==N_e, 'est_s_1'],
         pa[pa$a==0.6 & pa$s_1==0.035 & pa$initial_freq==f & pa$N_e==N_e, 'est_s_1'],
         pa[pa$a==0.6 & pa$s_1==0.05 & pa$initial_freq==f & pa$N_e==N_e, 'est_s_1'],
         pa[pa$a==0.8 & pa$s_1==0.02 & pa$initial_freq==f & pa$N_e==N_e, 'est_s_1'],
         pa[pa$a==0.8 & pa$s_1==0.035 & pa$initial_freq==f & pa$N_e==N_e, 'est_s_1'],
         pa[pa$a==0.8 & pa$s_1==0.05 & pa$initial_freq==f & pa$N_e==N_e, 'est_s_1']),
    pch=16, col='darkgrey',
    xlim = c(0.5, 3.5), 
    xaxt='n', 
    xlab='', ylab='selection coefficient s')
  #
  segments(c(0.6, 0.9, 1.2, 1.6, 1.9, 2.2, 2.6, 2.9, 3.2), rep(c(0.02, 0.035, 0.05), 3),
           c(0.8, 1.1, 1.4, 1.8, 2.1, 2.4, 2.8, 3.1, 3.4), rep(c(0.02, 0.035, 0.05), 3))
  #
  ys = c(mean(pa[pa$a==0.25 & pa$s_1==0.02 & pa$initial_freq==f & pa$N_e==N_e, 'est_s_1']),
         mean(pa[pa$a==0.25 & pa$s_1==0.035 & pa$initial_freq==f & pa$N_e==N_e, 'est_s_1']),
         mean(pa[pa$a==0.25 & pa$s_1==0.05 & pa$initial_freq==f & pa$N_e==N_e, 'est_s_1']),
         mean(pa[pa$a==0.6 & pa$s_1==0.02 & pa$initial_freq==f & pa$N_e==N_e, 'est_s_1']),
         mean(pa[pa$a==0.6 & pa$s_1==0.035 & pa$initial_freq==f & pa$N_e==N_e, 'est_s_1']),
         mean(pa[pa$a==0.6 & pa$s_1==0.05 & pa$initial_freq==f & pa$N_e==N_e, 'est_s_1']),
         mean(pa[pa$a==0.8 & pa$s_1==0.02 & pa$initial_freq==f & pa$N_e==N_e, 'est_s_1']),
         mean(pa[pa$a==0.8 & pa$s_1==0.035 & pa$initial_freq==f & pa$N_e==N_e, 'est_s_1']),
         mean(pa[pa$a==0.8 & pa$s_1==0.05 & pa$initial_freq==f & pa$N_e==N_e, 'est_s_1']))
  #
  segments(c(0.65, 0.95, 1.25, 1.65, 1.95, 2.25, 2.65, 2.95, 3.25), ys,
           c(0.75, 1.05, 1.35, 1.75, 2.05, 2.35, 2.75, 3.05, 3.35), ys, col='red')
  #
  axis(side = 1, at = 1:3, labels = c('a=0.25', 'a=0.6', 'a=0.8'))
  #
  legend('topright', c('simulated', 'estimated', 'mean estimated'),
         lty=c(1, NA, 1), pch=c(NA, 16, NA), col=c('black', 'darkgrey', 'red'))
  #
  #
  #
  #
  par(mar=c(3, 4, 1, 2))
  #
  plot(c(rep(0.7, 3), rep(1, 3), rep(1.3, 3), rep(1.7, 3), rep(2, 3), rep(2.3, 3), rep(2.7, 3), rep(3, 3), rep(3.3, 3)),
       c(pa[pa$s_1==0.02 & pa$a==0.25 & pa$initial_freq==f & pa$N_e==N_e, 'est_a'],
         pa[pa$s_1==0.02 & pa$a==0.6 & pa$initial_freq==f & pa$N_e==N_e, 'est_a'],
         pa[pa$s_1==0.02 & pa$a==0.8 & pa$initial_freq==f & pa$N_e==N_e, 'est_a'],
         pa[pa$s_1==0.035 & pa$a==0.25 & pa$initial_freq==f & pa$N_e==N_e, 'est_a'],
         pa[pa$s_1==0.035 & pa$a==0.6 & pa$initial_freq==f & pa$N_e==N_e, 'est_a'],
         pa[pa$s_1==0.035 & pa$a==0.8 & pa$initial_freq==f & pa$N_e==N_e, 'est_a'],
         pa[pa$s_1==0.05 & pa$a==0.25 & pa$initial_freq==f & pa$N_e==N_e, 'est_a'],
         pa[pa$s_1==0.05 & pa$a==0.6 & pa$initial_freq==f & pa$N_e==N_e, 'est_a'],
         pa[pa$s_1==0.05 & pa$a==0.8 & pa$initial_freq==f & pa$N_e==N_e, 'est_a']),
       pch=16, col='darkgrey',
       xlim = c(0.5, 3.5), ylim=c(0, 1),
       xaxt='n',
       xlab='', ylab='parameter a') # simulated parameter a
  #
  segments(c(0.6, 0.9, 1.2, 1.6, 1.9, 2.2, 2.6, 2.9, 3.2), rep(c(0.25, 0.6, 0.8), 3),
           c(0.8, 1.1, 1.4, 1.8, 2.1, 2.4, 2.8, 3.1, 3.4), rep(c(0.25, 0.6, 0.8), 3))
  #
  ys = c(mean(pa[pa$s_1==0.02 & pa$a==0.25 & pa$initial_freq==f & pa$N_e==N_e, 'est_a']),
         mean(pa[pa$s_1==0.02 & pa$a==0.6 & pa$initial_freq==f & pa$N_e==N_e, 'est_a']),
         mean(pa[pa$s_1==0.02 & pa$a==0.8 & pa$initial_freq==f & pa$N_e==N_e, 'est_a']),
         mean(pa[pa$s_1==0.035 & pa$a==0.25 & pa$initial_freq==f & pa$N_e==N_e, 'est_a']),
         mean(pa[pa$s_1==0.035 & pa$a==0.6 & pa$initial_freq==f & pa$N_e==N_e, 'est_a']),
         mean(pa[pa$s_1==0.035 & pa$a==0.8 & pa$initial_freq==f & pa$N_e==N_e, 'est_a']),
         mean(pa[pa$s_1==0.05 & pa$a==0.25 & pa$initial_freq==f & pa$N_e==N_e, 'est_a']),
         mean(pa[pa$s_1==0.05 & pa$a==0.6 & pa$initial_freq==f & pa$N_e==N_e, 'est_a']),
         mean(pa[pa$s_1==0.05 & pa$a==0.8 & pa$initial_freq==f & pa$N_e==N_e, 'est_a']))
  #
  segments(c(0.65, 0.95, 1.25, 1.65, 1.95, 2.25, 2.65, 2.95, 3.25), ys,
           c(0.75, 1.05, 1.35, 1.75, 2.05, 2.35, 2.75, 3.05, 3.35), ys, col='red')
  #
  axis(side = 1, at = 1:3, labels = c('s=0.02', 's=0.035', 's=0.05'))
  #
  mtext(text = paste0('initial freq.=', f, if (N_e==0) ', no drift' else paste0(', N_e=', N_e)),
        side = 3, line = 2, outer=TRUE)
  #
  dev.off()
}
#
#
#
for (N_e in c(0, 10000, 5000, 1000, 500)) {
  #
  pdf(paste0('results/power_analysis_s_1_a_vs_nb_der_al_', 
             if (N_e==0) 'no_drift' else paste0('Ne', N_e), '.pdf'), height=9)
  #
  par(mfrow=c(2, 1))
  #
  par(oma=c(1, 0, 4, 0))
  par(mar=c(3, 4, 1, 2))
  #
  xs = pa[pa$N_e==N_e, 'nb_der_allele']
  ys = abs(pa[pa$N_e==N_e, 'est_s_1'] - pa[pa$N_e==N_e, 's_1'])
  #
  nb_der_allele = str_count(paste0(c(genotypes_Britain, genotypes_NCentEur,
                                     genotypes_South, genotypes_Baltic), collapse=''), 'A')
  #
  plot(xs, ys, pch=16,
       xlab='', ylab='absolute estimation error in s')
  #  
  #
  lines(c(nb_der_allele, nb_der_allele), par('usr')[3:4], col='red')
  #
  legend('topright', c('simulation', 'observed nb. der. alleles'),
         pch=c(16, NA), lty=c(NA, 1), col=c('black', 'red'))
  #
  #
  #
  par(mar=c(4, 4, 0, 2))
  #
  ys = abs(pa[pa$N_e==N_e, 'est_a'] - pa[pa$N_e==N_e, 'a'])
  #
  plot(xs, ys, pch=16,
       xlab='nb. of derived LP alleles', ylab='absolute estimation error in a')
  #  
  lines(c(nb_der_allele, nb_der_allele), par('usr')[3:4], col='red')
  #
  mtext(if (N_e==0) 'no drift' else paste0('N_e = ', N_e), side = 3, line = 2, outer=TRUE)
  #
  dev.off()
}
#
#
#
pa_res = data.frame(row.names = c('N_e', 0, 10000, 5000, 1000, 500))
for (N_e in c(0, 10000, 5000, 1000, 500)) {
  #
  pdf(paste0('results/power_analysis_milk_lL_vs_nb_der_al_', 
             if (N_e==0) 'no_drift' else paste0('Ne', N_e), '.pdf'))
  #
  xs = pa[pa$N_e==N_e & pa$est_a < 0.95, 'nb_der_allele']
  ys = pa[pa$N_e==N_e & pa$est_a < 0.95, 'lL_const'] - pa[pa$N_e==N_e & pa$est_a < 0.95, 'lL_milk']
  xs_ = pa[pa$N_e==N_e & pa$est_a >= 0.95, 'nb_der_allele']
  ys_ = pa[pa$N_e==N_e & pa$est_a >= 0.95, 'lL_const'] - pa[pa$N_e==N_e & pa$est_a >= 0.95, 'lL_milk']
  nb_der_allele = str_count(paste0(c(genotypes_Britain, genotypes_NCentEur,
                                     genotypes_South, genotypes_Baltic), collapse=''), 'A')
  #
  above = sum(ys[xs>=nb_der_allele] >= (qchisq(0.05, df=1, lower.tail = FALSE) / 2)) / sum(xs>=nb_der_allele)
  total = sum(ys >= (qchisq(0.05, df=1, lower.tail = FALSE) / 2)) / length(xs)
  below = sum(ys[xs<nb_der_allele] >= (qchisq(0.05, df=1, lower.tail = FALSE) / 2)) / sum(xs<nb_der_allele)
  above_ = sum(ys_[xs_>=nb_der_allele] >= (qchisq(0.05, df=1, lower.tail = FALSE) / 2)) / sum(xs_>=nb_der_allele)
  total_ = sum(ys_ >= (qchisq(0.05, df=1, lower.tail = FALSE) / 2)) / length(xs_)
  below_ = sum(ys_[xs_<nb_der_allele] >= (qchisq(0.05, df=1, lower.tail = FALSE) / 2)) / sum(xs_<nb_der_allele)
  pa_res = rbind(pa_res, c(N_e, if (N_e>0) '0.005' else 'all', above, total, below, above_, total_, below_))
  #
  plot(xs, ys, pch=16,
       xlim = c(0, max(c(xs, xs_))), ylim = c(min(c(ys, ys_)), max(c(ys, ys_))), 
       xlab='nb. of derived LP alleles', ylab='log likelihood difference')
  #
  points(xs_, ys_, pch=16, col='darkgrey')
  #
  lines(c(nb_der_allele, nb_der_allele), par('usr')[3:4], col='blue', lty=3)
  #
  lines(par('usr')[1:2], c(qchisq(0.05, df=1, lower.tail = FALSE) / 2, qchisq(0.05, df=1, lower.tail = FALSE) / 2),
        col='blue', lty=2)
  #
  legend('top', c('simulation', 'simulation (a est. > 0.95)', 'observed nb. der. alleles', 'significance'),
         pch=c(16, 16, NA, NA), lty=c(NA, NA, 3, 2), lwd=c(NA, NA, 1, 1), col=c('black', 'darkgrey', 'blue', 'blue'))
  #
  mtext(paste0('likelihood difference milk over constant', if (N_e==0) ', no drift' else paste0(', N_e = ', N_e)),
        side = 3, line = 1)

  dev.off()
}
#
N_e = 0
f = 5e-3
xs = pa[pa$N_e==N_e & pa$est_a < 0.95 & pa$initial_freq==f, 'nb_der_allele'] # 
ys = pa[pa$N_e==N_e & pa$est_a < 0.95 & pa$initial_freq==f, 'lL_const'] - 
                                                      pa[pa$N_e==N_e & pa$est_a < 0.95 & pa$initial_freq==f, 'lL_milk']
xs_ = pa[pa$N_e==N_e & pa$est_a >= 0.95 & pa$initial_freq==f, 'nb_der_allele']
ys_ = pa[pa$N_e==N_e & pa$est_a >= 0.95 & pa$initial_freq==f, 'lL_const'] - 
                                                      pa[pa$N_e==N_e & pa$est_a >= 0.95 & pa$initial_freq==f, 'lL_milk']
nb_der_allele = str_count(paste0(c(genotypes_Britain, genotypes_NCentEur,
                                   genotypes_South, genotypes_Baltic), collapse=''), 'A')
above = sum(ys[xs>=nb_der_allele] >= (qchisq(0.05, df=1, lower.tail = FALSE) / 2)) / sum(xs>=nb_der_allele)
total = sum(ys >= (qchisq(0.05, df=1, lower.tail = FALSE) / 2)) / length(xs)
below = sum(ys[xs<nb_der_allele] >= (qchisq(0.05, df=1, lower.tail = FALSE) / 2)) / sum(xs<nb_der_allele)
above_ = sum(ys_[xs_>=nb_der_allele] >= (qchisq(0.05, df=1, lower.tail = FALSE) / 2)) / sum(xs_>=nb_der_allele)
total_ = sum(ys_ >= (qchisq(0.05, df=1, lower.tail = FALSE) / 2)) / length(xs_)
below_ = sum(ys_[xs_<nb_der_allele] >= (qchisq(0.05, df=1, lower.tail = FALSE) / 2)) / sum(xs_<nb_der_allele)
pa_res = rbind(pa_res, c(N_e, '0.005', above, total, below, above_, total_, below_))
#
colnames(pa_res) = c('N_e', 'initial freq.', 'above (a est. < 0.95)', 'total (a est. < 0.95)', 'below (a est. < 0.95)',
                     'above (a est. >= 0.95)', 'total (a est. >= 0.95)', 'below (a est. >= 0.95)')
write.csv(x = pa_res, file = 'results/power_analysis_milk_lL_vs_nb_der_al_res_table.csv')
#
#
#
# POWER ANALYSIS
#
#
#
for (a in c(0.25, 0.6, 0.8)) {
  for (s_1 in c(0.02, 0.035, 0.05)) {
    for (initial_freq_N_e in list(c(0.005, 0))) {
      N_e = initial_freq_N_e[2]
      #
      for (rep in 1:3) {
        out_infix = paste0('milk_noise_additive_f', clean_version(initial_freq_N_e[1]), '_s',
                           clean_version(s_1), '_a', clean_version(a), '_rep', rep,
                           if (N_e > 0) paste0('_Ne', clean_version(N_e)) else '')
        print(out_infix)
        #
        load(file = paste0('results/power_analysis_', out_infix, '.Robj'))
        #
        power_analysis[power_analysis$a==a & power_analysis$s_1==s_1 & power_analysis$N_e==N_e & 
                         power_analysis$initial_freq==initial_freq_N_e[1] & power_analysis$rep==rep, c(
                           'est_a_milk_10pc', 'est_s_1_milk_10pc', 'est_initial_freq_milk_10pc', 'lL_milk_milk_10pc',
                           'est_a_milk_20pc', 'est_s_1_milk_20pc', 'est_initial_freq_milk_20pc', 'lL_milk_milk_20pc',
                           'est_a_milk_30pc', 'est_s_1_milk_30pc', 'est_initial_freq_milk_30pc', 'lL_milk_milk_30pc',
                           'est_a_milk_50pc', 'est_s_1_milk_50pc', 'est_initial_freq_milk_50pc', 'lL_milk_milk_50pc')] =
          c(results_milk_noise$milk_10pc$par[3], results_milk_noise$milk_10pc$par[2],
            results_milk_noise$milk_10pc$par[1], results_milk_noise$milk_10pc$value,
            results_milk_noise$milk_20pc$par[3], results_milk_noise$milk_20pc$par[2], 
            results_milk_noise$milk_20pc$par[1], results_milk_noise$milk_20pc$value,
            results_milk_noise$milk_30pc$par[3], results_milk_noise$milk_30pc$par[2], 
            results_milk_noise$milk_30pc$par[1], results_milk_noise$milk_30pc$value,
            results_milk_noise$milk_50pc$par[3], results_milk_noise$milk_50pc$par[2],
            results_milk_noise$milk_50pc$par[1], results_milk_noise$milk_50pc$value)
      }
    }
  }
}
#
pa = power_analysis
f = 5e-03
#
noise_suffix = '_milk_50pc'
#
pdf('results/power_analysis_milk_noise_50pc_s_1_a_no_drift.pdf', height=9)
#
par(mfrow=c(3, 1))
par(oma=c(0, 0, 4, 0))
#
par(mar=c(4.5, 4, 1, 2))
#
xs = c(rep(0.7, 3), rep(1, 3), rep(1.3, 3), rep(1.7, 3), rep(2, 3), rep(2.3, 3), rep(2.7, 3), rep(3, 3), rep(3.3, 3))
ys1 = abs(pa[pa$a==0.25 & pa$s_1==0.02 & pa$initial_freq==f & pa$N_e==N_e, 'est_s_1'] - 
            pa[pa$a==0.25 & pa$s_1==0.02 & pa$initial_freq==f & pa$N_e==N_e, paste0('est_s_1', noise_suffix)])
ys2 = abs(pa[pa$a==0.25 & pa$s_1==0.035 & pa$initial_freq==f & pa$N_e==N_e, 'est_s_1'] - 
            pa[pa$a==0.25 & pa$s_1==0.035 & pa$initial_freq==f & pa$N_e==N_e, paste0('est_s_1', noise_suffix)])
ys3 = abs(pa[pa$a==0.25 & pa$s_1==0.05 & pa$initial_freq==f & pa$N_e==N_e, 'est_s_1'] -
            pa[pa$a==0.25 & pa$s_1==0.05 & pa$initial_freq==f & pa$N_e==N_e, paste0('est_s_1', noise_suffix)])
ys4 = abs(pa[pa$a==0.6 & pa$s_1==0.02 & pa$initial_freq==f & pa$N_e==N_e, 'est_s_1'] -
            pa[pa$a==0.6 & pa$s_1==0.02 & pa$initial_freq==f & pa$N_e==N_e, paste0('est_s_1', noise_suffix)])
ys5 = abs(pa[pa$a==0.6 & pa$s_1==0.035 & pa$initial_freq==f & pa$N_e==N_e, 'est_s_1'] -
            pa[pa$a==0.6 & pa$s_1==0.035 & pa$initial_freq==f & pa$N_e==N_e, paste0('est_s_1', noise_suffix)])
ys6 = abs(pa[pa$a==0.6 & pa$s_1==0.05 & pa$initial_freq==f & pa$N_e==N_e, 'est_s_1'] -
            pa[pa$a==0.6 & pa$s_1==0.05 & pa$initial_freq==f & pa$N_e==N_e, paste0('est_s_1', noise_suffix)])
ys7 = abs(pa[pa$a==0.8 & pa$s_1==0.02 & pa$initial_freq==f & pa$N_e==N_e, 'est_s_1'] -
            pa[pa$a==0.8 & pa$s_1==0.02 & pa$initial_freq==f & pa$N_e==N_e, paste0('est_s_1', noise_suffix)])
ys8 = abs(pa[pa$a==0.8 & pa$s_1==0.035 & pa$initial_freq==f & pa$N_e==N_e, 'est_s_1'] -
            pa[pa$a==0.8 & pa$s_1==0.035 & pa$initial_freq==f & pa$N_e==N_e, paste0('est_s_1', noise_suffix)])
ys9 = abs(pa[pa$a==0.8 & pa$s_1==0.05 & pa$initial_freq==f & pa$N_e==N_e, 'est_s_1'] - 
            pa[pa$a==0.8 & pa$s_1==0.05 & pa$initial_freq==f & pa$N_e==N_e, paste0('est_s_1', noise_suffix)])
plot(xs, c(ys1, ys2, ys3, ys4, ys5, ys6, ys7, ys8, ys9),
     pch=16, cex=1.2, col='darkgrey',
     xlim = c(0.5, 3.5), 
     xaxt='n', 
     xlab='', ylab='abs. estimation error in s')
#
ys = c(mean(ys1), mean(ys2), mean(ys3), mean(ys4), mean(ys5), mean(ys6), mean(ys7), mean(ys8), mean(ys9))
#
segments(c(0.65, 0.95, 1.25, 1.65, 1.95, 2.25, 2.65, 2.95, 3.25), ys,
         c(0.75, 1.05, 1.35, 1.75, 2.05, 2.35, 2.75, 3.05, 3.35), ys, col='red')
#
axis(side = 1, at = 1:3, labels = c('a=0.25', 'a=0.6', 'a=0.8'), line=1.5, tick=FALSE)
axis(side = 1, at = c(0.7, 1, 1.3, 1.7, 2, 2.3, 2.7, 3, 3.3), labels = rep(c('s=0.02', 's=0.0.035', 's=0.05'), 3),
     cex.axis=0.7)
#
legend('topright', c('estimated', 'mean estimated'),
       lty=c(NA, 1), pch=c(16, NA), col=c('darkgrey', 'red'))
#
#
#
par(mar=c(4.5, 4, 1, 2))
#
ys1 = abs(pa[pa$s_1==0.02 & pa$a==0.25 & pa$initial_freq==f & pa$N_e==N_e, 'est_a'] - 
            pa[pa$s_1==0.02 & pa$a==0.25 & pa$initial_freq==f & pa$N_e==N_e, paste0('est_a', noise_suffix)])
ys2 = abs(pa[pa$s_1==0.02 & pa$a==0.6 & pa$initial_freq==f & pa$N_e==N_e, 'est_a'] - 
            pa[pa$s_1==0.02 & pa$a==0.6 & pa$initial_freq==f & pa$N_e==N_e, paste0('est_a', noise_suffix)])
ys3 = abs(pa[pa$s_1==0.02 & pa$a==0.8 & pa$initial_freq==f & pa$N_e==N_e, 'est_a'] -
            pa[pa$s_1==0.02 & pa$a==0.8 & pa$initial_freq==f & pa$N_e==N_e, paste0('est_a', noise_suffix)])
ys4 = abs(pa[pa$s_1==0.035 & pa$a==0.25 & pa$initial_freq==f & pa$N_e==N_e, 'est_a'] -
            pa[pa$s_1==0.035 & pa$a==0.25 & pa$initial_freq==f & pa$N_e==N_e, paste0('est_a', noise_suffix)])
ys5 = abs(pa[pa$s_1==0.035 & pa$a==0.6 & pa$initial_freq==f & pa$N_e==N_e, 'est_a'] -
            pa[pa$s_1==0.035 & pa$a==0.6 & pa$initial_freq==f & pa$N_e==N_e, paste0('est_a', noise_suffix)])
ys6 = abs(pa[pa$s_1==0.035 & pa$a==0.8 & pa$initial_freq==f & pa$N_e==N_e, 'est_a'] -
            pa[pa$s_1==0.035 & pa$a==0.8 & pa$initial_freq==f & pa$N_e==N_e, paste0('est_a', noise_suffix)])
ys7 = abs(pa[pa$s_1==0.05 & pa$a==0.25 & pa$initial_freq==f & pa$N_e==N_e, 'est_a'] -
            pa[pa$s_1==0.05 & pa$a==0.25 & pa$initial_freq==f & pa$N_e==N_e, paste0('est_a', noise_suffix)])
ys8 = abs(pa[pa$s_1==0.05 & pa$a==0.6 & pa$initial_freq==f & pa$N_e==N_e, 'est_a'] -
            pa[pa$s_1==0.05 & pa$a==0.6 & pa$initial_freq==f & pa$N_e==N_e, paste0('est_a', noise_suffix)])
ys9 = abs(pa[pa$s_1==0.05 & pa$a==0.8 & pa$initial_freq==f & pa$N_e==N_e, 'est_a'] - 
            pa[pa$s_1==0.05 & pa$a==0.8 & pa$initial_freq==f & pa$N_e==N_e, paste0('est_a', noise_suffix)])
#
plot(xs, c(ys1, ys2, ys3, ys4, ys5, ys6, ys7, ys8, ys9),
     pch=16, cex=1.2, col='darkgrey',
     xlim = c(0.5, 3.5), ylim=c(0, 1),
     xaxt='n',
     xlab='', ylab='abs. estimation error in a') # simulated parameter a
#
ys = c(mean(ys1), mean(ys2), mean(ys3), mean(ys4), mean(ys5), mean(ys6), mean(ys7), mean(ys8), mean(ys9))
#
segments(c(0.65, 0.95, 1.25, 1.65, 1.95, 2.25, 2.65, 2.95, 3.25), ys,
         c(0.75, 1.05, 1.35, 1.75, 2.05, 2.35, 2.75, 3.05, 3.35), ys, col='red')
#
axis(side = 1, at = 1:3, labels = c('s=0.02', 's=0.035', 's=0.05'), line=1.5, tick=FALSE)
axis(side = 1, at = c(0.7, 1, 1.3, 1.7, 2, 2.3, 2.7, 3, 3.3), labels = rep(c('a=0.25', 'a=0.6', 'a=0.8'), 3),
     cex.axis=0.7)
#
diff_50pc = as.matrix(abs(input$milk[3:6] - input_milk_noise$milk_50pc[3:6]))
mtext(text = paste0('milk noise 50% (avg. abs. deviation=', round(mean(diff_50pc), digits = 3), ', max.=',
                    round(max(diff_50pc), digits = 2), '), initial freq.=', f, ', no drift'),
      side = 3, line = 2, outer=TRUE)
#
par(mar=c(4.5, 4, 1, 2))
#
xs = c(rep(0.7, 3), rep(1, 3), rep(1.3, 3), rep(1.7, 3), rep(2, 3), rep(2.3, 3), rep(2.7, 3), rep(3, 3), rep(3.3, 3))
#
.compute_vals <- function(dat, noise_suffix) {
  res = list('m_mn' = dat[, 'lL_milk'] - dat[, paste0('lL_milk', noise_suffix)], 
             'c_m' = dat[, 'lL_const'] - dat[, 'lL_milk'],
             'c_mn' = dat[, 'lL_const'] - dat[, paste0('lL_milk', noise_suffix)],
             'c_m_sig' = (dat[, 'lL_const'] - dat[, 'lL_milk']) >= (qchisq(0.05, df=1, lower.tail = FALSE) / 2),
             'c_mn_sig' = (dat[, 'lL_const'] - dat[, paste0('lL_milk', noise_suffix)]) >=
                                                                          (qchisq(0.05, df=1, lower.tail = FALSE) / 2))
  res[['c_m_c_mn_sig_diff']] = !xor(res$'c_m_sig', res$'c_mn_sig')#res$'c_m_sig' - res$'c_mn_sig'
  return(res)
}
ys1 = .compute_vals(pa[pa$s_1==0.02 & pa$a==0.25 & pa$initial_freq==f & pa$N_e==N_e, ], noise_suffix)
ys2 = .compute_vals(pa[pa$s_1==0.02 & pa$a==0.6 & pa$initial_freq==f & pa$N_e==N_e, ], noise_suffix)
ys3 = .compute_vals(pa[pa$s_1==0.02 & pa$a==0.8 & pa$initial_freq==f & pa$N_e==N_e, ], noise_suffix)
ys4 = .compute_vals(pa[pa$s_1==0.035 & pa$a==0.25 & pa$initial_freq==f & pa$N_e==N_e, ], noise_suffix)
ys5 = .compute_vals(pa[pa$s_1==0.035 & pa$a==0.6 & pa$initial_freq==f & pa$N_e==N_e, ], noise_suffix)
ys6 = .compute_vals(pa[pa$s_1==0.035 & pa$a==0.8 & pa$initial_freq==f & pa$N_e==N_e, ], noise_suffix)
ys7 = .compute_vals(pa[pa$s_1==0.05 & pa$a==0.25 & pa$initial_freq==f & pa$N_e==N_e, ], noise_suffix)
ys8 = .compute_vals(pa[pa$s_1==0.05 & pa$a==0.6 & pa$initial_freq==f & pa$N_e==N_e, ], noise_suffix)
ys9 = .compute_vals(pa[pa$s_1==0.05 & pa$a==0.8 & pa$initial_freq==f & pa$N_e==N_e, ], noise_suffix)
#
plot(xs, c(ys1$m_mn, ys2$m_mn, ys3$m_mn, ys4$m_mn, ys5$m_mn, ys6$m_mn, ys7$m_mn, ys8$m_mn, ys9$m_mn),
     pch=16, cex=1.2, col='darkgrey',
     xlim = c(0.5, 3.5),
     xaxt='n',
     xlab='', ylab='milk model log likelihood difference')
#
# just to make sure not to have to mark points in red
stopifnot(all(c(ys1$c_m_c_mn_sig_diff, ys2$c_m_c_mn_sig_diff, ys3$c_m_c_mn_sig_diff, ys4$c_m_c_mn_sig_diff,
                ys5$c_m_c_mn_sig_diff, ys6$c_m_c_mn_sig_diff, ys7$c_m_c_mn_sig_diff, ys8$c_m_c_mn_sig_diff, 
                ys9$c_m_c_mn_sig_diff)))
#
lines(par('usr')[1:2], c(0, 0), lwd=0.7)
#
axis(side = 1, at = 1:3, labels = c('s=0.02', 's=0.035', 's=0.05'), line=1.5, tick=FALSE)
axis(side = 1, at = c(0.7, 1, 1.3, 1.7, 2, 2.3, 2.7, 3, 3.3), labels = rep(c('a=0.25', 'a=0.6', 'a=0.8'), 3),
     cex.axis=0.7)
#
legend('bottomleft', legend = c('...unchanged', '...changed'), pch=16,
       col=c('darkgrey', 'red'), title='sig. of milk over constant model...')
#
dev.off()
#
#
#
plot_result('Settlement density', 'additive', logistic_fct_sel, res_cluster,
            time_points_Britain, genotypes_Britain, cluster_normalised[, c(1, 2, 3)],
            time_points_NCentEur, genotypes_NCentEur, cluster_normalised[, c(1, 2, 5)], 
            time_points_South, genotypes_South, cluster_normalised[, c(1, 2, 6)],
            time_points_Baltic, genotypes_Baltic, cluster_normalised[, c(1, 2, 4)])
#
plot_result('Cluster_inv', 'additive', logistic_fct_sel, res_cluster_inv,
            time_points_Britain, genotypes_Britain, cluster_normalised_inv[, c(1, 2, 3)],
            time_points_NCentEur, genotypes_NCentEur, cluster_normalised_inv[, c(1, 2, 5)], 
            time_points_South, genotypes_South, cluster_normalised_inv[, c(1, 2, 6)],
            time_points_Baltic, genotypes_Baltic, cluster_normalised_inv[, c(1, 2, 4)])
#
plot_result('Domestics', 'additive', logistic_fct_sel, res_domestics,
            time_points_Britain, genotypes_Britain, domestics_normalised[, c(1, 2, 3)],
            time_points_NCentEur, genotypes_NCentEur, domestics_normalised[, c(1, 2, 5)], 
            time_points_South, genotypes_South, domestics_normalised[, c(1, 2, 6)],
            time_points_Baltic, genotypes_Baltic, domestics_normalised[, c(1, 2, 4)])
#
plot_result('Wild animal consumption', 'additive', logistic_fct_sel, res_domestics_inv,
            time_points_Britain, genotypes_Britain, domestics_normalised_inv[, c(1, 2, 3)],
            time_points_NCentEur, genotypes_NCentEur, domestics_normalised_inv[, c(1, 2, 5)], 
            time_points_South, genotypes_South, domestics_normalised_inv[, c(1, 2, 6)],
            time_points_Baltic, genotypes_Baltic, domestics_normalised_inv[, c(1, 2, 4)])
#
plot_result('Milk', 'additive', logistic_fct_sel, res_milk,
            time_points_Britain, genotypes_Britain, milk_normalised[, c(1, 2, 3)],
            time_points_NCentEur, genotypes_NCentEur, milk_normalised[, c(1, 2, 5)], 
            time_points_South, genotypes_South, milk_normalised[, c(1, 2, 6)],
            time_points_Baltic, genotypes_Baltic, milk_normalised[, c(1, 2, 4)])
#
plot_result('Population fluctuation', 'additive', logistic_fct_sel, res_pop_fluc_globally,
            time_points_Britain, genotypes_Britain, pop_fluc_normalised[, c(1, 2, 3)],
            time_points_NCentEur, genotypes_NCentEur, pop_fluc_normalised[, c(1, 2, 5)], 
            time_points_South, genotypes_South, pop_fluc_normalised[, c(1, 2, 6)],
            time_points_Baltic, genotypes_Baltic, pop_fluc_normalised[, c(1, 2, 4)])
#
plot_result('Pop_fluc_inv', 'additive', logistic_fct_sel, res_pop_fluc_inv_globally,
            time_points_Britain, genotypes_Britain, pop_fluc_normalised_inv[, c(1, 2, 3)],
            time_points_NCentEur, genotypes_NCentEur, pop_fluc_normalised_inv[, c(1, 2, 5)], 
            time_points_South, genotypes_South, pop_fluc_normalised_inv[, c(1, 2, 6)],
            time_points_Baltic, genotypes_Baltic, pop_fluc_normalised_inv[, c(1, 2, 4)])
#
plot_result('Insolation', 'additive', logistic_fct_sel, res_insolation_globally,
            time_points_Britain, genotypes_Britain, insolation_normalised[, c(1, 2, 3)],
            time_points_NCentEur, genotypes_NCentEur, insolation_normalised[, c(1, 2, 5)], 
            time_points_South, genotypes_South, insolation_normalised[, c(1, 2, 6)],
            time_points_Baltic, genotypes_Baltic, insolation_normalised[, c(1, 2, 4)])
#
plot_result('Inverse insolation', 'additive', logistic_fct_sel, res_insolation_inv_globally,
            time_points_Britain, genotypes_Britain, insolation_normalised_inv[, c(1, 2, 3)],
            time_points_NCentEur, genotypes_NCentEur, insolation_normalised_inv[, c(1, 2, 5)], 
            time_points_South, genotypes_South, insolation_normalised_inv[, c(1, 2, 6)],
            time_points_Baltic, genotypes_Baltic, insolation_normalised_inv[, c(1, 2, 4)])
#
plot_result('Const', 'additive', logistic_fct_sel, res_const,
            time_points_Britain, genotypes_Britain, const,
            time_points_NCentEur, genotypes_NCentEur, const, 
            time_points_South, genotypes_South, const,
            time_points_Baltic, genotypes_Baltic, const)
#
#
#
plot_result('Cluster', 'dominant', dominant_sel, res_cluster_dominant,
            time_points_Britain, genotypes_Britain, cluster_normalised[, c(1, 2, 3)],
            time_points_NCentEur, genotypes_NCentEur, cluster_normalised[, c(1, 2, 5)], 
            time_points_South, genotypes_South, cluster_normalised[, c(1, 2, 6)],
            time_points_Baltic, genotypes_Baltic, cluster_normalised[, c(1, 2, 4)])
#
plot_result('Cluster_inv', 'dominant', dominant_sel, res_cluster_inv_dominant,
            time_points_Britain, genotypes_Britain, cluster_normalised_inv[, c(1, 2, 3)],
            time_points_NCentEur, genotypes_NCentEur, cluster_normalised_inv[, c(1, 2, 5)], 
            time_points_South, genotypes_South, cluster_normalised_inv[, c(1, 2, 6)],
            time_points_Baltic, genotypes_Baltic, cluster_normalised_inv[, c(1, 2, 4)])
#
plot_result('Domestics', 'dominant', dominant_sel, res_domestics_dominant,
            time_points_Britain, genotypes_Britain, domestics_normalised[, c(1, 2, 3)],
            time_points_NCentEur, genotypes_NCentEur, domestics_normalised[, c(1, 2, 5)], 
            time_points_South, genotypes_South, domestics_normalised[, c(1, 2, 6)],
            time_points_Baltic, genotypes_Baltic, domestics_normalised[, c(1, 2, 4)])
#
plot_result('Domestics_inv', 'dominant', dominant_sel, res_domestics_inv_dominant,
            time_points_Britain, genotypes_Britain, domestics_normalised_inv[, c(1, 2, 3)],
            time_points_NCentEur, genotypes_NCentEur, domestics_normalised_inv[, c(1, 2, 5)], 
            time_points_South, genotypes_South, domestics_normalised_inv[, c(1, 2, 6)],
            time_points_Baltic, genotypes_Baltic, domestics_normalised_inv[, c(1, 2, 4)])
#
plot_result('Milk', 'dominant', dominant_sel, res_milk_dominant,
            time_points_Britain, genotypes_Britain, milk_normalised[, c(1, 2, 3)],
            time_points_NCentEur, genotypes_NCentEur, milk_normalised[, c(1, 2, 5)], 
            time_points_South, genotypes_South, milk_normalised[, c(1, 2, 6)],
            time_points_Baltic, genotypes_Baltic, milk_normalised[, c(1, 2, 4)])
#
plot_result('Pop_fluc', 'dominant', dominant_sel, res_pop_fluc_dominant_globally,
            time_points_Britain, genotypes_Britain, pop_fluc_normalised[, c(1, 2, 3)],
            time_points_NCentEur, genotypes_NCentEur, pop_fluc_normalised[, c(1, 2, 5)], 
            time_points_South, genotypes_South, pop_fluc_normalised[, c(1, 2, 6)],
            time_points_Baltic, genotypes_Baltic, pop_fluc_normalised[, c(1, 2, 4)])
#
plot_result('Pop_fluc_inv', 'dominant', dominant_sel, res_pop_fluc_inv_dominant_globally,
            time_points_Britain, genotypes_Britain, pop_fluc_normalised_inv[, c(1, 2, 3)],
            time_points_NCentEur, genotypes_NCentEur, pop_fluc_normalised_inv[, c(1, 2, 5)], 
            time_points_South, genotypes_South, pop_fluc_normalised_inv[, c(1, 2, 6)],
            time_points_Baltic, genotypes_Baltic, pop_fluc_normalised_inv[, c(1, 2, 4)])
#
plot_result('Insolation', 'dominant', dominant_sel, res_insolation_dominant_globally,
            time_points_Britain, genotypes_Britain, insolation_normalised[, c(1, 2, 3)],
            time_points_NCentEur, genotypes_NCentEur, insolation_normalised[, c(1, 2, 5)], 
            time_points_South, genotypes_South, insolation_normalised[, c(1, 2, 6)],
            time_points_Baltic, genotypes_Baltic, insolation_normalised[, c(1, 2, 4)])
#
plot_result('Insolation_inv', 'dominant', dominant_sel, res_insolation_inv_dominant_globally,
            time_points_Britain, genotypes_Britain, insolation_normalised_inv[, c(1, 2, 3)],
            time_points_NCentEur, genotypes_NCentEur, insolation_normalised_inv[, c(1, 2, 5)], 
            time_points_South, genotypes_South, insolation_normalised_inv[, c(1, 2, 6)],
            time_points_Baltic, genotypes_Baltic, insolation_normalised_inv[, c(1, 2, 4)])
#
plot_result('Const', 'dominant', dominant_sel, res_const_dominant,
            time_points_Britain, genotypes_Britain, const,
            time_points_NCentEur, genotypes_NCentEur, const, 
            time_points_South, genotypes_South, const,
            time_points_Baltic, genotypes_Baltic, const)
#
#
#
pdf(file = 'results/log_likelihoods.pdf', width=5.9, height=5)
par(mar=c(9, 4, 3, 5))
#
ys_additive = -c(res_const$value,
                 #
                 res_milk$value,
                 #
                 res_insolation_globally$value, res_insolation_inv_globally$value,
                 #
                 res_pop_fluc_globally$value, res_pop_fluc_inv_globally$value,
                 #
                 res_cluster$value, res_cluster_inv$value,
                 #
                 res_domestics$value, res_domestics_inv$value)
#
ys_dominant = -c(res_const_dominant$value,
                 #
                 res_milk_dominant$value,
                 #
                 res_insolation_dominant_globally$value, res_insolation_inv_dominant_globally$value,
                 #
                 res_pop_fluc_dominant_globally$value, res_pop_fluc_inv_dominant_globally$value,
                 #
                 res_cluster_dominant$value, res_cluster_inv_dominant$value, 
                 #
                 res_domestics_dominant$value, res_domestics_inv_dominant$value)
#
#
plot(ys_dominant,
     pch='-', col='blue', cex=3,
     ylim = range(c(ys_additive, ys_dominant)), xaxt='n', xlab='', ylab='log likelihood')
#
points(ys_additive, pch='-', col='black', cex=3)
#
lines(c(1.5, 1.5), par('usr')[3:4], col='grey', lwd=0.7)
lines(c(2.5, 2.5), par('usr')[3:4], col='grey', lwd=0.7)
text(labels = 'Calcium assi.', x = 3.1, y = -131.6, cex = 0.8, pos = 2, srt = 90, col = 'grey')
lines(c(4.5, 4.5), par('usr')[3:4], col='grey', lwd=0.7)
text(labels = 'crisis', x = 6.45, y = -131.6, cex = 0.8, pos = 2, srt = 90, col = 'grey')
lines(c(6.5, 6.5), par('usr')[3:4], col='grey', lwd=0.7)
text(labels = 'chronic', x = 8.1, y = -131.6, cex = 0.8, pos = 2, srt = 90, col = 'grey')
lines(c(8.5, 8.5), par('usr')[3:4], col='grey', lwd=0.7)
#
lines(par('usr')[1:2], c(-res_const$value + qchisq(0.05, df=1, lower.tail = FALSE) / 2,
                         -res_const$value + qchisq(0.05, df=1, lower.tail = FALSE) / 2), col='red')
lines(par('usr')[1:2], c(-res_const_dominant$value + qchisq(0.05, df=1, lower.tail = FALSE) / 2,
                         -res_const_dominant$value + qchisq(0.05, df=1, lower.tail = FALSE) / 2),
      lty = 'dashed', col='red')
#
axis(side = 1, at = 1:10, labels = c('constant',
                                     #
                                     'milk',
                                     #
                                     'insolation', 'insolation inv.',
                                     #
                                     'population fluc.', 'population fluc. inv.',
                                     #
                                     'cluster', 'cluster inv.', 
                                     #
                                     'domesticates', 'wild'),
     las=2, cex=0.8)
#
par(xpd=TRUE)
legend('topright', c('additive', 'dominant', 'significance (additive)', 'significance (dominant)'), 
       pch=c('-', '-', NA, NA), pt.cex = c(3, 3, 1, 1), lty=c(NA, NA, 1, 2),
       col=c('black', 'blue', 'red', 'red'), inset = c(-0.23, -0.18), bg = 'white', cex = 0.8)
#
dev.off()
#
#
#
res_table = data.frame(model=c('constant, additive',
                               #
                               'milk, additive',
                               #
                               'insolation, additive', 'insolation inv., additive',
                               #
                               'population fluc., additive', 'population fluc. inv., additive',
                               #
                               'cluster, additive', 'cluster inv., additive', 
                               #
                               'domesticates, additive', 'wild, additive',
                               #
                               'constant, dominant',
                               #
                               'milk, dominant',
                               #
                               'insolation, dominant', 'insolation inv., dominant',
                               #
                               'population fluc., dominant', 'population fluc. inv., dominant',
                               #
                               'cluster, dominant', 'cluster inv., dominant', 
                               #
                               'domesticates, dominant', 'wild, dominant'),
                       nb_parameters=c(2, rep(3, 9), 2, rep(3, 9)),
                       log_likelihood=c(ys_additive, ys_dominant),
                       LRT_p_value=c(NA,
                         sapply(ys_additive, function(x) pchisq(2*(x + res_const$value), df=1, lower.tail = FALSE))[-1],
                             NA,
               sapply(ys_dominant, function(x) pchisq(2*(x + res_const_dominant$value), df=1, lower.tail = FALSE))[-1]))
#
write.csv(x = res_table, file = 'results/res_table.csv')
#
#
#
res_table
#
mt = res_table[2:10, c(1, 4)]
mt = mt[order(mt[2]), ]
mt = cbind(mt, ((1:nrow(mt)) / nrow(mt)) * 0.02)
mt = cbind(mt, mt[2] <= mt[3])
#
write.csv(x = mt, file = 'results/res_table_multiple_testing_additive.csv')
#
mt = res_table[12:20, c(1, 4)]
mt = mt[order(mt[2]), ]
mt = cbind(mt, ((1:nrow(mt)) / nrow(mt)) * 0.02)
mt = cbind(mt, mt[2] <= mt[3])
#
write.csv(x = mt, file = 'results/res_table_multiple_testing_dominant.csv')
#
#
#