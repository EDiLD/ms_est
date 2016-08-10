if (!exists("prj")) {
  stop("You need to create a object 'prj' that points to the top folder, 
       e.g. prj <- '/home/edisz/Documents/Uni/Projects/PHD/4BFG/Paper/ms_est' or
       prj <- '/home/user/Documents/projects_git/ms_est'!")
} else {
  source(file.path(prj, "src", "load.R"))
}

### ----------------------------------------------------------------------------
### Code for Modelling influence of catchment size and agriculture

source("https://gist.githubusercontent.com/EDiLD/c4b748eb57b71c027cbfe9714222a918/raw/4cf44afb401e715548e3b2fedb4cb7d30976c8df/Deriv.R")

# Load data ---------------------------------------------------------------
psm_sites <- fread(file.path(cachedir, 'psm_sites.csv'))
psm_samples <- fread(file.path(cachedir, 'psm_samples.csv'))
psm_sites_info <- fread(file.path(cachedir, 'psm_sites_info.csv'))
psm_variables <- fread(file.path(cachedir, 'psm_variables.csv'))
psm_maxtu <- fread(file.path(cachedir, 'psm_maxtu.csv'))
var_props <- fread(file.path(cachedir, 'var_props.csv'))


nrow(psm_sites_info[!(is.na(ezg_fin) | is.na(agri_fin))])
# 2376 with both data
# 265 completely missing
# 408 with either ezg or agri info missing

# join sites with info
setkey(psm_sites, site_id)
setkey(psm_sites_info, site_id)
psm_sites_wi <- psm_sites[ , list(site_id, state, easting, northing)][psm_sites_info[ , list(site_id, use, ezg_fin, agri_fin)]]



# Model number of RAK exceedances per sites -------------------------------

# select compounds with RAKS
raks <- psm_variables[!is.na(rak_uba), list(variable_id, name, cas, pgroup, rak_uba)]

# filter samples with RAKS
samples_raks <- psm_samples[variable_id %in% raks$variable]
# join samples & RAKS
setkey(samples_raks, variable_id)
setkey(raks, variable_id)
samples_raks <- raks[samples_raks]
# check if RAK is exceeded
samples_raks[ , rak_exceed := value_fin > rak_uba]

# calculate number of exceedances per sites
table(samples_raks$rak_exceed, useNA = 'always')
samples_exceed <- samples_raks[ , list(n_exceed = sum(rak_exceed), # = irrespective of compound
                   p_exceed = sum(rak_exceed) / length(rak_exceed),
                   n = length(rak_exceed)), by = site_id]


# join site_info with exceedances
setkey(samples_exceed, site_id)
setkey(psm_sites, site_id)
rak_exceed <- psm_sites_wi[samples_exceed]
rak_exceed
# 2970 sites
rm(samples_raks, samples_exceed)
rak_exceed <- rak_exceed[!(is.na(agri_fin) | is.na(ezg_fin)), ]
# 2343 sites with ezg and agri data


ggplot(rak_exceed, aes(x = ezg_fin, y = p_exceed)) +
  geom_point() +
  geom_smooth()

ggplot(rak_exceed, aes(x = agri_fin, y = p_exceed)) +
  geom_point() +
  geom_smooth()

ggplot(rak_exceed, aes(x = ezg_fin, y = n_exceed)) +
  geom_point() +
  geom_smooth()

ggplot(rak_exceed, aes(x = agri_fin, y = n_exceed)) +
  geom_point() +
  geom_smooth()

rak_exceed$logn <- log(rak_exceed$n)

p_raw <- ggplot(rak_exceed, aes(y = agri_fin, x = ezg_fin, col = log(n_exceed))) +
  geom_point(alpha = 0.5) +
  mytheme +
  scale_color_gradient(low = 'blue', high = 'red') +
  ylab("Agriculture [%]") +
  xlab("Catchment Size [km2]") +
  ggtitle('RAC Exceedances')
# p_raw
ggsave(file.path(prj, "supplement/ezgagrirac.pdf"),
       p_raw)

# model using gam 
rak_exceed$agri_fin <- rak_exceed$agri_fin*100
mod_p <- gam(n_exceed ~ s(agri_fin, bs = 'cr') + s(ezg_fin, bs = 'cr') + offset(logn), data = rak_exceed, 
           family = poisson, method = 'REML')
plot(mod_p, pages = 1)
plot(mod_p, pages = 1, residuals = TRUE) 
gam.check(mod_p)

# overdispersion
r <- resid(mod_p, type = "pearson")
sum(r^2) / (mod_p$df.res)
#! slightly present

# try negative binomial model
# with offset, automatic theta search and REML
mod_nb <- gam(n_exceed ~ s(agri_fin, bs = 'cr') + s(ezg_fin, bs = 'cr') + offset(logn), 
              data = rak_exceed,
              family = nb(),
              method = 'REML')

mod_nb_te <- gam(n_exceed ~ te(agri_fin, ezg_fin, bs = 'cr') + offset(logn), 
              data = rak_exceed,
              family = nb(),
              method = 'REML')
plot(mod_nb_te)
vis.gam(mod_nb_te, view = c('agri_fin', 'ezg_fin'))
anova(mod_nb, mod_nb_te, test = 'Chisq') 
# smoothing interaction not significant and can be omited



# offset out of formula: =ignored in predict?
#! Check
# mod_nb <- gam(n_exceed ~ s(agri_fin) + s(ezg_fin), offset = rak_exceed$logn, data = rak_exceed, 
#               family = nb())
plot(mod_nb, pages = 1)
plot(mod_nb, pages = 1, residuals = TRUE) 
gam.check(mod_nb)
# overdispersion
r <- resid(mod_nb, type = "pearson")
sum(r^2) / (mod_nb$df.res)
# OK


# calculate predictions for agri & ezg
# fix other variable and n at mean values
pdat_agri <- with(rak_exceed,
             data.frame(agri_fin = c(seq(min(agri_fin), max(agri_fin), length.out = 100)),
                        ezg_fin = rep(mean(ezg_fin), 100),
                        logn = rep(mean(logn), 100)))
pred_agri <- predict(mod_nb, newdata = pdat_agri, type = 'response', se.fit = TRUE)
pdat_agri <- transform(pdat_agri, 
                       fit_agri = pred_agri$fit)
pdat_agri <- transform(pdat_agri, 
                       up_agri = fit_agri + (1.96 * pred_agri$se.fit),
                       low_agri = fit_agri - (1.96 * pred_agri$se.fit))
pdat_agri$logn <- NULL
pdat_agri$ezg_fin <- NULL

pdat_ezg <- with(rak_exceed,
                  data.frame(ezg_fin = c(seq(min(ezg_fin), max(ezg_fin), length.out = 100)),
                             agri_fin = rep(mean(agri_fin), 100),
                             logn = rep(mean(logn), 100)))
pred_ezg <- predict(mod_nb, newdata = pdat_ezg, type = 'response', se.fit = TRUE)
pdat_ezg <- transform(pdat_ezg, 
                       fit_ezg = pred_ezg$fit)
pdat_ezg <- transform(pdat_ezg, 
                       up_ezg = fit_ezg + (1.96 * pred_ezg$se.fit),
                       low_ezg = fit_ezg - (1.96 * pred_ezg$se.fit))
# prepare data.frame for plotting
pdat_ezg$logn <- NULL
pdat_ezg$agri_fin <- NULL

pdat_agri <- melt(pdat_agri, measure.vars = 'agri_fin')
names(pdat_agri) <- c('fit', 'up', 'low', 'variable', 'value')
pdat_agri$variable <- as.character(pdat_agri$variable)

pdat_ezg <- melt(pdat_ezg, measure.vars = 'ezg_fin')
names(pdat_ezg) <- c('fit', 'up', 'low', 'variable', 'value')
pdat_ezg$variable <- as.character(pdat_ezg$variable)
pdat <- rbind(pdat_agri, pdat_ezg)

# calculate derivatives
# see http://www.fromthebottomoftheheap.net/2014/05/15/identifying-periods-of-change-with-gams/
# use the same pdat
mod_nb.d <- Deriv(mod_nb, newdata = data.frame(agri_fin = pdat$value[pdat$variable == 'agri_fin'],
                                               ezg_fin = pdat$value[pdat$variable == 'ezg_fin'],
                                               logn = mean(rak_exceed$logn)))

mod_nb.dci_agri <- confint(mod_nb.d, term = 'agri_fin')
mod_nb.dsig_agri <- signifD(pdat$value[pdat$variable == 'agri_fin'], 
                            d = mod_nb.d[['agri_fin']]$deriv,
                            mod_nb.dci_agri[['agri_fin']]$upper, 
                            mod_nb.dci_agri[['agri_fin']]$lower)

mod_nb.dci_ezg <- confint(mod_nb.d, term = 'ezg_fin')
mod_nb.dsig_ezg <- signifD(pdat$value[pdat$variable == 'ezg_fin'], 
                           d = mod_nb.d[['ezg_fin']]$deriv,
                            mod_nb.dci_ezg[['ezg_fin']]$upper, 
                           mod_nb.dci_ezg[['ezg_fin']]$lower)

# check if significant
pdat$sig[pdat$variable == 'agri_fin'] <- ifelse(!is.na(mod_nb.dsig_agri$incr) | 
                                                  !is.na(mod_nb.dsig_agri$decr), 
                                                TRUE, FALSE)
pdat$sig[pdat$variable == 'ezg_fin'] <- ifelse(!is.na(mod_nb.dsig_ezg$incr) | 
                                                 !is.na(mod_nb.dsig_ezg$decr), 
                                               TRUE, FALSE)
pdat <- transform(pdat,
                  sig_value = ifelse(sig, fit, NA))

mylabeller <- as_labeller(c(
  'agri_fin'="Agriculture [%]",
  'ezg_fin'="Catchment Size [km2]"
))


# plot
p <- ggplot(pdat, aes(x = value, y = fit, group = variable)) +
  geom_line() +
  geom_line(aes(y = up), lty ='dashed') +
  geom_line(aes(y = low), lty ='dashed') +
  geom_line(aes(y = sig_value), colour = 'red', lwd = 1.5) +
  facet_wrap(~variable, scales = 'free_x', labeller = mylabeller) +
  mytheme +
  xlab('Value') +
  ylab('No. RAC exceedances') +
  ylim(c(0, 1.3))
# p
ggsave(file.path(prj, "/fig/figrac.svg"),
       p, width = 8, height = 5)







# Model log(TUmax) --------------------------------------------------------
psm_maxtu # = view as calculated from database

# calculate withing R (for verification)
setkey(var_props, variable_id)
setkey(psm_samples, variable_id)

samples_lc50 <- var_props[ , list(variable_id, lc50_dm_fin)][psm_samples]
# rm samples without lc50 data
samples_lc50 <- samples_lc50[!is.na(lc50_dm_fin)]

# calculate log(tu)
# min value is -8.72
range(samples_lc50[value_fin > 0 , log10(value_fin / lc50_dm_fin)])

# logtu: set zero to -8.75
samples_lc50[ , logtu := ifelse(value_fin > 0, log10(value_fin / lc50_dm_fin), -8.75)]
# tu (use with log-normal dist)
samples_lc50[ , tu := value_fin / lc50_dm_fin]

# calculate per sample max(logtu)
logtumax <- samples_lc50[ , list(logtumax = max(logtu)), by = sample_id]


# join with data from database, to check.
setkey(logtumax, sample_id)
setkey(psm_maxtu, sample_id)
jj <- psm_maxtu[logtumax]
jj[ , diff := log_maxtu - logtumax]
range(jj$diff)
# Nice, perfect. Same output as with postgresql view


tumax <- samples_lc50[ , list(tumax = max(tu)), by = sample_id]


# join back site & date
setkey(tumax, sample_id)
setkey(samples_lc50, sample_id)
tumax <- tumax[unique(samples_lc50[ , list(sample_id, site_id, date)])]

# join site_info with logtumax
setkey(tumax, site_id)
setkey(psm_sites, site_id)
tumax_s <- psm_sites_wi[tumax]
tumax_s <- tumax_s[!(is.na(agri_fin) | is.na(ezg_fin)), ]
tumax_s

ggplot(tumax_s, aes(x = ezg_fin, y = tumax)) +
  geom_point() +
  geom_smooth() +
  scale_y_log10()

ggplot(tumax_s, aes(x = agri_fin, y = tumax)) +
  geom_point() +
  geom_smooth() +
  scale_y_log10()




# to avoid problems with independency we aggregate per site
# otherwise need a gamm with site as random effect?!

# # take maximum
# maxtumax <- tumax[ , list(maxtumax = max(tumax)), by = site_id]

# taker 95% quantile of all samples per sites
# take log
maxtumax <- tumax[ , list(log_95tumax = log10(quantile(tumax, 0.95))), by = site_id]


setkey(maxtumax, site_id)
setkey(psm_sites, site_id)
maxtumax_s <- psm_sites_wi[maxtumax]
maxtumax_s <- maxtumax_s[!(is.na(agri_fin) | is.na(ezg_fin)), ]

ggplot(maxtumax_s, aes(x = ezg_fin, y = log_95tumax)) +
  geom_point() +
  geom_smooth() 

ggplot(maxtumax_s, aes(x = agri_fin, y = log_95tumax)) +
  geom_point(alpha = 0.5) +
  geom_smooth() 



# use tobit gam (with normal errors) to model the effects

# censoring indicator (censored at -8.75)
maxtumax_s$log_95tumax <- ifelse(maxtumax_s$log_95tumax == -Inf, -8.75, maxtumax_s$log_95tumax)
maxtumax_s$is_cens <- ifelse(maxtumax_s$log_95tumax <= -8.75, TRUE, FALSE)

# log-gaussian Tobit GAM
library(gamlss.add)
dfa <- maxtumax_s[ , list(ezg_fin, agri_fin, log_95tumax, is_cens)]
setDF(dfa)
mod_logno_cens <- 
  gamlss(Surv(log_95tumax, is_cens == FALSE, type = 'left') ~
           pb(agri_fin) + pb(ezg_fin), 
                         data = dfa,
                         family = cens(NO, type = 'left'))
plot(mod_logno_cens)

pdf(file.path(prj, "supplement/ezgagritumodel.pdf"), width = 10)
term.plot(mod_logno_cens, rug = TRUE, pages = 1, scheme = 'lines',
          col.term = "black", 
          col.se = "black",
          main = c('Agriculture [%]', 'Catchment Size [km2]'),
          xlabs = '',
          ylabs = 'Effekt',
          ask = FALSE)
dev.off()


p_raw <- ggplot(dfa, aes(y = agri_fin, x = ezg_fin, col = log_95tumax)) +
  geom_point(alpha = 0.5) +
  mytheme +
  scale_color_gradient(low = 'blue', high = 'red') +
  ylab("Agriculture [%]") +
  xlab("Catchment Size [km2]") +
  ggtitle('95% percentile of log(TU)')
# p_raw
ggsave(file.path(prj, "supplement/ezgagritu.pdf"),
       p_raw)









# Model number of EQS exceedances per sites -------------------------------

# select compounds with EQS
eqs <- psm_variables[!is.na(wrrl_zhkuqn), list(variable_id, name, cas, pgroup, wrrl_zhkuqn)]

# filter samples with EQS
samples_eqs <- psm_samples[variable_id %in% eqs$variable]
# join samples & EQS
setkey(samples_eqs, variable_id)
setkey(eqs, variable_id)
samples_eqs <- eqs[samples_eqs]
# check if EQS is exceeded
samples_eqs[ , eqs_exceed := value_fin > wrrl_zhkuqn]

# calculate number of exceedances per sites
table(samples_eqs$eqs_exceed, useNA = 'always')
samples_exceed <- samples_eqs[ , list(n_exceed = sum(eqs_exceed), # = irrespective of compound
                                       p_exceed = sum(eqs_exceed) / length(eqs_exceed),
                                       n = length(eqs_exceed)), by = site_id]


# join site_info with exceedances
setkey(samples_exceed, site_id)
setkey(psm_sites, site_id)
eqs_exceed <- psm_sites_wi[samples_exceed]
eqs_exceed
# 2996 sites
rm(samples_eqs, samples_exceed)
eqs_exceed <- eqs_exceed[!(is.na(agri_fin) | is.na(ezg_fin)), ]
# 2344 sites with ezg and agri data
# offset
eqs_exceed$logn <- log(eqs_exceed$n)
eqs_exceed$agri_fin <- eqs_exceed$agri_fin*100


ggplot(eqs_exceed, aes(x = ezg_fin, y = p_exceed)) +
  geom_point() +
  geom_smooth()

ggplot(eqs_exceed, aes(x = agri_fin, y = p_exceed)) +
  geom_point() +
  geom_smooth()

ggplot(eqs_exceed, aes(x = ezg_fin, y = n_exceed)) +
  geom_point() +
  geom_smooth()

ggplot(eqs_exceed, aes(x = agri_fin, y = n_exceed)) +
  geom_point() +
  geom_smooth()



p_raw <- ggplot(eqs_exceed, aes(y = agri_fin, x = ezg_fin, col = log(n_exceed))) +
  geom_point(alpha = 0.5) +
  mytheme +
  scale_color_gradient(low = 'blue', high = 'red') +
  ylab("Agriculture [%]") +
  xlab("Catchment Size [km2]") +
  ggtitle('EQS Exceedances')
# p_raw
ggsave(file.path(prj, "supplement/ezgagrieqs.pdf"),
       p_raw)

# model using gam
mod_p <- gam(n_exceed ~ s(agri_fin, bs = 'cr') + s(ezg_fin, bs = 'cr') + offset(logn), 
             data = eqs_exceed,
             family = poisson, method = 'REML')
gam.check(mod_p)
#! increase k for ezg_fin
mod_p <- gam(n_exceed ~ s(agri_fin, bs = 'cr') + s(ezg_fin, bs = 'cr', k = 20) + offset(logn), 
             data = eqs_exceed,
             family = poisson, method = 'REML')
gam.check(mod_p) 
#! still problematic (also with 200)
plot(mod_p, pages = 1)

# overdispersion
r <- resid(mod_p, type = "pearson")
sum(r^2) / (mod_p$df.res)
#! slightly present
# 
# try negative binomial model
# with offset, automatic theta search and REML
mod_nb <- gam(n_exceed ~ s(agri_fin, bs = 'cr') + s(ezg_fin, bs = 'cr') + offset(logn),
              data = eqs_exceed,
              family = nb(),
              method = 'REML')
gam.check(mod_nb)
#! still problems with basis dimension?
#! maybe because of either: excess of zeros or ezg_gradient too long
plot(mod_nb, pages = 1)
# overdispersion
r <- resid(mod_nb, type = "pearson")
sum(r^2) / (mod_nb$df.res)
#! better

mod_nb_te <- gam(n_exceed ~ te(agri_fin, ezg_fin) + offset(logn),
                 data = eqs_exceed,
                 family = nb(),
                 method = 'REML')
plot(mod_nb_te)
vis.gam(mod_nb_te, view = c('agri_fin', 'ezg_fin'))
anova(mod_nb, mod_nb_te, test = 'Chisq')
# smoothing interaction not significant and can be omited


# calculate predictions for agri & ezg
# fix other variable and n at mean values
pdat_agri <- with(eqs_exceed,
                  data.frame(agri_fin = c(seq(min(agri_fin), max(agri_fin), length.out = 100)),
                             ezg_fin = rep(mean(ezg_fin), 100),
                             logn = rep(mean(logn), 100)))
pred_agri <- predict(mod_nb, newdata = pdat_agri, type = 'response', se.fit = TRUE)
pdat_agri <- transform(pdat_agri,
                       fit_agri = pred_agri$fit)
pdat_agri <- transform(pdat_agri,
                       up_agri = fit_agri + (1.96 * pred_agri$se.fit),
                       low_agri = fit_agri - (1.96 * pred_agri$se.fit))
pdat_agri$logn <- NULL
pdat_agri$ezg_fin <- NULL

pdat_ezg <- with(eqs_exceed,
                 data.frame(ezg_fin = c(seq(min(ezg_fin), max(ezg_fin), length.out = 100)),
                            agri_fin = rep(mean(agri_fin), 100),
                            logn = rep(mean(logn), 100)))
pred_ezg <- predict(mod_nb, newdata = pdat_ezg, type = 'response', se.fit = TRUE)
pdat_ezg <- transform(pdat_ezg,
                      fit_ezg = pred_ezg$fit)
pdat_ezg <- transform(pdat_ezg,
                      up_ezg = fit_ezg + (1.96 * pred_ezg$se.fit),
                      low_ezg = fit_ezg - (1.96 * pred_ezg$se.fit))
# prepare data.frame for plotting
pdat_ezg$logn <- NULL
pdat_ezg$agri_fin <- NULL

pdat_agri <- melt(pdat_agri, measure.vars = 'agri_fin')
names(pdat_agri) <- c('fit', 'up', 'low', 'variable', 'value')
pdat_agri$variable <- as.character(pdat_agri$variable)

pdat_ezg <- melt(pdat_ezg, measure.vars = 'ezg_fin')
names(pdat_ezg) <- c('fit', 'up', 'low', 'variable', 'value')
pdat_ezg$variable <- as.character(pdat_ezg$variable)
pdat <- rbind(pdat_agri, pdat_ezg)

# calculate derivatives
# see http://www.fromthebottomoftheheap.net/2014/05/15/identifying-periods-of-change-with-gams/
# use the same pdat
mod_nb.d <- Deriv(mod_nb, newdata = data.frame(agri_fin = pdat$value[pdat$variable == 'agri_fin'],
                                               ezg_fin = pdat$value[pdat$variable == 'ezg_fin'],
                                               logn = mean(eqs_exceed$logn)))

mod_nb.dci_agri <- confint(mod_nb.d, term = 'agri_fin')
mod_nb.dsig_agri <- signifD(pdat$value[pdat$variable == 'agri_fin'],
                            d = mod_nb.d[['agri_fin']]$deriv,
                            mod_nb.dci_agri[['agri_fin']]$upper,
                            mod_nb.dci_agri[['agri_fin']]$lower)

mod_nb.dci_ezg <- confint(mod_nb.d, term = 'ezg_fin')
mod_nb.dsig_ezg <- signifD(pdat$value[pdat$variable == 'ezg_fin'],
                           d = mod_nb.d[['ezg_fin']]$deriv,
                           mod_nb.dci_ezg[['ezg_fin']]$upper,
                           mod_nb.dci_ezg[['ezg_fin']]$lower)

# check if significant
pdat$sig[pdat$variable == 'agri_fin'] <- ifelse(!is.na(mod_nb.dsig_agri$incr) |
                                                  !is.na(mod_nb.dsig_agri$decr),
                                                TRUE, FALSE)
pdat$sig[pdat$variable == 'ezg_fin'] <- ifelse(!is.na(mod_nb.dsig_ezg$incr) |
                                                 !is.na(mod_nb.dsig_ezg$decr),
                                               TRUE, FALSE)
pdat <- transform(pdat,
                  sig_value = ifelse(sig, fit, NA))

mylabeller <- as_labeller(c(
  'agri_fin'="Agriculture [%]",
  'ezg_fin'="Catchment Size [km2]"
))


# plot
p <- ggplot(pdat, aes(x = value, y = fit, group = variable)) +
  geom_line() +
  geom_line(aes(y = up), lty ='dashed') +
  geom_line(aes(y = low), lty ='dashed') +
  geom_line(aes(y = sig_value), colour = 'red', lwd = 1.5) +
  facet_wrap(~variable, scales = 'free_x', labeller = mylabeller) +
  mytheme +
  xlab('Value') +
  ylab('No. EQS exceedances') 
# p
ggsave(file.path(prj, "/supplement/ezgagrieqsmodel.pdf"),
       p, width = 8, height = 5)
