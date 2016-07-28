if (!exists("prj")) {
  stop("You need to create a object 'prj' that points to the top folder, 
       e.g. prj <- '/home/edisz/Documents/Uni/Projects/PHD/4BFG/Paper/ms_est' or
       prj <- '/home/user/Documents/projects_git/ms_est'!")
} else {
  source(file.path(prj, "src", "load.R"))
}

### ----------------------------------------------------------------------------
### Code for Modelling influence of catchment size and agriculture

library(devtools)
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


# model using gam 
rak_exceed$agri_fin <- rak_exceed$agri_fin*100
library(mgcv)
mod_p <- gam(n_exceed ~ s(agri_fin) + s(ezg_fin) + offset(logn), data = rak_exceed, 
           family = poisson)
plot(mod_p, pages = 1)
plot(mod_p, pages = 1, residuals = TRUE) 
gam.check(mod_p)

# overdispersion
r <- resid(mod_p, type = "pearson")
sum(r^2) / (mod_p$df.res)
#! present

# try negative binomial model
# with offset, automatic theta search and REML
mod_nb <- gam(n_exceed ~ s(agri_fin) + s(ezg_fin) + offset(logn), data = rak_exceed,
           family = nb())
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
mod_nb.dsig_ezg <- signifD(pdat$value[pdat$variable == 'ezg_fin'], d = mod_nb.d[['ezg_fin']]$deriv,
                            mod_nb.dci_ezg[['ezg_fin']]$upper, mod_nb.dci_ezg[['ezg_fin']]$lower)

# check if significant
pdat$sig[pdat$variable == 'agri_fin'] <- ifelse(!is.na(mod_nb.dsig_agri$incr) | !is.na(mod_nb.dsig_agri$decr), TRUE, FALSE)
pdat$sig[pdat$variable == 'ezg_fin'] <- ifelse(!is.na(mod_nb.dsig_ezg$incr) | !is.na(mod_nb.dsig_ezg$decr), TRUE, FALSE)
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
  
ggsave(file.path(prj, "/fig/figrac.svg"),
       p, width = 8, height = 5)







# Model max(TUmax) --------------------------------------------------------
psm_maxtu # = as calculated from database

# calculate R (for verification)
setkey(var_props, variable_id)
setkey(psm_samples, variable_id)

samples_lc50 <- var_props[ , list(variable_id, lc50_dm_fin)][psm_samples]
# rm samples without lc50 data
samples_lc50 <- samples_lc50[!is.na(lc50_dm_fin)]

# calculate log(tu)
# min value is -8.72
range(samples_lc50[value_fin > 0 , log10(value_fin / lc50_dm_fin)])
# set zero to -8.75
samples_lc50[ , logtu := ifelse(value_fin > 0, log10(value_fin / lc50_dm_fin), -8.75)]

# calculate per sample max(logtu)
logtumax <- samples_lc50[ , list(logtumax = max(logtu)), by = sample_id]
# join back site & date
setkey(logtumax, sample_id)
setkey(samples_lc50, sample_id)
logtumax <- logtumax[unique(samples_lc50[ , list(sample_id, site_id, date)])]

# join site_info with logtumax
setkey(logtumax, site_id)
setkey(psm_sites, site_id)
logtumax_s <- psm_sites_wi[logtumax]
logtumax_s <- logtumax_s[!(is.na(agri_fin) | is.na(ezg_fin)), ]
logtumax_s

ggplot(logtumax_s, aes(x = ezg_fin, y = logtumax)) +
  geom_point() +
  geom_smooth()

ggplot(logtumax_s, aes(x = agri_fin, y = logtumax)) +
  geom_point() +
  geom_smooth()


# to avoid problems with independency I aggregate using the max per site
# otherwise need a gamm with site as random effect?!
# take maximum
maxlogtumax <- logtumax[ , list(maxlogtumax = max(logtumax)), by = site_id]
# taker 95% quantile of all samples per sites
maxlogtumax <- logtumax[ , list(maxlogtumax = quantile(logtumax, 0.95)), by = site_id]

setkey(maxlogtumax, site_id)
setkey(psm_sites, site_id)
maxlogtumax_s <- psm_sites_wi[maxlogtumax]
maxlogtumax_s <- maxlogtumax_s[!(is.na(agri_fin) | is.na(ezg_fin)), ]

ggplot(maxlogtumax_s, aes(x = ezg_fin, y = maxlogtumax)) +
  geom_point() +
  geom_smooth()

ggplot(maxlogtumax_s, aes(x = agri_fin, y = maxlogtumax)) +
  geom_point(alpha = 0.5) +
  geom_smooth()


#! use 95% percentile instead of max?

#! same pattern, how should this be modelled?
#! use censored gaussion

library(gamlss)
library(gamlss.add)
library(survival)
library(gamlss.cens)
