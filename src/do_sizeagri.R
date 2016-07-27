if (!exists("prj")) {
  stop("You need to create a object 'prj' that points to the top folder, 
       e.g. prj <- '/home/edisz/Documents/Uni/Projects/PHD/4BFG/Paper/ms_est' or
       prj <- '/home/user/Documents/projects_git/ms_est'!")
} else {
  source(file.path(prj, "src", "load.R"))
}

### ----------------------------------------------------------------------------
### Code for Modelling influence of catchment size and agriculture

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
plot(mod_nb, pages = 1)
plot(mod_nb, pages = 1, residuals = TRUE) 
gam.check(mod_nb)
# overdispersion
r <- resid(mod_nb, type = "pearson")
sum(r^2) / (mod_nb$df.res)

#! TODO: Model visualisation
#! derivation function?

rm(samples_raks, samples_exceed, raks, rak_exceed, nd1)



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
maxlogtumax <- logtumax[ , list(maxlogtumax = max(logtumax)), by = site_id]
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


#! same pattern, how should this be modelled?
#! zinf gaussian? tobin?
#! definitively use a two stage model!
#! where should be censored?
