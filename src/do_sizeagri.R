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
psm_sites_wi <- psm_sites[ , list(site_id, state, easting, 
                                  northing)][psm_sites_info[ , list(site_id, 
                                                                    use, ezg_fin, agri_fin)]]



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
# exclude NA (=sites wihtout ezg & agri info)
rak_exceed <- rak_exceed[!(is.na(agri_fin) | is.na(ezg_fin)), ]
# 2343 sites with ezg and agri data

# exclude sites > 150km2
rak_exceed <- rak_exceed[ezg_fin < 150]


# logn used as offset
rak_exceed$logn <- log(rak_exceed$n)


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
mod_p <- gam(n_exceed ~ s(agri_fin, bs = 'cr') + s(ezg_fin, bs = 'cr') + offset(logn), 
             data = rak_exceed, 
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

# model with interaction
mod_nb_ti <- gam(n_exceed ~ s(agri_fin, bs = 'cr') + s(ezg_fin, bs = 'cr') +  
                   ti(agri_fin, ezg_fin, bs = 'cr') + offset(logn), 
              data = rak_exceed,
              family = nb(),
              method = 'REML')
plot(mod_nb_ti)
vis.gam(mod_nb_ti, view = c('agri_fin', 'ezg_fin'))
anova(mod_nb, mod_nb_ti, test = 'Chisq') 
# smoothing interaction not of interest


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
odat <- melt(rak_exceed[ , list(ezg_fin, agri_fin)])
p <- ggplot() +
  geom_line(data = pdat, aes(x = value, y = fit, group = variable)) +
  geom_line(data = pdat, aes(x = value, y = up, group = variable), lty ='dashed') +
  geom_line(data = pdat, aes(x = value, y = low, group = variable), lty ='dashed') +
  geom_line(data = pdat, aes(x = value, y = sig_value, group = variable), 
            colour = 'red', lwd = 1.5) +
  facet_wrap(~variable, scales = 'free_x', labeller = mylabeller) +
  geom_rug(data = odat, aes(x = value), alpha = 0.05) +
  mytheme +
  xlab('') +
  ylab('No. RAC exceedances') +
  ylim(c(0, 1.3))
# p
ggsave(file.path(prj, "figure4.pdf"),
       p, width = 7, height = 7/1.6,
       units = 'in', dpi = 300, scale = 1)


