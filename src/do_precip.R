if (!exists("prj")) {
  stop("You need to create a object 'prj' that points to the top folder, 
       e.g. prj <- '/home/edisz/Documents/Uni/Projects/PHD/4BFG/Paper/ms_est' or
       prj <- '/home/user/Documents/projects_git/ms_est'!")
} else {
  source(file.path(prj, "src", "load.R"))
}

### ----------------------------------------------------------------------------
### Code to explore the effect of precipitation

# Load data ---------------------------------------------------------------
psm_sites <- fread(file.path(cachedir, 'psm_sites.csv'))
psm_samples <- fread(file.path(cachedir, 'psm_samples.csv'))
psm_sites_info <- fread(file.path(cachedir, 'psm_sites_info.csv'))
psm_variables <- fread(file.path(cachedir, 'psm_variables.csv'))

precip_dates <- readRDS(file.path(cachedir, 'precip_dates.rds'))
precip_dates1 <- readRDS(file.path(cachedir, 'precip_dates1.rds'))



# Risk Quotients ----------------------------------------------------------
# restrict samples to variables with rac
rac <- psm_variables[!is.na(rak_uba), list(variable_id, name, cas, pgroup, rak_uba)]
# join eqs
setkey(psm_samples, variable_id)
setkey(rac, variable_id)
samples_rac <- psm_samples[rac]
# rm left over from left join
samples_rac <- samples_rac[!is.na(site_id)]

# calculate Risk Quotient
samples_rac[ , rq := value_fin / rak_uba]
# rm unused cols
samples_rac <- samples_rac[ , list(sample_id, site_id, date, variable_id, rq)]


# join precipitation
setkey(precip_dates, sample_id)
setkey(samples_rac, sample_id)
samples_rac <- precip_dates[ , list(sample_id, val)][samples_rac]
setnames(samples_rac, 'val', 'precip0')

setkey(precip_dates1, sample_id)
setkey(samples_rac, sample_id)
samples_rac <- precip_dates1[ , list(sample_id, val)][samples_rac]
setnames(samples_rac, 'val', 'precip_1')


rm(precip_dates1, precip_dates, psm_samples)


# join catchment size
take_si <- psm_sites_info[!is.na(ezg_fin) | is.na(agri_fin) , list(site_id, agri_fin, ezg_fin)]
setkey(take_si, 'site_id')
setkey(samples_rac, 'site_id')
take <- take_si[samples_rac][!is.na(agri_fin)]



# rm obs without precip information
take <- take[!(is.na(precip_1) | is.na(precip0))]
# use log10(x + 0.05 transformation)
take[, log_precip_1 := log10(precip_1 + 0.05)]
take[, log_precip0 := log10(precip0 + 0.05)]
# make factors
take[ , v_id_fac := factor(variable_id)]
take[ , s_id_fac := factor(site_id)]
take[ , state_fac := factor(substr(site_id, 1, 2))]

# day1 and day0 are only sligthly correlated
with(take, cor(log_precip_1, log_precip0))
#! use them separetely in modelling


# split data into two groups (below and above 10mm)
take$precip_1_g <- factor(ifelse(take$precip_1 <= 15, 'low', 'high'), levels = c('low', 'high'))
take$precip0_g <- factor(ifelse(take$precip0 <= 15, 'low', 'high'), levels = c('low', 'high'))

# set season / quarters
take[ , season := cut(month(date), breaks = c(0.5, 3.5, 6.5, 9.5, 12.5), 
                      labels = c('Q1', 'Q2', 'Q3', 'Q4'))]

# add group
setkey(take, variable_id)
setkey(psm_variables, variable_id)
take <- psm_variables[ , list(variable_id, psm_type)][take]


# remove compounds with no risk
props <- take[  , list(prop = sum(rq>0) / length(rq),
              abs = sum(rq>0),
              tot = length(rq)), by = variable_id][order(prop, decreasing = TRUE)]
hist(props[prop < 0.2 , prop])
# keep compounds that were found in at least 5% of samples
# and with at least 500 observations
hist(props[tot < 1000, tot])

keep <- props[prop > 0.05 & tot > 1000]
# 25 compounds left

# export overview for supplement
keep
setkey(psm_variables, variable_id)
setkey(keep, variable_id)
keep_tab <- psm_variables[ , list(variable_id, name, cas, psm_type)][keep]
keep_tab[ , prop := round(prop*100, 2)]
setnames(keep_tab, c('id', 'Compound', 'CAS', 'Group', '\\%>LOQ', 'no. > LOQ', 'total no.'))
keep_tab$id <- NULL
keep_tab

keep_tab_x <- xtable(keep_tab, 
                    label = 'tab:var_model',
                    caption = '24 pesticides for which we modelled the relationship with precipitation and seasonality.
                    Order is the same as in Figure 5 of the articles. See Table \\ref{tab:var_model_coef} for model coefficients.',
                    align = 'lp{2.5cm}rlp{1.5cm}p{2cm}p{2cm}')

print(keep_tab_x, 
      file = file.path(prj, 'supplement/keeptab.tex'),
      tabular.environment="longtable",
      floating = FALSE,
      caption.placement = 'top',
      comment = FALSE,
      booktabs = TRUE,
      hline.after = c(-1, 0, nrow(keep_tab)),
      sanitize.text.function = identity,
      size="\\fontsize{8pt}{10pt}\\selectfont"
)



take <- take[variable_id %in% keep$variable_id]
psm_variables[variable_id %in% keep$variable_id, list(variable_id, name, psm_type)]




# Model -------------------------------------------------------------------
take_c <- take[variable_id == 457]

# gamma hurdle model with precipitation as groups, season and psm_type as predictors
# site within state as random effect
mod_l_m <- gamlss(rq ~ log_precip_1 + log_precip0 + season +
                    re(random=~1|state_fac/s_id_fac),
                  # model also pi with same predictors 
                  nu.formula =~precip_1_g + precip0_g + season + 
                    re(random=~1|state_fac/s_id_fac),
                  # sigma is constant
                  data = take_c,
                  family = ZAGA)
plot(mod_l_m)
summary(mod_l_m)
term.plot(mod_l_m)
res <- residuals(mod_l_m)
hist(res)


model_foo <- function(var){
  message('Running model on compound: ', var)
  take_c <<- take[variable_id == var]
  mod <- gamlss(rq ~ log_precip_1 + log_precip0 + season +
                      re(random=~1|state_fac/s_id_fac),
                    # model also pi with same predictors 
                    nu.formula =~log_precip_1 + log_precip0 + season + 
                      re(random=~1|state_fac/s_id_fac),
                    # sigma is constant
                    data = take_c,
                    family = ZAGA)
  # save model to cache
  # since would be to big to hold in RAM
  saveRDS(mod, file.path(cachedir, 'lmodels', paste0('mod_', var, '.rds')))
}

# model_foo(727)
# run model compounds
lapply(keep$variable_id, model_foo)


                  # function to extract the needed model components, 
model_extr <- function(file){
  mod <- readRDS(file)
  ss <- summary(mod)
  smod <- data.frame(ss)
  # remove NA coefficients
  mu.na <- which(is.na(mod$mu.coefficients))
  nu.na <- which(is.na(mod$nu.coefficients))
  
  # extract coefs
  smod$terms <- c(paste0('mu.', names(mod$mu.coefficients)[-mu.na]),
                'sigma',
                paste0('nu.', names(mod$nu.coefficients)[-nu.na]))
  smod$var <- gsub('mod_(.*)\\.rds', '\\1', basename(file))
  # CI
  mult <- qnorm(0.975)
  smod$Estimate[grepl('nu\\.', smod$terms)] <- -smod$Estimate[grepl('nu\\.', smod$terms)]
  smod$upci <- smod$Estimate + smod$Std..Error * mult
  smod$lowci <- smod$Estimate - smod$Std..Error * mult
  # calculate CIs
  return(smod)
}

res <- lapply(file.path(cachedir, 'lmodels', paste0('mod_', keep$variable_id, '.rds')), model_extr)


resdf <- rbindlist(res)
names(resdf) <- c('est', 'stderr', 'tval', 'pval', 'term', 'variable', 'upci', 'lowci')
resdf$term_type <- gsub('^(.*)\\..*', '\\1', resdf$term)
resdf$term2 <- gsub('^.*\\.(.*)$', '\\1', resdf$term)
resdf$pval <- round(resdf$pval, 4)
resdf$term_type[resdf$term_type == 'nu'] <- 'pi'

# indicatopr for significant based on p
resdf$estsig <- ifelse(resdf$pval < 0.05, resdf$est, NA)
# indicatopr for significant based on ci
resdf$cisig <- ifelse(sign(resdf$upci)  == sign(resdf$lowci), 'cisig', 'nocisig')

resdf$termind <- substr(resdf$term, 1, 2)

# export other table
keep_tab2 <- resdf[!grepl('Intercept', term) & term_type != 'sigma', 
                   list(variable, term2, term_type, est, lowci, upci)]
keep_tab2[ , variable := as.numeric(variable)]
setkey(keep_tab2, variable)
setkey(psm_variables, variable_id)
keep_tab2 <- psm_variables[ , list(variable_id, name)][keep_tab2]
keep_tab2 <- keep_tab2[order(term2, term_type)]
keep_tab2[ , est := round(est, 2)]
keep_tab2[ , lowci := round(lowci, 2)]
keep_tab2[ , upci := round(upci, 2)]
make_bold <- function(string) {
  paste0("\\textbf{", string, "}")
}
keep_tab2[ , est := as.character(est)]
keep_tab2[ sign(upci)  == sign(lowci), est := make_bold(est)]
keep_tab2[ , estci := paste0(est, '\\newline (', lowci, ' - ', upci, ')')]
keep_tab2[ , term2_p := mapvalues(term2, c("precip0", "precip_1", "seasonQ2", 
                                           "seasonQ3", "seasonQ4"),
                       c("$precip_0$", "$precip_{-1}$", "$season_{Q2}$", 
                         "$season_{Q3}$", "$season_{Q4}$"))]

keep_tab2[ , term_type_p := mapvalues(term_type, c('mu', 'pi'),
                                     c('$\\mu$', '$\\pi$'))]
setnames(keep_tab2, c("variable_id", "Compound", "term2", "term_type", "est", 
                      "lowci", "upci", "coefficient", 
                      "variable", "effect"))

keep_tab2_w <- keep_tab2[ , list(Compound, variable, effect, coefficient)]
keep_tab2_w <- dcast(keep_tab2_w, Compound + effect ~ variable)
keep_tab2_w <- keep_tab2_w[order(keep_tab2_w$effect), ]
rownames(keep_tab2_w)<-NULL

keep_tab2_x <- xtable(keep_tab2_w, 
                     label = 'tab:var_model_coef',
                     caption = 'Raw data for figure 5 in the main article. 
                     Bold values denote coefficients where the CI encompasses zero.',
                     align = 'lp{2cm}p{0.7cm}p{2cm}p{2cm}p{2cm}p{2cm}p{2cm}')

print(keep_tab2_x, 
      file = file.path(prj, 'supplement/keeptab2.tex'),
      tabular.environment="longtable",
      floating = FALSE,
      caption.placement = 'top',
      comment = FALSE,
      booktabs = TRUE,
      hline.after = c(-1, 0, 24, 48),
      sanitize.text.function = identity,
      size="\\fontsize{8pt}{10pt}\\selectfont"
)




# # display as tile plot
# ggplot() +
#   geom_tile(data = resdf[!grepl('Intercept|si|mu.pre|nu', resdf$term, )], aes(x = term, y = factor(variable), fill = estsig)) +
#   scale_fill_gradient2(na.value = "white") +
#   ggtitle('mu.season')
# 
# ggplot() +
#   geom_tile(data = resdf[!grepl('Intercept|si|Q|nu', resdf$term, )], aes(x = term, y = factor(variable), fill = estsig)) +
#   ggtitle('mu.precip') +
#   scale_fill_gradient2(na.value = "white") 
# 
# 
# ggplot() +
#   geom_tile(data = resdf[!grepl('Intercept|si|nu.pre|mu', resdf$term, )], aes(x = term, y = factor(variable), fill = -estsig)) +
#   scale_fill_gradient2(na.value = "white") +
#   ggtitle('nu.season')
# 
# ggplot() +
#   geom_tile(data = resdf[!grepl('Intercept|si|Q|mu', resdf$term, )], aes(x = term, y = factor(variable), fill = -estsig)) +
#   scale_fill_gradient2(na.value = "white") +
#   ggtitle('nu.precip')



# display estimates and errors
pdata <- resdf[!term2 %chin% c('(Intercept)', 'sigma') & term2 %in% c('log_precip0', 'log_precip_1')]
pdata[ , variable := factor(variable, levels = rev(sort(unique(pdata[ , variable]))))]
p_precip <- ggplot(data = pdata) +
  geom_pointrange(aes(x = term2, y = est, ymax = upci, ymin = lowci, fill = variable,
                      col = cisig), 
                  position=position_dodge(width = .6)) +
  geom_hline(aes(yintercept = 0), linetype = 'dotted') +
  facet_wrap(~term_type, labeller = label_parsed) +
  coord_flip() +
  mytheme +
  theme(legend.position="none") +
  scale_color_manual(values = c('black', 'grey70')) +
  labs(x = '', y = 'Coefficient') +
  scale_x_discrete(breaks = c('log_precip_1', 'log_precip0'), 
                   labels = c(expression('log'~prec[-1]), expression('log'~prec[0]))) +
  theme(strip.text.x = element_text(size = 22))

pdata2 <- resdf[!term2 %chin% c('(Intercept)', 'sigma') & !term2 %in% c('log_precip0', 'log_precip_1')]
pdata2[ , variable := factor(variable, levels = rev(sort(unique(pdata2[ , variable]))))]
p_season <- ggplot(data = pdata2) +
  geom_pointrange(aes(x = term2, y = est, ymax = upci, ymin = lowci, fill = variable,
                      col = cisig), 
                  position=position_dodge(width = .6)) +
  coord_flip() +
  facet_wrap(~term_type,labeller = label_parsed) +
  geom_hline(aes(yintercept = 0), linetype = 'dotted') +
  mytheme +
  theme(legend.position="none") +
  scale_color_manual(values = c('black', 'grey70')) +
  labs(x = '', y = 'Coefficient') +
  scale_x_discrete(breaks = c('seasonQ4', 'seasonQ3', 'seasonQ2'), 
                   labels = c(expression(seas[Q4]), expression(seas[Q3]),expression(seas[Q2]))) +
  theme(strip.text.x = element_text(size = 22))

p <- arrangeGrob(p_precip, p_season, ncol = 1)
# plot(p)
ggsave(file.path(prj, "figure5.pdf"), p, width = 10, height = 9)








# extract random effect variances

# function to extract the needed model components, 
ra_extr <- function(file){
  message('Working on file: ', file)
  mod <- readRDS(file)
  lme_out <- mod$mu.coefSmo[[1]]
  vc <- VarCorr(lme_out)
  suppressWarnings(storage.mode(vc) <- 'numeric')
  vc <- vc[c(2, 4), 'StdDev']
  names(vc) <- c('state/site', 'state')
  vc <- c(vc, var = gsub('mod_(.*)\\.rds', '\\1', basename(file)))
   # generalized rsq
  vc <- c(vc, rsq = Rsq(mod, type = "Cox Snell"))
  # calculate CIs
  return(vc)
}

ra_res <- lapply(file.path(cachedir, 'models', paste0('mod_', keep$variable_id, '.rds')), ra_extr)
ra_res <- do.call(rbind, ra_res)
ra_res <- apply(ra_res, 2, as.numeric)
# ramnge of rsq values
range(ra_res[,4])


# # try by compound models
# take_c <- take[variable_id == 349]
# 
# (tc <- take[ , list(ng0 = sum(rq > 0), ng0p = sum(rq > 0) / length(rq)), by = variable_id][order(ng0p,decreasing = TRUE )])
# tc[1:100]
# # # take_cc <- take[variable_id %in% tc[ng0p > 0.1, variable_id]]
# # take_cc <- take[variable_id %in% tc[1:5, variable_id]]
# 
# # start with mixed model (site nested within state) and linear effects
# # use only precipitation, as random effect takes most of ezg and agri?
# mod_l_m_g <- gamlss(rq ~ precip_1_g + precip0_g + 
#                       re(random=~1|state_fac/s_id_fac),
#                     # model also pi with same predictors 
#                     nu.formula =~precip_1 + precip0 +
#                       re(random=~1|state_fac/s_id_fac),
#                     # sigma is constant
#                     data = take_c,
#                     family = ZAGA)
# summary(mod_l_m_g)
# term.plot(mod_l_m_g)
# mod_l_m_g$mu.coefSmo
# res <- residuals(mod_l_m_g)
# pacf(res)
# # Ok can life with that...
# plot(coef(mod_l_m_g))




######
### Give up (for now) to model this as regression....

# # try ZAGA (=Gamma Hurdle model?)
# # to handle zeros correctly?
# # order according to date (for temporal autocorrelation)
# take2 <- take[order(take$date)]
# take2[ , state_fac := factor(substr(sample_id, 1, 2))]
# rm(take)
# 
# head(take2)
# ggplot(take2, aes(x = log_precip0, y = rq)) +
#   geom_point(alpha = 0.1) +
#   scale_y_log10()
# 
# # start with mixed model (site nested within state) and linear effectes
# mod_l_m <- gamlss(rq ~ log_precip_1 + log_precip0 + re(random=~1|state_fac/s_id_fac),
#                # model also pi 
#                nu.formula =~log_precip_1 + log_precip0 + re(random=~1|state_fac/s_id_fac),
#                # sigma is constant
#                data = take2,
#                family = ZAGA)
# res_l_m <- residuals(mod_l_m) # dunn-smyth residuals
# hist(res_l_m)    # OK
# acf(res_l_m)     #! temporal autocorrelation 
# pacf(res_l_m)
# term.plot(mod_l_m, what = 'mu') 
# term.plot(mod_l_m, what = 'nu') 
# summary(mod_l_m)
# mod_l_m$mu.coefSmo
# 
# # with simple random effect
# mod_l_ms <- gamlss(rq ~ log_precip_1 + log_precip0 + re(random=~1|s_id_fac),
#                   # model also pi 
#                   nu.formula =~log_precip_1 + log_precip0 + re(random=~1|state_fac),
#                   # sigma is constant
#                   data = take2,
#                   family = ZAGA)
# extractAIC(mod_l_m)
# extractAIC(mod_l_ms)
# #! nested random effect is much better
# 
# 
# # model with smoother from gamlss (but also linear (df = 0))
# # just to check if re and smoother work together
# mod_s_m <- gamlss(rq ~ cs(log_precip_1, df = 0) + cs(log_precip0, df = 0) + 
#                     re(random=~1|state_fac/s_id_fac),
#                 # model also pi 
#                 nu.formula =~cs(log_precip_1, df = 0) + cs(log_precip0, df = 0) + 
#                   re(random=~1|state_fac/s_id_fac),
#                 # sigma is constant
#                 data = take2,
#                 family = ZAGA)
# plot(fitted(mod_l_m), fitted(mod_s_m))
# # gives identical fits
# mod_l_m$mu.coefSmo[[1]]
# mod_s_m$mu.coefSmo[[3]]
# # identical random effects
# #! => works
# res_s_m <- residuals(mod_s_m) # dunn-smyth residuals
# hist(res_s_m)    # OK
# acf(res_s_m)     #! temporal autocorrelation 
# pacf(res_s_m)    #! temporal autocorrelation 
# 
# # model with smoother from mgcv (penalized cubic splince estimated via REML)
# gac <- ga.control(method = 'REML')
# mod_sr_m <- gamlss(rq ~ ga(~s(log_precip_1, bs = 'cr') + s(log_precip0, bs = 'cr'), control = gac) + 
#                          re(random=~1|state_fac/s_id_fac),
#                        # model also pi 
#                        nu.formula =~ga(~s(log_precip_1, bs = 'cr') + s(log_precip0, bs = 'cr'), control = gac) +
#                          re(random=~1|state_fac/s_id_fac),
#                        # sigma is constant
#                        data = take2,
#                        family = ZAGA)
# mod_sr_m$mu.coefSmo[[2]]
# #! phi is estimated and different from cubic smooth
# res_sr_m <- residuals(mod_sr_m) # dunn-smyth residuals
# hist(res_sr_m)    # OK
# acf(res_sr_m)     #! temporal autocorrelation 
# pacf(res_sr_m)    #! temporal autocorrelation 
# pacf(res_s_m) 
# term.plot(mod_sr_m, what = 'mu', pages = 1, ask = FALSE) # oversmooth?
# # CI seems no appropriate...
# term.plot(mod_sr_m, what = 'nu', pages = 1, ask = FALSE) # nearly linear
# 
# 
# 
# 
# # # model with smoother from gamlss (but also linear (df = 0))
# # # and ar1 process
# # mod_s_m_ar1 <- gamlss(rq ~ cs(log_precip_1, df = 0) + cs(log_precip0, df = 0) + 
# #                     re(random=~1|state_fac/s_id_fac, correlation = corARMA(form = ~ 1|state_fac/s_id_fac, p = 1)),
# #                   # model also pi 
# #                   nu.formula =~cs(log_precip_1, df = 0) + cs(log_precip0, df = 0) + 
# #                     re(random=~1|state_fac/s_id_fac, correlation = corARMA(form = ~ 1|state_fac/s_id_fac, p = 1)),
# #                   # sigma is constant
# #                   data = take2,
# #                   family = ZAGA)
# # mod_s_m$mu.coefSmo[[3]]
# # mod_s_m_ar1$mu.coefSmo[[3]]
# # #! phi is estimated! (and small)
# # #! works
# # res_s_m_ar1 <- residuals(mod_s_m_ar1) # dunn-smyth residuals
# # hist(res_s_m_ar1)    # OK
# # acf(res_s_m_ar1)     #! temporal autocorrelation 
# # pacf(res_s_m_ar1)    #! temporal autocorrelation 
# 
# 
# 
# # model with smoother from mgcv (penalized cubic splince estimated via REML)
# # and ar1 process
# gac <- ga.control(method = 'REML')
# mod_sr_m_ar1 <- gamlss(rq ~ ga(~s(log_precip_1, bs = 'cr') + s(log_precip0, bs = 'cr'), control = gac) + 
#                         re(random=~1|state_fac/s_id_fac, correlation = corAR1(form = ~ 1|state_fac/s_id_fac)),
#                       # model also pi 
#                       nu.formula =~ga(~s(log_precip_1, bs = 'cr') + s(log_precip0, bs = 'cr'), control = gac) +
#                         re(random=~1|state_fac/s_id_fac, correlation = corAR1(form = ~ 1|state_fac/s_id_fac)),
#                       # sigma is constant
#                       data = take2,
#                       family = ZAGA)
# #! take very long to run
# mod_s_m_ar1$mu.coefSmo[[3]]
# mod_sr_m_ar1$mu.coefSmo[[2]]
# #! phi is estimated and different from cubic smooth
# res_sr_m_ar1 <- residuals(mod_sr_m_ar1) # dunn-smyth residuals
# hist(res_sr_m_ar1)    # OK
# acf(res_sr_m_ar1)     #! temporal autocorrelation 
# pacf(res_sr_m_ar1)    #! temporal autocorrelation 
# pacf(res_s_m_ar1) 
# term.plot(mod_sr_m_ar1, what = 'mu', pages = 1, ask = FALSE) # oversmooth?
# # CI seems no appropriate...
# term.plot(mod_sr_m_ar1, what = 'nu', pages = 1, ask = FALSE) # nearly linear
# 
# 
# 
# # Try CAR1 (=continuos version of AR1)
# #! take forevever...
# gac <- ga.control(method = 'REML')
# # make numeric date
# take2$date_num <- as.numeric(factor(take2$date))
# mod_sr_m_car <- gamlss(rq ~ ga(~s(log_precip_1, bs = 'cr') + s(log_precip0, bs = 'cr'), control = gac) + 
#                          re(random=~1|state_fac/s_id_fac, correlation = corCAR1(form = ~ date_num|state_fac/s_id_fac)),
#                        # model also pi 
#                        nu.formula =~ga(~s(log_precip_1, bs = 'cr') + s(log_precip0, bs = 'cr'), control = gac) +
#                          re(random=~1|state_fac/s_id_fac, correlation = corCAR1(form = ~ date_num|state_fac/s_id_fac)),
#                        # sigma is constant
#                        data = take2,
#                        family = ZAGA)
# #! takes very long!!!
# mod_sr_m_car$mu.coefSmo[[2]]
# #! phi is estimated and different from cubic smooth
# res_sr_m_car <- residuals(mod_sr_m_car) # dunn-smyth residuals
# hist(res_sr_m_car)    # OK
# acf(res_sr_m_car)     #! temporal autocorrelation 
# acf(res_sr_m_ar2)     #! temporal autocorrelation slightly better
# pacf(res_sr_m_car)   
# pacf(res_sr_m_ar2) 
# term.plot(mod_sr_m_car, what = 'mu', pages = 1, ask = FALSE) # oversmooth? why this break?!
# term.plot(mod_sr_m_car, what = 'nu', pages = 1, ask = FALSE) 
























# # model with smoother from mgcv (penalized cubic splince estimated via REML)
# # and ar2 process
# gac <- ga.control(method = 'REML')
# mod_sr_m_ar2 <- gamlss(rq ~ ga(~s(log_precip_1, bs = 'cr') + s(log_precip0, bs = 'cr'), control = gac) + 
#                          re(random=~1|s_id_fac, correlation = corARMA(form = ~ 1|s_id_fac, p = 2)),
#                        # model also pi 
#                        nu.formula =~ga(~s(log_precip_1, bs = 'cr') + s(log_precip0, bs = 'cr'), control = gac) +
#                          re(random=~1|s_id_fac, correlation = corARMA(form = ~ 1|s_id_fac, p = 2)),
#                        # sigma is constant
#                        data = take2,
#                        family = ZAGA)
# mod_sr_m_ar1$mu.coefSmo[[2]]
# mod_sr_m_ar2$mu.coefSmo[[2]]
# #! phi is estimated and different from cubic smooth
# res_sr_m_ar2 <- residuals(mod_sr_m_ar2) # dunn-smyth residuals
# hist(res_sr_m_ar2)    # OK
# acf(res_sr_m_ar2)     #! temporal autocorrelation 
# acf(res_sr_m_ar1)     #! temporal autocorrelation slightly better
# pacf(res_sr_m_ar2)   
# pacf(res_sr_m_ar1) 
# term.plot(mod_sr_m_ar2, what = 'mu', pages = 1, ask = FALSE) # oversmooth? why this break?!
# term.plot(mod_sr_m_ar2, what = 'nu', pages = 1, ask = FALSE) 
# 
# 
# # LRT (correct?)
# a <- extractAIC(mod_sr_m_ar1)
# b <- extractAIC(mod_sr_m_ar2)
# 1-pchisq(b[2]-a[2], b[1]-a[1])
# 
# 
# # Try AR(3) process?!
# #! take forevever...
# gac <- ga.control(method = 'REML')
# mod_sr_m_ar3 <- gamlss(rq ~ ga(~s(log_precip_1, bs = 'cr') + s(log_precip0, bs = 'cr'), control = gac) + 
#                          re(random=~1|s_id_fac, correlation = corARMA(form = ~ 1|s_id_fac, p = 3)),
#                        # model also pi 
#                        nu.formula =~ga(~s(log_precip_1, bs = 'cr') + s(log_precip0, bs = 'cr'), control = gac) +
#                          re(random=~1|s_id_fac, correlation = corARMA(form = ~ 1|s_id_fac, p = 3)),
#                        # sigma is constant
#                        data = take2,
#                        family = ZAGA)
# #! takes forevever... abborted
# 
# 
# # Try CAR1 (=continuos version of AR1)
# #! take forevever...
# gac <- ga.control(method = 'REML')
# # make numeric date
# take2$date_num <- as.numeric(factor(take2$date))
# mod_sr_m_car <- gamlss(rq ~ ga(~s(log_precip_1, bs = 'cr') + s(log_precip0, bs = 'cr'), control = gac) + 
#                          re(random=~1|s_id_fac, correlation = corCAR1(form = ~ date_num|s_id_fac)),
#                        # model also pi 
#                        nu.formula =~ga(~s(log_precip_1, bs = 'cr') + s(log_precip0, bs = 'cr'), control = gac) +
#                          re(random=~1|s_id_fac, correlation = corCAR1(form = ~ date_num|s_id_fac)),
#                        # sigma is constant
#                        data = take2,
#                        family = ZAGA)
# #! take very long!!!
# 
# mod_sr_m_car$mu.coefSmo[[2]]
# #! phi is estimated and different from cubic smooth
# res_sr_m_car <- residuals(mod_sr_m_car) # dunn-smyth residuals
# hist(res_sr_m_car)    # OK
# acf(res_sr_m_car)     #! temporal autocorrelation 
# acf(res_sr_m_ar2)     #! temporal autocorrelation slightly better
# pacf(res_sr_m_car)   
# pacf(res_sr_m_ar2) 
# term.plot(mod_sr_m_car, what = 'mu', pages = 1, ask = FALSE) # oversmooth? why this break?!
# term.plot(mod_sr_m_car, what = 'nu', pages = 1, ask = FALSE) 















plot(rq ~ log_precip_1, data = take)
preds <- predictAll(mod2gc, type = 'response')
# observed vs fitted
plot(take$rq, preds$mu)
abline(0, 1)
plot(take$rq, preds$mu * (1-preds$nu))
abline(0, 1)


plot(take$log_precip_1, preds$mu * (1-preds$nu))
# same as
plot(take$log_precip_1, meanZAGA(mod2gc))


sort(table(take$site_id))
# plot only for one site
ts <- which(take$site_id == 'SL_189')
plot(take$log_precip_1[ts], meanZAGA(mod2gc)[ts])


# devel section.... -------------------------------------------------------


# psm_variables[grepl('Chlorp', name)]
ggplot(samples_rac[variable_id == 378], aes(x = precip_1, y = rq)) + 
  geom_point(alpha = 0.2) +
  scale_x_log10() +
  scale_y_log10() +
  geom_hline(aes(yintercept = 0.0015))


sub <- samples_rac[variable_id == 378]
# rm obs without precip information
sub <- sub[!(is.na(precip_1) | is.na(precip0))]
# assume all values belwo rq < 0.0015 as censored
hist(sub$rq[sub$rq < 0.01])
cens <- 0.0015
sub[, is_cens := ifelse(rq < cens, TRUE, FALSE)]
# set censored values to 1/2 censoring values 
sub[is_cens == TRUE, rq := cens / 2]
# log-transform (this is eqivalnet to a log10(x + 0.0075) transformation)
sub[, log_rq := log10(rq)]
# log-transform precipitation
hist(sub$precip_1[sub$precip_1 < 0.5])
# use log10(x + 0.05 transformation)
sub[, log_precip_1 := log10(precip_1 + 0.05)]

# devide into two groups (@10mm)
sub[ , precip_1_group := ifelse(precip_1 < 10, 'low', 'high')]
plot(log_rq ~ factor(precip_1_group), data = sub)



ggplot(sub, aes(x = log_precip_1, y = log_rq, col = is_cens)) + 
  geom_point(alpha = 0.2) 


# use mixed model
mod <- lmer(log_rq ~ precip_1_group + (1|site_id), data = sub)
summary(mod)
plot(mod)
plot(ranef(mod))
#! doesn't look nice
mod0<- lmer(log_rq ~ 1 + (1|site_id), data = sub)
anova(mod, mod0)

# gamma (cens set at 1/2 loq)
mod <- glmer(rq ~ precip_1 + (date|site_id), family = Gamma(link = "log"), data = sub)
summary(mod)
plot(mod)
plot(ranef(mod))
#! doesn't look nice

# normal with temporal ar1
require(nlme)
mod <- lme(log_rq ~ precip_1_group, random=~1|site_id, data = sub)
summary(mod)
plot(mod)
plot(ranef(mod))
plot(ACF(mod), alpha = 0.05)

modar <- lme(log_rq ~ precip_1_group, random=~1|site_id, data = sub, correlation = corAR1())
summary(modar)
plot(modar)
plot(ranef(modar))
plot(ACF(modar, resType="pearson"), alpha = 0.05)
# much better job, still sime negativ autocorrelation?




data(Pixel)
l1 <- lme(pixel~ day+I(day^2), data=Pixel, random=random=~1|Side,
          method="ML")
summary(l1)
t1<-gamlss(pixel~re(fixed=~day+I(day^2), random=~1|Side, 
                    opt="optim"), data=Pixel)
summary(t1)
getSmo(t1)
# fixed part in gamlss formula
t1a <- gamlss(pixel~day+I(day^2) + re(random=~1|Side), 
            data=Pixel, 
            opt="optim")
summary(t1a)
getSmo(t1a)

plot(fitted(l1), fitted(t1))
plot(fitted(l1), fitted(t1a))
# with smoother
t1b <- gamlss(pixel~cs(day, df = 1) + re(random=~1|Side), 
              data=Pixel, 
              opt="optim")
summary(t1b)
getSmo(t1b)

plot(fitted(t1a), fitted(t1b))
