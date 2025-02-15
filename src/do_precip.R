if (!exists("prj")) {
  stop("You need to create a object 'prj' that points to the top folder, 
       e.g. prj <- '/home/edisz/Documents/work/research/projects/2016/4BFG/Paper/ms_est' or
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


# restrict to sites < 100km² (or unknown)  and both data available (see do_overview.R)
take_site_id <- readRDS(file = file.path(cachedir, 'take_site_id.rds'))
psm_sites_info <- psm_sites_info[site_id %in% take_site_id]
psm_sites <- psm_sites[site_id %in% take_site_id]
psm_samples <- psm_samples[site_id %in% take_site_id]




# Risk Quotients ----------------------------------------------------------
# restrict samples to variables with rac
rac <- psm_variables[!is.na(rak_uba), list(variable_id, name, cas, pgroup, rak_uba)]
# join with racs
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
take_si <- psm_sites_info[!is.na(ezg_fin) | is.na(agri_fin) , 
                          list(site_id, agri_fin, ezg_fin)]
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


# set season / quarters
take[ , season := cut(month(date), breaks = c(0.5, 3.5, 6.5, 9.5, 12.5), 
                      labels = c('Q1', 'Q2', 'Q3', 'Q4'))]

# add group
setkey(take, variable_id)
setkey(psm_variables, variable_id)
take <- psm_variables[ , list(variable_id, psm_type)][take]


# remove compounds with no risk
props <- take[  , list(prop = sum(rq > 0) / length(rq),
              abs = sum(rq > 0),
              tot = length(rq)), by = variable_id][order(prop, 
                                                         decreasing = TRUE)]
hist(props[prop < 0.2 , prop])
# keep compounds that were found in at least 5% of samples
# and with at least 500 observations
hist(props[tot < 1000, tot])

(keep <- props[prop > 0.05 & tot > 1000])
# 23 compounds left





#  ------------------------------------------------------------------------
# export overview for supplement
keep
setkey(psm_variables, variable_id)
setkey(keep, variable_id)
keep_tab <- psm_variables[ , list(variable_id, name, cas, psm_type)][keep]
keep_tab[ , prop := round(prop*100, 2)]
setnames(keep_tab, c('id', 'Name', 'CAS', 'Group', '\\%>LOQ', 
                     'no. > LOQ', 'total no.'))
keep_tab$id <- NULL
keep_tab

keep_tab_x <- xtable(keep_tab, 
                    label = 'tab:var_model',
                    caption = c('23 pesticides for which we modelled the relationship between RQ and precipitation and seasonality, respectively.
                    Order is the same as in Figure 5 of the main text. See Table \\ref{tab:var_model_coef} for model coefficients.',
                                '23 pesticides for which we modelled the relationship with precipitation and seasonality.'),
                    align = 'lp{2.5cm}rlp{1.5cm}p{2cm}p{2cm}')

print(keep_tab_x,
      file = file.path(prj, 'supplement/keeptab.tex'),
      tabular.environment = "longtable",
      floating = FALSE,
      caption.placement = 'top',
      comment = FALSE,
      booktabs = TRUE,
      hline.after = c(-1, 0, nrow(keep_tab)),
      sanitize.text.function = identity#,
#      size="\\fontsize{8pt}{10pt}\\selectfont"
)






# subset
take <- take[variable_id %in% keep$variable_id]
psm_variables[variable_id %in% keep$variable_id, list(variable_id, name, 
                                                      psm_type)]

rm(samples_rac, psm_sites, psm_sites_info, take_si, rac, keep_tab, keep_tab_x)


# Model -------------------------------------------------------------------
# single compouds (for testing)
# take_c <- take[variable_id == 191]
take_c <- take[variable_id == 378]


gco <- glim.control(glm.trace = TRUE, bf.trace = FALSE)
# decrease number of inner cycles to speed up computations
glco <- gamlss.control(c.crit = 0.01, mu.step = 0.5) 
# increacse c.crit from 0.001 to speed up convergence for large data
# decrease step length, to omit improper parameters (slower convergence).



# # gamma hurdle model with precipitation as groups, season and psm_type as preds
# # state as random effect
# mod_l_m <- gamlss(rq ~ 0  + log_precip_1 + log_precip0 + season +
#                     re(random = ~1|state_fac/s_id_fac),
#                   nu.formula = ~ 0 + precip_1 + precip0 + season + 
#                     re(random = ~1|state_fac/s_id_fac),
#                   data = take_c,
#                   family = ZAGA,
#                   control = glco,
#                   i.control = gco)
# plot(mod_l_m)
# summary(mod_l_m)
# term.plot(mod_l_m)
# term.plot(mod_l_m, what = 'nu')
# res <- residuals(mod_l_m)
# hist(res)
# 
# foo <- function(x) log10(x + 0.05) # transformation
# newdata = data.frame(log_precip0 = foo(0.1),
#                      log_precip_1 = foo(0.1),
#                      season = 'Q2')
# predict(mod_l_m,  what = 'nu',  type = 'response')
# 
# 
# 
# # plot marginal data
# plot(log(rq) ~ log_precip_1 , data = take_c[take_c$rq > 0, ], pch = 16, col = rgb(0, 0, 0, 0.2))
# plot(I(ifelse(rq > 0, 1, 0)) ~ log_precip_1 , data = take_c, 
#      col = rgb(0, 0, 0, 0.1), pch = 16)
# plot(log(rq) ~ season , data = take_c[take_c$rq > 0, ])


# random slope model, does not converge!
# mod_l_m_ra <- gamlss(rq ~ log_precip_1 + log_precip0 + season +
#                     re(random=~1|state_fac/s_id_fac) +
#                       # random slope
#                       re(random=~precip_1 + precip0 + season|v_id_fac),
#                   # model also pi with same predictors
#                   nu.formula =~precip_1 + precip0 + season +
#                     re(random=~1|state_fac/s_id_fac) + re(random=~precip_1 + 
#                                                   precip0 + season|v_id_fac),
#                   # sigma is constant
#                   data = take,
#                   family = ZAGA,
#                   i.control = gco)



# fit model for each compound
model_foo <- function(var){
  message('Running model on compound: ', var)
  take_c <<- take[variable_id == var]
  mod <- gamlss(rq ~ 0 + log_precip_1 + log_precip0 + season +
                      re(random = ~1|state_fac/s_id_fac),
                    nu.formula = ~0+log_precip_1 + log_precip0 + season + 
                      re(random = ~1|state_fac/s_id_fac),
                    data = take_c,
                    family = ZAGA,
                    control = glco,
                    i.control = gco)
  # reduced models 
  mod_rmu <- update(mod, ~1, what = 'mu')
  mod_rnu <- update(mod, ~1, what = 'nu')
  mod_rfull <- update(mod, ~1, what = 'All')
  out <- list(mod = mod,
              mod_rmu = mod_rmu,
              mod_rnu = mod_rnu,
              mod_rfull = mod_rfull)
  
  saveRDS(out, file.path(cachedir, 'lmodels', paste0('mod_', var, '.rds')))
}

# model_foo(727)
# run model on compounds
run_model <- FALSE
# run_model <- TRUE
if (run_model) {
  file.remove(file.path(cachedir, 'lmodels', paste0('mod_', keep$variable_id, '.rds')))
  lapply(keep$variable_id, model_foo)
}




# function to extract the needed model components,
model_extr <- function(file){
  cont <- readRDS(file)
  mod <- cont[['mod']]
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
  # invert estimate for nu 
  # (so that a increase corresponds to a increase in p(x>LOQ))
  # standard is that a increase corresponse to p(x<LOQ)
  # nu = b * x  | nu = 1-pi
  # 1-pi = b*x
  # p = -b*x - 1
  smod$Estimate[grepl('nu\\.', smod$terms)] <- -smod$Estimate[grepl('nu\\.', smod$terms)]
  smod$upci <- smod$Estimate + smod$Std..Error * mult
  smod$lowci <- smod$Estimate - smod$Std..Error * mult
  
  # compare wih reduce models
  mod_rmu <- cont[['mod_rmu']]
  mod_rnu <- cont[['mod_rnu']]
  mod_rfull <- cont[['mod_rfull']]
  lliks <- c(mod = logLik(mod),
             mod_rmu = logLik(mod_rmu),
             mod_rnu = logLik(mod_rnu),
             mod_rfull = logLik(mod_rfull))
  
  out <- list(smod = smod, lliks = lliks)
  return(out)
}

# extract coefs from each model
res <- lapply(file.path(cachedir, 'lmodels', 
                        paste0('mod_', keep$variable_id, '.rds')), 
              model_extr)


smods <- lapply(res, '[[', 'smod')
resdf <- rbindlist(smods)
# nice formatting
names(resdf) <- c('est', 'stderr', 'tval', 'pval', 'term', 'variable', 
                  'upci', 'lowci')
resdf$term_type <- gsub('^(.*)\\..*', '\\1', resdf$term)
resdf$term2 <- gsub('^.*\\.(.*)$', '\\1', resdf$term)
resdf$pval <- round(resdf$pval, 4)

# indicator for significant based on p
resdf$estsig <- ifelse(resdf$pval < 0.05, resdf$est, NA)
# indicatopr for significant based on ci
resdf$cisig <- ifelse(sign(resdf$upci)  == sign(resdf$lowci), 
                      'cisig', 
                      'nocisig')
# parameter indictor
resdf$termind <- substr(resdf$term, 1, 2)





# Check deviances.
# ND = Sat - Null
# RD = Sat - Mod
# Explained deviance = 1 - ND / RD = (Mod - Null) / (Sat - Null)
# 
# But saturated model not available (computationally very intensive)
#  But difference from NUll gives indication (assuming that SAT - Null only scales)

lliks <- lapply(res, '[[', 'lliks')
lliks <- lapply(lliks, function(y) data.frame(t(y)))
lliks <- rbindlist(lliks)
lliks[ , variable_id := keep$variable_id]
# add names
lliks <- psm_variables[ , list(variable_id, name)][lliks, on = "variable_id"]
lliks
# # calculate proportional (to Null) increase
# lliks[ , c("dmod","dmu", "dnu") := list(dmod = (mod - mod_rfull),
#               dmu = (mod_rmu - mod_rfull),
#               dnu = (mod_rnu - mod_rfull)
#               )]
lliks
require(tidyr)
# to long format
lliks <- gather(lliks, key = model, value = value, -c(1:2))
lliks$value <- lliks$value - lliks$value[lliks$model == 'mod']
p_loglik <- ggplot(lliks, aes(y = value - 0.001, x = name, col = model, group = model)) +
  geom_line() +
  geom_point() +
  labs(y = 'log-Likelihood', x = '') +
  scale_color_brewer('Model',
                     palette = 'Dark2', 
                     breaks = c('mod', 'mod_rmu', 'mod_rnu', 'mod_rfull'),
                     labels = c('Full', 'Reduce mu', 'Reduce nu', 'Reduce all')) +
  coord_flip() +
  # scale_y_continuous(breaks = NULL) +
  labs(y = 'Likelihood') 
p_loglik
# ggsave('/home/edisz/Documents/work/research/projects/2016/1PHD/phd_defense/figs/logliks.pdf',
#        p_loglik,
#        width = 7, height = 8)

# example model verification:
# for glyphosate
takevar = 349
cont <- readRDS(file.path(cachedir, 'lmodels', 
                          paste0('mod_', takevar, '.rds')))
mod <- cont[['mod']]
plot(fitted(mod)[take[variable_id == takevar, rq] > 0], take[variable_id == takevar & rq > 0, rq], log = 'xy')
abline(0, 1)

pdat <- take[variable_id == takevar & rq > 0, ]
pdat$ftd <- fitted(mod)[take[variable_id == takevar, rq] > 0]
p_ftd_vs_obs <- ggplot(pdat, aes(x = rq, y = ftd)) +
  theme_edi() +
  geom_hex(binwidth = 0.1, alpha = 0.8) +
  # geom_point(alpha = 0.3) +
  geom_abline() +
  scale_x_log10(breaks = c(0.001, 0.01, 0.1, 1)) +
  scale_y_log10(breaks = c(0.001, 0.01, 0.1, 1)) +
  scale_fill_viridis('n') +
  labs(x = 'Observed RQ', y = 'Estimated RQ') +
  ggtitle(label = 'Observed vs. Fitted', subtitle = 'Glyphosate, only RQ > 0 shown, n = 1389')
p_ftd_vs_obs
# 
# ggsave('/home/edisz/Documents/work/research/projects/2016/1PHD/phd_defense/figs/pftdvsobs.pdf',
#        p_ftd_vs_obs,
#        width = 8, height = 7)


#  ------------------------------------------------------------------------
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
keep_tab2[ sign(upci)  == sign(lowci) & !grepl('season', term2), 
           est := make_bold(est)]
keep_tab2[ , estci := paste0(est, '\\newline (', lowci, ' - ', upci, ')')]
keep_tab2[ , term2_p := mapvalues(term2, c("log_precip0", "log_precip_1", 
                                           "seasonQ1", "seasonQ2", 
                                           "seasonQ3", "seasonQ4"),
                       c("$log~precip_0$", "$log~precip_{-1}$", "Quarter 1", 
                         "Quarter 2", "Quarter 3", "Quarter 4"))]

keep_tab2[ , term_type_p := mapvalues(term_type, c('mu', 'nu'),
                                     c('$\\mu$', '$\\nu$'))]
setnames(keep_tab2, c("variable_id", "Compound", "term2", "term_type", "est", 
                      "lowci", "upci", "coefficient", 
                      "variable", "effect"))

keep_tab2_w <- keep_tab2[ , list(Compound, variable, effect, coefficient)]
keep_tab2_w <- dcast(keep_tab2_w, Compound + effect ~ variable)
keep_tab2_w <- keep_tab2_w[order(keep_tab2_w$effect), ]
rownames(keep_tab2_w) <- NULL

keep_tab2_x <- xtable(keep_tab2_w, 
                     label = 'tab:var_model_coef',
                     caption = c('Coefficients and CI from per compound models. 
                     Bold values denote coefficients where the CI for precipitation encompasses zero.
                     Coefficients are on the link scale (log for $\\mu$ and logit for $\\nu$).',
                     'Coefficients and CI from per compound models.'),
                     align = 'lp{2cm}p{0.6cm}p{1.8cm}p{1.8cm}p{1.8cm}p{1.8cm}p{1.8cm}p{1.8cm}')

print(keep_tab2_x, 
      file = file.path(prj, 'supplement/keeptab2.tex'),
      tabular.environment = "longtable",
      floating = FALSE,
      caption.placement = 'top',
      comment = FALSE,
      booktabs = TRUE,
      hline.after = c(-1, 0, 22, 44),
      sanitize.text.function = identity,
      size = "\\fontsize{8pt}{10pt}\\selectfont"
)




#  ------------------------------------------------------------------------
# display estimates and errors
pdata <- resdf[!term2 %chin% c('(Intercept)', 'sigma') &
                 term2 %in% c('log_precip0', 'log_precip_1')]
pdata[ , variable := factor(variable, levels = rev(sort(
                                                unique(pdata[ , variable]))))]
p_precip <- ggplot(data = pdata) +
  geom_pointrange(aes(x = term2, y = est, ymax = upci, ymin = lowci, 
                      fill = variable,
                      col = cisig),
                  position = position_dodge(width = .6)) +
  geom_hline(aes(yintercept = 0), linetype = 'dotted') +
  facet_wrap(~term_type, labeller = label_parsed) +
  coord_flip() +
  mytheme +
  theme(legend.position = "none") +
  scale_color_manual(values = c('black', 'grey70')) +
  labs(x = '', y = '') +
  scale_x_discrete(breaks = c('log_precip_1', 'log_precip0'),
                   labels = c(expression('log'~prec[-1]), 
                              expression('log'~prec[0]))) +
  theme(strip.text.x = element_text(size = 22))

pdata2 <- resdf[!term2 %chin% c('(Intercept)', 'sigma') & 
                  !term2 %in% c('log_precip0', 'log_precip_1')]
pdata2[ , variable := factor(variable, levels = rev(sort(
                                                unique(pdata2[ , variable]))))]
p_season <- ggplot(data = pdata2) +
  geom_pointrange(aes(x = term2, y = est, ymax = upci, ymin = lowci, 
                      fill = variable,
                      col = cisig),
                  position = position_dodge(width = .6)) +
  coord_flip() +
  facet_wrap(~term_type,labeller = label_parsed) +
  geom_hline(aes(yintercept = 0), linetype = 'dotted') +
  mytheme +
  theme(legend.position = "none") +
  scale_color_manual(values = c('black', 'grey70')) +
  labs(x = '', y = 'Coefficient') +
  scale_x_discrete(breaks = c('seasonQ4', 'seasonQ3', 'seasonQ2', 'seasonQ1'),
                   labels = c(' Quarter Q4', ' Quarter Q3', 
                              ' Quarter Q2', ' Quarter Q1')) +
  theme(strip.text.x = element_text(size = 22))

p <- arrangeGrob(p_precip, p_season, ncol = 1)
# plot(p)
ggsave(file.path(prj, "supplement", "coefs.pdf"), p, width = 10, height = 9)

# phd <- '/home/edisz/Documents/work/research/projects/2016/1PHD/phd_thesis/appendix/smallstreams/one/'
# ggsave(file.path(phd,  "coefs.pdf"), p, width = 170, height = 170, unit = 'mm')


# # plot for talk2
# psm_variables[variable_id %in% keep$variable_id, list(variable_id, name, psm_type)]
# sort(table(take[variable_id == 349, site_id]))
# 
# df <- take[variable_id == 349 & site_id == 'RP_2375572100']
# take_df <- psm_samples[sample_id %chin% df$sample_id & variable_id == 349]
# p <- ggplot(take_df, aes(x = as.Date(date), y = value_fin)) + 
#   geom_point() +
#   geom_hline(aes(yintercept = 0.05), linetype = 'dashed', 
#                     col = 'darkorange', size = 1) +
#   mytheme +
#   labs(y = 'Glyphosate [ug/L]', x = 'Date') +
#   scale_x_date(date_breaks = '1 year', date_labels = '%Y') +
#   ggtitle('Erlenbach / Rheinzabern')
# 
# require(ggExtra)
# pout <- ggMarginal(p, type = 'histogram', margins = 'y')
# ggsave('/home/user/Documents/projects_git/talk_work2/fig/glyph.pdf', pout, 
#             width = 7, height = 5)




### ------------------------------------------------------------------------
# metaanalysis of coefficients

# test on one coef
unique(resdf$term)
# fixed effect meta analysis
mmod <- rma(est, sei = stderr, data = resdf[term == 'nu.log_precip_1'], 
            method = 'FE')
mmod
plot(mmod)


# by hand (See Harrison (2011), MEE, eqn 1+2)
dat <- resdf[term == 'nu.log_precip_1']
#+ weights
wi <- 1/dat$stderr^2
# weighted mean
sum(wi*dat$est) / sum(wi)
# CI
sum(wi*dat$est) / sum(wi) + sqrt(1/sum(wi)) * qnorm(0.975)
sum(wi*dat$est) / sum(wi) - sqrt(1/sum(wi)) * qnorm(0.975)


# random effect meta analysis
mmod <- rma(est, sei = stderr, data = resdf[term == 'nu.log_precip_1'], 
            method = 'REML')
mmod
plot(mmod)



# fit top each term
terms <- unique(resdf$term)
# terms <- terms[!(terms == 'sigma' | grepl('Intercept', terms))]
fit_meta <- function(tm){
  mmod <- rma(est, sei = stderr, data = resdf[term == tm], method = 'REML')
  out <- data.frame(term = tm, est = mmod$b[,1], upr = mmod$ci.ub, 
                    lwr = mmod$ci.lb)
  return(out)
}

resm <- lapply(terms, fit_meta)
resm
resmd <- rbindlist(resm)





# Absolute effect sizes ---------------------------------------------------


# LOQ ---------------------------------------------------------------------
foo <- function(x) log10(x + 0.05) # transformation
# log_precip_1: nu 
# intercept
int <- resmd$est[resmd$term ==  'nu.seasonQ2']
plogis(int) *100

# prop of exceeding loq in q2 at 0.1 mm
(zerop1 <- plogis(int + foo(0.1)*resmd$est[resmd$term ==  'nu.log_precip_1'])) * 100
(one <- plogis(int + foo(1)*resmd$est[resmd$term ==  'nu.log_precip_1'])) * 100
# prop of exceeding loq in q2 at 10 mm
(ten <- plogis(int + foo(15)*resmd$est[resmd$term ==  'nu.log_precip_1'])) * 100

(ten - zerop1) 
(ten - zerop1) / zerop1 * 100


# or individual models
zerop <- 100 * plogis(resdf[termind == 'nu' & term2 == 'seasonQ2', est] + foo(0.1) * resdf[termind == 'nu' & term2 == 'log_precip_1', est])
range(zerop)
tena <- 100 * plogis(resdf[termind == 'nu' & term2 == 'seasonQ2', est] + foo(10) * resdf[termind == 'nu' & term2 == 'log_precip_1', est])
range(tena)
(tena - zerop) / zerop * 100

# Q1 nu
plogis(resmd$est[resmd$term ==  'nu.seasonQ1'])
# Q2 nu
plogis(resmd$est[resmd$term ==  'nu.seasonQ2'])






# RQ ----------------------------------------------------------------------
# Q1 mu
exp(resmd$est[resmd$term ==  'mu.seasonQ1'])
# Q2 mu
exp(resmd$est[resmd$term ==  'mu.seasonQ2'])

int <- resmd$est[resmd$term ==  'mu.seasonQ2']
# predictions in q2 at 0.1 mm
(zerop1 <- exp(int + foo(0.1)*resmd$est[resmd$term ==  'mu.log_precip_1']))
# predictions in q2 at 15 mm
(ten <- exp(int + foo(15)*resmd$est[resmd$term ==  'nu.log_precip_1']))
(ten - zerop1) 
(ten - zerop1) / zerop1 * 100 


# or individual models
zerop <- exp(resdf[termind == 'mu' & term2 == 'seasonQ2', est] + foo(0.1) * resdf[termind == 'mu' & term2 == 'log_precip_1', est])
range(zerop)
tena <- exp(resdf[termind == 'mu' & term2 == 'seasonQ2', est] + foo(10) * resdf[termind == 'mu' & term2 == 'log_precip_1', est])
range(tena)
(tena - zerop)
range((tena - zerop))
(tena - zerop) / zerop * 100
resdf[termind == 'mu' & term2 == 'seasonQ2', variable]



# Plot of meta-analysis ---------------------------------------------------
# coefplot
resmd$type <- gsub('^(.*)\\.(.*)$', '\\1', resmd$term)
resmd$coeftype <- ifelse(grepl('season', resmd$term), 'season', 'precip')
# remove sigma
resmd <- resmd[!term == 'sigma']
resmd$type <- mapvalues(resmd$type, c('mu', 'nu'),
                        c('mu~(RQ)', 'nu~(LOQ)'))
resmd$term2 <- substr(resmd$term, nchar(resmd$term) - 1, nchar(resmd$term))


p_season <- ggplot(resmd[resmd$coeftype == 'season', ]) +
  geom_pointrange(aes(x = substr(term, nchar(term) - 1, nchar(term)), y = est, ymin = lwr, ymax = upr)) + 
  coord_flip() +
  facet_grid(.~type , labeller = label_parsed) +
  mytheme +
  labs(x = '', y = 'Coefficient') +
  scale_x_discrete(breaks = c('Q4', 'Q3', 
                              'Q2', 'Q1'),
                   labels = c(' Quarter Q4', ' Quarter Q3', 
                              ' Quarter Q2', ' Quarter Q1')) +
  theme(strip.text.x = element_blank())


p_precip <- ggplot(resmd[resmd$coeftype == 'precip', ]) +
  geom_pointrange(aes(x = gsub("^.*?\\.(.*)", "\\1", term), y = est, ymin = lwr, ymax = upr)) + 
  coord_flip() +
  geom_hline(aes(yintercept = 0), linetype = 'dotted') +
  facet_grid(.~type , labeller = label_parsed) +
  mytheme +
  labs(x = '', y = '') +
  scale_x_discrete(breaks = c('log_precip_1', 'log_precip0'),
                   labels = c(expression('log'~precip[-1]), 
                              expression('log'~precip[0]))) +
  theme(strip.text.x = element_text(size = 22))

p <- arrangeGrob(p_precip, p_season, ncol = 1, heights = c(1, 1.25))
plot(p)
# plot(p)
ggsave(file.path(prj, "supplement/mean_coef.pdf"), p, width = 7, height = 5.5)

# phd <- '/home/edisz/Documents/work/research/projects/2016/1PHD/phd_thesis/chapters/smallstreams/'
# ggsave(file.path(phd, "figure5.pdf"), p, width = 7, height = 5.5)




# plot model predictions.
pdata <- resmd[resmd$coeftype == 'season' & resmd$type == 'nu~(LOQ)', ]
# calculate for high precipitation
int <- pdata$est
pdata$high <- plogis(int + foo(15) * resmd[resmd$term == 'nu.log_precip_1' & resmd$type == 'nu~(LOQ)', ]$est) 

p <- ggplot() +
  geom_pointrange(data = pdata,
                  aes(x = substr(term, nchar(term) - 1, nchar(term)),
                      y = plogis(est),
                      ymin = plogis(lwr),
                      ymax = plogis(upr)),
                  size = 0.5) +
  geom_point(data = pdata,
             aes(x = substr(term, nchar(term) - 1, nchar(term)),
                 y = high),
             size = 3,
             col = rgb(0.847, 0.510, 0.106)) +
  scale_x_discrete(breaks = c('Q4', 'Q3',
                              'Q2', 'Q1'),
                   labels = c('Oct-Dec', 'Jul-Sep',
                              'Apr-Jun', 'Jan-Mar')) +
  mytheme + 
  labs(x = '', y = 'p(x > LOQ)') 
p
ggsave(file.path(prj, "figure5.pdf"), p, width = 4, height = 3.5)

# # plot for defense
# library(tikzDevice)
# p <- ggplot() +
#   geom_pointrange(data = pdata,
#                   aes(x = substr(term, nchar(term) - 1, nchar(term)),
#                       y = plogis(est),
#                       ymin = plogis(lwr),
#                       ymax = plogis(upr)),
#                   size = 1) +
#   geom_point(data = pdata,
#              aes(x = substr(term, nchar(term) - 1, nchar(term)),
#                  y = high),
#              size = 4,
#              col = rgb(0.847, 0.510, 0.106)) +
#   scale_x_discrete(breaks = c('Q4', 'Q3',
#                               'Q2', 'Q1'),
#                    labels = c('Oct-Dec', 'Jul-Sep',
#                               'Apr-Jun', 'Jan-Mar')) +
#   theme_edi(base_size = 18) +
#   labs(x = '', y = 'p(x > LOQ)') +
#   ggtitle('Model predictions',
#           subtitle = 'n = 23 compounds, orange: 15 mm precipitation.')
# p
# ggsave('/home/edisz/Documents/work/research/projects/2016/1PHD/phd_defense/figs/tikz/p_season.tikz',
#        p,
#        device = tikz,
#        width = 6)


p <- ggplot() +
  geom_pointrange(data = pdata,
                  aes(x = substr(term, nchar(term) - 1, nchar(term)),
                      y = plogis(est),
                      ymin = plogis(lwr),
                      ymax = plogis(upr)),
                  size = 1) +
  geom_point(data = pdata,
             aes(x = substr(term, nchar(term) - 1, nchar(term)),
                 y = high),
             size = 4,
             col = rgb(0.847, 0.510, 0.106)) +
  scale_x_discrete(breaks = c('Q4', 'Q3',
                              'Q2', 'Q1'),
                   labels = c('Oct-Dec', 'Jul-Sep',
                              'Apr-Jun', 'Jan-Mar')) +
  mytheme + 
  labs(x = '', y = 'p(x > LOQ)') 
  # ggtitle('Model predictions',
  #         subtitle = 'n = 23 compounds, orange: 15 mm precipitation.')
p
