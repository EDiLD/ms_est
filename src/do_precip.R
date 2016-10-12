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

(keep <- props[prop > 0.05 & tot > 1000])
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

rm(samples_rac, psm_sites, psm_sites_info, take_si, rac, keep_tab, keep_tab_x)


# Model -------------------------------------------------------------------
# one single compouns
take_c <- take[variable_id == 153]


gco <- glim.control(glm.trace = TRUE, bf.trace = FALSE)
# decrease number of inner cycles to speed up computations
glco <- gamlss.control(c.crit = 0.01, mu.step = 0.5) 
# increacse c.crit from 0.001 to speed up convergence for large data
# decrease step length, to omit improper parameters (slower convergence).



# gamma hurdle model with precipitation as groups, season and psm_type as predictors
# state as random effect
mod_l_m <- gamlss(rq ~ 0  + log_precip_1 + log_precip0 + season +
                    re(random=~1|state_fac/s_id_fac),
                  nu.formula =~ 0 + precip_1 + precip0 + season + 
                    re(random=~1|state_fac/s_id_fac),
                  data = take_c,
                  family = ZAGA,
                  control = glco,
                  i.control = gco)

take_c[season == 'Q4', rq]


plot(mod_l_m)
summary(mod_l_m)
term.plot(mod_l_m)
res <- residuals(mod_l_m)
hist(res)

# random slope model, does not converge!
# mod_l_m_ra <- gamlss(rq ~ log_precip_1 + log_precip0 + season +
#                     re(random=~1|state_fac/s_id_fac) +
#                       # random slope
#                       re(random=~precip_1 + precip0 + season|v_id_fac),
#                   # model also pi with same predictors
#                   nu.formula =~precip_1 + precip0 + season +
#                     re(random=~1|state_fac/s_id_fac) + re(random=~precip_1 + precip0 + season|v_id_fac),
#                   # sigma is constant
#                   data = take,
#                   family = ZAGA,
#                   i.control = gco)



# fit model for each compound
model_foo <- function(var){
  message('Running model on compound: ', var)
  take_c <<- take[variable_id == var]
  mod <- gamlss(rq ~ 0+log_precip_1 + log_precip0 + season +
                      re(random = ~1|state_fac/s_id_fac),
                    nu.formula =~0+log_precip_1 + log_precip0 + season + 
                      re(random = ~1|state_fac/s_id_fac),
                    data = take_c,
                    family = ZAGA,
                    control = glco,
                    i.control = gco)
  saveRDS(mod, file.path(cachedir, 'lmodels', paste0('mod_', var, '.rds')))
}

# model_foo(727)
# run model on compounds
run_model <- FALSE
if (run_model){
  lapply(keep$variable_id, model_foo)
}


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
  # invert estimate for nu 
  # (so that a increase corresponds to a increase in p(x>LOQ))
  # standard is that a increase corresponse to p(x<LOQ)
  # nu = b * x  | nu = 1-pi
  # 1-pi = b*x
  # p = -b*x - 1
  smod$Estimate[grepl('nu\\.', smod$terms)] <- -smod$Estimate[grepl('nu\\.', smod$terms)]
  smod$upci <- smod$Estimate + smod$Std..Error * mult
  smod$lowci <- smod$Estimate - smod$Std..Error * mult
  return(smod)
}

res <- lapply(file.path(cachedir, 'lmodels', paste0('mod_', keep$variable_id, '.rds')), model_extr)


resdf <- rbindlist(res)
names(resdf) <- c('est', 'stderr', 'tval', 'pval', 'term', 'variable', 'upci', 'lowci')
resdf$term_type <- gsub('^(.*)\\..*', '\\1', resdf$term)
resdf$term2 <- gsub('^.*\\.(.*)$', '\\1', resdf$term)
resdf$pval <- round(resdf$pval, 4)

# indicator for significant based on p
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
keep_tab2[ sign(upci)  == sign(lowci) & !grepl('season', term2), est := make_bold(est)]
keep_tab2[ , estci := paste0(est, '\\newline (', lowci, ' - ', upci, ')')]
keep_tab2[ , term2_p := mapvalues(term2, c("log_precip0", "log_precip_1", "seasonQ1", "seasonQ2", 
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
rownames(keep_tab2_w)<-NULL

keep_tab2_x <- xtable(keep_tab2_w, 
                     label = 'tab:var_model_coef',
                     caption = 'Coefficients and CI from per compound models. 
                     Bold values denote coefficients where the CI for precipitation encompasses zero.
                     Coefficients are on the link scale (log for $mu$ and logit for $\nu$).',
                     align = 'lp{2cm}p{0.6cm}p{1.8cm}p{1.8cm}p{1.8cm}p{1.8cm}p{1.8cm}p{1.8cm}')

print(keep_tab2_x, 
      file = file.path(prj, 'supplement/keeptab2.tex'),
      tabular.environment="longtable",
      floating = FALSE,
      caption.placement = 'top',
      comment = FALSE,
      booktabs = TRUE,
      hline.after = c(-1, 0, 23, 46),
      sanitize.text.function = identity,
      size="\\fontsize{8pt}{10pt}\\selectfont"
)




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
  labs(x = '', y = '') +
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
  scale_x_discrete(breaks = c('seasonQ4', 'seasonQ3', 'seasonQ2', 'seasonQ1'),
                   labels = c(' Quarter Q4', ' Quarter Q3', ' Quarter Q2', ' Quarter Q1')) +
  theme(strip.text.x = element_text(size = 22))

p <- arrangeGrob(p_precip, p_season, ncol = 1)
# plot(p)
ggsave(file.path(prj, "supplement", "coefs.pdf"), p, width = 10, height = 9)




# # plot for talk2
# psm_variables[variable_id %in% keep$variable_id, list(variable_id, name, psm_type)]
# sort(table(take[variable_id == 349, site_id]))
# 
# df <- take[variable_id == 349 & site_id == 'RP_2375572100']
# take_df <- psm_samples[sample_id %chin% df$sample_id & variable_id == 349]
# p <- ggplot(take_df, aes(x = as.Date(date), y = value_fin)) + 
#   geom_point() +
#   geom_hline(aes(yintercept = 0.05), linetype = 'dashed', col = 'darkorange', size = 1) +
#   mytheme +
#   labs(y = 'Glyphosate [ug/L]', x = 'Date') +
#   scale_x_date(date_breaks = '1 year', date_labels = '%Y') +
#   ggtitle('Erlenbach / Rheinzabern')
# 
# require(ggExtra)
# pout <- ggMarginal(p, type = 'histogram', margins = 'y')
# ggsave('/home/user/Documents/projects_git/talk_work2/fig/glyph.pdf', pout, width = 7, height = 5)




### ------------------------------------------------------------------------
# metaanalysis of coefficients

# test on one coef
unique(resdf$term)
# fixed effect meta analysis
mmod <- rma(est, sei=stderr, data = resdf[term == 'nu.log_precip_1'], method = 'FE')
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
mmod <- rma(est, sei=stderr, data = resdf[term == 'nu.log_precip_1'], method = 'REML')
mmod
plot(mmod)



# fit top each term
terms <- unique(resdf$term)
# terms <- terms[!(terms == 'sigma' | grepl('Intercept', terms))]
fit_meta <- function(tm){
  mmod <- rma(est, sei=stderr, data = resdf[term == tm], method = 'REML')
  out <- data.frame(term = tm, est = mmod$b[,1], upr = mmod$ci.ub, lwr = mmod$ci.lb)
  return(out)
}

resm <- lapply(terms, fit_meta)
resm
resmd <- rbindlist(resm)

# log_precip1 multiplicator for mu
exp(0.05582678) - 1
# one <- exp(-3.41532436) * exp(0.05582678 * 1)
# ten <- exp(-3.41532436) * exp(0.05582678 * 2)
# ten_one = exp(-3.41532436) * (exp(0.05582678 * 2) - exp(0.05582678 * 1))
# (ten-one) / one

# nu 
int <- -2.21435302
a <- exp(int + 0.27717549)
a / (1+a) # prob at 1mm)
(one <- plogis(int + 0.27717549))
(ten <- plogis(int + 0.27717549*2))
(hun <- plogis(int + 0.27717549*3))
(ten-one) 
(hun-ten)
(ten-one) / one 
(hun-ten) / ten

# Q1 mu
exp(-3.80443384)
# Q3 mu
exp(-3.40903988)
# Q1 nu
plogis(-3.15180443)
plogis(-2.21435302)



# coefplot
resmd$type <- gsub('^(.*)\\.(.*)$', '\\1', resmd$term)
resmd$coeftype <- ifelse(grepl('season', resmd$term), 'season', 'precip')
# remove sigma
resmd <- resmd[!term == 'sigma']
p_season <- ggplot(resmd[resmd$coeftype == 'season', ]) +
  geom_pointrange(aes(x = term, y = est, ymin = lwr, ymax = upr)) + 
  coord_flip() +
  facet_grid(.~type , scales ='free_x', labeller = label_parsed) +
  mytheme +
  labs(x = '', y = 'Coefficient') +
  scale_x_discrete(breaks = c('mu.seasonQ4', 'mu.seasonQ3', 'mu.seasonQ2', 'mu.seasonQ1'),
                   labels = c(' Quarter Q4', ' Quarter Q3', ' Quarter Q2', ' Quarter Q1')) +
  theme(strip.text.x = element_blank())


p_precip <- ggplot(resmd[resmd$coeftype == 'precip', ]) +
  geom_pointrange(aes(x = term, y = est, ymin = lwr, ymax = upr)) + 
  coord_flip() +
  geom_hline(aes(yintercept = 0), linetype = 'dotted') +
  facet_grid(.~type , scales ='free_x', labeller = label_parsed) +
  mytheme +
  labs(x = '', y = '') +
  scale_x_discrete(breaks = c('mu.log_precip_1', 'mu.log_precip0'),
                   labels = c(expression('log'~precip[-1]), expression('log'~precip[0]))) +
  theme(strip.text.x = element_text(size = 22))

p <- arrangeGrob(p_precip, p_season, ncol = 1, heights = c(1, 1.25))
# plot(p)
ggsave("figure5.pdf", p, width = 7, height = 5.5)


# # extract random effect variances
# # function to extract the needed model components, 
# ra_extr <- function(file){
#   message('Working on file: ', file)
#   mod <- readRDS(file)
#   lme_out <- mod$mu.coefSmo[[1]]
#   vc <- VarCorr(lme_out)
#   suppressWarnings(storage.mode(vc) <- 'numeric')
#   vc <- vc[c(2, 4), 'StdDev']
#   names(vc) <- c('state/site', 'state')
#   vc <- c(vc, var = gsub('mod_(.*)\\.rds', '\\1', basename(file)))
#   # generalized rsq
#   vc <- c(vc, rsq = Rsq(mod, type = "Cox Snell"))
#   # calculate CIs
#   return(vc)
# }
# 
# ra_res <- lapply(file.path(cachedir, 'models', paste0('mod_', keep$variable_id, '.rds')), ra_extr)
# ra_res <- do.call(rbind, ra_res)
# ra_res <- apply(ra_res, 2, as.numeric)
# # range of rsq values
# range(ra_res[,4])
