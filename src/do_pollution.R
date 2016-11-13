if (!exists("prj")) {
  stop("You need to create a object 'prj' that points to the top folder, 
       e.g. prj <- '/home/edisz/Documents/work/research/projects/2016/4BFG/Paper/ms_est' or
       prj <- '/home/user/Documents/projects_git/ms_est'!")
} else {
  source(file.path(prj, "src", "load.R"))
}

### ----------------------------------------------------------------------------
### Code to explore the pollution situation

# Load data ---------------------------------------------------------------
psm_sites <- fread(file.path(cachedir, 'psm_sites.csv'))
psm_samples <- fread(file.path(cachedir, 'psm_samples.csv'))
psm_sites_info <- fread(file.path(cachedir, 'psm_sites_info.csv'))
psm_variables <- fread(file.path(cachedir, 'psm_variables.csv'))
psm_maxtu <- fread(file.path(cachedir, 'psm_maxtu.csv'))
var_props <- fread(file.path(cachedir, 'var_props.csv'))


# restrict to sites < 100km² (or unknown)  and both data available (see do_overview.R)
take_site_id <- readRDS(file = file.path(cachedir, 'take_site_id.rds'))
psm_sites_info <- psm_sites_info[site_id %in% take_site_id]
psm_sites <- psm_sites[site_id %in% take_site_id]
psm_samples <- psm_samples[site_id %in% take_site_id]


## n neonics
psm_variables[name %chin% c('Thiacloprid', 'Imidacloprid', 'Clothianidin')]
neo_id <- psm_variables[name %chin% c('Thiacloprid', 'Imidacloprid', 'Clothianidin'), variable_id]
psm_samples[variable_id %in% neo_id, length(unique(sample_id)), by = variable_id]

# n measurements with RAC
psm_samples[variable_id %in% psm_variables[!is.na(rak_uba), variable_id]]


# Classify stream types ---------------------------------------------------

# # Analyze only small streams (catchment < 25km2)
# take_sites <- psm_sites_info[ezg_fin < 25, site_id]  
# take_samples <- psm_samples[site_id %chin% take_sites]

take_samples <- psm_samples

# # no. observations
# (tot <- nrow(take_samples))
# # no. sites
# length(unique(take_samples$site_id))
# # no. samples
# length(unique(take_samples$sample_id))
# # no. compounds
# length(unique(take_samples$variable_id))
# # values > loq
# length(take_samples[value_fin > 0, value_fin])
# 
# length(take_samples[variable_id == 457, value_fin])


# Risk Quotients ----------------------------------------------------------
# restrict samples to variables with rac
rac <- psm_variables[!is.na(rak_uba), list(variable_id, name, cas, pgroup, rak_uba)]
# join rac with samples
setkey(take_samples, variable_id)
setkey(rac, variable_id)
take_rac <- take_samples[rac]

# calculate Risk Quotient
take_rac[ , rq := value_fin / rak_uba]

# join variable name and type
setkey(take_rac, variable_id)
setkey(psm_variables, variable_id)
meas <- psm_variables[ , list(variable_id, name, psm_type)][take_rac]



# # remove variables with less then 1000 samples
# meas <- meas[variable_id %in% meas[ , list(n = length(value)), 
#                                   by = variable_id][n > 1000, variable_id]]
# 



# Compare with Stehle & Schulz ---------------------------------------
# restrict to insecticides compunds as Stehle & Schulz
rsvar <- c(290, # only alpha endosulfan, rs did no differenciate
           146, 213, 257, 397, 487, 488, 192, 1216, 167, 650,
           # beta cyfluthrin no in our datababse
           237, 872, 648, 645,
           # fenvalerate not in db
           905, 1195, 728, 365, 588, 589
           )
# compounds in common
com <- psm_variables[variable_id %in% rsvar & !is.na(rak_uba), list(variable_id, name, rak_uba)]
# add rac_sw from stehle
com$rs <- c(0.005, 0.1, 0.025, 0.3, 1.57, 2.8, 0.5)
# differences in RAC between studies
diff <- com$rak_uba - com$rs
mean(diff)
range(diff)

# restrict to only insecticides
insvar <- psm_variables[pgroup == 'organics, psm, insecticide', variable_id]
psm_variables[variable_id %in% insvar]
# meas_rvar <- meas[variable_id %in% rsvar]
meas_ins <- meas[variable_id %in% insvar]

# % of exceedances as proportion > loq
# nrow(meas_rvar[rq > 1]) / nrow(meas_rvar[rq > 0])
nrow(meas_ins[rq > 0])
length(unique(meas_ins$variable_id))
nrow(meas_ins[rq > 1]) / nrow(meas_ins[rq > 0]) * 100




# Maximum RQ --------------------------------------------------------------
# calculate maxRQ per variable and take only first 15 hits
take_var <- meas[ , list(media = max(rq)) , 
                  by = variable_id][order(media, 
                                          decreasing = TRUE)][1:15, variable_id]
take_meas <- meas[variable_id %in% take_var]

# highest RQs (15 subs, all samples)
take_meas[ , list(media = max(rq)) , by = name][order(media, decreasing = TRUE)]



## Exceedance as % of samples
# for each sample, take maximum RQ (from different compounds) and check exceedance
nrow(meas[ ,list(rq_max = max(rq)) , 
           by = sample_id][rq_max > 1]) /  length(unique(meas$sample_id)) * 100
# in 7.3\% of samples taken there was was RQ > 1 found
nrow(meas[ ,list(rq_max = max(rq)) , 
           by = sample_id][rq_max > 1]) /  length(unique(meas$sample_id[meas$rq > 0])) * 100
# in 14\% of samples with detecs was was RQ > 1 found


## Exceedances as % of measurements
mean(meas[ , rq > 1])*100 # 0.2 % of all measurements
nrow(meas[rq > 0])
mean(meas[rq > 0 , rq > 1]) * 100 # 5% of all detects

# % of detects
nrow(meas[rq > 1]) / nrow(meas[rq > 0]) * 100


## % number of sites
# number of SWB wit RAC exceedances
length(unique(meas[rq > 1, site_id]))
# 579 sites
length(unique(meas[, site_id]))
# form 2270
length(unique(meas[rq > 1, site_id])) / length(unique(meas[, site_id])) * 100
# = 25.5\%

length(unique(meas[rq > 1.12, site_id])) / length(unique(meas[, site_id])) * 100
# 24.8\% exceedance that are ecologically relevant (biodiversity 
#  reduciton of 30% at RQ= 1.12)

length(unique(meas[rq > 0.1, site_id])) / length(unique(meas[, site_id])) * 100
# 1/10 RQ == 11% reduction of biodiversity

# 23% of sites without detects
100 - length(unique(meas[rq  > 0, site_id])) / length(unique(meas[, site_id])) * 100



#  ------------------------------------------------------------------------

#show cummulative distribution of maxRQ
rqs <- c(0.1, 1, 1.12, logspace(log10(0.001), log10(100), length.out = 100))
prec <- numeric(length(rqs))
for (i in seq_along(rqs)) {
  message(i, '\n')
  prec[i] <- length(unique(meas[rq >= rqs[i], site_id])) 
}
prec <- prec / length(unique(meas[, site_id])) * 100

pdf(file.path(prj, "supplement/prac_ex.pdf"))
  plot(rqs[3:length(rqs)], prec[3:length(prec)], 
       log = 'x', 
       pch = 16,
       ylim = c(0, 100), 
       cex = 0.8,
       ylab = 'Fraction of sites',
       xlab = 'max(RQ)',
       xaxt = 'n') 
  axis(side = 1, at = c(0.001, 0.01, 0.1, 1,  10, 100))
  abline(v = 1, lty = 'dotted')
  abline(v = 0.1, lty = 'dotted')
  abline(h = prec[1], lty = 'dotted')
  abline(h = prec[2], lty = 'dotted')
dev.off()


# number of SWB sites with RAC 
length(unique(meas$site_id))   # 2270 (from 2301 in total)
# number of SWB samples
length(unique(meas$sample_id)) # 24344 samples with rac (from 24743 in total)



#  ------------------------------------------------------------------------
### table for comppunds with more the 1000 measurements
rac_dat <- psm_variables[ , list(variable_id, psm_type, name)][  # join with variable_names
  meas[ ,list(n_meas = .N,      # number of measurements
            n_detects = sum(value_fin > 0),    # number of detects
            p_detects = round(sum(value_fin > 0) / .N * 100, 1), # proportion of detects
            n_racex = sum(value_fin > rak_uba),  # number of exceedances
            p_racex = round(sum(value_fin > rak_uba) / .N * 100, 1), # prop of exceedances
            p_racex_d = round(sum(value_fin > rak_uba) / sum(value_fin > 0) * 100, 1) # prop of exceedances of detects
            ) , by = variable_id]
  ][n_meas > 1000]
rac_dat[ , variable_id := NULL]
rac_dat[ , psm_type := NULL]
rac_dat <- rac_dat[order(name)]
setnames(rac_dat, c('Name', 'No. ', 'No. \\textgreater LOQ', 
                    '\\% \\textgreater LOQ', 'No. RQ \\textgreater 1', 
                    '\\% RQ \\textgreater 1', '\\% RQ \\textgreater 1 | \\textgreater LOQ'))

rac_dat_x <- xtable(rac_dat, 
               label = 'tab:rac_dat',
               caption = c('Overview on RAC exceedances of the 78 compounds with more than 1000 measurements. No. = number of measurements;  \\% RQ \\textgreater 1 = RAC exceedances; \\% RQ \\textgreater 1 | \\textgreater LOQ= RAC exceedances as fraction of detects.',
                           'Overview on RAC exceedances of the 78 compounds with more than 1000 measurements.'),
               align = 'lp{3cm}rR{1.4cm}R{1.4cm}R{1.4cm}R{1.4cm}R{1.4cm}',
               digits = 1
)
print(rac_dat_x,
      file = file.path(prj, 'supplement/racdat.tex'),
      tabular.environment = "longtable",
      floating = FALSE,
      caption.placement = 'top',
      include.rownames = FALSE,
      comment = FALSE,
      booktabs = TRUE,
      hline.after = c(-1, 0, nrow(rac_dat)),
      sanitize.text.function = identity
)




#  ------------------------------------------------------------------------
# calculate percentage of non-detects
loqd <- take_meas[ , list(tot = length(rq),
                n_n0 = sum(rq > 0), # number of > LOQ
                n_0 = sum(rq == 0), # number of < LOQ
                p_n0 = round(sum(rq > 0) / length(rq) * 100, 1), # % 
                p_0 = round(sum(rq == 0) / length(rq) * 100, 1)), #%
                by = list(variable_id, name)][tot > 1000]
# same order as in plot
levs <- levels(reorder(take_meas[rq > 0, name], take_meas[rq > 0, rq], median))
loqd[ , name := factor(name, levels = levs[levs %in% loqd$name])]
loqd

prac <- ggplot() +
  geom_violin(data = take_meas[rq > 0 & name %in% loqd$name],
              aes(x = reorder(name, rq, FUN = median), y = rq, fill = psm_type)) +
  geom_hline(yintercept = 1, linetype = 'dotted') +
  coord_flip() +
  scale_y_log10(breaks = c(0.001, 0.01, 0.1, 1, 10, 100), 
                labels = c(0.001, 0.01, 0.1, 1, 10, 100),
                expand = c(0.05, 0.1)) +
  labs(x = 'Compound', y = 'Risk Quotient') +
  scale_fill_manual(name = 'Group',
                      values = c('#377eb8', '#4daf4a', '#B31010'),
                      labels = c('Fungicide', 'Herbicide', 'Insecticide')) +
  mytheme +
  theme(legend.position = 'bottom') +
  geom_text(data = loqd, aes(x = name, y = max(take_meas[rq > 0, rq]) + 2500, 
                             label = paste0(p_n0, '% (n=', tot, ')')),
            hjust = 'right', vjust = -0.3) 
# prac
ggsave(file.path(prj, "figure6.pdf"), prac, width = 7, height = 6.5)


# RQ exceedances for other compounds
take_samples[value_fin > 0, length(value_fin), by = variable_id][order(V1,decreasing = TRUE)]

# exceedances for Nicosulfuron, Diflufenican and Dimoxystrobin
psm_variables[name %like% c('Nicosu') | name %like% c('Diflufe') | name %like% c('Dimox')]
nrow(take_meas[variable_id %in% c(457) & rq > 1]) / nrow( take_meas[variable_id %in% c(457) & rq > 0]) * 100
nrow(take_meas[variable_id %in% c(272) & rq > 1]) / nrow( take_meas[variable_id %in% c(272) & rq > 0]) * 100
nrow(take_meas[variable_id %in% c(282) & rq > 1]) / nrow( take_meas[variable_id %in% c(282) & rq > 0]) * 100



# number of samples for Thiaclorpid
nrow(take_samples[variable_id == 588])



#  ------------------------------------------------------------------------
# detects
dt <- take_samples[ , list(tot = length(value_fin),
                     n_0 = sum(value_fin == 0),
                     n_d = sum(value_fin > 0),
                     p_0 = sum(value_fin == 0) / length(value_fin),
                     p_d = sum(value_fin > 0) / length(value_fin)), 
              by = variable_id]
setkey(dt, variable_id)
pdt <- psm_variables[ , list(variable_id, psm_type, name)][dt][order(p_d, decreasing = TRUE)]

# numbers for detection rates
pdt[order(p_d, decreasing = TRUE)][tot > 100]
pdt[order(p_d, decreasing = TRUE)][tot > 100 & psm_type != 'metabolite']


pdetects <- ggplot(pdt[p_d > 0.15 & tot > 100], aes(x = reorder(name, p_d), y = p_d * 100, col = psm_type, 
                            size = tot)) +
  geom_point() +
  coord_flip() + 
  mytheme +
  scale_color_manual(name = 'Group',
                    values = c('#377eb8', '#4daf4a', 'mediumorchid4'),
                    labels = c('Fungicide', 'Herbicide', 'Metabolite')) +
  scale_size_area(name = '# samples',
                  breaks = c(1000, 5000, 9000)) +
  labs(x = '', y = '% detects ')
# pdetects
ggsave(file.path(prj, "supplement/pdetects.pdf"), pdetects, width = 8, height = 6.5)




# Mixtures ----------------------------------------------------------------

# calculate the number of compounds per sample
mix <- take_samples[ , list(no_subs = sum(value_fin > 0)), by = list(sample_id, site_id, date)]

(tot <- nrow(mix))
# no detects
mix[no_subs == 0] # 3147 
nrow(mix[no_subs == 0]) / tot*100
# one compound
mix[no_subs == 1] # 1782
nrow(mix[no_subs == 1]) / tot * 100

# 2+ compound
mix[no_subs  >= 2] 
nrow(mix[no_subs >= 2]) / tot * 100

# max mixtures
max(mix$no_subs)


pmix <- ggplot(mix, aes(x = no_subs)) +
  geom_histogram(fill = 'grey50') +
  mytheme +
  labs(y = 'No. samples', x = 'No. compounds')
# pmix
ggsave(file.path(prj, "supplement/pmix.pdf"), pmix, width = 7, height = 6.5)

