if (!exists("prj")) {
  stop("You need to create a object 'prj' that points to the top folder, 
       e.g. prj <- '/home/edisz/Documents/Uni/Projects/PHD/4BFG/Paper/ms_est' or
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
dd <- psm_variables[ , list(variable_id, name, psm_type)][take_rac]

# remove variable with less then 1000 samples
dd <- dd[variable_id %in% dd[ , list(n = length(value)), by = variable_id][n > 1000, variable_id]]

# calculate max per variable and take only first 15 hits
take_var <- dd[ , list(media = max(rq)) , by = variable_id][order(media, decreasing = TRUE)][1:15, variable_id]
take_dd <- dd[variable_id %in% take_var]

# highest RQs
take_dd[ , list(media = max(rq)) , by = name][order(media, decreasing = TRUE)][1:15]

# percentage of RQ > 1
nrow(dd[rq > 1]) / nrow(dd) * 100

# take max per sample
rak_samples <- length(unique(dd$sample_id))
# percentake of samples with RQ > 1In 
nrow(dd[ ,list(rq_max = max(rq)) , by = sample_id][rq_max > 1]) / rak_samples * 100

# number of SWB sites
length(unique(dd$site_id))
# number of SWB samples
length(unique(dd$sample_id))
# number of RAC exceedances
nrow(dd[rq >1]) / nrow(dd) * 100
nrow(dd[rq >1]) / nrow(dd[rq>0]) * 100

# number of SW wit RAC exceedances
length(unique(dd[rq > 1, site_id]))
length(unique(dd[, site_id]))
length(unique(dd[rq > 1, site_id])) / length(unique(dd[, site_id])) * 100


length(unique(dd[rq > 1.12, site_id])) / length(unique(dd[, site_id]))* 100



# calculate percentage of non-detects
loqd <- take_dd[ , list(tot = length(rq),
                n_n0 = sum(rq > 0),
                n_0 = sum(rq == 0),
                p_n0 = round(sum(rq > 0) / length(rq) * 100, 1),
                p_0 = round(sum(rq == 0) / length(rq) * 100, 1)), by = list(variable_id, name)]
# same order as in plot
levs <- levels(reorder(take_dd[rq > 0, name], take_dd[rq > 0, rq], median))
loqd[ , name := factor(name, levels = levs)]
loqd

prac <- ggplot() +
  geom_violin(data = take_dd[rq > 0],
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
  geom_text(data = loqd, aes(x = name, y = max(take_dd[rq > 0, rq]) + 2500, 
                             label = paste0(p_n0, '% (n=', tot, ')')),
            hjust = 'right', vjust = -0.3) 
# prac
ggsave(file.path(prj, "figure6.pdf"), prac, width = 7, height = 6.5)


# RQ exceedances for others
take_samples[value_fin > 0, length(value_fin), by = variable_id]

psm_variables[name %like% c('Nicosu') | name %like% c('Diflufe') | name %like% c('Dimox')]
nrow(take_dd[variable_id %in% c(457) & rq > 1]) / nrow( take_dd[variable_id %in% c(457) & rq > 0]) * 100
nrow(take_dd[variable_id %in% c(272) & rq > 1]) / nrow( take_dd[variable_id %in% c(272) & rq > 0]) * 100
nrow(take_dd[variable_id %in% c(282) & rq > 1]) / nrow( take_dd[variable_id %in% c(282) & rq > 0]) * 100


# number of samples for Thiaclorpid
take_samples[variable_id == 588]


# detects
dt <- take_samples[ , list(tot = length(value_fin),
                     n_0 = sum(value_fin == 0),
                     n_d = sum(value_fin > 0),
                     p_0 = sum(value_fin == 0) / length(value_fin),
                     p_d = sum(value_fin > 0) / length(value_fin)), 
              by = variable_id]
setkey(dt, variable_id)
pdt <- psm_variables[ , list(variable_id, psm_type, name)][dt][order(p_d, decreasing = TRUE)]


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
