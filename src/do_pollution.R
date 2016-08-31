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


## n neonics
psm_variables[name %chin% c('Thiacloprid', 'Imidacloprid', 'Clothianidin')]
neo_id <- psm_variables[name %chin% c('Thiacloprid', 'Imidacloprid', 'Clothianidin'), variable_id]
psm_samples[variable_id %in% neo_id, length(unique(sample_id)), by = variable_id]

# n measurements with RCA
psm_samples[variable_id %in% psm_variables[!is.na(rak_uba), variable_id]]


# Classify stream types ---------------------------------------------------

# Analyze only small agricultural streams (catchment < 30km2 and >25% agriculture)
take_sites <- psm_sites_info[ezg_fin < 30 & agri_fin > 0.25, site_id]  
take_samples <- psm_samples[site_id %chin% take_sites]

# no. observations
(tot <- nrow(take_samples))
# no. sites
length(unique(take_samples$site_id))
# no. samples
length(unique(take_samples$sample_id))
# no. compounds
length(unique(take_samples$variable_id))
# values > loq
length(take_samples[value_fin > 0, value_fin])

length(take_samples[variable_id == 457, value_fin])


# Risk Quotients ----------------------------------------------------------
# restrict samples to variables with rac
rac <- psm_variables[!is.na(rak_uba), list(variable_id, name, cas, pgroup, rak_uba)]
# join eqs
setkey(take_samples, variable_id)
setkey(rac, variable_id)
take_rac <- take_samples[value_fin > 0][rac]
# rm left over from left join
take_rac <- take_rac[!is.na(site_id)]

# calculate Risk Quotient
take_rac[ , rq := value_fin / rak_uba]

# join variable name and type
setkey(take_rac, variable_id)
setkey(psm_variables, variable_id)
dd <- psm_variables[ , list(variable_id, name, psm_type)][take_rac]

# calculate media per variable and take only first 15 hits
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




prac <- ggplot(take_dd) +
  geom_violin(aes(x = reorder(name, rq, FUN = median), y = rq, fill = psm_type)) +
  geom_hline(yintercept = 1, linetype = 'dotted') +
  coord_flip() +
  scale_y_log10(breaks = c(0.001, 0.01, 0.1, 1, 10, 100), 
                labels = c(0.001, '', 0.1, 1, 10, 100)) +
  labs(x = 'Compound', y = 'Risk Quotient') +
  scale_fill_manual(name = 'Group',
                      values = c('#377eb8', '#4daf4a', '#B31010'),
                      labels = c('Fungicide', 'Herbicide', 'Insecticide')) +
  mytheme +
  theme(legend.position = 'bottom')



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
  geom_histogram(fill = 'grey75') +
  mytheme +
  labs(y = 'No. samples', x = 'No. compounds')


pall <- plot_grid(peqs, ptu, prac, pmix, labels = c('A', 'B', 'C', 'D'), 
                  label_size = 20)
ggsave(file.path(prj, "figure5.pdf"), pall, width = 7, height = 6.5, 
       units = 'in', dpi = 300, scale = 1.5)
