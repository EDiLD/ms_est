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


# Classify stream types ---------------------------------------------------

# Analyze only small agricultural streams (catchment < 30km2 and >25% agriculture)
take_sites <- psm_sites_info[ezg_fin < 30 & agri_fin > 0.25, site_id]  
take_samples <- psm_samples[site_id %chin% take_sites]

# no. observations
(tot <- nrow(take_samples))
# no. sites
length(unique(take_samples$site_id))
# no. compounds
length(unique(take_samples$variable_id))
# values > loq
length(take_samples[value_fin > 0, value_fin])

length(take_samples[variable_id == 457, value_fin])


# EQS as endpoint ---------------------------------------------------------
# restrict samples to variables with eqs
eqs <- psm_variables[!is.na(wrrl_zhkuqn), list(variable_id, name, cas, pgroup, wrrl_zhkuqn)]
# join eqs
setkey(take_samples, variable_id)
setkey(eqs, variable_id)
take_eqs <- take_samples[eqs]

# check if value is greater than EQS
take_eqs[ , g_wrrl_zhkuqn := value_fin > wrrl_zhkuqn]
# split up by compounds
deqs <- take_eqs[!is.na(g_wrrl_zhkuqn), list(abs = sum(g_wrrl_zhkuqn), 
                                               prop = sum(g_wrrl_zhkuqn) / length(g_wrrl_zhkuqn), 
                                               n = length(g_wrrl_zhkuqn)), 
                   by = list(variable_id)][order(prop, decreasing = TRUE)]
deqs

# join variable name and type
setkey(deqs, variable_id)
setkey(psm_variables, variable_id)
deqs <- psm_variables[ , list(variable_id, name, psm_type)][deqs][order(prop, decreasing = TRUE)]


# take only fist 10 compounds
deqs <- deqs[1:10]
deqs$name <- factor(deqs$name)
# # manually change type for irgarol
# dd[variable_id == 375, psm_type := 'biocide']

peqs <- ggplot() +
  geom_point(data = deqs,
             aes(x = prop * 100, y = reorder(name, prop), size = n, col = psm_type)) +
  scale_size_continuous(name = 'no.samples',
                        range = c(2, 6), breaks = c(min(deqs$n), max(deqs$n))) +
  labs(y = 'Compound', x = 'Proportion Samples > EQS [%]') +
  scale_color_manual(name = 'Group',
                     values = c('#4daf4a', '#B31010'),
                     labels = c('Herbicide', 'Insekticide')) +
  mytheme +
  theme(legend.key = element_blank(),
        legend.position = "bottom") +
  guides(colour = FALSE)




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
dd <- dd[variable_id %in% take_var]

prac <- ggplot(dd) +
  geom_violin(aes(x = reorder(name, rq, FUN = median), y = rq, fill = psm_type)) +
  geom_hline(yintercept = 1, linetype = 'dotted') +
  coord_flip() +
  scale_y_log10(breaks = c(0.001, 0.01, 0.1, 1, 10, 100), 
                labels = c(0.001, '', 0.1, 1, 10, 100)) +
  labs(x = 'Compound', y = 'Risk Quotient') +
  scale_fill_manual(name = 'Group',
                      values = c('#377eb8', '#4daf4a', '#B31010'),
                      labels = c('Fungizid', 'Herbizid', 'Insektizid')) +
  mytheme +
  theme(legend.position = 'bottom')




# Toxic Units -------------------------------------------------------------

# calculate withing R (for verification)
setkey(var_props, variable_id)
setkey(take_samples, variable_id)

take_lc50 <- var_props[ , list(variable_id, lc50_dm_fin)][take_samples]
# rm samples without lc50 data
take_lc50 <- take_lc50[!is.na(lc50_dm_fin)]

# calculate TU
take_lc50[ , tu := value_fin / lc50_dm_fin]

# calculate TUmax per sample
tumax <- take_lc50[ , list(tumax = max(tu)), by = sample_id]


ptu <- ggplot(data = tumax, aes(x = log10(tumax))) +
  geom_histogram(fill = 'grey75') +
  geom_rug() +
  mytheme +
  labs(x = expression(log[10]*'('~TU[max]~')'), y = 'No. samples')



# Mixtures ----------------------------------------------------------------

# calculate the number of compounds per sample
mix <- take_samples[ , list(no_subs = sum(value_fin > 0)), by = list(sample_id, site_id, date)]


pmix <- ggplot(mix, aes(x = no_subs)) +
  geom_histogram(fill = 'grey75') +
  mytheme +
  labs(y = 'No. samples', x = 'No. compounds')


pall <- plot_grid(peqs, ptu, prac, pmix, labels = c('A', 'C', 'B', 'D'), 
                  scale = c(1, 1, 1, 1))
ggsave(file.path(prj, "/fig/pall.svg"), pall, width = 11, height = 11)
