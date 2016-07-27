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


nrow(psm_sites_info[!(is.na(ezg_fin) | is.na(agri_fin))])
# 2376 with both data
# 265 completely missing
# 408 with either ezg or agri info missing



# Model number of RAK exceedances per sites -------------------------------


# join sites with info
setkey(psm_sites, site_id)
setkey(psm_sites_info, site_id)
psm_sites_wi <- psm_sites[psm_sites_info[ , list(site_id, use, ezg_fin, agri_fin)]]


### Endpoint: Number of RAK execeedances
raks <- psm_variables[!is.na(rak_uba), list(variable_id, name, cas, pgroup, rak_uba)]

samples_raks <- psm_samples[variable_id %in% raks$variable]
setkey(samples_raks, variable_id)
setkey(raks, variable_id)
samples_raks <- raks[samples_raks]
samples_raks[ , g_rak_uba := value_fin > rak_uba]
samples_raks <- samples_raks[month(date) %in% c(5, 6, 7)]

rak_thrs_sites <- samples_raks[!is.na(g_rak_uba) , list(sum_g_rak_uba_abs = sum(g_rak_uba),
                                  sum_g_rak_uba_prop = sum(g_rak_uba) / length(g_rak_uba),
                                  n_samp = length(g_rak_uba)), by = site_id][order(sum_g_rak_uba_prop, decreasing = TRUE)]
setkey(rak_thrs_sites, site_id)
setkey(psm_sites, site_id)

w3 <- psm_sites_wi[rak_thrs_sites]
w3 <- w3[ezg_fin < 100 & ezg_fin > 1 & n_samp > 10]

ggplot(w3, aes(x = ezg_fin, y = sum_g_rak_uba_prop)) +
  geom_point() +
  geom_smooth()


ggplot(w3, aes(x = agri_fin, y = sum_g_rak_uba_prop)) +
  geom_point() +
  geom_smooth()
