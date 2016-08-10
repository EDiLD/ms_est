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

# i) small, agricultural streams (catchment size < 30 km2 and  agriculture > 40%);
# ii) small, non-agricultural streams (catchment size < 30 km2 and agriculture < 10%);
# iii) medium sized, agricultural streams (30 km2 <= catchment size < 150 km2 and agriculture > 40%);
# iv) medium sized, non-agricultural streams (30 km2 <= catchment size < 150 km2 and agriculture < 10%).

# classify stream types
psm_sites_info[ezg_fin < 30 & agri_fin > 0.4, type := 'small, agricultural']  # small agricultural
psm_sites_info[ezg_fin < 30 & agri_fin < 0.1, type := 'small, non-agricultural'] # small non - agricultural
psm_sites_info[ezg_fin >= 30 & ezg_fin < 150 & agri_fin > 0.4, type := 'mid, agricultural'] # agricultural
psm_sites_info[ezg_fin >= 30 & ezg_fin < 150 & agri_fin < 0.1, type := 'mid, non-agricultural'] # non-agricultural

table(psm_sites_info$type, useNA = 'always') # 1040 sites ommitted.

