if (!exists("prj")) {
  stop("You need to create a object 'prj' that points to the top folder, 
       e.g. prj <- '/home/edisz/Documents/Uni/Projects/PHD/4BFG/Paper/ms_est' or
       prj <- '/home/user/Documents/projects_git/ms_est'!")
} else {
  source(file.path(prj, "src", "load.R"))
}

### ----------------------------------------------------------------------------
### Code to clean data from database


### ----------------------------------------------------------------------------
### Load data

# Load from server
drv <- dbDriver("PostgreSQL")
con <- dbConnect(drv, dbname = DBname, user = DBuser, host = DBhost, port = DBport)


# samples restricted to: newer than 2005 and pesticide
psm_samples <- dbGetQuery(con, "SELECT phch_samples.sample_id,
    phch_samples.site_id,
    phch_samples.date,
    phch_samples.variable_id,
    phch_samples.qualifier,
    phch_samples.value,
    phch_samples.unit,
    phch_samples.value_fin,
    phch_samples.loq,
    phch_samples.lod,
    phch_samples.sample_type,
    phch_samples.sample_duration,
    phch_samples.fraction,
    phch_samples.comment,
    phch_samples.source
   FROM phch.phch_samples,
    phch.phch_variables
  WHERE phch_variables.variable_id = phch_samples.variable_id 
                          AND date_part('year'::text, phch_samples.date) >= 2005::double precision 
                          AND phch_variables.psm = true")

psm_sites <- dbGetQuery(con, 'SELECT * FROM phch.phch_sites')
psm_sites_info <- dbGetQuery(con, 'SELECT * FROM phch.phch_sites_info')
psm_variables <- dbGetQuery(con, 'SELECT * FROM phch.phch_variables')
psm_maxtu <- dbGetQuery(con, 'SELECT * FROM views.maxtu')
var_props <- dbGetQuery(con, 'SELECT variable_id,
lc50_dm_fin, 
lc50_dm_fin_source
                        FROM phch.phch_variables_prop')

dbDisconnect(con)
dbUnloadDriver(drv)


# make data.table
setDT(psm_sites)
setDT(psm_samples)
setDT(psm_sites_info)
setDT(psm_variables)
setDT(psm_maxtu)
setDT(var_props)

psm_sites[ , geom := NULL]

# Clean data --------------------------------------------------------------

### filter data

## filter based on sample data



# sample type
table(unique(psm_samples[, list(sample_id, sample_type)])[ , sample_type], useNA = 'always')
23277 + 19624
100 - 4357 / (4357 + 23277 + 19624) * 100
#!~ 90% of samples are grabsamples.
#! NA assumed to be individual samples
length(unique(psm_samples[sample_type == 'composite', sample_id]))
length(unique(psm_samples[sample_type == 'composite', site_id]))

#! keep only grab samples
psm_samples <- psm_samples[sample_type == 'individual' | is.na(sample_type)]

# fraction
table(psm_samples$fraction, useNA = 'always')
#! assume NA as Gesamtprobe / unfiltr.
#! pesticides were all measured as whole fraction

# only psm-substances and no sum-paramters
keep_variables <- psm_variables[sum_parameter == '' & psm == TRUE, variable_id]
psm_samples <- psm_samples[variable_id %in% keep_variables]


# restrict sites to psm_samples
psm_sites <- psm_sites[site_id %chin% unique(psm_samples$site_id)]



## filter based on site data
# lakes
table(psm_sites$lake, useNA = 'always')
#! 58 sites form lakes (NA assumed as streams)
psm_sites <- psm_sites[is.na(lake) | lake == FALSE]
psm_samples <- psm_samples[site_id %in% unique(psm_sites$site_id)]


## restrict tables to samples
psm_sites <- psm_sites[site_id %chin% unique(psm_samples$site_id)]
psm_variables <- psm_variables[variable_id %in%  unique(psm_samples$variable_id)]
psm_sites_info <- psm_sites_info[site_id %chin% unique(psm_samples$site_id)]
psm_maxtu <- psm_maxtu[site_id %chin% unique(psm_samples$site_id)]



# Check data --------------------------------------------------------------

# check sites <-> samples
length(unique(psm_samples[ , site_id])) == length(unique(psm_sites[ , site_id]))
all(unique(psm_samples[ , site_id]) %in% unique(psm_sites[ , site_id]))
all(unique(psm_sites[ , site_id]) %in% unique(psm_samples[ , site_id]))


# check sites <-> sites_info
length(unique(psm_sites_info[ , site_id])) == length(unique(psm_sites[ , site_id]))
#! Info not for all sites available
length(unique(psm_sites[ , site_id])) - length(unique(psm_sites_info[ , site_id]))
#! for 265 missing
psm_sites[!psm_sites[ , site_id] %in% psm_sites_info[ , site_id]]
table(psm_sites[!psm_sites[ , site_id] %in% psm_sites_info[ , site_id], state])
#! These are the sites were none of the both algorithms worked & no data was supplied from the states


# check sites <-> maxtu
length(unique(psm_maxtu[ , site_id])) == length(unique(psm_sites[ , site_id]))
all(unique(psm_maxtu[ , site_id]) %in% unique(psm_sites[ , site_id]))
all(unique(psm_sites[ , site_id]) %in% unique(psm_maxtu[ , site_id]))

# check samples <-> variables
all(unique(psm_samples[ , variable_id]) %in% unique(psm_variables[ , variable_id]))


head(psm_variables)
psm_variables[ , rak_uba_comment := NULL]

# var_prop <-> variables
var_props <- var_props[variable_id %in% psm_variables$variable_id]



# Save to cache -----------------------------------------------------------

write.table(psm_sites, file.path(cachedir, 'psm_sites.csv'), 
            row.names = FALSE, sep = ';')
write.table(psm_samples, file.path(cachedir, 'psm_samples.csv'), 
            row.names = FALSE, sep = ';')
write.table(psm_sites_info, file.path(cachedir, 'psm_sites_info.csv'), 
            row.names = FALSE, sep = ';')
write.table(psm_variables, file.path(cachedir, 'psm_variables.csv'), 
            row.names = FALSE, sep = ';')
write.table(psm_maxtu, file.path(cachedir, 'psm_maxtu.csv'), 
            row.names = FALSE, sep = ';')
write.table(var_props, file.path(cachedir, 'var_props.csv'), 
            row.names = FALSE, sep = ';')

