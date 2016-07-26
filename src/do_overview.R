if (!exists("prj")) {
  stop("You need to create a object 'prj' that points to the top folder, 
       e.g. prj <- '/home/edisz/Documents/Uni/Projects/PHD/4BFG/Paper/ms_est' or
       prj <- '/home/user/Documents/projects_git/ms_est'!")
} else {
  source(file.path(prj, "src", "load.R"))
}

### ----------------------------------------------------------------------------
### Code for Overview on compiled dataset 


# Load data ---------------------------------------------------------------
psm_sites <- fread(file.path(cachedir, 'psm_sites.csv'))
psm_samples <- fread(file.path(cachedir, 'psm_samples.csv'))
psm_sites_info <- fread(file.path(cachedir, 'psm_sites_info.csv'))
psm_variables <- fread(file.path(cachedir, 'psm_variables.csv'))
psm_maxtu <- fread(file.path(cachedir, 'psm_maxtu.csv'))



# Spatial distribution ----------------------------------------------------
# (=Figure 1)
coordinates(psm_sites) <- ~easting + northing
proj4string(psm_sites) <- CRS("+init=epsg:31467")
psm_sites <- spTransform(psm_sites, CRS('+init=epsg:4326'))
psm_sites_bb <- psm_sites@bbox
psm_sites <- cbind(coordinates(psm_sites), psm_sites@data)
psm_sites$state_ab <- gsub('(.*?)_.*', '\\1', psm_sites$site_id)


library(raster)
adm1 <- raster::getData('GADM', country = 'DE', level = 1)
adm1 <- fortify(adm1)
p <- ggplot() +
  geom_polygon(data = adm1, aes(x = long, y = lat, group = group), fill = "grey90") +
  geom_path(data = adm1, aes(x = long, y = lat, group = group), size = .3) +
  geom_point(data = psm_sites, aes(x = easting, y = northing, col = state_ab), size = 1) +
  theme(legend.key = element_rect(fill = 'white')) +
  guides(colour = FALSE) +
  labs(x = 'Lon.', y = 'Lat.') +
  scale_color_hue(l = 50) +
  theme_bw() +
  coord_equal()
p
# ggsave(file.path(prj, "/fig/fig1.svg"),
#        p, width = 7, height = 7)




# Tabular overview --------------------------------------------------------
#  (= Tab. xx in Supplement)
psm_samples[ , state := gsub('(.*?)_.*', '\\1', site_id)]
psm_samples_tab <- psm_samples[, list(begin = min(date),
                                      end = max(date),
                                      no.sites = length(unique(site_id)),
                                      no.samples = length(unique(sample_id)),
                                      no.compounds = length(unique(variable_id))
),
by = state]
tot <- c('Total',
         as.character(min(psm_samples_tab$begin)),
         as.character(max(psm_samples_tab$end)),
         sum(psm_samples_tab$no.sites), 
         sum(psm_samples_tab$no.samples),
         length(unique(psm_samples$variable_id)))
psm_samples_tab <- rbind(data.frame(psm_samples_tab), tot)
psm_samples_tab <- psm_samples_tab[order(psm_samples_tab$state), ]
names(psm_samples_tab)[1] <- 'state \\textsuperscript{a}'
names(psm_samples_tab)[6] <- 'no.compounds \\textsuperscript{b}'

require(xtable)
xtab <- xtable(psm_samples_tab, 
       label = 'tab:phch_overview',
      caption = 'Overview on chemical samples. Only data from running waters and grab
sampling is shown. \\textsuperscript{a}: Abbreviations according to IS 3166-2:DE. \\textsuperscript{b}: Including metabolites',
      align = 'llllrrr'
)
print(xtab, 
      file = file.path(prj, 'supplement/phchoverview.tex'),
      caption.placement = 'top',
      include.rownames = FALSE,
      comment = FALSE,
      booktabs = TRUE,
      hline.after = c(-1, 0, 12, 13),
      sanitize.text.function = identity
      )


var_tab <- psm_variables[ , list(name, cas, pgroup, wrrl_zhkuqn, rak_uba, ger_auth_2015, eu_auth_2015)]
var_tab[ , ger_auth_2015 := as.character(ger_auth_2015)]
var_tab[ger_auth_2015 == 'FALSE', ger_auth_2015 := NA]
var_tab[ger_auth_2015 == 'TRUE', ger_auth_2015 := 'x']
var_tab[ , eu_auth_2015 := as.character(eu_auth_2015)]
var_tab[eu_auth_2015 == 'FALSE', eu_auth_2015 := NA]
var_tab[eu_auth_2015 == 'TRUE', eu_auth_2015 := 'x']
names(var_tab) <- c('Name', 'CAS', 'Group', 'MAC-EQS\\textsuperscript{a}', 'RAC \\textsuperscript{b}', 
                    'Auth. GER\\textsuperscript{c}', 'Auth. EU\\textsuperscript{d}')
var_tab$Group <- gsub('organics, psm, ', '', var_tab$Group)
# var_tab$Name <- gsub('ä', 'ae', var_tab$Name)
# fix bug with encoding (extra space or so....)
var_tab[Name == 'Benzoesäure', CAS := '65-85-0']

var_tab_x <- xtable(var_tab, 
                    label = 'tab:phch_var',
                    caption = 'Analysed chemical compounds. 
                    \\textsuperscript{a} Maximum Anual Concentration Environmental Quality Standard [ug/L].
                    \\textsuperscript{b} Regulatory Acceptable Concentration [ug/L] (Source: German EPA).
                    \\textsuperscript{c} Authorized in Germany (Source: BVL, 2015). 
                    \\textsuperscript{d} Authorized in the EU (Source: EU).',
                    align = 'lp{4cm}rlp{1cm}p{1.5cm}p{1.5cm}p{1cm}')
                    
print(var_tab_x, 
      file = file.path(prj, 'supplement/phchvar.tex'),
      tabular.environment="longtable",
      floating = FALSE,
      caption.placement = 'top',
      comment = FALSE,
      booktabs = TRUE,
      hline.after = c(-1, 0),
      sanitize.text.function = identity,
      size="\\fontsize{8pt}{10pt}\\selectfont"
)   




# Catchment size and landuse distribution ---------------------------------
options(stringsAsFactors = TRUE) # to fix bug in stat_density2d with polygons
ezg_lu <- ggplot(psm_sites_info[ezg_fin < 150 & !is.na(agri_fin) & !is.na(ezg_fin)], 
                 aes(x = ezg_fin, y = agri_fin * 100)) +
  stat_density2d(aes(alpha = ..level.., fill = ..level..), geom = "polygon") +
  geom_point(size = 0.5) +
  guides(alpha = FALSE, fill = FALSE) +
  mytheme +
  labs(x = expression('Catchment area ['~km^2~']'), y = expression('Agriculture [%]')) +
  scale_fill_gradient(low = "yellow", high = "red") +
  scale_x_continuous(breaks = c(0, 10, 25, 50, 100, 150))
ezg_lu
ezg_lu <- ggMarginal(ezg_lu, type = 'histogram', binwidth = 5)
ezg_lu
ggsave(file.path(prj, 'fig/ezg_lu.svg'), 
       width = 9, height = 7,
       device = grDevices::svg,
       ezg_lu)
options(stringsAsFactors = FALSE)




# Intersect samples with precipitation data -------------------------------

if (run_precip) {
  # path to regnie data
  # regpath <- '/home/edisz/Documents/Uni/Projects/PHD/4BFG/Project/data/regnie/'
  regpath <- '/home/user/Documents/projects_git/ms_est/data/regnie'
 
   # unique samplings
  samps <- unique(psm_samples[ , list(site_id, sample_id, date)])
  # join samps with coordinates
  setkey(samps, "site_id")
  setkey(psm_sites, "site_id")
  samps <- samps[psm_sites[ , list(site_id, easting, northing)]]
  
  # dates to search
  dates <- unique(samps[, date])
  
  
  container_dates <-  vector("list", length(dates)) 
  for (i in seq_along(dates)) {
    message('Processing date', i, ' out of ', length(dates))
    take <- dates[i]
    
    # read regnie data
    r <- try(esmisc::read_regnie(file.path(regpath, file_regnie(take))))
    if (inherits(r, 'try-error'))
      next
    # r <- projectRaster(from = r, crs = CRS('+init=epsg:31467'))
    # projectRaster() is slow - use gdalwarp instear
    tf <- tempfile(fileext = '.tif')
    tf2 <- tempfile(fileext = '.tif')
    writeRaster(r, tf)
    system(command = paste("gdalwarp -t_srs \'+init=epsg:31467\' -r near -overwrite", 
                           tf,
                           tf2))
    r <- raster(tf2)
    
    # get samplings at this date
    sl <- unique(samps[date == take, list(site_id, easting, northing)])
    # make spatial
    coordinates(sl) <- ~easting + northing
    proj4string(sl) <- CRS('+init=epsg:31467')
    # extract values from raster
    vals <- data.table(sl@data, val = extract(r, sl))
    # merge back and return data.frame
    out <- samps[date == take, list(sample_id, site_id, date)]
    setkey(out, site_id)
    setkey(vals, site_id)
    out <- out[vals]
    setkey(out, NULL)
    container_dates[[i]] <- out
  }
  
  # combine data from different dates
  precip_dates <- do.call(rbind, container_dates)
  saveRDS(precip_dates, file.path(cachedir, 'precip_dates.rds'))
  
  
  # precipitation at day before sampling
  dates1 <- (as.Date(dates) - 1) # day before sampling
  container_dates1 <-  vector("list", length(dates1)) 
  for (i in seq_along(dates1)) {
    message('Processing date', i, ' out of ', length(dates))
    # take = day before sampling
    take <- dates1[i]
    
    # read regnie data for day before
    r <- try(esmisc::read_regnie(file.path(regpath, file_regnie(take))))
    if (inherits(r, 'try-error'))
      next
    # r <- projectRaster(from = r, crs = CRS('+init=epsg:31467'))
    # projectRaster() is slow - use gdalwarp instear
    tf <- tempfile(fileext = '.tif')
    tf2 <- tempfile(fileext = '.tif')
    writeRaster(r, tf)
    system(command = paste("gdalwarp -t_srs \'+init=epsg:31467\' -r near -overwrite", 
                           tf,
                           tf2))
    r <- raster(tf2)
    
    # get samplings (=take + one date)
    sl <- unique(samps[date == take + 1, list(site_id, easting, northing)])
    # make spatial
    coordinates(sl) <- ~easting + northing
    proj4string(sl) <- CRS('+init=epsg:31467')
    # extract values from raster
    vals <- data.table(sl@data, val = extract(r, sl))
    # merge back and return data.frame
    out <- samps[date == take + 1, list(sample_id, site_id, date)]
    setkey(out, site_id)
    setkey(vals, site_id)
    out <- out[vals]
    setkey(out, NULL)
    container_dates1[[i]] <- out
  }
  precip_dates1 <- do.call(rbind, container_dates1)
  saveRDS(precip_dates1, file.path(cachedir, 'precip_dates1.rds'))
} else {
  precip_dates <- readRDS(file.path(cachedir, 'precip_dates.rds'))
  precip_dates1 <- readRDS(file.path(cachedir, 'precip_dates1.rds'))
}

precip_dates <- precip_dates[!is.na(val)]
precip_dates1 <- precip_dates1[!is.na(val)]
precp <- ggplot() +
  geom_histogram(data = precip_dates[val < 10], aes(x = val), fill = 'grey30', 
                 col = 'grey50', breaks = seq(0,50,2)) +
  geom_histogram(data = precip_dates[val >= 10], aes(x = val), fill = 'grey60', 
                 col = 'grey50', breaks = seq(0,50,2)) +
  mytheme +
  labs(x = 'Daily Precipitation [mm]', y = 'No. of samples') +
  ggtitle('At day of sampling (n = 41428)') +
  annotate('text', x = 4, y = 10000,
           label =  paste0(nrow(precip_dates[val < 10]), ' - ', 
                           round(nrow(precip_dates[val < 10]) / 
                                   nrow(precip_dates) * 100, 2), '%'
                           ), 
           hjust = 0.1, size = 5) +
  annotate('text', x = 20, y = 1000,
           label = paste0(nrow(precip_dates[val >= 10]), ' - ', 
                          round(nrow(precip_dates[val >= 10]) / 
                                  nrow(precip_dates) * 100, 2), '%'
                          ),
           hjust = 0, size = 5, col = 'grey30')
precp

precp1 <- ggplot() +
  geom_histogram(data = precip_dates1[val < 10], aes(x = val), fill = 'grey30', 
                 col = 'grey50', breaks = seq(0,50,2)) +
  geom_histogram(data = precip_dates1[val >= 10], aes(x = val), fill = 'grey60', 
                 col = 'grey50', breaks = seq(0,50,2)) +
  mytheme +
  labs(x = 'Daily Precipitation [mm]', y = 'No. of samples') +
  ggtitle('At day before of sampling (n = 42009)') +
  annotate('text', x = 4, y = 10000,
           label =  paste0(nrow(precip_dates1[val < 10]), ' - ', 
                           round(nrow(precip_dates1[val < 10]) / 
                                   nrow(precip_dates1) * 100, 2), '%'
           ), 
           hjust = 0.1, size = 5) +
  annotate('text', x = 20, y = 1000,
           label = paste0(nrow(precip_dates1[val >= 10]), ' - ', 
                          round(nrow(precip_dates1[val >= 10]) / 
                                  nrow(precip_dates1) * 100, 2), '%'
           ),
           hjust = 0, size = 5, col = 'grey30')
precp1
pp <- arrangeGrob(precp, precp1, nrow = 2)

ggsave(file.path(prj, 'supplement/precip.pdf'), 
       plot = pp, height = 12, width = 7)

