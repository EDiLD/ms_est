if (!exists("prj")) {
  stop("You need to create a object 'prj' that points to the top folder, 
       e.g. prj <- '/home/edisz/Documents/work/research/projects/2016/4BFG/Paper/ms_est' or
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
var_props <- fread(file.path(cachedir, 'var_props.csv'))


# restrict to sites < 100km² (or unknown) and > 2005
take_sites <- psm_sites_info[ezg_fin <= 100 | is.na(ezg_fin), site_id]

psm_sites_info <- psm_sites_info[site_id %in% take_sites]
psm_sites <- psm_sites[site_id %in% take_sites]
psm_samples <- psm_samples[site_id %in% take_sites]



# overview on Catchment derivation ----------------------------------------
# total number of sites
nrow(psm_sites_info)

table(psm_sites_info$ezg_fin_source, useNA = 'always')
# NA = no information available
# auth = from authorities direct
# drain = from river segments
# TRUE = from DEM
table(psm_sites_info$ezg_fin_source, useNA = 'always') / nrow(psm_sites_info) * 100
# 1401 (=59%) with data from authorities
# 372 (=16%) with data from stream segments
# 568 (=24%) with data from dem
# 28 (=1%) without any information

# sites with data from auth
nrow(psm_sites_info[!is.na(ezg_auth)]) / nrow(psm_sites_info)
# sites with good delienation
nrow(psm_sites_info[!is.na(ezg_gis) & use == TRUE]) / nrow(psm_sites_info)
# use == TRUE: catchment from DEM good
# use == drain: catchment from DEM bad, but river segment good
# use == NA: no good shapes available
nrow(psm_sites_info[!is.na(ezg_gis) & (use == 'drain')]) / nrow(psm_sites_info) * 100



table(psm_sites_info$agri_fin_source, useNA = 'always')
table(psm_sites_info$agri_fin_source, useNA = 'always') / nrow(psm_sites_info) * 100
# 571 (=24%) from authorities
# 522 (=22%) from drain
# 1231 (=52%) from DEM
# 45 (=2%) missing



nrow(psm_sites_info[!(is.na(ezg_fin) | is.na(agri_fin))])

nrow(psm_sites_info[!(is.na(ezg_fin) | is.na(agri_fin))]) / nrow(psm_sites_info)
# 2301 (=97%) with both informations


# restrict to sites < 100km² (or unknown)  and both data available
take_site_id <- psm_sites_info[!(is.na(ezg_fin) | is.na(agri_fin)), site_id]
saveRDS(take_site_id, file = file.path(cachedir, 'take_site_id.rds'))

# removed sites
psm_sites_info[!site_id %in% take_site_id]
psm_sites_info <- psm_sites_info[site_id %in% take_site_id]
psm_sites <- psm_sites[site_id %in% take_site_id]
psm_samples <- psm_samples[site_id %in% take_site_id]




# some numbers 
# no. observations
(tot <- nrow(psm_samples))
# no.of samples
length(unique(psm_samples$sample_id))
# no. sites
length(unique(psm_samples$site_id))
# no. compounds
length(unique(psm_samples$variable_id))
# values > loq
length(psm_samples[value_fin > 0, value_fin])
length(psm_samples[value_fin > 0, value_fin]) / nrow(psm_samples) * 100 




# Temporal overview of samples --------------------------------------------
# For Supplement
monthly <- psm_samples[ , length(unique(sample_id)), by = month(date)]
ggplot(monthly, aes(x = factor(month), y = V1)) + 
  geom_bar(stat = 'identity')
yearlymonthly <- psm_samples[ , length(unique(sample_id)), by = list(year(date), month(date))]
p_temp <-  ggplot(yearlymonthly, aes(x = factor(month), y = V1)) + 
  geom_bar(stat = 'identity') +
  facet_grid(~year) +
  labs(x = 'Month', y = 'No. samples') 
# p
ggsave(file.path(prj, 'supplement/temporal.pdf'), p_temp, width = 20, height = 7)

# # phd <- '/home/edisz/Documents/work/research/projects/2016/1PHD/phd_thesis/appendix/smallstreams/one/'
# p_temp <-  ggplot(yearlymonthly, aes(x = factor(month), y = V1)) + 
#   geom_bar(stat = 'identity') +
#   facet_grid(~year) +
#   labs(x = 'Month', y = 'No. samples') +
#   phdtheme +
#   theme(panel.grid.major = element_blank()) +
#   scale_x_discrete(labels = c('2' = '',
#                               '3' = '',
#                               '5' = '',
#                               '6' = '',
#                               '8' = '',
#                               '9' = '',
#                               '11' = '',
#                               '12' = ''))
# ggsave(file.path(phd, 'temporal.pdf'), p_temp, width = 250, units = 'mm')



# Spatial distribution ----------------------------------------------------
# (=Figure 1)
psm_sites_p <- as.data.frame(psm_sites)
coordinates(psm_sites_p) <- ~easting + northing
proj4string(psm_sites_p) <- CRS("+init=epsg:31467")
psm_sites_p <- spTransform(psm_sites_p, CRS('+init=epsg:4326'))
psm_sites_bb <- psm_sites_p@bbox
psm_sites_p <- cbind(coordinates(psm_sites_p), psm_sites_p@data)
psm_sites_p$state_ab <- gsub('(.*?)_.*', '\\1', psm_sites_p$site_id)
setDT(psm_sites_p)

# get administrative borders
adm1 <- raster::getData('GADM', country = 'DE', level = 1)
adm1 <- fortify(adm1)
p_map <- ggplot() +
  geom_polygon(data = adm1, aes(x = long, y = lat, group = group), fill = "grey90") +
  geom_path(data = adm1, aes(x = long, y = lat, group = group), size = .3) +
  geom_point(data = psm_sites_p, aes(x = easting, y = northing, col = state_ab), 
             size = 1) +
  theme(legend.key = element_rect(fill = 'white')) +
  # guides(colour = FALSE) +
  labs(x = 'Lon.', y = 'Lat.') +
  scale_color_hue(name = 'state', l = 50) +
  mytheme +
  coord_equal()
# p_map
ggsave(file.path(prj, "figure1.pdf"),
       p_map, width = 3.5, height = 3, units = 'in', dpi = 300, scale = 2)

# phd <- '/home/edisz/Documents/work/research/projects/2016/1PHD/phd_thesis/chapters/smallstreams/'
# ggsave(file.path(phd, "figure1.pdf"),
#        p_map, width = 176, height = 150, units = 'mm')


# Tabular overview --------------------------------------------------------
#  Sites and samples (= Tab. xx in Supplement)
psm_samples[ , state := gsub('(.*?)_.*', '\\1', site_id)]
psm_samples_tab <- psm_samples[, list(begin = min(date),
                                      end = max(date),
                                      no.sites = length(unique(site_id)),
                                      no.samples = length(unique(sample_id)),
                                      no.compounds = length(unique(variable_id))),
                               by = state]
tot <- c('Total',
         as.character(min(psm_samples_tab$begin)),
         as.character(max(psm_samples_tab$end)),
         sum(psm_samples_tab$no.sites), 
         sum(psm_samples_tab$no.samples),
         length(unique(psm_samples$variable_id)))
psm_samples_tab <- rbind(data.frame(psm_samples_tab), tot)
psm_samples_tab <- psm_samples_tab[order(psm_samples_tab$state), ]
psm_samples_tab$name <- c('Baden-Württemberg', 'Bavaria', 'Hesse', 
                          'Mecklenburg-Western Pomerania', 'Lower Saxony',
                          'North Rhine-Westphalia', 'Rhineland-Palatinate',
                          'Schleswig-Holstein', 'Saarland', 'Saxony', 
                          'Saxony-Anhalt', 'Thuringia', ''
                          )
setcolorder(psm_samples_tab, c(7, 1:6))
setnames(psm_samples_tab, c('name', 'abbrv.\\textsuperscript{a}',
                            'Begin', 'End', 'No. sites', 'No.samples',
                            'No. pesticides\\textsuperscript{b}'))


xtab <- xtable(psm_samples_tab, 
       label = 'tab:phch_overview',
      caption = c('Overview on chemical samples. Only data from running waters and grab
sampling is shown. \\textsuperscript{a}: Abbreviations according to ISO 3166-2:DE. 
      \\textsuperscript{b}: Including metabolites',
                  'Overview on chemical samples.'),
      align = 'lp{2.7cm}lllR{2cm}R{2cm}R{2cm}'
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


# Tabular overview --------------------------------------------------------
#  Variables (= Tab. xx in Supplement)
var_tab <- psm_variables[ , list(variable_id, name, cas, pgroup, wrrl_zhkuqn, rak_uba, 
                                 ger_auth_2015, eu_auth_2015)]
# join with properties-data
setkey(var_tab, variable_id)
setkey(var_props, variable_id)
var_tab <- var_props[var_tab] 
# restrict to variable found in psm_samples
take_var_id <- unique(psm_samples$variable_id)
var_tab <- var_tab[variable_id %in% take_var_id, ]
var_tab[!is.na(rak_uba)]
# 107 RAC values available

# prettify
var_tab[ , variable_id := NULL]
setcolorder(var_tab, c(names(var_tab)[3:5], names(var_tab)[8:9], names(var_tab)[1:2], names(var_tab)[6:7]))
var_tab <- var_tab[ , c(1,2,3,4, 5, 9), with = FALSE]

var_tab[ , ger_auth_2015 := as.character(ger_auth_2015)]
var_tab[ger_auth_2015 == 'FALSE', ger_auth_2015 := NA]
var_tab[ger_auth_2015 == 'TRUE', ger_auth_2015 := 'x']
var_tab[ , eu_auth_2015 := as.character(eu_auth_2015)]
var_tab[eu_auth_2015 == 'FALSE', eu_auth_2015 := NA]
var_tab[eu_auth_2015 == 'TRUE', eu_auth_2015 := 'x']

names(var_tab) <- c('Name', 'CAS', 'Group', 
                    'Auth. GER\\textsuperscript{a}', 
                    'Auth. EU\\textsuperscript{b}',
                    'RAC \\textsuperscript{c}'
                    )
var_tab$Group <- gsub('organics, psm, ', '', var_tab$Group)

# fix bug with encoding (extra space or so....)
var_tab[Name == 'Benzoesäure', CAS := '65-85-0']

# order by name
var_tab <- var_tab[order(Name)]

# round to four digits


var_tab_x <- xtable(var_tab, 
                    label = 'tab:phch_var',
                    caption = c('Overview on pesticides (and metabolites) in the database. \\
                    \\textsuperscript{a} Authorized in Germany (Source: German Federal Office of Consumer Protection and Food Safety (BVL) as at March 2015). 
                    \\textsuperscript{b} Authorized in the European union (Source: EU Pesticides database as at March 2015).
                    \\textsuperscript{c} Regulatory Acceptable Concentration [$\\mu g/L$] (Source: German Environment Agency (UBA) as at November 2015).',
                    'Overview on pesticides in the database.'),
                    align = 'lp{4cm}rlp{1.3cm}p{1.3cm}p{1.5cm}',
                    digits = 5)
                    
print(var_tab_x,
      file = file.path(prj, 'supplement/phchvar.tex'),
      tabular.environment = "longtable",
      floating = FALSE,
      caption.placement = 'top',
      comment = FALSE,
      booktabs = TRUE,
      hline.after = c(-1, 0),
      sanitize.text.function = identity#,
      # size="\\fontsize{8pt}{10pt}\\selectfont"
)



# Compound spectra --------------------------------------------------------

psm_variables[!is.na(rak_uba)]
# 107 compunds with RAKS


### Measured Spectra
## total number of substances
length(unique(psm_samples[ , variable_id]))
# no. of compounds per group
sort(table(psm_variables$pgroup, useNA = 'always'), decreasing = TRUE)
# 179 herbicides
# 117 insecticides
# 109 fungicides

# values > loq
length(psm_samples[value_fin > 0, value_fin])
length(psm_samples[value_fin > 0, value_fin]) / nrow(psm_samples) * 100 
# 4% above LOQ


var_bl <- psm_samples[ , list(length(unique(variable_id))) , by = substr(site_id, 1, 2)]
var_bl[substr %chin% c('ST', 'SL', 'TH')]
var_bl[substr %chin% c('RP', 'NI')]
var_bl[!substr %chin% c('RP', 'NI', 'ST', 'SL', 'TH')]

# bring to wide format (BL in rows, chemical in columns)
# aggregate using 'length' (=irrespective of value, how often has it been measured)
vw <- dcast.data.table(psm_samples, substr(site_id, 1, 2) ~ variable_id)
# make binary data.frame from data.table
# relies on a side-effect
makeone <-  function(DT) {
  for (i in names(DT)[-c(1)])
    DT[get(i) > 0, i := 1, with = FALSE]
}
makeone(vw)

# calculte jaccard distance between BL
rownames(vw) <- vw[ ,site_id]
dp <- vegdist(vw[ , -1, with = FALSE], method = 'jaccard')

# hierarchical clustering
hc <- hclust(dp, method = 'complete')
plot(hc, labels = vw$site_id)

# # fusion level plot
plot(hc$height, nrow(vw[ , -1, with = FALSE]):2, type = 'S')
text(hc$height, nrow(vw[ , -1, with = FALSE]):2, nrow(vw[ , -1, with = FALSE]):2, col = 'red')
# not easy to spot

#  #average silhoutte with to choose be number of clusters
nr <- nrow(vw[ , -1, with = FALSE])

# number of clusters to check
ks <- 2:(nr - 1)
asw <- numeric(length(ks))
for (i in seq_along(ks)) {
  sil <- silhouette(cutree(hc, k = ks[i]), dp)
  asw[i] <- summary(sil)$avg.width
}
ks[which.max(asw)]

pdf(file = file.path(prj, 'supplement/silhouette.pdf'))
  plot(ks, asw, type = 'h', ylab = 'Average Silhoutte Width', 
       xlab = 'Number of clusters')
  segments(ks[which.max(asw)], 0, ks[which.max(asw)], max(asw), col = 'red')
  points(ks[which.max(asw)], max(asw), pch = 16, col = 'red')
dev.off()

# pdf(file = file.path(phd, 'silhouette.pdf'))
#   par(mar = c(5,5,4,2))
#   plot(ks, asw, type = 'h', ylab = 'Average Silhoutte Width', 
#        xlab = 'Number of clusters', cex.lab = 2, cex.axis = 2)
#   segments(ks[which.max(asw)], 0, ks[which.max(asw)], max(asw), col = 'red', lwd = 2)
#   points(ks[which.max(asw)], max(asw), pch = 16, col = 'red', cex = 1.5)
# dev.off()

sil <- silhouette(cutree(hc, k = 2), dp)
plot(sil)


# two groups best

# # or maximise the correlation
# grpdist <- function(x) {
#   gr <- data.frame(factor(x))
#   distgr <- daisy(gr, "gower")
#   distgr
# }
# kt <- data.frame(k = 2:(nr-1), r = NA)
# for(i in seq_len(nrow(kt))) {
#   gr <- cutree(hc, k = kt$k[i])
#   distgr <- grpdist(gr)
#   kt$r[i] <- cor(dp, distgr, method = 'pearson')
# }
# kt
# plot(kt)


# 3 distinct groups: NI+RP; ST+SL+TH+BW; Rest
bl_groups <- cutree(hc, k = 2)

## nicer plot
# colors (red, blue or green)
# colo <- c("#E41A1C", "#4DAF4A", "#377EB8")
colo <- c("#E41A1C", "#377EB8")
hcd <- as.dendrogram(hc)

# function to get color labels
color_label <- function(n) {
  if (is.leaf(n)) {
    a <- attributes(n)
    col <- colo[bl_groups[a$label]]
    attr(n, "label") <- vw$site_id[a$label]
    attr(n, "nodePar") <- c(a$nodePar, lab.col = col)
  }
  n
}
# using dendrapply
clus = dendrapply(as.dendrogram(hc), color_label)

# make plot
pdf(file = file.path(prj, 'supplement/varclus.pdf'))
  plot(clus, ylab = 'Jaccard Distance')
dev.off()

# pdf(file = file.path(phd, 'varclus.pdf'), width = 5, height = 5)
#   plot(clus, ylab = 'Jaccard Distance', cex.lab = 1.5, cex.axis = 1.3)
# dev.off()


# PCO (principal coordinates)
pco1 <- wcmdscale(dp, k = 2, eig = TRUE)
rownames(pco1$points) <- unlist(vw[ , 1, with = FALSE])
pco_dat <- data.frame(scores(pco1), state = rownames(scores(pco1)), group = bl_groups)
evar <- round(pco1$eig / sum(pco1$eig), 2) * 100
xlab <- paste0('Axis 1 (', evar[1], '%)')
ylab <- paste0('Axis 2 (', evar[2], '%)')

p_mds <- ggplot(pco_dat, aes(x = Dim1, y = Dim2)) +
  mytheme +
  scale_colour_manual(guide = FALSE, values = colo) +
  xlab(xlab) +
  ylab(ylab) +
  geom_vline(xintercept = 0, linetype = 'dotted') +
  geom_hline(yintercept = 0, linetype = 'dotted') +
  geom_text(aes(label = state, col = as.factor(group)), size = 6) 
# p_mds


# tile plot / barcode
# bring binary data.frame to long format
# add 12 white colors
cols <- c(colo[bl_groups], rep('white', 12))
vwm <- melt(vw, id.vars = 'site_id')
# chaneg for colors
vwm$value <- ifelse(vwm$value == 1, 'ja', 'nein')
vwm$cols <- factor(paste0(vwm$value, vwm$site_id))
vwm$site_id <- factor(vwm$site_id, levels = c('TH', 'ST', 'SN', 'SL', 'SH',
                                              'NW', 'MV', 'HE', 'BW', 'BY', 'RP', 'NI'))

# plot
p_bar <- ggplot(vwm, aes(x = variable, y = site_id, fill = cols)) +
  geom_tile() +
  mytheme +
  scale_fill_manual(values = cols, name = 'gemessen') +
  theme(axis.text.x = element_blank(),
        axis.ticks.x = element_blank()) +
  labs(x = 'Compound', y = 'State') +
  guides(fill = FALSE) 
# p_bar
p <- plot_grid(p_bar, p_mds, rel_heights = c(3, 1))
# plot(p)

ggsave(file.path(prj, "figure2.pdf"),
       p, width = 7, height = 3.5, units = 'in', dpi = 300, scale = 1.5)


ggsave(file.path(phd, "figure2.pdf"),
       p, width = 7, height = 3.5, units = 'in', dpi = 300, scale = 1.5)




# Catchment size and landuse distribution ---------------------------------
options(stringsAsFactors = TRUE) # to fix bug in stat_density2d with polygons
ezg_lu <- ggplot(psm_sites_info[ezg_fin < 150 & !is.na(agri_fin) & !is.na(ezg_fin)], 
                 aes(x = ezg_fin, y = agri_fin * 100)) +
  stat_density2d(aes(alpha = ..level.., fill = ..level..), geom = "polygon") +
  geom_point(size = 0.5) +
  guides(alpha = FALSE, fill = FALSE) +
  mytheme +
  # phdtheme + 
  labs(x = expression('Catchment area ['~km^2~']'), y = expression('Agriculture [%]')) +
  # scale_fill_gradient(low = "yellow", high = "red") +
  scale_fill_viridis() +
  scale_x_continuous(breaks = c(0, 10, 25, 50, 100))
# ezg_lu
ezg_lu <- ggMarginal(ezg_lu, type = 'histogram', binwidth = 5)
# ezg_lu
ggsave(file.path(prj, 'figure3.pdf'), 
       ezg_lu,
       width = 3.3, height = 3.3/1.2,
       units = 'in', dpi = 300, scale = 2.5)

# ggsave(file.path(phd, 'figure3.pdf'),
#        ezg_lu,
#        width = 3, height = 3/1.2,
#        units = 'in', dpi = 300, scale = 2.5)

options(stringsAsFactors = FALSE)




# Intersect samples with precipitation data -------------------------------
run_precip <- FALSE
# if false, use precomouted data from cache 
if (run_precip) {
  # path to regnie data
  regpath <- '/home/edisz/Documents/work/research/projects/2016/4BFG/Project/data/regnie/'
  # regpath <- '/home/user/Documents/projects_git/ms_est/data/regnie'
 
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

# retrict to actual samples (were computedt for all)
precip_dates <- precip_dates[!is.na(val)]
precip_dates1 <- precip_dates1[!is.na(val)]
precip_dates <- precip_dates[sample_id %chin% psm_samples$sample_id]
precip_dates1 <- precip_dates1[sample_id %chin% psm_samples$sample_id]
precip_dates1 <- precip_dates1[sample_id %chin% precip_dates$sample_id]


precp <- ggplot() +
  geom_histogram(data = precip_dates[val < 10], aes(x = val), fill = 'grey30', 
                 col = 'grey50', breaks = seq(0,50,2)) +
  geom_histogram(data = precip_dates[val >= 10], aes(x = val), fill = 'grey60', 
                 col = 'grey50', breaks = seq(0,50,2)) +
  mytheme +
  labs(x = 'Daily Precipitation [mm]', y = 'No. of samples') +
  ggtitle('At day of sampling') +
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
# precp

precp1 <- ggplot() +
  geom_histogram(data = precip_dates1[val < 10], aes(x = val), fill = 'grey30', 
                 col = 'grey50', breaks = seq(0,50,2)) +
  geom_histogram(data = precip_dates1[val >= 10], aes(x = val), fill = 'grey60', 
                 col = 'grey50', breaks = seq(0,50,2)) +
  mytheme +
  labs(x = 'Daily Precipitation [mm]', y = 'No. of samples') +
  ggtitle('At day before of sampling') +
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
# precp1
pp <- arrangeGrob(precp, precp1, nrow = 2)

ggsave(file.path(prj, 'supplement/precip.pdf'),
       plot = pp, height = 12, width = 7)


# ggsave(file.path(phd, 'precip.pdf'),
#        plot = pp, height = 200, units = 'mm')


