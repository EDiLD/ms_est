###-----------------------------------------------------------------------------
### Defined functions

###-----------------------------------------------------------------------------
### Checking functions -
#' Small function to return unique levels if there aren't to much.
#' @param x a vector
#' @param cut how many levels should be printed?
chk_uniq <- function(x, cut=10) {
  unx <- unique(x)
  if (length(unx) <= cut) {
    return(unx)
  } else {
    return('manylevels')
  }
}


#' function to check sites
#' @param x a data.table of class 'physico_chem'
chk_sites <- function(x) {
  if(!is.data.table(x))
    stop('sites must be a data.table!')
  # check unique Ids
  setkey(x, NULL)
  unid_OK <- length(unique(x[ , site_id])) == nrow(x)
  if(unid_OK){
    print('Unique ID: PASS')
  } else {
    warning('Non unique ID!')
  }
  # check classes
  colclas <- sapply(x, class) 
  colclas_OK <- all(colclas == c(rep('character', 6), 'numeric', 'integer', 'integer'))
  if(colclas_OK){
    print('colClasses: PASS')
  } else {
    warning('colClasses not OK!')
  }
  # check names
  colnames <- names(x)
  colnames_OK <- all(colnames == c("site_id", "state", "site_nr", "site_name", 
                                   "stream", "lake", "ezg", "easting", "northing"))
  if(colnames_OK){
    print('colnames: PASS')
  } else {
    warning('colnames do not fit!')
  }
  # check duplicates
  dups_OK <- sum(duplicated(x)) == 0
  if(dups_OK){
    print('duplicates values: PASS')
  } else {
    message('duplicates in data! : NOPASS')
  }
  
  if(all(unid_OK, colclas_OK, colnames_OK, dups_OK)){
    print('overall: PASS')
    return(TRUE)
  } else {
    warning('overall: NOPASS')
    return(FALSE)
  }
}

# function to check chem-data
#' @param x a data.table
chk_chem <- function(x, variables){
  #! add check of value_fin (if qualifier is < then 0)
  if (!is.data.table(x))
    stop('chem must be a data.table!')
  setkey(x, NULL)
  # Check columns
  cols = c("sample_id", "site_id", "date", "variable_id", "fraction", 
           "qualifier", "value", "unit", "value_fin", "loq", "lod", "sample_type", 
           "sample_duration", "comment")
  cols_OK <- all(cols == names(x))
  if (cols_OK) {
    print('Colums: PASS')
  } else {
    message('Columns! : NOPASS!')
    print("Missing cols:",  cols[!cols %chin% names(x)])
  }
  # check class
  colclas <- sapply(x, class) 
  colclas_OK <- all(colclas == c(rep('character', 2), 'Date', 'integer', 
                                 rep('character', 2), 'numeric', 'character', 
                                 'numeric', 'numeric', 'numeric', 
                                 rep('character', 3)
                                 ))
  if (colclas_OK) {
    print('ColClasses: PASS')
  } else {
    message('Wrong Colclasses! : NOPASS')
  }
  
  # check NAs (no NA in "sample_id""site_id""date""variable_id" "unit""value_fin")
  NA_OK <- !any()
  if (NA_OK) {
    print('NA values: PASS')
  } else {
    message('NAs somewhere! : NOPASS')
    print(x[is.na(x[ , c(1:4, 8, 9), with = FALSE])])
  }
  
  # check ranges
  # ph
  ph_rg <- nrow(x[variable_id == 506 & (value_fin >= 14 | value_fin <= 0)]) == 0
  if (ph_rg) {
    print('pH Range: PASS')
  } else {
    message('pH Range : NOPASS')
    print(x[variable_id == 506 & (value_fin >= 14 | value_fin <= 0)])
  }
  # all vars except temp should be > 0
  nn <- nrow(x[!unit %chin% '°C' & value_fin < 0]) == 0
  if (nn) {
    print('Non-negative: PASS')
  } else {
    message('negative data!: NOPASS ')
    print(x[!unit %chin% '°C' & value_fin < 0])
  }
  
  # check units
 
  mm <- merge(unique(x[ , list(variable_id, unit)]), 
              unique(variables[ , list(variable_id, unit)]), by = 'variable_id')
  units_OK <- nrow(mm[!unit.x == unit.y & !variable_id %in% c(622, 623)]) == 0
  if (units_OK) {
    print('Units: PASS')
  } else {
    print(mm[!unit.x == unit.y & !variable_id %in% c(622, 623)])
    message('Units: NOPASS')
  }
  
  # full duplicates
  dups_OK <- !any(duplicated(x))
  if (dups_OK) {
    print('duplicates: PASS')
  } else {
    message('duplicates in data! : NOPASS')
  }
  
  # daily dups
  dups_OK2 <- !any(duplicated(x, by = c('site_id', 'date', 'variable_id')))
  if (dups_OK2) {
    print('daily duplicates: PASS')
  } else {
    message('daily duplicates in data! : NOPASS')
  }
  
  # no NA ins value_fin
  value_fin_OK <- !any(is.na(x[ , value_fin]))
  if (value_fin_OK) {
    print('value_finNA: PASS')
  } else {
    message('NAs in value_fin! : NOPASS')
  }
  
  qual_OK <- all(unique(x[ , qualifier]) %chin% c("<", "=", ">"))
  if (qual_OK) {
    print('qualifier: PASS')
  } else {
    message('Incorrect values in qualifier : NOPASS')
  }
  
  
  
  sampt_OK <- all(unique(x[ , sample_type]) %chin% c("individual", "composite", NA))
  if (sampt_OK) {
    print('sample_type: PASS')
  } else {
    message('Incorrect values in sample type : NOPASS')
  }
  
  # total
  if (all(cols_OK, colclas_OK, ph_rg, nn, units_OK, dups_OK, dups_OK2, value_fin_OK, qual_OK, sampt_OK)) {
    message('Total: PASS')
    return(TRUE)
  } else {
    message('Total: NOPASS')
    return(FALSE)
  }
}

# check if the two tables join smoothly
chk_join <- function(sites, chem){
  nm1 <- sum(is.na(match(sites$site_id, chem$site_id))) 
  if (nm1 > 0)
    message('sites without chemdata!')  
  nm2 <- sum(is.na(match(chem$site_id, sites$site_id))) 
  if (nm2 > 0)
    message('more sites in chem, than in sites!')
  out <- ifelse(nm1 == 0 & nm2 == 0, TRUE, FALSE)
  return(out)
}

# graphical check
chk_outl <- function(x){
  # check outliers
  vars <- c('Cl', 'GH', 'NH4', 'NO2', 'NO3', 'P_tot', 'SO4', 'TOC', 
            'OPO4', 'T', 'pH', 'O2', 'BSB5', 'EC', 'SAT')
  take <- x[variable %in% vars, ]
  take[,variable:= factor(variable)]
  ggplot(take, aes(x = value)) +
    geom_density() +
    facet_wrap(~variable, scales = 'free')
}

# take maximum of daily duplicates
clean_dmax <- function(chem){
  # check daily dups
  dups <- unique(chem[duplicated(chem, by = c('site_id', 'date', 'variable_id')), list(site_id, date, variable_id)])
  # which to aggregate
  toagg <- chem[sample_id %in% dups[ , sample_id] & date %in% dups[ , date] & variable_id %in% dups[ , variable_id]]
  noagg <- chem[!(sample_id %in% dups[ , sample_id] & date %in% dups[ , date] & variable_id %in% dups[ , variable_id])]
  # aggregate (take maximum of value_fin)
  agg <- toagg[toagg[ , .I[which.max(value_fin)], by = c('site_id', 'date', 'variable_id')]$V1]
  # add aggregated data
  return(rbindlist(list(noagg, agg)))
}

# returns string w/o leading or trailing whitespace
trim <- function (x) gsub("^\\s+|\\s+$", "", x)


# transform coordinates
# crs1, crs2 : crs of input and output
# example: mytrans(easting, northing, CRS("+init=epsg:31468"), CRS("+init=epsg:31467"))
# return a list of transformed coordinates
mytrans <- function(easting, northing, crs1, crs2){
  df <- data.frame(easting = easting, northing = northing)
  coordinates(df) <- c('easting', 'northing')
  proj4string(df) <- crs1
  df_t <- spTransform(df, crs2)
  return(list(coordinates(df_t)[,1], coordinates(df_t)[,2]))
}

# This function returns TRUE wherever elements are the same, including NA's,
# and FALSE everywhere else.
compareNA <- function(v1,v2) {
  same <- (v1 == v2) | (is.na(v1) & is.na(v2))
  same[is.na(same)] <- FALSE
  return(same)
}


###-----------------------------------------------------------------------------
### data cleaning
# function for bulk processing of saxonia data
extract_sax <- function(f){
  message("DB-processed :", f)
  d <- mdb.get(f, stringsAsFactors = FALSE) # need mdb-tools
  
  ### extract sites
  sites <- data.table(d[["LIMS_STAMMDAT_OW2"]])
  keep <- c('MKZ.NEU', 'NAME', 'GEWAESSER', 'RW', 'HW', 'STANDGEW')
  if(any(!keep %in% names(sites)))
    sites[ , keep[!keep %in% names(sites)] := 'NR']
  sites <- sites[ ,keep, with = FALSE]
  sites <- data.frame(lapply(sites, unclass), stringsAsFactors = FALSE)
  
  ### chem
  # get data tables
  ow <- d[grep('owdat_', names(d), value = TRUE)]
  # remove labelled class
  ow <- lapply(ow, function(y) data.frame(lapply(y, unclass), stringsAsFactors = FALSE))
  ow <- lapply(ow, setDT)
  
  # rm unused cols
  keep <- c('BEZUGP', 'EINHEIT', 'PRAEFIX', 'ERGEBNIS', 'MKZ.', 'MKZ', 'PBEMERK',
            'PARAM', 'PDATUM', 'BESTGR',  'NACHWGR', 'PGERAET', 'PMAT', 'PENTART')
  ow <- lapply(ow, function(x) x[, names(x) %in% keep, with = FALSE])
  ow <- lapply(ow, function(y) {
    if ('MKZ' %in% names(y))
      setnames(y, 'MKZ', 'MKZ.')
    return(y)
  })
  ow <- lapply(ow, setcolorder, c("BEZUGP", "EINHEIT", "PRAEFIX", "ERGEBNIS", "MKZ.", "PARAM", 
                                 "PBEMERK", "PDATUM", "PENTART", "PGERAET", "PMAT", "NACHWGR", 
                                 "BESTGR"))
  owd <- rbindlist(ow)
  
  # merge tables
  # units
  unit_tab <- data.table(d[["ZENTRAL_ZTEINH"]][ , c('SCHLUESSEL', 'KUERZEL')])
  setnames(unit_tab, c('EINHEIT', 'unit'))
  owd <- merge(owd, unit_tab, by = 'EINHEIT')
  rm(unit_tab)
  owd[ , EINHEIT := NULL]
  
  # bzp
  bzp_tab <- data.table(d[["ZENTRAL_ZTBEZUGP"]][ ,c("ID.BEZUGP", "PROBENBEZUG")])
  setnames(bzp_tab, c("BEZUGP", "PROBENBEZUG"))
  setkey(bzp_tab, BEZUGP)
  if ('integer' %in% class(owd$BEZUGP)) {
    # rm leading zeros
    bzp_tab[ , BEZUGP := gsub("^0+", "", BEZUGP)]
    owd[ , BEZUGP := as.character(BEZUGP)]
  }
  owd <- merge(owd, bzp_tab, by = 'BEZUGP')
  rm(bzp_tab)
  owd[ , BEZUGP := NULL]
  
  # param
  param_tab <- data.table(d[["ZENTRAL_ZTPARAM"]][ , c('SCHLUESSEL', "PARAMETER")])
  setnames(param_tab, c('PARAM', 'variable'))
  #   unique(owd[is.na(match(owd$PARAM, param_tab$PARAM)), PARAM])
  # lost some parameters
  owd <- merge(owd, param_tab, by = 'PARAM') 
  rm(param_tab)
  owd[ , PARAM := NULL]
  return(list(sites = sites, chem = owd))
}

###-----------------------------------------------------------------------------
#### Some function for memory management
#' garbage collection
cleanMem <- function(n=5) { for (i in 1:n) gc() }

# improved list of objects
.ls.objects <- function (pos = 1, pattern, order.by,
                         decreasing=FALSE, head=FALSE, n=5) {
  napply <- function(names, fn) sapply(names, function(x)
    fn(get(x, pos = pos)))
  names <- ls(pos = pos, pattern = pattern)
  obj.class <- napply(names, function(x) as.character(class(x))[1])
  obj.mode <- napply(names, mode)
  obj.type <- ifelse(is.na(obj.class), obj.mode, obj.class)
  obj.prettysize <- napply(names, function(x) {
    capture.output(print(object.size(x), units = "MB")) })
  obj.size <- napply(names, object.size)
  obj.dim <- t(napply(names, function(x)
    as.numeric(dim(x))[1:2]))
  vec <- is.na(obj.dim)[, 1] & (obj.type != "function")
  obj.dim[vec, 1] <- napply(names, length)[vec]
  out <- data.frame(obj.type, obj.size, obj.prettysize, obj.dim)
  names(out) <- c("Type", "Size", "PrettySize", "Rows", "Columns")
  if (!missing(order.by))
    out <- out[order(out[[order.by]], decreasing=decreasing), ]
  if (head)
    out <- head(out, n)
  out
}

# shorthand
lsos <- function(..., n=10) {
  .ls.objects(..., order.by="Size", decreasing=TRUE, head=TRUE, n=n)
}

#lsos()


### ggplot2 theme

mytheme <- theme_bw(base_size = 12, base_family = "Helvetica") + 
  theme(
    # panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(), 
        text = element_text(size = 14),
        axis.text = element_text(size = 12),
        axis.title.x = element_text(size = 14, face = "bold", vjust = 0),
        axis.title.y = element_text(size = 14, face = "bold", vjust = 1),
        legend.key = element_blank(),
        strip.background = element_blank(),
        strip.text = element_text(size = 14, face = 'bold'))

phdtheme <- theme_bw(base_size = 12, base_family = "Helvetica") + 
  theme(
    # panel.grid.major = element_blank(),
    panel.grid.minor = element_blank(), 
    text = element_text(size = 14),
    axis.text = element_text(size = 12),
    axis.title.x = element_text(size = 14, face = "bold", vjust = 0),
    axis.title.y = element_text(size = 14, face = "bold", vjust = 1),
    legend.key = element_blank(),
    strip.background = element_blank(),
    strip.text = element_text(size = 14, face = 'bold'))

theme_end <- theme_bw(base_size = 12, base_family = "Helvetica") + 
  theme(panel.grid.minor = element_blank(), 
        axis.title.x = element_text(face = "bold", vjust = 0),
        axis.title.y = element_text(face = "bold", vjust = 1),
        legend.key = element_blank(),
        strip.background = element_blank(),
        strip.text = element_text(face = 'bold'))

mytheme2 <- theme_bw(base_size = 20, base_family = "Helvetica") + 
  theme(
#     text = element_text(size = 14),
#     axis.text = element_text(size = 12),
#     axis.title.x = element_text(size = 14, face = "bold", vjust = 0),
#     axis.title.y = element_text(size = 14, face = "bold", vjust = 1),
    # legend.position = "bottom",
    legend.key = element_blank(),
    strip.background = element_blank(),
    strip.text = element_text(size = 14, face = 'bold'))


## geometric mean, from http://stackoverflow.com/questions/2602583/geometric-mean-is-there-a-built-in
geomean = function(x, na.rm=TRUE, zero.propagate = FALSE){
  if (any(x < 0, na.rm = TRUE)) {
    return(NaN)
  }
  if (zero.propagate) {
    if (any(x == 0, na.rm = TRUE)) {
      return(0)
    }
    exp(mean(log(x), na.rm = na.rm))
  } else {
    exp(sum(log(x[x > 0]), na.rm = na.rm) / length(x))
  }
}



###-----------------------------------------------------------------------------
#### Functions for REGNIE-data
#' get filename from date
#' @import data.table

file_regnie <- function(date){
  sf <- paste0('ra', year(date), 'm')
  y <- substr(year(date), 3, 4)
  m <- month(date)
  m <- ifelse(nchar(m) < 2, paste0(0, m), m)
  d <- mday(date)
  d <- ifelse(nchar(d) < 2, paste0(0, d), d)
  fi <- paste0('ra', y, m, d, '.gz')
  out <- file.path(sf, fi)
  return(out)
}



#' Read REGNIE data into R
#' @param file file path to gzip archive
#' @param crs reproject raster to specified CRS? If NULL the returned raster is WGS84.
#' @description This function reads REGNIE data and returns a projected raster object
#' @example read_regnie(file.path(regpath, sf, fi), crs = CRS('+init=epsg:31467'))
read_regnie <- function(file, crs = NULL){
  # read file
  zz <- gzfile(file)
  cont <- readLines(zz, n = 971)
  close(zz)
  # replace NA value
  cont <- gsub('-999', ' NA ', cont, fixed = TRUE)
  
  r <- raster(matrix(scan(textConnection(cont), what = integer(0), n = 971*611), 
                     nrow = 971, ncol = 611, byrow = TRUE))
  # r <-  raster(as.matrix(read.table(textConnection(cont), colClasses = 'integer', nrows = 971)))
  extent(r) <- c(5.833333, 16, 47, 55.083333)
  crs(r) <- "+proj=longlat"
  if (!is.null(crs)) {
    # reproject
    # r <- projectRaster(from = r, crs = CRS('+init=epsg:31467'))
    # projectRaster is slow - use gdalwarp
    tf <- tempfile(tmpdir = file.path(cachedir, 'tmp/regnie_raster'), fileext = '.tif')
    tf2 <- tempfile(tmpdir = file.path(cachedir, 'tmp/regnie_raster'), fileext = '.tif')
    writeRaster(r, tf)
    system(command = paste("gdalwarp -t_srs \'+init=epsg:31467\' -r near -overwrite", 
                           tf,
                           tf2))
    r <- raster(tf2)
  }
  # scale to mm
  r <- r/10
  if (!is.null(crs))
    unlink(c(tf, tf2))
  return(r)
}

# fix for ggsave [error: " plot should be a ggplot2 plot"]
ggsave_fix <- ggplot2::ggsave; body(ggsave_fix) <- body(ggplot2::ggsave)[-2]

#' log spaced sequence
#' @source http://grokbase.com/t/r/r-help/099nykfnnf/r-logarithmic-seq
#' 
logspace <- function( d1, d2, length.out) {
  exp(log(10)*seq(d1, d2, length.out=length.out))
}

