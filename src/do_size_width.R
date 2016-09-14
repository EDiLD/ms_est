if (!exists("prj")) {
  stop("You need to create a object 'prj' that points to the top folder, 
       e.g. prj <- '/home/edisz/Documents/Uni/Projects/PHD/4BFG/Paper/ms_est' or
       prj <- '/home/user/Documents/projects_git/ms_est'!")
} else {
  source(file.path(prj, "src", "load.R"))
}


### ----------------------------------------------------------------------------
### Comparison of catchment area and stream width


### ----------
# TH has good data 
th <- fread(file.path(cachedir, 'TH_2005-14.csv'))
th <- th[ , c(4, 5), with = FALSE]
setnames(th, c('ezg', 'width'))
th <- th[complete.cases(th)]
th[ , dataset := 'TH']

### ----------
### data from voss 2015
voss <- fread(file.path(cachedir, 'voss_2015.csv'))
voss <- voss[ , c(2, 3), with = FALSE]
setnames(voss, c('width', 'ezg'))
voss[ , c('width', 'ezg') := list(width = as.numeric(gsub(',', '.', width)),
           ezg = as.numeric(gsub(',', '.', ezg)))]
voss[ , dataset := 'Voss2015']


### ----------
### data from fernandez 2015
diego <- fread(file.path(cachedir, 'fernandez_2015.csv'))
diego  <- diego[ , c(3, 2), with = FALSE]
setnames(diego, c('width', 'ezg'))
diego[ , dataset := 'Fernandez2015']




# a function to fit power model and export predicts
foo <- function(data, n = 1000){
  # fit log-log linear model for initial values
  loglm <- lm(log(width) ~ log(ezg), data = data)
  # fit power model
  nlm <- nls(width ~ a*ezg^b, data = data, start = list(a = exp(coef(loglm)[1]), b = coef(loglm)[2]))
  # predict from fitted model
  pdat <- data.frame(ezg = seq(1, 10000, length.out = n))
  pdat$width <- predict(nlm, newdata = pdat)
  return(pdat)
}

ew <- rbind(th, voss, diego)
cols <- c('darkred', 'steelblue', 'darkgreen')


p_ezg_width <- ggplot() +
  scale_color_manual(values = cols, name = 'Dataset',
                     labels = c("TH", "Voss 2015", "Fernandez 2015"),
                     breaks = c("TH", "Voss2015", "Fernandez2015")) +
  geom_rect(aes(xmin = 0, xmax = 25, ymin = 0, ymax = 2.08), alpha = 0.3) +
  scale_x_log10(breaks = c(1, 25, 100, 1000)) +
  scale_y_log10(breaks = c(0.5, 1, 2, 10, 50)) +
  mytheme +
  geom_point(data = ew, aes(x = ezg, y = width, col = dataset), size = 2.5) +
  geom_line(data = foo(ew), aes(x = ezg, y = width), col = 'black') +
  geom_line(data = foo(ew[ew$dataset == 'TH', ]), aes(x = ezg, y = width), col = cols[2], linetype = 'dashed') +
  geom_line(data = foo(ew[ew$dataset == 'Voss2015', ]), aes(x = ezg, y = width), col = cols[3], linetype = 'dashed') +
  geom_line(data = foo(ew[ew$dataset == 'Fernandez2015', ]), aes(x = ezg, y = width), col = cols[1], linetype = 'dashed') +
  labs(y = 'Stream width [m]', x = expression(paste('Catchment Size [', km^2, ']'))) +
  annotate('text', x = 2, y = 50, label = 'Width = 0.33 Size^0.57', size = 5, hjust = 0) +
  theme(legend.position = 'right')
# p_ezg_width
# # export for presentation
ggsave(file.path(prj, 'supplement/width_size.pdf'), p_ezg_width, width = 7, height = 5)

# get prediction

loglm <- lm(log(width) ~ log(ezg), data = ew)
nlm <- nls(width ~ a*ezg^b, data = ew, start = list(a = exp(coef(loglm)[1]), b = coef(loglm)[2]))
coef(nlm)

# a catchment of 25km corresponds to a width of
coef(nlm)[1]*25^coef(nlm)[2]
# 2.1m  
