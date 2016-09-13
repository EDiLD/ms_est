if (!exists("prj")) {
  stop("You need to create a object 'prj' that points to the top folder, 
       e.g. prj <- '/home/edisz/Documents/Uni/Projects/PHD/4BFG/Project'!")
} else {
  source(file.path(prj, "src", "load.R"))
}


### ----------------------------------------------------------------------------
### Comparison of catchment area and stream width


### ----------
# TH has good data 
#! interessting data about ezg - width
th <- setDT(read_excel(file.path(datadir, 'physico_chem/TH/new/', 'Thueringen_PSM_FG_2005-14.xlsx')))
th <- th[ , c(8,10), with = FALSE]
setnames(th, c('ezg', 'width'))
th <- th[complete.cases(th)]
th[ , dataset := 'TH']

ggplot(th, aes(x = ezg, y = width)) +
  geom_point() +
  scale_x_log10() +
  scale_y_log10() 

# log-log linear model
loglm <- lm(log(width) ~ log(ezg), data = th)
coef(loglm)
summary(loglm)
# sum of squared error
sum((th$width - exp(fitted(loglm)))^2)

# power model
nlm <- nls(width ~ a*ezg^b, data = th, start = list(a = exp(coef(loglm)[1]), b = coef(loglm)[2]))
summary(nlm)
sum((th$width - fitted(nlm))^2)
# better odel

pdat <- data.frame(ezg = seq(1, 10000, length.out = 1000))
pdat$ploglm <- predict(loglm, newdata = pdat)
pdat$pnlm <- predict(nlm, newdata = pdat)

ggplot() +
  geom_point(data = th, aes(x = ezg, y = width)) +
  scale_x_log10() +
  scale_y_log10() +
  geom_line(data = pdat, aes(x = ezg, y = pnlm), col = 'blue') +
  geom_hline(aes(yintercept = 2)) +
  mytheme


### ----------
### data from voss 2015
voss <- fread(file.path(datadir, 'misc/AG/katha', 'totalData.csv'))
voss <- voss[ , c(5, 34), with = FALSE]
setnames(voss, c('width', 'ezg'))
voss[ , c('width', 'ezg') := list(width = as.numeric(gsub(',', '.', width)),
           ezg = as.numeric(gsub(',', '.', ezg)))]
voss[ , dataset := 'Voss2015']




### ----------
### data from fernandez 2015
diego <- fread(file.path(datadir, 'misc/AG/diego', 'diego.csv'))
diego  <- diego[ , c(3, 2), with = FALSE]
setnames(diego, c('width', 'ezg'))
diego[ , dataset := 'Fernandez2015']


### ----------
### data from ufz
ufz <- fread(file.path(datadir, 'misc/UFZ', '2_WIDTH_EZG.txt'))
ufz <- ufz[ , list(width = width/100, ezg = Shape_Area)]
ufz[ , dataset := 'ufz']
ufz <- ufz[!(duplicated(ufz) | is.na(width)), ]
# drop observation 15 (= outlier)
ufz <- ufz[-15, ]


# a function to fit power model and export predicts
foo <- function(data, n = 1000){
  loglm <- lm(log(width) ~ log(ezg), data = data)
  nlm <- nls(width ~ a*ezg^b, data = data, start = list(a = exp(coef(loglm)[1]), b = coef(loglm)[2]))
  pdat <- data.frame(ezg = seq(1, 10000, length.out = n))
  pdat$width <- predict(nlm, newdata = pdat)
  return(pdat)
}

ew <- rbind(th, voss, diego, ufz)
cols <- c('darkred', 'steelblue', 'darkgreen', 'darkorange3')


p_ezg_width <- ggplot() +
  scale_color_manual(values = cols, name = 'Datensatz',
                     labels = c("Thüringen", "Voss 2015", "Fernandez 2015", 'UFZ'),
                     breaks = c("TH", "Voss2015", "Fernandez2015", "ufz")) +
  geom_rect(aes(xmin = 0, xmax = 7, ymin = 0, ymax = 1), alpha = 0.3) +
  geom_rect(aes(xmin = 0, xmax = 23.5, ymin = 0, ymax = 2), alpha = 0.3) +
  scale_x_log10(breaks = c(1, 7, 23.5, 100, 1000)) +
  scale_y_log10(breaks = c(0.5, 1, 2, 10, 50)) +
  mytheme +
  geom_point(data = ew, aes(x = ezg, y = width, col = dataset), size = 2.5) +
  geom_line(data = foo(ew), aes(x = ezg, y = width), col = 'black') +
  geom_line(data = foo(ew[ew$dataset == 'TH', ]), aes(x = ezg, y = width), col = cols[2], linetype = 'dashed') +
  geom_line(data = foo(ew[ew$dataset == 'Voss2015', ]), aes(x = ezg, y = width), col = cols[4], linetype = 'dashed') +
  geom_line(data = foo(ew[ew$dataset == 'Fernandez2015', ]), aes(x = ezg, y = width), col = cols[1], linetype = 'dashed') +
  geom_line(data = foo(ew[ew$dataset == 'ufz', ]), aes(x = ezg, y = width), col = cols[3], linetype = 'dashed') +
  labs(y = 'mittl. Breite [m]', x = expression(paste('EZG [', m^2, ']'))) +
  annotate('text', x = 2, y = 50, label = 'Breite = 0.32 EZG^0.58', size = 5, hjust = 0) +
  theme(legend.position = 'right')
p_ezg_width
# # export for presentation
# ggsave("/home/edisz/Documents/Uni/Projects/PHD/4BFG/Project/reports/20151006_workshop_koblenz/pres/fig/ezg_width.pdf", p_ezg_width)


# log-log linear model
loglm <- lm(log(width) ~ log(ezg), data = ew)
# power model
nlm <- nls(width ~ a*ezg^b, data = ew, start = list(a = exp(coef(loglm)[1]), b = coef(loglm)[2]))
summary(nlm)

exp(log(2/coef(nlm)[1]) / coef(nlm)[2])
exp(log(1/coef(nlm)[1]) / coef(nlm)[2])

