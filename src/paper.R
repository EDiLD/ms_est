

# ### ----------------------------------------------------------------------------
# # Modelling influence of EZG & Landuse / threshold detection
# 
# # join sites with info
# setkey(psm_sites, site_id)
# setkey(psm_sites_info, site_id)
# psm_sites_wi <- psm_sites[psm_sites_info[ , list(site_id, use, ezg_fin, agri_fin)]]
# 
# 
# ### Endpoint: Number of RAK execeedances
# raks <- psm_variables[!is.na(rak_uba), list(variable_id, name, cas, pgroup, rak_uba)]
# 
# samples_raks <- psm_samples[variable_id %in% raks$variable]
# setkey(samples_raks, variable_id)
# setkey(raks, variable_id)
# samples_raks <- raks[samples_raks]
# samples_raks[ , g_rak_uba := value_fin > rak_uba]
# samples_raks <- samples_raks[month(date) %in% c(5, 6, 7)]
# 
# rak_thrs_sites <- samples_raks[!is.na(g_rak_uba) , list(sum_g_rak_uba_abs = sum(g_rak_uba),
#                                   sum_g_rak_uba_prop = sum(g_rak_uba) / length(g_rak_uba),
#                                   n_samp = length(g_rak_uba)), by = site_id][order(sum_g_rak_uba_prop, decreasing = TRUE)]
# setkey(rak_thrs_sites, site_id)
# setkey(psm_sites, site_id)
# 
# w3 <- psm_sites_wi[rak_thrs_sites]
# w3 <- w3[ezg_fin < 100 & ezg_fin > 1 & n_samp > 10]
# 
# ggplot(w3, aes(x = ezg_fin, y = sum_g_rak_uba_prop)) +
#   geom_point() +
#   geom_smooth()
# 
# 
# ggplot(w3, aes(x = agri_fin, y = sum_g_rak_uba_prop)) +
#   geom_point() +
#   geom_smooth()
# 
# 
# # model using zero inflated beta GAM
# 
# 
# mdata <- na.omit(w3[ , list(ezg_fin, agri_fin, sum_g_rak_uba_prop)])
# df <- 4
# method <- 'GAIC'
# mod <- gamlss(sum_g_rak_uba_prop ~ pb(agri_fin, df = df, method = method) + pb(ezg_fin, df = df, method = method), 
#               nu.formula = ~pb(agri_fin, df = df, method = method) + pb(ezg_fin, df = df, method = method), 
#               family = BEINF0,  # ZAGA
#               data = mdata
# )
# plot(mod)
# summary(mod)
# 
# # png('/home/edisz/Documents/Uni/Projects/PHD/4BFG/Project/reports/20151131_Endbericht/chapters/figs/lu_ezg_rak_1.png', 
# #     width = 2400, height = 2400, res = 300)
# term.plot(mod, what = 'mu', data = mdata, rug = TRUE, 
#           ylim = 'common')
# # dev.off()
# 
# term.plot(mod, what = 'nu', data = mdata, rug = TRUE, 
#           ylim = 'common')
# 
# 
# 
# xyz <- expand.grid(agri_fin = seq(0, 1, 0.02), ezg_fin = seq(0, 100, 2))
# xyz$y <- predict(mod, newdata = xyz, type = 'response', what = 'mu')
# 
# # pdf("/home/edisz/Documents/Uni/Projects/PHD/4BFG/Project/reports/20151006_workshop_koblenz/pres2/fig/agri_ezg_tumax_wireframe.pdf")
# wireframe(y ~ agri_fin * ezg_fin, xyz, 
#           scales = list(arrows=FALSE,cex=.5,tick.number = 10),
#           xlab = 'Landwirtschaft [%]', 
#           ylab = 'EZG [km2]', 
#           zlab = 'Anteil RAK Überschreitungen'
#           )
# # dev.off()
# 
# 
# 
# 
# ### Endpoint: TU
# # filter to only may + june
# mtu_mj <- psm_maxtu[, m := month(date)]
# setkey(mtu_mj, site_id)
# setkey(psm_sites_info, site_id)
# mtu_mj <- psm_sites_info[ , list(site_id, ezg_fin, agri_fin)][mtu_mj]
# mtu_mj[ , ltu := log10(maxtu + 0.00000001)]
# 
# ggplot(mtu_mj[ezg_fin < 150 & ezg_fin > 0.5 & m %in% c(5, 6, 7)], 
#        aes(x = ezg_fin, y = ltu)) +
#   geom_point(alpha = 0.2) +
#   geom_smooth() +
#   theme_bw() 
# 
# ggplot(mtu_mj[ezg_fin < 150 & ezg_fin > 0.5 & m %in% c(5, 6, 7)], 
#        aes(x = agri_fin, y = ltu)) +
#   geom_point(alpha = 0.2) +
#   geom_smooth() +
#   theme_bw()
# 
# 
# ggplot(mtu_mj[ezg_fin < 150 & ezg_fin > 0.5 & m %in% c(5, 6, 7)], 
#        aes( y = agri_fin, x = ezg_fin, col = ltu)) +
#   geom_point() +
#   theme_bw() +
#   scale_color_gradient(low = 'grey80', high = 'red')
# 
# 
# 
# 
# # model log10(TU) as censored gaussian GAM
# moddf <- mtu_mj[ , list(ltu, ezg_fin, agri_fin, site_id, date, m)]
#  # remove nas
# moddf <- na.omit(moddf)
# moddf <- moddf[ezg_fin < 100 & ezg_fin > 0.5 & m %in% c(5, 6, 7)]
# 
# 
# mod <- gamlss(Surv(moddf$ltu, moddf$ltu > -8, type = 'left') ~ pb(ezg_fin) + pb(agri_fin), 
#               data = moddf, 
#               family = cens(NO, type = 'left'))
# summary(mod)
# plot(mod)
# 
# 
# # png('/home/edisz/Documents/Uni/Projects/PHD/4BFG/Project/reports/20151131_Endbericht/chapters/figs/lu_ezg_tu_1.png', 
# #     width = 2400, height = 2400, res = 300)
# termplot(mod, daza = moddf, terms = 1, ylim = c(-2, 2), se = TRUE)
# # dev.off()
# 
# # png('/home/edisz/Documents/Uni/Projects/PHD/4BFG/Project/reports/20151131_Endbericht/chapters/figs/lu_ezg_tu_2.png', 
# #     width = 2400, height = 2400, res = 300)
# termplot(mod, daza = moddf, terms = 2, ylim = c(-2, 2), se = TRUE)
# # dev.off()
# 
# 
# xyz <- expand.grid(agri_fin = seq(0, 1, 0.02), ezg_fin = seq(0, 100, 2))
# xyz$y <- predict(mod, newdata = xyz)
# 
# 
# ggplot(xyz, aes(x = ezg_fin,  y = agri_fin, z = y, fill = y)) +
#   geom_tile() +
#   geom_contour() +
#   scale_fill_gradient(low = 'grey80', high = 'red')
# 
# 
# 
# 
# 
# wireframe(y ~ agri_fin * ezg_fin, xyz, 
#           scales = list(arrows = FALSE, cex = .5,tick.number = 10),
#           xlab = 'Landwirtschaft [%]', 
#           ylab = 'EZG [km2]', 
#           zlab = 'log10(maxTU)',
#           distance = .4)
# # dev.off()
# 
# 
# 
# 
# 
# ### Endpoint: Single substances
# 
# 
# ### What compounds are where reponsible for maxTU?
# 
# 
# 
# ### ---------------------------------------------------------------------------
# ### Chapter 5
# 
# 
# ### Beispielhafte Auwertung der Daten
# 
# # "Landwirtschaftliche Kleingewässer" definiert als EZG < 30km^2 und Landwirtschaftlicher Anteil > 50%
# # "Landwirtschaftliche Gewässer" 30 <= EZG < 100km^2 und Landwirtschaflicher Anteil > 50%
# # "Nicht Landwirtschaftliche Kleingewässer" definiert als EZG < 30km^2 und Landwirtschaftlicher Anteil < 25%
# # "Nicht landwirtschaftliche Gewässer" 30 <= EZG < 100km^2 und Landwirtschaflicher Anteil <25%
# 
# # klassify stream types
# psm_sites_info[ezg_fin < 30 & agri_fin > 0.4] # 843 MS
# psm_sites_info[ezg_fin < 30 & agri_fin > 0.4, type := 'klein, agrar']  # small agricultural
# psm_sites_info[ezg_fin >= 30 & ezg_fin <= 100 & agri_fin > 0.4] # 462 MS
# psm_sites_info[ezg_fin >= 30 & ezg_fin <= 100 & agri_fin > 0.4, type := 'mittel, agrar'] # agricultural
# psm_sites_info[ezg_fin < 30 & agri_fin < 0.1] # 254 MS
# psm_sites_info[ezg_fin < 30 & agri_fin < 0.1, type := 'klein, nicht-agrar'] # small non - agricultural
# psm_sites_info[ezg_fin >= 30 & ezg_fin <= 100 & agri_fin < 0.1] # 144 MS
# psm_sites_info[ezg_fin >= 30 & ezg_fin <= 100 & agri_fin < 0.1, type := 'mittel, nicht-agrar'] # non-agricultural
# 
# psm_sites_info[!is.na(type)]
# psm_sites_info[is.na(type)]
# 
# psm_sites_info[ezg_fin > 100] # 418
# psm_sites_info[agri_fin >= 0.1 & agri_fin <= 0.4] # 642
# 
# 
# ### ---------------------------------------
# ### Detections as Endpoint
# setkey(psm_sites_info, site_id)
# setkey(psm_samples, site_id)
# 
# samp <- psm_sites_info[,list(site_id, type)][psm_samples]
# samp <- samp[!is.na(type)]
# 
# # most often detected compounds
# dd <- samp[ , list(prop = sum(value_fin > 0) / length(value_fin), 
#                    abs = sum(value_fin > 0),
#                    n = length(value_fin)), 
#             by = list(variable_id, type)]
# 
# setkey(dd, variable_id)
# setkey(psm_variables, variable_id)
# dd <- psm_variables[ , list(variable_id, name, psm_type)][dd][order(prop, decreasing = TRUE)]
# 
# # restrict and order by appearance in sawb
# ord <- dd[type == 'klein, agrar'][order(prop, decreasing = TRUE), name]
# ord
# # take only fist 20 compounds
# ord <- ord[1:30]
# dd <- dd[name %chin% ord]
# dd$name <- factor(dd$name, levels = rev(ord))
# dd$type <- factor(dd$type, levels = unique(dd$type)[c(1, 4, 2, 3)])
# 
# ppos <- ggplot() +
#   geom_point(data = dd[prop > 0.10], 
#              aes(x = prop * 100, y = name, size = n, col = psm_type)) +
#   theme_bw() +
#   theme(legend.key = element_blank()) +
#   scale_size('# Proben', trans = "log10") +
#   labs(y = 'Wirkstoff', x = 'Anteil an Proben > LOQ [%]') +
#   scale_color_manual(name = 'Gruppe', 
#                      values = c('#377eb8', '#4daf4a', '#984ea3'),
#                      labels = c('Fungizid', 'Herbizid', 'Metabolit')) +
#   facet_wrap(~type) +
#   theme_end +
#   guides(colour = guide_legend(override.aes = list(size = 5)))
# ppos
# # ggsave("/home/edisz/Documents/Uni/Projects/PHD/4BFG/Project/reports/20151131_Endbericht/chapters/figs/ppos.png", width = 10, height = 7,
# #        ppos)
# 
# # anteil in kleingewässer immer höher als in größeren agri-gewässern
# # metabolite
# # unterschied agrar-nicht-agrar deutlich
# # geringer unterschied in der Abfolge der Wirkstoffe
# 
# 
# 
# 
# ### ---------------------------------------
# ### UQN as Endpoint
# 
# # restrict samples
# thrs <- psm_variables[!is.na(wrrl_zhkuqn), list(variable_id, name, cas, pgroup, wrrl_zhkuqn)]
# samples_thrs <- psm_samples[variable_id %in% thrs$variable]
# samples_thrs <- psm_sites_info[,list(site_id, type)][samples_thrs]
# samples_thrs <- samples_thrs[!is.na(type)]
# setkey(samples_thrs, variable_id)
# setkey(thrs, variable_id)
# samples_thrs <- thrs[samples_thrs]
# samples_thrs[ , g_wrrl_zhkuqn := value_fin > wrrl_zhkuqn]
# 
# dd <- samples_thrs[!is.na(g_wrrl_zhkuqn) , list(abs = sum(g_wrrl_zhkuqn),
#                                                             prop = sum(g_wrrl_zhkuqn) / length(g_wrrl_zhkuqn),
#                                                             n = length(g_wrrl_zhkuqn))
#                    , by = list(variable_id, type)][order(prop, decreasing = TRUE)]
# dd
# 
# setkey(dd, variable_id)
# dd <- psm_variables[ , list(variable_id, name, psm_type)][dd][order(prop, decreasing = TRUE)]
# 
# # restrict and order by appearance in sawb
# ord <- dd[type == 'klein, agrar'][order(prop, decreasing = TRUE), name]
# ord
# # take only fist 20 compounds
# ord <- ord[1:20]
# dd <- dd[name %chin% ord]
# dd$name <- factor(dd$name, levels = rev(ord))
# dd$type <- factor(dd$type, levels = unique(dd$type)[c(1, 4, 2, 3)])
# # # manually change type for irgarol
# # dd[variable_id == 375, psm_type := 'biocide']
# 
# puqn <- ggplot() +
#   geom_point(data = dd, 
#              aes(x = prop * 100, y = name, size = n, col = psm_type)) +
#   theme(legend.key = element_blank()) +
#   scale_size_continuous(name = '# Proben', trans = "log10", 
#                         range = c(2, 6), breaks = c(10,100, 1000)) +
#   labs(y = 'Wirkstoff', x = 'Anteil an Proben > UQN [%]') +
#   scale_color_manual(name = 'Gruppe', 
#                      values = c('#377eb8', '#4daf4a', '#B31010', '#984ea3'),
#                      labels = c('Fungizid', 'Herbizid', 'Insektizid', 'Metabolit')) +
#   theme_end +
#   facet_wrap(~type) +
#   guides(colour = guide_legend(override.aes = list(size = 5)))
# puqn
# # ggsave("/home/edisz/Documents/Uni/Projects/PHD/4BFG/Project/reports/20151131_Endbericht/chapters/figs/puqn.png", width = 10, height = 7,
# #        puqn)
# 
# 
# 
# ### ---------------------------------------
# ### RAK as Endpoint
# # restrict samples
# thrs <- psm_variables[!is.na(rak_uba), list(variable_id, name, cas, pgroup, rak_uba)]
# samples_thrs <- psm_samples[variable_id %in% thrs$variable]
# samples_thrs <- psm_sites_info[,list(site_id, type)][samples_thrs]
# samples_thrs <- samples_thrs[!is.na(type)]
# setkey(samples_thrs, variable_id)
# setkey(thrs, variable_id)
# samples_thrs <- thrs[samples_thrs]
# samples_thrs[ , g_rak := value_fin > rak_uba]
# 
# dd <- samples_thrs[!is.na(g_rak) , list(abs = sum(g_rak),
#                                                 prop = sum(g_rak) / length(g_rak),
#                                                 n = length(g_rak))
#                    , by = list(variable_id, type)][order(prop, decreasing = TRUE)]
# dd
# 
# setkey(dd, variable_id)
# dd <- psm_variables[ , list(variable_id, name, psm_type)][dd][order(prop, decreasing = TRUE)]
# 
# # restrict and order by appearance in sawb
# ord <- dd[type == 'klein, agrar'][order(prop, decreasing = TRUE), name]
# ord
# # take only fist 20 compounds
# ord <- ord[1:20]
# dd <- dd[name %chin% ord]
# dd$name <- factor(dd$name, levels = rev(ord))
# dd$type <- factor(dd$type, levels = unique(dd$type)[c(2, 3, 1, 4)])
# # # manually change type for irgarol
# # dd[variable_id == 375, psm_type := 'biocide']
# 
# prak <- ggplot() +
#   geom_point(data = dd, 
#              aes(x = prop * 100, y = name, size = n, col = psm_type)) +
#   theme(legend.key = element_blank()) +
#   scale_size_continuous(name = '# Proben', trans = "log10", 
#                         range = c(2, 6), breaks = c(10,100, 1000)) +
#   labs(y = 'Wirkstoff', x = 'Anteil an Proben > RAK [%]') +
#   scale_color_manual(name = 'Gruppe', 
#                        values = c('#377eb8', '#4daf4a', '#B31010'),
#                        labels = c('Fungizid', 'Herbizid', 'Insektizid')) +
#   theme_end +
#   facet_wrap(~type) +
#   guides(colour = guide_legend(override.aes = list(size = 5)))
# prak
# # ggsave("/home/edisz/Documents/Uni/Projects/PHD/4BFG/Project/reports/20151131_Endbericht/chapters/figs/prak.png", width = 10, height = 7,
# #        prak)
# 
# 
# 
# 
# ### how high are the RAK exceedances?
# rakü <- psm_samples[variable_id %in% dd[type == 'klein, agrar'][order(prop, decreasing = TRUE), variable_id][1:20] & value_fin > 0]
# setkey(rakü, variable_id)
# rakü <- psm_variables[ , list(variable_id, name, psm_type, rak_uba)][rakü]
# rakü[ , value_rak := value_fin / rak_uba]
# setkey(rakü, site_id)
# rakü <- psm_sites_info[,list(site_id, type)][rakü]
# rakü <- rakü[!is.na(type)]
# 
# prakü <- ggplot(rakü) +
#   geom_violin(aes(x = reorder(name, value_rak, FUN = median), y = value_rak, fill = psm_type)) +
#   scale_y_log10(breaks = c(0.01, 0.1, 1, 10, 100), limits = c(0.01, 150)) +
#   theme_end +
#   labs(x = '', y = 'Messwert / RAK') +
#   theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust = 1)) +
#   scale_fill_manual(name = 'Gruppe', 
#                       values = c('#377eb8', '#4daf4a', '#B31010'),
#                       labels = c('Fungizid', 'Herbizid', 'Insektizid')) +
#   geom_hline(yintercept = 1, linetype = 'dotted') +
#   coord_flip() +
#   guides(fill = guide_legend(override.aes = list(colour = NULL))) +  # remove diagonal line in legend
#   theme(legend.key = element_rect(colour = "black")) +
#   facet_wrap(~type) 
# prakü
# # ggsave("/home/edisz/Documents/Uni/Projects/PHD/4BFG/Project/reports/20151131_Endbericht/chapters/figs/prakü.png", width = 10, height = 7,
# #        prakü)
# 
# ### vorallem Neonics mit hohen Überschreitungen
# # kein so großer unterschied zwischen klein und groß...
# # Weils den Bach runterfließt?!
# 
# 
# 
# 
# ### ---------------------------------------
# ### Toxic Units as Endpoint
# 
# ### -------------------------------------------------------------------------
# ### TOXIC units as endpoint
# # add site & date information to maxtus
# setkey(psm_maxtu, variable_id)
# maxtu <- psm_variables[ , list(variable_id, name, psm_type)][psm_maxtu]
# setkey(maxtu, site_id)
# maxtu <- psm_sites_info[,list(site_id, type)][maxtu]
# maxtu <- maxtu[!is.na(type)]
# 
# 
# ### -------------------------
# # distribution of toxic units all sites
# # no_ex <- maxtu[list(abs = sum(log_maxtu > -2),
# #                     rel = sum(log_maxtu > -2) / length(log_maxtu)), 
# #                by = type]
# # no_ex
# 
# 
# tuhist <- ggplot() +
#   geom_histogram(data = maxtu[log_maxtu > -8], aes(x = log_maxtu), fill = 'grey75') +
#   geom_rug(data = maxtu[log_maxtu > -8], aes(x = log_maxtu)) +
#   facet_wrap(~type) +
#   theme_end +
#   labs(x = expression(log[10]*'('~TU[max]~')'), y = 'Anzahl Proben') +
#   geom_vline(data = maxtu[log_maxtu > -8], aes(x = log_maxtu), xintercept = -3, linetype = 'dotted') 
# tuhist
# # Proben mit Tumax < -8 nicht dargestellt
# # ggsave("/home/edisz/Documents/Uni/Projects/PHD/4BFG/Project/reports/20151131_Endbericht/chapters/figs/tuhist.png", width = 10, height = 7,
# #        tuhist)
# 
# 
# 
# freq <- as.data.table(table(maxtu[log_maxtu > -3, variable_id]))[order(N, decreasing = TRUE)][1:20, ]
# freq[ , variable_id := as.numeric(V1)]
# setkey(freq, variable_id)
# dd <- psm_variables[freq][ , list(variable_id, name, N, psm_type)][order(N, decreasing = TRUE)]
# 
# # Check DDT data
# tuc <- ggplot() +
#   geom_point(data = dd[1:10], 
#              aes(x = N, y = reorder(name, N), col = psm_type), size = 4) +
#   theme_end +
#   theme(legend.key = element_blank()) +
#   scale_size('# Proben', trans = "log10", range = c(2, 6)) +
#   labs(y = 'Wirkstoff', x = 'Anzahl maxTU') +
#   scale_color_manual(name = 'Gruppe', 
#                        values = c('#377eb8', '#4daf4a', '#B31010', '#984ea3'),
#                        labels = c('Fungizid', 'Herbizid', 'Insektizid', 'Metabolit'))
# tuc
# # wirkstoffe verantwortlich für tu > -3
# # ggsave("/home/edisz/Documents/Uni/Projects/PHD/4BFG/Project/reports/20151131_Endbericht/chapters/figs/tuc.png", 
# #        tuc, width = 5, height = 4.9)
# 
# 
# 
# ### ---------------------------------------
# ### Mischungen
# mm <- psm_samples[ , list(no_subs = sum(value_fin > 0)), by = list(sample_id, site_id, date)]
# setkey(mm, site_id)
# mm <- psm_sites_info[,list(site_id, type)][mm]
# mm <- mm[!is.na(type)]
# 
# mm_hist <- ggplot(mm, aes(x = no_subs)) +
#   geom_histogram(fill = 'grey75') +
#   facet_wrap(~type) +
#   theme_end +
#   labs(y = 'Anzahl Proben', x = 'Anzahl Stoffe in Mischung')
# mm_hist
# # ggsave("/home/edisz/Documents/Uni/Projects/PHD/4BFG/Project/reports/20151131_Endbericht/chapters/figs/mm_hist.png", 
# #        mm_hist, width = 7, height = 6.9)
# 
# 
# 
# ### ----------------------------------------------------------------------------
# ### Appendix tables
# ### --------------------------
# 
# 
# wlist <- fread(file.choose())
# 
# wlist$variable_id[!unique(wlist$variable_id) %in% unique(psm_samples$variable_id)]
# unique(psm_samples$variable_id)[!unique(psm_samples$variable_id) %in% unique(wlist$variable_id)]
# 
# w2 <- wlist[ , list(`HS-BfG`, Wirkstoff, bfg_id, variable_id)]
# w2 <- w2[2:nrow(w2)]
#   
# 
# 
# setkey(agg, variable_id)
# setkey(w2, variable_id)
# write.table(agg[w2], 'tmp.csv', sep = '|', row.names = FALSE)
# 
# 
# mlist <- fread(file.choose())
# 
# mlist$variable_id[!unique(mlist$variable_id) %in% unique(psm_samples$variable_id)]
# unique(psm_samples$variable_id)[!unique(psm_samples$variable_id) %in% unique(mlist$variable_id)]
# 
# m2 <- mlist[ , list(name, variable_id)]
# 
# agg <- psm_samples[ , list(nsamp = length(value_fin),
#                            gbg = sum(value_fin > 0),
#                            gbg_rel = sum(value_fin > 0) / length(value_fin),
#                            g01 = sum(value_fin > 0.1),
#                            g1 = sum(value_fin > 1)), by = variable_id]
# 
# setkey(agg, variable_id)
# setkey(m2, variable_id)
# write.table(agg[m2], 'tmp.csv', sep = '|', row.names = FALSE)
# 
# 
# setkey(psm_variables, variable_id)
# write.table(psm_variables[w2], 'tmp.csv', sep = '|', row.names = FALSE)
# 
# # Relate TU / RAK exceedance with number of compounds? (ala Malay, Moschet 2014)
# 
# 
# # number of substance measured
# no_subs_meas <- psm_samples[ , list(no_subs = length(value_fin)), by = list(sample_id, site_id, date)]
# range(no_subs_meas$no_subs)
# hist(no_subs_meas$no_subs)
# 
# psm_maxtu
# 
# 
# # merge maxtu & tox
# setkey(psm_maxtu, sample_id)
# setkey(no_subs_meas, sample_id)
# 
# nm <- no_subs_meas[psm_maxtu]
# plot(log_maxtu ~ no_subs, data = nm)
# 
# # what about LOQs and TUmax?
# 
# 
# 
# 
# # EZG -> number of compound / "dominanz of one compound" higher in small streams then in big streams
# vars <- unique(psm_samples$variable_id)
# 
# miss <- vars[!vars %in% vw$V1]
# 
# psm_variables[variable_id %in% miss]
# 
# # In wirkstoffliste but not in DB!
# vw$V1[!vw$V1 %in% vars]
# 
# # exclude lakes!
# vw <- read.table(text = '66
# 72
#                  76
#                  86
#                  116
#                  117
#                  131
#                  143
#                  247
#                  248
#                  249
#                  253
#                  275
#                  276
#                  278
#                  408
#                  409
#                  418
#                  423
#                  424
#                  471
#                  559
#                  676
#                  683
#                  739
#                  785
#                  825
#                  840
#                  841
#                  842
#                  843
#                  844
#                  845
#                  846
#                  847
#                  848
#                  849
#                  850
#                  1026
#                  1027
#                  1059
#                  1065
#                  1070
#                  1142
#                  1165
#                  1204
#                  1213
#                  1218
#                  118
# 120
#                  125
#                  142
#                  145
#                  146
#                  153
#                  166
#                  181
#                  186
#                  203
#                  597
#                  755
#                  206
#                  208
#                  213
#                  215
#                  232
#                  375
#                  1239
#                  123
#                  267
#                  292
#                  376
#                  237
#                  63
#                  1225
#                  485
#                  472
#                  483
#                  484
#                  242
#                  243
#                  244
#                  245
#                  246
#                  257
#                  263
#                  264
#                  266
#                  272
#                  279
#                  284
#                  285
#                  1224
#                  290
#                  291
#                  295
#                  306
#                  309
#                  314
#                  356
#                  357
#                  359
#                  642
#                  351
#                  352
#                  355
#                  353
#                  362
#                  378
#                  390
#                  397
#                  399
#                  401
#                  413
#                  407
#                  412
#                  422
#                  427
#                  429
#                  435
#                  476
#                  480
#                  487
#                  488
#                  499
#                  515
#                  517
#                  520
#                  524
#                  526
#                  528
#                  208
#                  540
#                  558
#                  70
#                  581
#                  580
#                  601
#                  608
#                  614
#                  867
#                  887
#                  868
#                  869
#                  870
#                  144
#                  871
#                  728
#                  872
#                  831
#                  690
#                  873
#                  874
#                  875
#                  147
#                  788
#                  151
#                  877
#                  878
#                  879
#                  167
#                  880
#                  180
#                  881
#                  702
#                  882
#                  191
#                  769
#                  192
#                  816
#                  703
#                  817
#                  651
#                  883
#                  884
#                  885
#                  851
#                  832
#                  886
#                  226
#                  227
#                  228
#                  234
#                  852
#                  706
#                  650
#                  888
#                  724
#                  238
#                  889
#                  648
#                  891
#                  647
#                  892
#                  893
#                  688
#                  686
#                  273
#                  274
#                  894
#                  280
#                  282
#                  890
#                  895
#                  896
#                  645
#                  897
#                  300
#                  299
#                  794
#                  708
#                  725
#                  819
#                  709
#                  308
#                  1120
#                  820
#                  312
#                  313
#                  733
#                  734
#                  899
#                  793
#                  900
#                  699
#                  681
#                  317
#                  901
#                  318
#                  838
#                  323
#                  821
#                  325
#                  789
#                  326
#                  327
#                  328
#                  839
#                  329
#                  861
#                  795
#                  822
#                  904
#                  710
#                  711
#                  905
#                  348
#                  349
#                  906
#                  713
#                  907
#                  908
#                  853
#                  823
#                  365
#                  714
#                  790
#                  374
#                  673
#                  726
#                  862
#                  379
#                  909
#                  384
#                  649
#                  387
#                  910
#                  715
#                  911
#                  1199
#                  834
#                  913
#                  824
#                  641
#                  915
#                  404
#                  916
#                  405
#                  406
#                  410
#                  695
#                  640
#                  738
#                  917
#                  414
#                  791
#                  716
#                  692
#                  919
#                  727
#                  452
#                  457
#                  920
#                  921
#                  496
#                  497
#                  498
#                  922
#                  505
#                  510
#                  923
#                  518
#                  924
#                  925
#                  523
#                  926
#                  525
#                  927
#                  836
#                  532
#                  718
#                  533
#                  928
#                  929
#                  930
#                  534
#                  931
#                  1220
#                  933
#                  536
#                  864
#                  539
#                  827
#                  693
#                  935
#                  937
#                  740
#                  565
#                  639
#                  828
#                  938
#                  638
#                  578
#                  742
#                  719
#                  747
#                  637
#                  939
#                  940
#                  720
#                  743
#                  588
#                  589
#                  837
#                  1215
#                  810
#                  941
#                  592
#                  594
#                  942
#                  600
#                  691
#                  829
#                  602
#                  792
#                  613
#                  694
#                  944
#                  746
#                  621
#                  860
#                  945
#                  721
#                  25
#                  809
#                  30
#                  64
#                  1124
#                  109
#                  115
#                  701
#                  678
#                  1216
#                  122
#                  960
#                  152
#                  1217
#                  185
#                  815
#                  1112
#                  1119
#                  729
#                  212
#                  1206
#                  1205
#                  1103
#                  233
#                  833
#                  705
#                  730
#                  252
#                  783
#                  672
#                  679
#                  277
#                  707
#                  675
#                  818
#                  646
#                  731
#                  732
#                  1212
#                  1121
#                  311
#                  644
#                  315
#                  865
#                  680
#                  643
#                  1022
#                  735
#                  813
#                  1023
#                  1025
#                  350
#                  358
#                  866
#                  712
#                  722
#                  854
#                  855
#                  856
#                  1122
#                  1102
#                  400
#                  736
#                  857
#                  737
#                  1024
#                  415
#                  426
#                  430
#                  766
#                  723
#                  1194
#                  1214
#                  991
#                  1054
#                  717
#                  478
#                  700
#                  863
#                  1195
#                  1219
#                  796
#                  835
#                  767
#                  527
#                  770
#                  784
#                  529
#                  1001
#                  1221
#                  687
#                  1002
#                  660
#                  677
#                  859
#                  552
#                  1004
#                  741
#                  1007
#                  1192
#                  1193
#                  656
#                  858
#                  749
#                  771
#                  744
#                  698
#                  745
#                  671
#                  898
#                  768
#                  316
#                  903
#                  1201
#                  1200
#                  1223
#                  918
#                  1198
#                  590
#                  830
#                  ')
# 
# 
# 
