cutoffs <- list(c(1,9), c(10,18), c(19,27), c(28,36),
             c(37,45), c(46,54), c(55,63), c(63,66))
for(i in 1:length(cutoffs)){
  str_list <- unique(tbl_plot$site)[cutoffs[[i]][1]:cutoffs[[i]][2]]
  print(tbl_plot[tbl_plot$site %in% str_list,] %>% 
    ggplot(aes(x = week, y = n, group = site))+
    geom_point()+
    xlim(15,32)+
    facet_wrap(~site)+
    theme_classic())
}