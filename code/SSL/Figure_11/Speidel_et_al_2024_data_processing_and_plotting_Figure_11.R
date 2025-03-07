library(readr)
library(plyr)
library(dplyr)
library(egg)
library(tidyverse)
library(igraph)
library(future)
library(future.apply)
library(magrittr)
library(readbulk)
library(ggpubr)
library(scales)
library(ggplot2)
library(ggh4x)
library(cowplot)
library(broom)
library(sf)
library(centr)
library(geosphere)
library(nngeo)

future::plan(multisession)

setwd('./Isolation_by_distance/stepping-stone_simulations/qpAdm/real_data/Speidel_et_al_2024_1_input_data/')

###reading the data, aggregated and slightly reformatted qpAdm results from the spreadsheets:
###https://www.biorxiv.org/content/biorxiv/early/2024/03/19/2024.03.15.585102/DC1/embed/media-1.xlsx?download=true
###the results are stratified by target group name, model label (model), proxy source name (source), 
###and model complexity (n_sources).
df <- read_delim("Speidel_results.txt", delim="\t", trim_ws=T)

df <- df[!duplicated(df[, c("target", "model", "source", "n_sources")]), ]
df1 <- filter(df, n_sources==1)
df2 <- filter(df, n_sources==2)
df3 <- filter(df, n_sources==3)

###1-way models:
outtab1 <- df1[,c(1,7,4,5,6,3)]
names(outtab1)[names(outtab1) == 'source'] <- 'source1'
names(outtab1)[names(outtab1) == 'EAF'] <- 'EAF1'
names(outtab1)[names(outtab1) == 'SE'] <- 'SE1'
write.table(outtab1, file="outtab1.txt", row.names=F, sep='\t')

###2-way models:
models <- (df2[!duplicated(df2[, c("target", "model")]), ])[,c(1:2)] %>% as.data.frame() %>% set_rownames(NULL)
combs <- c(1:length(models$target))
outtab2 = sapply(combs, function(x){
  test <- filter(df2, target==models[x,1], model==models[x,2])
  c("target"=as.vector(models[x,1]),
    "n_sources" = 2,
    "source1" = test$source[1],
    "source2" = test$source[2],
    "EAF1" = test$EAF[1],
    "EAF2" = test$EAF[2],
    "SE1" = test$SE[1],
    "SE2" = test$SE[2],
    "p" = test$p[1])
}) %>% t() %>% as.data.frame() %>% set_rownames(NULL)
outtab2 <- outtab2 %>% mutate(target = as.character(target))
outtab2 <- outtab2 %>% mutate(n_sources = as.numeric(n_sources))
outtab2 <- outtab2 %>% mutate(source1 = as.character(source1))
outtab2 <- outtab2 %>% mutate(source2 = as.character(source2))
outtab2 <- outtab2 %>% mutate(EAF1 = as.numeric(EAF1))
outtab2 <- outtab2 %>% mutate(EAF2 = as.numeric(EAF2))
outtab2 <- outtab2 %>% mutate(SE1 = as.numeric(SE1))
outtab2 <- outtab2 %>% mutate(SE2 = as.numeric(SE2))
outtab2 <- outtab2 %>% mutate(p = as.numeric(p))
write.table(outtab2, file="outtab2.txt", row.names=F, sep='\t')

###3-way models:
models <- (df3[!duplicated(df3[, c("target", "model")]), ])[,c(1:2)] %>% as.data.frame() %>% set_rownames(NULL)
combs <- c(1:length(models$target))
outtab3 = sapply(combs, function(x){
  test <- filter(df3, target==models[x,1], model==models[x,2])
  c("target"=as.vector(models[x,1]),
    "n_sources" = 3,
    "source1" = test$source[1],
    "source2" = test$source[2],
    "source3" = test$source[3],
    "EAF1" = test$EAF[1],
    "EAF2" = test$EAF[2],
    "EAF3" = test$EAF[3],
    "SE1" = test$SE[1],
    "SE2" = test$SE[2],
    "SE3" = test$SE[3],
    "p" = test$p[1])
}) %>% t() %>% as.data.frame() %>% set_rownames(NULL)
outtab3 <- outtab3 %>% mutate(target = as.character(target))
outtab3 <- outtab3 %>% mutate(n_sources = as.numeric(n_sources))
outtab3 <- outtab3 %>% mutate(source1 = as.character(source1))
outtab3 <- outtab3 %>% mutate(source2 = as.character(source2))
outtab3 <- outtab3 %>% mutate(source3 = as.character(source3))
outtab3 <- outtab3 %>% mutate(EAF1 = as.numeric(EAF1))
outtab3 <- outtab3 %>% mutate(EAF2 = as.numeric(EAF2))
outtab3 <- outtab3 %>% mutate(EAF3 = as.numeric(EAF3))
outtab3 <- outtab3 %>% mutate(SE1 = as.numeric(SE1))
outtab3 <- outtab3 %>% mutate(SE2 = as.numeric(SE2))
outtab3 <- outtab3 %>% mutate(SE3 = as.numeric(SE3))
outtab3 <- outtab3 %>% mutate(p = as.numeric(p))
write.table(outtab3, file="outtab3.txt", row.names=F, sep='\t')

rm(df1)
rm(df2)
rm(df3)


###calculation of group centroids on the landscape, average and median dates for groups###
###this is group composition information and geographic coordinates from:
###https://www.biorxiv.org/content/biorxiv/early/2024/03/19/2024.03.15.585102/DC1/embed/media-1.xlsx?download=true
df <- read_delim(file="Speidel_geographic_coordinates_dates.txt", delim="\t", trim_ws=T)

demes <- (df[!duplicated(df[, "qpAdm_label"]), ])[,1] %>% as.data.frame() %>% set_rownames(NULL)
combs <- c(1:length(demes$qpAdm_label))
outtab = sapply(combs, function(x){
  dates <- filter(df, qpAdm_label==demes[x,1])[,3]
  pts <- filter(df, qpAdm_label==demes[x,1])[,c(4:5)] %>% st_as_sf(coords = c("longitude", "latitude"), crs = 4326)
  center_mean <- mean_center(pts) %>% st_coordinates() %>% as.data.frame()
  c("qpAdm_label"=demes[x,1], 
    "average.date" = mean(dates$date_yBP),
    "median.date" = median(dates$date_yBP),
    "mean_center_latitude" = center_mean$Y[1],
    "mean_center_longitude" = center_mean$X[1])
}) %>% t() %>% as.data.frame() %>% set_rownames(NULL)
write.table(outtab, file="Speidel_geographic_centroids_dates.txt", row.names=F, sep='\t')


###calculating max|EAF-EF|###
outtab2 <- read_delim(file="outtab2.txt", delim='\t')
outtab3 <- read_delim(file="outtab3.txt", delim='\t')

###2-way models:
combs <- c(1:length(outtab2$target))
out2 = sapply(combs, function(x){
  max_deviation_fromEF = max(c(abs(outtab2$EAF1[x]-1/2), abs(outtab2$EAF2[x]-1/2)))
  max_se = outtab2$SE1[x]
  c("max_deviation_fromEF" = max_deviation_fromEF, "max_se" = max_se)
}) %>% t() %>% as.data.frame() %>% set_rownames(NULL)
colnames(out2) <- c("max_deviation_fromEF", "max_se")
outtab2 <- cbind(outtab2, out2)

###3-way models:
combs <- c(1:length(outtab3$target))
out3 = sapply(combs, function(x){
  max_deviation_fromEF = max(c(abs(outtab3$EAF1[x]-1/3), abs(outtab3$EAF2[x]-1/3), abs(outtab3$EAF3[x]-1/3)))
  max_se = max(c(outtab3$SE1[x], outtab3$SE2[x], outtab3$SE3[x]), na.rm=T)
  c("max_deviation_fromEF" = max_deviation_fromEF, "max_se" = max_se)
}) %>% t() %>% as.data.frame() %>% set_rownames(NULL)
colnames(out3) <- c("max_deviation_fromEF", "max_se")
outtab3 <- cbind(outtab3, out3)


###adding the geographic centroids and average/median dates to the qpAdm results###
###2-way models:
df <- read_delim(file="Speidel_geographic_centroids_dates.txt", delim="\t", trim_ws=T)
names(df)[names(df)=='qpAdm_label'] <- 'target'
outtab2 <- left_join(outtab2, df, by="target")
names(outtab2)[names(outtab2)=='average.date'] <- 'average.date.t'
names(outtab2)[names(outtab2)=='median.date'] <- 'median.date.t'
names(outtab2)[names(outtab2)=='mean_center_latitude'] <- 'mean_center_latitude.t'
names(outtab2)[names(outtab2)=='mean_center_longitude'] <- 'mean_center_longitude.t'
names(df)[names(df)=='target'] <- 'source1'
outtab2 <- left_join(outtab2, df, by="source1")
names(outtab2)[names(outtab2)=='average.date'] <- 'average.date.s1'
names(outtab2)[names(outtab2)=='median.date'] <- 'median.date.s1'
names(outtab2)[names(outtab2)=='mean_center_latitude'] <- 'mean_center_latitude.s1'
names(outtab2)[names(outtab2)=='mean_center_longitude'] <- 'mean_center_longitude.s1'
names(df)[names(df)=='source1'] <- 'source2'
outtab2 <- left_join(outtab2, df, by="source2")
names(outtab2)[names(outtab2)=='average.date'] <- 'average.date.s2'
names(outtab2)[names(outtab2)=='median.date'] <- 'median.date.s2'
names(outtab2)[names(outtab2)=='mean_center_latitude'] <- 'mean_center_latitude.s2'
names(outtab2)[names(outtab2)=='mean_center_longitude'] <- 'mean_center_longitude.s2'

###3-way models:
df <- read_delim(file="Speidel_geographic_centroids_dates.txt", delim="\t", trim_ws=T)
names(df)[names(df)=='qpAdm_label'] <- 'target'
outtab3 <- left_join(outtab3, df, by="target")
names(outtab3)[names(outtab3)=='average.date'] <- 'average.date.t'
names(outtab3)[names(outtab3)=='median.date'] <- 'median.date.t'
names(outtab3)[names(outtab3)=='mean_center_latitude'] <- 'mean_center_latitude.t'
names(outtab3)[names(outtab3)=='mean_center_longitude'] <- 'mean_center_longitude.t'
names(df)[names(df)=='target'] <- 'source1'
outtab3 <- left_join(outtab3, df, by="source1")
names(outtab3)[names(outtab3)=='average.date'] <- 'average.date.s1'
names(outtab3)[names(outtab3)=='median.date'] <- 'median.date.s1'
names(outtab3)[names(outtab3)=='mean_center_latitude'] <- 'mean_center_latitude.s1'
names(outtab3)[names(outtab3)=='mean_center_longitude'] <- 'mean_center_longitude.s1'
names(df)[names(df)=='source1'] <- 'source2'
outtab3 <- left_join(outtab3, df, by="source2")
names(outtab3)[names(outtab3)=='average.date'] <- 'average.date.s2'
names(outtab3)[names(outtab3)=='median.date'] <- 'median.date.s2'
names(outtab3)[names(outtab3)=='mean_center_latitude'] <- 'mean_center_latitude.s2'
names(outtab3)[names(outtab3)=='mean_center_longitude'] <- 'mean_center_longitude.s2'
names(df)[names(df)=='source2'] <- 'source3'
outtab3 <- left_join(outtab3, df, by="source3")
names(outtab3)[names(outtab3)=='average.date'] <- 'average.date.s3'
names(outtab3)[names(outtab3)=='median.date'] <- 'median.date.s3'
names(outtab3)[names(outtab3)=='mean_center_latitude'] <- 'mean_center_latitude.s3'
names(outtab3)[names(outtab3)=='mean_center_longitude'] <- 'mean_center_longitude.s3'


###classifying models into distal and proximal, all dates rounded to the nearest 100###
outtab2$temp_strat_average <- ifelse(round_any(outtab2$average.date.s1, 100)<round_any(outtab2$average.date.t, 100) | round_any(outtab2$average.date.s2, 100)<round_any(outtab2$average.date.t, 100), "proximal", "distal")
outtab2$temp_strat_median <- ifelse(round_any(outtab2$median.date.s1, 100)<round_any(outtab2$median.date.t, 100) | round_any(outtab2$median.date.s2, 100)<round_any(outtab2$median.date.t, 100), "proximal", "distal")
outtab3$temp_strat_average <- ifelse(round_any(outtab3$average.date.s1, 100)<round_any(outtab3$average.date.t, 100) | round_any(outtab3$average.date.s2, 100)<round_any(outtab3$average.date.t, 100) | round_any(outtab3$average.date.s3, 100)<round_any(outtab3$average.date.t, 100), "proximal", "distal")
outtab3$temp_strat_median <- ifelse(round_any(outtab3$median.date.s1, 100)<round_any(outtab3$median.date.t, 100) | round_any(outtab3$median.date.s2, 100)<round_any(outtab3$median.date.t, 100) | round_any(outtab3$median.date.s3, 100)<round_any(outtab3$median.date.t, 100), "proximal", "distal")


###calculating great circle distances between target and source groups in a model###
###2-way models:
models <- as.data.frame(cbind(outtab2$mean_center_longitude.t, outtab2$mean_center_latitude.t, outtab2$mean_center_longitude.s1, outtab2$mean_center_latitude.s1, outtab2$mean_center_longitude.s2, outtab2$mean_center_latitude.s2))
colnames(models) <- c("x.t", "y.t", "x.s1", "y.s1", "x.s2", "y.s2")

res1 = future_apply(models, 1, function(x) {
  out <- distCosine(c(x[1], x[2]), c(x[3], x[4]))
})
res2 = future_apply(models, 1, function(x) {
  out <- distCosine(c(x[1], x[2]), c(x[5], x[6]))
})

combs <- c(1:length(models$x.t))
outtab = sapply(combs, function(comb){
  t_s1_dist = res1[[comb]][[1]]/1000
  t_s2_dist = res2[[comb]][[1]]/1000
  c("distance_t.s1" = t_s1_dist, "distance_t.s2" = t_s2_dist, 
    "max.distance" = max(t_s1_dist, t_s2_dist), "avg.distance" = (t_s1_dist + t_s2_dist)/2)
}) %>% t() %>% as.data.frame() %>% set_rownames(NULL)
outtab2 <- cbind(outtab2, outtab)
write.table(outtab2, file="outtab2_merged_v2.txt", row.names=F, sep='\t')


###3-way models:
models <- as.data.frame(cbind(outtab3$mean_center_longitude.t, outtab3$mean_center_latitude.t, outtab3$mean_center_longitude.s1, outtab3$mean_center_latitude.s1, outtab3$mean_center_longitude.s2, outtab3$mean_center_latitude.s2, outtab3$mean_center_longitude.s3, outtab3$mean_center_latitude.s3))
colnames(models) <- c("x.t", "y.t", "x.s1", "y.s1", "x.s2", "y.s2", "x.s3", "y.s3")

res1 = future_apply(models, 1, function(x) {
  out <- distCosine(c(x[1], x[2]), c(x[3], x[4]))
})
res2 = future_apply(models, 1, function(x) {
  out <- distCosine(c(x[1], x[2]), c(x[5], x[6]))
})
res3 = future_apply(models, 1, function(x) {
  out <- distCosine(c(x[1], x[2]), c(x[7], x[8]))
})

combs <- c(1:length(models$x.t))
outtab = sapply(combs, function(comb){
  t_s1_dist = res1[[comb]][[1]]/1000
  t_s2_dist = res2[[comb]][[1]]/1000
  t_s3_dist = res3[[comb]][[1]]/1000
  c("distance_t.s1" = t_s1_dist, "distance_t.s2" = t_s2_dist, "distance_t.s3" = t_s3_dist, 
    "max.distance" = max(t_s1_dist, t_s2_dist, t_s3_dist), "avg.distance" = (t_s1_dist + t_s2_dist + t_s3_dist)/3)
}) %>% t() %>% as.data.frame() %>% set_rownames(NULL)
outtab3 <- cbind(outtab3, outtab)
write.table(outtab3, file="outtab3_merged_v2.txt", row.names=F, sep='\t')


###calculating azimuths between two courses: min. source-target-source angles (min. STS angles)###
angle_diff <- function(theta1, theta2){
  theta <- abs(theta1 - theta2) %% 360 
  return(ifelse(theta > 180, 360 - theta, theta))
}

###calculating angles (if min. source-target distance in a model >50 km!!!)###
###2-way models:
combs <- c(1:length(outtab2$target))
outtab = sapply(combs, function(x){
  if (outtab2$distance_t.s1[x]<=50 | outtab2$distance_t.s2[x]<=50 ) {
    "min_S1_T_S2_angle" = NA
  } else {
    A <- c(outtab2$mean_center_longitude.t[x], outtab2$mean_center_latitude.t[x])
    B <- c(outtab2$mean_center_longitude.s1[x], outtab2$mean_center_latitude.s1[x])
    C <- c(outtab2$mean_center_longitude.s2[x], outtab2$mean_center_latitude.s2[x])
    S1_T_S2_angle <- angle_diff(bearing(A, B, a=6378137, f=1/298.257223563), bearing(A, C, a=6378137, f=1/298.257223563))
    if (S1_T_S2_angle > 180) {"min_S1_T_S2_angle" = 360 - S1_T_S2_angle}
    else {"min_S1_T_S2_angle" = S1_T_S2_angle}
  } 
}) %>% as.data.frame() %>% set_rownames(NULL)
colnames(outtab) <- c("min_S1_T_S2_angle")
outtab2 <- cbind(outtab2, outtab)

###3-way models:
combs <- c(1:length(outtab3$target))
outtab = sapply(combs, function(x){
  if (outtab3$distance_t.s1[x]<=50 | outtab3$distance_t.s2[x]<=50 | outtab3$distance_t.s3[x]<=50 ) {
    "min_S1_T_S2_angle" = NA}
  else {
    A <- c(outtab3$mean_center_longitude.t[x], outtab3$mean_center_latitude.t[x])
    B <- c(outtab3$mean_center_longitude.s1[x], outtab3$mean_center_latitude.s1[x])
    C <- c(outtab3$mean_center_longitude.s2[x], outtab3$mean_center_latitude.s2[x])
    D <- c(outtab3$mean_center_longitude.s3[x], outtab3$mean_center_latitude.s3[x])
    S1_T_S2_angle1 <- angle_diff(bearing(A, B, a=6378137, f=1/298.257223563), bearing(A, C, a=6378137, f=1/298.257223563))
    if (S1_T_S2_angle1 < 180) {}
    else {S1_T_S2_angle1 = 360-S1_T_S2_angle1}
    S1_T_S2_angle2 <- angle_diff(bearing(A, B, a=6378137, f=1/298.257223563), bearing(A, D, a=6378137, f=1/298.257223563))
    if (S1_T_S2_angle2 < 180) {}
    else {S1_T_S2_angle2 = 360-S1_T_S2_angle2}
    S1_T_S2_angle3 <- angle_diff(bearing(A, C, a=6378137, f=1/298.257223563), bearing(A, D, a=6378137, f=1/298.257223563))
    if (S1_T_S2_angle3 < 180) {}
    else {S1_T_S2_angle3 = 360-S1_T_S2_angle3}
    "min_S1_T_S2_angle" = min(S1_T_S2_angle1, S1_T_S2_angle2, S1_T_S2_angle3)
  } 
}) %>% as.data.frame() %>% set_rownames(NULL)
colnames(outtab) <- c("min_S1_T_S2_angle")
outtab3 <- cbind(outtab3, outtab)


###calculate distance to the ideal symmetric models in the space of "model optimality metrics" 
###(max./average ST distances and min. STS angles)###
euclidean_dist <- function(x, y) sqrt(sum((x - y)^2))

###2-way models:
combs <- c(1:length(outtab2$target))
outtab = sapply(combs, function(x){
  if (is.na(outtab2$min_S1_T_S2_angle[x])) {
    c("dist_to_corner_maxdist" = NA, "dist_to_corner_avgdist" = NA)
  } else {
    c("dist_to_corner_maxdist" = euclidean_dist(c(outtab2$max.distance[x]/17.4, outtab2$min_S1_T_S2_angle[x]), c(50/17.4,180)),
      "dist_to_corner_avgdist" = euclidean_dist(c(outtab2$avg.distance[x]/17.4, outtab2$min_S1_T_S2_angle[x]), c(50/17.4,180)))
  } 
}) %>% t() %>% as.data.frame() %>% set_rownames(NULL)
colnames(outtab) <- c("dist_to_corner_maxdist", "dist_to_corner_avgdist")
outtab2 <- cbind(outtab2, outtab)

###3-way models:
combs <- c(1:length(outtab3$target))
outtab = sapply(combs, function(x){
  if (is.na(outtab3$min_S1_T_S2_angle[x])) {
    c("dist_to_corner_maxdist" = NA, "dist_to_corner_avgdist" = NA)
  } else {
    c("dist_to_corner_maxdist" = euclidean_dist(c(outtab3$max.distance[x]/17.4, outtab3$min_S1_T_S2_angle[x]*180/120), c(50/17.4,180)),
      "dist_to_corner_avgdist" = euclidean_dist(c(outtab3$avg.distance[x]/17.4, outtab3$min_S1_T_S2_angle[x]*180/120), c(50/17.4,180)))
  } 
}) %>% t() %>% as.data.frame() %>% set_rownames(NULL)
colnames(outtab) <- c("dist_to_corner_maxdist", "dist_to_corner_avgdist")
outtab3 <- cbind(outtab3, outtab)

rm(out2)
rm(out3)
rm(combs)
rm(models)
rm(outtab)


###adding model feasibility criteria###
###2-way models:
outtab2$positive22 <- ifelse(outtab2$EAF1-2*outtab2$SE1 >0 & outtab2$EAF1+2*outtab2$SE1 >0 & outtab2$EAF2-2*outtab2$SE2 >0 & outtab2$EAF2+2*outtab2$SE2 >0 & 
                               outtab2$EAF1-2*outtab2$SE1 <1 & outtab2$EAF1+2*outtab2$SE1 <1 & outtab2$EAF2-2*outtab2$SE2 <1 & outtab2$EAF2+2*outtab2$SE2 <1 & 
                               outtab2$p >=0.001, "fitting", "non-fitting")
outtab2$positive23 <- ifelse(outtab2$EAF1-2*outtab2$SE1 >0 & outtab2$EAF1+2*outtab2$SE1 >0 & outtab2$EAF2-2*outtab2$SE2 >0 & outtab2$EAF2+2*outtab2$SE2 >0 & 
                               outtab2$EAF1-2*outtab2$SE1 <1 & outtab2$EAF1+2*outtab2$SE1 <1 & outtab2$EAF2-2*outtab2$SE2 <1 & outtab2$EAF2+2*outtab2$SE2 <1 & 
                               outtab2$p >=0.01, "fitting", "non-fitting")
outtab2$positive24 <- ifelse(outtab2$EAF1-2*outtab2$SE1 >0 & outtab2$EAF1+2*outtab2$SE1 >0 & outtab2$EAF2-2*outtab2$SE2 >0 & outtab2$EAF2+2*outtab2$SE2 >0 & 
                               outtab2$EAF1-2*outtab2$SE1 <1 & outtab2$EAF1+2*outtab2$SE1 <1 & outtab2$EAF2-2*outtab2$SE2 <1 & outtab2$EAF2+2*outtab2$SE2 <1 & 
                               outtab2$p >=0.05, "fitting", "non-fitting")
outtab2$positive25 <- ifelse(outtab2$EAF1-2*outtab2$SE1 >0 & outtab2$EAF1+2*outtab2$SE1 >0 & outtab2$EAF2-2*outtab2$SE2 >0 & outtab2$EAF2+2*outtab2$SE2 >0 & 
                               outtab2$EAF1-2*outtab2$SE1 <1 & outtab2$EAF1+2*outtab2$SE1 <1 & outtab2$EAF2-2*outtab2$SE2 <1 & outtab2$EAF2+2*outtab2$SE2 <1 & 
                               outtab2$p >=0.1, "fitting", "non-fitting")
outtab2$positive26 <- ifelse(outtab2$EAF1-2*outtab2$SE1 >0 & outtab2$EAF1+2*outtab2$SE1 >0 & outtab2$EAF2-2*outtab2$SE2 >0 & outtab2$EAF2+2*outtab2$SE2 >0 & 
                               outtab2$EAF1-2*outtab2$SE1 <1 & outtab2$EAF1+2*outtab2$SE1 <1 & outtab2$EAF2-2*outtab2$SE2 <1 & outtab2$EAF2+2*outtab2$SE2 <1 & 
                               outtab2$p >=0.5, "fitting", "non-fitting")

outtab2$positive27 <- ifelse(outtab2$EAF1<1 & outtab2$EAF2<1 & outtab2$EAF1>0 & outtab2$EAF2>0 & outtab2$p >=0.001, "fitting", "non-fitting")
outtab2$positive28 <- ifelse(outtab2$EAF1<1 & outtab2$EAF2<1 & outtab2$EAF1>0 & outtab2$EAF2>0 & outtab2$p >=0.01, "fitting", "non-fitting")
outtab2$positive29 <- ifelse(outtab2$EAF1<1 & outtab2$EAF2<1 & outtab2$EAF1>0 & outtab2$EAF2>0 & outtab2$p >=0.05, "fitting", "non-fitting")
outtab2$positive30 <- ifelse(outtab2$EAF1<1 & outtab2$EAF2<1 & outtab2$EAF1>0 & outtab2$EAF2>0 & outtab2$p >=0.1, "fitting", "non-fitting")
outtab2$positive31 <- ifelse(outtab2$EAF1<1 & outtab2$EAF2<1 & outtab2$EAF1>0 & outtab2$EAF2>0 & outtab2$p >=0.5, "fitting", "non-fitting")

outtab2$positive32 <- ifelse(outtab2$EAF1<1.3 & outtab2$EAF2<1.3 & outtab2$EAF1>-0.3 & outtab2$EAF2>-0.3 & outtab2$p >=0.001, "fitting", "non-fitting")
outtab2$positive33 <- ifelse(outtab2$EAF1<1.3 & outtab2$EAF2<1.3 & outtab2$EAF1>-0.3 & outtab2$EAF2>-0.3 & outtab2$p >=0.01, "fitting", "non-fitting")
outtab2$positive34 <- ifelse(outtab2$EAF1<1.3 & outtab2$EAF2<1.3 & outtab2$EAF1>-0.3 & outtab2$EAF2>-0.3 & outtab2$p >=0.05, "fitting", "non-fitting")
outtab2$positive35 <- ifelse(outtab2$EAF1<1.3 & outtab2$EAF2<1.3 & outtab2$EAF1>-0.3 & outtab2$EAF2>-0.3 & outtab2$p >=0.1, "fitting", "non-fitting")
outtab2$positive36 <- ifelse(outtab2$EAF1<1.3 & outtab2$EAF2<1.3 & outtab2$EAF1>-0.3 & outtab2$EAF2>-0.3 & outtab2$p >=0.5, "fitting", "non-fitting")

###3-way models:
outtab3$positive22 <- ifelse(outtab3$EAF1-2*outtab3$SE1 >0 & outtab3$EAF1+2*outtab3$SE1 >0 & outtab3$EAF2-2*outtab3$SE2 >0 & outtab3$EAF2+2*outtab3$SE2 >0 & outtab3$EAF3-2*outtab3$SE3 >0 & outtab3$EAF3+2*outtab3$SE3 >0 & 
                               outtab3$EAF1-2*outtab3$SE1 <1 & outtab3$EAF1+2*outtab3$SE1 <1 & outtab3$EAF2-2*outtab3$SE2 <1 & outtab3$EAF2+2*outtab3$SE2 <1 & outtab3$EAF3-2*outtab3$SE3 <1 & outtab3$EAF3+2*outtab3$SE3 <1 & 
                               outtab3$p >=0.001, "fitting", "non-fitting")
outtab3$positive23 <- ifelse(outtab3$EAF1-2*outtab3$SE1 >0 & outtab3$EAF1+2*outtab3$SE1 >0 & outtab3$EAF2-2*outtab3$SE2 >0 & outtab3$EAF2+2*outtab3$SE2 >0 & outtab3$EAF3-2*outtab3$SE3 >0 & outtab3$EAF3+2*outtab3$SE3 >0 & 
                               outtab3$EAF1-2*outtab3$SE1 <1 & outtab3$EAF1+2*outtab3$SE1 <1 & outtab3$EAF2-2*outtab3$SE2 <1 & outtab3$EAF2+2*outtab3$SE2 <1 & outtab3$EAF3-2*outtab3$SE3 <1 & outtab3$EAF3+2*outtab3$SE3 <1 & 
                               outtab3$p >=0.01, "fitting", "non-fitting")
outtab3$positive24 <- ifelse(outtab3$EAF1-2*outtab3$SE1 >0 & outtab3$EAF1+2*outtab3$SE1 >0 & outtab3$EAF2-2*outtab3$SE2 >0 & outtab3$EAF2+2*outtab3$SE2 >0 & outtab3$EAF3-2*outtab3$SE3 >0 & outtab3$EAF3+2*outtab3$SE3 >0 & 
                               outtab3$EAF1-2*outtab3$SE1 <1 & outtab3$EAF1+2*outtab3$SE1 <1 & outtab3$EAF2-2*outtab3$SE2 <1 & outtab3$EAF2+2*outtab3$SE2 <1 & outtab3$EAF3-2*outtab3$SE3 <1 & outtab3$EAF3+2*outtab3$SE3 <1 & 
                               outtab3$p >=0.05, "fitting", "non-fitting")
outtab3$positive25 <- ifelse(outtab3$EAF1-2*outtab3$SE1 >0 & outtab3$EAF1+2*outtab3$SE1 >0 & outtab3$EAF2-2*outtab3$SE2 >0 & outtab3$EAF2+2*outtab3$SE2 >0 & outtab3$EAF3-2*outtab3$SE3 >0 & outtab3$EAF3+2*outtab3$SE3 >0 & 
                               outtab3$EAF1-2*outtab3$SE1 <1 & outtab3$EAF1+2*outtab3$SE1 <1 & outtab3$EAF2-2*outtab3$SE2 <1 & outtab3$EAF2+2*outtab3$SE2 <1 & outtab3$EAF3-2*outtab3$SE3 <1 & outtab3$EAF3+2*outtab3$SE3 <1 & 
                               outtab3$p >=0.1, "fitting", "non-fitting")
outtab3$positive26 <- ifelse(outtab3$EAF1-2*outtab3$SE1 >0 & outtab3$EAF1+2*outtab3$SE1 >0 & outtab3$EAF2-2*outtab3$SE2 >0 & outtab3$EAF2+2*outtab3$SE2 >0 & outtab3$EAF3-2*outtab3$SE3 >0 & outtab3$EAF3+2*outtab3$SE3 >0 & 
                               outtab3$EAF1-2*outtab3$SE1 <1 & outtab3$EAF1+2*outtab3$SE1 <1 & outtab3$EAF2-2*outtab3$SE2 <1 & outtab3$EAF2+2*outtab3$SE2 <1 & outtab3$EAF3-2*outtab3$SE3 <1 & outtab3$EAF3+2*outtab3$SE3 <1 & 
                               outtab3$p >=0.5, "fitting", "non-fitting")

outtab3$positive27 <- ifelse(outtab3$EAF1<1 & outtab3$EAF2<1 & outtab3$EAF3<1 & outtab3$EAF1>0 & outtab3$EAF2>0 & outtab3$EAF3>0 & outtab3$p >=0.001, "fitting", "non-fitting")
outtab3$positive28 <- ifelse(outtab3$EAF1<1 & outtab3$EAF2<1 & outtab3$EAF3<1 & outtab3$EAF1>0 & outtab3$EAF2>0 & outtab3$EAF3>0 & outtab3$p >=0.01, "fitting", "non-fitting")
outtab3$positive29 <- ifelse(outtab3$EAF1<1 & outtab3$EAF2<1 & outtab3$EAF3<1 & outtab3$EAF1>0 & outtab3$EAF2>0 & outtab3$EAF3>0 & outtab3$p >=0.05, "fitting", "non-fitting")
outtab3$positive30 <- ifelse(outtab3$EAF1<1 & outtab3$EAF2<1 & outtab3$EAF3<1 & outtab3$EAF1>0 & outtab3$EAF2>0 & outtab3$EAF3>0 & outtab3$p >=0.1, "fitting", "non-fitting")
outtab3$positive31 <- ifelse(outtab3$EAF1<1 & outtab3$EAF2<1 & outtab3$EAF3<1 & outtab3$EAF1>0 & outtab3$EAF2>0 & outtab3$EAF3>0 & outtab3$p >=0.5, "fitting", "non-fitting")

outtab3$positive32 <- ifelse(outtab3$EAF1<1.3 & outtab3$EAF2<1.3 & outtab3$EAF3<1.3 & outtab3$EAF1>-0.3 & outtab3$EAF2>-0.3 & outtab3$EAF3>-0.3 & outtab3$p >=0.001, "fitting", "non-fitting")
outtab3$positive33 <- ifelse(outtab3$EAF1<1.3 & outtab3$EAF2<1.3 & outtab3$EAF3<1.3 & outtab3$EAF1>-0.3 & outtab3$EAF2>-0.3 & outtab3$EAF3>-0.3 & outtab3$p >=0.01, "fitting", "non-fitting")
outtab3$positive34 <- ifelse(outtab3$EAF1<1.3 & outtab3$EAF2<1.3 & outtab3$EAF3<1.3 & outtab3$EAF1>-0.3 & outtab3$EAF2>-0.3 & outtab3$EAF3>-0.3 & outtab3$p >=0.05, "fitting", "non-fitting")
outtab3$positive35 <- ifelse(outtab3$EAF1<1.3 & outtab3$EAF2<1.3 & outtab3$EAF3<1.3 & outtab3$EAF1>-0.3 & outtab3$EAF2>-0.3 & outtab3$EAF3>-0.3 & outtab3$p >=0.1, "fitting", "non-fitting")
outtab3$positive36 <- ifelse(outtab3$EAF1<1.3 & outtab3$EAF2<1.3 & outtab3$EAF3<1.3 & outtab3$EAF1>-0.3 & outtab3$EAF2>-0.3 & outtab3$EAF3>-0.3 & outtab3$p >=0.5, "fitting", "non-fitting")

write.table(outtab2, file="outtab2_merged_v3.txt", row.names=F, sep='\t')
write.table(outtab3, file="outtab3_merged_v3.txt", row.names=F, sep='\t')


###merging results for all model complexities into one object###
outtab2 <- read_delim(file="./Isolation_by_distance/stepping-stone_simulations/qpAdm/real_data/Speidel_et_al_2024_1_input_data/outtab2_merged_v3.txt", delim='\t')
outtab3 <- read_delim(file="./Isolation_by_distance/stepping-stone_simulations/qpAdm/real_data/Speidel_et_al_2024_1_input_data/outtab3_merged_v3.txt", delim='\t')

speidel <- rbind.fill(outtab2, outtab3)
#min(filter(speidel, p!=0)$p)
speidel$log_p <- log10(speidel$p)
speidel$log_p[speidel$log_p == -Inf] <- -304
speidel2 <- speidel
#speidel2 <- speidel[,c(2,1,3,28,29,30,48,45,43,44,47,46,33,34,6,56,39,40,4,31,49,5,32,50,21,41,42,55,7,11,35,51,8,12,36,52,9,13,37,53,10,14,38,54,15,16,17,18,19,20,22,23,24)]
#speidel2 <- speidel2 %>% relocate(positive28, .after=positive27)
speidel2$n_sources <- gsub("2", "2-way", speidel2$n_sources, fixed=T)
speidel2$n_sources <- gsub("3", "3-way", speidel2$n_sources, fixed=T)


###binning min. STS angles###
speidel2<-filter(speidel2, !is.na(min_S1_T_S2_angle))
speidel2$min_S1_T_S2_angle_bin = cut_interval(speidel2$min_S1_T_S2_angle, length=10)

speidel2$min_S1_T_S2_angle_bin <- gsub("[0,10]", "10", speidel2$min_S1_T_S2_angle_bin, fixed=T)
speidel2$min_S1_T_S2_angle_bin <- gsub("(10,20]", "20", speidel2$min_S1_T_S2_angle_bin, fixed=T)
speidel2$min_S1_T_S2_angle_bin <- gsub("(20,30]", "30", speidel2$min_S1_T_S2_angle_bin, fixed=T)
speidel2$min_S1_T_S2_angle_bin <- gsub("(30,40]", "40", speidel2$min_S1_T_S2_angle_bin, fixed=T)
speidel2$min_S1_T_S2_angle_bin <- gsub("(40,50]", "50", speidel2$min_S1_T_S2_angle_bin, fixed=T)
speidel2$min_S1_T_S2_angle_bin <- gsub("(50,60]", "60", speidel2$min_S1_T_S2_angle_bin, fixed=T)
speidel2$min_S1_T_S2_angle_bin <- gsub("(60,70]", "70", speidel2$min_S1_T_S2_angle_bin, fixed=T)
speidel2$min_S1_T_S2_angle_bin <- gsub("(70,80]", "80", speidel2$min_S1_T_S2_angle_bin, fixed=T)
speidel2$min_S1_T_S2_angle_bin <- gsub("(80,90]", "90", speidel2$min_S1_T_S2_angle_bin, fixed=T)
speidel2$min_S1_T_S2_angle_bin <- gsub("(90,100]", "100", speidel2$min_S1_T_S2_angle_bin, fixed=T)
speidel2$min_S1_T_S2_angle_bin <- gsub("(100,110]", "110", speidel2$min_S1_T_S2_angle_bin, fixed=T)
speidel2$min_S1_T_S2_angle_bin <- gsub("(110,120]", "120", speidel2$min_S1_T_S2_angle_bin, fixed=T)
speidel2$min_S1_T_S2_angle_bin <- gsub("(120,130]", "130", speidel2$min_S1_T_S2_angle_bin, fixed=T)
speidel2$min_S1_T_S2_angle_bin <- gsub("(130,140]", "140", speidel2$min_S1_T_S2_angle_bin, fixed=T)
speidel2$min_S1_T_S2_angle_bin <- gsub("(140,150]", "150", speidel2$min_S1_T_S2_angle_bin, fixed=T)
speidel2$min_S1_T_S2_angle_bin <- gsub("(150,160]", "160", speidel2$min_S1_T_S2_angle_bin, fixed=T)
speidel2$min_S1_T_S2_angle_bin <- gsub("(160,170]", "170", speidel2$min_S1_T_S2_angle_bin, fixed=T)
speidel2$min_S1_T_S2_angle_bin <- gsub("(170,180]", "180", speidel2$min_S1_T_S2_angle_bin, fixed=T)

speidel2 <- speidel2 %>% relocate(min_S1_T_S2_angle, .after=avg.distance)
speidel2 <- speidel2 %>% relocate(min_S1_T_S2_angle_bin, .after=min_S1_T_S2_angle)
columns <-c("min_S1_T_S2_angle_bin")
speidel2[, columns] <- lapply(columns, function(x) as.numeric(speidel2[[x]]))


###binning average ST distances###
speidel2$avg.distance <- round_any(speidel2$avg.distance, 100, f = round) #!!!



###Visualizing fitting qpAdm models and pre-study model distributions in the space of "model optimality metrics"
###(average ST distances and min. STS angles). The results are shown for six selected composite model feasibility criteria, 
###with all qpAdm protocols and experiments from Speidel et al. (2024) combined:

###counting fitting and non-fitting models in 2D bins###
fit_criteria <- read_delim("fit_criteria2.txt", delim="\t", trim_ws=T)

crs = c(6,7,8,1,2,3)

###model feasibility criterion no. 1###
cr <- crs[1]
dfa <- speidel2[,c(2,29,31,33+cr)]
colnames(dfa) <- c("n_sources", "avg.distance", "min_S1_T_S2_angle_bin", "positive")
df01 <- as.data.frame(xtabs(~ n_sources+avg.distance+min_S1_T_S2_angle_bin+positive, data = dfa))
colnames(df01) <- c("n_sources", "avg.distance", "min_S1_T_S2_angle_bin", "positive", "model_count")
df01 <- df01[with(df01, order(n_sources, avg.distance, min_S1_T_S2_angle_bin, positive)), ]
df01f <- filter(df01, positive=="fitting")
df01nf <- filter(df01, positive=="non-fitting")
df01f$fit_criteria1 <- as.character(fit_criteria[cr,2])
df01f$fit_criteria2 <- as.character(fit_criteria[cr,3])
df01f$fit_criteria3 <- as.character(fit_criteria[cr,4])
df01nf$fit_criteria1 <- as.character(fit_criteria[cr,2])
df01nf$fit_criteria2 <- as.character(fit_criteria[cr,3])
df01nf$fit_criteria3 <- as.character(fit_criteria[cr,4])
df01fin <- cbind(df01f[,1:3], df01f$model_count/(df01f$model_count+df01nf$model_count))
colnames(df01fin) <- c("n_sources", "avg.distance", "min_S1_T_S2_angle_bin", "PR")
df01fin[df01fin == "NaN"] <- NA
columns <-c("avg.distance", "min_S1_T_S2_angle_bin", "PR")
df01fin[, columns] <- lapply(columns, function(x) as.numeric(as.character(df01fin[[x]])))
df01fin <- df01fin[with(df01fin, order(n_sources, avg.distance, min_S1_T_S2_angle_bin)), ]
df01fin$fit_criteria1 <- as.character(fit_criteria[cr,2])
df01fin$fit_criteria2 <- as.character(fit_criteria[cr,3])
df01fin$fit_criteria3 <- as.character(fit_criteria[cr,4])
counts01.f <- as.data.frame(xtabs(~ n_sources, data = filter(dfa, positive=="fitting")))
counts01.nf <- as.data.frame(xtabs(~ n_sources, data = filter(dfa, positive=="non-fitting")))
colnames(counts01.f) <- c("n_sources", "model_count")
colnames(counts01.nf) <- c("n_sources", "model_count")
counts01.f$fit_criteria1 <- as.character(fit_criteria[cr,2])
counts01.f$fit_criteria2 <- as.character(fit_criteria[cr,3])
counts01.f$fit_criteria3 <- as.character(fit_criteria[cr,4])
counts01.nf$fit_criteria1 <- as.character(fit_criteria[cr,2])
counts01.nf$fit_criteria2 <- as.character(fit_criteria[cr,3])
counts01.nf$fit_criteria3 <- as.character(fit_criteria[cr,4])
rm(dfa)
rm(df01)

###model feasibility criterion no. 2###
cr <- crs[2]
dfa <- speidel2[,c(2,29,31,33+cr)]
colnames(dfa) <- c("n_sources", "avg.distance", "min_S1_T_S2_angle_bin", "positive")
df02 <- as.data.frame(xtabs(~ n_sources+avg.distance+min_S1_T_S2_angle_bin+positive, data = dfa))
colnames(df02) <- c("n_sources", "avg.distance", "min_S1_T_S2_angle_bin", "positive", "model_count")
df02 <- df02[with(df02, order(n_sources, avg.distance, min_S1_T_S2_angle_bin, positive)), ]
df02f <- filter(df02, positive=="fitting")
df02nf <- filter(df02, positive=="non-fitting")
df02f$fit_criteria1 <- as.character(fit_criteria[cr,2])
df02f$fit_criteria2 <- as.character(fit_criteria[cr,3])
df02f$fit_criteria3 <- as.character(fit_criteria[cr,4])
df02nf$fit_criteria1 <- as.character(fit_criteria[cr,2])
df02nf$fit_criteria2 <- as.character(fit_criteria[cr,3])
df02nf$fit_criteria3 <- as.character(fit_criteria[cr,4])
df02fin <- cbind(df02f[,1:3], df02f$model_count/(df02f$model_count+df02nf$model_count))
colnames(df02fin) <- c("n_sources", "avg.distance", "min_S1_T_S2_angle_bin", "PR")
df02fin[df02fin == "NaN"] <- NA
columns <-c("avg.distance", "min_S1_T_S2_angle_bin", "PR")
df02fin[, columns] <- lapply(columns, function(x) as.numeric(as.character(df02fin[[x]])))
df02fin <- df02fin[with(df02fin, order(n_sources, avg.distance, min_S1_T_S2_angle_bin)), ]
df02fin$fit_criteria1 <- as.character(fit_criteria[cr,2])
df02fin$fit_criteria2 <- as.character(fit_criteria[cr,3])
df02fin$fit_criteria3 <- as.character(fit_criteria[cr,4])
counts02.f <- as.data.frame(xtabs(~ n_sources, data = filter(dfa, positive=="fitting")))
counts02.nf <- as.data.frame(xtabs(~ n_sources, data = filter(dfa, positive=="non-fitting")))
colnames(counts02.f) <- c("n_sources", "model_count")
colnames(counts02.nf) <- c("n_sources", "model_count")
counts02.f$fit_criteria1 <- as.character(fit_criteria[cr,2])
counts02.f$fit_criteria2 <- as.character(fit_criteria[cr,3])
counts02.f$fit_criteria3 <- as.character(fit_criteria[cr,4])
counts02.nf$fit_criteria1 <- as.character(fit_criteria[cr,2])
counts02.nf$fit_criteria2 <- as.character(fit_criteria[cr,3])
counts02.nf$fit_criteria3 <- as.character(fit_criteria[cr,4])
rm(dfa)
rm(df02)

###model feasibility criterion no. 3###
cr <- crs[3]
dfa <- speidel2[,c(2,29,31,33+cr)]
colnames(dfa) <- c("n_sources", "avg.distance", "min_S1_T_S2_angle_bin", "positive")
df03 <- as.data.frame(xtabs(~ n_sources+avg.distance+min_S1_T_S2_angle_bin+positive, data = dfa))
colnames(df03) <- c("n_sources", "avg.distance", "min_S1_T_S2_angle_bin", "positive", "model_count")
df03 <- df03[with(df03, order(n_sources, avg.distance, min_S1_T_S2_angle_bin, positive)), ]
df03f <- filter(df03, positive=="fitting")
df03nf <- filter(df03, positive=="non-fitting")
df03f$fit_criteria1 <- as.character(fit_criteria[cr,2])
df03f$fit_criteria2 <- as.character(fit_criteria[cr,3])
df03f$fit_criteria3 <- as.character(fit_criteria[cr,4])
df03nf$fit_criteria1 <- as.character(fit_criteria[cr,2])
df03nf$fit_criteria2 <- as.character(fit_criteria[cr,3])
df03nf$fit_criteria3 <- as.character(fit_criteria[cr,4])
df03fin <- cbind(df03f[,1:3], df03f$model_count/(df03f$model_count+df03nf$model_count))
colnames(df03fin) <- c("n_sources", "avg.distance", "min_S1_T_S2_angle_bin", "PR")
df03fin[df03fin == "NaN"] <- NA
columns <-c("avg.distance", "min_S1_T_S2_angle_bin", "PR")
df03fin[, columns] <- lapply(columns, function(x) as.numeric(as.character(df03fin[[x]])))
df03fin <- df03fin[with(df03fin, order(n_sources, avg.distance, min_S1_T_S2_angle_bin)), ]
df03fin$fit_criteria1 <- as.character(fit_criteria[cr,2])
df03fin$fit_criteria2 <- as.character(fit_criteria[cr,3])
df03fin$fit_criteria3 <- as.character(fit_criteria[cr,4])
counts03.f <- as.data.frame(xtabs(~ n_sources, data = filter(dfa, positive=="fitting")))
counts03.nf <- as.data.frame(xtabs(~ n_sources, data = filter(dfa, positive=="non-fitting")))
colnames(counts03.f) <- c("n_sources", "model_count")
colnames(counts03.nf) <- c("n_sources", "model_count")
counts03.f$fit_criteria1 <- as.character(fit_criteria[cr,2])
counts03.f$fit_criteria2 <- as.character(fit_criteria[cr,3])
counts03.f$fit_criteria3 <- as.character(fit_criteria[cr,4])
counts03.nf$fit_criteria1 <- as.character(fit_criteria[cr,2])
counts03.nf$fit_criteria2 <- as.character(fit_criteria[cr,3])
counts03.nf$fit_criteria3 <- as.character(fit_criteria[cr,4])
rm(dfa)
rm(df03)

###model feasibility criterion no. 4###
cr <- crs[4]
dfa <- speidel2[,c(2,29,31,33+cr)]
colnames(dfa) <- c("n_sources", "avg.distance", "min_S1_T_S2_angle_bin", "positive")
df04 <- as.data.frame(xtabs(~ n_sources+avg.distance+min_S1_T_S2_angle_bin+positive, data = dfa))
colnames(df04) <- c("n_sources", "avg.distance", "min_S1_T_S2_angle_bin", "positive", "model_count")
df04 <- df04[with(df04, order(n_sources, avg.distance, min_S1_T_S2_angle_bin, positive)), ]
df04f <- filter(df04, positive=="fitting")
df04nf <- filter(df04, positive=="non-fitting")
df04f$fit_criteria1 <- as.character(fit_criteria[cr,2])
df04f$fit_criteria2 <- as.character(fit_criteria[cr,3])
df04f$fit_criteria3 <- as.character(fit_criteria[cr,4])
df04nf$fit_criteria1 <- as.character(fit_criteria[cr,2])
df04nf$fit_criteria2 <- as.character(fit_criteria[cr,3])
df04nf$fit_criteria3 <- as.character(fit_criteria[cr,4])
df04fin <- cbind(df04f[,1:3], df04f$model_count/(df04f$model_count+df04nf$model_count))
colnames(df04fin) <- c("n_sources", "avg.distance", "min_S1_T_S2_angle_bin", "PR")
df04fin[df04fin == "NaN"] <- NA
columns <-c("avg.distance", "min_S1_T_S2_angle_bin", "PR")
df04fin[, columns] <- lapply(columns, function(x) as.numeric(as.character(df04fin[[x]])))
df04fin <- df04fin[with(df04fin, order(n_sources, avg.distance, min_S1_T_S2_angle_bin)), ]
df04fin$fit_criteria1 <- as.character(fit_criteria[cr,2])
df04fin$fit_criteria2 <- as.character(fit_criteria[cr,3])
df04fin$fit_criteria3 <- as.character(fit_criteria[cr,4])
counts04.f <- as.data.frame(xtabs(~ n_sources, data = filter(dfa, positive=="fitting")))
counts04.nf <- as.data.frame(xtabs(~ n_sources, data = filter(dfa, positive=="non-fitting")))
colnames(counts04.f) <- c("n_sources", "model_count")
colnames(counts04.nf) <- c("n_sources", "model_count")
counts04.f$fit_criteria1 <- as.character(fit_criteria[cr,2])
counts04.f$fit_criteria2 <- as.character(fit_criteria[cr,3])
counts04.f$fit_criteria3 <- as.character(fit_criteria[cr,4])
counts04.nf$fit_criteria1 <- as.character(fit_criteria[cr,2])
counts04.nf$fit_criteria2 <- as.character(fit_criteria[cr,3])
counts04.nf$fit_criteria3 <- as.character(fit_criteria[cr,4])
rm(dfa)
rm(df04)

###model feasibility criterion no. 5###
cr <- crs[5]
dfa <- speidel2[,c(2,29,31,33+cr)]
colnames(dfa) <- c("n_sources", "avg.distance", "min_S1_T_S2_angle_bin", "positive")
df05 <- as.data.frame(xtabs(~ n_sources+avg.distance+min_S1_T_S2_angle_bin+positive, data = dfa))
colnames(df05) <- c("n_sources", "avg.distance", "min_S1_T_S2_angle_bin", "positive", "model_count")
df05 <- df05[with(df05, order(n_sources, avg.distance, min_S1_T_S2_angle_bin, positive)), ]
df05f <- filter(df05, positive=="fitting")
df05nf <- filter(df05, positive=="non-fitting")
df05f$fit_criteria1 <- as.character(fit_criteria[cr,2])
df05f$fit_criteria2 <- as.character(fit_criteria[cr,3])
df05f$fit_criteria3 <- as.character(fit_criteria[cr,4])
df05nf$fit_criteria1 <- as.character(fit_criteria[cr,2])
df05nf$fit_criteria2 <- as.character(fit_criteria[cr,3])
df05nf$fit_criteria3 <- as.character(fit_criteria[cr,4])
df05fin <- cbind(df05f[,1:3], df05f$model_count/(df05f$model_count+df05nf$model_count))
colnames(df05fin) <- c("n_sources", "avg.distance", "min_S1_T_S2_angle_bin", "PR")
df05fin[df05fin == "NaN"] <- NA
columns <-c("avg.distance", "min_S1_T_S2_angle_bin", "PR")
df05fin[, columns] <- lapply(columns, function(x) as.numeric(as.character(df05fin[[x]])))
df05fin <- df05fin[with(df05fin, order(n_sources, avg.distance, min_S1_T_S2_angle_bin)), ]
df05fin$fit_criteria1 <- as.character(fit_criteria[cr,2])
df05fin$fit_criteria2 <- as.character(fit_criteria[cr,3])
df05fin$fit_criteria3 <- as.character(fit_criteria[cr,4])
counts05.f <- as.data.frame(xtabs(~ n_sources, data = filter(dfa, positive=="fitting")))
counts05.nf <- as.data.frame(xtabs(~ n_sources, data = filter(dfa, positive=="non-fitting")))
colnames(counts05.f) <- c("n_sources", "model_count")
colnames(counts05.nf) <- c("n_sources", "model_count")
counts05.f$fit_criteria1 <- as.character(fit_criteria[cr,2])
counts05.f$fit_criteria2 <- as.character(fit_criteria[cr,3])
counts05.f$fit_criteria3 <- as.character(fit_criteria[cr,4])
counts05.nf$fit_criteria1 <- as.character(fit_criteria[cr,2])
counts05.nf$fit_criteria2 <- as.character(fit_criteria[cr,3])
counts05.nf$fit_criteria3 <- as.character(fit_criteria[cr,4])
rm(dfa)
rm(df05)

###model feasibility criterion no. 6###
cr <- crs[6]
dfa <- speidel2[,c(2,29,31,33+cr)]
colnames(dfa) <- c("n_sources", "avg.distance", "min_S1_T_S2_angle_bin", "positive")
df06 <- as.data.frame(xtabs(~ n_sources+avg.distance+min_S1_T_S2_angle_bin+positive, data = dfa))
colnames(df06) <- c("n_sources", "avg.distance", "min_S1_T_S2_angle_bin", "positive", "model_count")
df06 <- df06[with(df06, order(n_sources, avg.distance, min_S1_T_S2_angle_bin, positive)), ]
df06f <- filter(df06, positive=="fitting")
df06nf <- filter(df06, positive=="non-fitting")
df06f$fit_criteria1 <- as.character(fit_criteria[cr,2])
df06f$fit_criteria2 <- as.character(fit_criteria[cr,3])
df06f$fit_criteria3 <- as.character(fit_criteria[cr,4])
df06nf$fit_criteria1 <- as.character(fit_criteria[cr,2])
df06nf$fit_criteria2 <- as.character(fit_criteria[cr,3])
df06nf$fit_criteria3 <- as.character(fit_criteria[cr,4])
df06fin <- cbind(df06f[,1:3], df06f$model_count/(df06f$model_count+df06nf$model_count))
colnames(df06fin) <- c("n_sources", "avg.distance", "min_S1_T_S2_angle_bin", "PR")
df06fin[df06fin == "NaN"] <- NA
columns <-c("avg.distance", "min_S1_T_S2_angle_bin", "PR")
df06fin[, columns] <- lapply(columns, function(x) as.numeric(as.character(df06fin[[x]])))
df06fin <- df06fin[with(df06fin, order(n_sources, avg.distance, min_S1_T_S2_angle_bin)), ]
df06fin$fit_criteria1 <- as.character(fit_criteria[cr,2])
df06fin$fit_criteria2 <- as.character(fit_criteria[cr,3])
df06fin$fit_criteria3 <- as.character(fit_criteria[cr,4])
counts06.f <- as.data.frame(xtabs(~ n_sources, data = filter(dfa, positive=="fitting")))
counts06.nf <- as.data.frame(xtabs(~ n_sources, data = filter(dfa, positive=="non-fitting")))
colnames(counts06.f) <- c("n_sources", "model_count")
colnames(counts06.nf) <- c("n_sources", "model_count")
counts06.f$fit_criteria1 <- as.character(fit_criteria[cr,2])
counts06.f$fit_criteria2 <- as.character(fit_criteria[cr,3])
counts06.f$fit_criteria3 <- as.character(fit_criteria[cr,4])
counts06.nf$fit_criteria1 <- as.character(fit_criteria[cr,2])
counts06.nf$fit_criteria2 <- as.character(fit_criteria[cr,3])
counts06.nf$fit_criteria3 <- as.character(fit_criteria[cr,4])
rm(dfa)
rm(df06)


dff <- rbind(df01f, df02f, df03f, df04f, df05f, df06f)
dfnf <- rbind(df01nf, df02nf, df03nf, df04nf, df05nf, df06nf)
dffin <- rbind(df01fin, df02fin, df03fin, df04fin, df05fin, df06fin)
counts.f <- rbind(counts01.f, counts02.f, counts03.f, counts04.f, counts05.f, counts06.f)
counts.nf <- rbind(counts01.nf, counts02.nf, counts03.nf, counts04.nf, counts05.nf, counts06.nf)
rm(df01f)
rm(df02f)
rm(df03f)
rm(df04f)
rm(df05f)
rm(df06f)
rm(df01nf)
rm(df02nf)
rm(df03nf)
rm(df04nf)
rm(df05nf)
rm(df06nf)
rm(df01fin)
rm(df02fin)
rm(df03fin)
rm(df04fin)
rm(df05fin)
rm(df06fin)
rm(counts01.f)
rm(counts02.f)
rm(counts03.f)
rm(counts04.f)
rm(counts05.f)
rm(counts06.f)
rm(counts01.nf)
rm(counts02.nf)
rm(counts03.nf)
rm(counts04.nf)
rm(counts05.nf)
rm(counts06.nf)

dffin$nf <- dff$model_count
dffin$nt <- dff$model_count+dfnf$model_count
dffin$nf <- as.double(dffin$nf)
dffin$nt <- as.double(dffin$nt)
dffin$nt[dffin$nt==0] <- NA ##IMPORTANT!
corr <- dffin  %>% group_by(n_sources, fit_criteria1, fit_criteria2, fit_criteria3)  %>% do(tidy(cor.test(.$nt, .$nf, alternative="two.sided", method = "spearman"))) 
#counts.f$PR_global <- counts.f$model_count/(counts.f$model_count+counts.nf$model_count)
#counts.f$n_total <- counts.f$model_count+counts.nf$model_count
#counts.f$n_total[13:84] <- NA


###plotting the results on real data:
p00<-ggplot(speidel2, aes(x=avg.distance, y=min_S1_T_S2_angle_bin))+
  geom_bin2d(aes(fill = after_stat(count)), binwidth = c(100, 10), origin = c(-50, -5))+
  scale_fill_gradientn(colours = pals::parula(20), trans = "log", breaks = c(1,3,10,30,100), labels = c(1,3,10,30,100), na.value='grey95')+
  ggtitle("Speidel et al. 2024: models to be tested (5,171 in total)")+
  scale_x_continuous(breaks=seq(200,2400,400), expand=c(0.05,0.05))+
  scale_y_continuous(breaks=c(30,60,90,120,150,180), expand=c(0.05,0.05))+
  coord_cartesian(xlim=c(200,2400), ylim=c(10,180), expand=T, default=F, clip="on")+
  #xlab("")+
  ylab("")+
  guides(fill = guide_colorbar(title = '', order=1, barwidth=4, barheight=18))+
  theme(plot.title = element_text(size=32, face='bold'),
        axis.text.x=element_blank(),
        axis.text.y=element_text(size=24, face='bold'),
        axis.title.x=element_blank(),
        axis.title.y=element_text(size=36, face='bold'),
        legend.position = "right",
        legend.key=element_rect(fill = 'white'),
        legend.title=element_text(size=36, face = "bold"),
        legend.text=element_text(size=30),
        panel.background=element_rect(fill = 'grey95'),
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        plot.margin = unit(c(10,110,10,10), "pt"),
        strip.text=element_text(size=26, face='bold'))+
  facet_nested(. ~ n_sources)

p01<-ggplot(dffin, aes(x=avg.distance, y=min_S1_T_S2_angle_bin))+
  geom_vline(aes(xintercept=1200), colour="red3", linewidth=1, alpha=0.5)+
  geom_tile(aes(fill = nf))+
  geom_text(data=counts.f, aes(x=400, y=200, label=scales::comma(model_count)), colour="darkgreen", size=10, fontface="bold", check_overlap=T)+
  geom_text(data=counts.nf, aes(x=2000, y=200, label=scales::comma(model_count)), colour="red3", size=10, fontface="bold", check_overlap=T)+
  geom_text(data=corr, aes(x=2400, y=150, label=round(estimate, 2)), colour="tomato4", size=9, fontface="bold", check_overlap=T, angle=90)+
  geom_hline(aes(yintercept=125), colour="red3", linewidth=1, alpha=0.5)+
  geom_hline(aes(yintercept=95), colour="red3", linewidth=1, alpha=0.5)+
  scale_fill_gradientn(colours = pals::parula(20), trans = "log", breaks = c(1,3,10,30), labels = c(1,3,10,30), na.value='grey95')+
  scale_x_continuous(breaks=seq(200,2400,400), expand=c(0.05,0.05))+
  scale_y_continuous(breaks=c(30,60,90,120,150,180), expand=c(0.05,0.05))+
  coord_cartesian(xlim=c(200,2400), ylim=c(10,200), expand=T, default=F, clip="on")+
  ggtitle('counts of feasible 2-3-way models\nin the space of optimality metrics, log. color scale,\nprogressively more stringent feasibility criteria')+
  xlab("average source-target distance, km")+
  ylab("min. S1-target-S2 angle")+
  guides(fill = guide_colorbar(title = '', order=1, barwidth=4, barheight=24))+
  theme(plot.title = element_text(color="darkmagenta", face="bold", size=28),
        axis.text.x=element_text(size=24, face='bold', angle=45, hjust=0.95, vjust=0.95),
        axis.text.y=element_text(size=24, face='bold'),
        axis.title.x=element_text(size=36, face='bold'),
        axis.title.y=element_text(size=36, face='bold'),
        legend.position = "right",
        legend.key=element_rect(fill = 'white'),
        legend.title=element_text(size=36, face = "bold"),
        legend.text=element_text(size=30),
        panel.background=element_rect(fill = 'white'),
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        plot.margin = unit(c(5,10,5,10), "pt"),
        strip.text=element_text(size=26, face='bold'))+
  facet_nested(fit_criteria3 + factor(fit_criteria1, level=c("% (0,1)", "%+/-2SE (0,1)")) + fit_criteria2 ~ n_sources)

pall <- plot_grid(
  plot_grid(p00+theme(legend.position="none"), p01+theme(legend.position="none"), ncol = 1, axis="l", align = "h", rel_heights = c(1.2,6.8)),
  plot_grid(get_legend(p00), get_legend(p01), ncol = 1, axis="l", align = "h", rel_heights = c(1.2,6.8)),
  rel_widths = c(24,6))



###adding simulated data, distal and proximal rotating 2-3-way models###
###these are intermediate outputs of the script "fitting_qpAdm_models_in_the_space_of_model_optimality_metrics_Figure_6.R":
dfrr <- read_bulk(directory = "./Isolation_by_distance/stepping-stone_simulations/qpAdm/histograms_model_metrics_basic_1-4way/", name_filter=".*_rotating_.*_high_for_2Doptimality_metrics.*", fun = readr::read_delim, delim="\t")[,c(1:4,6,10:48)]
dfrr <- filter(dfrr, model!="4-way")

###binning min. STS angles###
dfrr<-filter(dfrr, !is.na(min_S1_T_S2_angle))
dfrr$min_S1_T_S2_angle_bin = cut_interval(dfrr$min_S1_T_S2_angle, length=10)

dfrr$min_S1_T_S2_angle_bin <- gsub("[0,10]", "10", dfrr$min_S1_T_S2_angle_bin, fixed=T)
dfrr$min_S1_T_S2_angle_bin <- gsub("(10,20]", "20", dfrr$min_S1_T_S2_angle_bin, fixed=T)
dfrr$min_S1_T_S2_angle_bin <- gsub("(20,30]", "30", dfrr$min_S1_T_S2_angle_bin, fixed=T)
dfrr$min_S1_T_S2_angle_bin <- gsub("(30,40]", "40", dfrr$min_S1_T_S2_angle_bin, fixed=T)
dfrr$min_S1_T_S2_angle_bin <- gsub("(40,50]", "50", dfrr$min_S1_T_S2_angle_bin, fixed=T)
dfrr$min_S1_T_S2_angle_bin <- gsub("(50,60]", "60", dfrr$min_S1_T_S2_angle_bin, fixed=T)
dfrr$min_S1_T_S2_angle_bin <- gsub("(60,70]", "70", dfrr$min_S1_T_S2_angle_bin, fixed=T)
dfrr$min_S1_T_S2_angle_bin <- gsub("(70,80]", "80", dfrr$min_S1_T_S2_angle_bin, fixed=T)
dfrr$min_S1_T_S2_angle_bin <- gsub("(80,90]", "90", dfrr$min_S1_T_S2_angle_bin, fixed=T)
dfrr$min_S1_T_S2_angle_bin <- gsub("(90,100]", "100", dfrr$min_S1_T_S2_angle_bin, fixed=T)
dfrr$min_S1_T_S2_angle_bin <- gsub("(100,110]", "110", dfrr$min_S1_T_S2_angle_bin, fixed=T)
dfrr$min_S1_T_S2_angle_bin <- gsub("(110,120]", "120", dfrr$min_S1_T_S2_angle_bin, fixed=T)
dfrr$min_S1_T_S2_angle_bin <- gsub("(120,130]", "130", dfrr$min_S1_T_S2_angle_bin, fixed=T)
dfrr$min_S1_T_S2_angle_bin <- gsub("(130,140]", "140", dfrr$min_S1_T_S2_angle_bin, fixed=T)
dfrr$min_S1_T_S2_angle_bin <- gsub("(140,150]", "150", dfrr$min_S1_T_S2_angle_bin, fixed=T)
dfrr$min_S1_T_S2_angle_bin <- gsub("(150,160]", "160", dfrr$min_S1_T_S2_angle_bin, fixed=T)
dfrr$min_S1_T_S2_angle_bin <- gsub("(160,170]", "170", dfrr$min_S1_T_S2_angle_bin, fixed=T)
dfrr$min_S1_T_S2_angle_bin <- gsub("(170,180]", "180", dfrr$min_S1_T_S2_angle_bin, fixed=T)

dfrr <- dfrr %>% relocate(min_S1_T_S2_angle, .after=avg..distance)
dfrr <- dfrr %>% relocate(min_S1_T_S2_angle_bin, .after=min_S1_T_S2_angle)
columns <-c("min_S1_T_S2_angle_bin")
dfrr[, columns] <- lapply(columns, function(x) as.numeric(dfrr[[x]]))


###binning average ST distances###
dfrr$avg..distance <- round_any(dfrr$avg..distance, 0.5, f = round)


###Visualizing fitting qpAdm models and pre-study model distributions in the space of "model optimality metrics"
###(average ST distances and min. STS angles). The results are shown for six selected composite model feasibility criteria:

###counting fitting and non-fitting models in 2D bins###
fit_criteria <- read_delim("fit_criteria2.txt", delim="\t", trim_ws=T)

crs = c(6,7,8,1,2,3)

###model feasibility criterion no. 1###
cr <- crs[1]
dfa <- dfrr[,c(4:5,7,9,9+21+cr)]
colnames(dfa) <- c("simulation.replicate", "model", "avg..distance", "min_S1_T_S2_angle_bin", "positive")
df01 <- as.data.frame(xtabs(~ model+avg..distance+min_S1_T_S2_angle_bin+positive, data = dfa))
colnames(df01) <- c("model", "avg..distance", "min_S1_T_S2_angle_bin", "positive", "model_count")
df01 <- df01[with(df01, order(model, avg..distance, min_S1_T_S2_angle_bin, positive)), ]
df01f <- filter(df01, positive=="fitting")
df01nf <- filter(df01, positive=="non-fitting")
df01f$fit_criteria1 <- as.character(fit_criteria[cr,2])
df01f$fit_criteria2 <- as.character(fit_criteria[cr,3])
df01f$fit_criteria3 <- as.character(fit_criteria[cr,4])
df01nf$fit_criteria1 <- as.character(fit_criteria[cr,2])
df01nf$fit_criteria2 <- as.character(fit_criteria[cr,3])
df01nf$fit_criteria3 <- as.character(fit_criteria[cr,4])
df01fin <- cbind(df01f[,1:3], df01f$model_count/(df01f$model_count+df01nf$model_count))
colnames(df01fin) <- c("model", "avg..distance", "min_S1_T_S2_angle_bin", "PR")
df01fin[df01fin == "NaN"] <- NA
columns <-c("avg..distance", "min_S1_T_S2_angle_bin", "PR")
df01fin[, columns] <- lapply(columns, function(x) as.numeric(as.character(df01fin[[x]])))
df01fin <- df01fin[with(df01fin, order(model, avg..distance, min_S1_T_S2_angle_bin)), ]
df01fin$fit_criteria1 <- as.character(fit_criteria[cr,2])
df01fin$fit_criteria2 <- as.character(fit_criteria[cr,3])
df01fin$fit_criteria3 <- as.character(fit_criteria[cr,4])
counts01.f <- as.data.frame(xtabs(~ model, data = filter(dfa, positive=="fitting")))
counts01.nf <- as.data.frame(xtabs(~ model, data = filter(dfa, positive=="non-fitting")))
colnames(counts01.f) <- c("model", "model_count")
colnames(counts01.nf) <- c("model", "model_count")
counts01.f$fit_criteria1 <- as.character(fit_criteria[cr,2])
counts01.f$fit_criteria2 <- as.character(fit_criteria[cr,3])
counts01.f$fit_criteria3 <- as.character(fit_criteria[cr,4])
counts01.nf$fit_criteria1 <- as.character(fit_criteria[cr,2])
counts01.nf$fit_criteria2 <- as.character(fit_criteria[cr,3])
counts01.nf$fit_criteria3 <- as.character(fit_criteria[cr,4])
rm(dfa)
rm(df01)

###model feasibility criterion no. 2###
cr <- crs[2]
dfa <- dfrr[,c(4:5,7,9,9+21+cr)]
colnames(dfa) <- c("simulation.replicate", "model", "avg..distance", "min_S1_T_S2_angle_bin", "positive")
df02 <- as.data.frame(xtabs(~ model+avg..distance+min_S1_T_S2_angle_bin+positive, data = dfa))
colnames(df02) <- c("model", "avg..distance", "min_S1_T_S2_angle_bin", "positive", "model_count")
df02 <- df02[with(df02, order(model, avg..distance, min_S1_T_S2_angle_bin, positive)), ]
df02f <- filter(df02, positive=="fitting")
df02nf <- filter(df02, positive=="non-fitting")
df02f$fit_criteria1 <- as.character(fit_criteria[cr,2])
df02f$fit_criteria2 <- as.character(fit_criteria[cr,3])
df02f$fit_criteria3 <- as.character(fit_criteria[cr,4])
df02nf$fit_criteria1 <- as.character(fit_criteria[cr,2])
df02nf$fit_criteria2 <- as.character(fit_criteria[cr,3])
df02nf$fit_criteria3 <- as.character(fit_criteria[cr,4])
df02fin <- cbind(df02f[,1:3], df02f$model_count/(df02f$model_count+df02nf$model_count))
colnames(df02fin) <- c("model", "avg..distance", "min_S1_T_S2_angle_bin", "PR")
df02fin[df02fin == "NaN"] <- NA
columns <-c("avg..distance", "min_S1_T_S2_angle_bin", "PR")
df02fin[, columns] <- lapply(columns, function(x) as.numeric(as.character(df02fin[[x]])))
df02fin <- df02fin[with(df02fin, order(model, avg..distance, min_S1_T_S2_angle_bin)), ]
df02fin$fit_criteria1 <- as.character(fit_criteria[cr,2])
df02fin$fit_criteria2 <- as.character(fit_criteria[cr,3])
df02fin$fit_criteria3 <- as.character(fit_criteria[cr,4])
counts02.f <- as.data.frame(xtabs(~ model, data = filter(dfa, positive=="fitting")))
counts02.nf <- as.data.frame(xtabs(~ model, data = filter(dfa, positive=="non-fitting")))
colnames(counts02.f) <- c("model", "model_count")
colnames(counts02.nf) <- c("model", "model_count")
counts02.f$fit_criteria1 <- as.character(fit_criteria[cr,2])
counts02.f$fit_criteria2 <- as.character(fit_criteria[cr,3])
counts02.f$fit_criteria3 <- as.character(fit_criteria[cr,4])
counts02.nf$fit_criteria1 <- as.character(fit_criteria[cr,2])
counts02.nf$fit_criteria2 <- as.character(fit_criteria[cr,3])
counts02.nf$fit_criteria3 <- as.character(fit_criteria[cr,4])
rm(dfa)
rm(df02)

###model feasibility criterion no. 3###
cr <- crs[3]
dfa <- dfrr[,c(4:5,7,9,9+21+cr)]
colnames(dfa) <- c("simulation.replicate", "model", "avg..distance", "min_S1_T_S2_angle_bin", "positive")
df03 <- as.data.frame(xtabs(~ model+avg..distance+min_S1_T_S2_angle_bin+positive, data = dfa))
colnames(df03) <- c("model", "avg..distance", "min_S1_T_S2_angle_bin", "positive", "model_count")
df03 <- df03[with(df03, order(model, avg..distance, min_S1_T_S2_angle_bin, positive)), ]
df03f <- filter(df03, positive=="fitting")
df03nf <- filter(df03, positive=="non-fitting")
df03f$fit_criteria1 <- as.character(fit_criteria[cr,2])
df03f$fit_criteria2 <- as.character(fit_criteria[cr,3])
df03f$fit_criteria3 <- as.character(fit_criteria[cr,4])
df03nf$fit_criteria1 <- as.character(fit_criteria[cr,2])
df03nf$fit_criteria2 <- as.character(fit_criteria[cr,3])
df03nf$fit_criteria3 <- as.character(fit_criteria[cr,4])
df03fin <- cbind(df03f[,1:3], df03f$model_count/(df03f$model_count+df03nf$model_count))
colnames(df03fin) <- c("model", "avg..distance", "min_S1_T_S2_angle_bin", "PR")
df03fin[df03fin == "NaN"] <- NA
columns <-c("avg..distance", "min_S1_T_S2_angle_bin", "PR")
df03fin[, columns] <- lapply(columns, function(x) as.numeric(as.character(df03fin[[x]])))
df03fin <- df03fin[with(df03fin, order(model, avg..distance, min_S1_T_S2_angle_bin)), ]
df03fin$fit_criteria1 <- as.character(fit_criteria[cr,2])
df03fin$fit_criteria2 <- as.character(fit_criteria[cr,3])
df03fin$fit_criteria3 <- as.character(fit_criteria[cr,4])
counts03.f <- as.data.frame(xtabs(~ model, data = filter(dfa, positive=="fitting")))
counts03.nf <- as.data.frame(xtabs(~ model, data = filter(dfa, positive=="non-fitting")))
colnames(counts03.f) <- c("model", "model_count")
colnames(counts03.nf) <- c("model", "model_count")
counts03.f$fit_criteria1 <- as.character(fit_criteria[cr,2])
counts03.f$fit_criteria2 <- as.character(fit_criteria[cr,3])
counts03.f$fit_criteria3 <- as.character(fit_criteria[cr,4])
counts03.nf$fit_criteria1 <- as.character(fit_criteria[cr,2])
counts03.nf$fit_criteria2 <- as.character(fit_criteria[cr,3])
counts03.nf$fit_criteria3 <- as.character(fit_criteria[cr,4])
rm(dfa)
rm(df03)

###model feasibility criterion no. 4###
cr <- crs[4]
dfa <- dfrr[,c(4:5,7,9,9+21+cr)]
colnames(dfa) <- c("simulation.replicate", "model", "avg..distance", "min_S1_T_S2_angle_bin", "positive")
df04 <- as.data.frame(xtabs(~ model+avg..distance+min_S1_T_S2_angle_bin+positive, data = dfa))
colnames(df04) <- c("model", "avg..distance", "min_S1_T_S2_angle_bin", "positive", "model_count")
df04 <- df04[with(df04, order(model, avg..distance, min_S1_T_S2_angle_bin, positive)), ]
df04f <- filter(df04, positive=="fitting")
df04nf <- filter(df04, positive=="non-fitting")
df04f$fit_criteria1 <- as.character(fit_criteria[cr,2])
df04f$fit_criteria2 <- as.character(fit_criteria[cr,3])
df04f$fit_criteria3 <- as.character(fit_criteria[cr,4])
df04nf$fit_criteria1 <- as.character(fit_criteria[cr,2])
df04nf$fit_criteria2 <- as.character(fit_criteria[cr,3])
df04nf$fit_criteria3 <- as.character(fit_criteria[cr,4])
df04fin <- cbind(df04f[,1:3], df04f$model_count/(df04f$model_count+df04nf$model_count))
colnames(df04fin) <- c("model", "avg..distance", "min_S1_T_S2_angle_bin", "PR")
df04fin[df04fin == "NaN"] <- NA
columns <-c("avg..distance", "min_S1_T_S2_angle_bin", "PR")
df04fin[, columns] <- lapply(columns, function(x) as.numeric(as.character(df04fin[[x]])))
df04fin <- df04fin[with(df04fin, order(model, avg..distance, min_S1_T_S2_angle_bin)), ]
df04fin$fit_criteria1 <- as.character(fit_criteria[cr,2])
df04fin$fit_criteria2 <- as.character(fit_criteria[cr,3])
df04fin$fit_criteria3 <- as.character(fit_criteria[cr,4])
counts04.f <- as.data.frame(xtabs(~ model, data = filter(dfa, positive=="fitting")))
counts04.nf <- as.data.frame(xtabs(~ model, data = filter(dfa, positive=="non-fitting")))
colnames(counts04.f) <- c("model", "model_count")
colnames(counts04.nf) <- c("model", "model_count")
counts04.f$fit_criteria1 <- as.character(fit_criteria[cr,2])
counts04.f$fit_criteria2 <- as.character(fit_criteria[cr,3])
counts04.f$fit_criteria3 <- as.character(fit_criteria[cr,4])
counts04.nf$fit_criteria1 <- as.character(fit_criteria[cr,2])
counts04.nf$fit_criteria2 <- as.character(fit_criteria[cr,3])
counts04.nf$fit_criteria3 <- as.character(fit_criteria[cr,4])
rm(dfa)
rm(df04)

###model feasibility criterion no. 5###
cr <- crs[5]
dfa <- dfrr[,c(4:5,7,9,9+21+cr)]
colnames(dfa) <- c("simulation.replicate", "model", "avg..distance", "min_S1_T_S2_angle_bin", "positive")
df05 <- as.data.frame(xtabs(~ model+avg..distance+min_S1_T_S2_angle_bin+positive, data = dfa))
colnames(df05) <- c("model", "avg..distance", "min_S1_T_S2_angle_bin", "positive", "model_count")
df05 <- df05[with(df05, order(model, avg..distance, min_S1_T_S2_angle_bin, positive)), ]
df05f <- filter(df05, positive=="fitting")
df05nf <- filter(df05, positive=="non-fitting")
df05f$fit_criteria1 <- as.character(fit_criteria[cr,2])
df05f$fit_criteria2 <- as.character(fit_criteria[cr,3])
df05f$fit_criteria3 <- as.character(fit_criteria[cr,4])
df05nf$fit_criteria1 <- as.character(fit_criteria[cr,2])
df05nf$fit_criteria2 <- as.character(fit_criteria[cr,3])
df05nf$fit_criteria3 <- as.character(fit_criteria[cr,4])
df05fin <- cbind(df05f[,1:3], df05f$model_count/(df05f$model_count+df05nf$model_count))
colnames(df05fin) <- c("model", "avg..distance", "min_S1_T_S2_angle_bin", "PR")
df05fin[df05fin == "NaN"] <- NA
columns <-c("avg..distance", "min_S1_T_S2_angle_bin", "PR")
df05fin[, columns] <- lapply(columns, function(x) as.numeric(as.character(df05fin[[x]])))
df05fin <- df05fin[with(df05fin, order(model, avg..distance, min_S1_T_S2_angle_bin)), ]
df05fin$fit_criteria1 <- as.character(fit_criteria[cr,2])
df05fin$fit_criteria2 <- as.character(fit_criteria[cr,3])
df05fin$fit_criteria3 <- as.character(fit_criteria[cr,4])
counts05.f <- as.data.frame(xtabs(~ model, data = filter(dfa, positive=="fitting")))
counts05.nf <- as.data.frame(xtabs(~ model, data = filter(dfa, positive=="non-fitting")))
colnames(counts05.f) <- c("model", "model_count")
colnames(counts05.nf) <- c("model", "model_count")
counts05.f$fit_criteria1 <- as.character(fit_criteria[cr,2])
counts05.f$fit_criteria2 <- as.character(fit_criteria[cr,3])
counts05.f$fit_criteria3 <- as.character(fit_criteria[cr,4])
counts05.nf$fit_criteria1 <- as.character(fit_criteria[cr,2])
counts05.nf$fit_criteria2 <- as.character(fit_criteria[cr,3])
counts05.nf$fit_criteria3 <- as.character(fit_criteria[cr,4])
rm(dfa)
rm(df05)

###model feasibility criterion no. 6###
cr <- crs[6]
dfa <- dfrr[,c(4:5,7,9,9+21+cr)]
colnames(dfa) <- c("simulation.replicate", "model", "avg..distance", "min_S1_T_S2_angle_bin", "positive")
df06 <- as.data.frame(xtabs(~ model+avg..distance+min_S1_T_S2_angle_bin+positive, data = dfa))
colnames(df06) <- c("model", "avg..distance", "min_S1_T_S2_angle_bin", "positive", "model_count")
df06 <- df06[with(df06, order(model, avg..distance, min_S1_T_S2_angle_bin, positive)), ]
df06f <- filter(df06, positive=="fitting")
df06nf <- filter(df06, positive=="non-fitting")
df06f$fit_criteria1 <- as.character(fit_criteria[cr,2])
df06f$fit_criteria2 <- as.character(fit_criteria[cr,3])
df06f$fit_criteria3 <- as.character(fit_criteria[cr,4])
df06nf$fit_criteria1 <- as.character(fit_criteria[cr,2])
df06nf$fit_criteria2 <- as.character(fit_criteria[cr,3])
df06nf$fit_criteria3 <- as.character(fit_criteria[cr,4])
df06fin <- cbind(df06f[,1:3], df06f$model_count/(df06f$model_count+df06nf$model_count))
colnames(df06fin) <- c("model", "avg..distance", "min_S1_T_S2_angle_bin", "PR")
df06fin[df06fin == "NaN"] <- NA
columns <-c("avg..distance", "min_S1_T_S2_angle_bin", "PR")
df06fin[, columns] <- lapply(columns, function(x) as.numeric(as.character(df06fin[[x]])))
df06fin <- df06fin[with(df06fin, order(model, avg..distance, min_S1_T_S2_angle_bin)), ]
df06fin$fit_criteria1 <- as.character(fit_criteria[cr,2])
df06fin$fit_criteria2 <- as.character(fit_criteria[cr,3])
df06fin$fit_criteria3 <- as.character(fit_criteria[cr,4])
counts06.f <- as.data.frame(xtabs(~ model, data = filter(dfa, positive=="fitting")))
counts06.nf <- as.data.frame(xtabs(~ model, data = filter(dfa, positive=="non-fitting")))
colnames(counts06.f) <- c("model", "model_count")
colnames(counts06.nf) <- c("model", "model_count")
counts06.f$fit_criteria1 <- as.character(fit_criteria[cr,2])
counts06.f$fit_criteria2 <- as.character(fit_criteria[cr,3])
counts06.f$fit_criteria3 <- as.character(fit_criteria[cr,4])
counts06.nf$fit_criteria1 <- as.character(fit_criteria[cr,2])
counts06.nf$fit_criteria2 <- as.character(fit_criteria[cr,3])
counts06.nf$fit_criteria3 <- as.character(fit_criteria[cr,4])
rm(dfa)
rm(df06)


dff <- rbind(df01f, df02f, df03f, df04f, df05f, df06f)
dfnf <- rbind(df01nf, df02nf, df03nf, df04nf, df05nf, df06nf)
dffin <- rbind(df01fin, df02fin, df03fin, df04fin, df05fin, df06fin)
counts.f <- rbind(counts01.f, counts02.f, counts03.f, counts04.f, counts05.f, counts06.f)
counts.nf <- rbind(counts01.nf, counts02.nf, counts03.nf, counts04.nf, counts05.nf, counts06.nf)
rm(df01f)
rm(df02f)
rm(df03f)
rm(df04f)
rm(df05f)
rm(df06f)
rm(df01nf)
rm(df02nf)
rm(df03nf)
rm(df04nf)
rm(df05nf)
rm(df06nf)
rm(df01fin)
rm(df02fin)
rm(df03fin)
rm(df04fin)
rm(df05fin)
rm(df06fin)
rm(counts01.f)
rm(counts02.f)
rm(counts03.f)
rm(counts04.f)
rm(counts05.f)
rm(counts06.f)
rm(counts01.nf)
rm(counts02.nf)
rm(counts03.nf)
rm(counts04.nf)
rm(counts05.nf)
rm(counts06.nf)

dffin$nf <- dff$model_count
dffin$nt <- dff$model_count+dfnf$model_count
dffin$nf <- as.double(dffin$nf)
dffin$nt <- as.double(dffin$nt)
dffin$nt[dffin$nt==0] <- NA ##IMPORTANT!
corr <- dffin  %>% group_by(model, fit_criteria1, fit_criteria2, fit_criteria3)  %>% do(tidy(cor.test(.$nt, .$nf, alternative="two.sided", method = "spearman"))) 
#counts.f$PR_global <- counts.f$model_count/(counts.f$model_count+counts.nf$model_count)
#counts.f$n_total <- counts.f$model_count+counts.nf$model_count
#counts.f$n_total[13:84] <- NA


###plotting the results on simulated data:
p00<-ggplot(dfrr, aes(x=avg..distance, y=min_S1_T_S2_angle_bin))+
  geom_bin2d(aes(fill = after_stat(count)), binwidth = c(0.5, 10), origin = c(-0.25, -5))+
  scale_fill_gradientn(colours = pals::parula(20), trans = "log", breaks = c(1,10,100,1000,10000,100000), labels = c(1,10,100,1000,10000,100000), na.value='grey95')+
  ggtitle('"10-3 to 10-2" SSL: models to be tested (1,026,076 in total)')+
  scale_x_continuous(breaks=c(1,2,3,4,5,6,7,8,9), expand=c(0.05,0.05))+
  scale_y_continuous(breaks=c(30,60,90,120,150,180), expand=c(0.05,0.05))+
  coord_cartesian(xlim=c(1,9), ylim=c(10,180), expand=T, default=F, clip="on")+
  #xlab("")+
  ylab("")+
  guides(fill = guide_colorbar(title = '', order=1, barwidth=4, barheight=18))+
  theme(plot.title = element_text(size=32, face='bold'),
        axis.text.x=element_blank(),
        axis.text.y=element_text(size=24, face='bold'),
        axis.title.x=element_blank(),
        axis.title.y=element_text(size=36, face='bold'),
        legend.position = "right",
        legend.key=element_rect(fill = 'white'),
        legend.title=element_text(size=36, face = "bold"),
        legend.text=element_text(size=30),
        panel.background=element_rect(fill = 'grey95'),
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        plot.margin = unit(c(10,110,10,10), "pt"),
        strip.text=element_text(size=26, face='bold'))+
  facet_nested(. ~ model)

p01<-ggplot(dffin, aes(x=avg..distance, y=min_S1_T_S2_angle_bin))+
  geom_vline(aes(xintercept=4.75), colour="red3", linewidth=1, alpha=0.5)+
  geom_tile(aes(fill = nf))+
  geom_text(data=counts.f, aes(x=2.5, y=200, label=scales::comma(model_count)), colour="darkgreen", size=10, fontface="bold", check_overlap=T)+
  geom_text(data=counts.nf, aes(x=7, y=200, label=scales::comma(model_count)), colour="red3", size=10, fontface="bold", check_overlap=T)+
  geom_text(data=corr, aes(x=9, y=150, label=round(estimate, 2)), colour="tomato4", size=9, fontface="bold", check_overlap=T, angle=90)+
  geom_hline(aes(yintercept=125), colour="red3", linewidth=1, alpha=0.5)+
  geom_hline(aes(yintercept=95), colour="red3", linewidth=1, alpha=0.5)+
  scale_fill_gradientn(colours = pals::parula(20), trans = "log", breaks = c(1,10,100,1000,10000,100000), labels = c(1,10,100,1000,10000,100000), na.value='grey95')+
  scale_x_continuous(breaks=c(1,2,3,4,5,6,7,8,9), expand=c(0.05,0.05))+
  scale_y_continuous(breaks=c(30,60,90,120,150,180), expand=c(0.05,0.05))+
  coord_cartesian(xlim=c(1,9), ylim=c(10,200), expand=T, default=F, clip="on")+
  ggtitle('counts of feasible 2-3-way models\nin the space of optimality metrics, log. color scale,\nprogressively more stringent feasibility criteria')+
  xlab("average source-target distance")+
  ylab("min. S1-target-S2 angle")+
  guides(fill = guide_colorbar(title = '', order=1, barwidth=4, barheight=24))+
  theme(plot.title = element_text(color="darkmagenta", face="bold", size=28),
        axis.text.x=element_text(size=24, face='bold', angle=45, hjust=0.95, vjust=0.95),
        axis.text.y=element_text(size=24, face='bold'),
        axis.title.x=element_text(size=36, face='bold'),
        axis.title.y=element_text(size=36, face='bold'),
        legend.position = "right",
        legend.key=element_rect(fill = 'white'),
        legend.title=element_text(size=36, face = "bold"),
        legend.text=element_text(size=30),
        panel.background=element_rect(fill = 'white'),
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        plot.margin = unit(c(5,10,33,10), "pt"),
        strip.text=element_text(size=26, face='bold'))+
  facet_nested(fit_criteria3 + factor(fit_criteria1, level=c("% (0,1)", "%+/-2SE (0,1)")) + fit_criteria2 ~ model)

###arranging panels###
pall2 <- plot_grid(
  plot_grid(p00+theme(legend.position="none"), p01+theme(legend.position="none"), ncol = 1, axis="l", align = "h", rel_heights = c(1.2,6.8)),
  plot_grid(get_legend(p00), get_legend(p01), ncol = 1, axis="l", align = "hv", rel_heights = c(1.2,6.8)),
  rel_widths = c(24,6))

pall3 <- plot_grid(pall, pall2, ncol=2, align = "hv", axis="tb", rel_widths = c(1,1), rel_heights = c(1,1))

###This is the right-hand part of Figure 11###
pdf("distributions_avgdist_angles_2D_counts_2-3way_real_vs_sim.pdf", width=24, height=29)
print(pall3)
dev.off()
