# standard PCA of TraCE area-avaraged data

library(corrplot)
library(qgraph)
library(psych)
library(roll)
library(zoo)


# define functions for producing qgraphs of correlaitons and loadings
qgraph_loadings_plot <- function(loadings_in, title) {
  ld <- loadings(loadings_in)
  qg_pca <- qgraph(ld, title.cex=2.0, title=title, 
                   posCol = "darkgreen", negCol = "darkmagenta", arrows = FALSE, 
                   labels=attr(ld, "dimnames")[[1]])
}

# force-directed version
qgraph_loadings_spring_plot <- function(loadings_in, title) {
  ld <- loadings(loadings_in)
  qg_pca <- qgraph(ld, title.cex=2.0, title=title, 
                   posCol = "darkgreen", negCol = "darkmagenta", arrows = FALSE, 
                   labels=attr(ld, "dimnames")[[1]])

  qgraph(qg_pca, title.cex=2.0, title=title,
         layout = "spring",
         posCol = "darkgreen", negCol = "darkmagenta", arrows = FALSE, label.prop=1.0,
         node.height=0.5, node.width=0.5, vTrans=255, edge.width=0.75, label.cex=1.0,
         width=7, height=5, normalize=TRUE, edge.width=0.75 )
}

# W. Huber's "amend" function for restoring signs of flipping eigenvectors
# https://stats.stackexchange.com/questions/34396/
amend <- function(result) {
  result.m <- as.matrix(result)
  n <- dim(result.m)[1]
  delta <- apply(abs(result.m[-1,] - result.m[-n,]), 1, sum)
  delta.1 <- apply(abs(result.m[-1,] + result.m[-n,]), 1, sum)
  signs <- c(1, cumprod(rep(-1, n-1) ^ (delta.1 <= delta)))
  zoo(result * signs)
}

# path to save figures
fig_path = "../figs/"
# path to save animations
anim_path = "../animations/"

# read area-average anomalies
csv_path <- "../data/cav_files/"
csv_file <- "NAMidCont_WYNH_llsy_monlenadj_seas_anm_PI_hw30_aave_02.csv"

midw_anm <- read.csv(paste(csv_path, csv_file, sep=""))
names(midw_anm)
summary(midw_anm)
midw_anm$qdiv <- midw_anm$qdiv * 10e9
summary(midw_anm)

# rename variables for plotting
names(midw_anm) <- c("Age_ka", "Kex", "Qsw", "Qlw", "Qn", "Qe", "Qh", "Ef", "Pr", "Qdv", "UVQ", "Sh",  "Pw", "Vpd", "PET",      
                     "Eet", "P-E", "Al", "Cwd", "Sm", "R", "dS", "U5", "V5", "Vm5",  "Z5", "Om5", "Cld", "Z7", "Us", "Vs",
                     "SLP", "T85", "T2m", "Ts") 

# Holocene
midw_anm_hol <- midw_anm[midw_anm$Age_ka >= -11.7000, ]
head(midw_anm_hol); tail(midw_anm_hol)
n_anm <- dim(midw_anm_hol)[1]
print(n_anm)

# zoo plot
midw_anm_hol_zoo <- zoo(midw_anm_hol)
# replace index to get values to plot from older to younger
midw_anm_hol_zoo_index <- index(midw_anm_hol_zoo)
head(midw_anm_hol_zoo_index); tail(midw_anm_hol_zoo_index)
index(midw_anm_hol_zoo) <- -1.0 * 11700.0 + midw_anm_hol_zoo_index - 1.0
head(index(midw_anm_hol_zoo)); tail(index(midw_anm_hol_zoo))

# save a .png file
png_file <- paste(fig_path, "midw_anm_hol_zooplot", ".png", sep="")
png(file = png_file, width=1200, height= 1200)
  plot(midw_anm_hol_zoo[, 3:35], xlab="Age (-yr BP)", col="blue", cex.main = 3, font.main=1,
       main="North America Mid-Contient Anomalies (Relative to PI)")
dev.off()


# read area-average local anomalies
csv_path <- "../data/cav_files/"
csv_file <- "NAMidCont_WYNH_llsy_monlenadj_seas_locanm_hw30_aave_02.csv"

midw_locanm <- read.csv(paste(csv_path, csv_file, sep=""))
names(midw_locanm)
summary(midw_locanm)
midw_locanm$qdiv <- midw_locanm$qdiv * 10e9
summary(midw_locanm)

# rename variables for plotting
names(midw_locanm) <- c("Age_ka", "Kex", "Qsw", "Qlw", "Qn", "Qe", "Qh", "Ef", "Pr", "Qdv", "UVQ", "Sh",  "Pw", "Vpd", "PET",      
                     "Eet", "P-E", "Al", "Cwd", "Sm", "R", "dS", "U5", "V5", "Vm5",  "Z5", "Om5", "Cld", "Z7", "Us", "Vs",
                     "SLP", "T85", "T2m", "Ts") 

# Holocene
midw_locanm_hol <- midw_locanm[midw_locanm$Age_ka >= -11.7000, ]
head(midw_locanm_hol); tail(midw_locanm_hol)
n_locanm <- dim(midw_locanm_hol)[1]
print(n_locanm)

# zoo plot
midw_locanm_hol_zoo <- zoo(midw_locanm_hol)
midw_locanm_hol_zoo_index <- index(midw_locanm_hol_zoo)
head(midw_locanm_hol_zoo_index); tail(midw_locanm_hol_zoo_index)
index(midw_locanm_hol_zoo) <- -1.0 * 11700.0 + midw_locanm_hol_zoo_index - 1.0
head(index(midw_locanm_hol_zoo)); tail(index(midw_locanm_hol_zoo))

# save a .png file
png_file <- paste(fig_path, "midw_locanm_hol_zooplot", ".png", sep="")
png(file = png_file, width=1200, height= 1200)
  plot(midw_locanm_hol_zoo[, 3:35], xlab="Age (-yr BP)", col="blue", cex.main = 3.0, font.main=1,
       main="North America Mid-Contient Local Anomalies (Relative to PI)")
dev.off()

# global correlations

# corrplot of anomalies
png_file <- paste(fig_path, "corrplot_midw_anm_hol", ".png", sep="")
png(file = png_file, width=1200, height= 1200)
  cor_midw_anm_hol <- cor(midw_anm_hol[3:35])
  corrplot(cor_midw_anm_hol, method="color")
dev.off()

# corrplot of local anomalies
png_file <- paste(fig_path, "corrplot_midw_locanm_hol", ".png", sep="")
png(file = png_file, width=1200, height= 1200)
  cor_midw_locanm_hol <- cor(midw_locanm_hol[3:35])
  corrplot(cor_midw_locanm_hol, method="color")
dev.off()

# qgraph of anomalies
png_file <- paste(fig_path, "qgraph_midw_anm_hol", ".png", sep="")
png(file = png_file, width=1200, height= 1200)
  title = "Correlations, North American Mid-Continent Anomalies, Holocene"
  qgraph(cor_midw_anm_hol, title.cex=2, title=title,
         # layout = "spring", repulsion = 0.75,
         posCol = "darkgreen", negCol = "darkmagenta", arrows = FALSE, label.prop=1.0,
         node.height=0.5, node.width=0.5, vTrans=255, edge.width=0.75, label.cex=1.0,
         width=7, height=5, normalize=TRUE, edge.width=0.75) 
dev.off()

# qgraph of local anomalies
png_file <- paste(fig_path, "qgraph_midw_locanm_hol", ".png", sep="")
png(file = png_file, width=1200, height= 1200)
title = "Correlations, North American Mid-Continent Local Anomalies, Holocene"
qgraph(cor_midw_locanm_hol, title.cex=2, title=title,
       # layout = "spring", repulsion = 0.75,
       posCol = "darkgreen", negCol = "darkmagenta", arrows = FALSE, label.prop=1.0,
       node.height=0.5, node.width=0.5, vTrans=255, edge.width=0.75, label.cex=1.0,
       width=7, height=5, normalize=TRUE, edge.width=0.75 ) 
dev.off()

# force-directed version

# qgraph of anomalies
png_file <- paste(fig_path, "qgraph_spring_midw_anm_hol", ".png", sep="")
png(file = png_file, width=1200, height= 1200)
title = "Correlations, North American Mid-Continent Anomalies, Holocene"
qgraph(cor_midw_anm_hol, title.cex=2, title=title,
       layout = "spring", repulsion = 0.75,
       posCol = "darkgreen", negCol = "darkmagenta", arrows = FALSE, label.prop=1.0,
       node.height=0.5, node.width=0.5, vTrans=255, edge.width=0.75, label.cex=1.0,
       width=7, height=5, normalize=TRUE, edge.width=0.75) 
dev.off()

# qgraph of local anomalies
png_file <- paste(fig_path, "qgraph_spring_midw_locanm_hol", ".png", sep="")
png(file = png_file, width=1200, height= 1200)
title = "Correlations, North American Mid-Continent Local Anomalies, Holocene"
qgraph(cor_midw_locanm_hol, title.cex=2, title=title,
       layout = "spring", repulsion = 0.75,
       posCol = "darkgreen", negCol = "darkmagenta", arrows = FALSE, label.prop=1.0,
       node.height=0.5, node.width=0.5, vTrans=255, edge.width=0.75, label.cex=1.0,
       width=7, height=5, normalize=TRUE, edge.width=0.75 ) 
dev.off()


# rolling correlations

# window width and minimum number of observations
width <- 400
min_obs <- 200 # rolling correlations will start at 11.5 ka

# weights -- half Gaussian
x <- seq(-3, 0, by=1/width)
head(x); tail(x)
# weights <- rep(1, width) 
weights <- dnorm(x, 0, 1)
# weights <- rep(1.0, width)
plot(weights)

# rolling correlations -- anomalies

ptm <- proc.time()
midw_anm_hol_roll_corr <- roll_cor(as.matrix(midw_anm_hol[, 3:35]), width=width, weights=weights, min_obs = min_obs, online = FALSE)
proc.time() - ptm

# check
midw_anm_hol_roll_corr[, , 5000]
plot(midw_anm_hol_roll_corr["Pr", "T2m", ])
plot.ts(cbind(midw_anm_hol["Pr"], midw_anm_hol["T2m"]), type = "l")

# plot qgraphs of rolling correlations
ptm <- proc.time()
# i <- 10000
for (i in seq(201, n_anm, by=100)) {
  
  age <- sprintf("%5.3f", -1 * midw_anm_hol$Age_ka[i])
  if ((-1 * midw_anm_hol$Age_ka[i]) < 10) age <- paste("0", age, sep="")
  if ((-1 * midw_anm_hol$Age_ka[i]) == 0) age <- paste("0", sprintf("%5.3f", midw_anm_hol$Age_ka[i]), sep="")
  if ((-1 * midw_anm_hol$Age_ka[i]) < 0) age <- paste("+", sprintf("%5.3f", midw_anm_hol$Age_ka[i]), sep="")
  print(age)
  
  png_file <- paste(anim_path, "rc_anm/rc_anm_", age, ".png", sep="")
  png(file = png_file, width=1200, height= 1200)
  title = paste("Correlations, North American Mid-Continent Anomalies -- ", age, " ka", sep="")
  qgraph(midw_anm_hol_roll_corr[, , i], title.cex = 2.0, title=title,
         # layout = "spring", repulsion = 0.75,
         posCol = "darkgreen", negCol = "darkmagenta", arrows = FALSE, label.prop=1.0,
         node.height=0.5, node.width=0.5, vTrans=255, edge.width=0.75, label.cex=1.0,
         width=7, height=7, normalize=TRUE, edge.width=0.75, threshold=0.4 )
  dev.off()
}
proc.time() - ptm

# rolling correlations -- local anomalies

ptm <- proc.time()
midw_locanm_hol_roll_corr <- roll_cor(as.matrix(midw_locanm_hol[, 3:35]), width=width, weights=weights, min_obs = min_obs, online = FALSE)
proc.time() - ptm

# check
midw_locanm_hol_roll_corr[, , 5000]
plot(midw_locanm_hol_roll_corr["Pr", "T2m", ])
plot.ts(cbind(midw_locanm_hol["Pr"], midw_locanm_hol["T2m"]), type = "l")

# loc anomalies
ptm <- proc.time()
# i <- 10000
for (i in seq(201, n_locanm, by=100)) {
  
  age <- sprintf("%5.3f", -1 * midw_locanm_hol$Age_ka[i])
  if ((-1 * midw_locanm_hol$Age_ka[i]) < 10) age <- paste("0", age, sep="")
  if ((-1 * midw_locanm_hol$Age_ka[i]) == 0) age <- paste("0", sprintf("%5.3f", midw_locanm_hol$Age_ka[i]), sep="")
  if ((-1 * midw_locanm_hol$Age_ka[i]) < 0) age <- paste("+", sprintf("%5.3f", midw_locanm_hol$Age_ka[i]), sep="")
  print(age)
  
  png_file <- paste(anim_path, "rc_locanm/rc_locanm_", age, ".png", sep="")
  png(file = png_file, width=1200, height= 1200)
  title = paste("Correlations, North American Mid-Continent Local Anomalies -- ", age, " ka", sep="")
  qgraph(midw_locanm_hol_roll_corr[, , i], title.cex = 2.0, title=title,
         # layout = "spring", repulsion = 0.75,
         posCol = "darkgreen", negCol = "darkmagenta", arrows = FALSE, label.prop=1.0,
         node.height=0.5, node.width=0.5, vTrans=255, edge.width=0.75, label.cex=1.0,
         width=7, height=7, normalize=TRUE, edge.width=0.75, threshold=0.4 )
  dev.off()
}
proc.time() - ptm

# PCA

# global PCA of anomalies
nfactors <- 4
midw_anm_hol_pca_unrot <- principal(midw_anm_hol[, 3:35], nfactors = nfactors, rotate = "none")
midw_anm_hol_pca_unrot

# qgrapgh of of loadings
png_file <- paste(fig_path, "global_pca_midw_anm_hol", ".png", sep="")
png(file = png_file, width=1200, height= 1200)
  qgraph_loadings_plot(midw_anm_hol_pca_unrot, title="Unrotated component loadings, North American Mid-Continent Anomalies")
dev.off()

# global PCA of local anomalies
nfactors <- 4
midw_locanm_hol_pca_unrot <- principal(midw_locanm_hol[, 3:35], nfactors = nfactors, rotate = "none")
midw_locanm_hol_pca_unrot

# qgraph of loadings
png_file <- paste(fig_path, "global_pca_midw_locanm_hol", ".png", sep="")
png(file = png_file, width=1200, height= 1200)
qgraph_loadings_plot(midw_locanm_hol_pca_unrot, title="Unrotated component loadings, North American Mid-Continent Local Anomalies")
dev.off()


# rolling PCA with rollapply() and princomo()

# rolling PCA of anomalies

# zoo rollapply()
window <- 200
ptm <- proc.time()
roll_load_midw_anm_hol <- rollapply(zoo(midw_anm_hol[, 3:35]), window, 
                           function(x) (princomp(x))$loadings[, 1], by.column = FALSE, align = "right")
proc.time() - ptm

# check for flipping signs
plot(roll_load_midw_anm_hol)
roll_load_midw_anm_hol <- amend(roll_load_midw_anm_hol)
plot(roll_load_midw_anm_hol)

# better plot
roll_load_midw_anm_hol_index <- index(roll_load_midw_anm_hol)
index(roll_load_midw_anm_hol) <- -1.0 * 11700.0 + roll_load_midw_anm_hol_index - 1.0

# qgraph
png_file <- paste(fig_path, "roll_load_midw_anm_hol", ".png", sep="")
png(file = png_file, width=1200, height= 1200)
  plot(roll_load_midw_anm_hol, xlab="Age (-ka)", col="red", cex.main = 3, font.main=1,
       main="Rolling PCA -- Component 1 Loadings -- North America Mid-Contient Anomalies (Relative to PI)")
dev.off()

# rolling eigenvalues
window <- 200
rolling_eigen_midw_anm_hol <- rollapply(zoo(midw_anm_hol[, 3:35]), window, 
                    function(x) (princomp(x))$sdev[1:4], by.column = FALSE, align = "right")
plot(rolling_eigen_midw_anm_hol)

# PCA of PCA
nfactors <- 3
pca_roll_load_midw_anm_hol <- principal(roll_load_midw_anm_hol, nfactors = nfactors, rotate = "none")
pca_roll_load_midw_anm_hol
loadings(pca_roll_load_midw_anm_hol)

# qgrapgh
png_file <- paste(fig_path, "pca_roll_load_midw_anm_hol", ".png", sep="")
png(file = png_file, width=1200, height= 1200)
  qgraph_loadings_plot(pca_roll_load_midw_anm_hol, 
                       title="PCA of Rolling PCA Loadings[1] -- North America Mid-Contient Anomalies (Relative to PI)")
dev.off()

# rolling PCA of local anomalies

# zoo rollapply()
window <- 400
ptm <- proc.time()
roll_load_midw_locanm_hol <- rollapply(zoo(midw_locanm_hol[, 3:35]), window, 
                                    function(x) (princomp(x))$loadings[, 1], by.column = FALSE, align = "right")
proc.time() - ptm

# check for flipping signs
plot(roll_load_midw_locanm_hol)
roll_load_midw_locanm_hol <- amend(roll_load_midw_locanm_hol)
plot(roll_load_midw_locanm_hol)

# better plot
roll_load_midw_locanm_hol_index <- index(roll_load_midw_locanm_hol)
index(roll_load_midw_locanm_hol) <- -1.0 * 11700.0 + roll_load_midw_locanm_hol_index - 1.0

# qraph
png_file <- paste(fig_path, "pca_roll_load_midw_locanm_hol", ".png", sep="")
png(file = png_file, width=1200, height= 1200)
  plot(roll_load_midw_locanm_hol, xlab="Age (-ka)", col="red", cex.main = 3, font.main=1,
       main="Rolling PCA -- Component 1 Loadings -- North America Mid-Contient Local Anomalies (Relative to PI)")
dev.off()

# rolling eigenvalues
window <- 400
rolling_eigen_midw_locanm_hol <- rollapply(zoo(midw_locanm_hol[, 3:35]), window, 
                                        function(x) (princomp(x))$sdev[1:4], by.column = FALSE, align = "right")
plot(rolling_eigen_midw_locanm_hol)

# PCA of PCA
nfactors <- 4
pca_roll_load_midw_locanm_hol <- principal(roll_load_midw_locanm_hol, nfactors = nfactors, rotate = "none")
pca_roll_load_midw_locanm_hol
loadings(pca_roll_load_midw_locanm_hol)

# qgraph
png_file <- paste(fig_path, "pca_roll_load_midw_locanm_hol", ".png", sep="")
png(file = png_file, width=1200, height= 1200)
  qgraph_loadings_plot(pca_roll_load_midw_locanm_hol, 
    title="PCA of Rolling PCA Loadings[1] -- North America Mid-Contient Local Anomalies (Relative to PI)")
dev.off()