# SVD/PCA of TraCE data

# load packages and set paths
library(ncdf4)
library(RSpectra)
# library(RColorBrewer)
# library(lattice)
# library(rnaturalearth)

# # get a world outline
# world <- ne_coastline(scale = 110, returnclass = "sp")

# set paths
ncpath <- "../data/nc_files/"
file_label <- "TraCE_monthly"
out_path <- "../data/SVD/"

# TraCE-21k monthly time series"

# # actual values
# ncname <- "TraCE_tas_lln_monlenadj.nc"
# dname <- "tas"

# # anomalies
# ncname <- "TraCE_tas_lln_monlenadj_anm_PI_hw30.nc"
# dname <- "tas_anm"

# local anomalies
ncname <- "TraCE_tas_lln_monlenadj_locanm_hw30.nc"
dname <- "tas_locanm"

fname <- paste(file_label, "RSpectra", dname, sep="_")
scoresfile <- paste(out_path, fname, "_scores.csv", sep="")
statsfile <- paste(out_path, fname, "_stats.csv", sep="")
ncloadings <- paste(out_path, fname, "_loadings.nc", sep="")



# read data

# open the netCDF file
ncfname <- paste(ncpath, ncname, sep="")
ncin <- nc_open(ncfname)
print(ncin)

# lons and lats
lon <- ncvar_get(ncin, "lon")
nlon <- dim(lon)
head(lon)

lat <- ncvar_get(ncin, "lat", verbose = F)
nlat <- dim(lat)
head(lat)

# time variable
t <- ncvar_get(ncin, "time")
head(t); tail(t)

tunits <- ncatt_get(ncin, "time", "units")
tunits$value
nt <- dim(t)
nt
plt_xvals <- t

# get the data array
var_array <- ncvar_get(ncin, dname)
dlname <- ncatt_get(ncin, dname, "long_name")
dunits <- ncatt_get(ncin, dname, "units")
fillvalue <- ncatt_get(ncin, dname, "_FillValue")
dim(var_array)

# global attributes
title <- ncatt_get(ncin, 0, "title")
institution <- ncatt_get(ncin, 0, "institution")
datasource <- ncatt_get(ncin, 0, "source")
references <- ncatt_get(ncin, 0, "references")
history <- ncatt_get(ncin, 0, "history")
Conventions <- ncatt_get(ncin, 0, "Conventions")
nc_close(ncin)

# # plot a slice of data
# i <- 1
# var_slice <- matrix(NA, nlon, nlat)
# var_slice[1:(nlon/2), ] <- var_array[((nlon/2)+1):nlon, , i]
# var_slice[((nlon/2)+1):nlon, ] <- var_array[1:(nlon/2), , i]
# lon2 <- rep(0, nlon)
# lon2[1:(nlon/2)] <- lon[((nlon/2)+1):nlon] - 360.0
# lon2[((nlon/2)+1):nlon] <- lon[1:(nlon/2)]
# 
# # levelplot of the slice
# grid <- expand.grid(lon2=lon2, lat=lat)
# cutpts <- c(-60,-10,-2,-1,0,1,2,5,10,50)
# lp <- levelplot(var_slice ~ lon2 * lat, data=grid, at=cutpts, cuts=11, pretty=T, 
#           col.regions=(rev(brewer.pal(10,"RdBu"))))
# lp + latticeExtra::layer(sp.lines(world, col="black", lwd=0.5))

## # (not run) trim data to N.H.
## lat <- lat[47:91]
## nlat <- dim(lat)
## min(lat); max(lat)
## var_array <- var_array[,47:91,]

# reshape the 3-d array
var_vec_long <- as.vector(var_array)
length(var_vec_long)
X <- t(matrix(var_vec_long, nrow = nlon * nlat, ncol = nt))
dim(X)

# set n, p, and ncomp
n <- dim(X)[1]
p <- dim(X)[2]
ncomp <- 10
print (c(n, p, ncomp, n * p))

# SCD via RSpectra
time_RSpectra <- proc.time()
system.time( Z <- as.matrix(scale(X)) )
system.time( svd_Z <- svds(Z, k=ncomp, nu=ncomp, nv=ncomp) )
time_RSpectra <- proc.time() - time_RSpectra
time_RSpectra

# eigenvalues
eigenvalues <- svd_Z$d / sqrt(n-1)
eigenvalues

# loadings (v)
loadpca <- svd_Z$v * t(matrix(rep(eigenvalues, ncomp * p), nrow=ncomp, ncol=p))
head(loadpca); tail(loadpca)

# scores
scores <- (svd_Z$u*sqrt(n-1)) # svd_A$u #
colnames(scores) <- paste("S", as.character(seq(1, ncomp, by=1)), sep="")
head(scores); tail(scores)
apply(scores, 2, mean); apply(scores, 2, sd); apply(scores, 2, range)

# save results
eigenvalues_RSpectra <- eigenvalues
loadpca_RSpectra <- loadpca
scores_RSpectra <- scores
# # 
# # plot a slice of data
# i <- 1
# var_slice <- matrix(NA, nlon, nlat)
# dim(var_slice)
# var_slice[1:(nlon/2), ] <- var_array[((nlon/2)+1):nlon, , i]
# var_slice[((nlon/2)+1):nlon, ] <- var_array[1:(nlon/2), , i]
# lon2 <- rep(0, nlon)
# lon2[1:(nlon/2)] <- lon[((nlon/2)+1):nlon] - 360.0
# lon2[((nlon/2)+1):nlon] <- lon[1:(nlon/2)]
# 
# # levelplot of the slice
# grid <- expand.grid(lon2=lon2, lat=lat)
# cutpts <- c(-60,-10,-2,-1,0,1,2,5,10,50)
# lp <- levelplot(var_slice ~ lon2 * lat, data=grid, at=cutpts, cuts=11, pretty=T, 
#                 col.regions=(rev(brewer.pal(10,"RdBu"))))
# lp + latticeExtra::layer(sp.lines(world, col="black", lwd=0.5))
# 
# 
# # quick map of component loadings
# i <- 1
# var_slice <- matrix(loadpca_RSpectra[, i], ncol=nlon, nrow=nlat, byrow=TRUE)
# dim(var_slice)
# grid <- expand.grid(lon=lon, lat=lat)
# cutpts <- c(-1,-.8,-.6,-.4,-.2,0,.2,.4,.6,.8,1.0)
# lp <-levelplot(t(var_slice) ~ lon * lat, data=grid, at=cutpts, cuts=11, pretty=T,
#           col.regions=(rev(brewer.pal(10,"RdBu"))))
# lp
# 
# # quick map of component loadings
# i <- 1
# var_slice <- (matrix(loadpca_RSpectra[, i], ncol=nlon, nrow=nlat, byrow=TRUE))
# dim(var_slice)
# # var_slice[, 1:(nlon/2)] <- var_slice[, ((nlon/2)+1):nlon]
# # var_slice[, ((nlon/2)+1):nlon] <- var_slice[, 1:(nlon/2)]
# lon2 <- rep(0, nlon)
# lon2[1:(nlon/2)] <- lon[((nlon/2)+1):nlon] - 360.0
# lon2[((nlon/2)+1):nlon] <- lon[1:(nlon/2)]
# cutpts <- c(-1,-.8,-.6,-.4,-.2,0,.2,.4,.6,.8,1.0)
# lp <- levelplot(t(var_slice) ~ lon * lat, data=grid, at=cutpts, cuts=11, pretty=T,
#           col.regions=(rev(brewer.pal(10,"RdBu"))))
# lp
# lp + latticeExtra::layer(sp.lines(world, col="black", lwd=0.5))
# 
# # timeseries plot of scores
# plot(plt_xvals, scores[,5], type="l")

# write output

# RSpectra
# fname <- paste(file_label, "RSpectra", sep="_")
eigenvalues <- eigenvalues_RSpectra
loadpca <- loadpca_RSpectra
scores <- scores_RSpectra

# statistics
prop_EV <- (eigenvalues^2)/p
cum_prop_EV <- cumsum(prop_EV)
statsout <- cbind(eigenvalues,prop_EV,cum_prop_EV)
statsout
write.table(statsout, file=statsfile, row.names=TRUE, col.names=TRUE, sep=",")


# write scores
scoresout <- cbind(plt_xvals,scores)
write.table(scoresout, file=scoresfile, row.names=FALSE, col.names=TRUE, sep=",")

# netCDF file of loadings

# reshape loadings
dim_array <- c(nlon,nlat,ncomp)
loadpca <- array(loadpca, dim_array)

# define dimensions
londim <- ncdim_def("lon", "degrees_east", as.double(lon))
latdim <- ncdim_def("lat", "degrees_north", as.double(lat))
comp <- seq(1:ncomp)
ncompdim <- ncdim_def("comp", "SVD component", as.integer(comp))

# define variable
fillvalue <- 1e+32
dlname <- paste(dname, "loadings", sep=" ")
var_def <- ncvar_def(paste(dname,"loadings", sep="_"), "1", list(londim, latdim, ncompdim), fillvalue, 
  dlname, prec = "single")

# create netCDF file and put arrays
ncout <- nc_create(paste(ncloadings, sep=""), list(var_def), force_v4 = T)

# put loadings
ncvar_put(ncout, var_def, loadpca)

# put additional attributes into dimension and data variables
ncatt_put(ncout, "lon", "axis", "X")  #,verbose=FALSE) #,definemode=FALSE)
ncatt_put(ncout, "lat", "axis", "Y")
ncatt_put(ncout, "comp", "axis", "PC")

# add global attributes
title2 <- paste(title$value, "SVD component analysis using pcaMethods", sep="--")
ncatt_put(ncout, 0, "title", title2)
ncatt_put(ncout, 0, "institution", institution$value)
ncatt_put(ncout, 0, "source", datasource$value)
ncatt_put(ncout, 0, "references", references$value)
history <- paste("P.J. Bartlein", date(), sep = ", ")
ncatt_put(ncout, 0, "history", history)
ncatt_put(ncout, 0, "Conventions", Conventions$value)
nc_close(ncout)

