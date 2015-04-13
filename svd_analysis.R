### Name: svd_analysis.R
### Author: Jan Blanke
### Description: Single value decomposition analysis
#####################################################

library(sp)
library(raster)
library(rasterVis)
library(rgdal)

### Settings
data <- "npool_2069_2098" # natural_cmass, total_cpool, total_lai, avg_cflux, natural_avg_cmass, ...

setwd(paste("/media/jan/AddData/Simulations_processed/30y_avg_all/", data, sep=""))
system("ls -l")

### Load data
load("df_out.RData")
load("coords.RData")

df <- df.out
df <- as.matrix(df)

### Combinations
Q <- 2 # 2 rcps: 85 and 45
L <- 3 # 3 gcms: cnrm, mohc, ichec
J <- 2 # 2 projects: CORDEX, CMIP5
K <- 2 # 2 resolutions: 0.5 degree, 10 arcmin

nvars <- Q * L * J * K

# standardize (mean=0, std=1)
dfs <- scale(df)
dfs <- dfs/sqrt(nrow(dfs));

### Single value decomposition
pca <- svd(dfs)
#pca <- svd(df)
U <- pca$u # loadings

# sign correction
Ind <- t(U) %*% rowMeans(dfs) < 0 
U[, Ind] <- -U[, Ind]

### Regression on y with main components
mainc <- rep(1, nrow(U)) / sqrt(nrow(U))
U <- cbind(mainc, U) 

### Regression
library(pracma)
Ureg <- U[, 1:(nvars + 1)]
alpha <- mldivide(Ureg, df)
colnames(alpha) <- NULL 

alphaC <- rowMeans(alpha)
alpha_C <- repmat(a=cbind(alphaC), m=ncol(alpha), n=1)
a_t = alpha - alpha_C;

### Partition of regression coefficients
alpha_1 = matrix(0, nrow(alpha), ncol(alpha))
alpha_2 = matrix(0, nrow(alpha), ncol(alpha))
alpha_3 = matrix(0, nrow(alpha), ncol(alpha))
alpha_4 = matrix(0, nrow(alpha), ncol(alpha))

# 1. rcp
rcp85 <- grep("85", names(df.out))
rcp45 <- grep("45", names(df.out))

# 2. model
cnrm <- grep("CNRM", names(df.out))
mohc <- grep("MOHC", names(df.out))
ichec <- grep("ICHEC", names(df.out))

# 3. project
cmip5 <- grep("CMIP5", names(df.out))
cordex <- grep("CORDEX", names(df.out))

# 4. resolution
min10 <- grep("10min", names(df.out))
min30 <- grep("30min", names(df.out))

# group 1
x <- a_t[, rcp85]
alpha_1[, rcp85] <- repmat(a = cbind(rowMeans(x)), m = length(rcp85), n = 1)
x <- a_t[, rcp45]
alpha_1[, rcp45] <- repmat(a = cbind(rowMeans(x)), m = length(rcp45), n = 1)

# group 2
x <- a_t[, cnrm]
alpha_2[, cnrm] <- repmat(a = cbind(rowMeans(x)), m = length(cnrm), n = 1)
x <- a_t[, mohc]
alpha_2[, mohc] <- repmat(a = cbind(rowMeans(x)), m = length(mohc), n = 1)
x <- a_t[, ichec];
alpha_2[, ichec] <- repmat(a = cbind(rowMeans(x)), m = length(ichec), n = 1)

# group 3
x <- a_t[, cmip5];
alpha_3[, cmip5] <- repmat(a = cbind(rowMeans(x)), m = length(cmip5), n = 1)
x <- a_t[, cordex];
alpha_3[, cordex] <- repmat(a = cbind(rowMeans(x)), m = length(cordex), n = 1)

# group 4
x <- a_t[, min10];
alpha_4[, min10] <- repmat(a = cbind(rowMeans(x)), m = length(min10), n = 1)
x <- a_t[, min30];
alpha_4[, min30] <- repmat(a = cbind(rowMeans(x)), m = length(min30), n=1)

alpha_res <- alpha - alpha_1 - alpha_2 - alpha_3  - alpha_4 - alpha_C; 


### ANOVA Variance decomposition 

# mean over everything
df.mean <- mean(df)
# total sum of squares
SS.T <- sum((df - df.mean) ^ 2)

SS.com <- rep(0, nvars + 1)
SS.ind <- rep(0, nvars + 1)
SS.G1 <- rep(0, nvars + 1)
SS.G1 <- rep(0, nvars + 1)
SS.G2 <- rep(0, nvars + 1)
SS.G3 <- rep(0, nvars + 1)
SS.G4 <- rep(0, nvars + 1)

# common variance
SS.com[1] <- nvars * sum( ((Ureg[, 1] * alpha_C[1, 1]) - df.mean) ^ 2)
for (m in 2:(nvars + 1)) { 
  SS.com[m] <- nvars * sum((Ureg[, m] * alpha_C[m, 1]) ^ 2 )
}

# individual groups
for (m in 1:(nvars + 1)) {
  # G1: 2 vars
  SS.G1[m] <- L * K * J * sum(colSums((Ureg[, m] %o% alpha_1[m, c(rcp45[1], rcp85[1])]) ^ 2)) # rcp      
  # G2: 3 vars  
  SS.G2[m] <- Q * K * J * sum(colSums((Ureg[, m] %o% alpha_2[m, c(cnrm[1], mohc[1], ichec[1])]) ^ 2)) # model  
  # G3: 2 vars  
  SS.G3[m] <- Q * K * L * sum(colSums((Ureg[, m] %o% alpha_3[m, c(cmip5[1], cordex[1])]) ^ 2)) # project  
  # G4: 2 vars  
  SS.G4[m] <- Q * L * J * sum(colSums((Ureg[, m] %o% alpha_4[m, c(min10[1], min30[1])]) ^ 2)) # resolution  
  # Individual
  SS.ind[m] <- sum(colSums((Ureg[, m] %o% alpha_res[m, ]) ^ 2))  
}

# contribution of residuals 
SS.e <- sum(colSums((df - Ureg %*% alpha) ^ 2))

com.eff <- sum(SS.com)/SS.T 
G1.eff <- sum(SS.G1)/SS.T 
G2.eff <- sum(SS.G2)/SS.T 
G3.eff <- sum(SS.G3)/SS.T 
G4.eff <- sum(SS.G4)/SS.T 
ind.eff <- sum(SS.ind)/SS.T 
res.eff <- SS.e/SS.T 

all.eff <- c(com.eff, G1.eff, G2.eff, G3.eff, G4.eff, ind.eff, res.eff)
all.eff
sum(all.eff)
plot(all.eff, col="orange", pch=19)

#fileConn <- file("anova.txt")
#writeLines(as.character(all.eff), fileConn)
#close(fileConn)

#saveRDS(all.eff, file = paste("anova_", data, ".rds", sep=""))

###### Maps of mean approximation and group deviations
#out.stack <- stack()

### Estimation of mean
X_common <- matrix(0, nrow(df), 1) 
X_common <- Ureg %*% alpha_C[, 1]
X_common <- cbind(coords + 0.25, X_common)
X_common_sp <- SpatialPointsDataFrame(X_common[, 1:2], as.data.frame(X_common[, 3]), proj4string=CRS("+proj=longlat +datum=WGS84 +no_defs +ellps=WGS84 +towgs84=0,0,0"))

r <- raster(ncols=(bbox(X_common_sp)[1, 2] - bbox(X_common_sp)[1, 1]) / 0.171, nrows=(bbox(X_common_sp)[2, 2] - bbox(X_common_sp)[2, 1]) / 0.171, xmn=bbox(X_common_sp)[1, 1], xmx=bbox(X_common_sp)[1,2], ymn=bbox(X_common_sp)[2, 1], ymx=bbox(X_common_sp)[2, 2], crs=CRS("+proj=longlat +datum=WGS84 +no_defs +ellps=WGS84 +towgs84=0,0,0"))

out.stack <- stack(rasterize(as.data.frame(X_common)[, 1:2], r, as.data.frame(X_common)[, 3]))

### Estimation of variable RCP
X_ab1 <- matrix(0, nrow(df), 2) 

X_ab1[, 1] <- Ureg %*% alpha_1[, rcp45[1]]
X_ab1[, 2] <- Ureg %*% alpha_1[, rcp85[1]]

# mean of the absolute deviation from the mean
X_amean1 <- rowMeans(abs(X_ab1)) # mean is zero due to standardization
X_amean1 <- cbind(coords + 0.25, X_amean1)

#X_amean1_sp <- SpatialPointsDataFrame(X_amean1[, 1:2], as.data.frame(X_amean1[, 3]), proj4string=CRS("+proj=longlat +datum=WGS84 +no_defs +ellps=WGS84 +towgs84=0,0,0"))
out.stack[[2]] <- rasterize(as.data.frame(X_amean1)[, 1:2], r, as.data.frame(X_amean1)[, 3])


### Estimation of variable GCM
X_ab2 <- matrix(0, nrow(df), 3) 

X_ab2[, 1] <- Ureg %*% alpha_2[, cnrm[1]]
X_ab2[, 2] <- Ureg %*% alpha_2[, mohc[1]]
X_ab2[, 3] <- Ureg %*% alpha_2[, ichec[1]]

# mean of the absolute deviation from the mean
X_amean2 <- rowMeans(abs(X_ab2)) # mean is zero due to standardization
X_amean2 <- cbind(coords + 0.25, X_amean2)

#X_amean2_sp <- SpatialPointsDataFrame(X_amean2[, 1:2], as.data.frame(X_amean2[, 3]), proj4string=CRS("+proj=longlat +datum=WGS84 +no_defs +ellps=WGS84 +towgs84=0,0,0"))
out.stack[[3]] <- rasterize(as.data.frame(X_amean2)[, 1:2], r, as.data.frame(X_amean2)[, 3])

### Estimation of variable PROJECT
X_ab3 <- matrix(0, nrow(df), 2) 

X_ab3[, 1] <- Ureg %*% alpha_3[, cmip5[1]]
X_ab3[, 2] <- Ureg %*% alpha_3[, cordex[1]]

# mean of the absolute deviation from the mean
X_amean3 <- rowMeans(abs(X_ab3)) # mean is zero due to standardization
X_amean3 <- cbind(coords + 0.25, X_amean3)

#X_amean3_sp <- SpatialPointsDataFrame(X_amean3[, 1:2], as.data.frame(X_amean3[, 3]), proj4string=CRS("+proj=longlat +datum=WGS84 +no_defs +ellps=WGS84 +towgs84=0,0,0"))
out.stack[[4]] <- rasterize(as.data.frame(X_amean3)[, 1:2], r, as.data.frame(X_amean3)[, 3])

### Estimation of variable RESOLUTION
X_ab4 <- matrix(0, nrow(df), 2) 

X_ab4[, 1] <- Ureg %*% alpha_4[, min10[1]]
X_ab4[, 2] <- Ureg %*% alpha_4[, min30[1]]

# mean of the absolute deviation from the mean
X_amean4 <- rowMeans(abs(X_ab4)) 
X_amean4 <- cbind(coords + 0.25, X_amean4)

#X_amean4_sp <- SpatialPointsDataFrame(X_amean4[, 1:2], as.data.frame(X_amean4[, 3]), proj4string=CRS("+proj=longlat +datum=WGS84 +no_defs +ellps=WGS84 +towgs84=0,0,0"))
out.stack[[5]] <- rasterize(as.data.frame(X_amean4)[, 1:2], r, as.data.frame(X_amean4)[, 3])

### Estimation of INDIVIDUAL EFFECT
X_ind_no_mean <- matrix(0, nrow(df), ncol(alpha_res)) 
for(indx in 1:ncol(alpha_res)) {
  X_ind_no_mean[, indx] <- Ureg %*% alpha_res[, indx]
}
X_ind <- rowMeans(abs(X_ind_no_mean)) 
X_ind <- cbind(coords + 0.25, X_ind)

#X_ind_sp <- SpatialPointsDataFrame(X_ind[, 1:2], as.data.frame(X_ind[, 3]), proj4string=CRS("+proj=longlat +datum=WGS84 +no_defs +ellps=WGS84 +towgs84=0,0,0"))
out.stack[[6]] <- rasterize(as.data.frame(X_ind)[, 1:2], r, as.data.frame(X_ind)[, 3])

### Save
save(out.stack, file="out_stack_inkl_indiv_eff.RData")
 
