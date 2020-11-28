
## fit rectangle

# Notes -------------------------------------------------------------------

# PCA temelli
# PCA min ve max x ve y leri bul,
# PCA min deðerleri için fit lm
# PCA max deðerleri için fit lm


# Load libraries ----------------------------------------------------------

library(lidR)
library(dbscan)
library(ggplot2)
library(lmodel2)

library(concaveman)
library(rgdal)
library(raster)

# Read data ---------------------------------------------------------------
path <- "E:/DR_Sonrasi_Projeler/GeoSlam_Veri/Geomatics_geoslam_veriler/Mezar/las_Data/2020-11-07_04-58-39_seferihisar-rev-koordinatli-1.las"

#las <- readLAS(files = "test_2.las")
las <- readLAS(files = path)
#plot(las)


# Ground Segmentation -----------------------------------------------------

mycsf <- csf(FALSE, 0.1,0.5, time_step = 0.65)
las <- classify_ground(las, mycsf)

#plot(las, color = "Classification")

#dtm = grid_terrain(las, res = 0.5, algorithm = tin())

#plot(dtm)

las <- normalize_height(las, tin())

#plot(las)

# Cross section -----------------------------------------------------------

las.crs <- filter_poi(las, Z>=0.2 & Z<=0.3)

#plot(las.crs)
# Segmentation ------------------------------------------------------------
# DBSCAN

df <- las.crs@data
xy <- df[,1:2]
xyz <- df[,1:3]
## find suitable eps parameter using a k-NN plot for k = dim + 1
## Look for the knee!

kNNdistplot(xy, k = 5)
abline(h=.18, col = "red", lty=2)

res <- dbscan(xy, eps = .1, minPts = 5)
res

#pairs(xy, col = res$cluster + 1L)

## use precomputed frNN
fr <- frNN(xy, eps = .1)
res2 <- dbscan(xy, eps=0.1, minPts = 50)

res2
#pairs(xy, col = res2$cluster + 1L)

# Export clusters to file -------------------------------------------------
df.export <- data.frame(xyz)
las.cls.export <- LAS(df.export)

scalar.field <- add_lasattribute(las.cls.export, x=res2$cluster, name="cluster", desc = "Class")
print(scalar.field)
print(scalar.field@header)

writeLAS(scalar.field,"clusters.las")

n.cl <- max(unique(res2$cluster)) #cluster number

df <- list()

for (i in 1:n.cl) {
  
  cl <- xy[res2$cluster==i]
  df [[i]] <- cl
  
}
# PCA imply ---------------------------------------------------------------

data <- df[[4]]

data.pca <- prcomp(data, center = TRUE, scale. = FALSE)

rotated_data = as.data.frame(data.pca$x)

# Plot
gp <- ggplot(data=data ,aes(X,Y))+
  geom_point()+
  xlab("X")+ ylab("Y")+
  theme_classic(base_size = 20)
gp

gp.pca <- ggplot(data=rotated_data ,aes(PC1,PC2))+
  geom_point()+
  xlab("PC1")+ ylab("PC2")+
  theme_classic(base_size = 20)+ coord_fixed()
gp.pca

boxplot(rotated_data)

# Fit LM kesisim noktalari ------------------------------------------------

wd <- 0.2

q1.min <- quantile(rotated_data$PC1)
ind.1 <- rotated_data$PC1 >= q1.min[1] & rotated_data$PC1 <= (q1.min[1]+wd)
x <- rotated_data$PC1[ind.1]
y <- rotated_data$PC2[ind.1]
plot(x,y)
df.1 <- data.frame(x,y)


fit1 = lmodel2(y ~ x, data=df.1)
plot(fit1,method="MA", asp=1) # major axis regression

q1.max <- quantile(rotated_data$PC1)
ind.2 <- rotated_data$PC1 >= (q1.max[5]-wd) & rotated_data$PC1 <= (q1.max[5])
x <- rotated_data$PC1[ind.2]
y <- rotated_data$PC2[ind.2]
plot(x,y, asp=1)

df.2 <- data.frame(x,y)
fit2 = lmodel2(y ~ x, data=df.2)

plot(fit2,method="MA", asp=1) # major axis regression

q2.min <- quantile(rotated_data$PC2)
ind.1 <- rotated_data$PC2 >= q2.min[1] & rotated_data$PC2 <= (q2.min[1]+wd)
x <- rotated_data$PC1[ind.1]
y <- rotated_data$PC2[ind.1]
df.3 <- data.frame(x,y)

fit3 = lmodel2(y ~ x, data=df.3)

plot(fit3,method="MA", asp=1) # major axis regression

q2.max <- quantile(rotated_data$PC2)
ind.2 <- rotated_data$PC2 >= (q2.max[5]-wd) & rotated_data$PC2 <= q2.max[5]
x <- rotated_data$PC1[ind.2]
y <- rotated_data$PC2[ind.2]
plot(x,y, asp=1)

df.4 <- data.frame(x,y)

fit4 = lmodel2(y ~ x, data=df.4)

plot(fit4,method="MA", asp=1) # major axis regression

## Fit lm kesisim noktalarý tespit

cm  <- rbind(cbind(fit1$regression.results$Intercept[2],fit1$regression.results$Slope[2]),
             cbind(fit3$regression.results$Intercept[2],fit3$regression.results$Slope[2])) # Coefficient matrix
cor.1 <- c(-solve(cbind(cm[,2],-1)) %*% cm[,1])
cm  <- rbind(cbind(fit1$regression.results$Intercept[2],fit1$regression.results$Slope[2]),
             cbind(fit4$regression.results$Intercept[2],fit4$regression.results$Slope[2])) # Coefficient matrix
cor.2 <- c(-solve(cbind(cm[,2],-1)) %*% cm[,1])
cm  <- rbind(cbind(fit2$regression.results$Intercept[2],fit2$regression.results$Slope[2]),
             cbind(fit3$regression.results$Intercept[2],fit3$regression.results$Slope[2])) # Coefficient matrix
cor.3 <- c(-solve(cbind(cm[,2],-1)) %*% cm[,1])
cm  <- rbind(cbind(fit2$regression.results$Intercept[2],fit2$regression.results$Slope[2]),
             cbind(fit4$regression.results$Intercept[2],fit4$regression.results$Slope[2])) # Coefficient matrix
cor.4 <- c(-solve(cbind(cm[,2],-1)) %*% cm[,1])

corners <- data.frame(rbind(cor.1, cor.2, cor.3, cor.4))
colnames(corners) <- c("PC1", "PC2")

gp.pca <- ggplot(data=rotated_data ,aes(PC1,PC2))+
  geom_point()+
  xlab("X")+ ylab("Y")+
  theme_classic(base_size = 20)+ coord_fixed() +
  geom_point(data = corners, col="red")
gp.pca

# PCA convert to original -------------------------------------------------

corners2 <- as.matrix(corners)

orig <- t(t(corners2 %*% t(data.pca$rotation)) + data.pca$center)

# Concave hull algoritmasý ile çevirme ------------------------------------

polygon <- concaveman(orig)

plot(xy)
lines(polygon, type="l", add=T, col="red")

# Concave hull sonuclari shapefile ----------------------------------------

prj <- CRS("+init=epsg:5254")

sp <- SpatialPolygons(list(Polygons(list(Polygon(polygon)), ID=1)))

proj4string(sp)<- prj

class(sp)

# Write data to file ------------------------------------------------------


shapefile(sp,"mezarboundary",overwrite=TRUE )

# Export to map -----------------------------------------------------------



