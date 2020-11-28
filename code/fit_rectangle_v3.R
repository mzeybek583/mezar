
## M estimation


# Load Libraries ----------------------------------------------------------



## fit rectangle

# Notes -------------------------------------------------------------------

# PCA temelli
# PCA min ve max x ve y leri bul,
# PCA min deðerleri için fit lm
# PCA max deðerleri için fit lm


# Load libraries ----------------------------------------------------------
time <- proc.time()
library(lidR)
library(dbscan)
library(ggplot2)
library(lmodel2)

library(concaveman)
library(rgdal)
library(raster)
library(mclust)

# Read data ---------------------------------------------------------------
path <- "E:/DR_Sonrasi_Projeler/Mezar/R_test/cc_export/merged.las"

#las <- readLAS(files = "test_2.las")
las <- readLAS(files = path)
#plot(las)


# Segmentation ------------------------------------------------------------
# DBSCAN

data <- las@data
xy <- data[,1:2]
xyz <- data[,1:3]



# Export clusters to file -------------------------------------------------
cluster <- data$`Original cloud index`

n.cl <- max(unique(cluster)) +1 #cluster number

df <- list()
poly <- list()
sp <- list()
for (i in 1:n.cl) {
 #i<-1
  cl <- xy[cluster==i]
  df[[i]] <- cl

# PCA implementation ------------------------------------------------------

  data <- df[[i]]
  
  data.pca <- prcomp(data, center = TRUE, scale. = FALSE)
  
  rotated_data = as.data.frame(data.pca$x)
 # gp.pca <- ggplot(data=rotated_data ,aes(PC1,PC2))+
 #   geom_point()+
 #   xlab("PC1")+ ylab("PC2")+
 #   theme_classic(base_size = 20)+ coord_fixed()
 # gp.pca
  

# Model Clustering --------------------------------------------------------


  model.data <- matrix(unlist(rotated_data), ncol=2)
#  plot(model.data)
  mc_rect <- Mclust(model.data, G = 4)
  
#  plot(mc_rect, what = "density")
#  plot(mc_rect, what = "uncertainty")
  

# Order Fits --------------------------------------------------------------


  fit.lm <- list()
  
    j <-   which.min(mc_rect$parameters$mean[2,])
    hat <- model.data[mc_rect$classification==j,]
  #  plot(hat,asp=1)
    hat <- data.frame(hat)
    colnames(hat) <- c("X", "Y")
    fit = lmodel2(Y ~ X, data=hat, "interval", "interval", 99)
 #   plot(fit,method="RMA", asp=1) # major axis regression    
    fit.lm[[1]] <- fit
    
    j <- which.min(mc_rect$parameters$mean[1,])
    hat <- model.data[mc_rect$classification==j,]
#    plot(hat,asp=1)
    hat <- data.frame(hat)
    colnames(hat) <- c("X", "Y")
    fit = lmodel2(Y ~ X, data=hat, "interval", "interval", 99)
#    plot(fit,method="RMA", asp=1) # major axis regression    
    fit.lm[[2]] <- fit
    
    j <- which.max(mc_rect$parameters$mean[2,])
    hat <- model.data[mc_rect$classification==j,]
 #   plot(hat,asp=1)
    hat <- data.frame(hat)
    colnames(hat) <- c("X", "Y")
    fit = lmodel2(Y ~ X, data=hat, "interval", "interval", 99)
 #   plot(fit,method="RMA", asp=1) # major axis regression    
    fit.lm[[3]] <- fit
    
    j <- which.max(mc_rect$parameters$mean[1,])
    hat <- model.data[mc_rect$classification==j,]
  #  plot(hat,asp=1)
    hat <- data.frame(hat)
    colnames(hat) <- c("X", "Y")
    fit = lmodel2(Y ~ X, data=hat, "interval", "interval", 99)
#    plot(fit,method="RMA", asp=1) # major axis regression    
    fit.lm[[4]] <- fit
    
  
  ## Fit lm kesisim noktalarý tespit
  
  cm  <- rbind(cbind(fit.lm[[1]][["regression.results"]][["Intercept"]][2],fit.lm[[1]][["regression.results"]][["Slope"]][2]),
               cbind(fit.lm[[2]][["regression.results"]][["Intercept"]][2],fit.lm[[2]][["regression.results"]][["Slope"]][2])) # Coefficient matrix
  cor.1 <- c(-solve(cbind(cm[,2],-1)) %*% cm[,1])
  cm  <- rbind(cbind(fit.lm[[2]][["regression.results"]][["Intercept"]][2],fit.lm[[2]][["regression.results"]][["Slope"]][2]),
               cbind(fit.lm[[3]][["regression.results"]][["Intercept"]][2],fit.lm[[3]][["regression.results"]][["Slope"]][2])) # Coefficient matrix
  cor.2 <- c(-solve(cbind(cm[,2],-1)) %*% cm[,1])
  cm  <- rbind(cbind(fit.lm[[3]][["regression.results"]][["Intercept"]][2],fit.lm[[3]][["regression.results"]][["Slope"]][2]),
               cbind(fit.lm[[4]][["regression.results"]][["Intercept"]][2],fit.lm[[4]][["regression.results"]][["Slope"]][2])) # Coefficient matrix
  cor.3 <- c(-solve(cbind(cm[,2],-1)) %*% cm[,1])
  cm  <- rbind(cbind(fit.lm[[1]][["regression.results"]][["Intercept"]][2],fit.lm[[1]][["regression.results"]][["Slope"]][2]),
               cbind(fit.lm[[4]][["regression.results"]][["Intercept"]][2],fit.lm[[4]][["regression.results"]][["Slope"]][2])) # Coefficient matrix
  cor.4 <- c(-solve(cbind(cm[,2],-1)) %*% cm[,1])
  
  corners <- data.frame(rbind(cor.1, cor.2, cor.3, cor.4))
  colnames(corners) <- c("X", "Y")
  
  plot.corner <- data.frame(rotated_data)
   colnames(plot.corner) <- c("X", "Y")
 # gp.corner <- ggplot(data=plot.corner ,aes(X,Y))+
     # geom_point()+
     # xlab("X")+ ylab("Y")+
     # theme_classic(base_size = 20)+ coord_fixed() +
     # geom_point(data = corners, col="red")
 #  gp.corner
  
  # PCA convert to original -------------------------------------------------
  
  corners2 <- as.matrix(corners)
  
   orig <- t(t(corners2 %*% t(data.pca$rotation)) + data.pca$center)
   
  # Concave hull algoritmasý ile çevirme ------------------------------------
  
  polygon <- concaveman(orig)
  
  #plot(xy)
 # lines(polygon, type="l", col="red")
  
  # Concave hull sonuclari shapefile ----------------------------------------

    sp[[i]] <- SpatialPolygons(list(Polygons(list(Polygon(polygon)), ID=i)))
  
  
}

joined = SpatialPolygons(lapply(sp, function(x){x@polygons[[1]]}))
prj <- CRS("+init=epsg:5254")

proj4string(joined)<- prj

class(joined)

# Write data to file ------------------------------------------------------
#f.name <- paste("mezar",i,sep = "")
shapefile(joined,"mezars",overwrite=TRUE )
proc.time()-time
# Export to map -----------------------------------------------------------

sprintf("Program Bitti! Toplma sure %3.1f", (proc.time()-time)[3])

