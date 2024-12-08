---
title: spatial data visualization
author: Wilfried Dossou-Yovo
date: '2024-12-07'
slug: spatial-data-visualization
categories: []
tags: []
---

# Load libraries

``` r
#install.packages("blogdown")
#blogdown::new_site()
library(blogdown)
library(tidyverse)
```

```
## ── Attaching core tidyverse packages ──────────────────────── tidyverse 2.0.0 ──
## ✔ dplyr     1.1.4     ✔ readr     2.1.5
## ✔ forcats   1.0.0     ✔ stringr   1.5.1
## ✔ ggplot2   3.5.1     ✔ tibble    3.2.1
## ✔ lubridate 1.9.3     ✔ tidyr     1.3.1
## ✔ purrr     1.0.2     
## ── Conflicts ────────────────────────────────────────── tidyverse_conflicts() ──
## ✖ dplyr::filter() masks stats::filter()
## ✖ dplyr::lag()    masks stats::lag()
## ℹ Use the conflicted package (<http://conflicted.r-lib.org/>) to force all conflicts to become errors
```

``` r
library(sf)
```

```
## Linking to GEOS 3.10.2, GDAL 3.4.2, PROJ 8.2.1; sf_use_s2() is TRUE
```

``` r
library(terra)
```

```
## terra 1.7.65
## 
## Attaching package: 'terra'
## 
## The following object is masked from 'package:tidyr':
## 
##     extract
```

``` r
library(mapview)
library(geoR)
```

```
## --------------------------------------------------------------
##  Analysis of Geostatistical Data
##  For an Introduction to geoR go to http://www.leg.ufpr.br/geoR
##  geoR version 1.9-3 (built on 2023-12-11) is now loaded
## --------------------------------------------------------------
```

``` r
library(gstat)
library(spData)
```

```
## To access larger datasets in this package, install the spDataLarge
## package with: `install.packages('spDataLarge',
## repos='https://nowosad.github.io/drat/', type='source')`
```

``` r
library(tmap)
```

```
## Breaking News: tmap 3.x is retiring. Please test v4, e.g. with
## remotes::install_github('r-tmap/tmap')
```

``` r
library(sp)
library(viridis)
```

```
## Loading required package: viridisLite
```


``` r
data(meuse)
data(meuse.grid)
meuse <- st_as_sf(meuse, coords = c("x", "y"), crs = 28992)
meuse.grid <- st_as_sf(meuse.grid, coords = c("x", "y"),
                       crs = 28992)
```

# Kriging
## Anisotropy testing

``` r
# Compute directional variograms
v_dir <- variogram(log(zinc) ~ 1, meuse, alpha = c(0, 45, 90, 135))

# Plot the directional variograms
plot(v_dir)
```

<img src="{{< blogdown/postref >}}index_files/figure-html/unnamed-chunk-3-1.png" width="672" />

## Define anisotropic variogram model

``` r
# Define anisotropic variogram model
anisotropic_model <- vgm(psill = 0.8, model = "Sph", range = 1000, nugget = 0.1,
                         anis = c(45, 0.6))  # 45° major axis, anisotropy ratio 0.6
v <- variogram(log(zinc) ~ 1, data = meuse)
plot(v)
```

<img src="{{< blogdown/postref >}}index_files/figure-html/unnamed-chunk-4-1.png" width="672" />

``` r
# Fit the anisotropic variogram model
anisotropic_fitted <- fit.variogram(v, model = anisotropic_model)

# Plot the fitted anisotropic variogram
plot(v, anisotropic_fitted, cutoff = 1500, cex = 1.5)
```

<img src="{{< blogdown/postref >}}index_files/figure-html/unnamed-chunk-4-2.png" width="672" />

## kriging

``` r
k <- gstat(formula = log(zinc) ~ 1, data = meuse, model = anisotropic_model)
kpred <- predict(k, meuse.grid)
```

```
## [using ordinary kriging]
```


``` r
# Perform cross-validation
cv <- krige.cv(log(zinc) ~ 1, locations = meuse, model = anisotropic_fitted)

# Extract observed and predicted values
observed <- cv$observed
predicted <- cv$var1.pred
```


``` r
# Calculate R-squared
r2 <- 1 - (sum((observed - predicted)^2) / sum((observed - mean(observed))^2))

# Display R-squared
cat("R-squared:", r2, "\n")
```

```
## R-squared: 0.7154926
```


``` r
# Calculate RMSE
rmse <- sqrt(mean((observed - predicted)^2))
cat("RMSE:", rmse, "\n")
```

```
## RMSE: 0.3838017
```


``` r
# Calculate MAE
mae <- mean(abs(observed - predicted))
cat("MAE:", mae, "\n")
```

```
## MAE: 0.2933503
```


``` r
# Extract coordinates
kpred_coords <- cbind(kpred, st_coordinates(kpred))

# Check the structure of the new data
head(kpred_coords)
```

```
## Simple feature collection with 6 features and 4 fields
## Geometry type: POINT
## Dimension:     XY
## Bounding box:  xmin: 181100 ymin: 333660 xmax: 181220 ymax: 333740
## Projected CRS: Amersfoort / RD New
##   var1.pred  var1.var      X      Y              geometry
## 1  6.575652 0.4658882 181180 333740 POINT (181180 333740)
## 2  6.680143 0.3807574 181140 333700 POINT (181140 333700)
## 3  6.586995 0.4030249 181180 333700 POINT (181180 333700)
## 4  6.477784 0.4315752 181220 333700 POINT (181220 333700)
## 5  6.790864 0.2880628 181100 333660 POINT (181100 333660)
## 6  6.693545 0.3149947 181140 333660 POINT (181140 333660)
```


``` r
custom_levels <- c(4, 6, 7, 7.5)
custom_labels <- c("2 lbs/ac", "4 lbs/ac", "6 lbs/ac")
ggplot() + geom_sf(data = kpred, aes(color = var1.pred),
                   show.legend = FALSE) +
  geom_sf(data = meuse) +
  geom_contour_filled(data = kpred, aes(x = kpred_coords$X, y = kpred_coords$Y, z = var1.pred), #color = "gray100",
              # bins = 4
              breaks = custom_levels
               ) +
  scale_fill_viridis_d(labels = custom_labels) +  
  scale_fill_manual(values = c("green4", "red4", "yellow"))+
  labs(x = "Longitude", y = "Latitude", fill = "Zinc Prescription (lbs/ac)")+
  theme_bw()
```

```
## Scale for fill is already present.
## Adding another scale for fill, which will replace the existing scale.
```

<img src="{{< blogdown/postref >}}index_files/figure-html/unnamed-chunk-11-1.png" width="672" />


``` r
par(mar = c(0, 1, 0, 1))
pie(
  c(280, 60, 20),
  c('Sky', 'Sunny side of pyramid', 'Shady side of pyramid'),
  col = c('#0292D8', '#F7EA39', '#C4B632'),
  init.angle = -50, border = NA
)
```

<img src="{{< blogdown/postref >}}index_files/figure-html/unnamed-chunk-12-1.png" width="672" />

