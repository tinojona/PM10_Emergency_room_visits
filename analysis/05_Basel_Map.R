################################################################################
### Create a Map of Basel
################################################################################


library(terra);library(sf);library(tmap);library(osmdata);library(tmaptools)
library(easyclimate);library(ggplot2);library(tidyr);library(dplyr);library(osmdata)

rm(list= ls())


# Query Basel city administrative boundaries
basel_boundary <- opq("Basel, Switzerland") %>%
  add_osm_feature(key = "boundary", value = "administrative") %>%
  add_osm_feature(key = "name", value = "Basel") %>%
  osmdata_sf()

# Extract the polygons (check which is relevant)
boundary_polygons <- basel_boundary$osm_multipolygons


# Geocode addresses
addr1 <- geocode_OSM("Petersgraben, Basel, Switzerland")
addr2 <- geocode_OSM("Venusstrasse, Basel, Switzerland")

# Convert to sf points
points_sf <- st_as_sf(
  data.frame(
    name = c("Universitätsspital", "NABEL station"),
    lon = c(addr1$coords["x"], addr2$coords["x"]),
    lat = c(addr1$coords["y"], addr2$coords["y"])
  ),
  coords = c("lon", "lat"),
  crs = 4326
)


# convert coordinate reference system
boundary_polygons <- st_transform(boundary_polygons, 2056)
points_sf <- st_transform(points_sf, 2056)

# Get Rhine river data in Basel
rhein_data <- opq("Basel, Switzerland") %>%
  add_osm_feature(key = "waterway", value = "river") %>%
  add_osm_feature(key = "name", value = "Rhein") %>%
  osmdata_sf()
rhein_lines <- rhein_data$osm_lines


# save plot
png("plots/X_BASEL.png", width = 2500, height = 1700, res = 300)
par(mfrow = c(1,1),
    mar = c(3,3,.5,.5),
    mgp = c(1.8, .5, 0))

ggplot() +
  geom_sf(data = boundary_polygons, fill = "lightblue", color = "black", alpha = 0.2) +
  geom_sf(data = rhein_lines, color = "blue", size = 1) +
  geom_sf(data = points_sf, aes(color = name), size = 3) +
  coord_sf(xlim = c(basel_bbox["xmin"], basel_bbox["xmax"]),
           ylim = c(basel_bbox["ymin"], basel_bbox["ymax"]),
           expand = FALSE) +
  theme_minimal() +
  scale_color_discrete(name = "Location") +
  labs(title = "")

dev.off()





