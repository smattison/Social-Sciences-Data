library(igraph)
source("mosuo-functions.R")

locs <- c("lb", "yn")
ext <- "gml"
subdir <- "./Data/"

for(loc in locs) {
    exportfilename <- paste0(subdir, loc, "-friendship.", ext)

    g <- make_mosuo_networks(loc)

    write_graph(g, exportfilename, format = "gml")
        
}

