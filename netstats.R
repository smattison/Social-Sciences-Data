if(!(dir.exists("r-out"))) dir.create("r-out")

outputfile <- "./r-out/netstats.txt"
sink(outputfile, type = "output", split = TRUE)
sessionInfo()

library(igraph)
source("mosuo-functions.R")

## Names:
## ED density, aka edge density
## DC degree centralization
## BC betweenness centralization
## CC closeness centralization
## TR graph transitivity
## MD average shortest path length, aka mean distance

locs <- c("lb", "yn")

rownumber <- 1
for(loc in locs) {
    g <- make_mosuo_networks(loc)
    row <- unlist(summarize_network(g))
    row <- c(row, "loc" = loc)

    if(loc == locs[[1]]) {
        results <- matrix(ncol = length(row), nrow = length(locs))
        colnames(results) <- names(row)
    }

    results[rownumber, ] <- row
    rownumber <- rownumber + 1
}


results <- as.data.frame(results)

results

sink()
