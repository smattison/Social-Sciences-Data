if(!(dir.exists("r-out"))) dir.create("r-out")
if(!(dir.exists("img"))) dir.create("img")

outputfile <- "./r-out/analysis-pairwise.txt"
sink(outputfile, type = "output", split = TRUE)
sessionInfo()

library(igraph)
library(extrafont)
source("mosuo-functions.R")
source("palettes.R")

locs <- c("lb", "yn")
lb <- make_mosuo_networks("lb")
yn <- make_mosuo_networks("yn")

egosonly <- TRUE

lb_egos <- read.csv("./Data/LB_Ego2.csv", stringsAsFactors = FALSE)
yn_egos <- read.csv("./Data/YN_Ego2.csv", stringsAsFactors = FALSE)
selectcolumns <- c("SharingUnitID", "Ego1", "Ego1M", "Ego2", "Ego2M", "Ego2certainty")
lb_egos <- lb_egos[, selectcolumns]
yn_egos <- yn_egos[, selectcolumns]
lb_egos$loc <- "lb"
yn_egos$loc <- "yn"
egos <- rbind(lb_egos, yn_egos)

dc <- list(## degree
    lb_all = degree(lb, mode = "all"),
    lb_m = degree(lb, v = V(lb)[which(V(lb)$Gender == "M")], mode = "all"),
    lb_f = degree(lb, v = V(lb)[which(V(lb)$Gender == "F")], mode = "all"),
    yn_all = degree(yn, mode = "all"),
    yn_m = degree(yn, v = V(yn)[which(V(yn)$Gender == "M")], mode = "all"),
    yn_f = degree(yn, v = V(yn)[which(V(yn)$Gender == "F")], mode = "all")
)

bcn <- list(## normalized
    lb_all = betweenness(lb, directed = FALSE, normalized = TRUE),
    lb_m = betweenness(lb, v = V(lb)[which(V(lb)$Gender == "M")], directed = FALSE, normalized = TRUE),
    lb_f = betweenness(lb, v = V(lb)[which(V(lb)$Gender == "F")], directed = FALSE, normalized = TRUE),
    yn_all = betweenness(yn, directed = FALSE, normalized = TRUE),
    yn_m = betweenness(yn, v = V(yn)[which(V(yn)$Gender == "M")], directed = FALSE, normalized = TRUE),
    yn_f = betweenness(yn, v = V(yn)[which(V(yn)$Gender == "F")], directed = FALSE, normalized = TRUE)
)

lb_gcc <- get_gcc(lb)
yn_gcc <- get_gcc(yn)
ccn <- list(## normalized
    lb_all = closeness(lb_gcc, mode = "all", normalized = TRUE),
    lb_m = closeness(lb_gcc, v = V(lb_gcc)[which(V(lb_gcc)$Gender == "M")], mode = "all", normalized = TRUE),
    lb_f = closeness(lb_gcc, v = V(lb_gcc)[which(V(lb_gcc)$Gender == "F")], mode = "all", normalized = TRUE),
    yn_all = closeness(yn_gcc, mode = "all", normalized = TRUE),
    yn_m = closeness(yn_gcc, v = V(yn_gcc)[which(V(yn_gcc)$Gender == "M")], mode = "all", normalized = TRUE),
    yn_f = closeness(yn_gcc, v = V(yn_gcc)[which(V(yn_gcc)$Gender == "F")], mode = "all", normalized = TRUE)
)

tr <- list(## local (as opposed to network level) transitivity
    lb_all = transitivity(lb, type = "localundirected"),
    lb_m = transitivity(lb, vids = V(lb)[which(V(lb)$Gender == "M")], type = "localundirected"),
    lb_f = transitivity(lb, vids = V(lb)[which(V(lb)$Gender == "F")], type = "localundirected"),
    yn_all = transitivity(yn, type = "localundirected"),
    yn_m = transitivity(yn, vids = V(yn)[which(V(yn)$Gender == "M")], type = "localundirected"),
    yn_f = transitivity(yn, vids = V(yn)[which(V(yn)$Gender == "F")], type = "localundirected")
)

## distributions including all nodes
pdf(paste0("./img/degree-distributions-allnodes.pdf"),
    width = 7, height = 7, family = "Times New Roman")
four_panel_kde(dc, "Degree", .8)
dev.off()

pdf(paste0("./img/betweenness-distributions-allnodes.pdf"),
    width = 7, height = 7, family = "Times New Roman")
four_panel_kde(bcn, "Betweenness", .8)
dev.off()

pdf(paste0("./img/closeness-distributions-allnodes.pdf"),
    width = 7, height = 7, family = "Times New Roman")
four_panel_kde(ccn, "Closeness", .8)
dev.off()

pdf(paste0("./img/transitivity-distributions-allnodes.pdf"),
    width = 7, height = 7, family = "Times New Roman")
four_panel_kde(tr, "Transitivity", .8)
dev.off()

### Permutation Tests
print("First are the permutations for the full data set, for SI.")
### H_0: μ_X = μ_Y
### H_a: μ_X > μ_Y
### Output is: mean(X), mean(Y), p
### Order matters, with the larger mean first

## Degree
permutation_test(dc$yn_f, dc$yn_m, plot = FALSE)
permutation_test(dc$lb_m, dc$lb_f, plot = FALSE)
permutation_test(dc$lb_m, dc$yn_m, plot = FALSE)
permutation_test(dc$yn_f, dc$lb_f, plot = FALSE)

## Betweenness
permutation_test(bcn$yn_m, bcn$yn_f, plot = FALSE)
permutation_test(bcn$lb_f, bcn$lb_m, plot = FALSE)
permutation_test(bcn$lb_m, bcn$yn_m, plot = FALSE)
permutation_test(bcn$lb_f, bcn$yn_f, plot = FALSE)

## Closeness
permutation_test(ccn$yn_f, ccn$yn_m, plot = FALSE)
permutation_test(ccn$lb_f, ccn$lb_m, plot = FALSE)
permutation_test(ccn$yn_m, ccn$lb_m, plot = FALSE)
permutation_test(ccn$yn_f, ccn$lb_f, plot = FALSE)

## Transitivity
permutation_test(tr$yn_f, tr$yn_m, plot = FALSE)
permutation_test(tr$lb_f, tr$lb_m, plot = FALSE)
permutation_test(tr$lb_m, tr$yn_m, plot = FALSE)
permutation_test(tr$yn_f, tr$lb_f, plot = FALSE)



## This generates warnings about calculating p-values with exact ties for both Wilcoxon and K-S tests.
## Many nodes have the same degree, betweenness, etc.
## Results may be improved by adjusting the options on the two tests.

if(egosonly) {
    keepthese <- c(egos$Ego1, na.omit(egos$Ego2))
    keepidx <- lapply(dc, function(x) which(names(x) %in% keepthese))
    for(sl in c("lb_all", "lb_m", "lb_f", "yn_all", "yn_m", "yn_f")) {
        dc[[sl]] <- dc[[sl]][keepidx[[sl]]]
        bcn[[sl]] <- bcn[[sl]][keepidx[[sl]]]
        ccn[[sl]] <- ccn[[sl]][keepidx[[sl]]]
        tr[[sl]] <- tr[[sl]][keepidx[[sl]]]
    }
}

## localprint <- function(w) {
##     print(lapply(w, function(x) unlist(lapply(x, function(y) y$p.value))))
## }

## print("Degree")
## compare_dc <- suppressWarnings(comparisons(dc))
## localprint(compare_dc)
## print("Betweenness")
## compare_bc <- suppressWarnings(comparisons(bcn))
## localprint(compare_bc)
## print("Closeness")
## compare_cc <- suppressWarnings(comparisons(ccn))
## localprint(compare_cc)
## print("Transitivity")
## compare_tr <- suppressWarnings(comparisons(tr))
## localprint(compare_tr)

## ## Influential data point analysis
## egolist <- c(egos$Ego1, na.omit(egos$Ego2))
## rID <- sample(egolist, 1)#"YN2017"
## dcr <- lapply(dc, function(x) x[which(names(x) != rID)])
## compare_dcr <- suppressWarnings(comparisons(dcr))
## print(paste0("Removing ", rID, ", which has degree ", c(dc$yn_all, dc$lb_all)[rID]))
## print("Original LB F vs YN F WRS p-value:")
## print(compare_dc$medians$flbyn$p.value)
## print(paste0("With ", rID, " removed:"))
## print(compare_dcr$medians$flbyn$p.value)

pdf(paste0("./img/degree-distributions.pdf"),
    width = 7, height = 7, family = "Times New Roman")
four_panel_kde(dc, "Degree", .8)
dev.off()

pdf(paste0("./img/betweenness-distributions.pdf"),
    width = 7, height = 7, family = "Times New Roman")
four_panel_kde(bcn, "Betweenness", .8)
dev.off()

pdf(paste0("./img/closeness-distributions.pdf"),
    width = 7, height = 7, family = "Times New Roman")
four_panel_kde(ccn, "Closeness", .8)
dev.off()

pdf(paste0("./img/transitivity-distributions.pdf"),
    width = 7, height = 7, family = "Times New Roman")
four_panel_kde(tr, "Transitivity", .8)
dev.off()

### Permutation Tests
print("Second are the permutation tests for the respondents only, for the main text.")
### H_0: μ_X = μ_Y
### H_a: μ_X > μ_Y
### Output is: mean(X), mean(Y), p
### Order matters, with the larger mean first

## Degree
permutation_test(dc$yn_f, dc$yn_m, plot = FALSE)
permutation_test(dc$lb_m, dc$lb_f, plot = FALSE)
permutation_test(dc$lb_m, dc$yn_m, plot = FALSE)
permutation_test(dc$yn_f, dc$lb_f, plot = FALSE)

## Betweenness
permutation_test(bcn$yn_m, bcn$yn_f, plot = FALSE)
permutation_test(bcn$lb_f, bcn$lb_m, plot = FALSE)
permutation_test(bcn$lb_m, bcn$yn_m, plot = FALSE)
permutation_test(bcn$lb_f, bcn$yn_f, plot = FALSE)

## Closeness
permutation_test(ccn$yn_f, ccn$yn_m, plot = FALSE)
permutation_test(ccn$lb_f, ccn$lb_m, plot = FALSE)
permutation_test(ccn$yn_m, ccn$lb_m, plot = FALSE)
permutation_test(ccn$yn_f, ccn$lb_f, plot = FALSE)

## Transitivity
permutation_test(tr$yn_f, tr$yn_m, plot = FALSE)
permutation_test(tr$lb_f, tr$lb_m, plot = FALSE)
permutation_test(tr$lb_m, tr$yn_m, plot = FALSE)
permutation_test(tr$yn_f, tr$lb_f, plot = FALSE)


sink()
