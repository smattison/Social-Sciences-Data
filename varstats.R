if(!(dir.exists("r-out"))) dir.create("r-out")

outputfile <- "./r-out/varstats.txt"
sink(outputfile, type = "output", split = TRUE)
sessionInfo()

library(igraph)
source("mosuo-functions.R")

locs <- c("lb", "yn")
lb <- make_mosuo_networks("lb")
yn <- make_mosuo_networks("yn")

egosonly <- FALSE # TRUE

## For total gender and age distribution of all individuals discovered by the survey

lb_indiv <- subset(read.csv("./Data/LB_Indiv.csv", stringsAsFactors = FALSE), select = -X)
lb_indiv$loc <- "lb"
yn_indiv <- read.csv("./Data/YN_Indiv.csv", stringsAsFactors = FALSE)
colnames(yn_indiv)[which(colnames(yn_indiv) == "YearsofEdu")] <- "YearsOfEdu"
yn_indiv$YearsOfEdu <- suppressWarnings(as.numeric(yn_indiv$YearsOfEdu)) # because step above made it a character vector; this step also produces the warning message, "NAs introduced by coercion"
yn_indiv$loc <- "yn"
indiv <- rbind(lb_indiv, yn_indiv)
## Only consider those found by friendship edge
indiv <- subset(indiv, IndivID %in% c(V(lb)$name, V(yn)$name))

lb_egos <- read.csv("./Data/LB_Ego2.csv", stringsAsFactors = FALSE)
yn_egos <- read.csv("./Data/YN_Ego2.csv", stringsAsFactors = FALSE)
selectcolumns <- c("SharingUnitID", "Ego1", "Ego1M", "Ego2", "Ego2M", "Ego2certainty")
lb_egos <- lb_egos[, selectcolumns]
yn_egos <- yn_egos[, selectcolumns]
lb_egos$loc <- "lb"
yn_egos$loc <- "yn"
egos <- rbind(lb_egos, yn_egos)

indiv_red <- subset(indiv, IndivID %in% c(egos$Ego1, egos$Ego2))

print("Gender distribution for all individuals discovered by the survey:")
round(prop.table(table(subset(indiv, select = c(Gender, loc)), useNA = "always")), 4) * 100
table(subset(indiv, select = c(Gender, loc)))
print("Gender distribution for Ego1 and Ego2 only:")
round(prop.table(table(subset(indiv_red, select = c(Gender, loc)), useNA = "always")), 4) * 100
table(subset(indiv_red, select = c(Gender, loc)), useNA = "always")

print("Age distribution for all individuals discovered by the survey, by location and gender:")
print("YN")
unlist(summarize_variable(indiv$Age[indiv$loc == "yn" & indiv$Gender == "F"]))
unlist(summarize_variable(indiv$Age[indiv$loc == "yn" & indiv$Gender == "M"]))
print("LB")
unlist(summarize_variable(indiv$Age[indiv$loc == "lb" & indiv$Gender == "F"]))
unlist(summarize_variable(indiv$Age[indiv$loc == "lb" & indiv$Gender == "M"]))

print("Age distribution for Ego1 and Ego2, by location:")
print("YN")
print("women")
unlist(summarize_variable(indiv_red$Age[indiv_red$loc == "yn" & indiv_red$Gender == "F"]))
print("men")
unlist(summarize_variable(indiv_red$Age[indiv_red$loc == "yn" & indiv_red$Gender == "M"]))
print("LB")
print("Women")
unlist(summarize_variable(indiv_red$Age[indiv_red$loc == "lb" & indiv_red$Gender == "F"]))
print("Men")
unlist(summarize_variable(indiv_red$Age[indiv_red$loc == "lb" & indiv_red$Gender == "M"]))


print("Distribution of years of education for all individuals discovered by the survey, by location and gender:")
print("YN")
unlist(summarize_variable(indiv$YearsOfEdu[indiv$loc == "yn" & indiv$Gender == "F"]))
unlist(summarize_variable(indiv$YearsOfEdu[indiv$loc == "yn" & indiv$Gender == "M"]))
print("LB")
unlist(summarize_variable(indiv$YearsOfEdu[indiv$loc == "lb" & indiv$Gender == "F"]))
unlist(summarize_variable(indiv$YearsOfEdu[indiv$loc == "lb" & indiv$Gender == "M"]))

print("Distribution of years of education for Ego1 and Ego2, by location:")
print("YN")
print("women")
unlist(summarize_variable(indiv_red$YearsOfEdu[indiv_red$loc == "yn" & indiv_red$Gender == "F"]))
print("men")
unlist(summarize_variable(indiv_red$YearsOfEdu[indiv_red$loc == "yn" & indiv_red$Gender == "M"]))

print("LB")
print("Women")
unlist(summarize_variable(indiv_red$YearsOfEdu[indiv_red$loc == "lb" & indiv_red$Gender == "F"]))
print("men")
unlist(summarize_variable(indiv_red$YearsOfEdu[indiv_red$loc == "lb" & indiv_red$Gender == "M"]))

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

## for focusing on egos
keepthese <- c(egos$Ego1, na.omit(egos$Ego2))
keepidx <- lapply(dc, function(x) which(names(x) %in% keepthese))

listvars <- c("dc", "bcn", "ccn", "tr")
sublists <- c("yn_f", "yn_m", "lb_f", "lb_m")

print("Node-level statistics for the whole network, removing NA values:")
for(lst in listvars) {
    print(lst)
    for(sublst in sublists) {
        print(sublst)
        print(summary(get(lst)[[sublst]]))
    }
}

print("Node-level statistics for Ego1 and Ego2 only:")
for(lst in listvars) {
    print(lst)
    for(sublst in sublists) {
        print(sublst)
        print(summary(get(lst)[[sublst]][keepidx[[sublst]]]))
    }
}


sink()
                      
