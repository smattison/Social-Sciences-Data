library(igraph)
source("mosuo-functions.R")

lb <- make_mosuo_networks("lb")
yn <- make_mosuo_networks("yn")

## How many individuals could there be?
yn_indiv <- read.csv("./Data/YN_Indiv.csv")
nrow(subset(yn_indiv, Age >= 15))
lb_indiv <- read.csv("./Data/LB_Indiv.csv")
nrow(subset(lb_indiv, Age >= 15))

## How many cross-gender ties are there?
lbedf <- as.data.frame(ends(lb, E(lb)))
colnames(lbedf) <- c("i", "j")
ynedf <- as.data.frame(ends(yn, E(yn)))
colnames(ynedf) <- c("i", "j")

lbedf$igender <- sapply(lbedf$i, function(i) V(lb)[i]$Gender)
lbedf$jgender <- sapply(lbedf$j, function(j) V(lb)[j]$Gender)
ynedf$igender <- sapply(ynedf$i, function(i) V(yn)[i]$Gender)
ynedf$jgender <- sapply(ynedf$j, function(j) V(yn)[j]$Gender)

print("LB")
print("Number of cross gender edges")
print(nrow(subset(lbedf, igender != jgender)))
print("%")
print(nrow(subset(lbedf, igender != jgender))/ecount(lb))
print("YN")
print("#")
nrow(subset(ynedf, igender != jgender))
print("%")
nrow(subset(ynedf, igender != jgender))/ecount(yn)

## Which Ego1 and Ego2 have no friendship edges?
lb_fullel <- read.csv("./Data/LB_Edgelist.csv", stringsAsFactors = FALSE)
yn_fullel <- read.csv("./Data/YN_Edgelist.csv", stringsAsFactors = FALSE)

lb_egos <- read.csv("./Data/LB_Ego2.csv", stringsAsFactors = FALSE)
yn_egos <- read.csv("./Data/YN_Ego2.csv", stringsAsFactors = FALSE)
selectcolumns <- c("SharingUnitID", "Ego1", "Ego1M", "Ego2", "Ego2M", "Ego2certainty")
lb_egos <- lb_egos[, selectcolumns]
yn_egos <- yn_egos[, selectcolumns]
lb_egos$loc <- "lb"
yn_egos$loc <- "yn"
egos <- rbind(lb_egos, yn_egos)

yn_ego1_deg0 <- yn_egos$Ego1[which(!(yn_egos$Ego1 %in% V(yn)$name))]
yn_ego2_deg0 <- na.omit(yn_egos$Ego2)[which(!(na.omit(yn_egos$Ego2) %in% V(yn)$name))]
lb_ego1_deg0 <- lb_egos$Ego1[which(!(lb_egos$Ego1 %in% V(lb)$name))]
lb_ego2_deg0 <- na.omit(lb_egos$Ego2)[which(!(na.omit(lb_egos$Ego2) %in% V(lb)$name))]

## How many individuals of each gender have non-zero transitivity?
lbtr_m <- transitivity(lb,
                       type = "localundirected",
                       vids = V(lb)[
                           which(V(lb)$Gender == "M" #&
                                      #V(lb)$name %in% c(lb_egos$Ego1, lb_egos$Ego2)
                                 )
                       ])
lbtr_f <- transitivity(lb,
                       type = "localundirected",
                       vids = V(lb)[
                           which(V(lb)$Gender == "F" #&
                                      #V(lb)$name %in% c(lb_egos$Ego1, lb_egos$Ego2)
                                 )
                       ])
yntr_m <- transitivity(yn,
                       type = "localundirected",
                       vids = V(yn)[
                           which(V(yn)$Gender == "M" #&
                                      #V(yn)$name %in% c(yn_egos$Ego1, yn_egos$Ego2)
                                 )
                       ])
yntr_f <- transitivity(yn,
                       type = "localundirected",
                       vids = V(yn)[
                           which(V(yn)$Gender == "F" #&
                                      #V(yn)$name %in% c(yn_egos$Ego1, yn_egos$Ego2)
                                 )
                       ])

print("YN")
print("M")
print("#")
print(sum(yntr_m > 0, na.rm = TRUE))
print("%")
print(sum(yntr_m > 0, na.rm = TRUE)/length(yntr_m))
print("F")
print("#")
print(sum(yntr_f > 0, na.rm = TRUE))
print("%")
print(sum(yntr_f > 0, na.rm = TRUE)/length(yntr_f))

print("LB")
print("M")
print("#")
print(sum(lbtr_m > 0, na.rm = TRUE))
print("%")
print(sum(lbtr_m > 0, na.rm = TRUE)/length(lbtr_m))
print("F")
print("#")
print(sum(lbtr_f > 0, na.rm = TRUE))
print("%")
print(sum(lbtr_f > 0, na.rm = TRUE)/length(lbtr_f))

print("For all")
print("YN")
print(sum(yntr_m > 0, na.rm = TRUE) + sum(yntr_f > 0, na.rm = TRUE)) # non-zero transitivity
print(length(yntr_m) + length(yntr_f))
print("LB")
print(sum(lbtr_m > 0, na.rm = TRUE) + sum(lbtr_f > 0, na.rm = TRUE))
print(length(lbtr_m) + length(lbtr_f))

modes(
    transitivity(
        lb, "localundirected", isolates = "zero",
        vids = V(lb)[
            which(V(lb)$Gender == "M" & V(lb)$name %in% c(lb_egos$Ego1, lb_egos$Ego2))
        ]
    ))
modes(
    transitivity(
        lb, "localundirected", isolates = "zero",
        vids = V(lb)[
            which(V(lb)$Gender == "F" & V(lb)$name %in% c(lb_egos$Ego1, lb_egos$Ego2))
        ]
    ))
modes(
    transitivity(
        yn, "localundirected", isolates = "zero",
        vids = V(yn)[
            which(V(yn)$Gender == "M" & V(yn)$name %in% c(yn_egos$Ego1, yn_egos$Ego2))
        ]
    ))
modes(
    transitivity(
        yn, "localundirected", isolates = "zero",
        vids = V(yn)[
            which(V(yn)$Gender == "F" & V(yn)$name %in% c(yn_egos$Ego1, yn_egos$Ego2))
        ]
    ))
