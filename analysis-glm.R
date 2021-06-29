if(!(dir.exists("r-out"))) dir.create("r-out")
if(!(dir.exists("img"))) dir.create("img")

outputfile <- "./r-out/analysis-glm.txt"
sink(outputfile, type = "output", split = TRUE)
sessionInfo()

library(igraph)
library(extrafont)
source("mosuo-functions.R")
source("palettes.R")
## requires boot, lme4, and lmtest

lb <- make_mosuo_networks("lb")
yn<- make_mosuo_networks("yn")

egosonly <- FALSE

lb_indiv <- subset(read.csv("./Data/LB_Indiv.csv", stringsAsFactors = FALSE), select = -X)
lb_indiv$loc <- "lb"
yn_indiv <- read.csv("./Data/YN_Indiv.csv", stringsAsFactors = FALSE)
colnames(yn_indiv)[which(colnames(yn_indiv) == "YearsofEdu")] <- "YearsOfEdu"
yn_indiv$YearsOfEdu <- suppressWarnings(as.numeric(yn_indiv$YearsOfEdu)) # because step above made it a character vector; this step also produces the warning message, "NAs introduced by coercion"
yn_indiv$loc <- "yn"

dc <- list(
    lb = degree(lb),
    yn = degree(yn)
)

lb_indiv$dc <- sapply(lb_indiv$IndivID, function(x) if(x %in% names(dc$lb)) dc$lb[[x]] else NA)
yn_indiv$dc <- sapply(yn_indiv$IndivID, function(x) if(x %in% names(dc$yn)) dc$yn[[x]] else NA)

indiv <- rbind(lb_indiv, yn_indiv)

lb_egos <- read.csv("./Data/LB_Ego2.csv", stringsAsFactors = FALSE)
yn_egos <- read.csv("./Data/YN_Ego2.csv", stringsAsFactors = FALSE)
selectcolumns <- c("SharingUnitID", "Ego1", "Ego1M", "Ego2", "Ego2M", "Ego2certainty")
lb_egos <- lb_egos[, selectcolumns]
yn_egos <- yn_egos[, selectcolumns]
lb_egos$loc <- "lb"
yn_egos$loc <- "yn"
egos <- rbind(lb_egos, yn_egos)

indiv$ego <- ifelse(indiv$IndivID %in% c(egos$Ego1, egos$Ego2), TRUE, FALSE)
indiv$ego1 <- ifelse(indiv$IndivID %in% egos$Ego1, TRUE, FALSE)
indiv$ego2 <- ifelse(indiv$IndivID %in% egos$Ego2, TRUE, FALSE)
if(egosonly) {
    indiv <- subset(indiv, IndivID %in% c(egos$Ego1, egos$Ego2))
} else {
    indiv <- subset(indiv, indiv$IndivID %in% c(V(lb)$name, V(yn)$name))
}

data <- indiv
data$loc <- factor(data$loc)
data$Gender <- factor(data$Gender)

## Nested models, decision criteria is BIC
m0 <- glm(dc ~ ego2 + Gender, data = data, family = poisson(link = "log"))
m1 <- glm(dc ~ ego + ego2 + Age + Gender + loc, data = data, family = poisson(link = "log"))
m2 <- glm(dc ~ ego + ego2 + Age + YearsOfEdu + Gender + loc, data = data, family = poisson(link = "log")) # YearsOfEdu is confounded with location (higher in patriliny) and although it lowers the AIC it itself has no effect on the outcome and weakens the association with gender
BIC(m0)
BIC(m1)
BIC(m2)
m3 <- glm(dc ~ ego + ego2 + Age + Gender * loc, data = data, family = poisson(link = "log"))

## choose m1
print("Chosen model")
print(summary(m1))
confint(m1)

print("Alt model using interaction term")
print(summary(m3))
confint(m3)

m3u <- update(
    m3,
    data = within(data, {
        loc <- relevel(loc, ref = "yn")
        Gender <- relevel(Gender, ref = "F")
    }))
summary(m3u)
confint(m3u)
m3u2 <- update(
    m3,
    data = within(data, {
        loc <- relevel(loc, ref = "lb")
        Gender <- relevel(Gender, ref = "F")
    }))
summary(m3u2)
confint(m3u2)


diagnostics <- boot::glm.diag(m1)
diagnostics_filename <- paste0("./r-out/model-diagnostics.pdf")
pdf(diagnostics_filename, width = 7, height = 7)
boot::glm.diag.plots(m1, diagnostics)
dev.off()

## hm3 <- lme4::glmer(dc ~ Age + YearsOfEdu + Gender * loc + (1|SharingUnitID),
##                    data = data, family = poisson(link = "identity"), na.action = na.omit)
## hm3 <- lme4::glmer(dc ~ Gender * loc + (1 | IndivID), family = poisson,
##                    data = data, na.action = na.omit)
## print("Likelihood ratio test")
## print(lmtest::lrtest(m3, hm3))
## print("Hierarchical model")
## print(summary(hm3))

sink()

## Predicted Values
predicted <- exp(predict(m3)) # m1
data$predicted <- sapply(rownames(data),
                         function(x) if(x %in% names(predicted)) predicted[x] else NA)

## Boxplot
fcolor <- inauguration$Obama
mcolor <- inauguration$Biden

boxplot_filename <- paste0("./img/predicted-degree.pdf")
pdf(boxplot_filename, family = "Times New Roman", width = 7, height = 7)

par(bty = "n")
figdat <- subset(data, ego == TRUE)
figdat$Gender <- factor(figdat$Gender, c("F", "M"))
figdat$loc <- factor(figdat$loc, c("yn", "lb"))
ymin <- min(figdat$predicted)*.9
ymax <- max(figdat$predicted)*1.1
boxes <- boxplot(predicted ~ Gender + loc, data = figdat, plot = FALSE)
bxp(boxes, show.names = FALSE, xlab = "", ylab = "Predicted Degree",
    ylim = c(ymin, ymax),
    axes = FALSE, #log = "y",
    border = rep(c(fcolor, mcolor), times = 2),
    lwd = 2)
abline(v = mean(1:4), lty = 1, lwd = .5, col = "black")
text(x = 1.5, y = max(predicted)*1.05, pos = 3, label = "Matrilineal", cex = 1.5)
text(x = 3.5, y = max(predicted)*1.05, pos = 3, label = "Patrilineal", cex = 1.5)
axis(side = 2,
     at = seq(round(min(predicted)*0.95, 1), round(max(predicted)*1.05, 1), by = .2))
legend("bottomleft", bty = "n", col = c(fcolor, mcolor), lwd = 3, lty = 1,
       legend = c("Women", "Men"))

dev.off()

## Without location
boxplot_filename <- paste0("./img/predicted-degree-noloc.pdf")
pdf(boxplot_filename, family = "Times New Roman", width = 7, height = 7)

par(bty = "n")
figdat <- subset(data, ego == TRUE)
ymin <- min(figdat$predicted)*.9
ymax <- max(figdat$predicted)*1.1
boxes <- boxplot(predicted ~ Gender, data = figdat, plot = FALSE)
bxp(boxes, show.names = FALSE, xlab = "", ylab = "Predicted Degree",
    ylim = c(ymin, ymax),
    axes = FALSE, #log = "y",
    border = rep(c(fcolor, mcolor), times = 2),
    lwd = 2)

## abline(v = mean(1:4), lty = 1, lwd = .5, col = "black")
## text(x = 1.5, y = max(predicted)*1.05, pos = 3, label = "Patrilineal", cex = 1.5)
## text(x = 3.5, y = max(predicted)*1.05, pos = 3, label = "Matrilineal", cex = 1.5)
axis(side = 2,
     at = seq(round(min(predicted)*0.95, 1), round(max(predicted)*1.05, 1), by = .2))
legend("bottomleft", bty = "n", col = c(fcolor, mcolor), lwd = 3, lty = 1,
       legend = c("Women", "Men"))

dev.off()


## additional code for additional reporting
table(subset(indiv, select = c(ego1, Gender, loc)))
table(subset(indiv, select = c(ego2, Gender, loc)))
