## Provides functions for the "Market Transitions -- Mosuo" network study.

make_mosuo_networks <- function(loc) {
    "Given either 'lb' or 'yn', make the friendship network from that study site."
    require(igraph)

    if(loc == "lb") {
        indiv <- subset(read.csv("./Data/LB_Indiv.csv", stringsAsFactors = FALSE),
                        select = -X)
        el <- read.csv("./Data/LB_Edgelist.csv", stringsAsFactors = FALSE)
        egokey <- read.csv("./Data/LB_Ego2.csv", stringsAsFactors = FALSE)
    } else if(loc == "yn") {
        indiv <- read.csv("./Data/YN_Indiv.csv", stringsAsFactors = FALSE)
        colnames(indiv)[which(colnames(indiv) == "YearsofEdu")] <- "YearsOfEdu"
        indiv$YearsOfEdu <- suppressWarnings(as.numeric(indiv$YearsOfEdu)) # because step above made it a character vector; this step also produces the warning message, "NAs introduced by coercion"
        el <- read.csv("./Data/YN_Edgelist.csv", stringsAsFactors = FALSE)
        egokey <- read.csv("./Data/YN_Ego2.csv", stringsAsFactors = FALSE)
        
        ## one alter is missing from indiv that is present in el---it's YN9076, which should be 2000
        el$Alter[which(el$Alter == "YN9076")] <- "YN2000"
        ## one ego recorded as 2016 should be 2017
        el$Ego[which(el$Ego == "YN2016")] <- "YN2017"
        } else print("Specify a study site (either 'lb' or 'yn').")
    
    ties <- c("WomenHangOutAfterDinner", "MenHangOutAfterDinner")
    el <- el[el$Network %in% unlist(ties), ]

    egos <- c(egokey$Ego1, na.omit(egokey$Ego2))

    for(i in 1:nrow(el)) {
        q <- el[i, "Network"]
        u <- el[i, "Ego"]

        if(!(u %in% egokey$Ego1)) next

        alter <- el[i, "Alter"]
        keyrow <- subset(egokey, Ego1 == u)
        v <- keyrow$Ego2
        stopifnot(!is.na(v) | is.na(v))

        if(is.na(v)) next
        
        if(alter == v | alter == u) next

        if(q == "WomenHangOutAfterDinner") {
            if(keyrow$Ego1M == 1) el[i, "Ego"] <- v
        } else if(q == "MenHangOutAfterDinner") {
            if(keyrow$Ego1M == 0) el[i, "Ego"] <- v
        } else rowvec[[2]] <- NA
    }

    el$EdgeType <- "Friendship"

    el <- subset(el, subset = Ego %in% egos, select = -Network)
    indiv <- subset(indiv, IndivID %in% c(egos, unique(el$Alter)))
    
    g <- graph_from_data_frame(el, directed = FALSE, vertices = indiv)
    g <- simplify(g)

    return(g)
}

get_gcc <- function(g) {
    "Return largest (giant) connected component of a graph."
    require(igraph)
    
    components <- clusters(g, mode = "weak")
    gcc_id <- which.max(components$csize)
    v_ids <- V(g)[components$membership == gcc_id]
    g <- induced_subgraph(g, v_ids)

    return(g)
}

comparisons <- function(dl) {
    "Test for differences across major contrasts in the Mosuo data. Requires a list of data with the following subgroups specified: `lb_m`, `lb_f`, `yn_m`, and `yn_f`. Returns a list with a sublist for tests of medians (Wilcoxon) and a separate sublist for tests of distributions (Kolmogorov-Smirnov)."
    
    ## test for differences in medians and distributions four ways:
    ## lb male vs lb female
    ## yn male vs yn female
    ## lb male vs yn male
    ## lb female vs yn female

    ## tests of medians
    medians <- list(
        lbmf = wilcox.test(dl$lb_m, dl$lb_f),
        ynmf = wilcox.test(dl$yn_m, dl$yn_f),
        mlbyn = wilcox.test(dl$lb_m, dl$yn_m),
        flbyn = wilcox.test(dl$lb_f, dl$yn_f)
    )
    ## tests of distributions
    distributions <- list(
        lbmf = ks.test(dl$lb_m, dl$lb_f),
        ynmf = ks.test(dl$yn_m, dl$yn_f),
        mlbyn = ks.test(dl$lb_m, dl$yn_m),
        flbyn = ks.test(dl$lb_f, dl$yn_f)
    )

    results <- list(
        medians = medians,
        distributions = distributions
    )

    return(results)
}

four_panel_kde <- function(dl, plottitle = "", alpha = 1, bw = "nrd0",
                           mcolor = inauguration$Biden, fcolor = inauguration$Obama) {
    "Generate a plot with a KDE for each contrast."

    layout(matrix(c(1, 2, 3, 4), ncol = 2, byrow = TRUE),
           widths = c(1, 1), heights = c(1, 1), respect = TRUE)
    par(cex.lab = 1.5)

    gendercolors <- rep(c(fcolor, mcolor), times = 2)
    datas <- c("yn_f", "yn_m", "lb_f", "lb_m")

    kde <- list(
        lb_m = density(dl$lb_m[which(is.finite(dl$lb_m))], bw = bw, na.rm = TRUE),
        lb_f = density(dl$lb_f[which(is.finite(dl$lb_f))], bw = bw, na.rm = TRUE),
        yn_m = density(dl$yn_m[which(is.finite(dl$yn_m))], bw = bw, na.rm = TRUE),
        yn_f = density(dl$yn_f[which(is.finite(dl$yn_f))], bw = bw, na.rm = TRUE)
    )

    datmax <- max(unlist(dl)[which(is.finite(unlist(dl)))], na.rm = TRUE)
    kdemax <- max(unlist(lapply(kde, function(x) max(x$x))))
    
    xmin <- 0
    if(datmax <= 1 & kdemax > 1) {
        xmax <- 1
    } else {
        xmax <- kdemax
    }
    ymin <- 0
    ymax <- max(unlist(lapply(kde, function(x) max(x$y))))
    
    for(i in 1:4) {
        plot.new()
        if(i %in% c(1, 3)) par(mar = c(2.1, 3.1, .5, .1)) else par(mar = c(2.1, 2.1, .5, 1.1))
        plot.window(xlim = c(xmin, xmax),
                    ylim = c(ymin, ymax))

        dname <- datas[i]
        theseidx <- which(kde[[dname]][["x"]] >= xmin & kde[[dname]][["x"]] <= xmax)
        x <- kde[[dname]][["x"]][theseidx]
        y <- kde[[dname]][["y"]][theseidx]

        polygon(
            c(x, rev(x)),
            c(y, rep(0, length.out = length(y))),
            col = adjustcolor(gendercolors[i], alpha.f = alpha),
            border = gendercolors[i], lwd = 2)

        axis(1)
        if(i == 1) title(ylab = "Matrilineal", line = 1)
        if(i == 3) title(ylab = "Patrilineal", line = 1)

        if(i == 2) legend("topright", bty = "n", legend = c("Women", "Men"), pch = 21, cex = 1.5,
                          pt.bg = sapply(c(fcolor, mcolor),
                                         function(x) adjustcolor(x, alpha.f = alpha)),
                          col = c(fcolor, mcolor))
    }
}

centr_pagerank <- function(g, algo = "prpack") {
    "Calculate the PageRank centralization of a graph."
    require(igraph)
    
    scores <- page_rank(g, algo = algo)$vector
    sg <- make_star(vcount(tg), "in")
    sgscores <- page_rank(sg)$vector
    tmax <- sum(max(sgscores) - sgscores[which(sgscores != max(sgscores))])
    pc <- centralize(scores, theoretical.max = tmax)
    return(pc)
}

summarize_network <- function(g) {
    "Collect network-level summary statistics for Mosuo networks."
    require(igraph)

    calctypes <- c(# Capitolized to signify graph-level statistics
        "ED", # density, aka edge density
        "DC", # degree centralization
        "BC", # betweenness centralization
        "CC", # closeness centralization
        "TR", # graph transitivity
        "MD"  # average shortest path length, aka mean distance
    )
    netstats <- vector("list", length(calctypes))
    names(netstats) <- calctypes

    netstats[["ED"]] <- edge_density(g)
    netstats[["DC"]] <- centr_degree(g, "all")$centralization
    netstats[["BC"]] <- centr_betw(g)$centralization
    netstats[["CC"]] <- centr_clo(get_gcc(g), "all")$centralization
    netstats[["TR"]] <- transitivity(g)
    netstats[["MD"]] <- mean_distance(g)

    return(netstats)
}

geomean <- function(x, na.rm = TRUE) {
    "Calculate the geometric mean of a variable `x`."
    ## A less robust alternative: prod(x)^(1/length(x))
    
    return(exp(mean(log(x), na.rm = na.rm)))
}

modes <- function(x) {# https://stackoverflow.com/questions/2547402/how-to-find-the-statistical-mode/8189441#8189441
  ux <- unique(x)
  tab <- tabulate(match(x, ux))
  ux[tab == max(tab)]
}

summarize_variable <- function(x, na.rm = TRUE) {
    "Collect variable-level summary statistics."
    
    if(!(is.numeric(x))) return("x is not numeric.")

    varstats <- list(
        m = mean(x, na.rm = na.rm),
        sd = sd(x, na.rm = na.rm),
        xmin = min(x, na.rm = na.rm),
        xmax = max(x, na.rm = na.rm),
        gm = geomean(x, na.rm = na.rm),
        med = median(x, na.rm = na.rm),
        iqr = IQR(x, na.rm = na.rm),
        sumna = sum(is.na(x))
    )

    return(varstats)
}

permutation_test <- function(X, Y, nsims = 10000, plot = TRUE) {
    "Compare the difference in means between two samples. Asks if the mean of X is greater than the mean of Y."
    m1 <- mean(X, na.rm = TRUE)
    m2 <- mean(Y, na.rm = TRUE)
    ev <- m1 - m2

    dat <- c(na.omit(X), na.omit(Y))
    comps <- numeric(nsims)
    for(i in 1:nsims) {
        simdat <- sample(dat)
        s1 <- simdat[1:length(X)]
        s2 <- simdat[(length(X) + 1):length(simdat)]

        sv <- mean(s1) - mean(s2)
        comps[[i]] <- sv
    }

    p <- sum(comps > ev)/nsims

    if(plot) {
        hist(comps)
        abline(v = ev, col = "red", lwd = 3)
    }

    results <- c(m1 = m1, m2 = m2, p = p)
    return(results)
}

