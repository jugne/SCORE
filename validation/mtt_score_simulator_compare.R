library(ggplot2)
library(gridExtra)
library(ggthemes)

doCompare <- function(dsFileName, ds2FileName) {
    
    pdf(paste(strsplit(ds2FileName, "[.]")[[1]], ".pdf", sep=""), onefile = TRUE)
    
    dfs <- read.table(dsFileName, header=T)
    df <- read.table(ds2FileName, header=T)
    
    dts <- data.frame(dfs)
    dt <- data.frame(df)
    
    # Remove 10% for burnin
    N <- dim(df)[1]
    df <- df[-(1:ceiling(0.1*N)),]
    
    textSize = 10
    
    # Tree age/length plot
    
    p1 <- ggplot() + geom_density(data=dts, aes(x=network.height, colour="SCORE simulator"), size=1)
    p1 <- p1 + geom_density(data=dt, aes(x=tree.height, colour="MTT"), size=1)
    p1 <- p1 + theme_minimal(base_size = textSize) + labs(colour = "", linetype="", x="Height", y="Density")
    p1 <- p1 + theme_hc() + theme(legend.position = c(0.8,0.75))
    
    H <- c()
    L <- c()
    R <- c()
    l <- length(df$tree.height)
    # For comparison to run faster, calculations need to be done at fewer ponts.
    # It can be achieved by taking a fraction smaller than 0.001 below.
    nPoints <- round(l*0.01)
    by <- as.integer(l/nPoints)
    print(l)
    
    for (i in seq(l, 0, by=-by)){
        H <- append(H, log(ks.test(df$tree.height[1:(l-i)], dfs$network.height)$statistic))
        L <- append(L, log(ks.test(df$tree.length[1:(l-i)], dfs$network.totalLength)$statistic))
    }
    
    p3 <- ggplot() + geom_line(aes(x=seq(1, length(H)), y=H, colour="Tree height"), data.frame(H))
    p3 <- p3 + geom_line(aes(x=seq(1, length(L)), y=L, colour="Tree length"), data.frame(L), linetype="dotted")
    p3 <- p3 + theme_minimal(base_size = textSize) + labs(colour = "", linetype="", x=paste("Number of logged iterations (minus burnin) / ", by), y="Kolmogorov-Smirnov stat. (log scale)")
    
    
    grid.arrange(p1, p3)
    
    dev.off()
}

# Set the directory to the directory of the file (validation folder)
this.dir <- dirname(rstudioapi::getSourceEditorContext()$path)
setwd(this.dir)

doCompare("validation/test_xmls_for_MTT/simulate_noReassort_Score.log",
          "validation/test_xmls_for_MTT/simulate_noReassort_MTT.log")



