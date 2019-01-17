library(ggplot2)
library(gridExtra)

doCompare <- function(dsFileName, mcmcFileName) {
    
    pdf(paste(strsplit(mcmcFileName, "[.]")[[1]], ".pdf", sep=""), onefile = TRUE)
    
    dfs <- read.table(dsFileName, header=T)
    df <- read.table(mcmcFileName, header=T)
    
    dts <- data.frame(dfs)
    dt <- data.frame(df)
    
    # Remove 10% for burnin
    N <- dim(df)[1]
    df <- df[-(1:ceiling(0.1*N)),]
    
    textSize = 10
    
    # Tree age/length plot
    
    # theme_set(theme_minimal(base_size = 10)) 
    
    p1 <- ggplot() + geom_density(data=dt, aes(x=network.height, colour="MCMC", linetype="Network height")) 
    p1 <- p1 + geom_density(data=dt, aes(x=network.totalLength, colour="MCMC", linetype="Network length")) 
    p1 <- p1 + geom_density(data=dts, aes(x=network.height, colour="Direct simulation", linetype="Network height")) 
    p1 <- p1 + geom_density(data=dts, aes(x=network.totalLength, colour="Direct simulation", linetype="Network length"))
    p1 <- p1 + theme_minimal(base_size = textSize) + labs(colour = "", linetype="", x="Statistic", y="Density")
    
    
    # Network node count plot
    
    maxCount <- max(df$network.reassortmentNodeCount, dfs$network.reassortmentNodeCount)+1
    h <- hist(df$network.reassortmentNodeCount, plot=F,
              breaks=seq(-0.5,maxCount+0.5,by=1))
    hs <- hist(dfs$network.reassortmentNodeCount, plot=F,
               breaks=seq(-0.5,maxCount+0.5,by=1))
    
    dt_h <- data.frame(h$mids, h$density)
    dt_hs <- data.frame(hs$mids, hs$density)
    
    colnames(dt_h) <- c("mids", "density")
    colnames(dt_hs) <- c("mids", "density")
    
    p2 <- ggplot() + geom_point(aes(x=mids, y=density, colour="MCMC"), dt_h) + geom_line(aes(x=mids, y=density, colour="MCMC"), dt_h)
    p2 <- p2 + geom_point(aes(x=mids, y=density, colour="Direct simulation"), dt_hs) + geom_line(aes(x=mids, y=density, colour="Direct simulation"), dt_hs)
    p2 <- p2 + theme_minimal(base_size = textSize) + labs(colour = "", linetype="", x="Reassortment count", y="Posterior probability")
    
    
    # Kolmogorov-Smirnov test for network statistics
    
    # H <- c()
    # L <- c()
    # R <- c()
    # l <- length(df$network.height)
    # # For comparison to run faster, calculations need to be done at fewer ponts.
    # # It can be achieved by taking a fraction smaller than 0.001 below.
    # nPoints <- round(l*0.001)
    # by <- as.integer(l/nPoints)
    # 
    # for (i in seq(l, 0, by=-by)){
    #     H <- append(H, log(ks.test(df$network.height[1:(l-i)], dfs$network.height)$statistic))
    #     L <- append(L, log(ks.test(df$network.totalLength[1:(l-i)], dfs$network.totalLength)$statistic))
    #     R <- append(L, log(ks.test(df$network.reassortmentNodeCount[1:(l-i)], dfs$network.reassortmentNodeCount)$statistic))
    # }
    # 
    # p3 <- ggplot() + geom_line(aes(x=seq(1, length(H)), y=H, colour="Network height"), data.frame(H))
    # p3 <- p3 + geom_line(aes(x=seq(1, length(L)), y=L, colour="Network length"), data.frame(L), linetype="dotted")
    # p3 <- p3 + geom_line(aes(x=seq(1, length(R)), y=R, colour="Rassortment count"), data.frame(R), linetype = "dashed")
    # p3 <- p3 + theme_minimal(base_size = textSize) + labs(colour = "", linetype="", x=paste("Number of logged iterations (minus burnin) / ", by), y="Kolmogorov-Smirnov stat. (log scale)")
    # 
    # 
    # grid.arrange(p1, p2, p3)
    
    grid.arrange(p1, p2)
    
    dev.off()
}

setwd("/home/ugne/_18_Mokslai/MascotExtended/validation/simulator")

doCompare("/home/ugne/_18_Mokslai/CoalRe/validation/simulator/simulate_serial5taxon8seg.log",
          "/home/ugne/_18_Mokslai/MascotExtended/validation/simulator/simulate_notStructured.log")
