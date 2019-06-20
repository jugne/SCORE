library(ggplot2)
library(gridExtra)
library(ggthemes)

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
    
    
    p1 <- ggplot() + geom_density(data=dts, aes(x=network.height, colour="SCORE simulator", linetype="Network height")) 
    p1 <- p1 + geom_density(data=dts, aes(x=network.totalLength, colour="SCORE simulator", linetype="Network length")) 
    p1 <- p1 + geom_density(data=dt, aes(x=network.height, colour="CoalRe", linetype="Network height")) 
    p1 <- p1 + geom_density(data=dt, aes(x=network.totalLength, colour="CoalRe", linetype="Network length"))
    p1 <- p1 + theme_minimal(base_size = textSize) + labs(colour = "", linetype="", x="Statistic", y="Density")
    p1 <- p1 + theme_hc() + theme(legend.position = c(0.8,0.72)) + scale_color_calc()
    
    # Network node count plot
    
    maxCount <- max(df$network.reassortmentNodeCount, dfs$network.reassortmentNodeCount)+1
    l<-length(df$network.reassortmentNodeCount)
    h <- hist(df$network.reassortmentNodeCount, plot=F,
              breaks=seq(-0.5,maxCount+0.5,by=1))
    hs <- hist(dfs$network.reassortmentNodeCount, plot=F,
               breaks=seq(-0.5,maxCount+0.5,by=1))
    
    dt_h <- data.frame(h$mids, h$density)
    dt_hs <- data.frame(hs$mids, hs$density)
    
    colnames(dt_h) <- c("mids", "density")
    colnames(dt_hs) <- c("mids", "density")
    
    p2 <- ggplot() + geom_point(aes(x=mids, y=density, colour="SCORE simulator"), dt_hs) + geom_line(aes(x=mids, y=density, colour="SCORE simulator"), dt_hs)
    p2 <- p2 + geom_point(aes(x=mids, y=density, colour="CoalRe"), dt_h) + geom_line(aes(x=mids, y=density, colour="CoalRe"), dt_h)
    p2 <- p2 + theme_minimal(base_size = textSize) + labs(colour = "", linetype="", x="Reassortment count", y="Posterior probability")
    p2 <- p2 + theme_hc() + theme(legend.position = c(0.8,0.72)) + scale_color_calc()
    
    grid.arrange(p1, p2)
    
    dev.off()
}

doCompare("/home/ugne/_18_Mokslai/CoalRe/validation/simulator/simulate_serial5taxon8seg.log",
          "/home/ugne/_18_Mokslai/SCORE/validation/simulator/simulate_notStructured.log")



