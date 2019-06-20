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
    p1 <- p1 + theme_hc() + theme(legend.position = c(0.8,0.75)) + scale_color_calc()

    
    grid.arrange(p1)
    
    dev.off()
}

doCompare("/home/ugne/_18_Mokslai/SCORE/validation/simulator/simulate_noReassortment.log",
          "/home/ugne/_18_Mokslai/SCORE/validation/test_xmls_for_MTT/simulate_noReassort_MTT.log")



