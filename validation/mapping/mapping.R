library(ggplot2)
library(grid)
library(gridExtra)
library(coda)

# clear workspace
rm(list = ls())


doCompare <- function(dsFileName, mcmcFileName) {
    df1 <- read.table(dsFileName, header=T)
    df2 <- read.table(mcmcFileName, header=T)
    
    maxCount <- max(df1$network.count_0_to_1, df2$mappedNet.count_0_to_1)+1
    
    h1 <- hist(df1$network.count_0_to_1, breaks=seq(-0.5,maxCount+0.5,1), plot=F)
    h2 <- hist(df2$mappedNet.count_0_to_1, breaks=seq(-0.5,maxCount+0.5,1), plot=F)
    
    # Tree age/length plot
    pdf(gsub(".log", ".pdf", mcmcFileName), width = 10, height = 10) 
    
    par(mfcol=c(2,1))
    
    maxCount_y = max(max(h1$density), max(h2$density))
    
    plot(h1$mids, h1$density, 'o', col='blue',
         xlab='Count 0 -> 1', ylim = c(0, maxCount_y),, ylab='Relative Frequency')
    lines(h2$mids, h2$density, 'o', col='red')

    maxCount <- max(df1$network.count_1_to_0, df2$mappedNet.count_1_to_0)+1
    
    h1 <- hist(df1$network.count_1_to_0, breaks=seq(-0.5,maxCount+0.5,1), plot=F)
    h2 <- hist(df2$mappedNet.count_1_to_0, breaks=seq(-0.5,maxCount+0.5,1), plot=F)
    
    maxCount_y = max(max(h1$density), max(h2$density))
    
    plot(h1$mids, h1$density, 'o', col='blue',
         xlab='Count 1 -> 0', ylim = c(0, maxCount_y), ylab='Relative Frequency')
    lines(h2$mids, h2$density, 'o', col='red')
  dev.off()
}

doCompare("simulate_2tax_2seg_2type.typeStats.log",
          "test_approx_2tax_2seg_2type.stats.log")

doCompare("simulate_2tax_2seg_2type_diff_rates.typeStats.log",
          "test_approx_2tax_2seg_2type_diff_rates.stats.log")

doCompare("simulate_3tax_3seg_3type.typeStats.log",
          "test_approx_3tax_3seg_3type.stats.log")


