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
    
    n_cols <- ncol(df1)
    par(mfcol=c((n_cols-1)/3,3))
    
   for (i in 2:n_cols){
     maxCount <- max(df1[i], df2[i])+1
     
     h1 <- hist(df1[[i]], breaks=seq(-0.5,maxCount+0.5,1), plot=F)
     h2 <- hist(df2[[i]], breaks=seq(-0.5,maxCount+0.5,1), plot=F)
     
     maxCount_y = max(max(h1$density), max(h2$density))
     
     plot(h1$mids, h1$density, 'o', col='blue',
          xlab=paste('Count', names(df1[i])[1], sep = " ", collapse = NULL), ylim = c(0, 1), xlim= c(0,6), ylab='Relative Frequency')
     lines(h2$mids, h2$density, 'o', col='red')
     
     legend("topright", inset=0.04,
            c("MCMC", "Direct simulation"), lty=1, pch=1, lwd=2,
            col=c("red","blue"))
   }
    
  dev.off()
}

doCompare("/Users/jugne/Documents/Source/SCORE/validation/mapping/3tax_3seg_3types/high_migration/simulation/simulate_3tax_3seg_3type.typeStats.log",
          "/Users/jugne/Documents/Source/SCORE/validation/mapping/3tax_3seg_3types/high_migration/inference/test_approx_3tax_3seg_3type.stats.log")


doCompare("/Users/jugne/Documents/Source/SCORE/validation/mapping/3tax_3seg_3types/medium_mig/simulation/simulate_3tax_3seg_3type.typeStats.log",
          "/Users/jugne/Documents/Source/SCORE/validation/mapping/3tax_3seg_3types/medium_mig/inference/test_approx_3tax_3seg_3type.stats.log")

doCompare("/Users/jugne/Documents/Source/SCORE/validation/mapping/3tax_3seg_3types/slow_mig/simulation/simulate_3tax_3seg_3type.typeStats.log",
          "/Users/jugne/Documents/Source/SCORE/validation/mapping/3tax_3seg_3types/slow_mig/inference/test_approx_3tax_3seg_3type.stats.log")

