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
    l<-length(df$network.reassortmentNodeCount)
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
    
    H <- c()
    L <- c()
    R <- c()
    l <- length(df$network.height)
    # For comparison to run faster, calculations need to be done at fewer ponts.
    # It can be achieved by taking a fraction smaller than 0.001 below.
    step <- 0.001
    nPoints <- round(l*step)
    by <- as.integer(l/nPoints)
    print(by)

    for (i in seq(l, 0, by=-by)){
        H <- append(H, log(ks.test(df$network.height[1:(l-i)], dfs$network.height)$statistic))
        L <- append(L, log(ks.test(df$network.totalLength[1:(l-i)], dfs$network.totalLength)$statistic))
        R <- append(L, log(ks.test(df$network.reassortmentNodeCount[1:(l-i)], dfs$network.reassortmentNodeCount)$statistic))
    }

    p3 <- ggplot() + geom_line(aes(x=seq(1, length(H)), y=H, colour="Network height"), data.frame(H))
    p3 <- p3 + geom_line(aes(x=seq(1, length(L)), y=L, colour="Network length"), data.frame(L), linetype="dotted")
    p3 <- p3 + geom_line(aes(x=seq(1, length(R)), y=R, colour="Rassortment count"), data.frame(R), linetype = "dashed")
    p3 <- p3 + theme_minimal(base_size = textSize) + labs(colour = "", linetype="", x=paste("Number of logged iterations (minus burnin) / ", (by)), y="Kolmogorov-Smirnov stat. (log scale)")


    grid.arrange(p1, p2, p3)
    
    # grid.arrange(p1, p2)
    
    dev.off()
}

# Set the directory to the directory of the file (validation folder)
this.dir <- dirname(rstudioapi::getSourceEditorContext()$path)
setwd(this.dir)

## Structured + Reassortment

doCompare("exact_test/simulate_3tax_3seg_2type_diff_rates.log",
          "exact_test/exact_inf/test_exact_3tax_3seg_2type_NonCoupledMCMC_diff_rates_resimulate.log")

doCompare("exact_test/simulate_3tax_3seg_2type_diff_rates.log",
          "exact_test/exact_inf/test_exact_3tax_3seg_2type_NonCoupledMCMC_diff_rates.log")

doCompare("exact_test/simulate_3tax_3seg_2type.log",
          "exact_test/exact_inf/test_exact_3tax_3seg_2type_NonCoupledMCMC.log")


## Only Structured

doCompare("exact_test_no_reassort/simulate_3tax_3seg_no_reassort.log",
          "exact_test_no_reassort/test_exact_3tax_3seg_no_reassort.log")

## Only Reassortment

doCompare("exact_test_no_migration/simulate_5tax_8seg_no_migration.log",
          "exact_test_no_migration/test_exact_5tax_8seg_no_migration.log")

doCompare("exact_test_no_migration/simulate_5tax_2seg_no_migration.log",
          "exact_test_no_migration/test_exact_5tax_2seg_no_migration.log")

doCompare("exact_test_no_migration/simulate_5tax_3seg_no_migration.log",
          "exact_test_no_migration/test_exact_5tax_3seg_no_migration.log")

doCompare("exact_test_no_migration/simulate_2tax_3seg_no_migration.log",
          "exact_test_no_migration/test_exact_2tax_3seg_no_migration.log")

doCompare("exact_test_no_migration/simulate_2tax_3seg_no_migration.log",
          "exact_test_no_migration/test_exact_2tax_3seg_no_migration_1e8.log")


doCompare("exact_test_no_migration/simulate_2tax_3seg_no_migration.log",
          "../../CoalRe/validation/operators/test_compare_contemp2taxon3seg.log")



## Approximation 


doCompare("approximation/simulate_2tax_2seg_no_migration.log",
          "approximation/test_approx_2tax_2seg_no_migration.log")

doCompare("approximation/simulate_3tax_3seg_no_migration.log",
          "approximation/test_approx_3tax_3seg_no_migration.log")


doCompare("approximation/simulate_2tax_2seg_2type.log",
          "approximation/test_approx_2tax_2seg_2type.log")
