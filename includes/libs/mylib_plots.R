
# plots.plot_ScatterMatrix
# plots.ggplot_theme
# plots.ggplot_theme_casual
# plots.install_xkcd
# plots.ggplot_labels
# plots.hclust

############################################################################################
############################################################################################

plots.plot_ScatterMatrix <- function(df, columnsToPlot)
{
    library(GGally, quietly=T)
    p <- ggpairs(data=df, columns=columnsToPlot, upper="blank", lower=list(continuous="points"), legends=F, mapping=ggplot2::aes(colour=tag))
    
    for (i in 1:length(columnsToPlot))
    {
        # Address only the diagonal elements
        # Get plot out of matrix
        inner <- getPlot(p, i, i);
        
        # Add any ggplot2 settings you want (blank grid here)
        inner <- inner + theme(panel.grid=element_blank()) + theme(axis.text.x=element_blank())
        
        # Put it back into the matrix
        p <- putPlot(p, inner, i, i)
        
        for (j in 1:length(columnsToPlot))
        {
            if((i==1 & j==1))
            {
                # Move legend right
                inner <- getPlot(p, i, j)
                inner <- inner + theme(legend.position=c(length(columnsToPlot)-0.25,0.50))
                p <- putPlot(p, inner, i, j)
            }
            else{
                # Delete legend
                inner <- getPlot(p, i, j)
                inner <- inner + theme(legend.position="none")
                p <- putPlot(p, inner, i, j)
            }
        }
    }
    return(p)
}

############################################################################################

plots.ggplot_theme <- function(g)
{
    g <- g + theme_classic()
    g <- g + theme(axis.text=element_text(size=16, face="bold"), axis.title.y=element_text(vjust=0.9), axis.title.x=element_text(vjust=-0.6))
    g <- g + theme(axis.title=element_text(size=16, face="bold"))
    g <- g + theme(panel.border=element_blank())
    g <- g + theme(axis.line.x=element_line(color="black", size=0.7), axis.line.y=element_line(color="black", size=0.7))
    return(g)
}

############################################################################################

plots.ggplot_theme_casual <- function(g)
{
    if(!"xkcd" %in% rownames(installed.packages())) 
    {
      install.packages("xkcd")
    }
    library("xkcd")
    return (theme(axis.line.x=element_line(size=0.5, colour="black"), axis.line.y=element_line(size=0.5, colour="black"), axis.text.x=element_text(colour="black", size=10), axis.text.y=element_text(colour="black", size=10), panel.background=element_blank(), plot.title=element_text(family="xkcd"), text=element_text(family="xkcd")))
}

############################################################################################

plots.install_xkcd <- function()
{
  install.packages("xkcd")
  if(!"sysfonts" %in% rownames(installed.packages())) 
  {
    install.packages("sysfonts")
  }
  library(sysfonts)
  download.file("http://simonsoftware.se/other/xkcd.ttf", dest="xkcd.ttf", mode="wb")
  
  # The following steps are adapted from the xkcd package vignette
  # for a Mac OS X machine.
  system("cp xkcd.ttf ~/Library/Fonts")
  system("rm xkcd.ttf")
  stopifnot("xkcd.ttf" %in% font.files())
}

############################################################################################

plots.ggplot_labels <- function(g, xlabel, ylabel, title)
{
    g <- g + labs(x=paste0(xlabel), y=ylabel)
    g <- g + ggtitle(title)
    return(g)
}

############################################################################################

plots.hclust <- function(hc)
{
    library(ggplot2)
    # hc is the object from hclust()
    dhc <- as.dendrogram(hc)
    ddata <- dendro_data(dhc, type="rectangle")
    p <- ggplot(segment(ddata))
    p <- p + geom_segment(aes(x=x, y=y, xend=xend, yend=yend), size=0.7)
    p <- p + coord_flip()
    p <- p + scale_y_reverse(expand=c(1.3,0), breaks=seq(0,1,by=0.1))
    p <- p + scale_size("n")
    p <- p + geom_text(data=label(ddata), aes(x=x, y=y, label=label), , size=6, hjust=0, fontface="bold")
    p <- p + theme_classic()
    p <- p + theme(axis.ticks.y=element_blank(), axis.text.y=element_blank(), axis.title=element_text(size=16,face="bold",vjust=-0.65))
    p <- p + theme(axis.text.x=element_text(angle=90, hjust=0, size=16, vjust=0.4, face="bold"))
    p <- p + theme(axis.title=element_text(vjust=0.9,face="bold",size=16))
    p
}

############################################################################################

