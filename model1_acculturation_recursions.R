# R code for Model 1 in Mesoudi, A. (2018) Migration, acculturation, and the maintenance of between-group cultural variation. PLOS ONE

# multi-trait island model with conformist acculturation---------
# there are s subpopulations and s traits
# initially, perfect between-group structure
# subpop 1 has all trait 1, subpop 2 has all trait 2, etc.
# each generation, m individuals migrate
# each subpop then undergoes conformist acculturation
# pick n random demonstrators, disproportionately copy the most common
# track Fst

# load packages and color schemes-------------------------
library(gtools)  # permutations function needed for multinomial theorem
library(ggplot2)  # various graphing packages
library(viridis)
library(cowplot)
library(directlabels)

# colorblind pallette
cbPalette <- c("#999999", "#E69F00", "#56B4E9", "#009E73", "#D55E00", "#CC79A7")

# define main recursion model function--------------------------

Model1.rec <- function(s, m, t.max, n, a, r, showplot) {

  # create matrix to hold current frequences: each trait is a col, each subpop is a row
  freqs <- matrix(0, nrow = s, ncol = s)
  freqs[row(freqs) == col(freqs)] <- 1  # initially start with complete BGV
  
  # also create list of s matrices, for holding freqs in each subpop over timesteps
  all.freqs <- rep(list(matrix(0, nrow = t.max, ncol = s)),s)
    
    # put initial freqs into all.freqs
    for (i in 1:s) {
      all.freqs[[i]][1,] <- freqs[i,]
    }
  
  # function for getting Fst from set of frequencies
  getFstFreqs <- function(freqs) {
    
    total.var <- 1 - sum(colMeans(freqs)^2)  # 1 - sum of squared means of each trait
    
    within.var <- mean(1 - rowSums(freqs^2)) # mean of each group's (1 - the sum of squared freq of each trait)
    
    (total.var - within.var) / total.var  # return Fst
  }
  
  Fst <- rep(0,t.max)  # for holding Fst
  Fst[1] <- getFstFreqs(freqs)  # get Fst for 1st gen
  
  # function for applying conformity. x is frequencies for each trait
  doConformity <- function(n, a, r, x) {
    
    s <- length(x) # number of traits
    new.x <- rep(0,s)  # for holding new freqs
    
    # get all combinations of k values that sum to n
    k <- permutations(n+1, s, repeats.allowed = T) - 1  # all permutations
    k <- k[rowSums(k) == n,]  # only those that sum to n
    
    # define function for getting multinomial coefficients
    chooseMN <- function(n, k) {
      num <- factorial(n)
      den <- factorial(k[1])
      for (i in 2:length(k)) {
        den <- den * factorial(k[i])
      }
      
      num / den
    }
    
    # cycle through each combination of k values to get random formation probabilities
    for (k.row in 1:nrow(k)) {
      
      # get multinomial coefficient
      MN.coefficient <- chooseMN(n, k[k.row,])
      
      # get trait term and powers, i.e. probability this combination of models forms
      meeting.prob <- prod(x^k[k.row,])

      # get X_i, conformity modifier, for each trait (i=trait index)
      k.max <- max(k[k.row,])  # maximum k
      pi <- sum(k[k.row,] == k.max)  # pi, number of kmax traits (1 if only one, >1 if ties)
      X <- rep(0,s)  # X_i
      
      for (i in 1:s) {
        
        if (k[k.row,i] == k.max) { # if this k_i is sole or joint maximum
          X[i] <- k[k.row,i] / n + a * (1 / pi - k[k.row,i] / n)  # then increase proportional to a
        } else {
          X[i] <- k[k.row,i] / n - a * k[k.row,i] / n  # otherwise, decrease proportional to a
        }
      }
      
      # update trait values
      new.x <- new.x + MN.coefficient * meeting.prob * X
    }
    
    # add assortation: fraction r don't change frequency, due to homogenous demonstrators
    new.x <- (1-r)*new.x + r*x
    
    # normalise frequencies to sum to 1
    new.x <- new.x / sum(new.x)
    
    new.x  # return updated freqs as output
  }
  
  for (t in 2:t.max) {  # generation loop
    
    old.freqs <- freqs  # save old freqs
    
    # migration
    for (i in 1:s) {  # i denotes subpop
      for (j in 1:s) {  # j denotes trait
        freqs[i,j] <- old.freqs[i,j]*(1-m) + mean(old.freqs[,j])*m
      }
    }
    
    # conformity
    for (i in 1:s) {  # cycle thru subpops
        freqs[i,] <- doConformity(n, a, r, x = freqs[i,])
    }
    
    # rounding errors sometimes mess up the freqs; round them here
    freqs <- round(freqs, 5)
    
    Fst[t] <- getFstFreqs(freqs)  # get Fst for this timestep
    
    # put this timestep's freqs into all.freqs
      for (i in 1:s) {
        all.freqs[[i]][t,] <- freqs[i,]
      }
    
  }
  
  # plots, depending on showplot
  
  # Fst only
  if (showplot == 1) {
    par(mfrow=c(1,1), cex=1.1)
    plot(Fst, type = 'l', ylab = expression(F[ST]), xlab = "generation", ylim = c(0,1), lwd=2)
  }
  
  # all freqs plus Fst
  if (showplot == 2) {

    par(mfcol=c(ceiling((s+1)/3),3), cex = 0.9)  # 3 columns
    for (i in 1:s) {
      plot(all.freqs[[i]][,1], type = 'l', ylab = "frequency", xlab = "generation", ylim = c(0,1), col = cbPalette[1], lwd = 2)
      for (j in 2:s) {
        lines(all.freqs[[i]][,j], lwd = 2, col = cbPalette[j])
      }
    }
    plot(Fst, type = 'l', ylab = expression(F[ST]), xlab = "generation", ylim = c(0,1), lwd = 2)
    
  }
  
  Fst  # return Fst as output
}

# run basic model here---------------

# s subpopulations/traits, prob of migration m, t.max timesteps, n demonstrators for conformity of strength a and assortment r
# showplot: 0=no plots, 1=Fst only, 2=all freqs plus Fst
output.Model1.rec <- Model1.rec(s = 5, m = 0.3, t.max = 50, n = 5, a = 0.6, r = 0, showplot = 1)

# create Fig 1: multi-line time-series plot for a at m=0.01, m=0.1 and m=0.3-----------------------

# plot 1, at m=0.01
s.lines <- 5  # common parameters, change these once here
m.lines <- 0.01
t.lines <- 400
n.lines <- 5
r.lines <- 0

plot1.line1.rec <- Model1.rec(s = s.lines, m = m.lines, t.max = t.lines, n = n.lines, a = 0.0, r = r.lines, showplot = 0)
plot1.line2.rec <- Model1.rec(s = s.lines, m = m.lines, t.max = t.lines, n = n.lines, a = 0.02, r = r.lines, showplot = 0)
plot1.line3.rec <- Model1.rec(s = s.lines, m = m.lines, t.max = t.lines, n = n.lines, a = 0.05, r = r.lines, showplot = 0)
plot1.line4.rec <- Model1.rec(s = s.lines, m = m.lines, t.max = t.lines, n = n.lines, a = 0.1, r = r.lines, showplot = 0)
plot1.line5.rec <- Model1.rec(s = s.lines, m = m.lines, t.max = t.lines, n = n.lines, a = 1, r = r.lines, showplot = 0)
plot1.lines.rec <- c(plot1.line1.rec, plot1.line2.rec, plot1.line3.rec, plot1.line4.rec, plot1.line5.rec)
plot1.lines.rec <- data.frame(Fst = plot1.lines.rec, timestep = rep(1:t.lines, 5), line.id = c(rep("a=0",t.lines),rep("a=0.02",t.lines),rep("a=0.05",t.lines),rep("a=0.1",t.lines),rep("a=1",t.lines)))

plot1.timeseries.rec <- ggplot(data = plot1.lines.rec, aes(x = timestep, y = Fst, color = line.id)) + geom_line(size = 0.7) + labs(x = "timestep", y=expression(F[ST])) + scale_y_continuous(limits = c(0,1)) + scale_x_continuous(limits=c(0,t.lines+100), breaks = seq(0, t.lines, by = 100)) + theme_classic() + ggtitle("m = 0.01") + theme(axis.line.x = element_line(colour = "black", size = 0.5), axis.line.y = element_line(colour = "black", size = 0.5), legend.title=element_blank(), plot.title = element_text(size = 10, hjust = 0.5), axis.title.y = element_text(angle=0, vjust=0.5)) + scale_fill_manual(values=cbPalette) + scale_colour_manual(values=cbPalette) + guides(col = guide_legend(reverse = TRUE))

plot1.timeseries.rec <- direct.label(plot1.timeseries.rec, list(last.points, cex = 0.6, hjust = -0.1))

# plot 2, at m=0.1
m.lines <- 0.1  # changed for plot2, otherwise same as plot1
t.lines <- 50

plot2.line1.rec <- Model1.rec(s = s.lines, m = m.lines, t.max = t.lines, n = n.lines, a = 0.0, r = r.lines, showplot = 0)
plot2.line2.rec <- Model1.rec(s = s.lines, m = m.lines, t.max = t.lines, n = n.lines, a = 0.2, r = r.lines, showplot = 0)
plot2.line3.rec <- Model1.rec(s = s.lines, m = m.lines, t.max = t.lines, n = n.lines, a = 0.4, r = r.lines, showplot = 0)
plot2.line4.rec <- Model1.rec(s = s.lines, m = m.lines, t.max = t.lines, n = n.lines, a = 0.6, r = r.lines, showplot = 0)
plot2.line5.rec <- Model1.rec(s = s.lines, m = m.lines, t.max = t.lines, n = n.lines, a = 1, r = r.lines, showplot = 0)
plot2.lines.rec <- c(plot2.line1.rec, plot2.line2.rec, plot2.line3.rec, plot2.line4.rec, plot2.line5.rec)
plot2.lines.rec <- data.frame(Fst = plot2.lines.rec, timestep = rep(1:t.lines, 5), line.id = c(rep("a=0",t.lines),rep("a=0.2",t.lines),rep("a=0.4",t.lines),rep("a=0.6",t.lines),rep("a=1",t.lines)))

plot2.timeseries.rec <- ggplot(data = plot2.lines.rec, aes(x = timestep, y = Fst, color = line.id)) + geom_line(size = 0.7) + labs(x = "timestep", y=expression(F[ST])) + scale_y_continuous(limits = c(0,1)) + scale_x_continuous(limits=c(0,t.lines*1.2), breaks = seq(0, t.lines, by = t.lines/5)) + theme_classic() + ggtitle("m = 0.1") + theme(axis.line.x = element_line(colour = "black", size = 0.5), axis.line.y = element_line(colour = "black", size = 0.5), legend.title=element_blank(), plot.title = element_text(size = 10, hjust = 0.5), axis.title.y = element_text(angle=0, vjust=0.5)) + scale_fill_manual(values=cbPalette) + scale_colour_manual(values=cbPalette) + guides(col = guide_legend(reverse = TRUE))

plot2.timeseries.rec <- direct.label(plot2.timeseries.rec, list(last.points, cex = 0.6, hjust = -0.1))

# plot 3, at m=0.3
m.lines <- 0.3  # changed for plot2, otherwise same as plot1
t.lines <- 50

plot3.line1.rec <- Model1.rec(s = s.lines, m = m.lines, t.max = t.lines, n = n.lines, a = 0.0, r = r.lines, showplot = 0)
plot3.line2.rec <- Model1.rec(s = s.lines, m = m.lines, t.max = t.lines, n = n.lines, a = 0.6, r = r.lines, showplot = 0)
plot3.line3.rec <- Model1.rec(s = s.lines, m = m.lines, t.max = t.lines, n = n.lines, a = 0.8, r = r.lines, showplot = 0)
plot3.line4.rec <- Model1.rec(s = s.lines, m = m.lines, t.max = t.lines, n = n.lines, a = 0.9, r = r.lines, showplot = 0)
plot3.line5.rec <- Model1.rec(s = s.lines, m = m.lines, t.max = t.lines, n = n.lines, a = 1, r = r.lines, showplot = 0)
plot3.lines.rec <- c(plot3.line1.rec, plot3.line2.rec, plot3.line3.rec, plot3.line4.rec, plot3.line5.rec)
plot3.lines.rec <- data.frame(Fst = plot3.lines.rec, timestep = rep(1:t.lines, 5), line.id = c(rep("a=0",t.lines),rep("a=0.6",t.lines),rep("a=0.8",t.lines),rep("a=0.9",t.lines),rep("a=1",t.lines)))

plot3.timeseries.rec <- ggplot(data = plot3.lines.rec, aes(x = timestep, y = Fst, color = line.id)) + geom_line(size = 0.7) + labs(x = "timestep", y=expression(F[ST])) + scale_y_continuous(limits = c(0,1)) + scale_x_continuous(limits=c(0,t.lines*1.2), breaks = seq(0, t.lines, by = t.lines/5)) + theme_classic() + ggtitle("m = 0.3") + theme(axis.line.x = element_line(colour = "black", size = 0.5), axis.line.y = element_line(colour = "black", size = 0.5), legend.title = element_blank(), plot.title = element_text(size = 10, hjust = 0.5), axis.title.y = element_text(angle=0, vjust=0.5)) + scale_fill_manual(values=cbPalette) + scale_colour_manual(values=cbPalette) + guides(col = guide_legend(reverse = TRUE))

plot3.timeseries.rec <- direct.label(plot3.timeseries.rec, list(last.points, cex = 0.6, hjust = -0.1))

# use cowplot to create composite plot
all.plots.timeseries.rec <- plot_grid(plot1.timeseries.rec, plot2.timeseries.rec, plot3.timeseries.rec, labels=c("A", "B","C"), ncol = 3, nrow = 1, label_size = 16, rel_widths = c(1,1,1))

# show on screen
all.plots.timeseries.rec

# save to file (in the working directory)
save_plot("fig1.tif", all.plots.timeseries.rec, units = "cm", base_width = 18, base_height = 6, device = "tiff", dpi = 300)



# create Fig 2: plot Fst at varying a and n, for m=0.01, m=0.1 and m=0.3-------------------------------------------------------

# first make function to create one a-n plot
# NB n can be a single value, or several (e.g. "n = c(3,5,7)")
# 'gaps' specifies number of data points making up each line (default=10)
# a counter shows progress

plotan.model1.rec <- function (s, m, t.max, n, r, a.max=1, gaps=10, showplot=TRUE, showlegend=TRUE) {
  
  a.seq <- seq(0, a.max, by=1/gaps)  # a varies from 0 to a.max
  Fst.plotan <- data.frame(n = rep(n, each = (gaps+1)), a = rep(a.seq, length(n)), Fst = 0)  # for holding Fst values
  Fst.plotan$n <- as.factor(Fst.plotan$n)  # to plot lines need to specify n as factor
  
  counter <- 0  # for showing progress on screen
  
  for (i in 1:length(n)) {
    for (j in 1:(gaps+1)) {
      Fst <- Model1.rec(s, m, t.max, n = n[i], a = a.seq[j], r, showplot = FALSE)
      Fst.plotan$Fst[(i-1)*(gaps+1)+j] <- Fst[t.max]
      counter <- counter + 1
      cat(counter, "/", length(n)*(gaps+1), "\n")
    }
  }
  
  if (showlegend == FALSE) {
    leg <- "none"
  } else {
    leg <- "right"
  }
  
  plotan.plot <- ggplot(data = Fst.plotan, aes(x = a, y = Fst, color = n)) + geom_line(size = 0.7, alpha = 0.75) + labs(x = "a", y=expression(F[ST])) + theme_classic() + theme(legend.position = leg, axis.title.y = element_text(angle=0, vjust=0.5)) + scale_colour_manual(values=cbPalette) + scale_y_continuous(limits = c(0,1)) + annotate("text", y=0.1, x=0.8, label = paste("m =", m, sep=" "), size = 3)
  
  if (showplot == TRUE) {
    print(plotan.plot)
  }
  
  # save to file
  ggsave("plotan.png", width = 20, height = 20, units = "cm", dpi = 300)
  
  list(Fst.plotan, plotan.plot)  # output data and plot
}

# create one plot here
Fst.data <- plotan.model1.rec(s = 5, m = 0.2, t.max = 50, r = 0, gaps = 30, n = c(3,5,7))

# three-paneled plot of a and n, for 3 diff m values
n.lines <- c(3,5,7,11)  # change common values here
s.lines <- 5
t.max.lines <- 500
r.lines <- 0
gaps.lines <- 30

# get plots
plot1.an.rec <- plotan.model1.rec(s = s.lines, m = 0.01, t.max = t.max.lines, r = r.lines, n = n.lines, gaps = gaps.lines, showplot = FALSE, showlegend = FALSE)[[2]]

plot2.an.rec <- plotan.model1.rec(s = s.lines, m = 0.1, t.max = t.max.lines, r = r.lines, n = n.lines, gaps = gaps.lines, showplot = FALSE, showlegend = FALSE)[[2]]

plot3.an.rec <- plotan.model1.rec(s = s.lines, m = 0.3, t.max = t.max.lines, r = r.lines, n = n.lines, gaps = gaps.lines, showplot = FALSE)[[2]]

# use cowplot to create composite plot
legend <- get_legend(plot3.an.rec)
plot3.an.rec <- plot3.an.rec + theme(legend.position="none")
all.plots.an.rec <- plot_grid(plot1.an.rec, plot2.an.rec, plot3.an.rec, legend, labels=c("A", "B", "C"), ncol = 4, nrow = 1, label_size = 16, rel_widths = c(2,2,2,0.5))

# show on screen
all.plots.an.rec

# save to file
save_plot("fig2.tif", all.plots.an.rec, units = "cm", base_width = 18, base_height = 6, device = "tiff", dpi = 300)

# create Fig 3: 2x2 heat map for a and m, for different n and r -------------------------------------------------------

# first create heatmap function to create one heatmap
# for a grid.size x grid.size grid (e.g. if grid.size=10, a 10x10 grid)
# a counter tracks progress (can take a while)

drawHeatMap.model1.rec <- function(s, n, r, t.max, m.max, a.max, grid.size, showplot = TRUE, showlegend = TRUE) {
  
  # initialise grid.size x grid.size dataframe 'heat'
  m.heat <- rep(0:(grid.size-1), grid.size) # for now these are 0 to (grid.size-1) to use in loops, convert to 0 to m.max below
  a.heat <- rep(0:(grid.size-1), each = grid.size)
  value.heat <- rep(0,(grid.size^2))
  heat <- data.frame(m.heat, a.heat, value.heat)  
  
  counter <- 0
  
  # add values to heat matrix 
  for (m.loop in 0:(grid.size-1))
  {
    for (a.loop in 0:(grid.size-1))
    {
      Fst <- Model1.rec(s=s, n=n, r=r, m = m.loop*m.max/(grid.size-1), a = a.loop*a.max/(grid.size-1), t.max = t.max, showplot = FALSE)
      heat$value.heat[heat$a.heat == a.loop & heat$m.heat == m.loop] <- Fst[t.max]
      counter <- counter + 1
      cat(counter, "/", grid.size^2, "\n")
    }
  }
  
  # convert a.heat and m.heat values
  heat$a.heat <- heat$a.heat * a.max / (grid.size-1)
  heat$m.heat <- heat$m.heat * m.max / (grid.size-1)
  
  # create heatmap
  
  leg <- "right"
  if (showlegend == FALSE) { leg <- "none" }

  heat.plot <- ggplot(data = heat, aes(x = a.heat, y = m.heat)) + geom_tile(aes(fill = value.heat)) + labs(x = "a, acculturation", y = "m, migration") + ggtitle(paste("n = ", n, sep = "")) + theme_bw(base_size = 18) + theme(plot.title = element_text(size = 16, hjust = 0.5), legend.text = element_text(size = 10), legend.position = leg, axis.text = element_text(size = 12)) + scale_fill_viridis(name=expression(F[ST]))

  # show on screen
  if (showplot == TRUE) {
    print(heat.plot)
  }
  
  # save to file
  ggsave("model1_rec_am_heatmap.png", width = 20, height = 20, units = "cm", dpi = 300)
  
  list(heat, heat.plot)  # output heat dataframe and plot
}

# create one heat plot here
heat <- drawHeatMap.model1.rec(s = 5, n = 5, r = 0.05, t.max = 100, m.max = 0.6, a.max = 1, grid.size = 10)


# 2x2 heat plot in Fig 3 for different n and r

# use below commands to get data, then create plots with correct legend/axis titles
s.heats <- 5  # change common values here
t.max.heats <- 5
m.max.heats <- 0.6
a.max.heats <- 1
grid.size.heats <- 30

# turn Fst values into discrete categories
heat.breaks <- c(-1,0,0.1,0.2,0.3,0.4,0.5,0.6,0.7,0.8,0.9,1)
heat.labels <- c("0","0-0.1","0.1-0.2","0.2-0.3","0.3-0.4","0.4-0.5","0.5-0.6","0.6-0.7","0.7-0.8","0.8-0.9","0.9-1")


heat.n3r0.rec <- drawHeatMap.model1.rec(s = s.heats, n = 3, r = 0, t.max = t.max.heats, m.max = m.max.heats, a.max = a.max.heats, grid.size = grid.size.heats, showplot = F)[[1]]
# topleft: title n=3,r=0, y-axis, no x-axis
heat.n3r0.rec$value.heat.discrete <- cut(round(heat.n3r0.rec$value.heat,3), breaks = heat.breaks, labels = heat.labels)
topleft.n3r0.rec <- ggplot(data = heat.n3r0.rec, aes(x = a.heat, y = m.heat)) + geom_tile(aes(fill = value.heat.discrete)) + labs(title = "n = 3, r = 0", y = "m, migration") + theme(legend.position = "none", axis.title.x = element_blank(), axis.title = element_text(size = 11), axis.line = element_blank(), axis.text.x = element_text(size = 8), axis.text.y = element_text(size = 8), plot.title = element_text(size = 10)) + scale_fill_viridis(name=expression(F[ST]), discrete = TRUE) + scale_x_continuous(breaks = seq(0,1,by=0.2), expand = c(0,0)) + scale_y_continuous(expand = c(0,0))

heat.n7r0.rec <- drawHeatMap.model1.rec(s = s.heats, n = 7, r = 0, t.max = t.max.heats, m.max = m.max.heats, a.max = a.max.heats, grid.size = grid.size.heats, showplot = F)[[1]]
# top right: title n=7,r=0, no y-axis, no x-axis
heat.n7r0.rec$value.heat.discrete <- cut(round(heat.n7r0.rec$value.heat,3), breaks = heat.breaks, labels = heat.labels)
topright.n7r0.rec <- ggplot(data = heat.n7r0.rec, aes(x = a.heat, y = m.heat)) + geom_tile(aes(fill = value.heat.discrete)) + labs(title = "n = 7, r = 0") + theme(legend.position = "none", axis.title.x = element_blank(), axis.title.y = element_blank(), axis.line = element_blank(), axis.text.x = element_text(size = 8), axis.text.y = element_text(size = 8), plot.title = element_text(size = 10)) + scale_fill_viridis(name=expression(F[ST]), discrete = TRUE) + scale_x_continuous(breaks = seq(0,1,by=0.2), expand = c(0,0)) + scale_y_continuous(expand = c(0,0))

heat.n3r05.rec <- drawHeatMap.model1.rec(s = s.heats, n = 3, r = 0.5, t.max = t.max.heats, m.max = m.max.heats, a.max = a.max.heats, grid.size = grid.size.heats, showplot = F)[[1]]
# bottom left: title n=3,r=0.5, both axes
heat.n3r05.rec$value.heat.discrete <- cut(round(heat.n3r05.rec$value.heat,3), breaks = heat.breaks, labels = heat.labels)
botleft.n3r05.rec <- ggplot(data = heat.n3r05.rec, aes(x = a.heat, y = m.heat)) + geom_tile(aes(fill = value.heat.discrete)) + labs(title = "n = 3, r = 0.5", y = "m, migration", x = "a, acculturation") + theme(legend.position = "none", axis.line = element_blank(), axis.title = element_text(size = 11), axis.text.x = element_text(size = 8), axis.text.y = element_text(size = 8), plot.title = element_text(size = 10)) + scale_fill_viridis(name=expression(F[ST]), discrete = TRUE) + scale_x_continuous(breaks = seq(0,1,by=0.2), expand = c(0,0)) + scale_y_continuous(expand = c(0,0))

heat.n7r05.rec <- drawHeatMap.model1.rec(s = s.heats, n = 7, r = 0.5, t.max = t.max.heats, m.max = m.max.heats, a.max = a.max.heats, grid.size = grid.size.heats, showplot = F)[[1]]
# bottom right: title n=7,r=0.5, no y-axis, x-axis, include legend to extract later
heat.n7r05.rec$value.heat.discrete <- cut(round(heat.n7r05.rec$value.heat,3), breaks = heat.breaks, labels = heat.labels)
botright.n7r05.rec <- ggplot(data = heat.n7r05.rec, aes(x = a.heat, y = m.heat)) + geom_tile(aes(fill = value.heat.discrete)) + labs(title = "n = 7, r = 0.5", x = "a, acculturation") + theme(axis.title.y = element_blank(), axis.line = element_blank(), axis.title = element_text(size = 11), legend.title = element_text(size = 10), legend.text = element_text(size = 7), axis.text.x = element_text(size = 8), axis.text.y = element_text(size = 8), plot.title = element_text(size = 10), legend.key.width = unit(0.3,"cm"), legend.key.height = unit(0.3,"cm")) + scale_fill_viridis(name=expression(F[ST]), discrete = TRUE) + scale_x_continuous(breaks = seq(0,1,by=0.2), expand = c(0,0)) + scale_y_continuous(expand = c(0,0))

# create 4 graphs with correct combination of axis labels and titles, extract legend, and arrange all in grid using cowplot
legend <- get_legend(botright.n7r05.rec)
main_grid <- plot_grid(topleft.n3r0.rec, topright.n7r0.rec, botleft.n3r05.rec, botright.n7r05.rec + theme(legend.position="none"), align = "hv")
full_grid <- plot_grid(legend, main_grid, rel_widths = c(0.1,1))

# show on screen
full_grid

# save to file
save_plot("fig3.tif", full_grid, units = "cm", base_width = 14.4, base_height = 13.3, device = "tiff", dpi = 300)


