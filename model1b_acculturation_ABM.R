# code for agent-based version of Model 1 in Supplementary Material of Mesoudi, A. (in prep) Migration, acculturation, and the maintenance of between-group cultural variation

# semi-spatial migration model:
# s sub-populations each with N individuals
# each individual has cultural variant v (v from 1 to s)
# start with complete BGV: all N in s=1 has v=1, s=2 have v=2 etc.
# each individual migrates to random other sub-pop with prob m
# return migration to keep N constant
# acculturation a, prob of adopting most common variant amongst n randomly chosen agents in sub-group
# a applies to all, not just migrants, to reflect multi-generational acculturation. No mutation so migration is the only way of changing variant
# measure is degree of BGV/totalV, cultural Fst

# PACKAGES ---------------------------------
library(ggplot2)
library(viridis)
library(cowplot)
library(directlabels)

# colorblind pallette
cbPalette <- c("#999999", "#E69F00", "#56B4E9", "#009E73", "#D55E00", "#CC79A7")

# MODEL-------------------------------------

InitialiseAgents.Model1.ABM <- function (N, s) {
  # set up dataframe with subpop labels, and variants in each sub-pop
  subpop <- rep(1:s, each = N)
  variant <- rep(1:s, each = N) # start with complete BGV
  #variant <- sample(1:s, size = N*s, replace = TRUE) # alternative: initially random
  data.frame(subpop, variant, stringsAsFactors = FALSE)
}

GetFst.Model1.ABM <- function(d) {
  
  s <- max(d$subpop)  # get s, no. of traits and subpops; no need for N as all are equal
  
  X <- 0
  for (i in 1:s) {
    X <- X + mean(d$variant == i)^2  # sum up all the squared means of each variant i
  }
  total.var <- 1 - X  # 1 - sum of squared means
  
  X <- rep(0,s)  # variation per group
  for (j in 1:s) {  # for each group j
    for (i in 1:s) {  # for each variant i
      X[j] <- X[j] + mean(d$variant[d$subpop == j] == i)^2  # sum the squared means of each variant i
    }
    X[j] <- 1 - X[j]  # 1 - sum of squared means within group j
  }
  within.var <- mean(X)  # mean of all variances across all groups
  
  Fst <- (total.var - within.var) / total.var  # return Fst
  
  if (is.na(Fst)) {  # if Fst is NaN due to only 1 variant present, set to 0
    Fst <- 0
  }
  
  Fst  # return Fst
}

Migration.Model1.ABM <- function (m, d) {
  
  probs <- runif(1:nrow(d))  # probabilities for each agent to compare against m
  migrants <- d$variant[probs < m]  # with prob m, add an agent's variant to list of migrants
  d$variant[probs < m] <- sample(migrants, length(migrants))  # put migrants randomly back into empty slots
  
  d # return d
}

Acculturation.Model1.ABM <- function (a, n, r, d) {
  
  # for each agent, with prob a, set variant to most common variant amongst n randomly chosen agents in their subpop
  
  # get s and N from d
  s <- max(d$subpop)
  N <- nrow(d) / s
  
  success.a <- runif(nrow(d))  # get vector with a probs for each agent
  if (sum(success.a < a) > 0) {  # proceed only if there are any acculturating agents

    # create matrix, 1st n columns are demonstrators, next s columns are freqs
    demonstrators <- matrix(ncol = n+2*s, nrow = sum(success.a < a))  # list of demonstrators for <a agents

    # fill with random demonstrators from same subpop
    for (j in 1:n) {
      for (i in 1:s) {
        if (sum(success.a < a & d$subpop == i) > 0) {
          demonstrators[d$subpop[success.a < a] == i,j] <- sample(d$variant[success.a < a & d$subpop == i], length(d$variant[success.a < a & d$subpop == i]), replace = F)
        }
      }
    }

    # assortation: set each demonstrator trait to agent's trait with prob r
    for (i in 1:s) {  # for each trait
      success.r <- runif(nrow(demonstrators))  # r probs for each copying agent
      # for those rows where success.r<r, set to variant of that agent
      demonstrators[success.r < r,i] <- d$variant[success.a < a][success.r < r]
    }
    
    # now tally number of each variant in each s column
    for (i in 1:s) {
      demonstrators[,n+i] <- rowSums(demonstrators[,1:n] == i)
    }

    # convert s columns into 0=not highest, 1=highest (can be more than one 1, if joint-highest)
    demonstrators[,(n+1):(s+n)] <- demonstrators[,(n+1):(s+n)] == do.call(pmax, as.data.frame(demonstrators[,(n+1):(s+n)]))

    # put the single highest trait number in the final s columns (doesn't do ties yet)
    for (i in 1:s) {
      demonstrators[,n+s+i] <- ifelse(demonstrators[,n+i]==1 & rowSums(demonstrators[,(n+1):(s+n)]) == 1, i, 0)
    }

    # for ties (if any exist), put random variants into the final column
    if (0 %in% rowSums(demonstrators[,(n+s+1):(n+2*s)])) {
      demonstrators[rowSums(demonstrators[,(n+s+1):(n+2*s),drop=FALSE]) == 0, ncol(demonstrators)] <- apply(demonstrators[rowSums(demonstrators[,(n+s+1):(n+2*s),drop=FALSE]) == 0,1:n, drop=FALSE], 1, sample, 1)
    }

    # replace variants with most common variant for each agent, with probability a
    d$variant[success.a < a] <- rowSums(demonstrators[,(n+s+1):(n+2*s)])
  }

  d # return d
}

Model1.ABM <- function (N, s, m, a, n, r, t.max, runs, showplot) {
  
  Fst <- rep(0, t.max) # variable to keep Fst per generation
  Fst_min <- rep(0, t.max) # variable to keep min Fst per generation
  Fst_max <- rep(0, t.max) # variable to keep max Fst per generation
  
  for (i in 1:runs) {  # for-loop to cycle thru runs
    d <- InitialiseAgents.Model1.ABM(N, s) # initialise agent dataframe d
    
    for (t in 1:t.max) {  # cycle thru timesteps
      new.Fst <- GetFst.Model1.ABM(d)
      Fst[t] <- Fst[t] + new.Fst # add new Fst to running total
      
      if (i == 1) {  # 1st run, set max and min Fst to first value
        Fst_max[t] <- new.Fst
        Fst_min[t] <- new.Fst
      } else {
          if (new.Fst > Fst_max[t]) {Fst_max[t] <- new.Fst}  # otherwise update max
          if (new.Fst < Fst_min[t]) {Fst_min[t] <- new.Fst}  # and update min
      }
        
      d <- Migration.Model1.ABM(m,d) # migration
      d <- Acculturation.Model1.ABM(a,n,r,d)  # acculturation
    }
  }
  
  Fst <- Fst / runs  # average across all runs
  Fst.data <- data.frame(Fst, Fst_min, Fst_max)  # combine into data frame
  
  if (showplot == TRUE) {
  plot <- ggplot(data = Fst.data, aes(x = 1:nrow(Fst.data), y = Fst)) + geom_line() + geom_ribbon(aes(ymin = Fst_min, ymax = Fst_max), alpha = 0.2) + labs(x = "timestep", y=expression(F[ST])) + scale_y_continuous(limits = c(0,1.001)) + theme_classic(base_size = 14) + theme(axis.line.x = element_line(colour = "black", size = 0.5), axis.line.y = element_line(colour = "black", size = 0.5), axis.title.y = element_text(angle=0, vjust=0.5))
  print(plot)
    }
  
  Fst.data  # return Fst.data
}

# run model, storing Fst data---------------------------------
output.Model1.ABM <- Model1.ABM(N = 500, s = 5, m = 0.1, a = 0.2, n = 5, r = 0, t.max = 100, runs = 5, showplot = TRUE)


# PLOTS--------------------------------------------

# Figure S4: time series plots for different values of a, at two migration rates---------------------------------------

# run model above to get Fst.data for a set of parameters, then plot all together
# plot 1, m=0.01
N.lines <- 1000  # common parameters, change these once here
s.lines <- 5
m.lines <- 0.01
n.lines <- 5
r.lines <- 0
t.lines <- 400
runs.lines <- 10

plot1.line1.ABM <- Model1.ABM(N = N.lines, s = s.lines, m = m.lines, a = 0.0, n = n.lines, r = r.lines, t.max = t.lines, runs = runs.lines, showplot = FALSE)
plot1.line1.ABM$line.id <- "a=0"
plot1.line2.ABM <- Model1.ABM(N = N.lines, s = s.lines, m = m.lines, a = 0.02, n = n.lines, r = r.lines, t.max = t.lines, runs = runs.lines, showplot = FALSE)
plot1.line2.ABM$line.id <- "a=0.02"
plot1.line3.ABM <- Model1.ABM(N = N.lines, s = s.lines, m = m.lines, a = 0.05, n = n.lines, r = r.lines, t.max = t.lines, runs = runs.lines, showplot = FALSE)
plot1.line3.ABM$line.id <- "a=0.05"
plot1.line4.ABM <- Model1.ABM(N = N.lines, s = s.lines, m = m.lines, a = 0.1, n = n.lines, r = r.lines, t.max = t.lines, runs = runs.lines, showplot = FALSE)
plot1.line4.ABM$line.id <- "a=0.1"
plot1.line5.ABM <- Model1.ABM(N = N.lines, s = s.lines, m = m.lines, a = 1, n = n.lines, r = r.lines, t.max = t.lines, runs = runs.lines, showplot = FALSE)
plot1.line5.ABM$line.id <- "a=1"

plot1.lines.ABM <- rbind(plot1.line1.ABM, plot1.line2.ABM, plot1.line3.ABM, plot1.line4.ABM, plot1.line5.ABM)
plot1.lines.ABM$timestep <- rep(1:t.lines, 5)

plot1.timeseries.ABM <- ggplot(data = plot1.lines.ABM, aes(x = timestep, y = Fst, color = line.id)) + geom_line(size = 1) + geom_ribbon(aes(ymin = Fst_min, ymax = Fst_max, fill = line.id), alpha = 0.1, show.legend = FALSE) + labs(x = "timestep", y=expression(F[ST])) + scale_y_continuous(limits = c(0,1)) + scale_x_continuous(limits=c(0,t.lines*1.2), breaks = seq(0, t.lines, by = t.lines/4)) + theme_classic(base_size = 20) + ggtitle("m = 0.01") + theme(axis.line.x = element_line(colour = "black", size = 0.5), axis.line.y = element_line(colour = "black", size = 0.5), legend.title=element_blank(), plot.title = element_text(size = 20, hjust = 0.5), axis.title.y = element_text(size=25, angle=0, vjust=0.5)) + scale_fill_manual(values=cbPalette) + scale_colour_manual(values=cbPalette) + guides(col = guide_legend(reverse = TRUE))

plot1.timeseries.ABM <- direct.label(plot1.timeseries.ABM, list(last.points, cex = 1.4, hjust = -0.1))

# plot 2, m=0.1
m.lines <- 0.1  # these change, others from plot1
t.lines <- 50

plot2.line1.ABM <- Model1.ABM(N = N.lines, s = s.lines, m = m.lines, a = 0.0, n = n.lines, r = r.lines, t.max = t.lines, runs = runs.lines, showplot = FALSE)
plot2.line1.ABM$line.id <- "a=0"
plot2.line2.ABM <- Model1.ABM(N = N.lines, s = s.lines, m = m.lines, a = 0.2, n = n.lines, r = r.lines, t.max = t.lines, runs = runs.lines, showplot = FALSE)
plot2.line2.ABM$line.id <- "a=0.2"
plot2.line3.ABM <- Model1.ABM(N = N.lines, s = s.lines, m = m.lines, a = 0.4, n = n.lines, r = r.lines, t.max = t.lines, runs = runs.lines, showplot = FALSE)
plot2.line3.ABM$line.id <- "a=0.4"
plot2.line4.ABM <- Model1.ABM(N = N.lines, s = s.lines, m = m.lines, a = 0.6, n = n.lines, r = r.lines, t.max = t.lines, runs = runs.lines, showplot = FALSE)
plot2.line4.ABM$line.id <- "a=0.6"
plot2.line5.ABM <- Model1.ABM(N = N.lines, s = s.lines, m = m.lines, a = 1, n = n.lines, r = r.lines, t.max = t.lines, runs = runs.lines, showplot = FALSE)
plot2.line5.ABM$line.id <- "a=1"

plot2.lines.ABM <- rbind(plot2.line1.ABM, plot2.line2.ABM, plot2.line3.ABM, plot2.line4.ABM, plot2.line5.ABM)
plot2.lines.ABM$timestep <- rep(1:t.lines, 5)

plot2.timeseries.ABM <- ggplot(data = plot2.lines.ABM, aes(x = timestep, y = Fst, color = line.id)) + geom_line(size = 1) + geom_ribbon(aes(ymin = Fst_min, ymax = Fst_max, fill = line.id), alpha = 0.1, show.legend = FALSE) + labs(x = "timestep", y=expression(F[ST])) + scale_y_continuous(limits = c(0,1)) + scale_x_continuous(limits=c(0,t.lines*1.2), breaks = seq(0, t.lines, by = t.lines/5)) + theme_classic(base_size = 20) + ggtitle("m = 0.1") + theme(axis.line.x = element_line(colour = "black", size = 0.5), axis.line.y = element_line(colour = "black", size = 0.5), legend.title=element_blank(), plot.title = element_text(size = 20, hjust = 0.5), axis.title.y = element_text(size=25, angle=0, vjust=0.5)) + scale_fill_manual(values=cbPalette) + scale_colour_manual(values=cbPalette) + guides(col = guide_legend(reverse = TRUE))

plot2.timeseries.ABM <- direct.label(plot2.timeseries.ABM, list(last.points, cex = 1.4, hjust = -0.1))

# use cowplot to create composite plot
all.plots.timeseries.ABM <- plot_grid(plot1.timeseries.ABM, plot2.timeseries.ABM, labels=c("A", "B"), ncol = 2, nrow = 1, label_size = 30, rel_widths = c(1,1))

# show on screen
all.plots.timeseries.ABM

# save to file
save_plot("timeseries_model1_ABM.png", all.plots.timeseries.ABM, base_width = 30, base_height = 15, units = "cm")


# Figure S5: Plot values of Fst at varying a and n------------------------------

plotan.model1.ABM <- function (N, s, m, n, r, a.max = 1, t.max, runs, gaps = 10, showplot=TRUE, showlegend=TRUE) {
  
  a.seq <- seq(0, a.max, by=1/gaps)  # a varies from 0 to a.max
  Fst.plotan <- data.frame(n = rep(n, each = (gaps+1)), a = rep(a.seq, length(n)), Fst = 0)  # for holding Fst values
  Fst.plotan$n <- as.factor(Fst.plotan$n)  # to plot lines need to specify n as factor
  
  counter <- 0  # for showing progress on screen
  
  for (i in 1:length(n)) {
    for (j in 1:(gaps+1)) {
      Fst <- Model1.ABM(N, s, m, n = n[i], r, a = a.seq[j], t.max, runs, showplot = FALSE)
      Fst.plotan$Fst[(i-1)*(gaps+1)+j] <- Fst$Fst[t.max]
      counter <- counter + 1
      cat(counter, "/", length(n)*(gaps+1), "\n")
    }
  }
  
  if (showlegend == FALSE) {
    leg <- "none"
  } else {
    leg <- "right"
  }
  
  plotan.plot <- ggplot(data = Fst.plotan, aes(x = a, y = Fst, color = n)) + geom_line(size = 1.7) + labs(x = "a", y=expression(F[ST])) + theme_classic(base_size = 20) + theme(legend.position = leg, axis.title.y = element_text(angle=0, vjust=0.5)) + scale_colour_manual(values=cbPalette) + scale_y_continuous(limits = c(0,1)) + annotate("text", y=0.1, x=0.8, label = paste("m =", m, sep=" "), size = 6)
  
  if (showplot == TRUE) {
    print(plotan.plot)
  }
  
  # save to file
  ggsave("plotan_model1_ABM.png", width = 20, height = 20, units = "cm", dpi = 300)
  
  list(Fst.plotan, plotan.plot)  # output data and plot
}

# draw single plot
Fst.data <- plotan.model1.ABM(N = 300, s = 5, m = 0.1, n = c(3,5,7,11), r = 0, t.max = 500, runs = 10)


# three-paneled plot of a and n, for 3 diff m values
N.lines <- 1000  # change common values here
s.lines <- 5
r.lines <- 0
t.max.lines <- 500
runs.lines <- 10
n.lines <- c(13,30,100)
plot1.an <- plotan.model1.ABM(N=N.lines, s=s.lines, m = 0.01, r = r.lines, t.max=t.max.lines, n = n.lines, runs=runs.lines, showplot = FALSE, showlegend = FALSE)[[2]]
plot2.an <- plotan.model1.ABM(N=N.lines, s=s.lines, m = 0.1, r = r.lines, t.max=t.max.lines, n = n.lines, runs=runs.lines, showplot = FALSE, showlegend = FALSE)[[2]]
plot3.an <- plotan.model1.ABM(N=N.lines, s=s.lines, m = 0.2, r = r.lines, t.max=t.max.lines, n = n.lines, runs=runs.lines, showplot = FALSE)[[2]]

# use cowplot to create composite plot
legend <- get_legend(plot3.an)
plot3.an <- plot3.an + theme(legend.position="none")
all.plots.an <- plot_grid(plot1.an, plot2.an, plot3.an, legend, labels=c("A", "B", "C"), ncol = 4, nrow = 1, label_size = 24, rel_widths = c(2,2,2,0.5))

# show on screen
all.plots.an

# save to file
save_plot("plots_an_model1_ABM.png", all.plots.an, base_width = 35, base_height = 10, units = "cm")


# Figure S6: Heat map of m and a, at different n and r -------------------------------------------------------

# first make a function for a single heatplot
# for a grid.size x grid.size grid (e.g. if grid.size=10, a 10x10 grid)

drawHeatMap.model1.ABM <- function(N, n, s, r, t.max, m.max, a.max, runs, grid.size, showplot = TRUE, showlegend = TRUE) {
  
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
      Fst <- Model1.ABM(N=N, s=s, r=r, n=n, m = m.loop*m.max/(grid.size-1), a = a.loop*a.max/(grid.size-1), t.max = t.max, runs=runs, showplot = FALSE)
      heat$value.heat[heat$a.heat == a.loop & heat$m.heat == m.loop] <- Fst$Fst[t.max]
      counter <- counter + 1
      cat(counter, "/", grid.size^2, "\n")
    }
  }
  
  # convert c.heat and u.heat values
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
  ggsave("heatmap_model1_ABM.png", width = 20, height = 20, units = "cm", dpi = 300)
  
  list(heat, heat.plot)  # output heat dataframe and plot
}

drawHeatMap.model1.ABM(N = 300, s = 5, r = 0, n = 100, t.max = 500, runs = 10, m.max = 0.6, a.max = 1, grid.size = 30)

# now create 2x2 heat plot

# use below commands to get data, then create plots with correct legend/axis titles below
s.heats <- 5  # change common values here
t.max.heats <- 400
m.max.heats <- 0.6
a.max.heats <- 1
N.heats <- 300
runs.heats <- 5
grid.size.heats <- 30

heat.n3r0.ABM <- drawHeatMap.model1.ABM(runs=runs.heats, N=N.heats, s = s.heats, n = 3, r = 0, t.max = t.max.heats, m.max = m.max.heats, a.max = a.max.heats, grid.size = grid.size.heats, showplot = F)[[1]]
# top left:
topleft.n3r0.ABM <- ggplot(data = heat.n3r0.ABM, aes(x = a.heat, y = m.heat)) + geom_tile(aes(fill = value.heat)) + labs(title = "n = 3, r = 0", y = "m, migration") + theme(legend.position = "none", axis.title.x = element_blank(), axis.title = element_text(size = 18), axis.line = element_blank()) + scale_fill_viridis(name=expression(F[ST])) + scale_x_continuous(breaks = seq(0,1,by=0.2), expand = c(0,0)) + scale_y_continuous(expand = c(0,0))

heat.n7r0.ABM <- drawHeatMap.model1.ABM(runs=runs.heats, N=N.heats, s = s.heats, n = 7, r = 0, t.max = t.max.heats, m.max = m.max.heats, a.max = a.max.heats, grid.size = grid.size.heats, showplot = F)[[1]]
# top right:
topright.n7r0.ABM <- ggplot(data = heat.n7r0.ABM, aes(x = a.heat, y = m.heat)) + geom_tile(aes(fill = value.heat)) + labs(title = "n = 7, r = 0") + theme(legend.position = "none", axis.title.x = element_blank(), axis.title.y = element_blank(), axis.line = element_blank()) + scale_fill_viridis(name=expression(F[ST])) + scale_x_continuous(breaks = seq(0,1,by=0.2), expand = c(0,0)) + scale_y_continuous(expand = c(0,0))

heat.n3r05.ABM <- drawHeatMap.model1.ABM(runs=runs.heats, N=N.heats, s = s.heats, n = 3, r = 0.5, t.max = t.max.heats, m.max = m.max.heats, a.max = a.max.heats, grid.size = grid.size.heats, showplot = F)[[1]]
# bottom left:
botleft.n3r05.ABM <- ggplot(data = heat.n3r05.ABM, aes(x = a.heat, y = m.heat)) + geom_tile(aes(fill = value.heat)) + labs(title = "n = 3, r = 0.5", y = "m, migration", x = "a, acculturation") + theme(legend.position = "none", axis.line = element_blank(), axis.title = element_text(size = 18)) + scale_fill_viridis(name=expression(F[ST])) + scale_x_continuous(breaks = seq(0,1,by=0.2), expand = c(0,0)) + scale_y_continuous(expand = c(0,0))

heat.n7r05.ABM <- drawHeatMap.model1.ABM(runs=runs.heats, N=N.heats, s = s.heats, n = 7, r = 0.5, t.max = t.max.heats, m.max = m.max.heats, a.max = a.max.heats, grid.size = grid.size.heats, showplot = F)[[1]]
# bottom right:
botright.n7r05.ABM <- ggplot(data = heat.n7r05.ABM, aes(x = a.heat, y = m.heat)) + geom_tile(aes(fill = value.heat)) + labs(title = "n = 7, r = 0.5", x = "a, acculturation") + theme(axis.title.y = element_blank(), axis.line = element_blank(), axis.title = element_text(size = 18), legend.title = element_text(size = 18)) + scale_fill_viridis(name=expression(F[ST])) + scale_x_continuous(breaks = seq(0,1,by=0.2), expand = c(0,0)) + scale_y_continuous(expand = c(0,0))


# create 4 graphs with correct combination of axis labels and titles, extract legend, and arrange all in grid using cowplot
legend <- get_legend(botright.n7r05.ABM)
main_grid <- plot_grid(topleft.n3r0.ABM, topright.n7r0.ABM, botleft.n3r05.ABM, botright.n7r05.ABM + theme(legend.position="none"), align = "hv")
full_grid <- plot_grid(legend, main_grid, rel_widths = c(0.1,1))

save_plot("heat_grid_ABM.png", full_grid, base_height = 7)

