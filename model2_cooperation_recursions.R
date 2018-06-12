# code for Model 2 in Mesoudi, A. (in prep) Migration, acculturation, and the maintenance of between-group cultural variation

# multi-trait island model with conformist acculturation for a cooperative trait---------
# there is one focal subpopulations and 2 traits, cooperate or defect
# initially, focal subpopulation is 100% cooperate (p=1)
# assume that all other subpops are full of defectors (not simulated)
# each generation, m defectors migrate in with probability proportional to the difference between subpop mean fitness W and mean fitness of the all-defector meta-population (where p=0, and W=1)
# subpop then undergoes conformist acculturation as in Model 1
# pick n random demonstrators, disproportionately copy the most common trait with probability a
# then payoff-biased social learning within the sub-population with probability L
# track frequency of cooperation in the face of defecting immigrants

# load packages and color schemes-------------------------
library(gtools)  # permutations function needed for multinomial theorem
library(ggplot2)  # various graphing packages
library(viridis)
library(cowplot)
library(directlabels)

# colorblind palette
cbPalette <- c("#999999", "#E69F00", "#56B4E9", "#009E73", "#D55E00", "#CC79A7")

# define main recursion model function--------------------------

Model2.rec <- function(m, t.max, n, a, r, b, c, u, v, L, showplot) {

  # fitnesses checks:
  stopifnot(c+u<=1, c<b, v<b, (b-c)>=(u+v))
  
  # vector of frequencies; each row is a timestep
  p <- rep(0,t.max)
  p[1] <- 1  # initial freq of cooperation
  
  # function for getting fitnesses w_c, w_d and W from trait frequency p
  getFitness <- function(b, c, u, v, p) {
    w_c <- 1 + b*p - c - u*(1-p)
    w_d <- 1 + b*p - v*p
    W <- w_c*p + w_d*(1-p)
    list(w_c=w_c, w_d=w_d, W=W)
  }
  
  # function for payoff-biased migration
  # p is freq of cooperators at time t
  doMigration <- function(m, p) {
    
    w <- getFitness(b,c,u,v,p)  # fitnesses of subpop
    
    # get mu, payoff biased migration modifier, to ensure m ranges from 0-1
    # fitness assumptions ensure that mean fitness can never be smaller than mean fitness of all-defector population, W=1, so can never be negative by assumption
    # maximum difference is when p=1, where w_c=1+b-c; W-1 = b-c; mu is reciprocal of this
    mu <- 1/(b-c)
    
    new.p <- p*(1-m*mu*(w$W-1))  # migration rate weighted by fitness difference between this subpop mean fitness and all-D, which are always 1. Modifier mu keeps it 0-1
    
    new.p  # return updated frequency
  }
  
  # function for applying conformity within subpop. p is frequency of trait
  doConformity <- function(n, a, r, p) {
    
    new.p <- 0  # for holding new trait freq
    
    # get all combinations of k values that sum to n
    k <- permutations(n+1, 2, repeats.allowed = T) - 1  # all permutations
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
      # unlike Model 1, only one other trait with freq 1-p, so no need to cycle thru traits
      meeting.prob <- (p^k[k.row,1])*((1-p)^k[k.row,2])

      # get X, conformity modifier for this trait
      k.max <- max(k[k.row,])  # maximum k
      
      # first set X to unbiased transmission frequencies
      X <- k[k.row,1] / n
      
      # if there is a single maximum k, i.e. single most common trait, and 0 < k_i < n
      if (sum(k[k.row,] == k.max) == 1 && k[k.row,1] > 0 && k[k.row,1] < n) {  
        
        if (k[k.row,1] == k.max) {  # if this k_i is the sole maximum
          X <- X + a * (n - k[k.row,1]) / n  # then increase proportional to a
        } else {
          X <- X - a * k[k.row,1] / n  # otherwise, decrease proportional to a
        }
      }
      
      # update trait values
      new.p <- new.p + MN.coefficient * meeting.prob * X
    }
    
    # add assortation: fraction r don't change frequency, due to homogenous demonstrators
    new.p <- (1-r)*new.p + r*p
    
    new.p  # return updated trait freq as output
  }
  
  # function for payoff-biased learning within subpop. p is frequency of trait
  doPayoffLearning <- function(L, p) {
    
    w <- getFitness(b,c,u,v,p)  # get fitnesses
    
    # get gamma, payoff biased copying modifier
    # is reciprocal of largest of fitness difference at p=1 (v-c) or p=0 (c+u)
    gam <- 1/max(c+u,v-c)
    
    new.p <- p + L*p*(1-p)*gam*(w$w_c - w$w_d)  # update frequencies
    
    new.p  # return p
  }
  
  for (t in 2:t.max) {  # generation loop
    
    # payoff-biased migration
    p[t] <- doMigration(m, p = p[t-1])
    
    # conformity
    p[t] <- doConformity(n, a, r, p = p[t])
    
    # payoff-biased social learning
    p[t] <- doPayoffLearning(L, p = p[t])
    
    # rounding errors sometimes mess up the freqs; round them here
    #p[t] <- round(p[t], 5)
    
  }
  
  if (showplot == 1) {
    plot(p, type = 'l', ylim = c(0,1), ylab = "freq of cooperators", xlab = "time", col = "blue", lwd = 2)
  }
  
  p  # return trait frequencies
}

# run basic model here---------------

# prob of migration m, t.max timesteps, n demonstrators for conformity of strength a and assortment r
# fitness parameters b, c, u and v; L is strength of within-group payoff biased social learning
# showplot = 1 displays the plot
output.Model2.rec <- Model2.rec(m = 0.1, t.max = 150, n = 5, a = 0.0, r = 0.0, b = 1, c = 0.2, u = 0.1, v = 0.5, L = 0.7, showplot = 1)


# create Figure 4: multi-line time-series plot for a and L-----------------------

# plot 1, at m=0.1, a=0 only
m.lines <- 0.1
t.lines <- 6000
n.lines <- 5
r.lines <- 0
L.lines <- 0
b.lines <- 1
c.lines <- 0.2
u.lines <- 0.1
v.lines <- 0.5

plot1.line1.rec <- Model2.rec(m = m.lines, t.max = t.lines, n = n.lines, r = r.lines, b = b.lines, c = c.lines, u = u.lines, v = v.lines, L = L.lines, showplot = 0, a = 0.0)
plot1.lines.rec <- data.frame(p = plot1.line1.rec, timestep = rep(1:t.lines), line.id = c(rep("a=0",t.lines)))

plot1.timeseries.rec <- ggplot(data = plot1.lines.rec, aes(x = timestep, y = p, color = line.id)) + geom_line(size = 1.5) + labs(x = "timestep", y="freq of cooperators, p") + scale_y_continuous(limits = c(0,1)) + scale_x_continuous(limits=c(0,t.lines+t.lines/5), breaks = seq(0, t.lines, by = t.lines / 3)) + theme_classic(base_size = 20) + ggtitle("m = 0.1, a = 0, L = 0") + theme(axis.line.x = element_line(colour = "black", size = 0.5), axis.line.y = element_line(colour = "black", size = 0.5), legend.title=element_blank(), plot.title = element_text(size = 20, hjust = 0.5), axis.title.y = element_text(size=25, vjust=0.5)) + scale_fill_manual(values=cbPalette) + scale_colour_manual(values=cbPalette) + guides(col = guide_legend(reverse = TRUE))

plot1.timeseries.rec <- direct.label(plot1.timeseries.rec, list(last.points, cex = 1.4, hjust = -0.1))

# plot 2, varying a
t.lines <- 150
plot2.line1.rec <- Model2.rec(m = m.lines, t.max = t.lines, n = n.lines, r = r.lines, b = b.lines, c = c.lines, u = u.lines, v = v.lines, L = L.lines, showplot = 0, a = 0.1)
plot2.line2.rec <- Model2.rec(m = m.lines, t.max = t.lines, n = n.lines, r = r.lines, b = b.lines, c = c.lines, u = u.lines, v = v.lines, L = L.lines, showplot = 0, a = 0.3)
plot2.line3.rec <- Model2.rec(m = m.lines, t.max = t.lines, n = n.lines, r = r.lines, b = b.lines, c = c.lines, u = u.lines, v = v.lines, L = L.lines, showplot = 0, a = 0.4)
plot2.line4.rec <- Model2.rec(m = m.lines, t.max = t.lines, n = n.lines, r = r.lines, b = b.lines, c = c.lines, u = u.lines, v = v.lines, L = L.lines, showplot = 0, a = 1)

plot2.lines.rec <- c(plot2.line1.rec, plot2.line2.rec, plot2.line3.rec, plot2.line4.rec)
plot2.lines.rec <- data.frame(p = plot2.lines.rec, timestep = rep(1:t.lines, 4), line.id = c(rep("a=0.1",t.lines),rep("a=0.3",t.lines),rep("a=0.4",t.lines),rep("a=1",t.lines)))

plot2.timeseries.rec <- ggplot(data = plot2.lines.rec, aes(x = timestep, y = p, color = line.id)) + geom_line(size = 1.5) + labs(x = "timestep", y="freq of cooperators, p") + scale_y_continuous(limits = c(0,1)) + scale_x_continuous(limits=c(0,t.lines), breaks = seq(0, t.lines, by = t.lines / 3)) + theme_classic(base_size = 20) + ggtitle("m = 0.1, L = 0") + theme(axis.line.x = element_line(colour = "black", size = 0.5), axis.line.y = element_line(colour = "black", size = 0.5), legend.title=element_blank(), plot.title = element_text(size = 20, hjust = 0.5), axis.title.y = element_text(size=25, vjust=0.5)) + scale_fill_manual(values=cbPalette[-1]) + scale_colour_manual(values=cbPalette[-1]) + guides(col = guide_legend(reverse = TRUE))


# plot 3, at a=0 and varying L
t.lines <- 150
plot3.line1.rec <- Model2.rec(m = m.lines, t.max = t.lines, n = n.lines, r = r.lines, b = b.lines, c = c.lines, u = u.lines, v = v.lines, L = 0.1, showplot = 0, a = 0.0)
plot3.line2.rec <- Model2.rec(m = m.lines, t.max = t.lines, n = n.lines, r = r.lines, b = b.lines, c = c.lines, u = u.lines, v = v.lines, L = 0.5, showplot = 0, a = 0.0)
plot3.line3.rec <- Model2.rec(m = m.lines, t.max = t.lines, n = n.lines, r = r.lines, b = b.lines, c = c.lines, u = u.lines, v = v.lines, L = 0.6, showplot = 0, a = 0.0)
plot3.line4.rec <- Model2.rec(m = m.lines, t.max = t.lines, n = n.lines, r = r.lines, b = b.lines, c = c.lines, u = u.lines, v = v.lines, L = 1, showplot = 0, a = 0.0)

plot3.lines.rec <- c(plot3.line1.rec, plot3.line2.rec, plot3.line3.rec, plot3.line4.rec)
plot3.lines.rec <- data.frame(p = plot3.lines.rec, timestep = rep(1:t.lines, 4), line.id = c(rep("L=0.1",t.lines),rep("L=0.5",t.lines),rep("L=0.6",t.lines),rep("L=1",t.lines)))

plot3.timeseries.rec <- ggplot(data = plot3.lines.rec, aes(x = timestep, y = p, color = line.id)) + geom_line(size = 1.5) + labs(x = "timestep", y="freq of cooperators, p") + scale_y_continuous(limits = c(0,1)) + scale_x_continuous(limits=c(0,t.lines), breaks = seq(0, t.lines, by = t.lines / 3)) + theme_classic(base_size = 20) + ggtitle("m = 0.1, a = 0") + theme(axis.line.x = element_line(colour = "black", size = 0.5), axis.line.y = element_line(colour = "black", size = 0.5), legend.title=element_blank(), plot.title = element_text(size = 20, hjust = 0.5), axis.title.y = element_text(size=25, vjust=0.5)) + scale_fill_manual(values=cbPalette[-1]) + scale_colour_manual(values=cbPalette[-1]) + guides(col = guide_legend(reverse = TRUE))

# use cowplot to create composite plot
# use cowplot to create composite plot
all.plots.timeseries.rec <- plot_grid(plot1.timeseries.rec, plot2.timeseries.rec, plot3.timeseries.rec, labels=c("A", "B", "C"), ncol = 3, nrow = 1, label_size = 30, rel_widths = c(1,1.2,1.2))

# save to file
save_plot("timeseries_model2_rec.png", all.plots.timeseries.rec, base_width = 40, base_height = 15, units = "cm")


# create Figure 5: 2 heatmaps for m & a at L=0 and L=1------------------------

# first create function for a single a-m heatmap
# for a grid.size x grid.size grid (e.g. if grid.size=10, a 10x10 grid)
# a counter tracks progress (can take a while)

drawHeatMap.am.model2.rec <- function(m.max, t.max, n, a.max, r, b, c, u, v, L, grid.size, showplot = TRUE, showlegend = TRUE) {
  
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
      p <- Model2.rec(m = m.loop*m.max/(grid.size-1), t.max = t.max, n=n, a = a.loop*a.max/(grid.size-1), r=r, b=b, c=c, u=u, v=v, L=L, showplot = FALSE)
      heat$value.heat[heat$a.heat == a.loop & heat$m.heat == m.loop] <- p[t.max]
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

  # turn Fst values into discrete categories
  heat.breaks <- c(-1,0,0.1,0.2,0.3,0.4,0.5,0.6,0.7,0.8,0.9,1)
  heat.labels <- c("0","0-0.1","0.1-0.2","0.2-0.3","0.3-0.4","0.4-0.5","0.5-0.6","0.6-0.7","0.7-0.8","0.8-0.9","0.9-1")
  heat$value.heat.discrete <- cut(round(heat$value.heat,3), breaks = heat.breaks, labels = heat.labels)
  
  heat.plot <- ggplot(data = heat, aes(x = a.heat, y = m.heat)) + geom_tile(aes(fill = value.heat.discrete)) + labs(x = "a, acculturation", y = "m, migration") + ggtitle(paste("n = ", n, ", L = ", L, ", r = ", r, ", b = ", b, ", c = ", c, ", u = ", u, ", v = ", v, sep = "")) + theme_grey(base_size = 18) + theme(plot.title = element_text(size = 16, hjust = 0.5), legend.text = element_text(size = 10), legend.position = leg, axis.text = element_text(size = 12), axis.line = element_blank()) + scale_fill_viridis(name="p", discrete = TRUE, drop = FALSE) + scale_x_continuous(breaks = seq(0,1,by=0.2), expand = c(0,0)) + scale_y_continuous(expand = c(0,0))
  
  # show on screen
  if (showplot == TRUE) {
    print(heat.plot)
  }
  
  # save to file
  ggsave("model2_rec_am_heatmap.png", width = 20, height = 20, units = "cm", dpi = 300)
  
  list(heat, heat.plot)  # output heat dataframe and plot
}

# create one heatmap here
heat.am <- drawHeatMap.am.model2.rec(m.max = 0.6, t.max = 1000, n = 5, a.max = 1, r = 0.0, b = 1, c = 0.2, u = 0.1, v = 0.5, L = 1, grid.size = 30)

# create 2 heat plots in Fig 5 for different L

r.heats <- 0  # change common values across all heatmaps here
t.max.heats <- 1000
m.max.heats <- 0.6
a.max.heats <- 1
n.heats <- 5
grid.size.heats <- 30
b.heats <- 1
c.heats <- 0.2
u.heats <- 0.1
v.heats <- 0.5

heat1.rec <- drawHeatMap.am.model2.rec(n = n.heats, r = r.heats, t.max = t.max.heats, m.max = m.max.heats, a.max = a.max.heats, b = b.heats,  c = c.heats, u = u.heats, v = v.heats, L = 0, grid.size = grid.size.heats, showplot = F, showlegend = F)[[2]]
heat2.rec <- drawHeatMap.am.model2.rec(n = n.heats, r = r.heats, t.max = t.max.heats, m.max = m.max.heats, a.max = a.max.heats, b = b.heats,  c = c.heats, u = u.heats, v = v.heats, L = 1, grid.size = grid.size.heats, showplot = F, showlegend = T)[[2]]

# use cowplot to create composite plot
legend <- get_legend(heat2.rec)
heat2.rec <- heat2.rec + theme(legend.position="none")
heat1.rec <- heat1.rec + ggtitle("L = 0")
heat2.rec <- heat2.rec + ggtitle("L = 1")
all.heats.rec <- plot_grid(heat1.rec, heat2.rec, legend, labels=c("A", "B"), ncol = 3, nrow = 1, label_size = 24, rel_widths = c(2,2,0.5))

# save to file
save_plot("model2_rec_am_heats_varyL.png", all.heats.rec, base_width = 20, base_height = 10, units = "cm")

# create Figure 6: 3 heatmaps for m & L at v=0.3, v=0.5, v=0.7------------------------

# first create function for a single L-m heatmap

# for a grid.size x grid.size grid (e.g. if grid.size=10, a 10x10 grid)
# a counter tracks progress (can take a while)

drawHeatMap.Lm.model2.rec <- function(m.max, t.max, n, a, r, b, c, u, v, L.max, grid.size, showplot = TRUE, showlegend = TRUE) {
  
  # initialise grid.size x grid.size dataframe 'heat'
  m.heat <- rep(0:(grid.size-1), grid.size) # for now these are 0 to (grid.size-1) to use in loops, convert to 0 to m.max below
  L.heat <- rep(0:(grid.size-1), each = grid.size)
  value.heat <- rep(0,(grid.size^2))
  heat <- data.frame(m.heat, L.heat, value.heat)  
  
  counter <- 0
  
  # add values to heat matrix 
  for (m.loop in 0:(grid.size-1))
  {
    for (L.loop in 0:(grid.size-1))
    {
      p <- Model2.rec(m = m.loop*m.max/(grid.size-1), t.max = t.max, n=n, L = L.loop*L.max/(grid.size-1), r=r, b=b, c=c, u=u, v=v, a=a, showplot = FALSE)
      heat$value.heat[heat$L.heat == L.loop & heat$m.heat == m.loop] <- p[t.max]
      counter <- counter + 1
      cat(counter, "/", grid.size^2, "\n")
    }
  }
  
  # convert L.heat and m.heat values
  heat$L.heat <- heat$L.heat * L.max / (grid.size-1)
  heat$m.heat <- heat$m.heat * m.max / (grid.size-1)
  
  # create heatmap
  
  leg <- "right"
  if (showlegend == FALSE) { leg <- "none" }
  
  # turn Fst values into discrete categories
  heat.breaks <- c(-1,0,0.1,0.2,0.3,0.4,0.5,0.6,0.7,0.8,0.9,1)
  heat.labels <- c("0","0-0.1","0.1-0.2","0.2-0.3","0.3-0.4","0.4-0.5","0.5-0.6","0.6-0.7","0.7-0.8","0.8-0.9","0.9-1")
  heat$value.heat.discrete <- cut(round(heat$value.heat,3), breaks = heat.breaks, labels = heat.labels)
  
  heat.plot <- ggplot(data = heat, aes(x = L.heat, y = m.heat)) + geom_tile(aes(fill = value.heat.discrete)) + labs(x = "L, payoff-biased copying", y = "m, migration") + ggtitle(paste("a = ", a, ", n = ", n, ", r = ", r, ", b = ", b, ", c = ", c, ", u = ", u, ", v = ", v, sep = "")) + theme_grey(base_size = 18) + theme(plot.title = element_text(size = 16, hjust = 0.5), legend.text = element_text(size = 10), legend.position = leg, axis.text = element_text(size = 12), axis.line = element_blank()) + scale_fill_viridis(name="p", discrete = TRUE, drop = FALSE) + scale_x_continuous(breaks = seq(0,1,by=0.2), expand = c(0,0)) + scale_y_continuous(expand = c(0,0))
  
  # show on screen
  if (showplot == TRUE) {
    print(heat.plot)
  }
  
  # save to file
  ggsave("model2_rec_Lm_heatmap.png", width = 20, height = 20, units = "cm", dpi = 300)
  
  list(heat, heat.plot)  # output heat dataframe and plot
}

# create one L-m heatmap here
heat.Lm <- drawHeatMap.Lm.model2.rec(m.max = 0.6, t.max = 1000, n = 5, L.max = 1, r = 0.0, b = 1, c = 0.2, u = 0.1, v = 0.7, a = 0, grid.size = 20)

# create 3 heat plots in Fig 5 for different v

r.heats <- 0  # change common values here
t.max.heats <- 1000
m.max.heats <- 0.6
L.max.heats <- 1
n.heats <- 5
grid.size.heats <- 30
b.heats <- 1
c.heats <- 0.2
u.heats <- 0.1
a.heats <- 0

heat1.rec <- drawHeatMap.Lm.model2.rec(n = n.heats, r = r.heats, t.max = t.max.heats, m.max = m.max.heats, a = a.heats, b = b.heats,  c = c.heats, u = u.heats, v = 0.3, L.max = L.max.heats, grid.size = grid.size.heats, showplot = F, showlegend = F)[[2]]
heat2.rec <- drawHeatMap.Lm.model2.rec(n = n.heats, r = r.heats, t.max = t.max.heats, m.max = m.max.heats, a = a.heats, b = b.heats,  c = c.heats, u = u.heats, v = 0.5, L.max = L.max.heats, grid.size = grid.size.heats, showplot = F, showlegend = F)[[2]]
heat3.rec <- drawHeatMap.Lm.model2.rec(n = n.heats, r = r.heats, t.max = t.max.heats, m.max = m.max.heats, a = a.heats, b = b.heats,  c = c.heats, u = u.heats, v = 0.7, L.max = L.max.heats, grid.size = grid.size.heats, showplot = F, showlegend = T)[[2]]

# use cowplot to create composite plot
legend <- get_legend(heat3.rec)
heat3.rec <- heat3.rec + theme(legend.position="none")
heat1.rec <- heat1.rec + ggtitle("v = 0.3")
heat2.rec <- heat2.rec + ggtitle("v = 0.5")
heat3.rec <- heat3.rec + ggtitle("v = 0.7")
all.heats.rec <- plot_grid(heat1.rec, heat2.rec, heat3.rec, legend, labels=c("A", "B", "C"), ncol = 4, nrow = 1, label_size = 24, rel_widths = c(2,2,2,0.5))

# save to file
save_plot("model2_rec_Lm_heats_varyv.png", all.heats.rec, base_width = 30, base_height = 10, units = "cm")


# add 3 more heatmaps for m & a at v=0.3, v=0.5, v=0.7------------------------

# run the preceding code first to create 'all.heats.rec' showing m&L for v=0.3,0.5 and 0.7

# now create three more heatmaps for a-m at different values of v

# for a grid.size x grid.size grid (e.g. if grid.size=10, a 10x10 grid)
# a counter tracks progress (can take a while)

drawHeatMap.am.model2.rec <- function(m.max, t.max, n, L, r, b, c, u, v, a.max, grid.size, showplot = TRUE, showlegend = TRUE) {
  
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
      p <- Model2.rec(m = m.loop*m.max/(grid.size-1), t.max = t.max, n=n, a = a.loop*a.max/(grid.size-1), r=r, b=b, c=c, u=u, v=v, L=L, showplot = FALSE)
      heat$value.heat[heat$a.heat == a.loop & heat$m.heat == m.loop] <- p[t.max]
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
  
  # turn Fst values into discrete categories
  heat.breaks <- c(-1,0,0.1,0.2,0.3,0.4,0.5,0.6,0.7,0.8,0.9,1)
  heat.labels <- c("0","0-0.1","0.1-0.2","0.2-0.3","0.3-0.4","0.4-0.5","0.5-0.6","0.6-0.7","0.7-0.8","0.8-0.9","0.9-1")
  heat$value.heat.discrete <- cut(round(heat$value.heat,3), breaks = heat.breaks, labels = heat.labels)
  
  heat.plot <- ggplot(data = heat, aes(x = a.heat, y = m.heat)) + geom_tile(aes(fill = value.heat.discrete)) + labs(x = "a, acculturation", y = "m, migration") + ggtitle(paste("L = ", L, ", n = ", n, ", r = ", r, ", b = ", b, ", c = ", c, ", u = ", u, ", v = ", v, sep = "")) + theme_grey(base_size = 18) + theme(plot.title = element_text(size = 16, hjust = 0.5), legend.text = element_text(size = 10), legend.position = leg, axis.text = element_text(size = 12), axis.line = element_blank()) + scale_fill_viridis(name="p", discrete = TRUE, drop = FALSE) + scale_x_continuous(breaks = seq(0,1,by=0.2), expand = c(0,0)) + scale_y_continuous(expand = c(0,0))
  
  # show on screen
  if (showplot == TRUE) {
    print(heat.plot)
  }
  
  # save to file
  ggsave("model2_rec_am_heatmap.png", width = 20, height = 20, units = "cm", dpi = 300)
  
  list(heat, heat.plot)  # output heat dataframe and plot
}

# create one a-m heatmap here
heat.am <- drawHeatMap.am.model2.rec(m.max = 0.6, t.max = 1000, n = 5, a.max = 1, r = 0.0, b = 1, c = 0.2, u = 0.1, v = 0.7, L = 0, grid.size = 20)

# create 3 heat plots to add to Fig 6 for different v

r.heats <- 0  # change common values here
t.max.heats <- 1000
m.max.heats <- 0.6
a.max.heats <- 1
n.heats <- 5
grid.size.heats <- 30
b.heats <- 1
c.heats <- 0.2
u.heats <- 0.1
L.heats <- 0

heat1a.rec <- drawHeatMap.am.model2.rec(n = n.heats, r = r.heats, t.max = t.max.heats, m.max = m.max.heats, L = L.heats, b = b.heats,  c = c.heats, u = u.heats, v = 0.3, a.max = a.max.heats, grid.size = grid.size.heats, showplot = F, showlegend = F)[[2]]
heat2a.rec <- drawHeatMap.am.model2.rec(n = n.heats, r = r.heats, t.max = t.max.heats, m.max = m.max.heats, L = L.heats, b = b.heats,  c = c.heats, u = u.heats, v = 0.5, a.max = a.max.heats, grid.size = grid.size.heats, showplot = F, showlegend = F)[[2]]
heat3a.rec <- drawHeatMap.am.model2.rec(n = n.heats, r = r.heats, t.max = t.max.heats, m.max = m.max.heats, L = L.heats, b = b.heats,  c = c.heats, u = u.heats, v = 0.7, a.max = a.max.heats, grid.size = grid.size.heats, showplot = F, showlegend = T)[[2]]

# use cowplot to create composite plot
legend <- get_legend(heat3a.rec)
heat3a.rec <- heat3a.rec + theme(legend.position="none")
heat1a.rec <- heat1a.rec + ggtitle("v = 0.3")
heat2a.rec <- heat2a.rec + ggtitle("v = 0.5")
heat3a.rec <- heat3a.rec + ggtitle("v = 0.7")
all.heats.a.rec <- plot_grid(heat1a.rec, heat2a.rec, heat3a.rec, labels=c("D", "E", "F"), ncol = 4, nrow = 1, label_size = 24, rel_widths = c(2,2,2,0.5))


# combine the two sets of heatmaps
all.heats.a.L.rec <- plot_grid(all.heats.rec, all.heats.a.rec, nrow = 2)

# save to file
save_plot("model2_rec_Lam_heats_varyv.png", all.heats.a.L.rec, base_width = 30, base_height = 20, units = "cm")

