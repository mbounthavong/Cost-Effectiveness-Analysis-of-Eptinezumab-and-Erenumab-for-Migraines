#----------------------------#
#### Formatting Functions ####
#----------------------------#
gg_color_hue <- function(n) {
  hues = seq(15, 375, length = n + 1)
  hcl(h = hues, l = 65, c = 100)[1:n]
}

number_ticks <- function(n) {function(limits) pretty(limits, n)} #Function for number of ticks in ggplot


#----------------------------#
#### Scatter CEA frontier ####
#----------------------------#
ScatterCE <- function(strategies, m.c, m.e){
  library(reshape2)
  library(ellipse)
  library(plyr)
  library(ggplot2)
  library(scales)
  df.c  <- melt(m.c, variable.name="Strategy")
  levels(df.c$Strategy) <- strategies
  df.e <- melt(m.e, variable.name="Strategy")
  levels(df.e$Strategy) <- strategies
  CE <- cbind(df.c, df.e[,2])
  colnames(CE) <- c("Strategy","Cost","Effectiveness")
  
  # Ellipses code
  df_ell <- data.frame() #create an empty dataframe
  # for each level in df$groups 
  for(g in levels(CE$Strategy)){
    # create 100 points per variable around the mean of each group
    df_ell <- rbind(df_ell, 
                    cbind(as.data.frame(with(CE[CE$Strategy==g,], 
                                             ellipse::ellipse(cor(Effectiveness, Cost), 
                                                     scale=c(sd(Effectiveness),sd(Cost)), 
                                                     centre=c(mean(Effectiveness),mean(Cost)))
                    )),group=g))
  }
  Means <- ddply(CE,.(Strategy),summarise,
                 N = length(Cost),
                 Cost.mean = mean(Cost),
                 Eff.mean = mean(Effectiveness))  
  #Define ggplot object
  print(
  ggplot(CE, aes(x=Effectiveness, y=Cost, color=Strategy)) + 
    geom_point(size=0.7, alpha = 0.2, shape = 21) +
    geom_point(data = Means,aes(x = Eff.mean, y = Cost.mean, shape=Strategy),
               size=8,fill="white") +
    geom_text(data = Means,aes(x = Eff.mean, y = Cost.mean, label=c(1,2,3)),size=5,colour="gray",alpha=1) +
    geom_path(data=df_ell, aes(x=x, y=y,colour=group), size=1, linetype=2, alpha=1) + # draw ellipse lines
    ggtitle("Cost-Effectiveness Scatterplot") +
    scale_colour_discrete(l=50) +  # Use a slightly darker palette than normal
    scale_y_continuous(labels = dollar)+
    scale_x_continuous(breaks=number_ticks(6))+
    theme_bw() +
    theme(legend.position="bottom")
  )
}

#------------#
#### CEAF ####
#------------#
ceaf <- function(v.wtp, strategies, m.e, m.c, currency = "$"){
  library(reshape2)
  library(ggplot2)
  library(scales)
  n.sim <- nrow(m.e)
  n.str  <- ncol(m.e)
  # Matrix to store indicator of CEAC
  m.cea  <- array(0, dim = c(length(v.wtp), n.str))
  # Vector to store indicator of strategy at CEAF
  v.ceaf <- array(0, dim = c(length(v.wtp), 1))
  
  for(l in 1:length(v.wtp)){
    m.nmb <-  v.wtp[l]*m.e - m.c # Effectiveness minus Costs
    # Calculate point of CEAF, i.e., the strategy with the highest expected NMB
    v.ceaf[l, 1] <- which.max(colMeans(m.nmb))
    # Calculate points in CEAC, i.e, the probability that each strategy is cost-effective
    max.nmb <- max.col(m.nmb)
    opt <- table(max.nmb)
    m.cea[l, as.numeric(names(opt))] <- opt/n.sim
  }
  m.ceaf <- m.cea[cbind(1:length(v.wtp), v.ceaf)]
  
  df.cea <- data.frame(cbind(v.wtp, m.cea, m.ceaf))
  colnames(df.cea) <- c("WTP", strategies, "Frontier")
  
  df.ceac <- melt(df.cea, 
                  id.vars = "WTP", 
                  variable.name = "Strategy") 
  
  ## Plot CEAC & CEAF
  # Format to plot frontier
  strats <- 1:(length(unique(df.ceac$Strategy))-1)
  point.shapes <- c(strats+14, 0) # Shapes: http://sape.inf.usi.ch/quick-reference/ggplot2/shape
  colors <- c(gg_color_hue(length(strats)), "#696969")
  point.size <- c(rep(2, length(strats)), 4) # Trick consists on firts define size as aes then manually change it
  # Plot CEAC & CEAF
  print(
    ggplot(data = df.ceac, aes(x = WTP/1000, y = value)) +
      geom_point(aes(shape = Strategy, color = Strategy, size = Strategy)) +
      geom_line(aes(color = Strategy)) +
      ggtitle("Cost-Effectiveness Acceptability Curves and Frontier") + 
      # scale_colour_hue(l=50, values=colors) +
      scale_x_continuous(breaks=number_ticks(20))+
      scale_shape_manual(values=point.shapes) +
      #scale_shape(solid=TRUE) +
      scale_color_manual(values=colors) + 
      scale_size_manual(values = point.size) +
      #scale_alpha_manual(values=c(rep(0, length(strats)), 0.5)) + 
      xlab(paste("Willingness to pay (Thousand ", currency,"/QALY)", sep = "")) +
      ylab("Probability of Cost-Effectiveness") +
      theme_bw(base_size = 14) +
      theme(legend.position = "bottom")
  )
  
}


#----------------------------#
#### Expected Loss Curves ####
#----------------------------#
elc <- function(v.wtp, strategies, m.e, m.c){
  library(reshape2)
  library(ggplot2)
  library(scales)
  n.sim <- nrow(m.e)
  n.str  <- ncol(m.e)
  m.exp.loss <- matrix(0, nrow = length(v.wtp), ncol = n.str)
  for(l in 1:length(v.wtp)){
    m.nmb <- m.e*v.wtp[l] - m.c # Effectiveness minus Costs, with vector indexing
    max.str <- max.col(m.nmb)
    m.loss <- m.nmb - m.nmb[cbind(1:n.sim, max.str)]
    m.exp.loss[l, ] <- colMeans(m.loss)
  }
  # Optimal strategy based on lowest expected loss
  optimal.str <- max.col(m.exp.loss)
  # Expected loss of optimal strategy
  optimal.el <- m.exp.loss[cbind(1:length(v.wtp), optimal.str)]
  # Format expected loss for plotting
  df.exp.loss <- data.frame(cbind(v.wtp, m.exp.loss, optimal.el))
  colnames(df.exp.loss) <- c("WTP", strategies, "Frontier & EVPI")
  df.exp.loss.plot <- melt(df.exp.loss, 
                           id.vars = "WTP", 
                           variable.name = "Strategy")
  ## Plot expected losses
  # Format to plot frontier
  strats <- 1:(length(unique(df.exp.loss.plot$Strategy))-1)
  point.shapes <- c(strats+14, 0) # Shapes: http://sape.inf.usi.ch/quick-reference/ggplot2/shape
  colors <- c(gg_color_hue(length(strats)), "#696969")
  point.size <- c(rep(2, length(strats)), 4) # Trick consists on firts define size as aes then manually change it
  # Plot ELC
  print(
    ggplot(data = df.exp.loss.plot, aes(x = WTP/1000, y = -value)) +
      geom_point(aes(shape = Strategy, color = Strategy, size = Strategy)) +
      geom_line(aes(color = Strategy)) +
      ggtitle("Expected Loss Curves") + 
      #scale_colour_hue(l=50, values=colors) +
      scale_x_continuous(breaks=number_ticks(20))+
      scale_y_continuous(labels = comma, breaks = number_ticks(8)) +
      xlab("Willingness to Pay (Thousand $/QALY)") +
      ylab("Expected loss ($)") +
      scale_shape_manual(values=point.shapes) +
      scale_color_manual(values=colors) + 
      scale_size_manual(values = point.size) +
      theme_bw(base_size = 14) +
      theme(legend.position = "bottom")
  )
}

#----------------------------------------------------#
#### Expected Value of Perfect Information (EVPI) ####
#----------------------------------------------------#
evpi <- function(v.wtp, m.e, m.c){
  # Load required packages
  require(reshape2)
  require(matrixStats)
  # Create scalar with number of simulations
  n.sim <- nrow(m.e)
  # Create scalar with number of strategies (i.e. number of columns of m.e)
  n.str <- nrow(m.e)
  # Data frame to store EVPI for each WTP threshold
  evpi <- as.data.frame(array(0, dim = c(length(v.wtp), 2)))
  # Name data frame's columns
  colnames(evpi) <- c("WTP", "EVPI")
  # Assign vector with WTP thresholds to first column of `evpi`
  evpi$WTP <- v.wtp
  # Estimate the Loss matrix and EVPI at each WTP threshold
  for(l in 1:length(v.wtp)){
    # Compute NMB with vector indexing
    nmb <-  m.e*v.wtp[l]-m.c
    ## Find the optimal strategy with current info
    d.star <- which.max(colMeans(nmb))
    ## Calculate the opportunity loss from choosing d.star for each strategy
    m.loss <- nmb - nmb[, d.star]
    ## Compute EVPI
    evpi$EVPI[l] <- mean(rowMaxs(as.matrix(m.loss))) # needs to be a numeric matrix
  }
  # Load required packages
  require(ggplot2) 
  require(scales) 
  # Plot EVPI
  ggplot(data = evpi, aes(x = WTP/1000, y = EVPI)) +
    #geom_point() +
    geom_line() +  
    ggtitle("Expected Value of Perfect Information") +
    scale_x_continuous(labels = comma, breaks = number_ticks(20))+ 
    scale_y_continuous(labels = comma, breaks = number_ticks(6))+
    xlab("Willingness to Pay (Thousand $/QALY)") +
    ylab("EVPI ($)") +
    theme_bw(base_size = 14) +
    theme(legend.position="bottom")
}
