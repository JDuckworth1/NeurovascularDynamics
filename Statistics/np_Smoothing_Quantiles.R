# Loading data from csv file, plotting and performing kernel smoothing
# relevant functions from R np package: npcdistbw, npqreg.
# documentation: https://www.jstatsoft.org/article/view/v027i05
# https://cran.r-project.org/web/packages/np/np.pdf

library(ggplot2)
library(readr)
library(KernSmooth)
library(ks)
library(np)

saveloc = "PathToSaveLocation"
loadloc = "PathToDataLocation"
setwd(loadloc)

# LOAD AND VISUALIZE DATA
savestr <- "xyScatter_DataFileName"

data <- read_csv(paste(savestr,".csv",sep = ""))
setwd(saveloc)


# Filter Data as necessary:
data <- data[!is.nan(data[[2]]), ] # Removing NaN rows
data_plot <- data

dat_mat <- data.matrix(data) # Data matrix instead of data.frame
rawdata_plot <- ggplot(dat_mat, mapping = aes(x = xData, y = yData)) +
  geom_point(color = "blue") + 
  theme_minimal()
rawdata_plot

# Sort x data
# data <- data[order(data$xData), ]

# NO TRANSFORMATION #############
# Determine the BW
bw0 <- npcdistbw(formula = yData ~ xData,
                 data = data,
                 regtype = 'll',
                 bwmethod = "cv.ls",
                 bwtype = "adaptive_nn",
                 # bwtype = "generalized_nn",
                 # bwtype = "fixed",
                 cxkertype = "epanechnikov",
                 cxkerorder = 2,
                 cykertype = "epanechnikov",
                 cykerorder = 2
                 # bwtype = "fixed",  # Fixed bandwidth
                 # xbw = c(0.01),  # Fixed bandwidth for x
                 # ybw = c(0.01),  # Fixed bandwidth for y
                 # bandwidth.compute = FALSE,  # Don't compute bandwidth
                 # cxkertype = "uniform",  # Uniform kernel for x
                 # cykertype = "uniform"   # Uniform kernel for y
)
print(bw0)


# Perform quantile regression
model.q0.25 <- npqreg(bws=bw0, tau=0.25)
model.q0.50 <- npqreg(bws=bw0, tau=0.50)
model.q0.75 <- npqreg(bws=bw0, tau=0.75)

p<-ggplot(data_plot, aes(x = xData,y = yData)) +
  geom_point() +
  theme(
    panel.grid.major = element_line(color = "gray", linewidth = 0.5),
    panel.grid.minor = element_line(color = "lightgray", linewidth = 0.3)
  ) +
  xlim(-0.2, 0.75) +
  ylim(0, 1) +  
  theme(panel.grid = element_blank()) +  # Remove grid
  theme(panel.background = element_blank(), axis.line = element_line(color = "black"))# + # Keep only x and y axis

p + scale_x_continuous(breaks = c(-0.2, -0.15,-0.1, -0.05, 0, 0.05, 0.1, 0.15, 0.2, 0.25, 0.3, 0.35, 0.4, 0.45, 0.5, 0.55, 0.6, 0.65, 0.7),minor_breaks = seq(0))  # X-axis ticks every 10 units

data$q25 <- model.q0.25$quantile
data$q50 <- model.q0.50$quantile
data$q75 <- model.q0.75$quantile
p <- p + geom_line(data = data, aes(x = xData, y = q25), color = "red", linewidth = 0.5)
p <- p + geom_line(data = data, aes(x = xData, y = q50), color = "red", linewidth = 0.5)
p <- p + geom_line(data = data, aes(x = xData, y = q75), color = "red", linewidth = 0.5)
print(p)
ggsave(paste(savestr,"adaptive_nn.eps",sep = ""), plot = p, device = "eps", width = 6, height = 4, dpi = 300)
ggsave(paste(savestr,"adaptive_nn.png",sep = ""), plot = p, device = "png", width = 6, height = 4, dpi = 300)








