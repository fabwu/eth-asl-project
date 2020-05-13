
library(ggplot2)
library(svglite)
library(scales)
## library(viridis)

# read files
## df = read.csv("3c.csv")

# roof functions
bandwidth = function(x) 8*x
perf1 = function(x) 4
perf2 = function(x) 16 

# points

x = seq(0,4,.1)
df = data.frame(x = x)


# plot
y_breaks = c(1/32,1/16,1/8,1/4,1/2,1,2,4,8,16,32)
image = ggplot(data = df, mapping = aes(x = x)) +
    ggtitle("Roofline Plot") +
    stat_function(fun = bandwidth, colour="orange", lwd=2) +
    geom_label(aes(x=0.8, y=29, label = "P <= BI"), colour = "orange", hjust=-1/8, lwd=6) +
    stat_function(fun = perf1, colour="blue", lwd=2) +
    geom_label(aes(x=1, y=5, label = "P <= pi_1 = 5"), colour = "blue", hjust=-1/8, lwd=6) +
    stat_function(fun = perf2, colour="red", lwd=2) +
    geom_label(aes(x=1, y=20, label = "P <= pi_1 = 20"), colour = "red", hjust=-1/8, lwd=6) +
    geom_point(aes(x=10/32, y=3), colour="blue", fill= "white", shape=24, lwd=4) +
    geom_label(aes(x=10/32, y=3, label = "computation1 (3b)"),
               label.size = NA,
               colour = "blue", hjust=-1/16, vjust=1, lwd=6) +
    geom_point(aes(x=10/32, y=2), colour="blue", fill="white", shape=23, lwd=4) +
    geom_label(aes(x=10/32, y=2, label = "computation2 (3b)"),
               label.size = NA,
               colour = "blue", hjust=-1/16, vjust=1,  lwd=6) +
    geom_point(aes(x=10/32, y=5), color="blue", fill="white", shape=25, lwd=4) +
    geom_label(aes(x=10/32, y=5, label = "computation3 (3b)"),
               label.size = NA,
               color = "blue", hjust=-1/16, vjust=1.3,  lwd=6) +
    geom_point(aes(x=10/32, y=10), color="red", fill="white", shape=25, lwd=4) +
    geom_point(aes(x=10/32, y=10), color="red", fill="white", shape=24, lwd=4) +
    geom_point(aes(x=10/32, y=10), color="red", shape=25, lwd=4) +
    geom_label(aes(x=10/32, y=10, label = "computation1 & computation3 (3c, SIMD)"),
               label.size = NA,
               color = "red", hjust=-1/32, vjust=1, lwd=6) +
    geom_point(aes(x=10/32, y=8),
               color="red", fill="white", shape=23, lwd=4) +
    geom_label(aes(x=10/32, y=8, label = "computation2 (3c, SIMD)"),
               label.size = NA,
               color = "red", vjust=1, hjust=-1/16, lwd=6) +
    geom_point(aes(x=5/64, y=2), colour="purple", fill="white", shape=23, lwd=4) +
    geom_label(aes(x=5/64, y=2, label = "computation2 (3d)"),
               label.size = NA,
               color = "purple", vjust=1, hjust=-1/16, lwd=6) +
    ## geom_point(aes(x=5/64, y=2.5), colour="blue", fill="white", shape=23, lwd=4) +
    geom_point(aes(x=5/64, y=2.5), color="purple", fill="white", shape=25, lwd=4) +
    geom_point(aes(x=5/64, y=2.5), color="purple", fill="white", shape=24, lwd=4) +
    geom_point(aes(x=5/64, y=2.5), color="purple", shape=25, lwd=4) +
    geom_label(aes(x=5/64, y=2.5, label = "computation1 & computation3 (3d)"),
               label.size = NA,
               color = "purple", vjust=1, hjust=-1.5/32, lwd=6) +
    ## geom_label(aes(x=10/32, y=2, label = "computation2"),
    ##            label.size = NA,
    ##            colour = "purple", hjust=-1/16, vjust=1,  lwd=6) +
    scale_x_continuous(name=("I = W/Q [flops/byte]"),
                       trans='log2',breaks=y_breaks, labels=trans_format('log2', math_format(2^.x)), limits=c(1/32,2)) +
    scale_y_continuous(name=("P = W/T [flops/cycle]"),
                       trans='log2',breaks=y_breaks,labels=y_breaks, limits=c(1,32)) +
    theme_minimal() +
    theme(
        text = element_text(size=16),
        plot.title = element_text(size=24),
        ## axis.title.y = element_text(angle=0),
        axis.text=element_text(size=16),
        ## axis.title = element_text(size=24),
        legend.position = "none",
        legend.title = element_blank())
    ## xlim(1/8,32) +

image


## breaks = seq(0,23,4)
## minor_breaks = seq(0,23,1)
## labels = c("0",expression(2^{4}),expression(2^{8}), expression(2^{12}), expression(2^{16}), expression(2^{20}))


## # create plot
## theme_set(theme_light(base_size = 24))
## image = ggplot(data=df, aes(x=n)) +
##   ggtitle("Performance") +
##   geom_line(aes(y=mean, color="3c"), lwd=2) +
##   geom_point(aes(y=mean, color="3c"), lwd=4) +
##   geom_label(aes(x = 12, y = 7, label = "-O3 -march=native", color = "3c"), lwd=8) +
##   scale_x_continuous(breaks=breaks, labels=labels, expand = c(0, 2), minor_breaks=minor_breaks) +
##   scale_y_continuous(name=("[flops/cycle]"), limits=c(0,8)) +
##   theme(legend.position = c(0, 1),legend.justification = c(0, 1)) +
##   theme(
##         axis.title.y = element_text(angle=0),
##         legend.position = "none",
##         legend.title = element_blank())

## # show image
# image

## # save as svg
#ggsave(file="3-plot.pdf", plot=image, width=16, height=10)
