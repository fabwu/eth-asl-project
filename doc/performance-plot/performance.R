library(ggplot2)
library(svglite)

## read files
df64 = read.csv("data/lena_64.csv", sep=';')
df128 = read.csv("data/monkey_128.csv", sep=';')
df256 = read.csv("data/lena_256.csv", sep=';')
## df512 = read.csv("data/lena_512.csv", sep=';')

## create dataframe
df = data.frame(
    c(64,
      128,
      256),
    c(mean(df64$flops.cycle),
      mean(df128$flops.cycle),
      mean(df256$flops.cycle))
)
colnames(df) = c("size", "performance")

## colnames(df) = c("n", "mean")

## breaks = seq(64,128,256)
## minor_breaks = seq(0,23,1)
## labels = c("0",expression(2^{4}),expression(2^{8}), expression(2^{12}), expression(2^{16}), expression(2^{20}))


## Create plot
theme_set(theme_light(base_size = 24))
image = ggplot(data=df, aes(x=size)) +
    ggtitle("Performance") +
    geom_line(aes(y=performance), lwd=1) +
    geom_point(aes(y=performance), lwd=2) +
    scale_x_continuous(
        name="Input Size",
        expand = c(0, 20),
        ) +
    scale_y_continuous(name=("[flops/cycle]"), limits=c(0,1)) +
    theme(
        axis.title.y = element_text(angle=0),
        legend.position = "none",
        legend.title = element_blank())

## old homework plot
##     geom_line(aes(y=mean, color="3c"), lwd=2) +
##     geom_point(aes(y=mean, color="3c"), lwd=4) +
##     geom_label(aes(x = 12, y = 7, label = "-O3 -march=native", color = "3c"), lwd=8) +
##     scale_x_continuous(breaks=breaks, labels=labels, expand = c(0, 2), minor_breaks=minor_breaks) +
##     scale_y_continuous(name=("[flops/cycle]"), limits=c(0,8)) +
##     theme(legend.position = c(0, 1),legend.justification = c(0, 1)) +
##     theme(
##         axis.title.y = element_text(angle=0),
##         legend.position = "none",
##         legend.title = element_blank())

## display image
image

## save as svg
ggsave(file="performance.pdf", plot=image, width=16, height=10)
