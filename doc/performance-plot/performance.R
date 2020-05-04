library(ggplot2)
library(svglite)

## read files
df64 = read.csv("data/lena_64.csv", sep=';')
df128 = read.csv("data/monkey_128.csv", sep=';')
df256 = read.csv("data/lena_256.csv", sep=';')
df512 = read.csv("data/lena_512.csv", sep=';')
df1024 = read.csv("data/grey-parrot_1024.csv", sep=';')


## create dataframe
df = data.frame(
    c(64,
      128,
      256,
      512,
      1024),
    c(mean(df64$flops.cycle),
      mean(df128$flops.cycle),
      mean(df256$flops.cycle),
      mean(df512$flops.cycle),
      mean(df1024$flops.cycle))
)
colnames(df) = c("size", "performance")


## Create plot
theme_set(theme_light(base_size = 24))
image = ggplot(data=df, aes(x=size)) +
    ## ggtitle("Performance") +
    geom_line(aes(y=performance), lwd=2) +
    geom_point(aes(y=performance), lwd=4) +
    scale_x_continuous(
        name="Image Size (n x n)",
        expand = c(0, 20),
        ) +
    scale_y_continuous(name=("[flops/cycle]"), limits=c(0,4)) +
    theme(
        axis.title.y = element_text(angle=0),
        legend.position = "none",
        legend.title = element_blank())

## save plot
ggsave(file="performance.png",
       dpi = 300,
       plot=image,
       width=15,
       height=8)
