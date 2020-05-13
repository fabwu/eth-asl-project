library(ggplot2)
library(svglite)

## executables
executables = c(
    "ref_impl",
    "baseline",
    "01_precompute_indices"
)

## image sizes
## TODO: only read the csv defined here
image_sizes = c(
    64,
    128,
    256,
    512,
    1024
)

## collect data
df = data.frame()
for (exe in executables) {

    ## read files
    df64 = read.csv(paste("data/", exe, "/lena_64.csv", sep=""), sep=';')
    df128 = read.csv(paste("data/", exe, "/monkey_128.csv", sep=""), sep=';')
    df256 = read.csv(paste("data/", exe, "/lena_256.csv", sep=""), sep=';')
    df512 = read.csv(paste("data/", exe, "/lena_512.csv", sep=""), sep=';')
    df1024 = read.csv(paste("data/", exe, "/grey-parrot_1024.csv", sep=""), sep=';')

    df_temp = data.frame(
        rep(exe, length(image_sizes)),
        image_sizes,
        c(
            mean(df64$flops.cycle),
            mean(df128$flops.cycle),
            mean(df256$flops.cycle),
            mean(df512$flops.cycle),
            mean(df1024$flops.cycle)
        )
    )
    df = rbind(df, df_temp)

}
colnames(df) = c("executable", "size", "performance")


## Create plot
theme_set(theme_light(base_size = 24))
image = ggplot(data=df,
               aes(x=size,
                   y=performance,
                   group=executable,
                   color=executable)) +
    ## ggtitle("Performance") +
    geom_line(lwd=2) +
    geom_point(lwd=4) +
    scale_x_continuous(
        name="Image Size (n x n)",
        expand = c(0, 20)
    ) +
    scale_y_continuous(name=("[flops/cycle]"), limits=c(0,4)) +
    theme(
        axis.title.y = element_text(angle=0),
        legend.position = "bottom",
        legend.title = element_blank())

## save plot
ggsave(file="performance.png",
       dpi = 300,
       plot=image,
       width=15,
       height=8)
