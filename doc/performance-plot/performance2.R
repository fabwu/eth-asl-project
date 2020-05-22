library(ggplot2)
library(svglite)

## executables

    ## "baseline",
    ## "01_precompute_indices",
    ## "02_precomputations",
    ## "baseline",
executables = c(
    "25_ilp",
    "30_simd_precomp_rotations",
    "31_simd_precomp_rotations_no_bac_simd",
    "35_simd",
    "36_pascal"
    ## "ref_impl",
    ## "01_precompute_indices",
    ## "05_rtd_no_precompute",
    ## "08_remove_function_and_add_fma",
    ## "09_basel",
    ## "10_simd_brightness_and_contrast",
    ## "11_simd_rtd"
)

## for exe in baseline 01_precompute_indices

## image sizes
## TODO: only read the csv defined here
image_sizes = c(
    64,
    128,
    256,
    512,
    1024,
    2048,
    4096
)

## collect data
df = data.frame()
for (exe in executables) {

    ## read files
    files = c(
        paste("data/", exe, "/lena_64.csv", sep=""),
        paste("data/", exe, "/monkey_128.csv", sep=""),
        paste("data/", exe, "/lena_256.csv", sep=""),
        paste("data/", exe, "/lena_512.csv", sep=""),
        paste("data/", exe, "/grey-parrot_1024.csv", sep=""),
        paste("data/", exe, "/matterhorn_2048.csv", sep=""),
        paste("data/", exe, "/lion_4096.csv", sep="")
    )

    for(i in 1:7){
        if(file.exists(files[i])){
            csv = read.csv(files[i], sep=';')
            df_temp = data.frame(
                exe,
                image_sizes[i],
                mean(csv$flops.cycle),
                mean(csv$cycles),
                mean(csv$flops)
            )
            df = rbind(df, df_temp)
        }
    }

    ## df64 = read.csv(, sep=';')
    ## df128 = read.csv(, sep=';')
    ## df256 = read.csv(, sep=';')
    ## df512 = read.csv, sep=';')
    ## df1024 = read.csv, sep=';')
    ## df2048 = read.csv(, sep=';')

    ## df_temp = data.frame(
    ##     rep(exe, length(image_sizes)),
    ##     image_sizes,
    ##     c(
    ##         mean(df64$flops.cycle),
    ##         mean(df128$flops.cycle),
    ##         mean(df256$flops.cycle),
    ##         mean(df512$flops.cycle),
    ##         mean(df1024$flops.cycle),
    ##         mean(df2048$flops.cycle)
    ##     ),
    ##     c(
    ##         mean(df64$cycles),
    ##         mean(df128$cycles),
    ##         mean(df256$cycles),
    ##         mean(df512$cycles),
    ##         mean(df1024$cycles),
    ##         mean(df2048$cycles)
    ##     ),
    ##     c(
    ##         mean(df64$flops),
    ##         mean(df128$flops),
    ##         mean(df256$flops),
    ##         mean(df512$flops),
    ##         mean(df1024$flops),
    ##         mean(df2048$flops)
    ##     )
    ## )

}
colnames(df) = c("executable", "size", "performance", "cycles", "flops")


# baseplot
theme_set(theme_light(base_size = 24))
baseplot = ggplot(data=df,
                  aes(x=size,
                      group=executable,
                      color=executable)) +
    scale_x_continuous(
        name="Image Size (n x n)",
        expand = c(0, 20)
    ) +
    theme(
        axis.title.y = element_text(angle=0),
        legend.position = "bottom",
        legend.title = element_blank())


## Performance plot
performance_plot = baseplot +
    geom_line(aes(y=performance), lwd=2) +
    geom_point(aes(y=performance), lwd=4) +
    ggtitle("Performance") +
    scale_y_continuous(name=("[flops/cycle]"), limits=c(0,4))
ggsave(file="performance.png",
       dpi = 300,
       plot=performance_plot,
       width=15,
       height=8)

## Runtime plot
runtime_plot = baseplot +
    ggtitle("Runtime") +
    geom_line(aes(y=cycles), lwd=2) +
    geom_point(aes(y=cycles), lwd=4) +
    scale_y_continuous(name=("[cycles]"))

## Flops plot
flops_plot = baseplot +
    ggtitle("Flops") +
    geom_line(aes(y=flops), lwd=2) +
    geom_point(aes(y=flops), lwd=4) +
    scale_y_continuous(name=("[flops]"))

## Save plots
ggsave(file="performance.png",
       dpi = 300,
       plot=performance_plot,
       width=15,
       height=8)

ggsave(file="runtime.png",
       dpi = 300,
       plot=runtime_plot,
       width=15,
       height=8)

ggsave(file="flops.png",
       dpi = 300,
       plot=flops_plot,
       width=15,
       height=8)
