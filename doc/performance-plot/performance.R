library(ggplot2)
library(svglite)

executables = c(
    ## "0_baseline",
    ## "25_ilp",
    ## "35_simd"
    "0_baseline",
    "31_simd_precomp_rotations_no_bac_simd",
    "40_ilp_norot_90_270",
    "41_simd_norot_90_270"
)

## one
    ## "0_baseline",
    ## "25_ilp",
    ## "35_simd"

## two
    ## "0_baseline",
    ## "31_simd_precomp_rotations_no_bac_simd",
    ## "40_ilp_norot_90_270",
    ## "41_simd_norot_90_270"


## "35_simd",
## "31_simd_precomp_rotations_no_bac_simd",
## "40_ilp_norot_90_270",
## "41_simd_norot_90_270"

## 31=better-memory-access
## 40=ILP_without_90/270
## 41=SIMD_without_90/270

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
        ## new images
        ## paste("data/", exe, "/lion_64.csv", sep=""),
        ## paste("data/", exe, "/lion_128.csv", sep=""),
        ## paste("data/", exe, "/lion_256.csv", sep=""),
        ## paste("data/", exe, "/lion_512.csv", sep=""),
        ## paste("data/", exe, "/lion_1024.csv", sep=""),
        ## paste("data/", exe, "/lion_2048.csv", sep=""),
        ## paste("data/", exe, "/lion_4096.csv", sep="")

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
        axis.text.x = element_text(size=24),
        axis.text.y = element_text(size=24),
        legend.position = "none",
        legend.title = element_blank())


## Performance plot
performance_plot = baseplot +
    geom_line(aes(y=performance), lwd=2) +
    geom_point(aes(y=performance), lwd=4) +
    ggtitle("i7-8650U @ 1.9 GHz") +
    scale_y_continuous(name=("[flops/cycle]"), limits=c(0,4)) +
    ## geom_label(aes(x=1500, y=1.7, label = "Scalar"),
    ##            label.size = NA,
    ##            color = "#00B937", hjust=-1/32, vjust=1, lwd=8) +
    ## geom_label(aes(x=1500, y=2.4, label = "SIMD"),
    ##            label.size = NA,
    ##            color = "#609BFE", hjust=-1/32, vjust=1, lwd=8) +
    ## geom_label(aes(x=1500, y=0.8, label = "Baseline"),
    ##            label.size = NA,
    ##            color = "#F7756C", hjust=-1/32, vjust=1, lwd=8)
    geom_label(aes(x=1300, y=1.95, label = "Copy Memory"),
               label.size = NA,
               color = "#7BAD00", hjust=-1/32, vjust=1, lwd=8) +
    geom_label(aes(x=900, y=3.7, label = "SIMD (without 90/270)"),
               label.size = NA,
               color = "#C67BFE", hjust=-1/32, vjust=1, lwd=8) +
    geom_label(aes(x=900, y=1.3, label = "Scalar (without 90/270)"),
               label.size = NA,
               color = "#00BEC3", hjust=-1/32, vjust=1, lwd=8) +
    geom_label(aes(x=1000, y=0.45, label = "Baseline"),
               label.size = NA,
               color = "#F7756C", hjust=-1/32, vjust=1, lwd=8)
## ggsave(file="performance.png",
##        dpi = 300,
##        plot=performance_plot,
##        width=15,
##        height=8)

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
ggsave(file="performance.pdf",
       dpi = 300,
       plot=performance_plot,
       width=15,
       height=8)

ggsave(file="runtime.pdf",
       dpi = 300,
       plot=runtime_plot,
       width=15,
       height=8)

ggsave(file="flops.pdf",
       dpi = 300,
       plot=flops_plot,
       width=15,
       height=8)
