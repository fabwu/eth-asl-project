library(ggplot2)
library(ggthemes)
library(svglite)
library(latex2exp)

executables = c(
  ## Scalar O, GCCvsICC
  "25_ilp",
  ## Vectorized, GCCvsICC
  "51_simd_improved_rot"
)

flags = c(

  # changing main flags (with 51) -> Ofast/O3 wins
  ## "gcc",
  ## "gcc_O1",
  ## "gcc_O2",
  ## "gcc_O3",
  ## "gcc_Ofast"

  # icc vs gcc (with both)
  "icc_O3",
  "gcc_O3"

  # unroll (with 40/41) -> no effect
  ## "gcc_O3_unroll",
  ## "gcc_O3"
)

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
  for(flag in flags){

    ## read files
    files = c(
      paste("data/", flag, "-", exe, "/lion_64.csv", sep=""),
      paste("data/", flag, "-", exe, "/lion_128.csv", sep=""),
      paste("data/", flag, "-", exe, "/lion_256.csv", sep=""),
      paste("data/", flag, "-", exe, "/lion_512.csv", sep=""),
      paste("data/", flag, "-", exe, "/lion_1024.csv", sep=""),
      paste("data/", flag, "-", exe, "/lion_2048.csv", sep=""),
      paste("data/", flag, "-", exe, "/lion_4096.csv", sep="")
    )


    for(i in 1:7){
      if(file.exists(files[i])){
        csv = read.csv(files[i], sep=';')
        df_temp = data.frame(
          paste(flag, "-", exe),
          image_sizes[i],
          mean(csv$flops.cycle),
          mean(csv$cycles),
          mean(csv$flops)
        )
        df = rbind(df, df_temp)
      } else {
        print(paste("doesn't exist:", files[i]))
      }
    }
  }

}
colnames(df) = c("executable", "size", "performance", "cycles", "flops")


# baseplot

baseplot = ggplot(data=df,
                  aes(x=size,
                      group=executable,
                      color=executable)) +
    scale_x_continuous(
      minor_breaks=c(),
      breaks=c(64,512,1024,2048,4096),
      labels=c(64,512,1024,2048,4096),
      expand = c(0, 200)
    ) +
  theme_wsj() +
  ## scale_colour_wsj()+
    theme(
      plot.background = element_rect(fill="#eeeeee"),
      panel.background = element_rect(fill="#eeeeee"),
      ## panel.grid.major.x = element_line(colour = "black", linetype=3, size=rel(1)),
      axis.ticks = element_line(colour = "black", size = 1),
      axis.ticks.length=unit(.2, "cm"),
      text = element_text(size=22, family="mono"),
      plot.title = element_text(
        size=22, face="bold",
        # Scalar
        ## hjust=-0.26
        # SIMD
        ## hjust=-0.65
        # GCCvsICC
        hjust=-0.9
      ),
      plot.subtitle = element_text(
        size=22,
        # Vectorized
        ## hjust=-0.4,
        # SIMD
        ## hjust=-0.2
        # GCCvsICC
        hjust=-0.4
      ),
      ## axis.title.y = element_text(angle=0),
      axis.title.x = element_blank(),
      axis.title.y = element_blank(),
      axis.text.x = element_text(size=20, face="plain", family="mono"),
      axis.text.y = element_text(size=20, face="plain"),
      legend.position = "none",
      legend.title = element_blank())

scaleFUN <- function(x) sprintf("%.1f", x)

## Performance plot
performance_plot = baseplot +
  geom_line(aes(y=performance), lwd=2) +
  geom_point(aes(y=performance), lwd=4) +
  ggtitle(
    label = "GCC vs ICC for Scalar and SIMD Implementation",
    subtitle="Performance [flops/cycle] vs. input size") +
  scale_y_continuous(
    labels=scaleFUN,
    ## name=("[flops/cycle]"),
    # Scalar
    ## limits=c(0,2)
    # Vectorizied, GCCvsICC
    limits=c(0,6)
  ) +

  ## Scalar GCC Opts
  ## annotate("text", label = TeX("\\textbf{-march=native -O3 and -march=native -O2}"), parse=TRUE,
  ##          x=550, y=1.9, color="#00B0F6", size=7, hjust=0, family="mono") +
  ## annotate("text", label = TeX("\\textbf{-march=native -Ofast}"), parse = TRUE,
  ##          x=550, y=1.57, color="#E76BF3", size=7, hjust=0, family="mono") +
  ## annotate("text", label = TeX("\\textbf{-march=native -O1}"), parse=TRUE,
  ##          x=550, y=1.2, color="#A3A500", size=7, hjust=0, family="mono") +
  ## annotate("text", label = TeX("\\textbf{-march=native}"), parse=TRUE,
  ##          x=550, y=0.4, color="#F7756C", size=7, hjust=0, family="mono")

  ## Vectorized GCC Opts
  ## annotate("text", label = TeX("\\textbf{-march=native -Ofast}"), parse = TRUE,
  ##          x=1550, y=3.2, color="#E76BF3", size=7, hjust=0, family="mono") +
  ## annotate("text", label = TeX("\\textbf{-march=native -O3}"), parse=TRUE,
  ##          x=550, y=5.0, color="#00B0F6", size=7, hjust=0, family="mono") +
  ## annotate("text", label = TeX("\\textbf{-march=native -O2}"), parse = TRUE,
  ##          x=350, y=3.8, color="#00BF7D", size=7, hjust=0, family="mono") +
  ## annotate("text", label = TeX("\\textbf{-march=native -O1}"), parse=TRUE,
  ##          x=1750, y=4.2, color="#A3A500", size=7, hjust=0, family="mono") +
  ## annotate("text", label = TeX("\\textbf{-march=native}"), parse=TRUE,
  ##          x=550, y=1.0, color="#F7756C", size=7, hjust=0, family="mono")

  ## GCC vs ICC
  annotate("text", label = TeX("\\textbf{SIMD: GCC -march=native -O3}"), parse = TRUE,
           x=150, y=5.4, color="#C77CFF", size=7, hjust=0, family="mono") +
  annotate("text", label = TeX("\\textbf{SIMD: ICC -march=native -O3}"), parse=TRUE,
           x=150, y=3.3, color="#00BFC4", size=7, hjust=0, family="mono") +
  annotate("text", label = TeX("\\textbf{Scalar: GCC -march=native -O3}"), parse=TRUE,
           x=150, y=2.2, color="#7CAE00", size=7, hjust=0, family="mono") +
  annotate("text", label = TeX("\\textbf{Scalar: ICC -march=native -O3}"), parse=TRUE,
           x=150, y=0.7, color="#F7756C", size=7, hjust=0, family="mono")

    ##
    ## geom_label(aes(x=1800, y=2.4, label = "-O3"),
    ##            label.size = NA,
    ##            color = "#00B937", hjust=-1/32, vjust=1, lwd=8) +
    ## geom_label(aes(x=1600, y=3.2, label = "-Ofast -funsafe-math -ffast-math"),
    ##            label.size = NA,
    ##            color = "#609BFE", hjust=-1/32, vjust=1, lwd=8) +
    ## geom_label(aes(x=1500, y=0.8, label = "No flags"),
    ##            label.size = NA,
    ##            color = "#F7756C", hjust=-1/32, vjust=1, lwd=8)


    ## geom_label(aes(x=1500, y=1.7, label = "Scalar"),
    ##            label.size = NA,
    ##            color = "#00B937", hjust=-1/32, vjust=1, lwd=8) +
    ## geom_label(aes(x=1500, y=2.4, label = "SIMD"),
    ##            label.size = NA,
    ##            color = "#609BFE", hjust=-1/32, vjust=1, lwd=8) +
    ## geom_label(aes(x=1500, y=0.8, label = "Baseline"),
    ##            label.size = NA,
    ##            color = "#F7756C", hjust=-1/32, vjust=1, lwd=8)



    ## geom_label(aes(x=1300, y=1.95, label = "Copy Memory"),
    ##            label.size = NA,
    ##            color = "#7BAD00", hjust=-1/32, vjust=1, lwd=8) +
    ## geom_label(aes(x=900, y=3.7, label = "SIMD (without 90/270)"),
    ##            label.size = NA,
    ##            color = "#C67BFE", hjust=-1/32, vjust=1, lwd=8) +
    ## geom_label(aes(x=900, y=1.3, label = "Scalar (without 90/270)"),
    ##            label.size = NA,
    ##            color = "#00BEC3", hjust=-1/32, vjust=1, lwd=8) +
    ## geom_label(aes(x=1000, y=0.45, label = "Baseline"),
    ##            label.size = NA,
    ##            color = "#F7756C", hjust=-1/32, vjust=1, lwd=8)
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
       width=10,
       height=6)

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
