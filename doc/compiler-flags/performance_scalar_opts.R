library(ggplot2)
library(ggthemes)
library(svglite)
library(latex2exp)

executables = c(
  "25_ilp"
)

flags = c(
  "gcc",
  "gcc_O1",
  "gcc_O2",
  "gcc_O3",
  "gcc_Ofast"
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
        hjust=-0.26
      ),
      plot.subtitle = element_text(
        size=22,
        hjust=-0.4,
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
    label = "GCC Flags for Scalar Implementation",
    subtitle="Performance [flops/cycle] vs. input size") +
  scale_y_continuous(
    labels=scaleFUN,
    ## name=("[flops/cycle]"),
    # Scalar
    ## limits=c(0,2)
    # Vectorizied, GCCvsICC
    limits=c(0,2)
  ) +

  ## Scalar GCC Opts
  annotate("text", label = TeX("\\textbf{-march=native -O3 and -march=native -O2}"), parse=TRUE,
           x=550, y=1.9, color="#00B0F6", size=7, hjust=0, family="mono") +
  annotate("text", label = TeX("\\textbf{-march=native -Ofast}"), parse = TRUE,
           x=550, y=1.57, color="#E76BF3", size=7, hjust=0, family="mono") +
  annotate("text", label = TeX("\\textbf{-march=native -O1}"), parse=TRUE,
           x=550, y=1.2, color="#A3A500", size=7, hjust=0, family="mono") +
  annotate("text", label = TeX("\\textbf{-march=native}"), parse=TRUE,
           x=550, y=0.4, color="#F7756C", size=7, hjust=0, family="mono")

## Save plots
ggsave(file="performance_scalar_opts.pdf",
       dpi = 300,
       plot=performance_plot,
       width=10,
       height=6)
