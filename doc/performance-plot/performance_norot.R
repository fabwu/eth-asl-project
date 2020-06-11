library(ggplot2)
library(ggthemes)
library(svglite)
library(latex2exp)

executables = c(
    "40_ilp_norot_90_270",
    "41_simd_norot_90_270"
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

    ## read files
    files = c(
        ## new images
        paste("data/", exe, "/lion_64.csv", sep=""),
        paste("data/", exe, "/lion_128.csv", sep=""),
        paste("data/", exe, "/lion_256.csv", sep=""),
        paste("data/", exe, "/lion_512.csv", sep=""),
        paste("data/", exe, "/lion_1024.csv", sep=""),
        paste("data/", exe, "/lion_2048.csv", sep=""),
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
}
colnames(df) = c("executable", "size", "performance", "cycles", "flops")

scaleFUN <- function(x) sprintf("%.1f", x)

## baseplot
baseplot = ggplot(data=df,
                  aes(x=size,
                      group=executable,
                      color=executable)) +
  scale_x_continuous(
    minor_breaks=c(),
    breaks=c(64,512,1024,2048,4096),
    labels=c(64,512,1024,2048,4096),
    expand = c(0, 200)
  )

## Performance plot
performance_plot = baseplot +
  theme_wsj() +
  theme(
    plot.background = element_rect(fill="#eeeeee"),
    panel.background = element_rect(fill="#eeeeee"),
    ## panel.grid.major.x = element_line(colour = "black", linetype=3, size=rel(1)),
    axis.ticks = element_line(colour = "black", size = 1),
    axis.ticks.length=unit(.2, "cm"),
    text = element_text(size=22, family="mono"),
    plot.title = element_text(
      size=22, face="bold",
      hjust=-0.13
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
    legend.title = element_blank()
  ) +
  geom_line(aes(y=performance), lwd=2) +
  geom_point(aes(y=performance), lwd=4) +
  ggtitle(
    label = "Test Implementations",
    subtitle="Performance [flops/cycle] vs. input size") +
  scale_y_continuous(
    labels=scaleFUN,
    limits=c(-.5,6)
  ) +
  annotate("text", label = TeX("\\textbf{SIMD: without 90/270 deg. rotation}"), parse=TRUE,
           x=150, y=4.2, color="#00BFC4", size=7, hjust=0, family="mono") +
  annotate("text", label = TeX("\\textbf{Scalar: without 90/270 deg. rotation}"), parse=TRUE,
           x=150, y=1, color="#F7756C", size=7, hjust=0, family="mono")

## Save plots
ggsave(file="performance_norot.pdf",
       dpi = 300,
       plot=performance_plot,
       width=10,
       height=6)
