## needed?
## library(tidyverse)
library(ggplot2)
library(ggthemes)
library(latex2exp)


pi_scalar <- 4
pi_simd <- 16
beta <- 9

log2 <- function(x) log(x)/log(2)

roofline <- function(x) ifelse(x*beta<pi_scalar, log2(x*beta), log2(pi_scalar))

executables = c(
  "0_baseline",
  "25_ilp",
  ## "35_simd",
  ## "40_ilp_norot_90_270",
  ## "41_simd_norot_90_270",
  "51_simd_improved_rot"
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

df = data.frame()
for (exe in executables) {

  ## read files
  files = c(
    paste("../performance-plot/data/", exe, "/lion_64.csv", sep=""),
    paste("../performance-plot/data/", exe, "/lion_128.csv", sep=""),
    paste("../performance-plot/data/", exe, "/lion_256.csv", sep=""),
    paste("../performance-plot/data/", exe, "/lion_512.csv", sep=""),
    paste("../performance-plot/data/", exe, "/lion_1024.csv", sep=""),
    paste("../performance-plot/data/", exe, "/lion_2048.csv", sep=""),
    paste("../performance-plot/data/", exe, "/lion_4096.csv", sep="")
  )

  for(i in 1:7){
    if(file.exists(files[i])){
      csv = read.csv(files[i], sep=';')
      cycles = mean(csv$cycles)
      flops = mean(csv$flops)
      bytes = 8*image_sizes[i]^2
      df_temp = data.frame(
        exe,
        image_sizes[i],
        cycles,
        flops,
        bytes,
        flops/cycles,
        flops/bytes
      )
      df = rbind(df, df_temp)
    }
  }
}

colnames(df) = c("opt", "n", "cycles", "flops", "bytes", "perf", "op_int")

roofline <- ggplot(mapping = aes(x = op_int, y = perf, color = opt), data = df) +

  theme_wsj() +
  theme(
    plot.background = element_rect(fill="#eeeeee"),
    panel.background = element_rect(fill="#eeeeee"),
    panel.grid.major.y = element_line(colour = "darkgray", linetype=3, size=rel(1)),
    axis.ticks = element_line(colour = "black", size = 1),
    axis.ticks.length=unit(.2, "cm"),
    text = element_text(size=22, family="mono"),
    plot.title = element_text(
      size=22, face="bold",
      hjust=-.07
    ),
    plot.subtitle = element_text(
      size=22,
      hjust=-0.12,
      ),
    ## axis.title.y = element_text(angle=0),
    axis.title.x = element_text(size=20, face="plain", family="mono",
                                margin = margin(t = 15, r = 0, b = 0, l = 0)
                                ),
    axis.title.y = element_blank(),
    axis.text.x = element_text(size=20, face="plain", family="mono"),
    axis.text.y = element_text(size=20, face="plain"),
    legend.position = "none",
    legend.title = element_blank()
  ) +

  # start/end size

  # roofline with dashed helper
  geom_vline(xintercept = pi_scalar / beta, linetype='dashed', color = 'darkgrey') +
  geom_vline(xintercept = pi_simd / beta, linetype='dashed', color = 'darkgrey') +
  geom_hline(yintercept = pi_scalar, color = 'purple', lwd=1) +
  geom_hline(yintercept = pi_simd, color = 'purple', lwd=1) +
  geom_abline(slope = 1, intercept = 3.2, color = 'purple', lwd=1) +

  # uncomment to see exact roofline function
  ## stat_function(fun = roofline) +

  scale_y_continuous(trans='log2', limits = c(NA, 20), breaks=c(1,4,16)) +
  scale_x_continuous(trans='log2', limits = c(0.08, NA)) +

  ## geom_text(data = subset(df, opt == '0_baseline' & n == 2048),
  ##           aes(x=1400, y=pi_scalar,
  ##               label = "Peak scal. (4 flops/cycle)"), colour = 'darkgrey', vjust=-1,
  ##           family="mono",
  ##           fontface="bold",
  ##           size=7
  ##           ) +

  annotate("text", label = expression(bold("Peak Scalar")), parse=TRUE,
           x=1.2, y=4.6, color="black", size=7, hjust=0, family="mono") +
  annotate("text", label = expression(bold("Peak SIMD")), parse=TRUE,
           x=3, y=18.6, color="black", size=7, hjust=0, family="mono") +
  annotate("text", label = expression(bold("Bandwidth")), parse=TRUE,
           angle=61.5,
           x=0.08, y=1, color="black", size=7, hjust=0, family="mono") +
  annotate("text", label = expression(bold("9 bytes/cycle")), parse=TRUE,
           angle=61.5,
           x=0.08, y=.5, color="black", size=7, hjust=0, family="mono") +

  ## geom_text(data = subset(df, opt == '0_baseline' & n == 2048), aes(x=1300, y=pi_simd, label = "Peak vec. (16 flops/cycle)"), colour = 'darkgrey', vjust=-1) +
  ## geom_text(data = subset(df, opt == '0_baseline' & n == 2048), aes(x=0.17, y=1.4, angle = 64, label = "Bandwidth (9 bytes/cycle)"), colour = 'darkgrey', vjust=-1) +

  geom_label(
    data = subset(df, n == 64 & opt != '0_baseline'),
    aes(label=n),
    vjust = 1,
    hjust = 1.1,
    label.size = NA,
    alpha = 0.0,
    family="mono",
    fontface="bold",
    size=7
  ) +
  geom_label(
    data = subset(df, n == 4096 & opt != '0_baseline'),
    aes(label=n),
    vjust = 1,
    hjust = -0.1,
    label.size = NA,
    alpha = 0.0,
    family="mono",
    fontface="bold",
    size=7
  ) +
  geom_label(
    data = subset(df, n == 64 & opt == '0_baseline'),
    aes(label=n),
    vjust = 0.3,
    hjust = 1.1,
    label.size = NA,
    alpha = 0.0,
    family="mono",
    fontface="bold",
    size=7
  ) +
  geom_label(
    data = subset(df, n == 2048 & opt == '0_baseline'),
    aes(label=n),
    vjust = -0.6,
    hjust = 0.6,
    label.size = NA,
    alpha = 0.0,
    family="mono",
    fontface="bold",
    size=7
  ) +


  geom_point(lwd=4) +
  geom_line(lwd=2) +


  annotate("text", label = TeX("\\textbf{SIMD}"), parse = TRUE,
           x=430, y=6.4, color="#619CFF", size=7, hjust=0, family="mono") +
  annotate("text", label = TeX("\\textbf{Scalar}"), parse=TRUE,
           x=430, y=2.3, color="#00BA38", size=7, hjust=0, family="mono") +
  annotate("text", label = TeX("\\textbf{Baseline}"), parse=TRUE,
           x=430, y=.7, color="#F8766D", size=7, hjust=0, family="mono") +

  ggtitle(
    label = "Roofline",
    subtitle="Performance [flops/cycle]") +

  labs(
    ## subtitle = 'Performance [flops/cycle]',
    x = 'Operational Intensity [flops/byte]'
  )

## plot(roofline)

ggsave(file="roofline.pdf",
       dpi = 300,
       plot=roofline,
       width=10,
       height=6)
