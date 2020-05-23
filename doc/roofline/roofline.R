library(tidyverse)

pi_scalar <- 4
pi_simd <- 16
beta <- 9

log2 <- function(x) log(x)/log(2)

roofline <- function(x) ifelse(x*beta<pi_scalar, log2(x*beta), log2(pi_scalar))

executables = c(
  "0_baseline",
  "25_ilp",
  "31_simd_precomp_rotations_no_bac_simd",
  "35_simd",
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

df = data.frame()
for (exe in executables) {
  
  ## read files
  files = c(
    paste("../performance-plot/data/", exe, "/lena_64.csv", sep=""),
    paste("../performance-plot/data/", exe, "/monkey_128.csv", sep=""),
    paste("../performance-plot/data/", exe, "/lena_256.csv", sep=""),
    paste("../performance-plot/data/", exe, "/lena_512.csv", sep=""),
    paste("../performance-plot/data/", exe, "/grey-parrot_1024.csv", sep=""),
    paste("../performance-plot/data/", exe, "/matterhorn_2048.csv", sep=""),
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
  # uncomment to see exact roofline function
  #stat_function(fun = roofline) +
  
  geom_point() +
  geom_line() +
  
  # start/end size
  geom_label(data = subset(df, n == 64 | n == 4096 | (opt == '0_baseline' & n == 2048)),
                            aes(label=n), vjust = -0.4, label.size = NA) +
  # roofline with dashed helper
  geom_vline(xintercept = pi_scalar / beta, linetype='dashed', color = 'darkgrey') +
  geom_hline(yintercept = pi_scalar, color = 'darkgrey') +
  geom_hline(yintercept = pi_simd, color = 'darkgrey') +
  geom_vline(xintercept = pi_simd / beta, linetype='dashed', color = 'darkgrey') +
  geom_abline(slope = 1, intercept = 3.2, color = 'darkgrey') +
  
  scale_y_continuous(trans='log2', limits = c(NA, 20)) +
  scale_x_continuous(trans='log2', limits = c(0.08, NA)) +
  
  geom_text(data = subset(df, opt == '0_baseline' & n == 2048), aes(x=1400, y=pi_scalar, label = "Peak scal. (4 flops/cycle)"), colour = 'darkgrey', vjust=-1) +
  geom_text(data = subset(df, opt == '0_baseline' & n == 2048), aes(x=1300, y=pi_simd, label = "Peak vec. (16 flops/cycle)"), colour = 'darkgrey', vjust=-1) +
  geom_text(data = subset(df, opt == '0_baseline' & n == 2048), aes(x=0.17, y=1.4, angle = 64, label = "Bandwidth (9 bytes/cycle)"), colour = 'darkgrey', vjust=-1) +
  
  labs(
    y = 'Performance [flops/cycle]',
    x = 'Operational Intensity [flops/byte]'
  ) +
  
  ggtitle('i7-8650U @ 1.9 GHz') +
  
  theme_light() +
  theme(legend.position = "none"
)

plot(roofline)

ggsave(file="roofline.pdf",
       dpi = 300,
       plot=roofline,
       width=8,
       height=5)
