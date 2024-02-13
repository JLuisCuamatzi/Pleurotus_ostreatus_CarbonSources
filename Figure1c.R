# author: jcuamatzi
# date:   2023-09-31

# set as working directory the path where this script is saved
setwd(dirname(rstudioapi::getActiveDocumentContext()$path))

getwd()
rm(list = ls()) # Clean env


libraries <- c("ggplot2", "data.table", "dplyr", "ggh4x", "reshape2", "ggtext", "readxl", "tidyverse", 
               "ComplexHeatmap", "dendextend")

for (lib in libraries) {
  if (!requireNamespace(lib, quietly = TRUE)) {
    suppressPackageStartupMessages(install.packages(lib, dependencies = TRUE))
  }
  suppressPackageStartupMessages(library(lib, character.only = TRUE))
}


rm(lib, libraries)


df.1 <- read_xlsx(path = "Data_MycelialGrowthRate.xlsx", sheet = "Data2Plot")

df.1$Condition <- factor(df.1$Condition, levels = c("Control", "AYG", "RBBR", "AB129"))

plot.growthRate <- df.1 %>% 
  ggplot(aes(x = Condition, y = MycelialGrowthRate, fill = Condition)) +
  geom_col(aes(x = Condition, y = MycelialGrowthRate)) +
  facet_wrap(~CarbonSource) +
  scale_fill_manual(values = c("gray", "yellow3", "royalblue4", "deepskyblue2" )) +
  scale_y_continuous(limits = c(0, 6), breaks = seq(0, 6, 1),
                     minor_breaks = seq(0, 6, by = 0.5),
                     expand = c(0,0),
                     guide = "axis_minor")+
  labs(x = "\nDye", y = "Mycelial Growth Rate (mm/day)\n") +
  theme_classic() +
  theme(
    legend.title = element_text(face = "bold", color = "black", size = 12),
    axis.title = element_text(face = "bold", size = 12, color = "black"),
    ggh4x.axis.ticks.length.minor = rel(0.75)
    
  )
  
plot.growthRate

ggsave(filename = "Figure_1c.png", plot = plot.growthRate, dpi = 300, width = 8, height = 5, units = "in")  


# with numbers
plot.growthRate.b <- df.1 %>% 
  ggplot(aes(x = Condition, y = MycelialGrowthRate, fill = Condition)) +
  geom_col(aes(x = Condition, y = MycelialGrowthRate)) +
  geom_text(aes(y = MycelialGrowthRate + 0.5, label = sprintf("%.4f", MycelialGrowthRate))) +
  facet_wrap(~CarbonSource) +
  scale_fill_manual(values = c("gray", "yellow3", "royalblue4", "deepskyblue2" )) +
  scale_y_continuous(limits = c(0, 6), breaks = seq(0, 6, 1),
                     minor_breaks = seq(0, 6, by = 0.5),
                     expand = c(0,0),
                     guide = "axis_minor")+
  labs(x = "\nDye", y = "Mycelial Growth Rate (mm/day)\n") +
  theme_classic() +
  theme(
    legend.title = element_text(face = "bold", color = "black", size = 12),
    axis.title = element_text(face = "bold", size = 12, color = "black"),
    ggh4x.axis.ticks.length.minor = rel(0.75)
    
  )


# ggsave(filename = "Plot.RadialGrowthRate.b.png", plot = plot.growthRate.b, dpi = 300, width = 8, height = 5, units = "in")  

## Stat.test

kruskal.test(MycelialGrowthRate ~ CarbonSource, data = df.1)

rm(list = ls())



