# author: jcuamatzi
# date:   2023-09-31

# set as working directory the path where this script is saved
setwd(dirname(rstudioapi::getActiveDocumentContext()$path))


rm(list = ls()) # Clean env


# Load necessary libraries


libraries <- c("ggplot2", "data.table", "dplyr", "ggh4x", "reshape2", "ggtext", "readxl")

for (lib in libraries) {
  if (!requireNamespace(lib, quietly = TRUE)) {
    suppressPackageStartupMessages(install.packages(lib, dependencies = TRUE))
  }
  suppressPackageStartupMessages(library(lib, character.only = TRUE))
}


rm(lib, libraries)


# 1) Read Data
df <- read_xlsx(path = "Data_GeneExpression_Pleos-DyeP.xlsx", sheet = "GeneExpressionData_2heatmap")

# Set levels for each gene in the heat map
df$Gene <- factor(df$Gene, levels = c("Pleos-Dyp4", "Pleos-Dyp2", "Pleos-Dyp1"))

# Fix typo in "glicerol"
df$Condition <- gsub("Glicerol", "Glycerol", df$Condition)

# Set levels for each fermentation
df$Condition <- factor(df$Condition, levels = c("Glucose", "Glycerol", "Glycerol+AYG"))

# Change to factor the column position
df$Position <- as.factor(df$Position)

# Plot using geom tile from ggplot2
heatmap.dyp <- df %>% 
  ggplot() +
  geom_tile(aes(y = Gene,  x = Position, fill = log2))+
  geom_hline(yintercept =0.5, color = "white", linewidth = 5)+
  geom_hline(yintercept =1.5, color = "white", linewidth = 5)+
  geom_hline(yintercept =2.5, color = "white", linewidth = 5)+
  #scale_fill_gradient(low = "white", high = "red", limits=c(-14, 14), breaks = seq(-14, 14, 4))+
  scale_fill_gradientn(colours = c("blue", "gray95", "red" ), limits=c(-14, 14), breaks = seq(-14, 14, 4))+
  scale_y_discrete(labels = c("Pleos-Dyp1" = "<i>Pleos</i>-DyeP1",
                              "Pleos-Dyp2" = "<i>Pleos</i>-DyeP2",
                              "Pleos-Dyp4" = "<i>Pleos</i>-DyeP4"))+
  scale_x_discrete(labels = c("1" = "120",
                              "2" = "144",
                              "3" = "168",
                              "4" = "216",
                              "5" = "240",
                              "6" = "312",
                              "7" = "360",
                              "8" = "504"))+
  #scale_fill_gradient(low="red", high="blue") +
  facet_grid(~Condition)+
  labs(x = "<br>Time (h)", y = "Gene<br>", fill = "Fold Change<br>(log<sub>2</sub>)")+
  theme_classic()+
  theme(
    # axis aesthetics
    axis.title.y = element_markdown(face = "bold", size = 14, color = "black"),
    axis.title.x = element_markdown(face = "bold", size = 14, color = "black"),
    # axis aesthetics
    axis.text.y = element_markdown(size = 13, color = "black"),
    axis.text.x = element_markdown(size = 13, color = "black"),
    # strip aesthetics
    strip.text = element_text(size = 15, color = "black", face = "bold"),
    strip.background = element_blank(),
    #
    legend.title = element_markdown(size = 12, face = "bold", color = "black"),
    #
    axis.line.y = element_blank(),
    
    axis.ticks.length.x = unit(0.2, "cm"),
    axis.ticks.length.y = unit(0.2, "cm")
  )

heatmap.dyp

ggsave(filename = "Figure_6.png", plot = heatmap.dyp, width = 13, height = 4, units = "in", dpi = 300)


rm(list = ls())


