# author:           jcuamatzi
# date(created):    2023-09-31
# date(update):     2024-02-07
# set as working directory the path where this script is saved
setwd(dirname(rstudioapi::getActiveDocumentContext()$path))
getwd()

rm(list = ls()) # Clean env


# Load necessary libraries


libraries <- c("ggplot2", "data.table", "dplyr", "ggh4x", "reshape2", "ggtext", "readxl", "tidyverse", 
               "ComplexHeatmap", "dendextend")

for (lib in libraries) {
  if (!requireNamespace(lib, quietly = TRUE)) {
    suppressPackageStartupMessages(install.packages(lib, dependencies = TRUE))
  }
  suppressPackageStartupMessages(library(lib, character.only = TRUE))
}


rm(lib, libraries)

## Plot an exploratory heatmap
# 1) Read Data
df.raw <- read_xlsx(path = "Data_GeneExpression_Pleos-DyeP.xlsx", sheet = "GeneExpressionData_2Heatmap_Raw")

df.raw <- df.raw %>% 
  pivot_longer(cols = R1:R3,
               names_to = "Replicate",
               values_to = "DeltaCt")

df.raw$log2 <- log2(df.raw$DeltaCt)

df.mean <- df.raw %>% 
  group_by(Condition, Gene, Time) %>% 
  reframe(MeanLog2 = mean(log2),
          SDLog2 = sd(log2))

df.mean <- df.mean %>% 
  mutate(MeanLog2 = ifelse(MeanLog2 == -Inf, 0, MeanLog2),
         SDLog2 = ifelse(SDLog2 == "NaN", 0, SDLog2))

# Add italics to gene name using html style
df.mean$Gene <- gsub("Pleos-DyP", "<i>Pleos-dyp</i>", df.mean$Gene)

# define the order for each gene
#levels(factor(df.mean$Gene))

df.mean$Gene <- factor(df.mean$Gene, levels = c("<i>Pleos-dyp</i>4", "<i>Pleos-dyp</i>2", "<i>Pleos-dyp</i>1"))


# Plot using geom tile from ggplot2
plot.heatmap.a <- df.mean %>% 
  ggplot() +
  geom_tile(aes(y = Gene,  x = as.factor(Time), fill = MeanLog2))+
  geom_hline(yintercept =0.5, color = "white", linewidth = 5)+
  geom_hline(yintercept =1.5, color = "white", linewidth = 5)+
  geom_hline(yintercept =2.5, color = "white", linewidth = 5)+
   geom_text(aes(y = Gene, x = as.factor(Time), label = sprintf("%.2f", MeanLog2) ), size = 2.5)+
  scale_fill_gradientn(colours = c("blue", "gray95", "red" ), limits=c(-14, 14), breaks = seq(-14, 14, 4))+
  facet_grid(~Condition) +
  labs(x = "<br>Time (h)", y = "Gene<br> <br> <br>", fill = "Relative Expression <br> Levels (log<sub>2</sub>)")+
  theme_classic()+
  theme(
    # axis aesthetics
    axis.title.y = element_markdown(face = "bold", size = 14, color = "black"),
    axis.title.x = element_markdown(face = "bold", size = 14, color = "black"),
    # axis aesthetics
    axis.text.y = element_markdown(size = 11, color = "black"),
    axis.text.x = element_markdown(size = 11, color = "black"),
    # strip aesthetics
    strip.text = element_text(size = 14, color = "black", face = "bold"),
    strip.background = element_blank(),
    #
    legend.title = element_markdown(size = 12, face = "bold", color = "black"),
    #
    axis.line.y = element_blank(),
    
    axis.ticks.length.x = unit(0.2, "cm"),
    axis.ticks.length.y = unit(0.2, "cm")
  )

plot.heatmap.b <- df.mean %>% 
  ggplot() +
  geom_tile(aes(y = Gene,  x = as.factor(Time), fill = MeanLog2))+
  geom_hline(yintercept =0.5, color = "white", linewidth = 5)+
  geom_hline(yintercept =1.5, color = "white", linewidth = 5)+
  geom_hline(yintercept =2.5, color = "white", linewidth = 5)+
  #geom_text(aes(y = Gene, x = as.factor(Time), label = sprintf("%.2f", MeanLog2) ), size = 2.5)+
  scale_fill_gradientn(colours = c("blue", "gray95", "red" ), limits=c(-14, 14), breaks = seq(-14, 14, 4))+
  facet_grid(~Condition) +
  labs(x = "<br>Time (h)", y = "Gene<br> <br> <br>", fill = "Relative Expression <br> Levels (log<sub>2</sub>)")+
  theme_classic()+
  theme(
    # axis aesthetics
    axis.title.y = element_markdown(face = "bold", size = 14, color = "black"),
    axis.title.x = element_markdown(face = "bold", size = 14, color = "black"),
    # axis aesthetics
    axis.text.y = element_markdown(size = 11, color = "black"),
    axis.text.x = element_markdown(size = 11, color = "black"),
    # strip aesthetics
    strip.text = element_text(size = 14, color = "black", face = "bold"),
    strip.background = element_blank(),
    #
    legend.title = element_markdown(size = 12, face = "bold", color = "black"),
    #
    axis.line.y = element_blank(),
    
    axis.ticks.length.x = unit(0.2, "cm"),
    axis.ticks.length.y = unit(0.2, "cm")
  )

plot.heatmap.a
plot.heatmap.b

# ggsave(filename = "Figure_6a.png", plot = plot.heatmap.a,
#         width = 13, height = 4, units = "in", dpi = 300) # Numbers
# 
getwd()
ggsave(filename = "Figure_6.png", plot = plot.heatmap.b,
        width = 13, height = 4, units = "in", dpi = 300) # No Numbers

# ### NO LOG2
# df.raw
# 
# df.mean <- df.raw %>% 
#   group_by(Condition, Gene, Time) %>% 
#   reframe(Mean = mean(DeltaCt),
#           SD = sd(DeltaCt))
# 
# df.mean <- df.mean %>% 
#   mutate(Mean = ifelse(Mean == -Inf, 0, Mean),
#          SD = ifelse(SD == "NaN", 0, SD))
# 
# # Add italics to gene name using html style
# df.mean$Gene <- gsub("Pleos-DyP", "<i>Pleos-dyp</i>", df.mean$Gene)
# 
# # define the order for each gene
# #levels(factor(df.mean$Gene))
# 
# df.mean$Gene <- factor(df.mean$Gene, levels = c("<i>Pleos-dyp</i>4", "<i>Pleos-dyp</i>2", "<i>Pleos-dyp</i>1"))
# 
# 
# # Plot using geom tile from ggplot2
# plot.heatmap.c <- df.mean %>% 
#   ggplot() +
#   geom_tile(aes(y = Gene,  x = as.factor(Time), fill = Mean))+
#   geom_hline(yintercept =0.5, color = "white", linewidth = 5)+
#   geom_hline(yintercept =1.5, color = "white", linewidth = 5)+
#   geom_hline(yintercept =2.5, color = "white", linewidth = 5)+
#   geom_text(aes(y = Gene, x = as.factor(Time), label = sprintf("%.2f", Mean) ), size = 2.5)+
#   #scale_fill_gradientn(colours = c("blue", "gray95", "red" ), limits=c(-14, 14), breaks = seq(-14, 14, 4))+
#  scale_fill_gradientn(colours = c("gray90", "red" ))+
#   facet_grid(~Condition) +
#   labs(x = "<br>Time (h)", y = "Gene<br> <br> <br>", fill = "Relative Expression <br> Levels")+
#   theme_classic()+
#   theme(
#     # axis aesthetics
#     axis.title.y = element_markdown(face = "bold", size = 14, color = "black"),
#     axis.title.x = element_markdown(face = "bold", size = 14, color = "black"),
#     # axis aesthetics
#     axis.text.y = element_markdown(size = 11, color = "black"),
#     axis.text.x = element_markdown(size = 11, color = "black"),
#     # strip aesthetics
#     strip.text = element_text(size = 14, color = "black", face = "bold"),
#     strip.background = element_blank(),
#     #
#     legend.title = element_markdown(size = 12, face = "bold", color = "black"),
#     #
#     axis.line.y = element_blank(),
#     
#     axis.ticks.length.x = unit(0.2, "cm"),
#     axis.ticks.length.y = unit(0.2, "cm")
#   ); plot.heatmap.c
# 
# plot.heatmap.d <- df.mean %>% 
#   ggplot() +
#   geom_tile(aes(y = Gene,  x = as.factor(Time), fill = Mean))+
#   geom_hline(yintercept =0.5, color = "white", linewidth = 5)+
#   geom_hline(yintercept =1.5, color = "white", linewidth = 5)+
#   geom_hline(yintercept =2.5, color = "white", linewidth = 5)+
#   #geom_text(aes(y = Gene, x = as.factor(Time), label = sprintf("%.2f", MeanLog2) ), size = 2.5)+
#   #scale_fill_gradientn(colours = c("blue", "gray95", "red" ), limits=c(-14, 14), breaks = seq(-14, 14, 4))+
#   scale_fill_gradientn(colours = c("gray90", "red" ))+
#   facet_grid(~Condition) +
#   labs(x = "<br>Time (h)", y = "Gene<br> <br> <br>", fill = "Relative Expression <br> Levels (log<sub>2</sub>)")+
#   theme_classic()+
#   theme(
#     # axis aesthetics
#     axis.title.y = element_markdown(face = "bold", size = 14, color = "black"),
#     axis.title.x = element_markdown(face = "bold", size = 14, color = "black"),
#     # axis aesthetics
#     axis.text.y = element_markdown(size = 11, color = "black"),
#     axis.text.x = element_markdown(size = 11, color = "black"),
#     # strip aesthetics
#     strip.text = element_text(size = 14, color = "black", face = "bold"),
#     strip.background = element_blank(),
#     #
#     legend.title = element_markdown(size = 12, face = "bold", color = "black"),
#     #
#     axis.line.y = element_blank(),
#     
#     axis.ticks.length.x = unit(0.2, "cm"),
#     axis.ticks.length.y = unit(0.2, "cm")
#   );plot.heatmap.d
# 
# plot.heatmap.c
# 
# 
# ggsave(filename = "Figure_6c.png", plot = plot.heatmap.c,
#         width = 13, height = 4, units = "in", dpi = 300) # Numbers
# # 
# ggsave(filename = "Figure_6d.png", plot = plot.heatmap.d,
#         width = 13, height = 4, units = "in", dpi = 300) # No Numbers




# 
rm(list = ls())



# Bar plot of delta delta ct
# Reading the data without log2 transformation
df.ddCt <- read_xlsx(path = "Data_GeneExpression_Pleos-DyeP.xlsx", sheet = "Data_DeltaCt")

df.ddCt <- df.ddCt %>% 
  pivot_longer(cols = R1:R3,
               names_to = "Replicates",
               values_to = "ddCt")
##

df.ddCt <- df.ddCt %>% group_by(Gene, Time, Condition) %>% 
  mutate(Log2ddCt = log2(ddCt)) %>% 
  reframe(MeanLog2ddCt = mean(Log2ddCt),
          SDLog2ddCt = sd(Log2ddCt))

# plot.bar.geneExpression <- df.mean %>% 

df.ddCt$Gene <- gsub("Pleos-dyp", "<i>Pleos-dyp</i>", df.ddCt$Gene)

df.ddCt <- df.ddCt %>% 
  mutate(MeanLog2ddCt = ifelse(Gene == "<i>Pleos-dyp</i>4" & 
                                 (Condition == "Glycerol" | Condition == "Glycerol + AYG"), 0, MeanLog2ddCt ),
         SDLog2ddCt = ifelse(Gene == "<i>Pleos-dyp</i>4" & 
                               (Condition == "Glycerol" | Condition == "Glycerol + AYG"), 0, SDLog2ddCt ))

plot.bar.geneExpression.c <- df.ddCt %>% 
  #filter(Condition == "Glucose") %>% 
  ggplot(aes(x = as.factor(Time), y = MeanLog2ddCt, fill = Gene)) +
  geom_col(color = "black",
           position = position_dodge(), 
           alpha = 0.75) +
  geom_hline(yintercept = 0) +
  scale_y_continuous(limits = c(-15, 15), 
                     breaks = seq(-15, 15, 3),
                     minor_breaks = seq(-15, 15, by = 1),
                     expand = c(0,0),
                     guide = "axis_minor") +
  geom_errorbar(aes(ymin = MeanLog2ddCt - SDLog2ddCt,
                    ymax = MeanLog2ddCt + SDLog2ddCt),
                position = position_dodge(width = 0.8), # Adjust width as needed
                width = 0.25, # Adjust the width of the error bars
                color = "black") +
  facet_wrap( ~Condition, nrow = 3, scales = "free") +
  labs(x = "Time (h)", y = paste0("Relative Expression Levels (log<sub>2</sub>)") )+
  
  scale_fill_manual(values = c("darkblue", "darkorange", "darkgreen"))+
  theme_classic()+
  theme(legend.text = element_markdown(),
        legend.title = element_markdown(face = "bold"),
        axis.title.y = element_markdown(face = "bold", color = "black"),
        axis.title.x = element_markdown(face = "bold", color = "black"),
        axis.text = element_text(color = "black", size = 10),
        ggh4x.axis.ticks.length.minor = rel(0.75)); plot.bar.geneExpression.c

getwd()
ggsave(filename = "Figure_7.png", 
       plot = plot.bar.geneExpression.c, width = 8, height = 10, units = "in", dpi = 300)


# # with numbers
# plot.bar.geneExpression.d <- df.ddCt %>% 
#   #filter(Condition == "Glucose") %>% 
#   ggplot(aes(x = as.factor(Time), y = MeanLog2ddCt, fill = Gene)) +
#   geom_col(color = "black",
#            position = position_dodge(), 
#            alpha = 0.75) +
#   geom_hline(yintercept = 0) +
#   geom_errorbar(aes(ymin = MeanLog2ddCt - SDLog2ddCt,
#                     ymax = MeanLog2ddCt + SDLog2ddCt),
#                 position = position_dodge(width = 0.85), # Adjust width as needed
#                 width = 0.25, # Adjust the width of the error bars
#                 color = "black") +
#   scale_y_continuous(limits = c(-15, 15), 
#                      breaks = seq(-15, 15, 3),
#                      minor_breaks = seq(-15, 15, by = 1),
#                      expand = c(0,0),
#                      guide = "axis_minor") +
#   geom_text(aes(y = 12, label = sprintf("%.2f", MeanLog2ddCt)), 
#             position = position_dodge(width = 0.85),
#             size = 2) +
#   facet_wrap( ~Condition, nrow = 3, scales = "free") +
#   labs(x = "Time (h)", y = paste0("Relative Expression Levels (log<sub>2</sub>)") )+
#   
#   scale_fill_manual(values = c("darkblue", "darkorange", "darkgreen"))+
#   theme_classic()+
#   theme(legend.text = element_markdown(),
#         legend.title = element_markdown(face = "bold"),
#         axis.title.y = element_markdown(face = "bold", color = "black"),
#         axis.title.x = element_markdown(face = "bold", color = "black"),
#         axis.text = element_text(color = "black", size = 10),
#         strip.text = element_text(color = "black", face = "bold", size = 12),
#         ggh4x.axis.ticks.length.minor = rel(0.75)); plot.bar.geneExpression.d
# 
# ggsave(filename = "Figure_GeneExpression.DDCt.d.20240209.png",
#        plot = plot.bar.geneExpression.d, width = 8, height = 10, units = "in", dpi = 300)
# 
# 

# Bar plot of delta delta ct
# Reading the data without log2 transformation
# Analysis : glu + ayg using glu as control; gly + ayg using gly as control
# rm(list = ls())
# df.ddCt <- read_xlsx(path = "Data_GeneExpression_Pleos-DyeP.xlsx", sheet = "Data_DeltaCt_2")
# 
# df.ddCt <- df.ddCt %>% 
#   pivot_longer(cols = R1:R3,
#                names_to = "Replicates",
#                values_to = "ddCt")
# ##
# 
# df.ddCt <- df.ddCt %>% group_by(Gene, Time, Condition) %>% 
#   mutate(Log2ddCt = log2(ddCt)) %>% 
#   reframe(MeanLog2ddCt = mean(Log2ddCt),
#           SDLog2ddCt = sd(Log2ddCt))
# 
# # plot.bar.geneExpression <- df.mean %>% 
# 
# df.ddCt$Gene <- gsub("Pleos-dyp", "<i>Pleos-dyp</i>", df.ddCt$Gene)
# 
# df.ddCt <- df.ddCt %>% 
#   mutate(MeanLog2ddCt = ifelse(Gene == "<i>Pleos-dyp</i>4" & 
#                                  (Condition == "Glycerol + AYG"), 0, MeanLog2ddCt ),
#          SDLog2ddCt = ifelse(Gene == "<i>Pleos-dyp</i>4" & 
#                                (Condition == "Glycerol + AYG"), 0, SDLog2ddCt ))
# 
# plot.bar.geneExpression.e <- df.ddCt %>% 
#   #filter(Condition == "Glucose") %>% 
#   ggplot(aes(x = as.factor(Time), y = MeanLog2ddCt, fill = Gene)) +
#   geom_col(color = "black",
#            position = position_dodge(), 
#            alpha = 0.75) +
#   geom_hline(yintercept = 0) +
#   scale_y_continuous(limits = c(-15, 15), 
#                      breaks = seq(-15, 15, 3),
#                      minor_breaks = seq(-15, 15, by = 1),
#                      expand = c(0,0),
#                      guide = "axis_minor") +
#   geom_errorbar(aes(ymin = MeanLog2ddCt - SDLog2ddCt,
#                     ymax = MeanLog2ddCt + SDLog2ddCt),
#                 position = position_dodge(width = 0.8), # Adjust width as needed
#                 width = 0.25, # Adjust the width of the error bars
#                 color = "black") +
#   facet_wrap( ~Condition, nrow = 3, scales = "free") +
#   labs(x = "Time (h)", y = paste0("Relative Expression Levels (log<sub>2</sub>)") )+
#   
#   scale_fill_manual(values = c("darkblue", "darkorange", "darkgreen"))+
#   theme_classic()+
#   theme(legend.text = element_markdown(),
#         legend.title = element_markdown(face = "bold"),
#         axis.title.y = element_markdown(face = "bold", color = "black"),
#         axis.title.x = element_markdown(face = "bold", color = "black"),
#         axis.text = element_text(color = "black", size = 10),
#         ggh4x.axis.ticks.length.minor = rel(0.75)); plot.bar.geneExpression.e
# 
# 
# ggsave(filename = "Figure_GeneExpression.DDCt.e.20240209.png", 
#        plot = plot.bar.geneExpression.e, width = 8, height = 10, units = "in", dpi = 300)
# 
# 
# # with numbers
# plot.bar.geneExpression.f <- df.ddCt %>% 
#   #filter(Condition == "Glucose") %>% 
#   ggplot(aes(x = as.factor(Time), y = MeanLog2ddCt, fill = Gene)) +
#   geom_col(color = "black",
#            position = position_dodge(), 
#            alpha = 0.75) +
#   geom_hline(yintercept = 0) +
#   geom_errorbar(aes(ymin = MeanLog2ddCt - SDLog2ddCt,
#                     ymax = MeanLog2ddCt + SDLog2ddCt),
#                 position = position_dodge(width = 0.85), # Adjust width as needed
#                 width = 0.25, # Adjust the width of the error bars
#                 color = "black") +
#   scale_y_continuous(limits = c(-15, 15), 
#                      breaks = seq(-15, 15, 3),
#                      minor_breaks = seq(-15, 15, by = 1),
#                      expand = c(0,0),
#                      guide = "axis_minor") +
#   geom_text(aes(y = 12, label = sprintf("%.2f", MeanLog2ddCt)), 
#             position = position_dodge(width = 0.85),
#             size = 2) +
#   facet_wrap( ~Condition, nrow = 3, scales = "free") +
#   labs(x = "Time (h)", y = paste0("Relative Expression Levels (log<sub>2</sub>)") )+
#   
#   scale_fill_manual(values = c("darkblue", "darkorange", "darkgreen"))+
#   theme_classic()+
#   theme(legend.text = element_markdown(),
#         legend.title = element_markdown(face = "bold"),
#         axis.title.y = element_markdown(face = "bold", color = "black"),
#         axis.title.x = element_markdown(face = "bold", color = "black"),
#         axis.text = element_text(color = "black", size = 10),
#         strip.text = element_text(color = "black", face = "bold", size = 12),
#         ggh4x.axis.ticks.length.minor = rel(0.75)); plot.bar.geneExpression.f
# 
# ggsave(filename = "Figure_GeneExpression.DDCt.f.20240209.png",
#        plot = plot.bar.geneExpression.f, width = 8, height = 10, units = "in", dpi = 300)
# 
# 
# # without log2
# # Analysis : glu + ayg using glu as control; gly + ayg using gly as control
# rm(list = ls())
# df.ddCt <- read_xlsx(path = "Data_GeneExpression_Pleos-DyeP.xlsx", sheet = "Data_DeltaCt_2")
# 
# df.ddCt <- df.ddCt %>% 
#   pivot_longer(cols = R1:R3,
#                names_to = "Replicates",
#                values_to = "ddCt")
# ##
# 
# df.ddCt <- df.ddCt %>% group_by(Gene, Time, Condition) %>% 
#   #mutate(Log2ddCt = log2(ddCt)) %>% 
#   reframe(MeanddCt = mean(ddCt),
#           SDddCt = sd(ddCt))
# 
# # plot.bar.geneExpression <- df.mean %>% 
# 
# df.ddCt$Gene <- gsub("Pleos-dyp", "<i>Pleos-dyp</i>", df.ddCt$Gene)
# 
# df.ddCt <- df.ddCt %>% 
#   mutate(MeanLog2ddCt = ifelse(Gene == "<i>Pleos-dyp</i>4" & 
#                                  (Condition == "Glycerol + AYG"), 0, MeanddCt ),
#          SDLog2ddCt = ifelse(Gene == "<i>Pleos-dyp</i>4" & 
#                                (Condition == "Glycerol + AYG"), 0, SDddCt ))
# 
# plot.bar.geneExpression.g <- df.ddCt %>% 
#   #filter(Condition == "Glucose") %>% 
#   ggplot(aes(x = as.factor(Time), y = MeanddCt, fill = Gene)) +
#   geom_col(color = "black",
#            position = position_dodge(), 
#            alpha = 0.75) +
#   geom_hline(yintercept = 0) +
#   #scale_y_continuous(limits = c(-15, 15), breaks = seq(-15, 15, 3), minor_breaks = seq(-15, 15, by = 1), expand = c(0,0), guide = "axis_minor") +
#   geom_errorbar(aes(ymin = MeanddCt - SDddCt,
#                     ymax = MeanddCt + SDddCt),
#                 position = position_dodge(width = 0.8), # Adjust width as needed
#                 width = 0.25, # Adjust the width of the error bars
#                 color = "black") +
#   facet_wrap( ~Condition, nrow = 3, scales = "free") +
#   labs(x = "Time (h)", y = paste0("Relative Expression Levels") )+
#   
#   scale_fill_manual(values = c("darkblue", "darkorange", "darkgreen"))+
#   theme_classic()+
#   theme(legend.text = element_markdown(),
#         legend.title = element_markdown(face = "bold"),
#         axis.title.y = element_markdown(face = "bold", color = "black"),
#         axis.title.x = element_markdown(face = "bold", color = "black"),
#         axis.text = element_text(color = "black", size = 10),
#         ggh4x.axis.ticks.length.minor = rel(0.75)); plot.bar.geneExpression.g
# 
# 
# ggsave(filename = "Figure_GeneExpression.DDCt.g.NoLog2.20240209.png", 
#        plot = plot.bar.geneExpression.g, width = 8, height = 10, units = "in", dpi = 300)
# 
# 
# # with numbers
# plot.bar.geneExpression.h <- df.ddCt %>% 
#   #filter(Condition == "Glucose") %>% 
#   ggplot(aes(x = as.factor(Time), y = MeanddCt, fill = Gene)) +
#   geom_col(color = "black",
#            position = position_dodge(), 
#            alpha = 0.75) +
#   geom_hline(yintercept = 0) +
#   geom_errorbar(aes(ymin = MeanddCt - SDddCt,
#                     ymax = MeanddCt + SDddCt),
#                 position = position_dodge(width = 0.85), # Adjust width as needed
#                 width = 0.25, # Adjust the width of the error bars
#                 color = "black") +
#   # scale_y_continuous(limits = c(-15, 15), 
#   #                    breaks = seq(-15, 15, 3),
#   #                    minor_breaks = seq(-15, 15, by = 1),
#   #                    expand = c(0,0),
#   #                    guide = "axis_minor") +
#   geom_text(aes(y =  MeanddCt + (MeanddCt*0.20), label = sprintf("%.2f", MeanddCt)), 
#             position = position_dodge(width = 0.85),
#             size = 2) +
#   facet_wrap( ~Condition, nrow = 3, scales = "free") +
#   labs(x = "Time (h)", y = paste0("Relative Expression Levels") )+
#   
#   scale_fill_manual(values = c("darkblue", "darkorange", "darkgreen"))+
#   theme_classic()+
#   theme(legend.text = element_markdown(),
#         legend.title = element_markdown(face = "bold"),
#         axis.title.y = element_markdown(face = "bold", color = "black"),
#         axis.title.x = element_markdown(face = "bold", color = "black"),
#         axis.text = element_text(color = "black", size = 10),
#         strip.text = element_text(color = "black", face = "bold", size = 12),
#         ggh4x.axis.ticks.length.minor = rel(0.75)); plot.bar.geneExpression.h
# 
# ggsave(filename = "Figure_GeneExpression.DDCt.h.NoLog2.20240209.png",
#        plot = plot.bar.geneExpression.h, width = 8, height = 10, units = "in", dpi = 300)

rm(list = ls())
  
