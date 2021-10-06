library(readxl)
library(ggplot2)
library(dplyr)

ImageDataCombined <- read_excel("ImageAnalysis/ImageDataCombined.xlsx")

ImageDataCombined <- mutate(ImageDataCombined, Density = NumCells / Area)

mu <- group_by(ImageDataCombined, Type, Patient) %>%
  dplyr::summarise(medianNumCells = median(log10(NumCells)), 
                   medianDensity = median(Density), 
                   medianArea = median(log10(Area)))

quantiles <- 10^quantile(log10(ImageDataCombined$NumCells), 
         probs = c(0.25, 0.5, 0.75))

g1 <- ggplot(ImageDataCombined, aes(x = Type, y = log10(NumCells), 
                              color = factor(Patient), shape = factor(Patient))) + 
  geom_jitter(width = 0.2) + 
  geom_hline(aes(yintercept = log10(quantiles[1])),
             linetype = "dashed", col = "red") + 
  geom_hline(aes(yintercept = log10(quantiles[2])),
             linetype = "dashed", col = "black") + 
  geom_hline(aes(yintercept = log10(quantiles[3])),
             linetype = "dashed", col = "red") + 
  labs(color = "patient", shape = "patient") + 
  scale_x_discrete(name = "cancer type") + 
  scale_y_continuous(breaks = 0:4, labels = 10^(0:4), 
                     name = "number of cells") + 
  theme_bw()

g1 <- arrangeGrob(g1, top = textGrob("a", x = unit(0, "npc"), 
                                     y = unit(0.5, "npc"), 
                                     just = c("left", "top"),
                                     gp = gpar(fontsize = 18)))

g2 <- ggplot(ImageDataCombined, aes(x = Type, y = Density, 
                              color = factor(Patient), shape = factor(Patient))) + 
  geom_jitter(width = 0.2) + 
  labs(color = "patient", shape = "patient") + 
  scale_x_discrete(name = "cancer type") + 
  scale_y_continuous(name = expression('cell density (per mm'^2*')'), 
                     limits = c(0, 0.015)) + 
  theme_bw()

g2 <- arrangeGrob(g2, top = textGrob("b", x = unit(0, "npc"), 
                                     y = unit(0.5, "npc"), 
                                     just = c("left", "top"),
                                     gp = gpar(fontsize = 18)))

grid.arrange(g1, g2, ncol = 2)
