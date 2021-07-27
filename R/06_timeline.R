## HEADER ####
## Florence Galliers
## BMSB Distribution Diagram
## Last Edited: 2020-05-11

# Install packages
library(openxlsx)
library(ggplot2)
library(hrbrthemes)
library(tidyverse)
library(directlabels)
library(viridis)


## Import Data
bmsb <- read.xlsx("data/06_timeline_distribution.xlsx")

head(bmsb)
unique(bmsb$continent)
bmsb$continent <- factor(bmsb$continent, levels = c("North America", "South America", "Africa", "Europe"))

bmsb_long <- pivot_longer(bmsb,
             cols = 3:6,
             names_to = "invasion",
             values_to = "date")

head(bmsb_long)

bmsb_long$invasion <- factor(bmsb_long$invasion, levels = c("end",
                                                            "widespread",
                                                            "established",
                                                            "first"))


ggplot(bmsb_long, aes(fill = invasion,
                      y = date,
                      x = region)) + 
  geom_bar(position = "stack", stat = "identity", alpha = 1) +
  ggtitle("Timeline of BMSB Invasion") +
  xlab("") +
  coord_flip(ylim = c(1995, 2021)) 
  





  ggplot(bmsb) +
    # first seen line
  geom_segment(aes(x = reorder(region, -first),
                   xend = region, 
                   y = first,
                   yend = end),
               color = "grey60", size = 0.7, linetype = 2) +
    # label lines
  geom_dl(aes(region, end, label = region), 
          method = list(dl.trans(x = x + 0.3), "last.points", cex = 0.65, rot = 0),
          size = 0.3) +
    # first seen points
 geom_point(aes(x = region, y = first), color = "grey60", size = 0.8 ) +
 # geom_point(aes(x = region, y = first), color = "white", size = 1.6 ) +
    # established populations line
  geom_segment(aes(x = region, xend = region, 
                   y = established, yend = end), color = "grey60", size = 1) +
    # established population points
  geom_point(aes(x = region, y = established), color = "grey60", size = 1.2) +
  #geom_point(aes(x = region, y = established), color = "#16A085", size = 2) +
    # widespread population line
  geom_segment(aes(x = region, xend = region,
                   y = widespread, yend = end), color = "black", size = 1) +
    # widespread population points
  geom_point(aes(x = region, y = widespread), color = "black", size = 1.2) +
  #geom_point(aes(x = region, y = widespread), color = "#2980B9", size = 2) +
    # facet by region
    facet_grid((continent) ~ ., 
               scales = "free", 
               space = "free",
               switch = "both") +
  coord_flip() +
  theme_void() +
  theme(
    axis.text.x = element_text(angle = 0),
    axis.title = element_text(),
    axis.line.x = element_line(),
    axis.ticks.x = element_line(),
    axis.ticks.length.x = unit(3, "pt"),
    plot.margin = margin(1, 1, 1, 1, "cm"),
    strip.text.y.left = element_text(angle = 0)
    ) +
    scale_y_continuous(limits=c(1995, 2023), breaks = seq(1995, 2021, 1)) +
  xlab("") +
  ylab("Year")

  ggsave("output/timeline-plot.png", last_plot(),
         height = 21, width = 35, units = "cm")

  

today <- as.character(Sys.Date())
bmsb$today <- today
  
bmsb$first <- as.Date(paste(bmsb$first, 1, 1, sep = "-"))

vistime(bmsb, 
        col.start = "first", 
        col.end = "today", 
        col.group = "continent", 
        col.event = "region", 
        col.color = "color", 
        col.fontcolor = "fontcolor", 
        optimize_y = TRUE, 
        linewidth = NULL, 
        title = "From Asia, to North America, to Europe: BMSB Invasion Timeline", 
        show_labels = TRUE, 
        background_lines = 10)
