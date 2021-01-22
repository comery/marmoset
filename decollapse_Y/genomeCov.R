#++++     Credit: Giulio Formenti giulio.formenti@gmail.com     ++++

setwd(dirname(rstudioapi::getSourceEditorContext()$path))

library(ggplot2)
library(dplyr)
library(grid)
library(gridExtra)
library(readxl)
library(gtable)

read_excel_allsheets <- function(filename, tibble = FALSE) {
  sheets <- readxl::excel_sheets(filename)
  x <- lapply(sheets, function(X) readxl::read_excel(filename, sheet = X))
  if(!tibble) x <- lapply(x, as.data.frame)
  names(x) <- sheets
  x
}

plot<-function(data,xlim,yvalues,title){
  ggplot(data, aes(x = Coverage, y = get(yvalues), color=group)) + 
    geom_line(size=5) +
    geom_vline(aes(xintercept=32), color="black", linetype="dashed", size=5) +
    geom_vline(aes(xintercept=64), color="black", linetype="dashed", size=5) +
    theme(
      panel.background = element_blank(),
      axis.title = element_text(family = "Arial",
                                size = rel(10), colour = "black"),
      axis.title.y = element_text(margin = margin(t = 0, r = 50, b = 0, l = 50)),
      axis.title.x = element_text(margin = margin(t = 50, r = 0, b = 50, l = 0)),
      axis.text = element_text(family = "Arial",
                               size = rel(8), colour = "black"),
      axis.ticks.length = unit(2, "cm"),
      axis.ticks = element_line(colour = "black", size=(2)),
      plot.title = element_text(size = rel(15), face = "bold", hjust = 0.5,
                                margin = margin(t = 50, r = 0, b = 50, l = 0)),
      axis.line = element_line(colour = 'black', size = 2),
    )+ggtitle(title)+ylab(yvalues)+
      scale_x_continuous(expand = c(0, 0), limits = c(0,xlim)) + scale_y_continuous(expand = c(0, 0))
}

grid_arrange_shared_legend <- function(...) {
  plots <- list(...)
  g <- ggplotGrob(plots[[1]] + theme(
    legend.position="bottom", 
    legend.text = element_text(size = rel(8), face = "bold"),
    legend.key.size = unit.c(unit(10, "cm")),
    legend.title = element_blank(),
    legend.background = element_blank(),
    legend.key = element_rect(fill = NA, color = NA)))$grobs
  legend <- g[[which(sapply(g, function(x) x$name) == "guide-box")]]
  lheight <- sum(legend$height)
  grid.arrange(
    do.call(arrangeGrob, lapply(plots, function(x)
      x + theme(legend.position="none", plot.margin = unit(c(4,4,4,4), "cm")))),
    legend,
    ncol = 1,
    heights = unit.c(unit(1, "npc") - lheight, lheight))
}

data <- mysheets <- read_excel_allsheets("genomeCov.xlsx")

Y_unplaced<-bind_rows("Collapsed" = data$`original_Y+unplaced`, "Decollapsed (SDA)" = data$`SDA_Y+unplaced`, .id = "group")
Y_only<-bind_rows("original_Yonly" = data$original_Yonly, "SDA_Yonly" = data$SDA_Yonly, .id = "group")

png("chrY_SDA.png", width = 8000, height = 5000) 

gr1<-plot(Y_unplaced,1000,'Frequency','Y-linked and unplaced scaffolds')
gr2<-plot(Y_only,400,'Count', 'Y-linked only')

grid_arrange_shared_legend(gr1,gr2)

dev.off()

