#!/usr/bin/env Rscript

library(ggplot2)
library(scatterplot3d)
library(plyr)

arg <- commandArgs(trailingOnly=TRUE)

data <- read.csv(arg[1], header=FALSE, sep="\t")
names(data) <- c("species", "position", "coverage")

pdf("myplot.pdf", width=20, height=4)

ggplot(data, aes(x=position, y=coverage, colour=species, group=species)) + geom_line() + ggtitle(arg[2])

data2 <- ddply(data, .(species), summarise, med = median(coverage))
ggplot(data, aes(x = species, y = coverage)) + geom_boxplot() + geom_text(data = data2, aes(x = species, y = med, label = med), size =3, vjust = -1.0)

#ggplot (data, aes(x=species, y=coverage)) + geom_boxplot(outlier.size=0.2) + ggtitle(arg[2])
ggplot(data, aes(x=position, y=coverage, colour=species, group=species)) + geom_point(size=0.5) + stat_smooth(method="auto", se=FALSE, size=3) + ggtitle(arg[2])

scatterplot3d(data$species, data$position, data$coverage, highlight.3d = TRUE)


dev.off()



#arg <- commandArgs(trailingOnly=TRUE)
#print(arg[1])
#print(arg[2])
#mydata <- read.table(arg[1], header=FALSE, sep="\t")
#fit <- 
#qplot(V2, V3, data = mydata, colour = V1, geom = "line", ylim=c(0,arg[2]))
#ggsave(file="Routfile.pdf")
#png(filename="test")
#plot(fit)
#dev.off()



