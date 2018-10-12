require(ggplot2)
require(reshape2)

setwd("~/Documents/Postdoc/assigning/revisions/")

# plot should be 2x1 panel, memory (Gb) and CPU (minutes) vs sketch size
# for each species, own coloured line (colourblind pallette)
# in key, give N and L for the species
# note that run on 8 cores (increasing memory)

sketch_times <- read.delim("~/Documents/Postdoc/assigning/revisions/sketch_size/sketch_times.txt", 
                           header=FALSE, stringsAsFactors=FALSE)
colnames(sketch_times) = c("Species", "Sketch size", "CPU", "Memory")
sketch_times$CPU = sketch_times$CPU/3600
sketch_times$Memory = sketch_times$Memory/(1024^2)
sketch_times = melt(sketch_times, id.vars = c("Species", "Sketch size"))

sketch_cpu = sketch_times[sketch_times$variable=="CPU",]
colnames(sketch_cpu) = c("Species", "sketch_size", "Resource", "time")
sketch_memory = sketch_times[sketch_times$variable=="Memory",]
colnames(sketch_memory) = c("Species", "sketch_size", "Resource", "mem_use")

ggplot(data = sketch_cpu, mapping=aes(x=sketch_size, y=time, colour=factor(Species))) + 
  geom_segment(aes(x=1E4, xend=1E4, y=0, yend=8), colour="red", linetype = "longdash") + 
  geom_point(size = 3) +
  geom_path(size = 1.5) +
  #scale_x_continuous() +
  #scale_x_log10(breaks=c(3000,10000,30000,100000,300000)) +
  #scale_y_log10(breaks=c(0.1, 0.5, 1, 1.5, 10)) +
  xlab("Sketch size") +
  ylab("CPU time (hours)") + 
  scale_colour_manual(values=c("#999999", "#E69F00", "#56B4E9"), 
                      name="Species",
                      breaks=c("pneumo", "staph", "listeria"),
                      labels=c("S. pneumoniae\n(N = 616; L = 2Mb)", 
                               "S. aureus\n(N = 284; L = 3Mb)",
                               "L. moncytogenes\n(N = 128; L = 3Mb)")) +
  theme_bw(base_size = 18) + 
  theme(legend.key.size = unit(3, 'lines'))

ggplot(data = sketch_memory, mapping=aes(x=sketch_size, y=mem_use, colour=factor(Species))) + 
  geom_segment(aes(x=1E4, xend=1E4, y=0, yend=18), colour="red", linetype = "longdash") + 
  geom_point(size = 3) +
  geom_path(size = 1.5) +
  #scale_x_continuous() +
  #scale_x_log10(breaks=c(3000,10000,30000,100000,300000)) +
  #scale_y_log10(breaks=c(0.1, 0.5, 1, 1.5, 10)) +
  xlab("Sketch size") +
  ylab("Max memory use (Gb)") + 
  scale_colour_manual(values=c("#999999", "#E69F00", "#56B4E9"), 
                      name="Species",
                      breaks=c("pneumo", "staph", "listeria"),
                      labels=c("S. pneumoniae\n(N = 616; L = 2Mb)", 
                               "S. aureus\n(N = 284; L = 3Mb)",
                               "L. moncytogenes\n(N = 128; L = 3Mb)")) +
  theme_bw(base_size = 18) + 
  theme(legend.key.size = unit(3, 'lines'))
