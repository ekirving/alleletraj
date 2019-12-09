#!/usr/bin/env Rscript
require(ggplot2, quietly=T, warn.conflicts=F)
require(gganimate, quietly=T, warn.conflicts=F)
require(R.utils, quietly=T, warn.conflicts=F)
require(data.table, quietly=T, warn.conflicts=F)
require(reshape2, quietly=T, warn.conflicts=F)
require(viridis, quietly = T, warn.conflicts=F)
require(ggExtra, quietly = T, warn.conflicts=F)

# get the command line arguments
args <- commandArgs(trailingOnly = TRUE)
output_prefix <- args[1]

# TODO remove when done testing
# output_prefix <- 'test-chain1'

time_file <- paste0(output_prefix, '.time.gz')
traj_file <- paste0(output_prefix, '.traj.gz')
param_file <- paste0(output_prefix, '.param.gz')

# find the length of the longest trajectory
max_length <- max(count.fields(time_file, sep = " "))

# load the time and trajectory data (and set the number of columns)
time <- read.table(time_file, sep=" ", fill=TRUE, header = F, col.names=c("gen", 1:(max_length-1)))
traj <- read.table(traj_file, sep=" ", fill=TRUE, header = F, col.names=c("gen", 1:(max_length-1)))

# melt the data using "gen" as the key and join the two dataframes
path <- merge(
  melt(time, id.vars = "gen", value.name = "time", na.rm=T), 
  melt(traj, id.vars = "gen", value.name = "traj", na.rm=T), 
  by = c('gen', 'variable'))

# load the other parameter samples
mcmc.params <- read.table(param_file, sep="\t", header=T)

# untransform the end frequency
mcmc.params$end_freq <- (1 - cos(mcmc.params$end_freq)) / 2

p <- ggplot() +
  geom_line(data=path, aes(x=time, y=traj, group=gen, color=gen)) +
  scale_colour_viridis(direction = -1, option="viridis") +
  geom_point(data=mcmc.params, aes(x=age, y=0, fill=gen), shape=21) +
  geom_point(data=mcmc.params, aes(x=0, y=end_freq, fill=gen), shape=21) +
  # scale_fill_viridis(direction = -1, option="magma") +
  scale_fill_gradient(low = "#ffffff", high = "#ff0000") +
  labs(title = '', x = 'Time', y = 'Frequency', color = "Path", fill = "Age") +
  theme_bw()

png(filename = paste0(output_prefix, '.png'), width=14, height=7, units='in', res=300)

# marginal density of age
ggMarginal(p, margins = 'x', x=age, y=0, type="density")

dev.off()

# # animated line plot with colours and fade
# a <- ggplot(data=data.merged, aes(x=time, y=traj, group=gen, alpha=gen, color=gen)) +
#   geom_line() +
#   transition_time(gen) +
#   scale_colour_viridis(direction = -1) +
#   labs(title = 'Generations: {frame_time} -', x = 'Time', y = 'Frequency') +
#   theme_bw()
# 
# animate(a, fps=1)
