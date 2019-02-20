library(ggpubr)
library(ggplot2)
library(gridExtra)

df = read.table(snakemake@input[[1]], header=1, stringsAsFactors=FALSE, row.names=NULL)
print(head(df))

xlab = ""
bottom = "Log10 number of intervals."
desc = snakemake@params[["description"]]
func = snakemake@params[["function"]]

ylab_time = "Log10 seconds"
ylab_memory = "GB Memory"

print(desc)
print(func)

f = ggplot(df, aes(x=Log10NBIntervals, y=Seconds)) + geom_line(aes(colour=Library)) + facet_wrap(~Filetype, ncol=1) + geom_errorbar(aes(x=Log10NBIntervals, ymin=Seconds - SecondsSD, ymax=Seconds + SecondsSD, width=0.05)) + xlab(xlab) + ylab("") + theme(plot.title = element_text(size=8), axis.text.x = element_text(angle = 45, hjust = 1, vjust=1)) + ggtitle(ylab_time)
                                        # theme(legend.position="none")
f$layout$l[f$layout$name == "title"] <- 0

g = ggplot(df, aes(x=Log10NBIntervals, y=MaxRSSGB)) + geom_line(aes(colour=Library)) + facet_wrap(~Filetype, ncol=1) + geom_errorbar(aes(x=Log10NBIntervals, ymin=MaxRSSGB - MemorySD, ymax=MaxRSSGB + MemorySD, width=0.05)) + xlab(xlab) + ylab("") + theme(plot.title = element_text(size=8), axis.text.x = element_text(angle = 45, hjust = 1, vjust=1)) + ggtitle(ylab_memory)

g$layout$l[g$layout$name == "title"] <- 1


fig = ggarrange(f, g, ncol=2, common.legend=TRUE, legend="right")

fig = annotate_figure(fig, bottom=text_grob(bottom, vjust=-2, size=8))

fig = grid.arrange(fig, top=text_grob(desc, size=10, vjust=0.66))

fig = grid.arrange(fig, top=text_grob(func, size=18))

## desc_grob = text_grob()

## fig = arrangeGrob(fig, desc_grob)

ggsave(snakemake@output[[1]], height=4, width=6, fig)
