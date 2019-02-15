library(ggpubr)
library(ggplot2)
library(gridExtra)

df = read.table(snakemake@input[[1]], header=1, stringsAsFactors=FALSE)
print(head(df))

f = ggplot(df, aes(x=Log10NBIntervals, y=Seconds)) + geom_line(aes(colour=Library)) + facet_wrap(~Function) + geom_errorbar(aes(x=Log10NBIntervals, ymin=Seconds - SecondsSD, ymax=Seconds + SecondsSD, width=0.05)) + xlab("") + ylab("Running time in log10 seconds.") + theme(axis.text.x = element_text(angle = 45, hjust = 1, vjust=1)) + ggtitle("Time usage")
                                        # theme(legend.position="none")
f$layout$l[f$layout$name == "title"] <- 0

g = ggplot(df, aes(x=Log10NBIntervals, y=MaxRSSGB)) + geom_line(aes(colour=Library)) + facet_wrap(~Function) + geom_errorbar(aes(x=Log10NBIntervals, ymin=MaxRSSGB - MemorySD, ymax=MaxRSSGB + MemorySD, width=0.05)) + xlab("") + ylab("MaxRSS use in GB.") + theme(axis.text.x = element_text(angle = 45, hjust = 1, vjust=1)) + ggtitle("Memory usage")
g$layout$l[g$layout$name == "title"] <- 1

fig = ggarrange(f, g, ncol=2, common.legend=TRUE, legend="right")

fig = annotate_figure(fig, bottom=text_grob("Log10 number of intervals."))

ggsave(snakemake@output[[1]], height=4, width=6)
