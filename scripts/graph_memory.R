library(ggplot2)

df = read.table(snakemake@input[[1]], header=1)

## df$NBIntervals = as.factor(df$NBIntervals)

p = ggplot(data=df, aes(x=NBIntervals, y=MaxRSSGB, color=Library)) + geom_point() + geom_smooth(method=lm) + facet_wrap(~Function) + ggtitle("Memory usage: PyRanges vs. R GenomicRanges") + xlab("Number of intervals in millions") + ylab("GB")

ggsave(snakemake@output[[1]], p)
