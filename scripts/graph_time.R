library(ggplot2)

df = read.table(snakemake@input[[1]], header=1)

## df$NBIntervals = as.factor(df$NBIntervals)
print(head(df))

p = ggplot(data=df, aes(x=Log10NBIntervals, y=Seconds, color=Library)) + geom_point() +  geom_line() + facet_wrap(~Function) + ggtitle("Time usage: PyRanges vs. R GenomicRanges") + xlab("Number of intervals in millions") + ylab("Time in seconds")

ggsave(snakemake@output[[1]], p)
