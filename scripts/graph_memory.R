library(ggplot2)


library(gtools)
df = read.table(snakemake@input[[1]], header=1, stringsAsFactors=FALSE)

## df$Library = as.factor(df$Library, labels=mixedsort(df$Library))
print(mixedsort(unique(df$Library)))
df$Library = factor(df$Library, labels=mixedsort(unique(df$Library)))

p = ggplot(data=df, aes(x=Log10NBIntervals, y=MaxRSSGB, color=Library)) + geom_line() + facet_wrap(~Function) + ggtitle("Memory usage: PyRanges vs. R GenomicRanges") + xlab("Log10 nb intervals") + ylab("GB") + geom_errorbar(aes(x = Log10NBIntervals, ymin = MaxRSSGB - MemorySD, ymax = MaxRSSGB + MemorySD, width=0.05))


if (snakemake@wildcards["subset"] == "_subset"){
  ggsave(snakemake@output[[1]], p, height=5, width=5)

} else {
  ggsave(snakemake@output[[1]], p)
  }
