library(ggplot2)
library(gtools)


df = read.table(snakemake@input[[1]], header=1, stringsAsFactors=FALSE)

## df$NBIntervals = as.factor(df$NBIntervals)
## df$Library = as.factor(df$Library, mixedsort(df$Library))
df$Library = factor(df$Library, levels=mixedsort(unique(df$Library)))
print(head(df))
print(tail(df))


if (snakemake@params[["subset"]]){
  ncol = 4
} else {
  ncol = NULL
}

p = ggplot(data=df, aes(x=Log10NBIntervals, y=Seconds, color=Library))  + geom_line() + facet_wrap(~Function, ncol=ncol) + ggtitle("Time usage: PyRanges vs. R GenomicRanges") + xlab("Log10 nb intervals") + ylab("Log10 seconds") + geom_errorbar(aes(x=Log10NBIntervals, ymin= Seconds - SecondsSD, ymax=Seconds + SecondsSD, width=0.05))

ggsave(snakemake@output[[1]], p)
