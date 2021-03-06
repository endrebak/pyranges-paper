library(ggplot2)
library(gtools)


df = read.table(snakemake@input[[1]], header=1, stringsAsFactors=FALSE)

## df$NBIntervals = as.factor(df$NBIntervals)
## df$Library = as.factor(df$Library, mixedsort(df$Library))
df$Library = factor(df$Library, levels=mixedsort(unique(df$Library)))
print(head(df))
print(tail(df))


if (snakemake@params[["subset"]]){
  ncol = NULL
} else {
  ncol = NULL
}

p = ggplot(data=df, aes(x=Log10NBIntervals, y=Seconds, color=Library))  + geom_line() + facet_wrap(~Function, ncol=ncol) + labs(title = snakemake@params[["title"]]) + xlab("Log10 nb intervals") + ylab("Log10 seconds") + geom_errorbar(aes(x=Log10NBIntervals, ymin=Seconds - SecondsSD, ymax=Seconds + SecondsSD, width=0.05)) # , subtitle=snakemake@params[["subtitle"]]

ggsave(snakemake@output[[1]], p)
