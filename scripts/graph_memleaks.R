library(ggplot2)

df = read.table(snakemake@input[[1]], header=1)

df$Iterations = as.factor(df$Iterations)
df = df[df$LeakType != "Indirectly",] # Indirectly was always zero, so is on top of Definitely

p = ggplot(data=df, aes(x=Iterations, y=MegaBytes, color=LeakType)) + geom_point() + facet_wrap(~Function) + ggtitle("MemoryLeaks PyRanges") + xlab("Number of iterations") + ylab("Megabytes")

ggsave(snakemake@output[[1]], p)
