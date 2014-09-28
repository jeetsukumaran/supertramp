library(adegenet)
summary.df = read.table("processed/summary.txt", header=T)
summary.df = na.omit(summary.df)
groups = summary.df$dispersal.model
cols.to.drop <- c(
                  "dispersal.model",
                  "birth.rate",
                  "death.rate",
                  "dispersal.rate",
                  "niche.evolution.prob",
                  "edges",
                  "est.birth.rate",
                  "length",
                  "size"
                  )
predictors = summary.df[,!(names(summary.df) %in% cols.to.drop)]
dapc.result = dapc(predictors, groups)
model.prefs = data.frame(dispersal.model=groups,
                      pp.model=dapc.result$posterior)
var.contr = data.frame(var=rownames(dapc.result$var.contr),
                       LD1=as.vector(dapc.result$var.contr)
                       )
var.contr = var.contr[order(-var.contr$LD1),]

