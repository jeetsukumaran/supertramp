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

analyze.dapc = function(n.pc, n.da) {
    dapc.result = dapc(predictors, groups, n.pc=n.pc, n.da=n.da)
    var.contr = data.frame(var=rownames(dapc.result$var.contr),
                        LD1=as.vector(dapc.result$var.contr)
                        )
    var.contr = var.contr[order(-var.contr$LD1),]
    model.prefs = data.frame(dispersal.model=groups,
                        pp.model=dapc.result$posterior)
    mean.pp.for.correct.model = mean(c(model.prefs[model.prefs$dispersal.model == "constrained",]$pp.model.constrained,
                                    model.prefs[model.prefs$dispersal.model == "unconstrained",]$pp.model.unconstrained))
    debug1 = nrow(model.prefs[model.prefs$dispersal.model=="constrained" & model.prefs$pp.model.constrained>0.5,])
    debug2 = nrow(model.prefs[model.prefs$dispersal.model=="unconstrained" & model.prefs$pp.model.unconstrained>0.5,])
    mean.count.correct.model.preferred = nrow(model.prefs[model.prefs$dispersal.model=="constrained" & model.prefs$pp.model.constrained>0.5,]) + nrow(model.prefs[model.prefs$dispersal.model=="unconstrained" & model.prefs$pp.model.unconstrained>0.5,])
    mean.prop.correct.model.preferred = mean.count.correct.model.preferred / nrow(model.prefs)
    rv = list(
              dapc.result=dapc.result,
              var.contr=var.contr,
              model.prefs=model.prefs,
              mean.pp.for.correct.model=mean.pp.for.correct.model,
              mean.count.correct.model.preferred=mean.count.correct.model.preferred,
              mean.prop.correct.model.preferred=mean.prop.correct.model.preferred,
              debug1=debug1,
              debug2=debug2
              )
    rv
}

assess.vars = function() {
    result = data.frame()
    for (n.pc in 2:nrow(predictors)) {
        for (n.da in 1:nrow(predictors)) {
            x = analyze.dapc(n.pc, n.da)
            cat(paste(n.pc,
                        n.da,
                        x$mean.pp.for.correct.model,
                        x$mean.prop.correct.model.preferred,
                        "\n",
                        sep="\t\t"
                        ))
            subresult = data.frame(
                        n.pc = n.pc,
                        n.da = n.da,
                        mean.pp.for.correct.model=x$mean.pp.for.correct.model,
                        mean.prop.correct.model.preferred=x$mean.prop.correct.model.preferred
                        )
            result = rbind(result, subresult)
        }
    }
    result
}
