library(plyr)
library(MASS)
library(glmnet) 
library(rms)
library(glinternet)
library(survival)
library(butcher)

rm(list = ls())

info <- sessionInfo()

load(file="Papers/sim9154-sup-0001-supinfo/allfunctions.RData")

## Simulations

ws <- NULL # to include "ws" in ws
ws <- ls()

nsim <- 250

# Grid of simulation scenarios
N <- c(400, 1200, 3600)
TE <- c(0,-1) 
n.main <- 12 
n.int <- c(0,12)
eventRate <- 0.25 
n.noise <- c(0)  

simGrid <- expand.grid(n.main=n.main, TE=TE,
                       n.int=n.int, eventRate=eventRate,
                       N=N, n.noise=n.noise)
simGrid
simGrid$p <- 12
simGrid <- simGrid[order(simGrid$N, simGrid$n.int, simGrid$TE), ]
row.names(simGrid) <- NULL

set.seed(1)
main.coef <- 2^(-seq(0,5.5,.5))
dev1 <- c(runif(9, .9,1.1)*main.coef[1:9] - main.coef[1:9], -2^(c(-1,-2,-3)))

ws <- c(ws, "nsim", "simGrid", "main.coef", "dev1")
rm(list = ls()[!ls() %in% ws])

ptm.overall <- proc.time()

methodList <- c("overall", "hom", "hom.pen", "het", "het.ridge", "het.lasso",  
                "hgl", "het.ck", "risk.mod", "risk.mod.pen", "test.mod.eff", 
                "sessionInfo")

resultsNew <- lapply(methodList, function(x){vector(mode = "list", length = nrow(simGrid))})
names(resultsNew) <- methodList

resultsNew <- lapply(resultsNew, function(x){
  lapply(x, function(xx){
    xx <- vector(mode = "list", length = nsim)})})
rm(methodList)

betaList <- list()

ws <- ls()

i=j=1
for(i in 1:nsim){
  
  set.seed(i)
  
  for(j in 1:nrow(simGrid)){
    
    a <- 0
    
    data <- simData(n=simGrid$N[j],
                    p=simGrid$p[j],
                    rho=0.1, 
                    trt.perc=1/2, 
                    trt.coef = simGrid$TE[j]*log(1/.6),
                    main.coef = main.coef,
                    int.coef = if(simGrid$n.int[j] == 0){
                      rep(0, simGrid$p[j])
                    } else {dev1},
                    eventRate = simGrid$eventRate[j],
                    symm = FALSE, test.set = TRUE, n.test.sets = 1,
                    n.test = 10000) #!#! was 25000 or 20000
    train <- data$train
    test <- data$test
    beta <- data$true.beta
    true.test <- data$true.test
    
    if(i == 1){betaList[[j]] <- beta}
    
    hom.form <- as.formula(paste("y ~ trt +",
                                 paste(grep("trt", grep("v", names(train), value=TRUE), invert=TRUE, value=TRUE),
                                       collapse = "+"), sep=""))
    
    het.form <- as.formula(paste("y ~ trt +", 
                                 paste(grep("v", names(train), value=TRUE), collapse = "+"), sep=""))
    
    risk.mod.form <- as.formula(paste("y ~ ",
                                      paste(grep("trt", grep("v", names(train), value=TRUE), invert=TRUE, value=TRUE),
                                            collapse = "+"), sep=""))
    
    nd <- y01(data=test, trt.var="trt", ctrl.group=0, trt.group=1, 
              main.cols=grep("trt", grep("v", names(test), value=TRUE), invert=TRUE, value=TRUE))
    
    true.y0 <- invlogit(as.matrix(nd$ctrl) %*% beta)
    true.y1 <- invlogit(as.matrix(nd$trt) %*% beta)
    
    betahat <- rep(0, length(beta))
    names(betahat) <- names(beta)
    names(betahat)[1] <- "(Intercept)"
    betahat.template <- betahat
    
    ## overall
    naive <- jhTryCatch(
      overall(train.data=train,outcome="y",
              trt.var="trt", ctrl.group=0, trt.group=1))
    
    if(!is.null(naive$value)){
      betahat <- betahat.template
      coefs <- coef(naive$value$mod)
      betahat[names(coefs)] <- coefs
      
      # for numerical reasons, treat treatment effect < 10e-12 as zero (to prevent val.prob() from failing)
      naive$value$mod$coefficients["trt"] <- ifelse(abs(coef(naive$value$mod)["trt"]) < 0.005, 
                                                    0, coef(naive$value$mod)["trt"])
      pred <- predict(naive$value$mod, newdata = test, type = "response")
      pred.y0 <- predict(naive$value$mod, newdata = nd$ctrl, type = "response")
      pred.y1 <- predict(naive$value$mod, newdata = nd$trt, type = "response")
      resultsNew$overall[[j]][[i]]$risk.diag <- risk.diagnostics(pred=pred, obs=test$y, trt=test$trt, true.test=true.test)
      resultsNew$overall[[j]][[i]]$deltai.diag <- deltai.diagnostics(pred.y0, pred.y1, true.y0, true.y1, pred, test$y, test$trt)
      
      naive$value$mod <- betahat
    }
    
    resultsNew$overall[[j]][[i]]$train <- naive
    
    rm(naive, betahat, coefs, pred, pred.y0, pred.y1)
    
    a <- a+1
    print(c(j,i,a))
    
    ## hom
    hom1 <- jhTryCatch(
      hom(mod.form=hom.form, train.data=train,
          trt.var = "trt",ctrl.group = 0, trt.group = 1, family = "binomial"))
    
    if(!is.null(hom1$value)){
      betahat <- betahat.template
      coefs <- coef(hom1$value$mod)
      betahat[names(coefs)] <- coefs
      
      pred <- predict(hom1$value$mod, newdata = test, type = "response")
      pred.y0 <- predict(hom1$value$mod, newdata = nd$ctrl, type = "response")
      pred.y1 <- predict(hom1$value$mod, newdata = nd$trt, type = "response")
      resultsNew$hom[[j]][[i]]$risk.diag <- risk.diagnostics(pred=pred, obs=test$y, trt=test$trt, true.test=true.test)
      resultsNew$hom[[j]][[i]]$deltai.diag <- deltai.diagnostics(pred.y0, pred.y1, true.y0, true.y1, pred, test$y, test$trt)
      
      hom1$value$mod <- betahat
    }
    
    resultsNew$hom[[j]][[i]]$train <- hom1
    
    rm(hom1, betahat, coefs, pred, pred.y0, pred.y1)
    
    a <- a+1
    print(c(j,i,a))
    
    ## hom.pen
    hom2 <- jhTryCatch(
      hom.pen(mod.form=hom.form, train.data=train,
              trt.var = "trt",ctrl.group = 0, trt.group = 1, family = "binomial",
              penalty.type="ridge",
              unpenalized=NULL,
              alpha.grid=NULL,
              lambda.criterion="lambda.min"))
    
    if(!is.null(hom2$value)){
      betahat <- betahat.template
      coefs <- coef(hom2$value$mod, s=as.numeric(hom2$value$lambda))
      betahat[dimnames(coefs)[[1]]] <- coefs
      
      pred <- predict(hom2$value$mod,
                      newx = as.matrix(test[ ,all.vars(hom.form)[-1]]),
                      s = as.numeric(hom2$value$lambda), type = "response")
      pred.y0 <- predict(hom2$value$mod,
                         newx = as.matrix(nd$ctrl[ ,all.vars(hom.form)[-1]]),
                         s = as.numeric(hom2$value$lambda), type = "response")
      pred.y1 <- predict(hom2$value$mod,
                         newx = as.matrix(nd$trt[ ,all.vars(hom.form)[-1]]),
                         s = as.numeric(hom2$value$lambda), type = "response")
      resultsNew$hom.pen[[j]][[i]]$risk.diag <- risk.diagnostics(pred=pred, obs=test$y, trt=test$trt, true.test=true.test)
      resultsNew$hom.pen[[j]][[i]]$deltai.diag <- deltai.diagnostics(pred.y0, pred.y1, true.y0, true.y1, pred, test$y, test$trt)
      
      hom2$value$mod <- betahat
    }
    
    resultsNew$hom.pen[[j]][[i]]$train <- hom2
    
    rm(hom2, betahat, coefs, pred, pred.y0, pred.y1)
    
    a <- a+1
    print(c(j,i,a))
    
    ## het
    het1 <- jhTryCatch(
      het(mod.form=het.form, train.data=train,
          trt.var = "trt", ctrl.group = 0, trt.group = 1, family = "binomial"))
    
    if(!is.null(het1$value)){
      betahat <- betahat.template
      coefs <- coef(het1$value$mod)
      betahat[names(coefs)] <- coefs
      
      pred <- predict(het1$value$mod, newdata = test, type = "response")
      pred.y0 <- predict(het1$value$mod, newdata = nd$ctrl, type = "response")
      pred.y1 <- predict(het1$value$mod, newdata = nd$trt, type = "response")
      resultsNew$het[[j]][[i]]$risk.diag <- risk.diagnostics(pred=pred, obs=test$y, trt=test$trt, true.test=true.test)
      resultsNew$het[[j]][[i]]$deltai.diag <- deltai.diagnostics(pred.y0, pred.y1, true.y0, true.y1, pred, test$y, test$trt)
      
      het1$value$mod <- betahat
    }
    
    resultsNew$het[[j]][[i]]$train <- het1
    
    rm(het1, betahat, coefs, pred, pred.y0, pred.y1)
    
    a <- a+1
    print(c(j,i,a))
    
    ## het.pen
    het2 <- jhTryCatch(
      het.pen(mod.form=het.form, train.data=train,
              trt.var = "trt",ctrl.group = 0, trt.group = 1, family = "binomial",
              penalty.type="ridge",
              unpenalized=NULL,
              alpha.grid=NULL,
              lambda.criterion="lambda.min"))
    
    if(!is.null(het2$value)){
      betahat <- betahat.template
      coefs <- coef(het2$value$mod, s=as.numeric(het2$value$lambda))
      betahat[dimnames(coefs)[[1]]] <- coefs
      
      pred <- predict(het2$value$mod,
                      newx = as.matrix(test[ ,all.vars(het.form)[-1]]),
                      s = as.numeric(het2$value$lambda), type = "response")
      pred.y0 <- predict(het2$value$mod,
                         newx = as.matrix(nd$ctrl[ ,all.vars(het.form)[-1]]),
                         s = as.numeric(het2$value$lambda), type = "response")
      pred.y1 <- predict(het2$value$mod,
                         newx = as.matrix(nd$trt[ ,all.vars(het.form)[-1]]),
                         s = as.numeric(het2$value$lambda), type = "response")
      resultsNew$het.ridge[[j]][[i]]$risk.diag <- risk.diagnostics(pred=pred, obs=test$y, trt=test$trt, true.test=true.test)
      resultsNew$het.ridge[[j]][[i]]$deltai.diag <- deltai.diagnostics(pred.y0, pred.y1, true.y0, true.y1, pred, test$y, test$trt)
      
      het2$value$mod <- betahat
    }
    
    resultsNew$het.ridge[[j]][[i]]$train <- het2
    
    rm(het2, betahat, coefs, pred, pred.y0, pred.y1)
    
    a <- a+1
    print(c(j,i,a))
    
    ## het.lasso
    het3 <- jhTryCatch(
      het.pen(mod.form=het.form, train.data=train,
              trt.var = "trt",ctrl.group = 0, trt.group = 1, family = "binomial",
              penalty.type="lasso",
              unpenalized=NULL,
              alpha.grid=NULL,
              lambda.criterion="lambda.min"))
    
    if(!is.null(het3$value)){
      betahat <- betahat.template
      coefs <- coef(het3$value$mod, s=as.numeric(het3$value$lambda))
      betahat[dimnames(coefs)[[1]]] <- coefs
      
      pred <- predict(het3$value$mod,
                      newx = as.matrix(test[ ,all.vars(het.form)[-1]]),
                      s = as.numeric(het3$value$lambda), type = "response")
      pred.y0 <- predict(het3$value$mod,
                         newx = as.matrix(nd$ctrl[ ,all.vars(het.form)[-1]]),
                         s = as.numeric(het3$value$lambda), type = "response")
      pred.y1 <- predict(het3$value$mod,
                         newx = as.matrix(nd$trt[ ,all.vars(het.form)[-1]]),
                         s = as.numeric(het3$value$lambda), type = "response")
      resultsNew$het.lasso[[j]][[i]]$risk.diag <- risk.diagnostics(pred=pred, obs=test$y, trt=test$trt, true.test=true.test)
      resultsNew$het.lasso[[j]][[i]]$deltai.diag <- deltai.diagnostics(pred.y0, pred.y1, true.y0, true.y1, pred, test$y, test$trt)
      
      het3$value$mod <- betahat
    }
    
    resultsNew$het.lasso[[j]][[i]]$train <- het3
    
    rm(het3, betahat, coefs, pred, pred.y0, pred.y1)
    
    a <- a+1
    print(c(j,i,a))
    
    ## hgl
    hgl.fit <- jhTryCatch(
      hgl(mod.form=het.form, train.data=train,
          trt.var = "trt", ctrl.group = 0, trt.group = 1,
          family = "binomial"))
    
    if(!is.null(hgl.fit$value)){
      betahat <- betahat.template
      tmp <- ls()
      
      # reparametrization of glinternet solution (see Lim and Hastie 2015)
      index <- with(hgl.fit$value, which(lambda == mod$lambda))
      as <- hgl.fit$value$mod$activeSet[[index]]
      bh <- hgl.fit$value$mod$betahat[[index]]
      pCont <- length(grep("v", all.vars(hom.form)))
      intercept.index <- 1
      if(!is.null(as$cat)){
        cat.indices <- 2:3
        if(!is.null(as$cont)){
          cont.indices <- 4:(4+length(as$cont)-1)
        } else {cont.indices <- NULL}
        if(!is.null(as$catcont)){
          x.number <- nrow(as$catcont)
          x.indices <- matrix(NA, nrow=x.number, ncol=2)
          main.cat.indices <- matrix(NA, nrow=x.number, ncol=2)
          current <- max(cont.indices)
          for(xn in 1:x.number){
            main.cat.indices[xn, ] <- (current + 1):(current + 2)
            x.indices[xn, ] <- (current + 3):(current + 4)
            current <- current+4
          }
        } else {main.cat.indices <- NULL; x.indices <- NULL}
      } else {
        cat.indices <- NULL
        if(!is.null(as$cont)){
          cont.indices <- 2:(2+length(as$cont)-1)
        } else {cont.indices <- NULL}
        if(!is.null(as$catcont)){
          x.number <- nrow(as$catcont)
          x.indices <- matrix(NA, nrow=x.number, ncol=2)
          main.cat.indices <- matrix(NA, nrow=x.number, ncol=2)
          current <- max(cont.indices)
          for(xn in 1:x.number){
            main.cat.indices[xn, ] <- (current + 1):(current + 2)
            x.indices[xn, ] <- (current + 3):(current + 4)
            current <- current+4
          }
        } else {main.cat.indices <- NULL; x.indices <- NULL}
      }
      betahat["(Intercept)"] <- bh[intercept.index] + sum(bh[c(cat.indices[1], main.cat.indices[,1])])
      betahat["trt"] <- sum(bh[c(cat.indices[2], main.cat.indices[,2])]) - sum(bh[c(cat.indices[1], main.cat.indices[,1])])
      betahat[paste0("v", 1:pCont)][as$cont[ ,1]] <- bh[cont.indices]
      betahat[paste0("v", 1:pCont)][as$catcont[ ,2]] <- betahat[paste0("v", 1:pCont)][as$catcont[ ,2]] + bh[x.indices[ ,1]]
      betahat[paste0("v", 1:pCont, "trt")][as$catcont[ ,2]] <- bh[x.indices[ ,2]] - bh[x.indices[ ,1]]
      rm(list = ls()[!ls() %in% tmp])
      
      pred <- predict(hgl.fit$value$mod,
                      X=as.matrix(test[ ,all.vars(hom.form)[-1]]),
                      lambda = hgl.fit$value$lambda, type = "response")
      pred.y0 <- predict(hgl.fit$value$mod,
                         X=as.matrix(nd$ctrl[ ,all.vars(hom.form)[-1]]),
                         lambda = hgl.fit$value$lambda, type = "response")
      pred.y1 <- predict(hgl.fit$value$mod,
                         X=as.matrix(nd$trt[ ,all.vars(hom.form)[-1]]),
                         lambda = hgl.fit$value$lambda, type = "response")
      resultsNew$hgl[[j]][[i]]$risk.diag <- risk.diagnostics(pred=pred, obs=test$y, trt=test$trt, true.test=true.test)
      resultsNew$hgl[[j]][[i]]$deltai.diag <- deltai.diagnostics(pred.y0, pred.y1, true.y0, true.y1, pred, test$y, test$trt)
      
      hgl.fit$value$mod <- betahat
    }
    
    resultsNew$hgl[[j]][[i]]$train <- hgl.fit
    
    rm(hgl.fit, betahat, pred, pred.y0, pred.y1)
    
    a <- a+1
    print(c(j,i,a))
    
    ## Content knowledge 1
    ck.form <- as.formula(paste("y ~ trt +", paste0("v", c(1:12), collapse = "+"),"+", paste0("v", 9:12,"trt", collapse = "+")))
    ck <- jhTryCatch(
      het.pen(mod.form=ck.form, train.data=train,
              trt.var = "trt",ctrl.group = 0, trt.group = 1, family = "binomial",
              penalty.type="ridge",
              unpenalized=NULL,
              alpha.grid=NULL,
              lambda.criterion="lambda.min"))
    
    if(!is.null(ck$value)){
      betahat <- betahat.template
      coefs <- coef(ck$value$mod, s=as.numeric(ck$value$lambda))
      betahat[dimnames(coefs)[[1]]] <- coefs
      
      pred <- predict(ck$value$mod,
                      newx = as.matrix(test[ ,all.vars(ck.form)[-1]]),
                      s = as.numeric(ck$value$lambda), type = "response")
      pred.y0 <- predict(ck$value$mod,
                         newx = as.matrix(nd$ctrl[ ,all.vars(ck.form)[-1]]),
                         s = as.numeric(ck$value$lambda), type = "response")
      pred.y1 <- predict(ck$value$mod,
                         newx = as.matrix(nd$trt[ ,all.vars(ck.form)[-1]]),
                         s = as.numeric(ck$value$lambda), type = "response")
      resultsNew$het.ck[[j]][[i]]$risk.diag <- risk.diagnostics(pred=pred, obs=test$y, trt=test$trt, true.test=true.test)
      resultsNew$het.ck[[j]][[i]]$deltai.diag <- deltai.diagnostics(pred.y0, pred.y1, true.y0, true.y1, pred, test$y, test$trt)
      
      ck$value$mod <- betahat
    }
    
    resultsNew$het.ck[[j]][[i]]$train <- ck
    
    rm(ck, betahat, coefs, pred, pred.y0, pred.y1)
    
    a <- a+1
    print(c(j,i,a))
    
    ## risk.modeling
    rm <- jhTryCatch(
      risk.mod(mod.form=risk.mod.form, train.data=train,
               trt.var="trt", ctrl.group=0,
               trt.group=1, family="binomial"))
    
    if(!is.null(rm$value)){
      betahat.y0mod <- coef(rm$value$y0mod)
      betahat.deltamod <- coef(rm$value$deltamod)
      
      test$lp <- predict(rm$value$y0mod, newdata = test, type = "link")
      pred <- predict(rm$value$deltamod, newdata = test, type = "response")
      test$lp <- NULL
      nd$ctrl$lp <- predict(rm$value$y0mod, newdata = nd$ctrl, type = "link")
      pred.y0 <- predict(rm$value$deltamod, newdata = nd$ctrl, type = "response")
      nd$ctrl$lp <- NULL
      nd$trt$lp <- predict(rm$value$y0mod, newdata = nd$trt, type = "link")
      pred.y1 <- predict(rm$value$deltamod, newdata = nd$trt, type = "response")
      nd$trt$lp <- NULL
      resultsNew$risk.mod[[j]][[i]]$risk.diag <- risk.diagnostics(pred=pred, obs=test$y, trt=test$trt, true.test=true.test)
      resultsNew$risk.mod[[j]][[i]]$deltai.diag <- deltai.diagnostics(pred.y0, pred.y1, true.y0, true.y1, pred, test$y, test$trt)
      
      rm$value$y0mod <- betahat.y0mod
      rm$value$deltamod <- betahat.deltamod
    }
    
    resultsNew$risk.mod[[j]][[i]]$train <- rm
    
    rm(rm, betahat.y0mod, betahat.deltamod, pred, pred.y0, pred.y1)
    
    a <- a+1
    print(c(j,i,a))
    
    # risk.mod.pen
    rmpen <- jhTryCatch(
      risk.mod.pen(mod.form=risk.mod.form, train.data=train,
                   trt.var = "trt",ctrl.group = 0, trt.group = 1,
                   family = "binomial",
                   penalty.type="ridge",
                   unpenalized=NULL,
                   alpha.grid=NULL,
                   lambda.criterion="lambda.min"))
    
    if(!is.null(rmpen$value)){
      coefs <- coef(rmpen$value$y0mod, s=as.numeric(rmpen$value$y0lambda))
      betahat.y0mod <- as.numeric(coefs)
      names(betahat.y0mod) <- dimnames(coefs)[[1]]
      betahat.deltamod <- coef(rmpen$value$deltamod)
      
      test$lp <- predict(rmpen$value$y0mod,
                         newx = as.matrix(test[ ,all.vars(risk.mod.form)[-1]]),
                         s = as.numeric(rmpen$value$y0lambda), type = "link")
      pred <- predict(rmpen$value$deltamod,
                      newdata = test, type="response")
      test$lp <- NULL
      nd$ctrl$lp <- predict(rmpen$value$y0mod,
                            newx = as.matrix(nd$ctrl[ ,all.vars(risk.mod.form)[-1]]),
                            s = as.numeric(rmpen$value$y0lambda), type = "link")
      pred.y0 <- predict(rmpen$value$deltamod,
                         newdata = nd$ctrl, type="response")
      nd$ctrl$lp <- NULL
      nd$trt$lp <- predict(rmpen$value$y0mod,
                           newx = as.matrix(nd$trt[ ,all.vars(risk.mod.form)[-1]]),
                           s = as.numeric(rmpen$value$y0lambda), type = "link")
      pred.y1 <- predict(rmpen$value$deltamod,
                         newdata = nd$trt, type="response")
      nd$trt$lp <- NULL
      resultsNew$risk.mod.pen[[j]][[i]]$risk.diag <- risk.diagnostics(pred=pred, obs=test$y, trt=test$trt, true.test=true.test)
      resultsNew$risk.mod.pen[[j]][[i]]$deltai.diag <- deltai.diagnostics(pred.y0, pred.y1, true.y0, true.y1, pred, test$y, test$trt)
      
      rmpen$value$y0mod <- betahat.y0mod
      rmpen$value$deltamod <- betahat.deltamod
    }
    
    resultsNew$risk.mod.pen[[j]][[i]]$train <- rmpen
    
    rm(rmpen, betahat.y0mod, betahat.deltamod, coefs, pred, pred.y0, pred.y1)
    
    a <- a+1
    print(c(j,i,a))
    
    # test.mod effect
    tm <- jhTryCatch(
      test.mod(hom.form=hom.form, het.form=het.form,
               risk.mod.form=risk.mod.form, train.data=train,
               trt.var="trt", ctrl.group=0,
               trt.group=1, family="binomial", type="effect"))
    
    if(!is.null(tm$value)){
      betahat <- betahat.template
      coefs <- coef(tm$value$mod)
      betahat[names(coefs)] <- coefs
      
      pred <- predict(tm$value$mod, newdata = test, type = "response")
      pred.y0 <- predict(tm$value$mod, newdata = nd$ctrl, type = "response")
      pred.y1 <- predict(tm$value$mod, newdata = nd$trt, type = "response")
      resultsNew$test.mod.eff[[j]][[i]]$risk.diag <- risk.diagnostics(pred=pred, obs=test$y, trt=test$trt, true.test=true.test)
      resultsNew$test.mod.eff[[j]][[i]]$deltai.diag <- deltai.diagnostics(pred.y0, pred.y1, true.y0, true.y1, pred, test$y, test$trt)
      
      tm$value$mod <- betahat
    }
    
    resultsNew$test.mod.eff[[j]][[i]]$train <- tm
    
    rm(tm, betahat, coefs, pred, pred.y0, pred.y1)
    
    a <- a+1
    print(c(j,i,a))
    
    print(c("scenario"=j, "scenario total"=nrow(simGrid), "simulation"=i, "simulation total"=nsim))
    
    ws <- c(ws, "i", "betaList")
    rm(list=ls()[!ls() %in% ws])
  }
}

resultsNew$sessionInfo <- list(info=info, betaList=betaList, time.elapsed=(proc.time() - ptm.overall)["elapsed"])

elapsed <- proc.time() - ptm.overall 
elapsed / nsim # estimated seconds for 1 complete simulation (over all scenarios)

