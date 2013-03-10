COBRA <-
function(train.design,
         train.responses,
         split,
         test,
         machines,
         machines.names,
         grid = 100,
         parallel = FALSE,
         nb.cpus = 2,
         plots = FALSE,
         logs = FALSE,
         progress = TRUE,
         path = "")
  {
    COBRA.knn <- function(x) ## knn
      {
        a <- knn.reg(train=train.design.1,y=train.responses.1,test=x)
        return(a$pred)
      }
    COBRA.lars <- function(x) ## lasso / lars
      {
        a <- lars(train.design.1,train.responses.1,intercept=TRUE)
        res <- predict.lars(a,x,s=length(a$df))$fit
        return(res)
      }
    COBRA.tree <- function(x) ## cart
      {
        a <- tree(as.formula(paste("train.responses.1~",paste("V",c(1:d),sep="",collapse="+"),collapse="",sep="")), data = as.data.frame(train.design.1))
        res <- as.vector(predict(a,as.data.frame(x)))
        return(res)
      }
    COBRA.ridge <- function(x) ## ridge
      {
        a <- linearRidge(as.formula(paste("train.responses.1~",paste("V",c(1:d),sep="",collapse="+"),collapse="",sep="")), data = as.data.frame(train.design.1))
        res <- as.vector(predict(a,as.data.frame(x)))
        return(res)
      }
    COBRA.randomForest <- function(x) ## random forests
      { 
        a <- randomForest(x=train.design.1,y=train.responses.1)
        b <- as.vector(predict(a,x))
        return(b)
      }
    wrapper.COBRA <- function(e)
      {
        poids <- matrix(nrow = n.eval, ncol = n.2, data = 0)
        res <- numeric(n.eval)
        poids <- .C("COBRA",
                    as.integer(n.2),
                    as.integer(n.eval),
                    as.double(test.machines.2),
                    as.double(test.eval),
                    as.integer(n.machines),
                    as.double(alpha),
                    as.double(poids),
                    as.double(e))[[7]]
        dim(poids) <- c(n.eval,n.2)
        poids <- sweep(poids,1,rowSums(poids),FUN="/")
        poids[which(is.nan(poids))] <- 0
        res <- poids%*%train.responses.2
        return(res)
      }
    epsBuild <- function(x)
      {
        return(max(x)-min(x))
      }
    func.risk <- function(x)
      {
        res <- t(train.responses.3 - x)%*%(train.responses.3 - x)/(n.3)
        return(res)
      }
    n.train <- dim(train.design)[1]
    d <- dim(train.design)[2]
    n.test <- dim(test)[1]
    if(missing(machines))
      {
        if(missing(split))
          {
            split <- numeric(2)
            split[1] <- floor(n.train/3)
            split[2] <- n.train - split[1]
            n.2 <- split[2] - split[1]
            n.3 <- n.train - split[2]
          }
        train.design.1 <- train.design[1:split[1],]
        train.design.2 <- train.design[(split[1]+1):split[2],]
        train.design.3 <- train.design[(split[2]+1):n.train,]
        train.responses.1 <- train.responses[1:split[1]]
        train.responses.2 <- train.responses[(split[1]+1):split[2]]
        train.responses.3 <- train.responses[(split[2]+1):n.train]
        if(progress) cat('No machines provided. Loading default machines...\n')
        machines <- c(COBRA.lars,COBRA.ridge,COBRA.knn,COBRA.tree,COBRA.randomForest)
        machines.names <- c('Lars','Ridge','KNN','CART','RF')
        library(FNN)
        library(lars)
        library(tree)
        library(ridge)
        library(randomForest)
        n.machines <- length(machines)
        if(progress) cat('Processing datasets')
        test.machines.2 <- matrix(nrow = n.2, ncol = n.machines, data = 0)
        test.machines.3 <- matrix(nrow = n.3, ncol = n.machines, data = 0)
        test.machines <- matrix(nrow = n.test, ncol = n.machines, data = 0)
        for (imachines in 1:n.machines)
          {
            test.machines.2[,imachines] <- machines[[imachines]](train.design.2)
            test.machines.3[,imachines] <- machines[[imachines]](train.design.3)
            test.machines[,imachines] <- machines[[imachines]](test)
            if(progress) cat('... ',floor(100*imachines/n.machines),'%',sep='')
          }
        cat('.\n')
      } else {
        if(missing(split))
          {
            split <- floor(n.train/2)
            n.2 <- split
            n.3 <- n.train - split
          }
        train.responses.2 <- train.responses[1:split]
        train.responses.3 <- train.responses[(split+1):n.train]
        n.machines <- dim(machines)[2]
        test.machines.2 <- machines[1:split,]
        test.machines.3 <- machines[(split+1):n.train,]
        test.machines <- machines[(n.train+1):n.test,]
      }
    if(progress) cat('Calibrating parameters')
    emin <- 1e-9
    emax <- max(apply(test.machines.2, 2, epsBuild),apply(test.machines.3, 2, epsBuild),apply(test.machines, 2, epsBuild))/2
    estep <- (emax - emin)/grid
    evect <- seq(emin, emax, estep)
    n.e <- length(evect)
    risks <- matrix(nrow = n.e, ncol = n.machines, data = 0)
    alphaseq <- seq(1, n.machines)/n.machines
    epsoptvec <- numeric(n.machines)
    n.eval <- n.3
    test.eval <- test.machines.3
    if(!parallel)
      {
        for(j in alphaseq)
          {
            alpha <- j
            imachines <- j*n.machines
            junk <- sapply(X = evect, FUN = wrapper.COBRA)
            risks[, imachines] <- apply(X = junk, MARGIN = 2, FUN = func.risk)
            epsoptvec[imachines] <- evect[which.min(risks[,imachines])]
            if(progress) cat('... ',floor(100*imachines/n.machines),'%',sep='')
          }
      } else {
        library(snowfall)
        sfInit(parallel = TRUE, cpus = nb.cpus)
        sfClusterEval(dyn.load("COBRA.so"))
        for(j in alphaseq)
          {
            alpha <- j
            imachines <- j*n.machines
            sfExportAll()
            junk <- sfSapply(x = as.matrix(evect), fun = wrapper.COBRA)
            risks[, imachines] <- apply(X = junk, MARGIN = 2, FUN = func.risk)
            epsoptvec[imachines] <- evect[which.min(risks[,imachines])]
            if(progress) cat('... ',floor(100*imachines/n.machines),'%',sep='')
          }
        sfStop()
      }
    alphaopt <- which.min(risks)%/%n.e+1
    alpha <- alphaseq[alphaopt]
    eps <- epsoptvec[alphaopt]
    train.COBRA <- wrapper.COBRA(eps)
    if(progress) cat('.\nRetaining ',alphaopt,' machine(s) out of ',n.machines,', with epsilon = ',eps,'.\n',sep='')
    n.eval <- n.test
    test.eval <- test.machines
    predictor <- wrapper.COBRA(eps)
    res <- list(predict = predictor)
    if(logs)
      {
        risks.machines <- numeric(n.machines + 1)
        risks.machines[1:n.machines] <- apply(X = test.machines.3, MARGIN = 2, FUN = func.risk)
        risks.machines[n.machines + 1] <- func.risk(train.COBRA)
        write.table(t(risks.machines), file = paste(path,"risks.txt",sep=""), append = FALSE, row.names  = FALSE, col.names = FALSE)
      }
    if(plots)
      {
        scenarii <- c("1 machine",paste(2:n.machines,"machines"))
        pdf(paste(path,"epscal.pdf",sep=""))
        plot(evect,risks[,1],ylim=c(0,max(risks)),xlab="Epsilon",ylab="Quadratic Risk",type="l")
        #title("Testing different scenarii to adjust parameters")
        legend(x="bottomright",legend=scenarii,col=1:n.machines,cex=0.8,lty=1,bty="n")
        abline(v = epsoptvec[1],lty=3)
        for(j in 2:n.machines)
          {
            lines(evect,risks[,j],col = j)
            abline(v = epsoptvec[j],col=j,lty=3)
            if(j == alphaopt)
              {
                points(eps,min(risks),pch = 19,col = j,cex = 1.5)
              }
          }
        dev.off()
        pdf(paste(path,"predict-machines.pdf",sep=""))
        lim1 <- min(min(train.responses.3),min(train.COBRA))
        lim2 <- max(max(train.responses.3),max(train.COBRA))
        plot(train.responses.3,train.COBRA,col=1,pch=19,xlab="Responses",ylab="Predictions",xlim=c(lim1,lim2),ylim=c(lim1,lim2))
        #title("Training dataset. Prediction versus true responses")
        abline(0,1,lty=2)
        for(imachines in 1:n.machines)
          {
            points(train.responses.3,test.machines.3[,imachines],col=imachines+1)
          }
        legend(x="bottomright",legend=c('COBRA',machines.names),col=1:(n.machines+1),cex=0.8,pch=c(19,rep(21,n.machines)),bty="n")
        dev.off()
      }
    return(res)
  }
