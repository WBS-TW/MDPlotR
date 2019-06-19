

findhalo_rt <- function(mz,
                     intensity,
                     sf = 79/78.917789,
                     step = 0.001,
                     stepsd1=0.003,
                     stepsd2=0.005,
                     mzc=700,
                     cutoffint = 1000,
                     cutoffr=0.4,
                     rt,
                     clustercf =10){
  mzr <- round(mz)
  sm <- mz*sf
  sd <- ceiling(sm)-sm
  smsd <- ifelse(mz<=mzc,stepsd1,stepsd2)
  smstep <- seq(0,1,step)
  rt <- rt
  
  data <- cbind.data.frame(mz=mz,
                           mzr = mzr,
                           sm = sm,
                           sd =sd,
                           intensity=intensity,
                           rt = rt)
  data2 <<- data
  
  result <- NULL
  for(i in 1:length(smstep)){
    maxi = smstep[i]+smsd
    mini = smstep[i]-smsd
    index = sd<maxi & sd>mini
    
    li <- data[index & intensity > cutoffint,]
    mzt <- mzr[index & intensity > cutoffint]
    rtt <- rt[index & intensity > cutoffint]
    
    if(length(mzt)>=2){
      c <- cutree(hclust(dist(mzt)),h=clustercf)
      t <- cutree(hclust(dist(rtt)), h = clustercf)
      u <- paste0(c, t)
      cn <- length(unique(u))
      lit <- cbind.data.frame(li,u,i)
      for (j in 1:cn){
        li2 <- lit[lit[,7]==j,]
        mzt2 <- lit$mzr[lit[,7]==j]
        
        if(length(mzt2)>=2){
          if(length(unique(li2$intensity))>1){
            ratio <- max(li2$intensity[li2$intensity != max(li2$intensity)]) / max(li2$intensity)
            diff <- abs(li2$mzr[round(li2$intensity) == round(max(li2$intensity[li2$intensity != max(li2$intensity)]))] - li2$mzr[which.max(li2$intensity)])
          }else{
            ratio <- 1
            diff <- abs(li2$mzr[1]-li2$mzr[2])
          }
          
          if(ratio>cutoffr&round(diff)==2){
            result <- rbind.data.frame(result,li2)
          }
        }
      }
    }
  }
  return(result[!duplicated(result$mz), ])
}


findhalo_no_rt <- function(mz,
                      intensity,
                      sf = 79/78.917789,
                      step = 0.001,
                      stepsd1=0.003,
                      stepsd2=0.005,
                      mzc=700,cutoffint = 1000,
                      cutoffr=0.4,
                      clustercf = 10){
  mzr <- round(mz)
  sm <- mz*sf
  sd <- ceiling(sm)-sm
  smsd <- ifelse(mz<=mzc,stepsd1,stepsd2)
  smstep <- seq(0,1,step)

  data <- cbind.data.frame(mz=mz, mzr = mzr, sm = sm,sd =sd, intensity=intensity)
  data2 <<- data

  result <- NULL
  for(i in 1:length(smstep)){
    maxi = smstep[i]+smsd
    mini = smstep[i]-smsd
    index = sd<maxi & sd>mini
    
    li <- data[index&intensity>cutoffint,]
    mzt <- mzr[index&intensity>cutoffint]
    
    if(length(mzt)>=2){
      c <- cutree(hclust(dist(mzt)),h= clustercf)
      cn <- length(unique(c))
      lit <- cbind.data.frame(li,c,i)
      for (j in 1:cn){
        li2 <- lit[lit[,6]==j,]
        mzt2 <- lit$mzr[lit[,6]==j]
        
        if(length(mzt2)>=2){
          if(length(unique(li2$intensity))>1){
            ratio <- max(li2$intensity[li2$intensity != max(li2$intensity)]) / max(li2$intensity)
            diff <- abs(li2$mzr[round(li2$intensity) == round(max(li2$intensity[li2$intensity != max(li2$intensity)]))] - li2$mzr[which.max(li2$intensity)])
          }else{
            ratio <- 1
            diff <- abs(li2$mzr[1]-li2$mzr[2])
          }
          
          if(ratio>cutoffr&round(diff)==2){
            result <- rbind.data.frame(result,li2)
          }
        }
      }
    }
  }
  return(result[!duplicated(result$mz), ])
}
