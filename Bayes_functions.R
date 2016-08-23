   # Fitted parameter for the bayesian framework
BAYESIAN_FITTED<-c(0.407277142798302, 0.554007336744485, 0.63777155771234, 0.693989162719009, 0.735450014674917, 0.767972534429806, 0.794557287143399, 0.816906816601605, 0.83606796225341, 0.852729446430296, 0.867370424541641, 0.880339760590323, 0.891900995024999, 0.902259181289864, 0.911577919359,0.919990301665853, 0.927606458124537, 0.934518806350661, 0.940805863754375, 0.946534836475715, 0.951763691199255, 0.95654428191308, 0.960920179487397, 0.964930893680829, 0.968611312149038, 0.971992459313836, 0.975102110004818, 0.977964943023096, 0.980603428208439, 0.983037660179428, 0.985285800977406, 0.987364285326685, 0.989288037855441, 0.991070478823525, 0.992723699729969, 0.994259575477392, 0.995687688867975, 0.997017365051493, 0.998257085153047, 0.999414558305388, 1.00049681357804, 1.00151036237481, 1.00246080204981, 1.00335370751909, 1.0041939329768, 1.0049859393417, 1.00573382091263, 1.00644127217376, 1.00711179729107, 1.00774845526417, 1.00835412715854, 1.00893143010366, 1.00948275846309, 1.01001030293661, 1.01051606798079, 1.01100188771288, 1.01146944044216, 1.01192026195449, 1.01235575766094, 1.01277721370986)
 
  CONST_i <- sort(c(((2^(seq(-39,0,length.out=201)))/2)[1:200],(c(0:11,13:99)+0.5)/100,1-(2^(seq(-39,0,length.out=201)))/2))
  
# Given x, M & p, returns a pdf 
   calculate_bayes <- function ( x=3, N=10, p=0.33,
                                i=CONST_i,
                                max_sigma=20,length_sigma=4001
                              ){
    if(!0%in%N){
      G <- max(length(x),length(N),length(p))
      x=array(x,dim=G)
      N=array(N,dim=G)
      p=array(p,dim=G)
      sigma_s<-seq(-max_sigma,max_sigma,length.out=length_sigma)
      sigma_1<-log({i/{1-i}}/{p/{1-p}})
      index<-min(N,60)
      y<-dbeta(i,x+BAYESIAN_FITTED[index],N+BAYESIAN_FITTED[index]-x)*(1-p)*p*exp(sigma_1)/({1-p}^2+2*p*{1-p}*exp(sigma_1)+{p^2}*exp(2*sigma_1))
      y<-y/mean(y)
      if(!is.na(unique(y))){
        tmp<-approx(sigma_1,y,sigma_s)$y
        tmp/sum(tmp)/{2*max_sigma/{length_sigma-1}}
      }else{
        return(NA)
      }
    }else{
      return(NA)
    }
  }  

   
  # Computes the 95% CI for a pdf
  calcBayesCI <- function(Pdf,low=0.025,up=0.975,max_sigma=20, length_sigma=4001){
    if(length(Pdf)!=length_sigma) return(NA)
    sigma_s=seq(-max_sigma,max_sigma,length.out=length_sigma)
    cdf = cumsum(Pdf)
    cdf = cdf/cdf[length(cdf)]  
    return( c(sigma_s[findInterval(low,cdf)-1] , sigma_s[findInterval(up,cdf)]) ) 
  }
  
  # Computes a mean for a pdf
  calcBayesMean <- function(Pdf,max_sigma=20,length_sigma=4001){
    if(length(Pdf)!=length_sigma) return(NA)
    Pdf=Pdf/sum(Pdf)
    sigma_s=seq(-max_sigma,max_sigma,length.out=length_sigma)
 #   norm = {length_sigma-1}/2/max_sigma
    return( (Pdf%*%sigma_s)  ) 
  }
  
  # Returns the mean, and the 95% CI for a pdf
  calcBayesOutputInfo <- function(Pdf,low=0.025,up=0.975,max_sigma=20, length_sigma=4001){
    if(is.na(Pdf)) 
     return(rep(NA,3))  
    bCI = calcBayesCI(Pdf=Pdf,low=low,up=up,max_sigma=max_sigma,length_sigma=length_sigma)
    bMean = calcBayesMean(Pdf=Pdf,max_sigma=max_sigma,length_sigma=length_sigma)
    return(c(bMean, bCI))
  }   
  calcBayesSD <- function(Pdf,max_sigma=20,length_sigma=4001){
    if(length(Pdf)!=length_sigma) return(NA)
    sigma_s=seq(-max_sigma,max_sigma,length.out=length_sigma)
    norm = {length_sigma-1}/2/max_sigma
    return( sqrt((((Pdf)%*%((sigma_s^2/norm)))-(Pdf%*%sigma_s/norm)^2)  ) )
  }

  # Computes the p-value of a pdf
  computeSigmaP <- function(Pdf, length_sigma=4001, max_sigma=20){
    if(length(Pdf)>1){
    Pdf<-Pdf/sum(Pdf)
      pVal = {sum(Pdf[1:{{length_sigma-1}/2}]) + Pdf[{{length_sigma+1}/2}]/2}
      if(pVal>0.5){
        pVal = pVal-1
      }
      return(pVal)
    }else{
      return(NA)
    }
  }    
  
  # Compute p-value of two distributions
  compareTwoDistsFaster <-function(sigma_S=seq(-20,20,length.out=4001), N=10000, dens1=runif(4001,0,1), dens2=runif(4001,0,1)){
  #print(c(length(dens1),length(dens2)))
  if(length(dens1)>1 & length(dens2)>1 ){
    dens1<-dens1/sum(dens1)
    dens2<-dens2/sum(dens2)
    cum2 <- cumsum(dens2)-dens2/2
    tmp<- sum(sapply(1:length(dens1),function(i)return(dens1[i]*cum2[i])))
    #print(tmp)
    if(tmp>0.5)tmp<-tmp-1
    return( tmp )
    }
    else {
    return(NA)
    }
    #return (sum(sapply(1:N,function(i)(sample(sigma_S,1,prob=dens1)>sample(sigma_S,1,prob=dens2))))/N)
  }  
  
 length_sigma<-4001
 max_sigma<-20

sigma_s<-seq(-max_sigma,max_sigma,length.out=length_sigma)


  ##Covolution
  break2chunks<-function(G=1000){
  base<-2^round(log(sqrt(G),2),0)
  return(c(rep(base,floor(G/base)-1),base+G-(floor(G/base)*base)))
  }  
  
  PowersOfTwo <- function(G=100){
    exponents <- array()
    i = 0
    while(G > 0){
      i=i+1
      exponents[i] <- floor( log2(G) )
      G <- G-2^exponents[i]
    }
    return(exponents)
  }
  
  convolutionPowersOfTwo <- function( cons, length_sigma=4001 ){
    G = ncol(cons)
    if(G>1){
      for(gen in log(G,2):1){
        ll<-seq(from=2,to=2^gen,by=2)
        sapply(ll,function(l){cons[,l/2]<<-weighted_conv(cons[,l],cons[,l-1],length_sigma=length_sigma)})
      }
    }
    return( cons[,1] )
  }
  
  convolutionPowersOfTwoByTwos <- function( cons, length_sigma=4001,G=1 ){
    if(length(ncol(cons))) G<-ncol(cons)
    groups <- PowersOfTwo(G)
    matG <- matrix(NA, ncol=length(groups), nrow=length(cons)/G )
    startIndex = 1
    for( i in 1:length(groups) ){
      stopIndex <- 2^groups[i] + startIndex - 1
      if(stopIndex!=startIndex){
        matG[,i] <- convolutionPowersOfTwo( cons[,startIndex:stopIndex], length_sigma=length_sigma )
        startIndex = stopIndex + 1
      }
      else {
        if(G>1) matG[,i] <- cons[,startIndex:stopIndex]
        else matG[,i] <- cons
        #startIndex = stopIndex + 1
      }
    }
    return( list( matG, groups ) )
  }
  
  weighted_conv<-function(x,y,w=1,m=100,length_sigma=4001){
    lx<-length(x)
    ly<-length(y)
    if({lx<m}| {{lx*w}<m}| {{ly}<m}| {{ly*w}<m}){
      if(w<1){
        y1<-approx(1:ly,y,seq(1,ly,length.out=m))$y
        x1<-approx(1:lx,x,seq(1,lx,length.out=m/w))$y
        lx<-length(x1)
        ly<-length(y1)
      }
      else {
        y1<-approx(1:ly,y,seq(1,ly,length.out=m*w))$y
        x1<-approx(1:lx,x,seq(1,lx,length.out=m))$y
        lx<-length(x1)
        ly<-length(y1)
      }
    }
    else{
      x1<-x
      y1<-approx(1:ly,y,seq(1,ly,length.out=floor(lx*w)))$y
      ly<-length(y1)
    }
    tmp<-approx(x=1:(lx+ly-1),y=convolve(x1,rev(y1),type="open"),xout=seq(1,lx+ly-1,length.out=length_sigma))$y
    tmp[tmp<=0] = 0
    return(tmp/sum(tmp))
  }
  
  calculate_bayesGHelper <- function( listMatG,length_sigma=4001 ){
    matG <- listMatG[[1]]
    groups <- listMatG[[2]]
    i = 1
    resConv <- matG[,i]
    denom <- 2^groups[i]
    if(length(groups)>1){
      while( i<length(groups) ){
        i = i + 1
        resConv <- weighted_conv(resConv, matG[,i], w= {{2^groups[i]}/denom} ,length_sigma=length_sigma)
        #cat({{2^groups[i]}/denom},"\n")
        denom <- denom + 2^groups[i]
      }
    }
    return(resConv)
  }
  
  # Given a list of PDFs, returns a convoluted PDF    
  groupPosteriors <- function( listPosteriors, max_sigma=20, length_sigma=4001 ,Threshold=2 ){  
    listPosteriors = listPosteriors[ !is.na(listPosteriors) ]
    Length_Postrior<-length(listPosteriors)
    if(Length_Postrior>1 & Length_Postrior<=Threshold){
      cons = matrix(unlist(listPosteriors),length(listPosteriors[[1]]),length(listPosteriors))
      listMatG <- convolutionPowersOfTwoByTwos(cons,length_sigma=length_sigma)
      y<-calculate_bayesGHelper(listMatG,length_sigma=length_sigma)
      return( y/sum(y)/(2*max_sigma/(length_sigma-1)) )
    }else if(Length_Postrior==1) return(listPosteriors[[1]])
    else  if(Length_Postrior==0) return(NA)
    else {
      cons = matrix(unlist(listPosteriors),length(listPosteriors[[1]]),length(listPosteriors))
      y = fastConv(cons,max_sigma=max_sigma, length_sigma=length_sigma )
      return( y/sum(y)/(2*max_sigma/(length_sigma-1)) )
    }
  }

  fastConv<-function(cons, max_sigma=20, length_sigma=4001){
    chunks<-break2chunks(G=ncol(cons))
    if(ncol(cons)==3) chunks<-2:1
    index_chunks_end <- cumsum(chunks)
    index_chunks_start <- c(1,index_chunks_end[-length(index_chunks_end)]+1)
    index_chunks <- cbind(index_chunks_start,index_chunks_end)
    
    case <- sum(chunks!=chunks[1])
    if(case==1) End <- max(1,((length(index_chunks)/2)-1))
    else End <- max(1,((length(index_chunks)/2)))
    
    firsts <- sapply(1:End,function(i){
          	    indexes<-index_chunks[i,1]:index_chunks[i,2]
          	    convolutionPowersOfTwoByTwos(cons[ ,indexes])[[1]]
          	  })
    if(case==0){
    	result<-calculate_bayesGHelper( convolutionPowersOfTwoByTwos(firsts) )
    }else if(case==1){
      last<-list(calculate_bayesGHelper(
      convolutionPowersOfTwoByTwos( cons[ ,index_chunks[length(index_chunks)/2,1]:index_chunks[length(index_chunks)/2,2]] )
                                      ),0)
      result_first<-calculate_bayesGHelper(convolutionPowersOfTwoByTwos(firsts))
      result<-calculate_bayesGHelper(
        list(
          cbind(
          result_first,last[[1]]),
          c(log(index_chunks_end[length(index_chunks)/2-1],2),log(index_chunks[length(index_chunks)/2,2]-index_chunks[length(index_chunks)/2,1]+1,2))
        )
      )
    }
    return(as.vector(result))
  }
    

 length_sigma<-4001
 max_sigma<-20

 sigma_s<-seq(-max_sigma,max_sigma,length.out=length_sigma)




