# zdepth(x,s) computes the univariate zonoid depth of x w.r.t. sample s
# can handle any discrete distribution with probability vector prob

zdepth <- function(x,s,prob=rep(1/length(s),length(s)),ord=FALSE) {
  prob=prob/sum(prob)
  if(ord==FALSE) {or<-order(s);prob<-prob[or];s<-s[or]}
  n <- length(s)
  if((x<s[1])|(x>s[n])) {return(0)}
  me <- sum(s*prob)
  if(x==me) {return(1)}
  if(x>me) {x <- -x ; s <- -s[n:1] ; prob <- prob[n:1]}
  i <- 1
  ac <- s[i]*prob[i]
  ap <- prob[i]
  while(x*ap>=ac & i<n) {i <- i+1 ; ap <- ap+prob[i] ; ac <- ac+s[i]*prob[i]}
  return((ap*s[i]-ac)/(s[i]-x))
}


# m1m2obj is the objective function to maximize for the interval depth

m1m2obj<-function(m,x,lb){
  n <- length(x)
  if(lb==0 | lb==1) return(0)
  r <- floor(n*lb)+1
  z1 <- zdepth(m[1],x,prob=c(rep(1/(n*lb),r-1),(n*lb-floor(n*lb))/(n*lb),rep(0,n-r)),ord=TRUE)
  z2 <- zdepth(m[2],x,prob=c(rep(0,r-1),(r-n*lb)/(n*(1-lb)),rep(1/(n*(1-lb)),n-r)),ord=TRUE)
  return(2*min(lb*z1,(1-lb)*z2))
}


# m1m2depth(m,x) computes the interval depth of m w.r.t. sample x

m1m2depth<-function(m,x,ord=FALSE) {
  if(m[2]<m[1]) return(0)
  if(ord==FALSE) x <- sort(x)
  lo <- mean(x<m[1])
  up <- mean(x<=m[2])
  if(lo==1 | up==0 | up==lo) {return(0)}
  return(optimize(m1m2obj,m=m,x=x,lower=lo,upper=up,maximum=TRUE)$objective)
}


# m1m2(samp) computes the zonoid interval of level 1/2 of sample samp

m1m2<-function(samp) {
  mu<-mean(samp)
  desv<-mean(abs(samp-median(samp)))
  return(c(mu-desv,mu+desv))
}

# sm1m2depth(samp,x) computes the interval depth of the zonoid interval 
# of level 1/2 of sample samp w.r.t. sample x

sm1m2depth<-function(samp,x,ord=FALSE) {
  return(m1m2depth(m=m1m2(samp),x,ord))
}


# objective function of 2-dimensional region depth

m1m2r2obj<- function(ang,samp,X) {
  return(sm1m2depth(samp=cos(ang)*samp[,1]+sin(ang)*samp[,2],
                   x=cos(ang)*X[,1]+sin(ang)*X[,2]))
}


# 2-dimensional region depth of zonoid region of level 1/2 of samp wrt X

m1m2depthr2<-function(samp,X,nang=10) {
  depth <- NULL
  for(i in 1:nang) {
    depth<- min(depth,
                optimize(m1m2r2obj,samp=samp,X=X,lower=2*pi*(i-1)/nang,upper=2*pi*i/nang,maximum=FALSE,tol=.Machine$double.eps^0.25)$objective)
  }
  return(depth)
}

# objective function of 3-dimensional region depth

m1m2r3obj <- function(ang,samp,X) {
  return(sm1m2depth(samp=sin(ang[1])*cos(ang[2])*samp[,1]+sin(ang[1])*sin(ang[2])*samp[,2]+cos(ang[1])*samp[,3],
                   x=sin(ang[1])*cos(ang[2])*X[,1]+sin(ang[1])*sin(ang[2])*X[,2]+cos(ang[1])*X[,3]))
}


# 3-dimensional region depth of zonoid region of level 1/2 of samp wrt X

m1m2depthr3<-function(samp,X,nang=10) {
  depth <- 1
  for(th in pi*seq(from=0,to=1,length.out=nang)[-1]) {
    for(phi in 2*pi*seq(from=0,to=1,length.out=nang)[-1]) {
      depth <- min(depth,optim(c(th,phi),m1m2r3obj,samp=samp,X=X)$value)
    }
  }
  return(depth)
}

# extreme points of the zonoid region of level 1/2 of the bivariate dataset data

zonoidreg <- function(data,nlines=1000){
  n <- nrow(data)
  angle <- seq(0,2*pi,length.out=nlines) 
  extreme.points <- NULL
  zon1 <- mean(sort(data[,1]*cos(angle[1])+data[,2]*sin(angle[1]))[(n/2+1):n])
  
  for(i in 2:nlines){
    zon2 <-  mean(sort(data[,1]*cos(angle[i])+data[,2]*sin(angle[i]))[(n/2+1):n])
    extreme.points <- rbind(extreme.points,solve(matrix(c(cos(angle[i-1]),sin(angle[i-1]),
                                                          cos(angle[i]),sin(angle[i])),byrow=TRUE,ncol=2),
                                                 c(zon1,zon2)))
    zon1 <- zon2
  }
  extreme.points <- extreme.points[chull(extreme.points),]
  extreme.points <- rbind(extreme.points,extreme.points[1,])
  return(extreme.points)
}
