# intervaldepth
computes the (univariate) interval depth and the (bi- and trivariate) region depths 


# interval depth
x.data=rnorm(1000)
y.data=rnorm(100)
# interval depth of zonoid interval of level 1/2 of y.data ,computed as m1m2(y.data), w.r.t. x.data
m1m2depth(m1m2(y.data),x.data)
# interval depth of zonoid region of w.r.t. x.data
sm1m2depth(y.data,x.data)

# Bivariate region depth
require(MASS)
x.data=mvrnorm(1000,mu=c(0,0),Sigma=matrix(c(1,0,0,1),nrow=2))
y.data=mvrnorm(100,mu=c(0,0),Sigma=matrix(c(1,0,0,1),nrow=2))
# Bivariate region depth of zonoid region of level y.data w.r.t. x.data
m1m2depthr2(y.data,x.data)
