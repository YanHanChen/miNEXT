
DetAbu = function(x, zero=FALSE){
  x <- unlist(x)
  n <- sum(x)
  f1 <- sum(x==1)
  f2 <- sum(x==2)
  f3 <- sum(x==3)
  if(f2==0){
    f1 <- max(f1 - 1, 0)
    f2 <- 1
  }
  A1 <- f1 / n * ((n-1)*f1 / ((n-1)*f1 + 2*max(f2,1)))
  A2 <- f2 / choose(n, 2) * ((n-2)*f2 / ((n-2)*f2 + 3*max(f3,1)))^2
  if(zero==FALSE) x <- x[x>0]
  q.solve <- function(q){
    e <- A1 / sum(x/n*exp(-q*x))
    out <- sum((x/n * (1 - e * exp(-q*x)))^2) - sum(choose(x,2)/choose(n,2)) + A2
    abs(out)
  }
  q <- tryCatch(optimize(q.solve, c(0,1))$min, error = function(e) {1})
  e <- A1 / sum(x/n*exp(-q*x))
  o <- x/n * (1 - e * exp(-q*x))
  o
}

#
###########################################
#' Estimating undetected species relative abundance
#' 
#' \code{UndAbu} Estimating undetected species relative abundance
#' @param x a vector of species abundance frequency
#' @return a numerical vector
UndAbu <- function(x){
  x <- unlist(x)
  n <- sum(x)  
  f1 <- sum(x==1)
  f2 <- sum(x==2)
  f3 <- sum(x==3)
  f4 <- max(sum(x == 4), 1)
  f0.hat <- ceiling(ifelse(f2 == 0, (n - 1) / n * f1 * (f1 - 1) / 2, (n - 1) / n * f1 ^ 2/ 2 / f2))  #estimation of unseen species via Chao1
  if(f0.hat < 0.5){
    return(NULL)
  } 
  
  if(f1>0 & f2>0){
    A1 <- f1 / n * ((n-1)*f1 / ((n-1)*f1 + 2*f2))    
  }else if(f1>1 & f2==0){
    A1 <- (f1-1) / n * ((n-1)*f1 / ((n-1)*f1 + 2))
  }else{
    return(NULL)
  }
  
  if(f2>0 & f3>0){
    A2 <- f2 / choose(n, 2) * ((n-2)*f2 / ((n-2)*f2 + 3*f3))^2  
  }else if(f2>1 & f3==0){
    A2 <- (f2-1) / choose(n, 2) * ((n-2)*f2 / ((n-2)*f2 + 3))^2  
  }else{
    A2 <- 0
  }
  
  
  R <- ifelse(A2>0, A1^2/A2, 1)
  j <- 1:f0.hat
  f.solve <- function(x){ 
    out <- sum(x^j)^2 / sum((x^j)^2) - R
    abs(out)
  }
  b <-  tryCatch(optimize(f.solve, lower=(R-1)/(R+1), upper=1, tol=1e-5)$min, error = function(e) {(R-1)/(R+1)})
  a <- A1 / sum(b^j)
  p <- a * b^j
  if(f0.hat == 1) p <- A1
  p
}



Dq_in<-function(x1,x2,mm,q){
  out = apply(mm,1,function(mm) {
    mm1 = mm[1]
    mm2 = mm[2]
    D_share(x1, x2, mm1, mm2, q)
  })
}



Dq0_2_ext_cpp<-function(p1_hat,p2_hat,n1,n2,mmext){
  out = apply(mmext,1,function(mm) {
    mm1 = mm[1]
    mm2s = mm[2]
    h0hat<-h0_hat_cpp(p1_hat, p2_hat, mm1, mm2s,n1,n2)
  })
}


Dq0_mix_ext_cpp<-function(p1_hat,p2_hat,n1,n2,mmext){
  out = apply(mmext,1,function(mm) {
    mm1 = mm[1]
    mm2s = mm[2]
    h0hat<-h0_hat_cpp(p1_hat, p2_hat, mm1, mm2s,n1,n2)
  })
}

h0hat<-function( pi1, pi2, m1,  m2s, n2 ){
  output = 0
  output = (1-pi1)^m1*(1-pi2)^n2*(1-(1-pi2)^m2s)
  return (output)
}


Dq0_2_ext_old<-function(x1,x2,p1_hat,p2_hat,n1,n2,mmext){
  out = apply(mmext,1,function(mm) {
    mm1 = mm[1]
    mm2s = mm[2]
    
    p1a = p1_hat[(x1>=1)&(x2>=1)]
    p2a = p2_hat[(x1>=1)&(x2>=1)]
    p1b = p1_hat[(x1==0)&(x2>=1)]
    p2b = p2_hat[(x1==0)&(x2>=1)]
    p2a = p2_hat[(x2>=1)]
    sub1 = sum((1-p1a)^mm1*(1-p2a)^n2*(1-(1-p2a)^(mm2s))/((1-(1-p1a)^n1)*(1-(1-p2a)^n2)))
    # sub1 = sum(h0hat(p1a,p2a,mm1,mm2s,n2)/((1-(1-p1a)^n1)*(1-(1-p2a)^n2)))
  })
}


Dq0_2_ext<-function(x1,x2,p1_hat,p2_hat,n1,n2,mmext){
  out = apply(mmext,1,function(mm) {
    mm1 = mm[1]
    mm2s = mm[2]
    
    p1a = p1_hat[(x1>=1)&(x2>=1)]
    p2a = p2_hat[(x1>=1)&(x2>=1)]
    p1b = p1_hat[(x1==0)&(x2>=1)]
    p2b = p2_hat[(x1==0)&(x2>=1)]
    #p2a = p2_hat[(x2>=1)]
    p1anew<<-p1a
    p2anew<<-p2a
    sub1=sum((1-p1a)^mm1*(1-p2a)^n2*(1-(1-p2a)^(mm2s))/((1-(1-p1a)^n1)*(1-(1-p2a)^n2)))
    
    # sub1 = sum((1-p1a)^mm1*(1-p2a)^n2*(1-(1-p2a)^(mm2s))/((1-(1-p1a)^n1)*(1-(1-p2a)^n2)))
    # sub1 = sum(h0hat(p1a,p2a,mm1,mm2s,n2)/((1-(1-p1a)^n1)*(1-(1-p2a)^n2)))
  })
}


Dq0_2_ext_new<-function(x1,x2,p1_hat,p2_hat,n1,n2,mmext){
  out = apply(mmext,1,function(mm) {
    mm1 = mm[1]
    mm2s = mm[2]
    
    p1a = p1_hat[(x1>=1)&(x2>=1)]
    p2a = p2_hat[(x1>=1)&(x2>=1)]
    p1b = p1_hat[(x1==0)&(x2>=1)]
    p2b = p2_hat[(x1==0)&(x2>=1)]
    p2a = p2_hat[(x2>=1)]
    
    p1anew<<-p1a
    p2anew<<-p2a
    sub1=sum((1-p2a)^n2*(1-(1-p2a)^(mm2s))/((1-(1-p2a)^n2)))
    
    # sub1 = sum((1-p1a)^mm1*(1-p2a)^n2*(1-(1-p2a)^(mm2s))/((1-(1-p1a)^n1)*(1-(1-p2a)^n2)))
    # sub1 = sum(h0hat(p1a,p2a,mm1,mm2s,n2)/((1-(1-p1a)^n1)*(1-(1-p2a)^n2)))
  })
}



Dq2<-function(x1,x2,mm,n1,n2){
  
  out = apply(mm,1,function(mm) {
    m1 = mm[1]
    m2 = mm[2]
    
    if(m1 !=0 || m2!=0){
      
      D2 = 1/(m1+m2)+sum(x1*(x1-1)/(n1*(n1-1)))*m1*(m1-1)/(m1+m2)^2+sum(x2*(x2-1)/(n2*(n2-1)))*m2*(m2-1)/(m1+m2)^2+2*sum(x1*x2/(n1*n2))*m1*m2/(m1+m2)^2
      D2<-1/D2
      
    }
    else if(m1==0 && m2==0)
      D2<-0
  })
}




Dq1_ext<-function(x1,x2,p1_hat,p2_hat,n1,n2,mm){
  
  out = apply(mm,1,function(mm) {
    m1 = mm[1]
    m2 = mm[2]
    # print(paste("m1=",m1,"m2=",m2,Sys.time()))
    h1hat<-h1_hat(p1_hat, p2_hat, x1,x2,m1, m2,n1,n2)
    # print(paste("h1hat=",h1hat,Sys.time()))
    dq1_ext<-exp(h1hat)
  })
}


Dq1_ext_choose<-function(x1,x2,p1_hat,p2_hat,n1,n2,mm){
  maxx1<-max(x1>0)
  maxx2<-max(x2>0)
  out = apply(mm,1,function(mm) {
    m1 = mm[1]
    m2 = mm[2]
    loopm1<-min(m1,maxx1)
    loopm2<-min(m2,maxx2)
    h1hat<-h1_hat_choose(p1_hat, p2_hat, x1,x2,m1, m2,n1,n2)
    print(paste("h1_hat_choose=",h1hat,Sys.time()))
    dq1_ext<-exp(h1hat)
  })
}






