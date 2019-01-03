#include <Rcpp.h>
using namespace Rcpp;


// You can include R code blocks in C++ files processed with sourceCpp
// (useful for testing and development). The R code will be automatically 
// run after the compilation.
//

// [[Rcpp::export]]
double Dm(int m, NumericVector Xi) {
  int D = sum(Xi>0);
  int n = sum(Xi);
  double out;
  NumericVector X = Xi[Xi>0];
  if(m <= n){
    int index1 = X.size();
    NumericVector fkhat(index1);
    for (int k = 0; k < index1; k++){
      fkhat[k] =  exp(Rf_lchoose(n - X[k], m) - Rf_lchoose(n, m)) ;
    }
    out = D - sum(fkhat);
  }else{
    
    int f1 = sum(X==1);
    
    if(f1 == 0){
      out = D;
    }else{
      int f2 = sum(X==2);
      double f0_hat;
      if(f2 > 0){
        f0_hat = (n-1)*pow(f1, 2)/(n*2*f2);
      }else{
        f0_hat = (n-1)*f1*(f1-1)/(2*n);
      }
      out = D + round(f0_hat)*(1-pow((1-f1/(n*f0_hat+f1)),(m-n)));
    }
    
  }
  return out ;
}


// [[Rcpp::export]]
double Dmm(int m1,int m2, NumericVector X1i, NumericVector X2i, NumericVector p1_hat, NumericVector p2_hat) {
  int n1 = sum(X1i);
  int n2 = sum(X2i);
  NumericVector X_share = X1i + X2i;
  NumericVector X_sharen = X_share[X_share>0];
  NumericVector X1i_share = X1i[X_share>0];
  NumericVector X2i_share = X2i[X_share>0];
  // NumericVector X1i_share = X1i[1-(X1i==0)-(X2i==0)+((X1i==0)&(X2i==0))];
  // NumericVector X2i_share = X2i[1-(X1i==0)-(X2i==0)+((X1i==0)&(X2i==0))];
  int index1 = X_sharen.size();
  // int index1 =sum(1-(X1i==0)-(X2i==0)+((X1i==0)&(X2i==0)));
  double output;
  NumericVector fkhat(index1);
  if((m1 <= n1) & (m2 <= n2)){
    for (int k = 0; k < index1; k++){
      fkhat[k] =  (1-exp(Rf_lchoose(n1 - X1i_share[k], m1) - Rf_lchoose(n1, m1)))*(1-exp(Rf_lchoose(n2 - X2i_share[k], m2) - Rf_lchoose(n2, m2))) ;
    }
    output = sum(fkhat);
  }else{
    p1_hat = p1_hat[(X_share>0)];
    p2_hat = p2_hat[(X_share>0)];
    // p1_hat = p1_hat[1-(X1i==0)-(X2i==0)+((X1i==0)&(X2i==0))];
    // p2_hat = p2_hat[1-(X1i==0)-(X2i==0)+((X1i==0)&(X2i==0))];
    if((m1 <= n1) & (m2 > n2)){
      
      for (int k = 0; k < index1; k++){
        fkhat[k] =  (1-exp(Rf_lchoose(n1 - X1i_share[k], m1) - Rf_lchoose(n1, m1)))*(1-exp(Rf_lchoose(n2 - X2i_share[k], n2) - Rf_lchoose(n2, n2))) ;
      }
      output = sum(fkhat) + sum(pow(1-p1_hat, m1)*pow(1-p2_hat, n2)*(1-pow((1-p2_hat), (m2-n2)))/(1-pow(1-p1_hat, n1)*pow(1-p2_hat, n2)));
    
    }else{
      for (int k = 0; k < index1; k++){
        fkhat[k] =  (1-exp(Rf_lchoose(n1 - X1i_share[k], n1) - Rf_lchoose(n1, n1)))*(1-exp(Rf_lchoose(n2 - X2i_share[k], m2) - Rf_lchoose(n2, m2))) ;
      }
      output = sum(fkhat) +  sum(pow(1-p1_hat, n1)*pow(1-p2_hat, m2)*(1-pow((1-p1_hat), (m1-n1)))/(1-pow(1-p1_hat, n1)*pow(1-p2_hat, n2)));
    
    }
  }
  return(output);
}

// [[Rcpp::export]]
double Dmm_share(int m1,int m2, NumericVector X1i, NumericVector X2i, NumericVector p1_hat, NumericVector p2_hat) {
  int n1 = sum(X1i);
  int n2 = sum(X2i);
  LogicalVector X_share = (X1i>=1) & (X2i>=1);
  NumericVector X1i_share = X1i[X_share];
  NumericVector X2i_share = X2i[X_share];
  // NumericVector X1i_share = X1i[1-(X1i==0)-(X2i==0)+((X1i==0)&(X2i==0))];
  // NumericVector X2i_share = X2i[1-(X1i==0)-(X2i==0)+((X1i==0)&(X2i==0))];
  int index1 = sum(X_share);
  // int index1 =sum(1-(X1i==0)-(X2i==0)+((X1i==0)&(X2i==0)));
  double output;
  NumericVector fkhat(index1);
  if((m1 <= n1) & (m2 <= n2)){
    for (int k = 0; k < index1; k++){
      fkhat[k] =  (1-exp(Rf_lchoose(n1 - X1i_share[k], m1) - Rf_lchoose(n1, m1)))*(1-exp(Rf_lchoose(n2 - X2i_share[k], m2) - Rf_lchoose(n2, m2))) ;
    }
    output = sum(fkhat);
  }else{
    p1_hat = p1_hat[(X_share>0)];
    p2_hat = p2_hat[(X_share>0)];
    
    if((m1 <= n1) & (m2 > n2)){
      
      for (int k = 0; k < index1; k++){
        fkhat[k] =  (1-exp(Rf_lchoose(n1 - X1i_share[k], m1) - Rf_lchoose(n1, m1)))*(1-exp(Rf_lchoose(n2 - X2i_share[k], n2) - Rf_lchoose(n2, n2))) ;
      }
      output = sum(fkhat) + sum(pow(1-p1_hat, m1)*pow(1-p2_hat, n2)*(1-pow((1-p2_hat), (m2-n2)))/(1-pow(1-p1_hat, n1)*pow(1-p2_hat, n2)));
      
    }else{
      for (int k = 0; k < index1; k++){
        fkhat[k] =  (1-exp(Rf_lchoose(n1 - X1i_share[k], n1) - Rf_lchoose(n1, n1)))*(1-exp(Rf_lchoose(n2 - X2i_share[k], m2) - Rf_lchoose(n2, m2))) ;
      }
      output = sum(fkhat) +  sum(pow(1-p1_hat, n1)*pow(1-p2_hat, m2)*(1-pow((1-p1_hat), (m1-n1)))/(1-pow(1-p1_hat, n1)*pow(1-p2_hat, n2)));
      
    }
  }
  return(output);
}

double Hypergeometric(int K, int k, int N,  int n) {
  return exp(Rf_lchoose(K,k)+Rf_lchoose(N-K,n-k)-Rf_lchoose(N,n));
  //return Rf_choose(K,k)*Rf_choose(N-K,n-k)/Rf_choose(N,n);
}

double fk(int k1, int k2, int m1, int m2, NumericVector x1, NumericVector y1){
  int n1 = sum(x1);
  int n2 = sum(y1);
  NumericVector x = x1[(x1>=k1) & (y1>=k2)];
  NumericVector y = y1[(x1>=k1) & (y1>=k2)];
  double output = 0;
  if(sum((x1>=k1) & (y1>=k2))==0){
    output = output;
  }else{
    for(int i = 0; i<x.size(); i++){
      output = output + Hypergeometric(x[i], k1, n1, m1)*Hypergeometric(y[i], k2, n2, m2);
    }
  }
  return output;
}


//for specially q
// [[Rcpp::export]]
double D_share(NumericVector xi,NumericVector yi,double m1, double m2,double q){
  //NumericVector xi = X(_,0);
  //NumericVector yi = X(_,1);
  double fk1k2 = 0;
  double output = 0;
  if(q==1){
    for(int k1 = 0; k1<(m1+1) ; k1++){
      for(int k2 = 0; k2<(m2+1) ; k2++){
        if((k1 == 0) & (k2 == 0)){
          fk1k2 = fk1k2; 
        }else{
          fk1k2 = fk1k2 -((k1+k2)/(m1+m2))*log((k1+k2)/(m1+m2))*fk(k1, k2, m1, m2, xi, yi);
        }
      }
      output = exp(fk1k2);
    }
  }else{
    for(int k1 = 0; k1<(m1+1) ; k1++){
      for(int k2 = 0; k2<(m2+1) ; k2++){
        if((k1 == 0) & (k2 == 0)){
          fk1k2 = fk1k2; 
        }else{
          fk1k2 = fk1k2 + pow(((k1+k2)/(m1+m2)), q)*fk(k1, k2, m1, m2, xi, yi);
        }
      }
    }
    output = pow(fk1k2, 1/(1-q));
  }
  return output;
}


double D_share_share(NumericVector xi,NumericVector yi,double m1, double m2,double q){
  //NumericVector xi = X(_,0);
  //NumericVector yi = X(_,1);
  double fk1k2 = 0;
  double output = 0;
  if(q==1){
    for(int k1 = 0; k1<(m1+1) ; k1++){
      for(int k2 = 0; k2<(m2+1) ; k2++){
        if(!((k1 >= 1) & (k2 >= 1))){
          fk1k2 = fk1k2; 
        }else{
          fk1k2 = fk1k2 -((k1+k2)/(m1+m2))*log((k1+k2)/(m1+m2))*fk(k1, k2, m1, m2, xi, yi);
        }
      }
      output = exp(fk1k2);
    }
  }else{
    for(int k1 = 0; k1<(m1+1) ; k1++){
      for(int k2 = 0; k2<(m2+1) ; k2++){
        if(!((k1 >= 1) & (k2 >= 1))){
          fk1k2 = fk1k2; 
        }else{
          fk1k2 = fk1k2 + pow(((k1+k2)/(m1+m2)), q)*fk(k1, k2, m1, m2, xi, yi);
        }
      }
    }
    output = pow(fk1k2, 1/(1-q));
  }
  return output;
}

// [[Rcpp::export]]
double Efk_q1(int m1, int m2, int k1, int k2,NumericVector p1_hat, NumericVector p2_hat,NumericVector xi, NumericVector yi){
  double  output; 
  int n1 = sum(xi);
  int n2 = sum(yi);
  
  NumericVector sub11;
  NumericVector sub12;
  NumericVector sub13;
  
  double sub1;
  double sub2;
  double sub3;
  
  sub11 = exp(Rf_lchoose(m1, k1)+k1*log(p1_hat)+(m1-k1)*log(1-p1_hat)+Rf_lchoose(m2, k2)+k2*log(p2_hat)+(m2-k2)*log(1-p2_hat)-log((1-pow(1-p1_hat, n1))*(1-pow(1-p2_hat, n2))));
  sub12 = exp(Rf_lchoose(m2, k2)+k2*log(p2_hat)+(m2-k2)*log(1-p2_hat)-log((1-pow(1-p2_hat, n2))));
  sub13 = exp(Rf_lchoose(m1, k1)+k1*log(p1_hat)+(m1-k1)*log(1-p1_hat)-log((1-pow(1-p1_hat, n1))));
  sub11 = sub11[((xi>=1)&(yi>=1))];
  if(k1==0){
    if(m1==0){
      sub12 = sub12[((yi>=1))];
    }else{sub12 = sub12[((xi==0)&(yi>=1))];}
    }else{sub12 = 0;}
  if(k2==0){
    if(m2==0){
      sub13 = sub13[((xi>=1))];
    }else{
      sub13 = sub13[((xi>=1)&(yi==0))];
    }
  }else{sub13 = 0;}
  sub1 = sum(sub11);
  sub2 = sum(sub12);
  sub3 = sum(sub13);
  if(m1==0){
    output = sub2;  
  }else if(m2==0){
    output = sub3;
  }else{
    sub1 = 0;
    output = sub1+sub2+sub3;
  }
  
  return output;
}
  
// [[Rcpp::export]]
double D1_f(NumericVector xi,NumericVector yi,double m1, double m2, NumericVector p1_hat, NumericVector p2_hat){
  int n1 = sum(xi);
  int n2 = sum(yi);
  double output;
  NumericVector X_share = xi + yi;
  double sub1 = 0;
  double sub2 = 0;
  
  if((m1 <= n1) & (m2 <= n2)){
    output = D_share(xi, yi, m1, m2, 1);
  }else{
    if((m1 <= n1) & (m2 > n2)){
      for(int k1 = 0; k1 < (m1+1); k1++){
        for(int k2 = 0; k2 < (m2+1); k2++){
          if((k1==0)&(k2==0)){
            sub1 = sub1;
          }else{
            sub1 = sub1 -(k1+k2)/(m1+m2)*log((k1+k2)/ (m1+m2))*Efk_q1(m1, m2, k1, k2, p1_hat, p2_hat, xi, yi);
          }
        }
      }
      for(int k1 = 0; k1 < (m1+1); k1++){
        for(int k2 = 0; k2 < (n2+1); k2++){
          if((k1==0)&(k2==0)){
            sub2 = sub2;
          }else{
            sub2 = sub2 -(k1+k2)/(m1+n2)*log((k1+k2)/ (m1+n2))*Efk_q1(m1, n2, k1, k2, p1_hat, p2_hat, xi, yi); 
          }
        }
      }
      output = D_share(xi, yi, m1, n2, 1)*exp(sub1-sub2);
    }else{
      for(int k2 = 0; k2 < (m2+1); k2++){
        for(int k1 = 0; k1 < (m1+1); k1++){
          if((k1==0)&(k2==0)){
            sub1 = sub1;
          }else{
            sub1 = sub1 -(k1+k2)/(m1+m2)*log((k1+k2)/(m1+m2))*Efk_q1(m1, m2, k1, k2, p1_hat, p2_hat, xi, yi);
              }
        }
      }
      for(int k2 = 0; k2 < (m2+1); k2++){
        for(int k1 = 0; k1 < (n1+1); k1++){
          if((k1==0)&(k2==0)){
            sub2 = sub2;
          }else{
            sub2 = sub2 -(k1+k2)/(n1+m2)*log((k1+k2)/(n1+m2))*Efk_q1(n1, m2, k1, k2, p1_hat, p2_hat, xi, yi);  
          }
        }
      }
      output = D_share(xi, yi, n1, m2, 1) * exp(sub1)/exp(sub2);
      //output = sub1;
    }
  }
    
  
  return output;
}

// [[Rcpp::export]]
NumericVector D1_f_share(NumericVector xi,NumericVector yi,double m1, double m2, NumericVector p1_hat, NumericVector p2_hat){
  int n1 = sum(xi);
  int n2 = sum(yi);
  double output_1;
  double output_2;
  NumericVector X_share = xi + yi;
  double sub1_1 = 0;
  double sub2_1 = 0;
  double sub1_2 = 0;
  double sub2_2 = 0;
  NumericVector sub3 = 0;
  if((m1 <= n1) & (m2 <= n2)){
    output_1 = D_share_share(xi, yi, m1, m2, 1);
    output_2 = D_share_share(xi, yi, m1, m2, 2);
  }else{
    if((m1 <= n1) & (m2 > n2)){
      for(int k1 = 0; k1 < (m1+1); k1++){
        for(int k2 = 0; k2 < (m2+1); k2++){
          if(!((k1>=1) & (k2>=1))){
            sub1_1 = sub1_1;
            sub1_2 = sub1_2;
          }else{
            sub3 = exp(Rf_lchoose(m1, k1)+k1*log(p1_hat)+(m1-k1)*log(1-p1_hat)+Rf_lchoose(m2, k2)+k2*log(p2_hat)+(m2-k2)*log(1-p2_hat)-log((1-pow(1-p1_hat, n1))*(1-pow(1-p2_hat, n2))));
            sub3 = sub3[((xi>=1)&(yi>=1))];
            sub1_1 = sub1_1 -(k1+k2)/(m1+m2)*log((k1+k2)/ (m1+m2))*sum(sub3);
            sub1_2 = sub1_2 + pow(((k1+k2)/(m1+m2)),2)*sum(sub3);
            
          }
        }
        
      }
      
      for(int k1 = 0; k1 < (m1+1); k1++){
        for(int k2 = 0; k2 < (n2+1); k2++){
          if(!((k1>=1) & (k2>=1))){
            sub2_1 = sub2_1;
            sub2_2 = sub2_2;
          }else{
            sub3 = exp(Rf_lchoose(m1, k1)+k1*log(p1_hat)+(m1-k1)*log(1-p1_hat)+Rf_lchoose(n2, k2)+k2*log(p2_hat)+(n2-k2)*log(1-p2_hat)-log((1-pow(1-p1_hat, n1))*(1-pow(1-p2_hat, n2))));
            sub3 = sub3[((xi>=1)&(yi>=1))];
            sub2_1 = sub2_1 -(k1+k2)/(m1+n2)*log((k1+k2)/ (m1+n2))*sum(sub3); 
            sub2_2 = sub2_2 + pow(((k1+k2)/(m1+n2)), 2)*sum(sub3);
          }
        }
      }
      
        output_1 = D_share_share(xi, yi, m1, n2, 1)*exp(sub1_1-sub2_1);
        output_2 = D_share_share(xi, yi, m1, n2, 2)*pow(sub1_2, 1/(1-2))/pow(sub2_2, 1/(1-2));
    }else{
      for(int k2 = 0; k2 < (m2+1); k2++){
        for(int k1 = 0; k1 < (m1+1); k1++){
          if(!((k1>=1) & (k2>=1))){
            sub1_1 = sub1_1;
            sub1_2 = sub1_2;
          }else{
            sub3 = exp(Rf_lchoose(m1, k1)+k1*log(p1_hat)+(m1-k1)*log(1-p1_hat)+Rf_lchoose(m2, k2)+k2*log(p2_hat)+(m2-k2)*log(1-p2_hat)-log((1-pow(1-p1_hat, n1))*(1-pow(1-p2_hat, n2))));
            sub3 = sub3[((xi>=1)&(yi>=1))];
            sub1_1 = sub1_1 -(k1+k2)/(m1+m2)*log((k1+k2)/ (m1+m2))*sum(sub3); 
            sub1_2 = sub1_2 + pow(((k1+k2)/(m1+m2)), 2)*sum(sub3) ;
            
          }
        }
      }
      for(int k2 = 0; k2 < (m2+1); k2++){
        for(int k1 = 0; k1 < (n1+1); k1++){
          if(!((k1>=1) & (k2>=1))){
            sub2_1 = sub2_1;
            sub2_2 = sub2_2;
          }else{
            sub3 = exp(Rf_lchoose(n1, k1)+k1*log(p1_hat)+(n1-k1)*log(1-p1_hat)+Rf_lchoose(m2, k2)+k2*log(p2_hat)+(m2-k2)*log(1-p2_hat)-log((1-pow(1-p1_hat, n1))*(1-pow(1-p2_hat, n2))));
            sub3 = sub3[((xi>=1)&(yi>=1))];
            sub2_1 = sub2_1 -(k1+k2)/(n1+m2)*log((k1+k2)/ (n1+m2))*sum(sub3); 
            sub2_2 = sub2_2 + pow(((k1+k2)/(n1+m2)), 2)*sum(sub3) ;
            
          }
        }
      }
        output_1 = D_share_share(xi, yi, n1, m2, 1) * exp(sub1_1)/exp(sub2_1);
        output_2 = D_share_share(xi, yi, n1, m2, 2) * pow(sub1_2, 1/(1-2))/pow(sub2_2, 1/(1-2));
      //output = sub1;
    }
  }
  
  
  return NumericVector::create(output_1,output_2);
}


// [[Rcpp::export]]
double D0_f(NumericVector xi,NumericVector yi,double m1, double m2, NumericVector p1_hat, NumericVector p2_hat){
  int n1 = sum(xi);
  int n2 = sum(yi);
  double output;
  // NumericVector X_share = xi + yi;
  double sub1 = 0;
  
  if((m1 <= n1) & (m2 <= n2)){
    output = D_share(xi, yi, m1, m2, 0);
  }else{
    if((m1 <= n1) & (m2 > n2)){
      for(int k1 = 0; k1 < (m1+1); k1++){
        for(int k2 = (n2+1); k2 < (m2+1); k2++){
          sub1 = sub1 + sum(exp(Rf_lchoose(m1, k1)+k1*log(p1_hat)+(m1-k1)*log(1-p1_hat)+Rf_lchoose(m2, k2)+k2*log(p2_hat)+(m2-k2)*log(1-p2_hat)-log(1-pow(1-p1_hat, n1)*pow(1-p2_hat, n2)))); 
        }
      }
      output = D_share(xi, yi, m1, n2, 0) + sub1;
    }else{
      for(int k2 = 0; k2 < (m2+1); k2++){
        for(int k1 = (n1+1); k1 < (m1+1); k1++){
          sub1 = sub1 + sum(exp(Rf_lchoose(m1, k1)+k1*log(p1_hat)+(m1-k1)*log(1-p1_hat)+Rf_lchoose(m2, k2)+k2*log(p2_hat)+(m2-k2)*log(1-p2_hat)-log(1-pow(1-p1_hat, n1)*pow(1-p2_hat, n2)))); 
        }
      }
      output = D_share(xi, yi, n1, m2, 0) + sub1;
    }
  }
  return output;
}

// [[Rcpp::export]]
NumericVector sh_un_abun(NumericVector xi,int n, int m){
  int s = xi.size();
  NumericVector out(s, 0.0);
  
  for(int i=0;i<s;i++){
    out[i] = 1-exp(Rf_lchoose((n-xi[i]),m)-Rf_lchoose(n,m));
  }
  return(out);
}

// [[Rcpp::export]]
NumericVector sh_sh_abun(NumericVector xi1, NumericVector xi2, int n1, int m1, int n2, int m2){
  int s = xi1.size();
  NumericVector out(s, 0.0);
  
  for(int i=0;i<s;i++){
    out[i] = (1-exp(Rf_lchoose((n1-xi1[i]),m1)-Rf_lchoose(n1,m1)))*(1-exp(Rf_lchoose((n2-xi2[i]),m2)-Rf_lchoose(n2,m2)));
  }
  return(out);
}