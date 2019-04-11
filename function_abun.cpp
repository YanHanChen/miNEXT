#include <Rcpp.h>
using namespace Rcpp;


// You can include R code blocks in C++ files processed with sourceCpp
// (useful for testing and development). The R code will be automatically 
// run after the compilation.
// this is created by hsiaotung at 201903


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

// this is D_share by HsiaoTung to improve speed
// [[Rcpp::export]]
double D_share(NumericVector xi,NumericVector yi,double m1, double m2,double q){
  //NumericVector xi = X(_,0);
  //NumericVector yi = X(_,1);
  double fk1k2 = 0;
  double output = 0;
  
  double max1 = max(xi);
  double max2 = max(yi);
  double m1loop = std::min(m1,max1);
  double m2loop = std::min(m2,max2);
  
  
  if(q==1){
    for(int k1 = 0; k1<(m1loop+1) ; k1++){
      for(int k2 = 0; k2<(m2loop+1) ; k2++){
        if((k1 == 0) & (k2 == 0)){
          fk1k2 = fk1k2; 
        }else{
          fk1k2 = fk1k2 -((k1+k2)/(m1+m2))*log((k1+k2)/(m1+m2))*fk(k1, k2, m1, m2, xi, yi);
        }
      }
      
      output = exp(fk1k2);
      
    }
    if ((m1loop==0) & (m2loop==0)) output=1;
  }else{
    for(int k1 = 0; k1<(m1loop+1) ; k1++){
      for(int k2 = 0; k2<(m2loop+1) ; k2++){
        if((k1 == 0) & (k2 == 0)){
          fk1k2 = fk1k2; 
        }else{
          fk1k2 = fk1k2 + pow(((k1+k2)/(m1+m2)), q)*fk(k1, k2, m1, m2, xi, yi);
        }
      }
    }
    output = pow(fk1k2, 1/(1-q));
    if ((m1loop==0) & (m2loop==0)) output=0;
  }
  return output;
}



double h0_cpp(double pi1, double pi2, int m1,  int m2s, int n2 ){
  double output = 0;
  output = pow((1-pi1),m1)*pow((1-pi2),n2)*(1-pow((1-pi2),m2s));
  return output;
}

// [[Rcpp::export]]
double h0_hat_cpp(NumericVector pi1, NumericVector pi2, int m1, int m2s, int n1, int n2){
  double output_all= 0; 
  //  double output_sh = 0;
  if(m1>0){
    NumericVector pi1_tmp = pi1[(pi1>0) & (pi2>0)];
    NumericVector pi2_tmp = pi2[(pi1>0) & (pi2>0)];
    double sumsh = 0;
    for(int i=0; i < pi1_tmp.size(); i++){
      // sumsh = sumsh +  h0(pi1_tmp[i],pi2_tmp[i],m1,m2s,n2)/(1-pow(1-pi1_tmp[i], n1))/(1-pow(1-pi2_tmp[i], n2));
      sumsh = sumsh +  h0_cpp(pi1_tmp[i],pi2_tmp[i],m1,m2s,n2)/((1-pow(1-pi1_tmp[i], n1))*(1-pow(1-pi2_tmp[i], n2)));
    }
    
    pi1_tmp = pi1[(pi1==0) & (pi2>0)];
    pi2_tmp = pi2[(pi1==0) & (pi2>0)];
    double sumx0 = 0;
    for(int i=0; i < pi1_tmp.size(); i++){
      sumx0 = sumx0 +  h0_cpp(0,pi2_tmp[i],m1,m2s,n2)/(1-pow(1-pi2_tmp[i], n2));
    }
    output_all = sumsh+sumx0;
    //    output_sh = sumsh;
  }
  else if(m1 == 0){
    NumericVector pi2_tmp = pi2[(pi2>0)];
    double sum2 = 0;
    for(int i=0; i < pi2_tmp.size(); i++){
      sum2 = sum2 + h0_cpp(0,pi2_tmp[i],0,m2s,n2)/(1-pow(1-pi2_tmp[i], n2));
    }
    output_all = sum2;
    //    output_sh = sum2;
  }
  //  NumericVector output = NumericVector::create(output_all, output_sh);
  double output = output_all;
  return output;
}





//20190326 new

double Efk_q1_new(double pi1, double pi2,int m1,int m2,int k1,int k2){
  double result;
  if((pi1 == 0) & (k1==0)){
    result = Rf_dbinom(k2, m2 , pi2, 0);
    
  }
  else if ((pi2 == 0) & (k2==0)){
    result = Rf_dbinom(k1, m1 , pi1, 0);
  }
  else{
    result = Rf_dbinom(k1, m1 , pi1, 0)*Rf_dbinom(k2, m2 , pi2, 0);
  }
  if((k1 == 0) & (k2==0)){
    result = 0;
  }
  
  
  return result;
}


//20190326 to test differenct subfunction for choose and dbinom

double Efk_q1_new_choose(double pi1, double pi2,int m1,int m2,int k1,int k2){
  double result;
  if((pi1 == 0) & (k1==0)){
    //result = Rf_dbinom(k2, m2 , pi2, 0);
    result=exp(Rf_lchoose(m2,k2)+k1*log(pi2)+(m2-k2)*(1-pi2));
    
  }
  else if ((pi2 == 0) & (k2==0)){
    //result = Rf_dbinom(k1, m1 , pi1, 0);
    result=exp(Rf_lchoose(m1,k1)+k1*log(pi1)+(m1-k1)*(1-pi1));
  }
  else{
    //result = Rf_dbinom(k1, m1 , pi1, 0)*Rf_dbinom(k2, m2 , pi2, 0);
    result=exp(Rf_lchoose(m1,k1)+k1*log(pi1)+(m1-k1)*(1-pi1)+Rf_lchoose(m2,k2)+k1*log(pi2)+(m2-k2)*(1-pi2));
  }
  if((k1 == 0) & (k2==0)){
    result = 0;
  }
  
  
  return result;
}





double h1(double pi1, double pi2, double m1, double m2, double n1, double n2,int s){
  
  double tmp1 = 0;
  double tmp2 = 0;
  
  for(int k2=s; k2 <= m2; k2++){
    for(int k1=s; k1 <= m1; k1++){
      if((k1 == 0) & (k2 == 0)){ tmp1 = 0; }
      else{ tmp1 = tmp1 + (k1+k2)/(m1+m2)*log((k1+k2)/(m1+m2))*Efk_q1_new(pi1,pi2,m1,m2,k1,k2); }
    }
  }
  tmp1 = -tmp1;
  
  for(int k2=s; k2 <= n2; k2++){
    for(int k1=s; k1 <= m1; k1++){
      if((k1 == 0) & (k2 == 0)){ tmp2 = 0; }
      else{ tmp2 = tmp2 + (k1+k2)/(m1+n2)*log((k1+k2)/(m1+n2))*Efk_q1_new(pi1,pi2,m1,n2,k1,k2); }
    }
  }
  tmp2 = -tmp2;
  // Rcout << "h1 function, tmp1= " << tmp1 << std::endl;
  // Rcout << "h1 function, tmp2= " << tmp2 << std::endl;
  double result = tmp1-tmp2;
  return result;
}



double h1_choose(double pi1, double pi2, double m1, double m2, double n1, double n2,int s){
  
  double tmp1 = 0;
  double tmp2 = 0;
  
  for(int k2=s; k2 <= m2; k2++){
    for(int k1=s; k1 <= m1; k1++){
      if((k1 == 0) & (k2 == 0)){ tmp1 = 0; }
      else{ tmp1 = tmp1 + (k1+k2)/(m1+m2)*log((k1+k2)/(m1+m2))*Efk_q1_new_choose(pi1,pi2,m1,m2,k1,k2); }
    }
  }
  tmp1 = -tmp1;
  
  for(int k2=s; k2 <= n2; k2++){
    for(int k1=s; k1 <= m1; k1++){
      if((k1 == 0) & (k2 == 0)){ tmp2 = 0; }
      else{ tmp2 = tmp2 + (k1+k2)/(m1+n2)*log((k1+k2)/(m1+n2))*Efk_q1_new_choose(pi1,pi2,m1,n2,k1,k2); }
    }
  }
  tmp2 = -tmp2;
  // Rcout << "h1 function, tmp1= " << tmp1 << std::endl;
  // Rcout << "h1 function, tmp2= " << tmp2 << std::endl;
  double result = tmp1-tmp2;
  return result;
}



double h1_assem2( double pi2,double m2, double n2){
  
  
  double tmp1 = 0;
  double tmp2 = 0;
  
  for(int k2=1; k2 <= m2; k2++){
    
    tmp1 = tmp1 + (k2)/(m2)*log((k2)/(m2))*Efk_q1_new(0,pi2,0,m2,0,k2);
  }
  tmp1 = -tmp1;
  
  for(int k2=1; k2 <= n2; k2++){
    
    tmp2 = tmp2 + (k2)/(n2)*log((k2)/(n2))*Efk_q1_new(0,pi2,0,n2,0,k2); 
    
  }
  tmp2 = -tmp2;
  double result = tmp1-tmp2;
  //Rcout << "tmp1  is " << tmp1 << std::endl;
  // Rcout << "tmp2  is " << tmp2 << std::endl;
  //Rcout << "h1_assem2  is " << result << std::endl;
  return result;
}





// [[Rcpp::export]]
double h1_hat(NumericVector pi1, NumericVector pi2, NumericVector xi1, NumericVector xi2, int m1, int m2, int n1, int n2){
  double output_all = 0;
  //  double output_sh = 0;
  if(m1>0){
    NumericVector pi1_tmp = pi1[(xi1>0) & (xi2>0)];
    NumericVector pi2_tmp = pi2[(xi1>0) & (xi2>0)];
    double n1 = sum(xi1);
    double n2 = sum(xi2);
    double sumsh = 0;
    //Rcout << "h1_hat size=" << pi1_tmp.size() << std::endl;
    //    double sumsh_p = 0;
    for(int i=0; i < pi1_tmp.size(); i++){
      sumsh = sumsh +  h1(pi1_tmp[i],pi2_tmp[i],m1,m2,n1,n2,0)/(1-pow(1-pi1_tmp[i], n1))/(1-pow(1-pi2_tmp[i], n2));
      //      sumsh_p = sumsh_p +  h1(pi1_tmp[i],pi2_tmp[i],t1,t2,UT1,UT2,T1,T2,1)/(1-pow(1-pi1_tmp[i], T1))/(1-pow(1-pi2_tmp[i], T2));
    }
    
    pi1_tmp = pi1[(xi1==0) & (xi2>0)];
    pi2_tmp = pi2[(xi1==0) & (xi2>0)];
    double sumx0 = 0;
    for(int i=0; i < pi1_tmp.size(); i++){
      sumx0 = sumx0 +  h1(0,pi2_tmp[i],m1,m2,n1,n2,0)/(1-pow(1-pi2_tmp[i], n2));
    }
    
    pi1_tmp = pi1[(xi1>0) & (xi2==0)];
    pi2_tmp = pi2[(xi1>0) & (xi2==0)];
    double sumy0 = 0;
    for(int i=0; i < pi1_tmp.size(); i++){
      sumy0 = sumy0 +  h1(pi1_tmp[i],0,m1,m2,n1,n2,0)/(1-pow(1-pi1_tmp[i], n1));
    }
    output_all = sumsh+sumx0+sumy0;
    //    output_sh = sumsh_p;
    
  }
  else if(m1 == 0){
    NumericVector pi2_tmp = pi2[(xi2>0)];
    //Rcout << "The value pi2_tmp is " << pi2_tmp << std::endl;
    double sum2 = 0;
    
    for(int i=0; i < pi2_tmp.size(); i++){
      
      sum2 = sum2 +  h1_assem2(pi2_tmp[i],m2,n2)/(1-pow(1-pi2_tmp[i], n2));
      //Rcout << "The value sum2 z is " << z << std::endl;
    }
    output_all = sum2;
    //    output_sh = sum2;
  }
  //  NumericVector output = NumericVector::create(output_all, output_sh);
  double output = output_all;
  return output;
}




// [[Rcpp::export]]
double h1_hat_choose(NumericVector pi1, NumericVector pi2, NumericVector xi1, NumericVector xi2, int m1, int m2, int n1, int n2){
  double output_all = 0;
  //  double output_sh = 0;
  if(m1>0){
    NumericVector pi1_tmp = pi1[(xi1>0) & (xi2>0)];
    NumericVector pi2_tmp = pi2[(xi1>0) & (xi2>0)];
    double n1 = sum(xi1);
    double n2 = sum(xi2);
    double sumsh = 0;
    Rcout << "pi1_tmp.size " << pi1_tmp.size() << std::endl;
    //    double sumsh_p = 0;
    for(int i=0; i < pi1_tmp.size(); i++){
      sumsh = sumsh +  h1_choose(pi1_tmp[i],pi2_tmp[i],m1,m2,n1,n2,0)/(1-pow(1-pi1_tmp[i], n1))/(1-pow(1-pi2_tmp[i], n2));
      //      sumsh_p = sumsh_p +  h1(pi1_tmp[i],pi2_tmp[i],t1,t2,UT1,UT2,T1,T2,1)/(1-pow(1-pi1_tmp[i], T1))/(1-pow(1-pi2_tmp[i], T2));
    }
    
    pi1_tmp = pi1[(xi1==0) & (xi2>0)];
    pi2_tmp = pi2[(xi1==0) & (xi2>0)];
    double sumx0 = 0;
    for(int i=0; i < pi1_tmp.size(); i++){
      sumx0 = sumx0 +  h1_choose(0,pi2_tmp[i],m1,m2,n1,n2,0)/(1-pow(1-pi2_tmp[i], n2));
    }
    
    pi1_tmp = pi1[(xi1>0) & (xi2==0)];
    pi2_tmp = pi2[(xi1>0) & (xi2==0)];
    double sumy0 = 0;
    for(int i=0; i < pi1_tmp.size(); i++){
      sumy0 = sumy0 +  h1_choose(pi1_tmp[i],0,m1,m2,n1,n2,0)/(1-pow(1-pi1_tmp[i], n1));
    }
    output_all = sumsh+sumx0+sumy0;
    //    output_sh = sumsh_p;
    
  }
  else if(m1 == 0){
    NumericVector pi2_tmp = pi2[(xi2>0)];
    //Rcout << "The value pi2_tmp is " << pi2_tmp << std::endl;
    double sum2 = 0;
    
    for(int i=0; i < pi2_tmp.size(); i++){
      
      sum2 = sum2 +  h1_assem2(pi2_tmp[i],m2,n2)/(1-pow(1-pi2_tmp[i], n2));
      //Rcout << "The value sum2 z is " << z << std::endl;
    }
    output_all = sum2;
    //    output_sh = sum2;
  }
  //  NumericVector output = NumericVector::create(output_all, output_sh);
  double output = output_all;
  return output;
}


// [[Rcpp::export]]
NumericVector un_abun(NumericVector xi,int n, int m){
  int s = xi.size();
  NumericVector out(s, 0.0);
  
  for(int i=0;i<s;i++){
    out[i] = 1-exp(Rf_lchoose((n-xi[i]),m)-Rf_lchoose(n,m));
  }
  return(out);
}

// [[Rcpp::export]]
NumericVector sh_abun(NumericVector xi1, NumericVector xi2, int n1, int m1, int n2, int m2){
  int s = xi1.size();
  NumericVector out(s, 0.0);
  
  for(int i=0;i<s;i++){
    out[i] = (1-exp(Rf_lchoose((n1-xi1[i]),m1)-Rf_lchoose(n1,m1)) * exp(Rf_lchoose((n2-xi2[i]),m2)-Rf_lchoose(n2,m2)) );
  }
  return(out);
}


