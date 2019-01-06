#*******************************************************************************
#********************************************************************************
#
## R scipts "Abundance" for Chao et al. (2019) paper on proportional mixture of rarefaction/extrapolation. 
## This R code for abundance data to compute and plot proportional mixture of two rarefaction/extrapolation curves, and
# the number of shared and unique species in any two mixture of rarefaction curves.
# This code includes three parts:
# (1) Comupute composite diversity of any sample and species composition (shared and unique species);
# (2) Plot individual rarefaction/extrapolation curves and mixture curves;
# (3) Example.
#
# NOTE: The packages "ggplot2", "dplyr", "reshape2", "ggpubr", "Rcpp", must be 
# installed and loaded before running the scripts. 
# 
#
#*******************************************************************************
#*******************************************************************************


####################################################################################
#
# (1). Comupute composite diversity of any sample and species composition (shared and unique species)
#
####################################################################################
library(Rcpp)
library(ggplot2)
library(dplyr)
library(ggpubr)
library(reshape2)
sourceCpp('function_abun.cpp')

#' abundance(data, knots = 20) for abundance data, comupute composite diversity of any sample and species composition (shared and unique species)
#' @param data a Sx2 dataframe, the intact assemblage (main) assemblage should be the first column.
#' @param allpts specifying whether to compute all combinations of sampling units of two assemblages. Default is FALSE.
#' @param size a vector specifying the smapling units of intact (main) assemblage. Default is NULL.
#' @param knots the number of points that the mixture diveristy will be computed. Default is 20.
#' @return a list containing 4 tables. Th first 3 are diversities of the two assemblages and the mixed one. The 4th table is the species composition of the mixed assemblage. 
abundance <- function(data1, knots = 10){
  x1 = data1[, 1]
  x2 = data1[, 2]
  n1 = sum(x1)
  n2 = sum(x2)
  D1 = sum(data1[, 1]>0)
  D2 = sum(data1[, 2]>0)
  n = c(n1, n2)
  D = c(D1, D2)
  tm_n = which.max(n)
  tm_D = which.max(D)
  if(tm_D==1){
    if(n1<=n2){
      m2 = round(seq(0, (n1), length.out = knots))
      m1 = n1 - m2
      m = cbind(m1, m2)
      m2 = round(seq(0, (n2), length.out = knots))
    }else{
      m2 = round(c(seq(0, n2, length.out = knots),
                   seq(n2, (n1), length.out = knots)))
      m1 = n1 - m2
      m = cbind(m1, m2)
    }
  }else{
    if(n1<n2){
      m1 = round(c(seq(0, n1, length.out = knots),
                   seq(n1, (n2), length.out = knots)))
      m2 = n2- m1 
      m = cbind(m1, m2)
    }else{
      m1 = round(seq(0, (n2), length.out = knots))
      m2 = m2 - m1 
      m = cbind(m1, m2)
      m1 = round(seq(0, (n2), length.out = knots))
    }
  }
  
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
  
  p1 = DetAbu(data1[,1], zero = T)
  p2 = DetAbu(data1[,2], zero = T)
  
  p = cbind(p1, p2)
  p1_hat = p[, 1]
  p2_hat = p[, 2]
  D0_1 =apply(cbind(m1, 0), 1, function(mm){
    m1 = mm[1]
    m2 = mm[2]
    if(m1<=n1&m2<=n2){
      out = D_share(x1, x2, m1, m2, 0)
    }else if(m1<=n1&m2>n2){
      p1a = p1_hat[(x1>=1)&(x2>=1)]
      p2a = p2_hat[(x1>=1)&(x2>=1)]
      p1b = p1_hat[(x1==0)&(x2>=1)]
      p2b = p2_hat[(x1==0)&(x2>=1)]
      sub1 = ifelse((m1==0)|(m2==0), sum((1-p1a)^m1*(1-p2a)^n2*(1-(1-p2a)^(m2-n2))/((1-(1-p1a)^n1)*(1-(1-p2a)^n2))), 0)
      sub2 = ifelse(m2==0, sum((1-p2b)^n2*(1-(1-p2b)^(m2-n2))/(1-(1-p2b)^n2)), 0)
      out = D_share(x1, x2, m1, n2, 0)+sub1+sub2
    }else{
      p1a = p1_hat[(x1>=1)&(x2>=1)]
      p2a = p2_hat[(x1>=1)&(x2>=1)]
      p1b = p1_hat[(x2==0)&(x1>=1)]
      p2b = p2_hat[(x2==0)&(x1>=1)]
      sub1 = ifelse((m1==0)|(m2==0), sum((1-p1a)^n1*(1-p2a)^m2*(1-(1-p1a)^(m1-n1))/((1-(1-p1a)^n1)*(1-(1-p2a)^n2))), 0)
      sub2 = ifelse(m1==0, sum((1-p1b)^n1*(1-(1-p1b)^(m1-n1))/(1-(1-p1b)^n1)), 0)
      out = D_share(x1, x2, n1, m2, 0)+sub1+sub2
    }
  })
  D0_2 =apply(cbind(0, m2), 1, function(mm){
    m1 = mm[1]
    m2 = mm[2]
    if(m1<=n1&m2<=n2){
      out = D_share(x1, x2, m1, m2, 0)
    }else if(m1<=n1&m2>n2){
      p1a = p1_hat[(x1>=1)&(x2>=1)]
      p2a = p2_hat[(x1>=1)&(x2>=1)]
      p1b = p1_hat[(x1==0)&(x2>=1)]
      p2b = p2_hat[(x1==0)&(x2>=1)]
      sub1 = ifelse((m1==0)|(m2==0), sum((1-p1a)^m1*(1-p2a)^n2*(1-(1-p2a)^(m2-n2))/((1-(1-p1a)^n1)*(1-(1-p2a)^n2))), 0)
      sub2 = ifelse(m2==0, sum((1-p2b)^n2*(1-(1-p2b)^(m2-n2))/(1-(1-p2b)^n2)), 0)
      out = D_share(x1, x2, m1, n2, 0)+sub1+sub2
    }else{
      p1a = p1_hat[(x1>=1)&(x2>=1)]
      p2a = p2_hat[(x1>=1)&(x2>=1)]
      p1b = p1_hat[(x2==0)&(x1>=1)]
      p2b = p2_hat[(x2==0)&(x1>=1)]
      sub1 = ifelse((m1==0)|(m2==0), sum((1-p1a)^n1*(1-p2a)^m2*(1-(1-p1a)^(m1-n1))/((1-(1-p1a)^n1)*(1-(1-p2a)^n2))), 0)
      sub2 = ifelse(m1==0, sum((1-p1b)^n1*(1-(1-p1b)^(m1-n1))/(1-(1-p1b)^n1)), 0)
      out = D_share(x1, x2, n1, m2, 0)+sub1+sub2
    }
  })
  D0 =apply(m, 1, function(mm){
    m1 = mm[1]
    m2 = mm[2]
    if(m1<=n1&m2<=n2){
      out = D_share(x1, x2, m1, m2, 0)
    }else if(m1<=n1&m2>n2){
      p1a = p1_hat[(x1>=1)&(x2>=1)]
      p2a = p2_hat[(x1>=1)&(x2>=1)]
      p1b = p1_hat[(x1==0)&(x2>=1)]
      p2b = p2_hat[(x1==0)&(x2>=1)]
      sub1 = ifelse((m1==0)|(m2==0), sum((1-p1a)^m1*(1-p2a)^n2*(1-(1-p2a)^(m2-n2))/((1-(1-p1a)^n1)*(1-(1-p2a)^n2))), 0)
      sub2 = ifelse(m2==0, sum((1-p2b)^n2*(1-(1-p2b)^(m2-n2))/(1-(1-p2b)^n2)), 0)
      out = D_share(x1, x2, m1, n2, 0)+sub1+sub2
    }else{
      p1a = p1_hat[(x1>=1)&(x2>=1)]
      p2a = p2_hat[(x1>=1)&(x2>=1)]
      p1b = p1_hat[(x2==0)&(x1>=1)]
      p2b = p2_hat[(x2==0)&(x1>=1)]
      sub1 = ifelse((m1==0)|(m2==0), sum((1-p1a)^n1*(1-p2a)^m2*(1-(1-p1a)^(m1-n1))/((1-(1-p1a)^n1)*(1-(1-p2a)^n2))), 0)
      sub2 = ifelse(m1==0, sum((1-p1b)^n1*(1-(1-p1b)^(m1-n1))/(1-(1-p1b)^n1)), 0)
      out = D_share(x1, x2, n1, m2, 0)+sub1+sub2
    }
  })

  m1_i = unique(m1[m1 >= n1-n2])
  m2_i = unique(n1 - m1_i)
  datash = data1[(data1[,1]>0 & data1[,2]>0), , drop=F]
  dataun1 = data1[(data1[,1]>0 & data1[,2]==0), , drop=F]
  dataun2 = data1[(data1[,1]==0 & data1[,2]>0), , drop=F]
  un1 = sapply(m1_i,function(i){
    sum(sh_un_abun(xi = dataun1[,1], n = n1, m = i))
  })
  un2 = sapply(m2_i,function(i){
    sum(sh_un_abun(xi = dataun2[,2], n = n2, m = i))
  })
  sh1 = sapply(m1_i,function(i){
    sum(sh_un_abun(xi = datash[,1], n = n1, m = i))
  })
  sh2 = sapply(m2_i,function(i){
    sum(sh_un_abun(xi = datash[,2], n = n2, m = i))
  })
  sh12 = sapply(1:length(m1_i),function(i){
    sum(sh_sh_abun(xi1 = datash[,1],xi2 = datash[,2], n1 = n1,m1 =  m1_i[i],
                   n2 = n2, m2 = m2_i[i]))
  })
  sh = sh1 + sh2 -sh12
  total = un1 + un2 + sh
  q0_ana = data.frame(m1 = m1_i, m2 = m2_i,
                      q0_un1 = un1, q0_un2 = un2,
                      q0_sh = sh)
  colnames(q0_ana)[3:5] = c(paste("Unique to",colnames(data1)[1]),
                            paste("Unique to",colnames(data1)[2]),
                            "Share")
  q0_ana = as.matrix(q0_ana)

  D1 = apply(m, 1, function(mi)D1_f(xi = x1, yi = x2, m1 = mi[1], m2 = mi[2], p1_hat = p1_hat, p2_hat = p2_hat))
  D1_1 = apply(cbind(m1, 0), 1, function(mi)D1_f(xi = x1, yi = x2, m1 = mi[1], m2 = 0, p1_hat = p1_hat, p2_hat = p2_hat))
  D1_2 = apply(cbind(0, m2), 1, function(mi)D1_f(xi = x1, yi = x2, m1 = 0, m2 = mi[2], p1_hat = p1_hat, p2_hat = p2_hat))
  x11 =x1
  x22 =x2
  x1 = x1[rowSums(data1)>0]
  x2 = x2[rowSums(data1)>0]
  D2.f  <- function(m){
    m1 = m[1]
    m2 = m[2]
    D2 = 1/(m1+m2)+sum(x1*(x1-1)/(n1*(n1-1)))*m1*(m1-1)/(m1+m2)^2+sum(x2*(x2-1)/(n2*(n2-1)))*m2*(m2-1)/(m1+m2)^2+2*sum(x1*x2/(n1*n2))*m1*m2/(m1+m2)^2
    ifelse(m1==0&m2==0, 0, 1/D2)
  }
  D2 = apply(m, 1, D2.f)
  D2_1 = apply(cbind(m1, 0), 1, D2.f)
  D2_2 = apply(cbind(0, m2), 1, D2.f)
  
  output0 = cbind(m1, m2, D0_1, D0_2, D0)
  output1 = cbind(m1, m2, D1_1, D1_2, D1)
  output2 = cbind(m1, m2, D2_1, D2_2, D2)
  
  colnames(output0) = c("m1", "m2", paste(colnames(data1)), "Mixed")
  colnames(output1) = c("m1", "m2", paste(colnames(data1)), "Mixed")
  colnames(output2) = c("m1", "m2", paste(colnames(data1)), "Mixed")
  
  list(q0 = output0, q1 = output1, q2 = output2, q0_ana = q0_ana)

}

####################################################################################
#
# (2). Plot individual rarefaction/extrapolation curves and mixture curves
#
####################################################################################
               
#' multi.plot(data, ans1, type) plot the outcome of abundance(abundance data) or incidence(incidence data).
#' @param data the Sx2 data used in incidence.
#' @param ans1 the outcome of incidence.
#' @param type datatype of data, "abundance" or "incidence". Default is "abundance".
#' @return a list containing two plots: $div for diversity of q = 0, 1, 2 and $comp for species composition.
multi.plot <- function(data, ans1, type = "abundance"){
  if(type == "abundance"){text1 = "individuals"}else{text1 = "sampling units"}
  if ( sum(data[,1] > 0) < sum(data[,2] > 0) ){data1 = data[,c(2,1)]}else{data1 = data}
  draw.f_Div <- function(data1, output){
    if(type == "abundance"){
      x1 = data1[, 1]
      x2 = data1[, 2]
      n1 =  sum(x1)
      n2 = sum(x2)
      n = c(n1, n2)
      D1 = sum(data1[, 1]>0)
      D2 = sum(data1[, 2]>0)
      D = c(D1, D2) 
    }else{
      x1 = data1[, 1]
      x2 = data1[, 2]
      n1 =  x1[1]
      n2 =  x2[1]
      n = c(n1, n2)
      D1 = sum(data1[, 1]>0)-1
      D2 = sum(data1[, 2]>0)-1
      D = c(D1, D2) 
    }
    
    lty1 = cbind(as.character(output[,1]<=n1),as.character(output[,2]<=n2))
    index1 = 1
    output = as.data.frame(output)
    output$m_main = output[, index1]
    output1 = data.frame(m = c(output$m1, output$m2, output$m_main), value = c(output[, 3], output[, 4], output[, 5]))
    nn = length(output[, 3])
    name1 = names(output)[c(3, 4)]
    output1$type1 = factor(c(rep(name1[1], nn), rep(name1[2], nn), rep("Mixed", nn)), levels = c(name1[1],name1[2],"Mixed"), labels = c(name1[1],name1[2],"Mixed"))
    
    if(n1 == n2){
      output1$type2 = factor(c(lty1[, 1], lty1[, 2], lty1[, which.min(n)]), levels = c("TRUE"), labels = c("Rarefaction"))
      output2 = data.frame(nn = c((n1), (n2)), va = c(output[, 3][output[, 1]==(n1)], output[, 4][output[, 2]==(n2)]))
    }else{
      if(which.max(D)==which.min(n)){
        output1$type2 = factor(c(lty1[, 1], lty1[, 2], lty1[, which.min(n)]), levels = c("TRUE"), labels = c("Rarefaction"))
        if(which.min(n)==1){
          output2 = data.frame(nn = c(n1, n2, n1), va = c(output[, 3][output[, 1]==(n1)][1], output[, 4][output[, 2]==(n2)], min(output[,5])))
        }else{
          output2 = data.frame(nn = c(n1, n2, n2), va = c(output[, 3][output[, 1]==(n1)], output[, 4][output[, 2]==(n2)][1], output[, 3][output[, 1]==(n2)]))
        }    
      }else{
        lty1[,which.min(n)][output[,which.min(n)]==n[which.min(n)]][2] = FALSE
        output1$type2 = factor(c(lty1[, 1], lty1[, 2], lty1[, which.min(n)]), levels = c("TRUE", "FALSE"), labels = c("Rarefaction", "Extrapolation"))  
        output2 = data.frame(nn = c(n1, (n2)), va = c(output[, 3][output[, 1]==(n1)][1], output[, 4][output[, 2]==(n2)][1]))
      }
    }
    output1 = output1[!is.na(output1$value),]
    output2 = output2[!is.na(output2$va),]
    
    pp = ggplot(output1)+
      geom_hline(aes(yintercept = output[,4][output[, 2]==max(output$m2)]), col = "darkgray", linetype = 3, size = 1.25)+
      geom_hline(aes(yintercept = output[,3][output[, 1]==max(output$m_main)][1]), col = "darkgray", linetype = 3, size = 1.25)+
      geom_line(aes(x = m, y = value, col = type1, size = type1,linetype = type2))+
      scale_size_manual(breaks=c(name1[1],name1[2],"Mixed"), values=c(1, 1, 3),
                        labels = c(name1[1],name1[2],"Mixed"))+
      scale_color_manual(values = c("black", "#3498DB","#C0392B"),
                         breaks = c(name1[1],name1[2],"Mixed"),
                         labels = c(name1[1],name1[2],"Mixed"))+
      guides(linetype = FALSE)+
      xlim(c(0, max(output1$m)+5))+
      ylab("Diversity")+
      theme_bw()+
      theme(panel.grid.major=element_blank(),panel.grid.minor=element_blank())+
      theme(legend.position="bottom",legend.title=element_blank(),
            legend.text=element_text(size=15),legend.key.width  = unit(1.5,"cm"))+
      theme(plot.title = element_text(size=14, face="bold.italic",hjust = 0))+
      theme(axis.text.x= element_text(size = 18,colour = "black",margin=unit(c(0.2,0.2,0.5,0.5), "cm")),
            axis.text.y= element_text(size = 18,colour = "black",margin=unit(c(0.2,0.2,0.2,0.2), "cm")),
            axis.title.x=element_text(size = 20,lineheight = 1.1),axis.title.y=element_text(size = 16),
            axis.ticks = element_line(size = 1.5,colour = "black"),axis.ticks.length = unit(-0.25,"line"))+
      scale_x_continuous(name = paste0("Number of ", text1, "\n(Proportion % in intact)"),
                         breaks = round(seq(0,max(output$m_main),length.out = 6)),
                         labels = paste0(round(seq(0,max(output$m_main),length.out = 6)), "\n(", round(seq(0,max(output$m_main),length.out = 6)/(max(output$m_main))*100, 1), "%)"))
    if(n1==n2){
      pp + geom_point(aes(output2$nn[1], output2$va[1]),colour="black",size=8)+
        geom_point(aes(output2$nn[2], output2$va[2]),colour="#3498DB",size=8)
    }else{
      if(which.max(D)==which.min(n)){
        if(which.min(n)==1){
          pp + geom_point(aes(output2$nn[1], output2$va[1]),colour="black",size=8) +
            geom_point(aes(output2$nn[2], output2$va[2]),colour="#3498DB",size=8) +
            geom_point(aes(output2$nn[3], output2$va[3]),colour="#3498DB",size=8, shape= 1)
          
        }else{
          pp + geom_point(aes(output2$nn[1], output2$va[1]),colour="black",size=8) +
            geom_point(aes(output2$nn[2], output2$va[2]),colour="#3498DB",size=8) +
            geom_point(aes(output2$nn[3], output2$va[3]),colour="#black",size=8, shape= 1)
        }
      }else{
        if(which.min(D)==1){color1 = 1}else{color1 ="#3498DB"}
        pp + geom_point(aes(output2$nn[1], output2$va[1]),colour="black",size=8)+
          geom_point(aes(output2$nn[2], output2$va[2]),colour="#3498DB",size=8)+
          geom_point(aes(max(n), max(output[, which.min(n)+2], na.rm = T)),colour=color1,size=8, shape=1)
      }
    }
  }
  plot_comb1 <- function(data1, output_d, output_p){
    if(type == "abundance"){
      x1 = data1[, 1]
      x2 = data1[, 2]
      n1 =  sum(x1)
      n2 = sum(x2)
      n = c(n1, n2)
      D1 = sum(data1[, 1]>0)
      D2 = sum(data1[, 2]>0)
      D = c(D1, D2) 
    }else{
      x1 = data1[, 1]
      x2 = data1[, 2]
      n1 =  x1[1]
      n2 =  x2[1]
      n = c(n1, n2)
      D1 = sum(data1[, 1]>0)-1
      D2 = sum(data1[, 2]>0)-1
      D = c(D1, D2) 
    }
    sites = colnames(data1)
    
    p_names = colnames(output_p)[3:5]
    output_p = output_p[,c(-2)]
    output_p = data.frame(output_p)
    output_p = melt(data = output_p,id.vars = c("m1"),variable.name = "col")
    output_p$col = as.character(output_p$col)
    
    output_d = data.frame(output_d)
    output_d = output_d[,-c(2,4)]
    output_d <- output_d %>% data.frame() %>% mutate(lty={
      (output_d$m1 >= n1-n2) %>% as.numeric()
    })
    if(sum(output_d$m1==(n1-n2))==2 ){
      output_d$lty[which(output_d$lty==1) %>% max]=0
    }
    output_d$lty[output_d$lty==1]=c("Rarefaction")
    output_d$lty[output_d$lty==0]=c("Extrapolation")
    
    output_d = melt(data = output_d,id.vars = c("m1","lty"),variable.name = "col")
    output_d[output_d$col== sites[1],2] = "Rarefaction"
    
    output_d$lty <- factor(x = output_d$lty,levels = c("Rarefaction","Extrapolation"))
    output_d$col <- factor(x = output_d$col,levels = c(sites[1],"Mixed"))
    output_p$col <- factor(x = output_p$col,levels = unique(output_p$col))
    output_p$lty <- "portion"
    output = rbind(output_d,output_p)
    pp = ggplot(output)+
      theme_bw()+
      geom_hline(aes(yintercept = output[output$col=="Mixed",4][1]), col = "darkgray", linetype = 3, size = 1.25)+
      geom_line(data = output, aes(x = m1, y = value, linetype = lty, col = col, size = col))+
      scale_size_manual( breaks = c(sites[1],"Mixed",p_names[3],p_names[1],p_names[2]),
                         values=c(1.5, 3, 2, 2, 2))+
      scale_color_manual(values=c("black", "#C0392B","black","#3498DB","purple"),
                         breaks = c(sites[1],"Mixed",p_names[3],p_names[1],p_names[2]))+
      scale_linetype_manual(values=c("solid", "2121","1111"))+
      xlim(c(0, max(output_d$m1)+5))+
      ylab("Diversity")+
      theme(panel.grid.major=element_blank(),panel.grid.minor=element_blank())+
      theme(legend.position="none",legend.title = element_blank(),
            legend.text=element_text(size=15),legend.key.width  = unit(1.5,"cm"))+
      theme(axis.text.x= element_text(size = 18,colour = "black",margin=unit(c(0.2,0.2,0.2,0.2), "cm")),
            axis.text.y= element_text(size = 18,colour = "black",margin=unit(c(0.2,0.2,0.2,0.2), "cm")),
            axis.title.x=element_text(size = 20),axis.title.y=element_text(size = 20),
            axis.ticks = element_line(size = 1.5),axis.ticks.length = unit(-0.25,"line"))+
      scale_x_continuous(name = paste0("Number of ", text1, "\n(Proportion % in ", sites[1], ")" ),
                         breaks = round(seq(0,max(output_d$m1),length.out = 6)),
                         labels = paste0(round(seq(0,max(output_d$m1),length.out = 6)), "\n(", round(seq(0,max(output_d$m1),length.out = 6)/(max(output_d$m1))*100, 1), "%)"))
  }
  
  p0 = draw.f_Div(data1,ans1[[1]])+ggtitle("q=0, species richness")+theme(plot.title = element_text(size = 20, face = "bold"))
  p1 = draw.f_Div(data1,ans1[[2]])+ggtitle("q=1")+theme(plot.title = element_text(size = 20, face = "bold"))
  p2 = draw.f_Div(data1,ans1[[3]])+ggtitle("q=2")+theme(plot.title = element_text(size = 20, face = "bold"))
  p0_ana = plot_comb1(data1 = data1,output_d = ans1[[1]],output_p = ans1[[4]])+ylim(c((min(unlist(ans1[[1]][, -(1:2)]))), ceiling(max(unlist(ans1[[1]][, -(1:2)])))))+
    ggtitle("q=0, species composition")+theme(plot.title = element_text(size = 20, face = "bold"))
  g1 = ggarrange(p0, p1, p2, ncol=1, nrow=3, common.legend = TRUE, legend="bottom",labels = c("(a)", "(b)", "(c)"),font.label = list(size = 20))
  return(list(div = g1, comp = p0_ana))
}
####################################################################################
#
# (2). Example
#
####################################################################################
spider = read.table("spider.txt")
result_abun_spider = abundance(spider, knots = 10)
multi.plot(data = spider, ans1 = result_abun_spider, type = "abundance")
