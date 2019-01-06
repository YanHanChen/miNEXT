#*******************************************************************************
#********************************************************************************
#
## R scipts "Incidence" for Chao et al. (2019) paper on proportional mixture of rarefaction/extrapolation. 
## This R code for incidence data to compute and plot proportional mixture of two rarefaction/extrapolation curves, and
# the number of shared and unique species in any two mixture of rarefaction curves.
# This code includes three parts:
# (1) Comupute composite diversity of any sample and species composition (shared and unique species);
# (2) Plot individual rarefaction/extrapolation curves and mixture curves;
# (3) Example.
#
# NOTE: The packages "ggplot2", "dplyr", "reshape2", "ggpubr", "Rcpp", "vegan" must be 
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
library(vegan)
library(ggpubr)
library(reshape2)
sourceCpp('function_inci.cpp')

#' dq_hat(data, allpts = FALSE, size = NULL, knots = 20 ) for incidence data, estimates the diversity and species composition of mixed assemblage consisting of various portions of two assemblages.
#' @param data a Sx2 dataframe, the intact assemblage (main) assemblage should be the first column.
#' @param allpts specifying whether to compute all combinations of sampling units of two assemblages. Default is FALSE.
#' @param size a vector specifying the smapling units of intact (main) assemblage. Default is NULL.
#' @param knots the number of points that the mixture diveristy will be computed. Default is 20.
#' @return a list containing 4 tables. Th first 3 are diversities of the two assemblages and the mixed one. The 4th table is the species composition of the mixed assemblage. 
dq_hat <- function(data, allpts = FALSE, size = NULL, knots = 20 ){
  d2_hat  <- function(x1,x2,t1,t2){
    if(t1!=0 || t2!=0){
      u1 = sum(x1)
      u2 = sum(x2)
      u1_hat = t1*u1/T1;
      u2_hat = t2*u2/T2;
      D2 = 1/(u1_hat+u2_hat)+sum(x1*(x1-1)/(T1*(T1-1)))*t1*(t1-1)/(u1_hat+u2_hat)^2+sum(x2*(x2-1)/(T2*(T2-1)))*t2*(t2-1)/(u1_hat+u2_hat)^2+2*sum(x1*x2)*t1*t2/(T1*T2)/(u1_hat+u2_hat)^2
      1/D2
    }
    else if (t1==0 && t2==0){
      D2 = 0
      D2
    }
    
  }
  DetInc <- function(y, zero=FALSE){
    y <- unlist(y)
    nT <- max(y)  
    y <- y[-1]
    Q1 <- sum(y==1)
    Q2 <- sum(y==2)
    Q3 <- sum(y==3)
    if(Q2==0){
      Q1 <- max(Q1 - 1, 0)
      Q2 <- 1
    }
    A1 <- Q1 / nT * ((nT-1)*Q1 / ((nT-1)*Q1 + 2*max(Q2,1)))
    A2 <- Q2 / choose(nT, 2) * ((nT-2)*Q2/((nT-2)*Q2 + 3*max(Q3,1)))^2
    if(zero==FALSE) y <- y[y>0]
    q.solve <- function(q){
      e <- A1 / sum(y/T*exp(-q*y))
      out <- sum((y/nT * (1 - e * exp(-q*y)))^2) - sum(choose(y,2)/choose(nT,2)) + A2
      abs(out)
    }
    q <- tryCatch(optimize(q.solve, c(0,1))$min, error = function(e) {1})
    e <- A1 / sum(y/nT*exp(-q*y))
    o <- y/nT * (1 - e * exp(-q*y))
    o
  }
  ts = as.numeric(data[1,])
  T1 = ts[1]
  T2 = ts[2]
  d0ob = c(sum(data[,1]>0),sum(data[,2]>0))
  if(T1 == T2){
    extr = 0
    if(d0ob[1]<d0ob[2]){
      data = data[,c(2,1)]
      # tmp = data[,1]
      # data[,1] = data[,2]
      # data[,2] = tmp
      # colnames(data) = rev(colnames(data))
    }
  }else if(((T1-T2)*(d0ob[1]-d0ob[2])) > 0){
    extr = 1
    if(T1<T2){
      data = data[,c(2,1)]
      # tmp = data[,1]
      # data[,1] = data[,2]
      # data[,2] = tmp
      # colnames(data) = rev(colnames(data))
    }
  }else if (((T1-T2) * (d0ob[1]-d0ob[2])) < 0 ){
    extr = 0
    if(T1>T2){
      data = data[,c(2,1)]
      # tmp = data[,1]
      # data[,1] = data[,2]
      # data[,2] = tmp
      # colnames(data) = rev(colnames(data))
    }
  }
  ts = as.numeric(data[1,])
  T1 = ts[1]
  T2 = ts[2]
  if(allpts == TRUE){
    t1 = seq(T1,0)
    t2 = T1-t1
  }else if (allpts == FALSE & !is.null(size)){
    if(max(size)>T1)
      stop(paste0("Because the community ",colnames(data)[1]," has higher divrsity, its sampling units should be viewed main."))
    
    t1 = sort(c(size,(T1-T2),T1),decreasing = T)
    t1 = t1[!duplicated(t1[t1>=0])]
    t2 = T1-t1
  }else if (allpts == FALSE & is.null(size)){
    t1 = round(seq(T1,0,length.out = knots),0)
    t1 = sort(c(t1,(T1-T2)),decreasing = T)
    t1 = t1[t1>=0] %>% unique
    t2 = T1-t1
  }
  
  datap = data[-1,]
  datap[,1] = DetInc(data[,1],zero = T)
  datap[,2] = DetInc(data[,2],zero = T)
  
  data = data[-1,]
  D = c(sum(data[,1]>0),sum(data[,2]>0))

  t1_in = sort(t1[t1>=(T1-T2)],decreasing = T) 
  t2_in = t1_in[1] - t1_in
  d01_int = t(sapply(1:length(t1_in),function(i) D0_rare(xi = c(T1,data[,1]),yi = c(T2,data[,2]),t1 = t1_in[i],t2 = t2_in[i])))
  
  d01 = d01_int

  sites = colnames(data)
  datash = data[(data[,1]>0 & data[,2]>0),]
  dataun1 = data[(data[,1]>0 & data[,2]==0),]
  dataun2 = data[(data[,1]==0 & data[,2]>0),]
  
  un1 = sapply(t1_in,function(i){
    sum(sh_un_inci(yi = dataun1[,1], T = T1, t = i))
  })
  un2 = sapply(t2_in,function(i){
    sum(sh_un_inci(yi = dataun2[,2], T = T2, t = i))
  })
  sh1 = sapply(t1_in,function(i){
    sum(sh_un_inci(yi = datash[,1], T = T1, t = i))
  })
  sh2 = sapply(t2_in,function(i){
    sum(sh_un_inci(yi = datash[,2], T = T2, t = i))
  })
  sh12 = sapply(1:length(t1_in),function(i){
    sum(sh_sh_inci(yi1 = datash[,1],yi2 = datash[,2], T1 = T1, t1 =  t1_in[i],
                   T2 = T2, t2 = t2_in[i]))
  })
  sh = sh1 + sh2 -sh12
  q0_ana = data.frame(m1 = t1_in, m2 = t2_in,
                      q0_un1 = un1, q0_un2 = un2,
                      q0_sh = sh)
  colnames(q0_ana)[3:5] = c(paste("Unique to",sites[1]),
                            paste("Unique to",sites[2]),
                            "Share")
  
  a1_q01 =  t(sapply(1:length(t1[t1<=T1]),function(i) D0_rare(xi = c(T1,data[,1]),yi = c(T2,data[,2]),t1 = t1[t1<=T1][i],t2 = 0))) 
  a1_q2 =  sapply(1:length(t1),function(i) d2_hat(x1 = data[,1],x2 = data[,2],sort(t1,decreasing = T)[i],0))
  
  
  if(extr == 1){ 
    t1_ex = sort(t1[t1<(T1-T2)],decreasing = T) 
    t2_ex = T1-t1_ex
    
    d01_ext = t(sapply(1:length(t1_ex),function(i) D0_rare(xi = c(T1,data[,1]),yi = c(T2,data[,2]),t1 = t1_ex[i],t2 = T2)))
    d0_ext_2 = sapply(1:length(t1_ex),function(i) h0_hat(datap[,1],datap[,2],t1_ex[i],t2_ex[i],T1,T2))
    d1_ext_2 = sapply(1:length(t1_ex),function(i) h1_hat(datap[,1],datap[,2],data[,1],data[,2],t1_ex[i],t2_ex[i],T1,T2))
    d1_ext_2 = exp(d1_ext_2)
    d01_ext[,1] = d01_ext[,1] + d0_ext_2
    d01_ext[,2] = d01_ext[,2] * d1_ext_2
    d01 = rbind(d01_int,d01_ext)
    
    a2_q01_ext = t(sapply(1:length(t1_ex),function(i) D0_rare(xi = c(T1,data[,1]),yi = c(T2,data[,2]),t1 =0,t2 = T2)))
    a2_q0_ext_2 = t(sapply(1:length(t2_ex),function(i) h0_hat(datap[,1],datap[,2],0,t2_ex[i],T1,T2)))[,1]
    a2_q1_ext_2 = t(sapply(1:length(t1_ex),function(i) h1_hat(datap[,1],datap[,2],data[,1],data[,2],0,t2_ex[i],T1,T2)))[,1]
    a2_q1_ext_2 = exp(a2_q1_ext_2)
    a2_q01_ext[,1] = a2_q01_ext[,1] + a2_q0_ext_2
    a2_q01_ext[,2] = a2_q01_ext[,2] * a2_q1_ext_2
  }
  d2_in_ex = sapply(1:length(t1),function(i) d2_hat(x1 = data[,1],x2 = data[,2],sort(t1,decreasing = T)[i],sort(t2)[i]))
  d012 = cbind(t1 = t1, d01,d2_in_ex)
  
  if(T2 > T1) {
    t2_exceed = seq(T1,T2,length.out = 5) %>% round 
    t2_in = c(t2_in, t2_exceed) %>% unique
    t2 = c(t2, t2_exceed) %>% unique
  }
  a2_q01 = t(sapply(1:length(t2_in),function(i) D0_rare(xi = c(T1,data[,1]),yi = c(T2,data[,2]),t1 = 0,t2 = t2_in[i])))
  a2_q2 =  sapply(1:length(t2),
                  function(i) d2_hat(x1 = data[,1],x2 = data[,2],0,t2[i])) 
  if(extr == 1){ a2_q01 = rbind(a2_q01,a2_q01_ext) }
  
  
  d0 = cbind(t1 = t1, t2 = t2, a1_q01[,1], a2_q01[,1], d012[,2]);colnames(d0) = c("t1","t2",colnames(data)[1],colnames(data)[2],"Mixed")
  d1 = cbind(t1 = t1, t2 = t2, a1_q01[,2], a2_q01[,2], d012[,3]);colnames(d1) = c("t1","t2",colnames(data)[1],colnames(data)[2],"Mixed")
  d2 = cbind(t1 = t1, t2 = t2, a1_q2, a2_q2, d012[,4]);colnames(d2) = c("t1","t2",colnames(data)[1],colnames(data)[2],"Mixed")
  
  if(extr == 1){ 
    cross = which(d0[,2] == T2)
    d0 = rbind(d0,d0[cross,]);d0 = d0[order(d0[,2]),];rownames(d0) = NULL
    d1 = rbind(d1,d1[cross,]);d1 = d1[order(d1[,2]),];rownames(d1) = NULL
    d2 = rbind(d2,d2[cross,]);d2 = d2[order(d2[,2]),];rownames(d2) = NULL
  }
  
  colnames(d0)[c(1,2)] = c("m1","m2")
  colnames(d1)[c(1,2)] = c("m1","m2")
  colnames(d2)[c(1,2)] = c("m1","m2")
  if(T2 > T1){
    d0[d0[,2]>T1,c(3,5)]=NA
    d1[d1[,2]>T1,c(3,5)]=NA
    d2[d2[,2]>T1,c(3,5)]=NA
  }
  return(list(q0 = d0, q1 = d1, q2 = d2, q0_ana = q0_ana))
}

####################################################################################
#
# (2). Plot individual rarefaction/extrapolation curves and mixture curves
#
####################################################################################
#' multi.plot(data, ans1, type) plot the outcome of S_abundance(abundance data) or dq_hat(incidence data).
#' @param data the Sx2 data used in dq_hat.
#' @param ans1 the outcome of dq_hat.
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
# (3). Example
#
####################################################################################
bird = read.table("bird.txt")
result_inci_bird = dq_hat(bird, knots = 10)
multi.plot(data = bird, ans1 = result_inci_bird, type = "incidence")
