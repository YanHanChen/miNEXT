####################################################################################
#
#Plot individual rarefaction/extrapolation curves and mixture curves
#
####################################################################################
#' multi.plot(data, ans1, type) plot the outcome of Abundance(abundance data) or Incidence(incidence data).
#' @param data the Sx2 data used in Abundance/Incidence.
#' @param ans1 the outcome of Abundance/Incidence.
#' @param type datatype of data, "abundance" or "incidence". Default is "abundance".
#' @return a list containing two plots: the first one for diversity of q = 0, 1, 2 and the second one for species composition.
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
    if(ncol(output) != 7){
      output = output %>% cbind(.,UCL = .[,5] ) %>% cbind(., LCL = .[,5])
    }
    lty1 = cbind(as.character(output[,1]<=n1),as.character(output[,2]<=n2))
    index1 = 1
    output = as.data.frame(output)
    output$m_main = output[, index1]
    output1 = data.frame(m = c(output$m1, output$m2, output$m_main), value = c(output[, 3], output[, 4], output[, 5]))
    nn = length(output[, 3])
    name1 = names(output)[c(3, 4)]
    output1$type1 = factor(c(rep(name1[1], nn), rep(name1[2], nn), rep("Mixture", nn)), levels = c(name1[1],name1[2],"Mixture"), labels = c(name1[1],name1[2],"Mixture"))
    
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
    output1 = output1[!is.na(output1$value),] %>% cbind(.,LCL = output$UCL) %>% cbind(.,UCL = output$LCL)
    output1[output1$type1!="Mixture",c(5,6)] =  output1[output1$type1!="Mixture",2]
    output2 = output2[!is.na(output2$va),]
    
    pp = ggplot(output1)+
      geom_hline(aes(yintercept = output[,4][output[, 2]==max(output$m2)]), col = "darkgray", linetype = 3, size = 1.25)+
      geom_hline(aes(yintercept = output[,3][output[, 1]==max(output$m_main)][1]), col = "darkgray", linetype = 3, size = 1.25)+
      geom_line(aes(x = m, y = value, col = type1, size = type1,linetype = type2))+
      scale_size_manual(breaks=c(name1[1],name1[2],"Mixture"), values=c(1, 1, 3),
                        labels = c(name1[1],name1[2],"Mixture"))+
      scale_color_manual(values = c("black", "#3498DB","#C0392B"),
                         breaks = c(name1[1],name1[2],"Mixture"),
                         labels = c(name1[1],name1[2],"Mixture"))+
      guides(linetype = FALSE)+
      xlim(c(0, max(output1$m)+5))+
      ylab("Diversity")+
      theme_bw()+
      theme(panel.grid.major=element_blank(),panel.grid.minor=element_blank())+
      theme(legend.position="none",legend.title=element_blank(),
            legend.text=element_text(size=15),legend.key.width  = unit(1.5,"cm"))+
      theme(plot.title = element_text(size=14, face="bold.italic",hjust = 0))+
      theme(axis.text.x= element_text(size = 18,colour = "black",margin=unit(c(0.2,0.2,0.5,0.5), "cm")),
            axis.text.y= element_text(size = 18,colour = "black",margin=unit(c(0.2,0.2,0.2,0.2), "cm")),
            axis.title.x=element_text(size = 20,lineheight = 1.1),axis.title.y=element_text(size = 16),
            axis.ticks = element_line(size = 1.5,colour = "black"),axis.ticks.length = unit(-0.25,"line"))+
      scale_x_continuous(name = paste0("Number of ", text1, "\n(Proportion % in original habitat)"),
                         breaks = round(seq(0,max(output$m_main),length.out = 6)),
                         labels = paste0(round(seq(0,max(output$m_main),length.out = 6)), "\n(", round(seq(0,max(output$m_main),length.out = 6)/(max(output$m_main))*100, 1), "%)")) + 
      geom_ribbon(aes(x = m, ymin = LCL, ymax = UCL, fill = type1),alpha = 0.2)+
      scale_fill_manual(values = c("black", "#3498DB","#C0392B"),
                        breaks = c(name1[1],name1[2],"Mixture"),
                        labels = c(name1[1],name1[2],"Mixture"))
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
    
    if(ncol(output_d) != 7){
      output_d = output_d %>% cbind(.,UCL = .[,5] ) %>% cbind(., LCL = .[,5])
      tmp <- output_p[,rep(3:5,2)]
      colnames(tmp) <- c(paste0("LCL.",p_names),paste0("UCL.",p_names))
      output_p = cbind(output_p,tmp)
    }
    
    output_p = output_p[,c(-2)] %>% mutate(lty={
      (output_p$m1 >= n1-n2) %>% as.numeric()
    })
    
    output_d <- data.frame(output_d[,-c(2,4)]) %>% mutate(lty={
      (.$m1 >= n1-n2) %>% as.numeric()
    })
    if(sum(output_d$m1==(n1-n2))==2 ){
      corss <- which(output_d$lty==1) %>% max
      output_d$lty[corss]=0
      # output_p$lty[corss]=0
    }
    
    
    output_p <- data.frame(m1 = rep(output_p[,1],3), lty = rep(output_p[,"lty"],3),
                           col = rep(colnames(output_p)[2:4],each = nrow(output_p)),
                           est = c(output_p[,2],output_p[,3],output_p[,4]),
                           LCL = c(output_p[,5],output_p[,6],output_p[,7]),
                           UCL = c(output_p[,8],output_p[,9],output_p[,10]))
    # output_p <- reshape(data = output_p, direction = "long", varying = colnames(output_p)[-c(1,11)],
    #         timevar='col', times= c("U_forest","U_logged","Share"), v.names=c("EST",'LCL', 'UCL'),
    #         idvar=c('m1',"lty"), new.row.names = NULL)
    # rownames(output_p) = NULL
    
    output_d = melt(data = output_d[,c(1:3,6)],id.vars = c("m1","lty"),variable.name = "col", value.name = "est")
    output_d[output_d$col== sites[1],2] = 1
    output_d$LCL <- output_d$est
    output_d$UCL <- output_d$est
    
    output_d$col <- factor(x = output_d$col,levels = c(sites[1],"Mixture"))
    output_d$lty <- paste0(output_d$lty, output_d$col) %>% factor(.,levels =c(paste0("1",sites[1]),
                                                                              "1Mixture","0Mixture"))
    
    output_p$col <- factor(x = output_p$col,levels = unique(output_p$col))
    output_p$lty <- paste0(output_p$lty, output_p$col) %>% factor(.,levels = unique(.))
    
    output = rbind(output_d,output_p)
    pp = ggplot(output)+
      theme_bw()+
      geom_hline(aes(yintercept = output[output$col=="Mixture",4][1]), col = "darkgray", linetype = 3, size = 1.25)+
      geom_line(data = output, aes(x = m1, y = est, linetype = lty, col = col, size = lty))+
      scale_size_manual( breaks = levels(output$lty),
                         values=c(1.5, 3, 3, 2, 2, 2))+
      scale_color_manual(values=c("black", "#C0392B","black","#3498DB","purple"),
                         breaks = c(sites[1],"Mixture",p_names[3],p_names[1],p_names[2]))+
      scale_linetype_manual(breaks = levels(output$lty),
                            values=c("solid","solid","2121","4111","311111","1111"))+
      xlim(c(0, max(output_d$m1)+5))+
      ylab("Diversity")+
      theme(panel.grid.major=element_blank(),panel.grid.minor=element_blank())+
      theme(legend.position="none",legend.title = element_blank(),
            legend.text=element_text(size=15),legend.key.width  = unit(1.5,"cm"))+
      theme(axis.text.x= element_text(size = 18,colour = "black",margin=unit(c(0.2,0.2,0.2,0.2), "cm")),
            axis.text.y= element_text(size = 18,colour = "black",margin=unit(c(0.2,0.2,0.2,0.2), "cm")),
            axis.title.x=element_text(size = 20),axis.title.y=element_text(size = 20),
            axis.ticks = element_line(size = 1.5),axis.ticks.length = unit(-0.25,"line"))+
      scale_x_continuous(name = paste0("Number of ", text1, "\n(Proportion % in original habitat)" ),
                         breaks = round(seq(0,max(output_d$m1),length.out = 6)),
                         labels = paste0(round(seq(0,max(output_d$m1),length.out = 6)), "\n(", round(seq(0,max(output_d$m1),length.out = 6)/(max(output_d$m1))*100, 1), "%)"))+
      geom_ribbon(aes(x = m1, ymin = LCL, ymax = UCL, fill = col), alpha = 0.1)+
      scale_fill_manual(values=c("black", "#C0392B","black","#3498DB","purple"),
                        breaks = c(sites[1],"Mixture",p_names[3],p_names[1],p_names[2]))
  }
  
  samey = max(ans1[[2]][,-c(1,2)])
  p0 = draw.f_Div(data1,ans1[[1]])+ggtitle("q=0, species richness")+theme(plot.title = element_text(size = 20, face = "bold"))
  p1 = draw.f_Div(data1,ans1[[2]])+ggtitle("q=1")+theme(plot.title = element_text(size = 20, face = "bold"))
  p2 = draw.f_Div(data1,ans1[[3]])+ggtitle("q=2")+theme(plot.title = element_text(size = 20, face = "bold"))+ylim(c(0,samey))
  p0_ana = plot_comb1(data1 = data1,output_d = ans1[[1]],output_p = ans1[[4]])+
    # ylim(c((min(unlist(ans1[[1]][, -(1:2)]))), ceiling(max(unlist(ans1[[1]][, -(1:2)])))))+
    ggtitle("q=0, species composition")+theme(plot.title = element_text(size = 20, face = "bold"))
  g1 = ggarrange(p0, p1, p2, ncol=1, nrow=3, common.legend = TRUE, legend="bottom",labels = c("(a)", "(b)", "(c)"),font.label = list(size = 20))
  print(g1)
  print(p0_ana)
}