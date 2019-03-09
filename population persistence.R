##########################################################################

library("gplots")
ptime<-function(x1){
  a<-lapply(x1,function(x) x[1])
  b<-sqrt(length(a))
  
  cc<-matrix(unlist(a),b,b)
  return(cc)
  
  #heatmap.2(c,Rowv = FALSE,Colv = FALSE,col=heat.colors(10),scale = "none",key = TRUE,
  #   symkey = FALSE,density.info = "none",trace="none",cexRow = 0.9,xlab ="ii",ylab = "jj" )
  
}


rtran<-function(x){
  n<-ncol(x)
  b<-matrix(,n,n)
  for(i in 1:n){
    b[i,]<-x[(n+1-i),]
  }
  return(b)
}

ctran<-function(x){
  n<-ncol(x)
  b<-matrix(,n,n)
  for(i in 1:n){
    b[,i]<-x[,(n+1-i)]
  }
  return(b)
}
library(RColorBrewer)
#RColorBrewer

display.brewer.all(n=11)
display.brewer.all()

cols<-brewer.pal(n=9,name="Spectral")
##################################################################################
#####################################################################################
#########################################
#####################################################################################
#########################################
#####################################################################################
#########################################0

no_cores <- detectCores() - 1
cl <- makeCluster(no_cores)
registerDoParallel(cl)
mm<-rep(1:10,10)
jjj<-function(x){
  b<-0
  for (i in 1:x) {
    a<-rep(i,x)
    b<-c(b,a)
  }
  c<-b[-1]
  return(c)
}
nn<-jjj(10)

v_d<-foreach(ii=nn,jj=mm,.multicombine = TRUE)  %dopar% persistence(t=1000,beta =0.6*jj,size=50,f_mean = 0.5,f_sd = 0.3,h_mean = 0.6,h_sd = 0.3,velocity=0.001*ii)
write.csv(v_d,file = "velocity-dispersal.csv")

stopCluster(cl)

hm0<-ptime(v_d)
cname<-seq(1,10,1)
rname<-seq(0.6,6,0.6)
colnames(hm0)<-cname
rownames(hm0)<-rname



heatmap.2(hm0,Rowv = FALSE,Colv = FALSE,col=cols,scale = "none",key = TRUE,
          symkey = FALSE,density.info = "none",trace="none",key.xlab ="Persistence probability",
          key.title = TRUE, cexRow = 0.9,xlab =expression(Velocity(10^-4)),ylab = "Dispersal cost" )


write.csv(hm0,file = "hm0.csv")


hm0<-read.csv(file = "C:/Users/dell/Desktop/ÍË»¯ËÙÂÊ/hm0.csv",header = T,row.names = 1)
hm0<-as.matrix(unlist(hm0))
dim(hm0)<-c(10,10)
cname<-seq(1,10,1)
rname<-seq(0.6,6,0.6)
colnames(hm0)<-cname
rownames(hm0)<-rname
jpeg(file="velocity_dispersal.jpeg",width=200, height=200, units='mm',res=400)


heatmap.2(hm0/1000,Rowv = FALSE,Colv = FALSE,col=cols,scale = "none",key = TRUE,
          symkey = FALSE,density.info = "none",trace="none",key.xlab ="Persistence probability",
          key.title = TRUE, cexRow = 0.9,xlab =expression(Velocity(10^-4)),ylab = "Dispersal cost" )
dev.off()

#################################################

vd<-t(hm0)
extinction<-array(vd)
velocity<-rep(seq(0.0001,0.001,0.0001),10)

velp_ext<-as.data.frame(cbind(velocity,extinction))

library(ggplot2)

jpeg(file="velocity_persistence.jpeg",width=100, height=100, units='mm',res=300)


qplot(velocity,extinction/1000,data=velp_ext,geom="smooth",ylab = "Persistence probability")
dev.off()


################################################################################




velocity_1<-v_d[[51:60]]

v_d1<-as.data.frame(unlist(v_d))

dim(v_d1)
class(v_d1)


###################################################################################
###############1.dispersal-heterogeneity
no_cores <- detectCores() - 1
cl <- makeCluster(no_cores)
registerDoParallel(cl)
mm<-rep(1:50,50)
jjj<-function(x){
  b<-0
  for (i in 1:x) {
    a<-rep(i,x)
    b<-c(b,a)
  }
  c<-b[-1]
  return(c)
}
nn<-jjj(50)

d_h<-foreach(ii=nn,jj=mm,.multicombine = TRUE)  %dopar% persistence(t=1000,beta =0.12*ii,size=50,f_mean = 0.5,f_sd = 0.3,h_mean = 0.6,h_sd = jj*0.006,velocity=0.005)
write.csv(d_h,file = "dispersal-heterogeneity.csv")

stopCluster(cl)

hm1<-ptime(d_h)

heatmap.2(hm1,Rowv = FALSE,Colv = FALSE,col=heat.colors(10),scale = "none",key = TRUE,
          symkey = FALSE,density.info = "none",trace="none",cexRow = 0.9,xlab ="dispersal",ylab = "heterogeneity" )
write.csv(hm1,file = "hm1.csv")

##################################################################################
##########2.heterogeneity-dv

no_cores <- detectCores() - 1
cl <- makeCluster(no_cores)
registerDoParallel(cl)
mm<-rep(1:50,50)
jjj<-function(x){
  b<-0
  for (i in 1:x) {
    a<-rep(i,x)
    b<-c(b,a)
  }
  c<-b[-1]
  return(c)
}
nn<-jjj(50)

h_dv<-foreach(ii=nn,jj=mm,.multicombine = TRUE)  %dopar% persistence(t=1000,beta =4,size=50,f_mean = 0.5,f_sd =ii*0.006,h_mean = 0.6,h_sd = jj*0.006,velocity=0.005)
write.csv(h_dv,file = "heterogeneity-dv.csv")
stopCluster(cl)

hm2<-ptime(h_dv)
heatmap.2(hm2,Rowv = FALSE,Colv = FALSE,col=heat.colors(10),scale = "none",key = TRUE,
          symkey = FALSE,density.info = "none",trace="none",cexRow = 0.9,xlab ="dv",ylab = "heterogeneity" )
write.csv(hm2,file = "hm2.csv")
##################################################################################
##########3.size-dv
no_cores <- detectCores() - 1
cl <- makeCluster(no_cores)
registerDoParallel(cl)
mm<-rep(1:50,50)
jjj<-function(x){
  b<-0
  for (i in 1:x) {
    a<-rep(i,x)
    b<-c(b,a)
  }
  c<-b[-1]
  return(c)
}
nn<-jjj(50)

s_dv<-foreach(ii=nn,jj=mm,.multicombine = TRUE)  %dopar% persistence(t=1000,beta =4,size=1*jj,f_mean = 0.5,f_sd = ii*0.006,h_mean = 0.6,h_sd = 0.3,velocity=0.005)
write.csv(s_dv,file = "size-dv.csv")
stopCluster(cl)

hm3<-ptime(s_dv)
heatmap.2(hm3,Rowv = FALSE,Colv = FALSE,col=heat.colors(10),scale = "none",key = TRUE,
          symkey = FALSE,density.info = "none",trace="none",cexRow = 0.9,xlab ="dv",ylab = "size" )
write.csv(hm3,file = "hm3.csv")


##################################################################################
##########4.size-dispersal
no_cores <- detectCores() - 1
cl <- makeCluster(no_cores)
registerDoParallel(cl)
mm<-rep(1:50,50)
jjj<-function(x){
  b<-0
  for (i in 1:x) {
    a<-rep(i,x)
    b<-c(b,a)
  }
  c<-b[-1]
  return(c)
}
nn<-jjj(50)

s_d<-foreach(ii=nn,jj=mm,.multicombine = TRUE)  %dopar% persistence(t=1000,beta =0.12*ii,size=1*jj,f_mean = 0.5,f_sd = 0.3,h_mean = 0.6,h_sd = 0.3,velocity=0.005)
write.csv(s_d,file = "size-dispersal.csv")
stopCluster(cl)
hm4<-ptime(s_d)
heatmap.2(hm4,Rowv = FALSE,Colv = FALSE,col=heat.colors(10),scale = "none",key = TRUE,
          symkey = FALSE,density.info = "none",trace="none",cexRow = 0.9,xlab ="dispersal",ylab = "size" )

write.csv(hm4,file = "hm4.csv")



m1<-ctran(hm1)
m2<-hm2
m3<-hm3
m4<-ctran(hm4)
mm<-rbind(cbind(m4,m3),cbind(m1,m2))

heatmap.2(mm,Rowv = FALSE,Colv = FALSE,col=heat.colors(8),scale = "none",key = TRUE,
          symkey = FALSE,density.info = "none",trace="none",cexRow = 0.9,xlab ="Population demography",ylab = "Landscape heterogeneity" )

write.csv(mm,file = "heatmapmatrix.csv")


#####################################

no_cores <- detectCores() - 1
cl <- makeCluster(no_cores)
registerDoParallel(cl)
mm<-rep(1:10,10)
jjj<-function(x){
  b<-0
  for (i in 1:x) {
    a<-rep(i,x)
    b<-c(b,a)
  }
  c<-b[-1]
  return(c)
}
nn<-jjj(10)

v_dif<-foreach(ii=nn,jj=mm,.multicombine = TRUE)  %dopar% persistence(t=1000,beta =4,size=50,f_mean = 0.5,f_sd = 0.03*jj,h_mean = 0.6,h_sd = 0.3,velocity=0.001*ii)
write.csv(v_dif,file = "velocity-difference.csv")

stopCluster(cl)

hm5<-ptime(v_dif)

heatmap.2(hm5,Rowv = FALSE,Colv = FALSE,col=cols,scale = "none",key = TRUE,
          symkey = FALSE,density.info = "none",trace="none",cexRow = 0.9,xlab ="velocity",ylab = "individual difference" )
c4write.csv(hm5,file = "hm5.csv")


