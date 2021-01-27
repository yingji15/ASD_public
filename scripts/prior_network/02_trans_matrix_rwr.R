# Code adapted from Quan Wang
 require("igraph") 

##################################
# load adjacency matrix from GO
####################################
 network<-read.delim("GO.matrix",as.is=T)
 
 
##################################
# get transition matrix from adjacency matrix
####################################
 
 network<-network[,-1]
 rownames(network)<-colnames(network)

 for(i in 1:ncol(network))
  network[i,i]<-0
 
 for(i in 1:ncol(network))
  if(sum(network[,i])>0)
   network[,i]<-network[,i]/sum(network[,i]) 

 adj_m<-as.matrix(network)
 
##################################
# derive the reaching probabilities between genes via sampling strategies
##################################

 restart_p<-0.3
 pro_p<-numeric()
 candidate<-colnames(adj_m)

 t0<-proc.time()
 for(i in 1:length(candidate))
 #for(i in 1:10)
 {
  cat("i:",i,"\n",sep="")
  s<-rep(0,nrow(adj_m))
  s[which(colnames(adj_m)==candidate[i])]<-1

  s0<-s
  s1<-(1-restart_p)*as.numeric(adj_m%*%s0)+restart_p*s
  #j<-1
  while( (sum((s1-s0)^2))>1e-6 )
  {
   s0<-s1;
   s1<-(1-restart_p)*as.numeric(adj_m%*%s0)+restart_p*s
   #j<-j+1
  }
  pro_p<-cbind(pro_p,s1)
 }
 proc.time()-t0

 colnames(pro_p)<-candidate
 rownames(pro_p)<-rownames(adj_m)
 #pro_p<-pro_p[,match(unique(colnames(pro_p)),colnames(pro_p))]
 save(pro_p,file=paste("go_propogation_probality_rp_",restart_p,".RData",sep=""))

 ####============== end of generating propagation probilities ======####

 q("no")


