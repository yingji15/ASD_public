# add average network distance to positive genes to BFs

# Code adapted from Quan Wang

library(dplyr)

# positive seed gene 

load("tada65.RData")
 


# get average network distance to positive seed genes

# this is the matrix generated from step 2
load("go_propogation_probality_rp_0.3.RData")
colnum<-which(colnames(pro_p) %in% tada65)
tadadis<-pro_p[,colnum] #only those selected columns
tadadissum<-rowSums(tadadis)


# load new scores (combined binary and continuous BFs)

load("bin_ct_comb_092120_out.RData")

for (i in 1:nrow(comb) ){
 if ( comb[i,1] %in% tada65 ){ 
	other<- setdiff(tada65,comb[i,1] )
	a = pro_p[which( rownames(pro_p) ==comb[i,1]) ,which(colnames(pro_p)%in% other) ]
	comb[i,c("netdis")]<- mean(a)
 }else{
 
 	a = pro_p[which( rownames(pro_p) ==comb[i,1]) ,which(colnames(pro_p)%in% tada65) ]
  comb[i,c("netdis")]<- mean(a)
 }
 }
 
 

#standarize network distance

range01 <- function(x){(x-min(x))/(max(x)-min(x))}
comb2 = comb[!is.nan(comb$netdis),]
comb2$stdis<-range01(comb2$netdis)



#multiply distance with other feature products
comb2[,"bf_dis"] = comb2$bin2_ct2*comb2$stdis


# save
save(comb2, file="net_bf.RData")
