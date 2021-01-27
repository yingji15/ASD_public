# sfari enrich
# 12/10/2020
# validation by enrichment in sfari
# take argument: gene-score list
# return sfari enrichment values
.libPaths("/home/jiy1/R/rlib-3.6.0")
library(optparse)

option_list <- list(
    ## Type : logical, integer, souble complex, character
    make_option(c("--list"), action="store", #dimension
                type='character', help="gene scoring list: 1st column gene name, 2nd column score")
)

library(readxl)
library(dplyr)

opt <- parse_args(OptionParser(option_list=option_list))
print(opt)

score.file = opt$list

score = read.table(score.file,head=T,as.is=T)


fname = sapply( strsplit(score.file,"[.]"),'[',1)
print( fname)




# get predicted genes + results
load("/home/jiy1/ASD/ASD_paper/other_methods/sfari_lists.RData")

enasd2<-function(l,r, score=score, sel=netdis_bf, testasd = testasd){
    #test<-arrange(comb2,desc(bf_dis))[l:r,"gene"]
    
    test<-arrange(score,desc({{sel}}))[l:r,"Gene"]
    x = r-l+1
    a<-sum(test %in% testasd )
    pop<-setdiff(score$Gene,test )
    b<-sum( pop %in% testasd )
    d<-length(pop)
    t <- fisher.test( matrix(c(a,b,x-a,nrow(score)),2,2),alternative='greater')
    print ( fisher.test( matrix(c(a,b,x-a,nrow(score)-x-b),2,2),alternative='greater')$p.value )
    print ( fisher.test( matrix(c(a,b,x-a,nrow(score)-x-b),2,2),alternative='greater'))
    print(c(a,b,x-a,nrow(score)-x-b ))
    list = list("t"=t,"a"=a,'b'=b,'x_a'=x-a,'n_x_b'=nrow(score)-x-b )
    return(list)
}

enasd2(1,100,score=score,sel=Score,s10)

cuts2 <- data.frame(l=c(1,1,1,1,1001,2001,3001),r=c(100,200,500,1000,2000,3000,4000))
out2 <- matrix(0,nrow=nrow(cuts2), ncol=17)
for (i in 1:nrow(cuts2)){
    out2[i,1] <- lb <- cuts2[i,1]
    out2[i,2] <- rb <- cuts2[i,2]
    #num = cuts[i]
    # tier 1 
    o1 <- enasd2(l=lb,r=rb,score=score,sel=Score,s10)
    t <- o1$t
    out2[i,3] <- t$p.value
    out2[i,4] <- t$estimate
    out2[i,5] <- t$conf.int[1]
    out2[i,6] <- t$conf.int[2]
    out2[i,7] <- o1$a
    # tier 2
    
    o1 <- enasd2(l=lb,r=rb,score=score,sel=Score,s20)
    t <- o1$t
    out2[i,8] <- t$p.value
    out2[i,9] <- t$estimate
    out2[i,10] <- t$conf.int[1]
    out2[i,11] <- t$conf.int[2]
    out2[i,12] <- o1$a
    
    # tier 3
    o1 <- enasd2(l=lb,r=rb,score=score,sel=Score,s30)
    t <- o1$t
    out2[i,13] <- t$p.value
    out2[i,14] <- t$estimate
    out2[i,15] <- t$conf.int[1]
    out2[i,16] <- t$conf.int[2]
    out2[i,17] <- o1$a
    
    
}

colnames(out2) <- c("lb",'rb',"t1_p-value","t1_OR","t1_OR_ci_l","t1_OR_ci_r",'t1_hits', "t2_p-value","t2_OR","t2_OR_ci_l","t2_OR_ci_r",'t2_hits',"t3_p-value","t3_OR","t3_OR_ci_l","t3_OR_ci_r",'t3_hits')


sel = c("lb",'rb',"t1_p-value","t1_OR", "t2_p-value","t2_OR", "t3_p-value","t3_OR")
out3 = out2[4:7,sel]



print(out3)

colnames(out3)[3:8] <- paste(sel[3:8],fname,sep="_") 

print(out3)


save(out3,file=paste0(fname,".sfari.RData") )
# 

# THEN edit the latex output in overleaf
