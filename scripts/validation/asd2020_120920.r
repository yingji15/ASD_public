# validation by enrichment in asd2020
# take argument: gene list
# return asd2020 enrichment

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


asds = read.table("asd2020.txt",as.is=T,head=T)
score.file = opt$list
score = read.table(score.file,head=T,as.is=T)



fname = sapply( strsplit(score.file,"[.]"),'[',1)
print( fname)



test <- function(l,r,score,sel){
    t1<-arrange(score,desc({{sel}}))[l:r,"Gene"]
    return(t1)
}

test(l=1,r=10,score,Score)

tasd2020<-function(l,r, score=score, sel=netdis_bf, evi=asds){
    # use embracing when wrapping in a function;
    test<-arrange(score,desc({{sel}}))[l:r,"Gene"]
    
    #
    pop<-setdiff(score$Gene,test )
    
    x = r-l+1
    # enrich in asd genes
    a1 = sum(asds[asds$gene %in% test,"ASC102_2018"])
    
    b1<-sum( asds[asds$gene %in% pop,"ASC102_2018"] )
    x = r-l+1
    d<-length(pop)
    t1 <- fisher.test( matrix(c(a1,b1,x-a1,d-b1),2,2),alternative='greater')
    
    
    
    
    list = list( 'asd102'=t1, 'a1'=a1,'b1'=b1,'x'=x,'d'=d)
    return(list)
}






# summarize those into table
cuts2 <- data.frame(l=c(1,101,201,501,1001,2001,3001,1),r=c(100,500,500,1000,2000,3000,4000,1000))
out3 <- matrix(0,nrow=nrow(cuts2), ncol=7)
for (i in 1:nrow(cuts2)){
    out3[i,1] <- lb <- cuts2[i,1]
    out3[i,2] <- rb <- cuts2[i,2]
    # tier 1 
    o1 <- tasd2020(l=lb,r=rb,score=score,sel=Score,evi=asds)
    
    
    asd102 = o1$asd102
    out3[i,3] <- asd102$p.value
    out3[i,4] <- asd102$estimate
    out3[i,5] <- asd102$conf.int[1]
    out3[i,6] <- o1$a1
    out3[i,7] <- o1$b1
    

    
    
    
    
    
}


part = c("asd102_p","asd102_or","asd102_conf","asd102_cand","asd102_bg")

colnames(out3) <- c("lb","ub",paste(part,fname,sep="_") )

print(out3)


save(out3,file=paste0(fname,".asd2020.RData") )


