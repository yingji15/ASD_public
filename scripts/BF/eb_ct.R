# get BFs using parameters estimated by Empirical bayes 
# include both binary and continuous features
# give output: genome wide score

# compute scores only
# Rscript eb_ct.R

# compute and write score results (include all genes)
# Rscript eb_ct.R  --write
suppressMessages(require(optparse))

option_list <- list(
    
    make_option(c("-v", "--verbose"), action="store_true", default=FALSE,
                help="Print extra output [default]"),
    make_option(c("--write"), action="store_true", 
                dest="verbose", help="Print little output")
    
)


opt <- parse_args(OptionParser(option_list=option_list))

print(opt$verbose)



library(dplyr)


#####################
# load data
#######################


# processed data

tab = read.table("aut_feature_0628.txt",as.is=T,head=T)
dim(tab)


# get positive training genes
load("tada65.RData")



# get random negative genes

neg = sample(tab$GeneSymbol[!tab$GeneSymbol %in% tada65],500)

# now form a 'other' tab for all testing genes (not training)
other = tab[!tab$GeneSymbol %in% c(tada65,neg),]


######################################
# continuous features
######################################



# selected features

ctname = c("HIS","GHIS","pLI","RVIS.pop_maf_0.05.","ncRVIS")
print(ctname)

# preprocess data, center and scale use negative training data
# resulting data: testc, otherc, posc, negc


otherc = other[,colnames(other) %in% c("GeneSymbol",ctname)]
#train
posc = tab[tab$GeneSymbol %in% tada65, colnames(tab) %in% c("GeneSymbol",ctname)]
negc = tab[tab$GeneSymbol %in% neg,colnames(tab) %in% c("GeneSymbol",ctname)]

# method of moments estimator
# posc: matrix with training genes
eb <- function(i=i, posc = posc){
    
    # i: indicator for which column to process
    a1 = mean(posc[,i])
    a2 = mean(posc[,i]**2)
    a3 = mean(posc[,i]**3)
    a4 = mean(posc[,i]**4)
    
    # hyperparameters
    m0 = a1
    k0 = 1
    up0 = (-(14/3)*(a1**4)+4*(a1**2)*a2 + 2*(a2**2)-(4/3)*a4)/(-(2/3)*(a1**4)+(a2**2)-(1/3)*a4)
    s0 = (-(1/6)*(a2-(a1**2))*(5*(a1**4)-6*(a1**2)*a2+a4 ))/(-(7/3)*(a1**4)+2*(a1**2)*a2+(a2**2)-(2/3)*a4)
    
    
    # parameters based on hyperparameters
    loc1 <- m0
    scale1 <- sqrt(s0*2)
    
    list = list("loc1"=loc1,"scale1"=scale1,"df"=up0,'m0'=m0,'k0'=k0,'up0'=up0,'s0'=s0 )
    return(list)
    
    
}

# example of one column
hisp = eb(i=2,posc=posc)




out <- data.frame(matrix(NA,nrow=nrow(otherc),ncol=length(ctname)))

for ( c in 2:ncol(posc)){
    
    
    print(ctname[c-1])
    # m1 posterior
    n1 = eb(i=c,posc=posc)
    
    # m0 posterior
    n2 = eb(i=c,posc=negc)
    
    # any gene without label:
    for (i in 1:nrow(otherc)){
        thisx = otherc[i,c]
        # probability of observing this data under m1
        # t distribution: "loc1"=loc1,"scale1"=scale1,"df"=up0
        stdVals1 <- (thisx - n1$loc1)/(n1$scale1)
        n1f <- dt(stdVals1,df=n1$df)/(n1$scale1)
        
        # probability of observing this data under m1
        # t distribution: "loc1"=loc1,"scale1"=scale1,"df"=up0
        stdVals1 <- (thisx - n2$loc1)/n2$scale1
        n0f <- dt(stdVals1,df=n2$df)/n2$scale1
        
        # BF: ratio of the two probabilities
        out[i,c] = n1f/n0f
    }
    
}


ci_other1 = out[,tokeep]
ci_other2 = data.frame(ci_other1) %>% mutate(prod=Reduce('*',.)) %>% as.data.frame()

ci_other2[,"gene"] = other$GeneSymbol


colnames(ci_other2)[1:(ncol(ci_other2)-2)]=ctname


######################################
# binary features
######################################

# selected binary features
sel = read.table("lr_0825_colnames.txt",as.is=T,head=F)

ctname2 = c("HIS","GHIS","pLI","RVIS.pop_maf_0.05.","ncRVIS","pcGERP_percentile")
biname = setdiff(sel[,1], c(ctname2,"gene"))
biname = setdiff(biname,"prod13")



# selected binary features



otherb = other[,colnames(other) %in% c("GeneSymbol",biname)]
#train
posb = tab[tab$GeneSymbol %in% tada65, colnames(tab) %in% c("GeneSymbol",biname)]
negb = tab[tab$GeneSymbol %in% neg,colnames(tab) %in% c("GeneSymbol",biname)]

ebb <- function(i=i, posb = posb){
    
    a1 = sum(posb[,i])
    
    
    n = length(posb[,i])
    
    
    # parameter
    #theta = astar/(astar+bstar)
    astar = (a1/n)/(1-a1/n)
    bstar = 1

    list = list("astar"=astar,"bstar"=bstar)
    return(list)
    
    
}



# helper function for beta-binomial

# (n choose x) beta(x+u,n-x+v)/beta(u,v)
dbb <- function(x,a,b){
    
    o=beta(x+a,1-x+b)/beta(a,b)
    return(o)
}
    



outb <- data.frame(matrix(NA,nrow=nrow(otherb),ncol=length(biname)))

for ( c in 2:ncol(posb)){
    
    
    print(biname[c-1])
    # m1 posterior
    n1 = ebb(i=c,posb=posb)
    
    # m0 posterior
    n2 = ebb(i=c,posb=negb)
    
    
    for (i in 1:nrow(otherb)){
        thisx = otherb[i,c]
        
        
        n1f = dbb(x=thisx,a=n1$astar,b=n1$bstar)
        n0f = dbb(x=thisx,a=n2$astar,b=n2$bstar)
        # BF
        outb[i,c] = n1f/n0f
    }
    
}

# only keep those with BF not equal to NA 
tokeep = which(sapply(outb, function(x)all(!is.na(x))))


bi_other1 = outb[,tokeep]
bi_other2 = data.frame(bi_other1) %>% mutate(prod=Reduce('*',.)) %>% as.data.frame()

bi_other2[,"gene"] = otherb$GeneSymbol




colnames(bi_other2)[1:(ncol(bi_other2)-2)]=biname[tokeep-1]


######################################
# output
######################################

# merge output scores
cmb = merge(bi_other2,ci_other3,by="gene",suffixes = c(".bi",".ct"))


if (opt$verbose ){
    
    save(cmb,tada65,file="output/cmb_bi_ct_102220.RData")
    

}
