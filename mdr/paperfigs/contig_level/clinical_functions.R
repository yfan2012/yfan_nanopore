findMethFreq <- function(x) {
    ##collapses called meth into methfreq.
    ##x is meth df above, grouped by chrom and motif

    motiflen=nchar(x$motif[1])
    
    motifpos=NULL
    for (i in x$pos) {
        diffs=abs(x$pos-i)
        motifgroup=x[diffs<=motiflen,] %>%
            arrange(-methfrac) %>%
            mutate(totcalls=methnum+umethnum) %>%
            arrange(-totcalls)
        motifpos=bind_rows(motifpos, motifgroup[1,])
    }

    motifpos=unique(motifpos) %>%
        select(-totcalls)
    return(motifpos)
}


findpeaks <- function(test) {
    ##normalized coverage frequency tibble
    ##return height and time (highest point of bimodality, and time spent in bimodality)
    
    steps=rev(seq(.01,1,.01))
    covrange=diff(range(test$cov))
    
    height=0
    time=0
    spread=0
    for (i in steps) {
        over=which(test$normfreq>=i)
        if (length(over)>1) {
            peaks=1+sum(diff(over)>1)
            if (peaks>1) {
                time=time+.01
                starts=c(over[1], over[which(diff(over)>1)]+1)
                ends=c(over[which(diff(over)>1)+1], over[length(over)])
                mids=starts+(ends-starts)/2
                ispread=max(diff(mids))/covrange
                if (spread<ispread && ispread>.05) {
                    spread=ispread
                    if (height<i){
                        height=i
                    }
                }
                
            }
        }
    }
    return(tibble(tig=test$tig[1], h=height, t=time, s=spread))
}


qqmse <- function(test){
    ##poisson coverage 
    mean=sum(test$cov*test$freq)/sum(test$freq)
    std=sqrt(mean)
    
    ordered_freq=sort(test$freq)
    ptile=test %>%
        rowwise() %>%
        mutate(percentile=sum(test$freq[test$cov<=cov])/sum(test$freq)) %>%
        mutate(theoretical=pnorm(cov, mean, std)) %>%
        mutate(res=(percentile=theoretical)^2)

    mse=sqrt(sum(ptile$res))
    return(tibble(tig=test$tig[1], mse=mse))
}
