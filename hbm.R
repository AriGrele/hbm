library(stringr);library(jagsUI);library(ggplot2)
#################################################################################
model_data=setClass('model_data',slots=list(
  'input'     ='call',
  'dist'      ='character',
  'model_dir' ='character',
  'model_name'='character',
  'variables' ='vector',
  'vars'      ='list',
  'model'     ='character',
  'model_data'='list',
  'scales'    ='vector',
  'save'      ='vector',
  'filter'    ='data.frame'))
run_model=setClass('run_model',slots=list(
  'n.adapt'   ='numeric',
  'n.burnin'  ='numeric',
  'n.iter'    ='numeric',
  'n.thin'    ='numeric',
  'n.chains'  ='numeric',
  'jags_model'='jagsUI',
  'output'    ='data.frame'))
hbm_data=setClass('hbm_object',slots=list(
  'name'        ='character',
  'dir'         ='character',
  'data'        ='data.frame',
  'source_model'='character'),
  contains=c('model_data','run_model'))
semb_data=setClass('semb_object',slots=list(
  'name'        ='character',
  'dir'         ='character',
  'data'        ='data.frame',
  'source_model'='character'),
  contains=c('model_data','run_model'))
#################################################################################
pd=function(x){                                                                 #calculates percent of vector with same sign as vector median
  side=sign(median(x))                                                          #sign of median
  return(sum(side*x>0)/length(x))}                                              #return probability of direction
#################################################################################
setGeneric("summary", function(hbm) standardGeneric("summary"))
setMethod("summary","hbm_object",function(hbm){                                 #creates summary aggregation of model output
  l=1+length(hbm@vars[[2]])
  p=list()
  for(i in 1:l){p[[i]]=hbm@output[,i]}
  a=with(hbm@output,aggregate(response,p,pd))
  a$mean=with(hbm@output,aggregate(response,p,mean))$x
  a$CIl=with(hbm@output,aggregate(response,p,function(x) quantile(x,0.025)))$x
  a$CIu=with(hbm@output,aggregate(response,p,function(x) quantile(x,0.975)))$x
  a=setNames(a,c(names(hbm@output)[-c(length(hbm@output))],
                 'PD','mean','CI_lower','CI_upper'))
  return(a)})
#################################################################################
setGeneric("traces", function(hbm,cull) standardGeneric("traces"))
setMethod("traces","hbm_object",function(hbm,cull=0){                           #create traceplots from hbm output, takes arguments of model, directory to save to, amount to cull from output
  samples=hbm@jags_model$samples                                                #select samples and columns
  cols=colnames(hbm@jags_model$samples[[1]])
  allucols=unique(gsub('\\[.+]$','',cols));ucols=allucols[1]                    #correct colnames
  for(u in allucols){                                                           #for each col:
    stop=0
    for(q in ucols){if(regexpr(u,q)>0|regexpr(q,u)>0){stop=1}}                  #   stop when non-unique column reached
    if(stop==0){ucols=c(ucols,u)}}                                              #   append col after stop
  
  nchains=hbm@jags_model$mcmc.info$n.chains                                     #extract mcmc info from model
  xmax=hbm@jags_model$mcmc.info$n.samples/nchains
  plots=list();k=0                                                              #empty vars for output
  for(u in ucols){                                                              #for each column:
    out=list()
    for(i in cols[regexpr(u,cols)>0]){                                          
      k=k+1
      e=stringr::str_extract_all(i,'(\\[\\d+,\\d+\\])|(\\[\\d+\\])')
      pos=eval(parse(text=gsub('\\]','\\)',gsub('\\[','c\\(',e))))
      rhat=ifelse(length(pos)==1,hbm@jags_model$Rhat[[gsub('\\[.+]$','',i)]][pos[1]],
                  ifelse(length(pos)==0,hbm@jags_model$Rhat[[gsub('\\[.+]$','',i)]],
                         hbm@jags_model$Rhat[[gsub('\\[.+]$','',i)]][pos[1],pos[2]]))
      m=suppressMessages(reshape2::melt(as.data.frame(
        coda::as.array.mcmc.list(samples[,i])))[c(T,rep(F,cull)),])
      out[[i]]=ggplot(m)+
        geom_line(aes(x=rep(1:xmax,nchains)[c(T,rep(F,cull))],
                      y=value,color=variable),alpha=1/nchains*2,size=1)+
        theme_classic()+
        theme(text=element_text(size=65),
              panel.border=element_rect(fill=F,size=5),
              axis.text = element_blank())+
        ggtitle(i,subtitle =paste('Rhat =',round(rhat,4)))+
        guides(color='none')+
        xlab('')+ylab('')+
        scale_x_continuous(expand = c(0,0))+
        scale_y_continuous(expand=c(0,0))
      if(rhat<1.1){out[[i]]=out[[i]]+
        scale_color_manual(values=pals::ocean.ice(nchains))}
      else{out[[i]]=out[[i]]+
        scale_color_manual(values=pals::ocean.matter(nchains))}
      progress(k/length(cols),50)}
    n=floor(sqrt(length(out)))
    if(mean(regexpr('traces',list.dirs(hbm@dir))<0)==1){
      suppressMessages(dir.create(paste0(hbm@dir,'/traces/')))}
    png(paste0(hbm@dir,'/traces/',hbm@name,u,'.png'),1000*n,1000*length(out)/n)
    do.call(gridExtra::grid.arrange,c(out,ncol=n))
    dev.off()}})
#################################################################################
progress=function(percent,len,char='='){                                        #print progress bar to console, takes arguments of percent progress, nchar to print at 100%, character to form the bar from
  k=round(percent*len)                                                          #nchar to print at current percent
  cat('|',                                                                      #print
      rep(char,k),
      rep(' ',len-k),'| ',percent*100,'%                  \r',sep='')
  if(percent==1){cat('\n')}}                                                    #new line at 100%
#################################################################################
expand=function(x){o=c();for(i in x){o=c(o,1:i)};return(o)}
#################################################################################
gap=function(l,u,count){                                                        #returns missing values, takes argument of lower range vector, upper value, count
  o=c(0);n=u[1];i=1
  for(q in 2:length(l)){
    if(l[q]==0&(l[q-1]==count[i]|l[q-1]==0)){i=i+1;n=u[i];o=c(o,0)}
    else{o=c(o,n)}}
  return(o)}
#################################################################################
repframe=function(d,n){                                                         #acts like rep() for data.frames, takes arguments of data.frame, count
  out=setNames(as.data.frame(matrix(nrow=n,ncol=length(d))),names(d))
  for(i in 1:length(d)){out[,i]=rep(d[,i],n)}
  return(out)}
#################################################################################
chunk=function(ar,len){                                                         #breaks array into equal chunks, takes arguments of array, chunk length
  o=list();c=0                                                                  #define vars
  for(i in len*0:(length(ar)/len-1)){i=i+1;c=c+1;o[[c]]=ar[(i):(i+len-1)]}      #split into groups of set length, append to list
  return(o)}
#################################################################################
response=function(hbm){
  if(!(hbm@dist%in%c('dnorm'))){return()}
  return(switch(hbm@dist,
         'dnorm'=str_interp('${hbm@dist}(mu[i],tau)')))}
#################################################################################
ilen=function(data,v,x=1){                                                      #returns counts of factors in each level of hierarchy, take argument of data.frame, vector of col names
  o=c()                                                                         #define var
  d=as.numeric(as.factor(data[,v[x]]))                                          #characters to numbers
  if(x==length(v)){return(length(unique(d)))}                                   #return on final value
  else{for(i in unique(d)){o[i]=(ilen(data[d==i,],v,x+1))};return(o)}}          #otherwise repeat for sub-hierarchies 
#################################################################################
flip=function(x){                                                               #takes a named vector and returns vector of names, now named for the values they represent in original vector
  if(is.list(x)){x=x[[1]]}                                                      #edge use case
  o=data.frame('i'=NA,'name'=NA)                                                #empty data.frame
  for(i in 1:length(x)){o=rbind(o,data.frame('i'=x[i],'name'=names(x[i])))}     #iterate through items and append names to frame
  return(with(unique(na.omit(o)),setNames(name,i)))}                            #return renamed names
#################################################################################
items=function(data,v,x=1){                                                     #recursive function for generating named list from data.frame columns, take argument of dataframe, vector of col names
  o=c()                                                                         #define var
  d=as.numeric(as.factor(data[,v[x]]))                                          #characters to numbers
  d=setNames(d,data[,v[x]])                                                     #rename based on original values
  n=setNames(unique(data[,v[x]]),unique(d))                                     
  if(x==length(v)){return(d)}                                                   #return on final value
  else{for(i in unique(d)){o=c(o,items(data[d==i,],v,x+1))};return(o)}}         #otherwise repeat for sub-hierarchies
#################################################################################
semb=function(data,model,...){                                                   #write and compute Bayesian structural equation models, takes arguments of model string, data, directory to save to
  output=new('semb_object')
  param=list(...)
  defaults=list('name'        ='unnamed_semb',
                'dir'         ='.',
                'model_dir'   ='.',
                'n.adapt'     =2000,
                'n.burnin'    =1000,
                'n.iter'      =10000,
                'n.thin'      =5,  
                'n.chains'    =4)
  
  default=c(list('data' =data,
                 'input'=match.call(expand.dots=T)),
            lapply(names(defaults),\(x)ifelse(is.null(param[[x]]),
                                              defaults[[x]],
                                              param[[x]]))|>
              setNames(names(defaults)))
  for(d in names(default)){
    if(class(default[[d]])==class(slot(output,d))){slot(output,d)=default[[d]]}
    else{
      cat('Invalid input for variable "',d,'"\n',sep='')
      return()}}

  rows=str_split(model,'\n')[[1]]
  row=c();lefts=c();rights=c();old=c()
  for(r in 1:length(rows)){
    item=c(LETTERS,letters)[r]
    rights=c(rights,paste(item,'[1]',sep=''))
    groups=str_split(rows[r],'~')[[1]]
    left=paste(groups[1],'[i]',sep='')
    lefts=c(lefts,left)
    old=c(old,groups[1])
    Right=str_split(groups[2],'\\+')[[1]]
    new=c(paste(item,'[1]',sep=''))
    
    for(i in 1:length(Right)){
      old=c(old,Right[i])
      new=c(new,paste(item,'[',i+1,']*y.',Right[i],'[i]',sep=''))
      rights=c(rights,paste(item,'[',i+1,']',sep=''))}
    right=paste(new,collapse=' + ')
    row=c(row,paste(left,right,sep=' = '))}
  
  row=paste(row,collapse='\n')
  upper=paste('model {\nfor (i in 1:N){',row,'}',sep='\n')
  new=c()
  for(l in 1:length(lefts)){
    base=gsub('\\[i\\]','',lefts[l])
    new=c(new,paste('y.',lefts[l],'~dnorm(',lefts[l],',tau[',l,'])',sep=''),
          paste('y.',base,'.sim[i]~dnorm(',lefts[l],',tau[',l,'])',sep=''),
          paste('y.',base,'.res[i]=y.',lefts[l],'-',lefts[l],sep=''),
          paste('y.',base,'.sres[i]=y.',base,'.sim[i]-',lefts[l],sep=''))}
  
  row=paste(new,collapse='\n')
  middle=paste('for (i in 1:N){',row,'}',sep='\n')
  new=c()
  for(r in rights){new=c(new,paste(r,'dnorm(0.0,0.01)',sep=' ~ '))}
  
  lower=paste(new,collapse='\n')
  bottom=paste('for (j in 1:',
               length(lefts),
               '){\nsigma[j] ~ dgamma(1,1)\ntau[j] <- pow(sigma[j], -2)\n}\n}',
               sep='')
  new=c();parameters=c()
  for(l in 1:length(lefts)){
    base=gsub('\\[i\\]','',lefts[l])
    new=c(new,paste('y.',
                    base,'.fit=sum(pow(y.',base,'.res[],2)/tau[',l,'])',sep=''),
          paste('y.',
                base,'.sfit=sum(pow(y.',base,'.sres[],2)/tau[',l,'])',sep=''))
    parameters=c(parameters,
                 paste('y.',base,'.sfit',sep=''),
                 paste('y.',base,'.fit',sep=''))}
  tail=paste(new,collapse='\n')
  formula=paste(upper,middle,lower,tail,bottom,sep='\n')
  
  sink(paste0(output@dir,'/',output@name,'.jags'))
  cat(formula,fill=T)
  sink()
  
  output@model_data=list('N'=nrow(output@data))
  for(o in unique(old)){output@model_data=c(output@model_data,
                                     setNames(list(output@data[,o]),
                                              paste('y.',o,sep='')))}
  
  inits=function(){list(sigma=rep(1,length(lefts)))}  
  parameters=c(parameters,c(LETTERS,letters)[1:length(lefts)])
  output@jags_model=tryCatch(jags(output@model_data,
                                  inits,
                                  parameters,
                                  paste0(output@dir,'/',output@name,'.jags'),
                                  n.chains=output@n.chains,
                                  n.thin=output@n.thin,
                                  n.iter=output@n.iter,
                                  n.burnin=output@n.burnin,
                                  parallel=F),
                             error=function(e){print(e);return(NULL)})
  
  cat('\nMaking trace plots\n')
  traces(output)
  return(output)}
#################################################################################
setGeneric("fits", function(hbm,...) standardGeneric("fits"))
setMethod("fits","semb_object",function(hbm,...){                               #summaries fit measures for Bayesian models, takes arguments of one or more model outputs
  mc=match.call(expand.dots=T)
  mods=c(list(hbm),list(...))
  out=data.frame('response'=NA,'ppp'=NA,'DIC'=NA,
                 'names'=NA,'intercept'=NA,'slope'=NA,'r'=NA)
  if(length(mods)>1){
    for(m in 1:length(mods)){
      mc=gsub('^.+ = ','c(',mc)                                                 #bizarre functionality, if this function breaks it's probably this line
      mc=mc[mc!='fits']
      new=fits(mods[[m]])
      new$names=rep(mc[m],nrow(new))
      out=rbind(out,new)}
    out=na.omit(out)}
  else{
    sims=hbm@jags_model$sims.list
    groups=c()
    for(n in names(sims)){
      if(regexpr('\\.s?fit$',n)<=0){sims[[n]]=NULL}
      else{groups[n]=gsub('\\.s?fit$','',n)}}
    DIC=ifelse(hbm@jags_model$calc.DIC,hbm@jags_model$DIC,'NA')
    n=names(groups)
    par(mfrow=c(1,ceiling(length(n)/2)))
    g=1;for(i in 1:(length(n)/2)){
      sims=hbm@jags_model$sims.list
      x=sims[[n[g+1]]];y=sims[[n[g]]]
      l=lm(y~x)
      out=rbind(out,data.frame('response'=groups[n[g]],
                               'ppp'=pp.check(hbm@jags_model,n[g+1],n[g]),
                               'DIC'=DIC,'names'='model','intercept'=coef(l)[1],
                               'slope'=coef(l)[2],r=summary(l)$r.squared))
      g=g+2}
    par(mfrow=c(1,1))}
  return(na.omit(out))})
#################################################################################
ddic=function(mod){                                                             #compares model fits, DIC for bsems, takes argument of frame of fit outputs
  if(class(mod)!='data.frame'){cat('Not output of fits()\n')}
  else{
    if(length(unique(mod$names))<2){cat('Only one model input\n')}
    dics=unique(mod$DIC)
    n=unique(mod$names)
    out=data.frame('name'=NA,'delta_DIC'=NA)
    used=c()
    for(i in 1:length(n)){
      used=c(used,i)
      for(j in 1:length(dics)){
        if(i!=j&!(j %in% used)){
          out=rbind(out,data.frame('name'=paste(n[i],n[j],sep=' - '),
                                   'delta_DIC'=dics[i]-dics[j]))}}}
    return(na.omit(out))}}
#################################################################################
setGeneric("format_data", function(hbm) standardGeneric("format_data"))
setMethod("format_data","hbm_object",function(hbm){
  hbm@model_data=list('N'=nrow(hbm@data))
  hbm@model_data[as.character(hbm@input$model)[2]]=
    list(hbm@data[,as.character(hbm@input$model)[2]])
  full_names=c()
  
  for(n in hbm@variables){full_names=c(full_names,strsplit(n,':')[[1]])}
  
  hbm@variables=unique(full_names)
  
  for(n in unique(c(hbm@vars[[1]],hbm@variables))){
    if(!is.numeric(hbm@data[,n])){
      hbm@model_data[n]=list(as.numeric(as.factor(hbm@data[,n])))}
    else{hbm@model_data[n]=list(hbm@data[,n])}}
  
  convert=list()
  if(length(hbm@vars[[2]])>1){for(n in 2:length(hbm@vars[[2]])){
    hbm@model_data[hbm@vars[[2]][n]]=list(items(hbm@data,hbm@vars[[2]][2:n]))
    convert[[hbm@vars[[2]][n]]]=flip(hbm@model_data[[hbm@vars[[2]][n]]])
    hbm@model_data[paste('N',hbm@vars[[2]][n],sep='')]=
      list(ilen(hbm@data,hbm@vars[[2]][2:n]))}}
  hbm@scales=list()
  
  for(m in hbm@variables){
    if(length(hbm@model_data[[m]])>1&!(m%in%hbm@vars[[2]])){
      if(m!='error'){
        hbm@scales[[m]]=attributes(scale(hbm@model_data[[m]],center=F))$scale
        hbm@model_data[[m]]=as.numeric(scale(hbm@model_data[[m]],center=F))}
      else{
        hbm@scales[[m]]=1
        hbm@model_data[[m]]=as.numeric(hbm@model_data[[m]])}}}
  return(hbm)})
#################################################################################
setGeneric("format_model", function(hbm) standardGeneric("format_model"))
setMethod("format_model","hbm_object",function(hbm){
  vars=hbm@vars[[2]][-1]
  o=as.data.frame(matrix(nrow=hbm@model_data$N,ncol=2*length(vars)))|>
    setNames(rep(vars,2))
  for(i in 1:length(vars)){
    o[i]=hbm@model_data[[vars[i]]]
    o[i+2]=names(hbm@model_data[[vars[i]]])}
  filter=unique(o)
  hbm@filter=filter
  sims=hbm@jags_model$sims.list
  for(n in names(sims)){if(regexpr('\\.s?fit$',n)>0){sims[[n]]=NULL}}
  name=names(sims)
  A=name[regexec('[aA]lpha',name)>0]
  B=name[regexec('beta',name)>0]
  l=list()
  for(n in vars){l[[n]]=(hbm@model_data[paste('N',n,sep='')])[[1]]}
  
  left=setNames(as.data.frame(matrix(nrow=sum(l[[length(l)]]),ncol=length(l))),
                vars)
  if(length(l)>1){
    for(i in length(l):2){
      if(sum(!is.na(left[[i]]))==0){
        left[i]=expand(l[[i]])
        left[i][left[i]==1]=0}
      left[i-1]=gap(left[,i],expand(l[[i-1]]),as.numeric(l[i][[1]]))}
    left[length(left)][left[length(left)]==0]=1
    for(i in 1:(length(l)-1)){
      left[,i][left[,i]==0]=expand(l[[i]])}
    for(i in 1:(length(left)-1)){
      first=unique(left[1:i])
      for(q in 1:(length(left)-i)){first=cbind(first,rep('all',length(first)))}
      left=rbind(setNames(first,names(left)),left)}
    left=rbind(setNames(rep('all',length(left)),names(left)),unique(left))
    for(i in names(left)){
      left=rbind(left[left[,i]=='all',],
                 left[left[,i]!='all',][order(
                   as.numeric(left[,i][left[,i]!='all'])),])}}
  if(length(l)==1){
    left[1]=expand(l[[1]])
    left=rbind(setNames(data.frame('all'),names(left)),left)}
  p=list()
  
  out=setNames(as.data.frame(matrix(ncol=3+length(left))),
               c(names(left),c('upper','lower','response')))
  counts=list();medians=list()
  betas=B[regexpr(paste(paste('_',names(left),'$',sep=''),collapse='|'),B)<0]
  cat('\nformating slopes\n')
  for(b in betas){
    Bs=B[regexpr(b,B)>0]
    for(i in 1:nrow(left)){
      col=names(left)[regexpr('all',left[i,])<0]
      if(length(col)>0){
        lup=tryCatch(left[i,col[length(col)-1]],error=function(e) return('all'))
        if(length(lup)==0){lup='all'}
        llo=tryCatch(left[i,col[length(col)]],error=function(e) return('all'))
        upp=tryCatch(Bs[regexpr(col[length(col)-1],Bs)>0],
                     error=function(e) return('beta_all'))
        col=paste('_',col,'$',sep='')
        lab=Bs[regexpr(col[length(col)],Bs)>0]
        
        for(slab in lab){
          if(!(slab %in% names(counts))){counts[[slab]]=0}
          counts[[slab]]=counts[[slab]]+1
          right=data.frame('upper'=upp,
                           'lower'=lab,
                           'response'=medians[[lup]]+
                             as.numeric(as.data.frame(sims[slab])[counts[[slab]]][,1]))
          if(length(col)<length(left)){
            medians[[llo]]=right$response}
          if(length(left)>1){out=rbind(out,
                                       cbind(repframe(left[i,],nrow(right)),right))}
          else{out=rbind(out,cbind(
            setNames(as.data.frame(rep(left[i,],nrow(right))),
                     names(left)),right))}}}
      else{
        right=data.frame('upper'='beta_all','lower'=Bs[1],'response'=sims[[Bs[[1]]]])
        if(length(left)>1){out=rbind(out,cbind(repframe(left[i,],nrow(right)),right))}
        else{out=rbind(out,cbind(setNames(as.data.frame(rep(left[i,],nrow(right))),names(left)),right))}
        medians[['all']]=right$response[right$upper=='beta_all']}}
    progress(match(b,betas)/length(betas),50)}
  counts=list();medians=list()
  alphas=A[regexpr(paste(paste('_',names(left),'$',sep=''),collapse='|'),A)<0]
  cat('\nformating intercepts\n')
  for(a in alphas){
    Bs=A[regexpr(a,A)>0]
    for(i in 1:nrow(left)){
      col=names(left)[regexpr('all',left[i,])<0]
      if(length(col)>0){
        lup=tryCatch(left[i,col[length(col)-1]],error=function(e) return('all'))
        if(length(lup)==0){lup='all'}
        llo=tryCatch(left[i,col[length(col)]],error=function(e) return('all'))
        upp=tryCatch(Bs[regexpr(col[length(col)-1],Bs)>0],error=function(e) return('alpha_all'))
        col=paste('_',col,'$',sep='')
        lab=Bs[regexpr(col[length(col)],Bs)>0]
        for(slab in lab){
          
          if(!(slab %in% names(counts))){counts[[slab]]=0}
          counts[[slab]]=counts[[slab]]+1
          
          right=data.frame('upper'=upp,
                           'lower'=lab,
                           'response'=medians[[lup]]+
                             as.numeric(as.data.frame(sims[slab])[counts[[slab]]][,1]))
          if(length(col)<length(left)){
            medians[[llo]]=right$response}
          if(length(left)>1){out=rbind(out,cbind(repframe(left[i,],nrow(right)),right))}
          else{out=rbind(out,cbind(setNames(as.data.frame(rep(left[i,],nrow(right))),names(left)),right))}}
      }
      else{
        right=data.frame('upper'='alpha_all','lower'=Bs[1],'response'=sims[[Bs[[1]]]])
        if(length(left)>1){out=rbind(out,cbind(repframe(left[i,],nrow(right)),right))}
        else{out=rbind(out,cbind(setNames(as.data.frame(rep(left[i,],nrow(right))),names(left)),right))}
        medians[['all']]=right$response[right$upper=='alpha_all']}}
    progress(match(a,alphas)/length(alphas),50)}
  out=na.omit(out)#[c(T,rep(F,99)),] #culling
  out$response[out$lower %in% A]=out$response[out$lower %in% A]*as.numeric(hbm@scales[1])
  cat('\nbacktransforming data\n')
  for(sb in names(hbm@scales[-1])){
    out$response[regexpr(sb,out$lower)>0]=out$response[regexpr(sb,out$lower)>0]*as.numeric(hbm@scales[1])/as.numeric(hbm@scales[sb])
    progress(match(sb,names(hbm@scales[-1]))/length(names(hbm@scales[-1])),50)}
  hbm@output=out
  count=0
  cat('\nfixing hierarchy names\n')
  for(i in unique(names(filter))){
    group=unique(filter[,regexpr(i,unique(names(filter)))>0])
    group_filter=setNames(c('all',group[,2]),c('all',group[,1]))
    hbm@output[,i]=sapply(hbm@output[,i],\(x)group_filter[x])
    count=count+1
    progress(count/length(unique(names(filter))),50)}
  return(hbm)})
#################################################################################
setGeneric("run_model", function(hbm) standardGeneric("run_model"))
setMethod("run_model","hbm_object",function(hbm){
  hbm@jags_model=jagsUI::jags(data=hbm@model_data,
                              n.adapt=hbm@n.adapt,
                              n.burnin=hbm@n.burnin,
                              n.iter=hbm@n.iter,
                              n.chains=hbm@n.chains,
                              modules="glm",
                              model.file=paste(hbm@model_dir,'/',
                                               hbm@model_name,
                                               '.txt',sep=''),
                              parameters.to.save=hbm@save,
                              verbose=T,
                              DIC=F)
  return(hbm)})
#################################################################################
setGeneric("write_model", function(hbm) standardGeneric("write_model"))
setMethod("write_model","hbm_object",function(hbm){
  slot(hbm,'variables')=c(as.character(hbm@input$model)[2],
                          strsplit(gsub('\\((.+)\\)','\\1',
                                        as.character(hbm@input$model)[3]),
                                   ' \\+ ')[[1]])
  hbm@vars=strsplit(as.character(hbm@input$model)[3],'\\(')[[1]]|>
    (\(x)gsub(')','',x))()|>
    (\(x)gsub(' \\+ ',',',x))()|>
    (\(x)gsub(',$','',x))()|>
    strsplit(',')|>
    setNames(c('Variables','Groups'))

  var1_terms=c()
  hbm@vars[['Interaction terms']]=list()
  for(q in hbm@vars[[1]]){
    if(regexpr(':',q)>0){
      split=strsplit(q,':')[[1]]
      var1_terms=c(var1_terms,split[1])
      if(split[1]%in%names(hbm@vars[[3]])){
        hbm@vars[[3]][[split[1]]]=c(hbm@vars[[3]][[split[1]]],split[-1])}
      else{hbm@vars[[3]]=c(hbm@vars[[3]],setNames(list(split[-1]),split[1]))}}
    else{var1_terms=c(var1_terms,q)}}
  hbm@vars[[1]]=var1_terms
  
  new_i='i'
  while(T){                                                                     #if i in variable names generate replacement variable for iterations in model
    if(new_i%in%hbm@variables){new_i=paste(new_i,sample(letters,1),sep='')}
    else{break}}
  
  new=c(str_interp('${as.character(hbm@input$model)[2]}[i]~${response(hbm)}'),
        str_interp('${as.character(hbm@input$model)[2]}.sim[i]~${response(hbm)}'),
        str_interp('${as.character(hbm@input$model)[2]}.res[i]=${as.character(hbm@input$model)[2]}[i]-mu[i]'),
        str_interp('${as.character(hbm@input$model)[2]}.sres[i]=${as.character(hbm@input$model)[2]}.sim[i]-mu[i]'))
  left=paste(new,collapse='\n\t')
  
  for(q in hbm@vars[[2]]){
    if(q!=''){hbm@data[q]=factor(hbm@data[,q],levels=unique(hbm@data[,q]))}}
  
  hbm@data=hbm@data[order(hbm@data[,hbm@vars[[2]][1]]),]                        
  rv=is=c();c=0
  hbm@vars[[2]]=tryCatch(c('',hbm@vars[[2]]),error=function(e){c('')})
  x=paste(c('',rep('[',length(hbm@vars[[2]])-1)),hbm@vars[[2]],sep='')|>
    (\(x)paste(x,'[i]',sep=''))()|>
    (\(x)paste(x,c('',rep(']',length(hbm@vars[[2]])-1)),sep=''))()
  
  for(q in hbm@vars[[2]]){
    c=c+1;is=c(is,paste(c('',x[-1])[1:c],collapse=''))
    rv=c(rv,paste('alpha',q,sep='_'))}
  
  for(w in hbm@vars[[1]]){
    for(q in hbm@vars[[2]]){rv=c(rv,paste('beta',w,q,sep='_'))}}
  
  rv=chunk(gsub('_$','[i]',paste(rv,is,sep='')),length(hbm@vars[[2]]))
  right=c()
  rv=unique(rv)
  c=0
  
  for(q in rv[-1]){
    term=gsub('beta_','',gsub('\\[.+\\]','',q))[1]
    if(term%in%names(hbm@vars[[3]])){
      ints=c()
      for(item in hbm@vars[[3]][[term]]){
        int_q=gsub('beta_',paste('beta_int_',item,sep=''),q)
        rv=c(rv,list(int_q))
        ints=c(ints,str_interp("(${paste(int_q,collapse='+')})*${item}[i][i]"))}
      c=c+1;right=c(right,str_interp("(${paste(paste(q,collapse='+'),paste(ints,collapse='+'),sep='+')})*${hbm@vars[[1]][c]}[i]"))}
    else{c=c+1;right=c(right,str_interp("(${paste(q,collapse='+')})*${hbm@vars[[1]][c]}[i]"))}}
  
  rv=unique(rv)
  right=gsub('\\[i\\]\\)',')',gsub('\\[i\\]\\+','+',right))
  rn=gsub('\\[i\\]$','',rv[[1]])
  rn_ints=list()
  simple=rn
  
  for(term in hbm@vars[[3]]){
    for(item in unique(term)){
      int=gsub('^alpha',paste('Alpha_int_',item,sep=''),simple)
      rn_ints=c(rn_ints,list(int))
      int=paste('(',paste(int,collapse='+'),')*',item,'[i]',sep='')
      rn=c(rn,paste(int))}}
  
  top=str_interp("for(i in 1:N){\n\t${left}\n\tmu[i]=${paste(rn,collapse='+')}+\n\t${paste(right,collapse='+\n\t')}}")
  s=t=c()
  rv=c(rv,rn_ints)
  
  for(q in hbm@vars[[2]]){
    snam='sigma';tnam='tau'
    s=c(s,paste(snam,q,sep='_'))
    t=c(t,paste(tnam,q,sep='_'))}
  
  s=gsub('_$','',s)
  t=gsub('_$','',t)
  bottom=paste(s,'~dunif(0,100)\n',t,'=1/(',s,'*',s,')',sep='')|>
    (\(x)paste(x,collapse='\n'))()
  if(q!=''){
    bottom=paste('a',s,'~dunif(0,100)\n','a',t,'=1/(',s,'*',s,')',sep='')|>
      (\(x)paste(x,collapse='\n'))()|>
      (\(x)paste(bottom,x,sep='\n'))()}
  bottom=str_interp("${bottom}\n${paste(paste(gsub('\\\\[.+\\\\]','',as.data.frame(rv)[1,]),'~dnorm(0,0.01)',sep=''),collapse='\n')}")
  middle=c()
  
  if(length(hbm@vars[[2]])>1){
    for(q in 1:(length(hbm@vars[[2]])-1)){
      loop=str_interp('for(${letters[q]} in 1:N${hbm@vars[[2]][-1][q]}${c("",paste("[",letters,"]",sep=""))[q]}){')
      types=c()
      for(r in rv){
        if(regexpr('^[aA]lpha',r[1])>0){
          types=c(types,str_interp('~dnorm(0,a${t[q+1]})'))}
        else{types=c(types,str_interp('~dnorm(0,${t[q+1]})'))}}
      
      lines=gsub(is[q+1],paste(paste('[',letters[1:q],']',sep=''),collapse=''),
                 paste(
                   paste(as.data.frame(rv)[q+1,],types,sep=''),
                   collapse='\n\t'),
                 fixed=T)
      middle=c(middle,paste(loop,lines,sep='\n\t'))}}
    
  middle=paste(middle,collapse='\n')|>
    (\(x)paste(c(x,rep('}',length(hbm@vars[[2]])-1)),collapse=''))()
  tail=paste(c(str_interp('${as.character(hbm@input$model)[2]}.fit=sum(pow(${as.character(hbm@input$model)[2]}.res[],2)/tau)'),
               str_interp('${as.character(hbm@input$model)[2]}.sfit=sum(pow(${as.character(hbm@input$model)[2]}.sres[],2)/tau)')),collapse='\n')
  hbm@model=str_interp("model{\n${paste(top,middle,bottom,tail,sep='\n')}\n}")|>
    (\(x)gsub('\\]\\[',',',x))()|>
    (\(x)gsub('\\[i\\]',paste(c('\\[',new_i,'\\]'),collapse=''),x))()|>
    (\(x)gsub('for\\(i',paste(c('for\\(',new_i),collapse=''),x))()
  if(hbm@source_model==''){
    writeLines(hbm@model,paste(hbm@model_dir,'/',hbm@model_name,'.txt',sep=''))}
  
  hbm@save=c(str_interp('${as.character(hbm@input$model)[2]}.sfit'),
            str_interp('${as.character(hbm@input$model)[2]}.fit'))
  for(i in rv){hbm@save=c(hbm@save,i)}
  hbm@save=gsub('\\[.+\\]','',hbm@save)
  return(hbm)})
#################################################################################
hbm=function(data,model,...){
  output=new('hbm_object')
  param=list(...)
  defaults=list('name'        ='unnamed_hbm',
                'model_name'  ='unnamed_model',
                'dir'         ='.',
                'dist'        ='dnorm',
                'model_dir'   ='.',
                'n.adapt'     =2000,
                'n.burnin'    =1000,
                'n.iter'      =10000,
                'n.chains'    =4,
                'source_model'='')
  
  default=c(list('data' =data,
                 'input'=match.call(expand.dots=T)),
            lapply(names(defaults),\(x)ifelse(is.null(param[[x]]),
                                                       defaults[[x]],
                                                       param[[x]]))|>
            setNames(names(defaults)))
  for(d in names(default)){
    if(class(default[[d]])==class(slot(output,d))){slot(output,d)=default[[d]]}
    else{
      cat('Invalid input for variable "',d,'"\n',sep='')
      return()}}
  
  if(!(output@dist%in%c('dnorm'))){
    cat('Distribution "',output@dist,'" not yet implemented',sep='')
    return()}
  
  if(output@source_model==''){
    output=tryCatch(write_model(output),error=\(x){print(x);return(output)})}
  output=tryCatch(format_data(output),error=\(x){print(x);return(output)})
  output=tryCatch(run_model(output),error=\(x){print(x);return(output)})
  cat('\nFormatting data\n')
  output=tryCatch(format_model(output),error=\(x){print(x);return(output)})
  cat('\nmaking trace plots\n')
  traces(output,cull=9)
  
  write.csv(summary(output),
            paste0(output@dir,'/',output@name,'.csv'),row.names=F)
  for(n in names(output@jags_model$Rhat)){
    if(max(output@jags_model$Rhat[[n]],na.rm=T)>1.1){
      cat(rep('#',50),'\nRhat Greater than 1.1 in\n',
          paste(output@input,collapse = ' '),'\n',rep('#',50),'\n',sep='')
      break}}
  return(output)}
#################################################################################
wave=function(data,var,s){                                                      #create single distribution plots, takes argument of model output, focal variables, scale
  g=ggplot()+
    geom_vline(xintercept = 0,size=2*s,linetype='dashed')+
    geom_density(data=data[data[,var$var[1]]=='all',],
                 aes_string(x='response',group=var$var[1]),size=1.5*s,
                 fill='black',color='black',alpha=.5)+
    geom_density(data=data[data[,var$var[1]]!='all',],
                 aes_string(x='response',color=var$var[1]),size=1.5*s)
  return(list(g))}

wave2=function(data,var,s){                                                     #create multi-distribution plots, takes argument of model output, focal variables, scale
  data=subset(data,!(data[,var$var[1]]=='all'&data[,var$var[2]]=='all'))
  out=ggplot()+
      geom_vline(xintercept = 0,size=2*s,linetype='dashed')+
      geom_density(data=data[data[,var$var[1]]=='all',],
                   aes_string(x='response',group=var$var[2]),size=1.5*s,
                   fill='black',color='black',alpha=.5)+
      geom_density(data=data[data[,var$var[1]]!='all',],
                   aes_string(x='response',color=var$var[1]),size=1.5*s)+
      facet_grid(reformulate(var$var[2]))
  return(list(out))}

wave3=function(data,var,s,int){                                                     #create interaction + distribution plots, takes argument of model output, focal variables, scale
  data=subset(data,regexpr(int,lower)>0)
  data$interact=as.factor(regexpr(int,data$lower)>0)
  out=ggplot()+
    geom_vline(xintercept = 0,size=2*s,linetype='dashed')+
    geom_density(data=data[data[,var$var[1]]=='all',],
                 aes_string(x='response',y='..scaled..',group=var$var[1]),size=1.5*s,
                 fill='black',color='black',alpha=.5)+
    geom_density(data=data[data[,var$var[1]]!='all',],
                 aes_string(x='response',y='..scaled..',color=var$var[1]),size=1.5*s)+
    xlab('Difference in slopes')
  return(list(out))}

setGeneric("ocean", function(hbm,...) standardGeneric("ocean"))
setMethod("ocean","hbm_object",function(hbm,vars,fill='lower',s=1,interaction='none'){             #create groups of waveplots, takes arguments of model output, focal groups, fill groups, scale
  data=hbm@output
  out=list()
  for(g in vars){
    h=data[regexpr(g,data$lower)>0,]
    v=hbmvar(h,c(fill))
    if(interaction!='none'){out[[g]]=wave3(h,v,s,interaction)}
    else{if(length(v$var)==1){out[[g]]=wave(h,v,s)}
      else{if(v$var[2]=='upper'){out[[g]]=wave(h,v,s)}
        else{out[[g]]=wave2(h,v,s)}}}}
  return(unlist(out, recursive = FALSE))})
#################################################################################
dotplot=function(data,var,s){                                                   #create single dotplot, takes argument of model output, focal variables, scale
  df=setNames(as.data.frame(matrix(ncol=4)),c('variable','y2.5','y97.5','y50'))
  for(i in unique(data[,var$var[1]])){
    q=data$response[data[,var$var[1]]==i]
    df=rbind(df,data.frame('variable'=i,'y2.5'=quantile(q,0.05),
                           'y97.5'=quantile(q,0.95),'y50'=quantile(q,.5)))}
  g=ggplot(data)+
    geom_hline(yintercept=0,size=2*s)+
    geom_boxplot(data=na.omit(df),width=.05*s,lwd=s,
                 aes(x=as.factor(variable),ymin = y2.5, lower = y2.5, 
                     middle = y50, upper = y97.5, ymax = y97.5),
                 stat = "identity",fill='black')+
    geom_point(data=na.omit(df),aes(x=as.factor(variable),y=y50),size=5*s)
  
  return(g)}

dotplot2=function(data,var,s){                                                  #create multi-dotplots (new version), takes argument of model output, focal variables, scale
  df=setNames(as.data.frame(matrix(ncol=5)),
              c(var$var[2],'v2','y2.5','y97.5','y50'))
  for(i in unique(data[,var$var[2]])){
    for(j in unique(data[,var$var[1]])){
      q=data$response[data[,var$var[2]]==i&data[,var$var[1]]==j]
      df=rbind(df,setNames(data.frame('v'=i,'v2'=j,
                                      'y2.5'=quantile(q,0.05),
                                      'y97.5'=quantile(q,0.95),
                                      'y50'=quantile(q,.5)),
                           c(var$var[2],'v2','y2.5','y97.5','y50')))}}
  g=ggplot(data)+
    geom_hline(yintercept=0,size=2*s)+
    facet_grid(reformulate(var$var[2]))+
    geom_boxplot(data=na.omit(df),
                 width=.05*s,lwd=s,position = position_dodge(.9),
                 aes(x=as.factor(v2),group=as.factor(v2),
                     ymin = y2.5, lower = y2.5, middle = y50, 
                     upper = y97.5, ymax = y97.5),
                 stat = "identity",fill='black')+
    geom_point(data=na.omit(df),aes(x=as.factor(v2),
                                    group=as.factor(v2),y=y50),size=5*s)
  return(g)}

dotplot3=function(data,var,s,int){                                                   #create single dotplot, takes argument of model output, focal variables, scale
  data=subset(data,regexpr(int,lower)>0)
  data$interact=as.factor(regexpr(int,data$lower)>0)
  df=setNames(as.data.frame(matrix(ncol=4)),c('variable','y2.5','y97.5','y50'))
  for(i in unique(data[,var$var[1]])){
    q=data$response[data[,var$var[1]]==i]
    df=rbind(df,data.frame('variable'=i,'y2.5'=quantile(q,0.05),
                           'y97.5'=quantile(q,0.95),'y50'=quantile(q,.5)))}
  g=ggplot(data)+
    geom_hline(yintercept=0,size=2*s)+
    geom_boxplot(data=na.omit(df),width=.05*s,lwd=s,
                 aes(x=as.factor(variable),ymin = y2.5, lower = y2.5, 
                     middle = y50, upper = y97.5, ymax = y97.5),
                 stat = "identity",fill='black')+
    geom_point(data=na.omit(df),aes(x=as.factor(variable),y=y50),size=5*s)
  
  return(g)}

setGeneric("polka", function(hbm,...) standardGeneric("polka"))
setMethod("polka","hbm_object",function(hbm,vars,fill='lower',s=1,interaction='none'){             #create groups of dotplots, takes arguments of model output, focal groups, fill group, scale
  data=hbm@output
  out=list()
  for(g in vars){
    h=data[regexpr(g,data$lower)>0,]
    v=hbmvar(h,c(fill))
    if(interaction!='none'){out[[g]]=dotplot3(h,v,s,interaction)}
    else{if(length(v$var)==1){out[[g]]=dotplot(h,v,s)}
      else{if(v$var[2]=='upper'){out[[g]]=dotplot(h,v,s)}
        else{out[[g]]=dotplot2(h,v,s)}}}}
  return(out)})
#################################################################################
setGeneric("hbmgroup", function(hbm,...) standardGeneric("hbmgroup"))
setMethod("hbmgroup","hbm_object",function(hbm,groups){                         #Selects specific subsets of data based on focal group, takes arguments of model output, focal groups
  v=hbm@vars[[2]]
  out=setNames(as.data.frame(matrix(ncol=length(hbm@output))),names(hbm@output))
  for(g in groups){
    if(!(g %in% v)){out=rbind(out,hbm@output[regexpr(g,hbm@output$lower)>0,])}
    else{out=rbind(out,hbm@output[regexpr('all',hbm@output[,g])>0,])}}
  for(g in groups){
    if(!(g %in% v)){out=out[regexpr(g,out$lower)>0,]}
    else{out=out[regexpr('all',out[,g])>0,]}}
  return(na.omit(out))})

hbmvar=function(data,groups){                                                   #selects specific variable names based on focal groups, takes arguments of model output, focal groups
  n=names(data)
  filter=c()
  for(g in groups){
    if(g %in% n){
      if(match(g,n)==1){Var=c(g)}
      else{Var=c(g,n[match(g,n)-1])}}
    else{filter=c(filter,g)}}
  return(list('var'=Var,'filter'=filter))}

mgroup=function(data,groups){                                                   #creates data groupings based on focal group names, takes arguments of model output, focal groups
  k=0
  for(g in groups){
    if(k==0){h=hbmgroup(data,c(g));k=1}
    else{h=rbind(h,hbmgroup(data,c(g)))}}
  return(h)}

#################################################################################
cello=function(data,var,s,label='none',lsize=1){                                #create single violin plot with CI bars, takes argument of model output, focal variables, scale
  df=setNames(as.data.frame(matrix(ncol=5)),c('variable','y2.5','y97.5','y50','lab'))
  for(i in unique(data[,var$var[1]])){
    q=data$response[data[,var$var[1]]==i]
    df=rbind(df,data.frame('variable'=i,'y2.5'=quantile(q,0.05),
                           'y97.5'=quantile(q,0.95),'y50'=quantile(q,.5),
                           'lab'=paste0(sign(median(q))*round(100*pd(q),1),'%')))}
  g=ggplot(data)+
    geom_hline(yintercept=0,size=2*s)+
    geom_violin(aes_string(x=var$var[1],y='response',
                           fill=var$var[1]),size=s,scale='width')+
    geom_boxplot(data=na.omit(df),width=.05*s,lwd=s,
                 aes(x=as.factor(variable),
                     ymin = y2.5, lower = y2.5, middle = y50, 
                     upper = y97.5, ymax = y97.5),stat = "identity")
  if(label!='none'){g=g+geom_text(data=na.omit(df),aes(label=lab,
                                                       x=as.factor(variable),
                                                       y=label),size=lsize)}
  
  return(g)}

cello2=function(data,var,s){                                                    #create multi-violin plots with CI bars, takes argument of model output, focal variables, scale
  df=setNames(as.data.frame(matrix(ncol=5)),
              c(var$var[2],'v2','y2.5','y97.5','y50'))
  for(i in unique(data[,var$var[2]])){
    for(j in unique(data[,var$var[1]])){
      q=data$response[data[,var$var[2]]==i&data[,var$var[1]]==j]
      df=rbind(df,setNames(data.frame('v'=i,'v2'=j,'y2.5'=quantile(q,0.05),
                                      'y97.5'=quantile(q,0.95),
                                      'y50'=quantile(q,.5)),
                           c(var$var[2],'v2','y2.5','y97.5','y50')))}}
  g=ggplot(data)+
    geom_hline(yintercept=0,size=2*s)+
    facet_grid(reformulate(var$var[2]))+
    geom_violin(aes_string(x=var$var[1],y='response',
                           fill=var$var[1]),size=s,scale='width')+
    geom_boxplot(data=na.omit(df),width=.05*s,lwd=s,position=position_dodge(.9),
                 aes(x=as.factor(v2),group=as.factor(v2),
                     ymin = y2.5, lower = y2.5, middle = y50, 
                     upper = y97.5, ymax = y97.5),stat = "identity")
  return(g)}

cello3=function(data,var,s,label='none',lsize=1,int){                                #create single violin plot with CI bars, takes argument of model output, focal variables, scale
  data=subset(data,regexpr(int,lower)>0)
  data$interact=as.factor(regexpr(int,data$lower)>0)
  df=setNames(as.data.frame(matrix(ncol=5)),c('variable','y2.5','y97.5','y50','lab'))
  for(i in unique(data[,var$var[1]])){
    q=data$response[data[,var$var[1]]==i]
    df=rbind(df,data.frame('variable'=i,'y2.5'=quantile(q,0.05),
                           'y97.5'=quantile(q,0.95),'y50'=quantile(q,.5),
                           'lab'=paste0(sign(median(q))*round(100*pd(q),1),'%')))}
  g=ggplot(data)+
    geom_hline(yintercept=0,size=2*s)+
    geom_violin(aes_string(x=var$var[1],y='response',
                           fill=var$var[1]),size=s,scale='width')+
    geom_boxplot(data=na.omit(df),width=.05*s,lwd=s,
                 aes(x=as.factor(variable),
                     ymin = y2.5, lower = y2.5, middle = y50, 
                     upper = y97.5, ymax = y97.5),stat = "identity")
  if(label!='none'){g=g+geom_text(data=na.omit(df),aes(label=lab,
                                                       x=as.factor(variable),
                                                       y=label),size=lsize)}
  
  return(g)}

setGeneric("bass", function(hbm,...) standardGeneric("bass"))
setMethod("bass","hbm_object",function(hbm,groups,fill='lower',s=1,label='none',lsize=1,interaction='none'){               #create groups cello plots, takes arguments of model output, focal groups, fill group, scale
  data=hbm@output
  out=list()
  for(g in groups){
    h=data[regexpr(g,data$lower)>0,]
    v=hbmvar(h,c(fill))
    if(interaction!='none'){out[[g]]=cello3(h,v,s,label=label,lsize=lsize,int=interaction)}
    else{if(length(v$var)==1){out[[g]]=cello(h,v,s,label=label,lsize=lsize)}
      else{if(v$var[2]=='upper'){out[[g]]=cello(h,v,s,label=label,lsize=lsize)}
        else{out[[g]]=cello2(h,v,s)}}}}
  return(out)})
#################################################################################
clamp=function(x,minimum,maximum){return(ifelse(x<minimum,minimum,ifelse(x>maximum,maximum,x)))}

randseq=function(x,size){out=c();for(i in 1:size){out=c(out,sample(x,1))};return(out)}

fixedmean=function(n,m,d){
  out=c(rnorm(n-1,m,d))
  return(c(out,(m-sum(out)/n)*n))}

# x => list(distrubition, mean, sd or other parameter)
# y => list(x,alpha,slope)
# g => list(groups,sd)

mixed=function(n,params){
  data=setNames(as.data.frame(matrix(ncol=length(names(params)),nrow=n)),names(params))
  rslope=list();ralpha=list()
  resp=names(params)[sapply(names(params),\(x)length(params[[x]])>=3&is.character(params[[x]][[1]]))]
  for(p in setdiff(names(params),resp)){
    if(length(params[[p]])==2){
      data[p]=randseq(params[[p]][[1]],n)
      rslope[[p]]=setNames(fixedmean(length(params[[p]][[1]]),params[[resp[1]]][[3]],params[[p]][[2]]),params[[p]][[1]])
      ralpha[[p]]=setNames(fixedmean(length(params[[p]][[1]]),params[[resp[1]]][[2]],params[[p]][[2]]),params[[p]][[1]])}}
  for(p in setdiff(names(params),names(rslope))){
    if(length(params[[p]])>=3){
      if(is.function(params[[p]][[1]])){data[p]=params[[p]][[1]](n,params[[p]][[2]],params[[p]][[3]])}
      else{
        slope=rep(1/(mean(rslope[[1]])^(length(rslope)-1)),n)
        alpha=rep(0,n)
        for(i in names(rslope)){
          slope=slope*rslope[[i]][data[,i]]
          alpha=alpha+ralpha[[i]][data[,i]]/length(names(ralpha))}
        data[p]=data[params[[p]][[1]]]*slope+alpha+rnorm(n,0,params[[p]][[2]])}}}
  return(list(data,ralpha,rslope))}

####examples#####################################################################
#generate data
data=mixed(200,list('mass'=list(runif,10,20),
                    'length'=list('mass',1,2),
                    'site'=list(LETTERS[1:5],1),
                    'species'=list(c('s1','s2','s3'),1),
                    'sex'=list(c('M','F'),1)))[[1]]
  
ggplot(data,aes(x=mass,y=abs(length),color=as.factor(sex)))+
  geom_point()+
  geom_smooth(method=lm,alpha=0)+
  facet_grid(cols=vars(site),rows=vars(species))+
  theme_classic()
#single level
o=hbm(rbind(data,data),length~mass+(site))
summary(o) 
ocean(o,'mass','site')[[1]]
polka(o,'mass','site')[[1]]
bass(o,'mass','site')[[1]]

#multiple levels
o=hbm(data,length~mass+(site+species))
summary(o) 
ocean(o,'mass','species')[[1]]
polka(o,'mass','species')[[1]]
bass(o,'mass','species')[[1]]

#interaction
o=hbm(data,length~mass:sex+(site))
summary(o) 
ocean(o,'mass','site',interaction='sex')[[1]]
polka(o,'mass','site',interaction='sex')[[1]]
bass(o,'mass','site',interaction='sex')[[1]]


o=hbm(data,length~mass+(site+species))
summary(o)




# 
# data$sex=as.numeric(data$sex)
# model='mass~length+sex
# length~sex'
# b=semb(data,model)
