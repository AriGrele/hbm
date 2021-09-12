source('hbm.R')
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
mix=mixed(200,list('mass'=list(runif,10,20),
                   'length'=list('mass',1,2),
                   'site'=list(LETTERS[1:5],1),
                   'species'=list(c('s1','s2','s3'),1),
                   'sex'=list(c('M','F'),1),
                   'color'=list(c('red','blue','yellow'),1)))[[1]]

data=data.frame(mass=mix$mass)
data$slope=1:4
data$site=letters[1:4]
data$length=rnorm(200,data$mass*data$slope,1)
data$offspring=rgamma(200,shape=1,rate=1/data$mass*data$slope)
data$survival=(data$mass>15)+0
cowplot::plot_grid(
  ggplot(data,aes(x=mass,y=length,color=site))+
    geom_point()+
    geom_smooth(method=lm),
  ggplot(data,aes(x=mass,y=offspring,color=site))+
    geom_point()+
    geom_smooth(method=lm),
  ggplot(data,aes(x=mass,y=survival,color=site))+
    geom_point()+
    geom_smooth(method=lm),
  ggplot(data,aes(x=length,color=site))+geom_density(),
  ggplot(data,aes(x=offspring,color=site))+geom_density(),
  ggplot(data,aes(x=survival,color=site))+geom_density(bw=.01),ncol=3)
mix$species[mix$site=='A']=c('fish','moose')
mix$species[mix$site=='B']=c('s1','s4')
ggplot(mix,aes(x=mass,y=length,color=paste(species,color)))+
  geom_point()+
  geom_smooth(method=lm)+
  facet_grid(cols=vars(site))

set.seed(1)
source('hbm.R')

o1=hbm(mix,length~mass)

bass(o1,'mass','site')
resid(o1)
summary(o1)
fits(o1)

set.seed(1)
o2=hbm(data,offspring~mass+(site),dist='dgamma')
resid(o2)
fits(o2)
set.seed(1)

o3=hbm(data,survival~mass+(site),dist='dbern')
resid(o3)
o3fits(o3)
summary(o3)

mix$sex=as.numeric(as.factor(mix$sex))
model='mass~length+sex
length~sex'
b=bsem(data,model)
