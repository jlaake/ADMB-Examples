# Number of histories
n=500
# Number of occasions
nocc=10
# Matrix of survival events - fixed survival
alive=matrix(rbinom(n*(nocc-1),1,.9),ncol=nocc-1)
# Cummulative survival events (once dead stays dead)
alive=t(apply(alive,1,cumprod))
# Individual variation in capture probability
p=plogis(rnorm(n,0,1))
# Capture events
seen=matrix(rbinom(n*(nocc-1),1,rep(p,each=nocc-1)),byrow=T,ncol=nocc-1)
# Construct capture history matrix for a single release cohort
chmat=cbind(rep(1,n),seen*alive)
ch=apply(chmat,1,paste,collapse="")
# Create dataframe with ch values
test=data.frame(ch=ch,stringsAsFactors=FALSE)
# attach marked package
library(marked)
# process data frame but don't accumulate same capture histories
test.proc=process.data(test,model="cjs",begin.time=1,accumulate=F)
# create design data
test.ddl=make.design.data(test.proc)
# Fit model
admb.mod=crm(test.proc,test.ddl,use.admb=TRUE,re=TRUE,extra.args="-gh 10")
# Show results
admb.mod
# Now fit the same model with MARK via RMark
# detach marked and attach RMark
detach("package:marked")
library(RMark)
# process same dataframe for RMark
test.proc=process.data(test,model="CJSRandom",begin.time=1)
# make design data (different format from marked)
test.ddl=make.design.data(test.proc)
# fit model
mark.mod=mark(test.proc,test.ddl)
# show results
mark.mod

