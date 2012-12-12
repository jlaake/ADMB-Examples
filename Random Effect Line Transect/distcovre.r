fitds=function(obs,width,detfct="hn",scale.formula=~1,exponent.formula=~1,extra.args="-est -gh 10")
{
  # create scale design matrix with formula and data
  scale_mat=model.matrix(scale.formula,obs)
  if(detfct=="haz")exponent_mat=model.matrix(exponent.formula,obs)
  # write out data file
  con=file(paste(tplfile,".dat",sep=""),open="wt")
  writeLines(as.character(nrow(obs)),con)
  writeLines(as.character(width),con)
  write(obs$distance,con,ncolumns=1)
#  if(detfct=="haz")
#  {
#     writeLines("2",con)
#     writeLines("2",con)
#     writeLines("2  2",con)
#     writeLines(paste(ncol(scale_mat)," ",ncol(exponent_mat),sep=""),con)
#     write(t(scale_mat),con,ncolumns=ncol(scale_mat))
#     write(t(exponent_mat),con,ncolumns=ncol(exponent_mat))
#  }else
#  {
#     writeLines("1",con)
#     writeLines("1",con)
#     writeLines("2",con)
#     writeLines(paste(ncol(scale_mat),sep=""),con)
#     write(t(scale_mat),con,ncolumns=ncol(scale_mat))
#   }   
  close(con)
  run_admb(tplfile,verbose=TRUE,extra.args=extra.args)
  results=read_admb(tplfile)
#  cnames=paste("scale:",colnames(scale_mat),sep="")
#  if(detfct=="haz")
#    cnames=c(cnames,paste("exponent:",colnames(exponent_mat),sep=""))
#  names(results$coefficients)=cnames
#  rownames(results$vcov)=cnames
#  colnames(results$vcov)=cnames
  return(results)
}

tplfile="inthnre_f2"
# compile tpl file
library(RandomScale)
prepare_admb()
compile_admb(tplfile,re=TRUE,verbose=TRUE)
# simulate some data 
x=simdata(n=200,w=20,beta=1,beta_eps=-.5)
obs=data.frame(distance=x,object=1:n,seen=rep(1,n))
# fit for different models
detfct=fitds(obs,20,"hn",~1,extra.args="-est -gh 10")
detfct$loglik
flike=fitadmb(x,w=20,likelihood="f2")
flike$loglik
flike=fitadmb(x,w=20,likelihood="f1")
flike$loglik
glike=fitadmb(x,w=20,likelihood="g")
glike$loglik

