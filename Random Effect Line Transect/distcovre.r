fitds=function(obs,width,detfct="hn",scale.formula=~1,exponent.formula=~1)
{
  # create scale design matrix with formula and data
  scale_mat=model.matrix(scale.formula,obs)
  if(detfct=="haz")exponent_mat=model.matrix(exponent.formula,obs)
  # write out data file
  con=file(paste(tplfile,".dat",sep=""),open="wt")
  writeLines(as.character(nrow(obs)),con)
  write(obs$distance,con,ncolumns=1)
  writeLines(as.character(width),con)
  if(detfct=="haz")
  {
     writeLines("2",con)
     writeLines("2",con)
     writeLines("2  2",con)
     writeLines(paste(ncol(scale_mat)," ",ncol(exponent_mat),sep=""),con)
     write(t(scale_mat),con,ncolumns=ncol(scale_mat))
     write(t(exponent_mat),con,ncolumns=ncol(exponent_mat))
  }else
  {
     writeLines("1",con)
     writeLines("1",con)
     writeLines("2",con)
     writeLines(paste(ncol(scale_mat),sep=""),con)
     write(t(scale_mat),con,ncolumns=ncol(scale_mat))
   }   
  close(con)
  run_admb(tplfile)
  results=read_admb(tplfile)
  cnames=paste("scale:",colnames(scale_mat),sep="")
  if(detfct=="haz")
    cnames=c(cnames,paste("exponent:",colnames(exponent_mat),sep=""))
  names(results$coefficients)=cnames
  rownames(results$vcov)=cnames
  colnames(results$vcov)=cnames
  return(results)
}

tplfile="distcovre"
# compile tpl file
prepare_admb()
compile_admb(tplfile,re=TRUE,verbose=TRUE)
# simulate some data 
n=100
x=abs(rnorm(n,0,rnorm(n,4,1)))
obs=data.frame(distance=x,object=1:n,seen=rep(1,n))
# fit for different models
detfct=fitds(obs,10,"hn",~1)



