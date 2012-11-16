#include <admodel.h>

#include <df1b2fun.h>

#include <adrndeff.h>

  extern "C"  {
    void ad_boundf(int i);
  }
#include <distcovre.htp>

  df1b2_parameters * df1b2_parameters::df1b2_parameters_ptr=0;
  model_parameters * model_parameters::model_parameters_ptr=0;
model_data::model_data(int argc,char * argv[]) : ad_comm(argc,argv)
{
  n.allocate("n");
  xs.allocate(1,n,"xs");
  width.allocate("width");
  ifct.allocate("ifct");
  pt.allocate("pt");
  links.allocate(1,pt,"links");
  k.allocate(1,pt,"k");
  dm.allocate(1,pt,1,n,1,k,"dm");
}

void model_parameters::initializationfunction(void)
{
  sigma.set_initial_value(0.2);
  beta.set_initial_value(0.6666667);
}

model_parameters::model_parameters(int sz,int argc,char * argv[]) : 
 model_data(argc,argv) , function_minimizer(sz)
{
  model_parameters_ptr=this;
  initializationfunction();
  beta.allocate(1,pt,1,k,2,"beta");
  sigma.allocate(1,pt,0.000001,5,1,"sigma");
  parmat.allocate(1,pt,1,n,"parmat");
  #ifndef NO_AD_INITIALIZE
    parmat.initialize();
  #endif
  u.allocate(1,pt,1,n,"u");
  prior_function_value.allocate("prior_function_value");
  likelihood_function_value.allocate("likelihood_function_value");
  f.allocate("f");  /* ADOBJECTIVEFUNCTION */
}
void model_parameters::userfunction(void)
{
  f =0.0;
   int i, j;
   dvariable mu;
   for (i=1;i<=pt;i++)
   {
    par_index=i;
    parmat(i)=reals(dm(i),beta(i),links(i));
   }
   f=0;
   for (j=1;j<=n;j++)
   {
      obs_index=j;
      mu=adromb(0,width,8);
	  f-= log(h(xs(j))) - log(mu);
	  f-= -0.5*square(u(1,j))-0.9189385332046727;
   }  
   cout << "beta = " << beta << endl;
   cout << "sigma = " << sigma << endl;   
   cout << "f = " << f << endl;   
}

dvar_vector model_parameters::reals(dmatrix& dm, dvar_vector& beta, int ilink)
{
    dvar_vector tmp;
	if(ilink==1)
        tmp=dm*beta+u(par_index)*sigma(par_index);
    if(ilink==2)
        tmp=exp(dm*beta+u(par_index)*sigma(par_index));
    if(ilink==3)
        tmp=1/(1+exp(-dm*beta-u(par_index)*sigma(par_index)));
    return tmp;
}

dvariable model_parameters::h(const double x)
{
   dvariable tmp;
   if(ifct==1)
     tmp=exp(-.5*square(x/parmat(1,obs_index)));
   if(ifct==2)
   {
     if(x<0.0000001)
	   tmp=1;
	 else
       tmp=1-exp(-(pow(x/parmat(1,obs_index),-(1+parmat(2,obs_index)))));
   }
   return tmp;
}

dvariable model_parameters::adromb(double a,double b,int ns)
{
  const double base = 4;
  const int JMAX=50;
  dvariable ss=0.0;
  int MAXN = min(JMAX, ns);
  dvar_vector s(1,MAXN+1);
  for(int j=1; j<=MAXN+1; j++)
  {
    s[j] = trapzd(a,b,j,ss);
  }
  for(int iter=1; iter<=MAXN+1; iter++)
  {
    for(int j=1; j<=MAXN+1-iter; j++)
    {
      s[j] = (pow(base,iter)*s[j+1]-s[j])/(pow(base,iter)-1);
    }
  }
  return s[1];
}

dvariable model_parameters::trapzd(double a,double b,int n,const dvariable& s)
{
  double x,num_interval,hn;
  dvariable sum;
 // static dvariable s;
  static int interval;
  int j;
  model_parameters * ptr= (model_parameters *) mycast();
  if (n == 1) {
    interval=1;
    return (s=0.5*(b-a)*(h(a)+h(b)));
  } else {
    num_interval=interval;
    hn=(b-a)/num_interval;
    x=a+0.5*hn;
    for (sum=0.0,j=1;j<=interval;j++,x+=hn) sum += h(x);
    interval *= 2;
    s=0.5*(s+(b-a)*sum/num_interval);
    return s;
  }
}
  long int arrmblsize=0;

int main(int argc,char * argv[])
{
  ad_set_new_handler();
  ad_exit=&ad_boundf;
  gradient_structure::set_MAX_NVAR_OFFSET(50502); 
  gradient_structure::set_NUM_DEPENDENT_VARIABLES(800);
  gradient_structure::set_GRADSTACK_BUFFER_SIZE(100000);
  gradient_structure::set_CMPDIF_BUFFER_SIZE(1000000);
  arrmblsize=500000;
    gradient_structure::set_NO_DERIVATIVES();
    gradient_structure::set_YES_SAVE_VARIABLES_VALUES();
  #if defined(__GNUDOS__) || defined(DOS386) || defined(__DPMI32__)  || \
     defined(__MSVC32__)
      if (!arrmblsize) arrmblsize=150000;
  #else
      if (!arrmblsize) arrmblsize=25000;
  #endif
    df1b2variable::noallocate=1;
    df1b2_parameters mp(arrmblsize,argc,argv);
    mp.iprint=10;

    function_minimizer::random_effects_flag=1;
    df1b2variable::noallocate=0;
    mp.preliminary_calculations();
    mp.computations(argc,argv);
    return 0;
}

extern "C"  {
  void ad_boundf(int i)
  {
    /* so we can stop here */
    exit(i);
  }
}

void model_parameters::preliminary_calculations(void){
  admaster_slave_variable_interface(*this);
  }

model_data::~model_data()
{}

model_parameters::~model_parameters()
{}

void model_parameters::report(void){}

void model_parameters::final_calcs(void){}

void model_parameters::set_runtime(void){}

#ifdef _BORLANDC_
  extern unsigned _stklen=10000U;
#endif


#ifdef __ZTC__
  extern unsigned int _stack=10000U;
#endif

void df1b2_parameters::user_function(void)
{
  f =0.0;
   int i, j;
   df1b2variable mu;
   for (i=1;i<=pt;i++)
   {
    par_index=i;
    parmat(i)=reals(dm(i),beta(i),links(i));
   }
   f=0;
   for (j=1;j<=n;j++)
   {
      obs_index=j;
      mu=adromb(0,width,8);
	  f-= log(h(xs(j))) - log(mu);
	  f-= -0.5*square(u(1,j))-0.9189385332046727;
   }  
   cout << "beta = " << beta << endl;
   cout << "sigma = " << sigma << endl;   
   cout << "f = " << f << endl;   
}

df1b2vector df1b2_parameters::reals(dmatrix& dm, df1b2vector& beta, int ilink)
{
    df1b2vector tmp;
	if(ilink==1)
        tmp=dm*beta+u(par_index)*sigma(par_index);
    if(ilink==2)
        tmp=exp(dm*beta+u(par_index)*sigma(par_index));
    if(ilink==3)
        tmp=1/(1+exp(-dm*beta-u(par_index)*sigma(par_index)));
    return tmp;
}

df1b2variable df1b2_parameters::h(const double x)
{
   df1b2variable tmp;
   if(ifct==1)
     tmp=exp(-.5*square(x/parmat(1,obs_index)));
   if(ifct==2)
   {
     if(x<0.0000001)
	   tmp=1;
	 else
       tmp=1-exp(-(pow(x/parmat(1,obs_index),-(1+parmat(2,obs_index)))));
   }
   return tmp;
}

df1b2variable df1b2_parameters::adromb(double a,double b,int ns)
{
  const double base = 4;
  const int JMAX=50;
  df1b2variable ss=0.0;
  int MAXN = min(JMAX, ns);
  df1b2vector s(1,MAXN+1);
  for(int j=1; j<=MAXN+1; j++)
  {
    s[j] = trapzd(a,b,j,ss);
  }
  for(int iter=1; iter<=MAXN+1; iter++)
  {
    for(int j=1; j<=MAXN+1-iter; j++)
    {
      s[j] = (pow(base,iter)*s[j+1]-s[j])/(pow(base,iter)-1);
    }
  }
  return s[1];
}

df1b2variable df1b2_parameters::trapzd(double a,double b,int n,const df1b2variable& s)
{
  double x,num_interval,hn;
  df1b2variable sum;
 // static df1b2variable s;
  static int interval;
  int j;
  model_parameters * ptr= (model_parameters *) mycast();
  if (n == 1) {
    interval=1;
    return (s=0.5*(b-a)*(h(a)+h(b)));
  } else {
    num_interval=interval;
    hn=(b-a)/num_interval;
    x=a+0.5*hn;
    for (sum=0.0,j=1;j<=interval;j++,x+=hn) sum += h(x);
    interval *= 2;
    s=0.5*(s+(b-a)*sum/num_interval);
    return s;
  }
}
   
void df1b2_pre_parameters::setup_quadprior_calcs(void) 
{ 
  df1b2_gradlist::set_no_derivatives(); 
  quadratic_prior::in_qp_calculations=1; 
}  
  
void df1b2_pre_parameters::begin_df1b2_funnel(void) 
{ 
  (*re_objective_function_value::pobjfun)=0; 
  other_separable_stuff_begin(); 
  f1b2gradlist->reset();  
  if (!quadratic_prior::in_qp_calculations) 
  { 
    df1b2_gradlist::set_yes_derivatives();  
  } 
  funnel_init_var::allocate_all();  
}  
 
void df1b2_pre_parameters::end_df1b2_funnel(void) 
{  
  lapprox->do_separable_stuff(); 
  other_separable_stuff_end(); 
} 
  
void model_parameters::begin_df1b2_funnel(void) 
{ 
  if (lapprox)  
  {  
    {  
      begin_funnel_stuff();  
    }  
  }  
}  
 
void model_parameters::end_df1b2_funnel(void) 
{  
  if (lapprox)  
  {  
    end_df1b2_funnel_stuff();  
  }  
} 

void df1b2_parameters::allocate(void) 
{
  beta.allocate(1,pt,1,k,2,"beta");
  sigma.allocate(1,pt,0.000001,5,1,"sigma");
  parmat.allocate(1,pt,1,n,"parmat");
  #ifndef NO_AD_INITIALIZE
    parmat.initialize();
  #endif
  u.allocate(1,pt,1,n,"u");
  prior_function_value.allocate("prior_function_value");
  likelihood_function_value.allocate("likelihood_function_value");
  f.allocate("f");  /* ADOBJECTIVEFUNCTION */
}
