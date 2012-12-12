DATA_SECTION
   init_int n;                        // number of distances
   init_vector xs(1,n);               // distances
   init_number width;                 // truncation half-width of transect
   init_int ifct;                     // type of detection function 1=hn, 2=haz
   init_int pt;                       // number of parameter types
   init_ivector links(1,pt);          // link number but would rather use strings; not sure how 1=identity, 2=log, 3=logit
   init_ivector k(1,pt);              // vector of number of parameters for each type; cols in design matrix
   init_3darray dm(1,pt,1,n,1,k);     // design matrices - one for each parameter type
   int obs_index;                     // index to observation
   int par_index;                     // index to parameter
INITIALIZATION_SECTION
   sigma 0.2
   beta 1.38
PARAMETER_SECTION
   init_matrix beta(1,pt,1,k,2);        // beta parameters for each parameter type
   init_bounded_vector sigma(1,pt,0.000001,5,1);
   matrix parmat(1,pt,1,n);           // matrix of parameter values; 1 to n and 1 to pt types of parameters
   random_effects_matrix u(1,pt,1,n);
   objective_function_value f;        // negative log-likelihood
PROCEDURE_SECTION
   int i, j;
   dvariable mu;
   dvariable g;
// funnel_dvariable mu;
// Create matrix of real parameter values which is pt rows and n columns
   for (i=1;i<=pt;i++)
   {
    par_index=i;
    parmat(i)=reals(dm(i),beta(i),links(i));
   }
   cout << "parmat = " << parmat << endl;
// loop over each observation computing sum of log-likelihood values
   f=0;
   g=0;
   for (j=1;j<=n;j++)
   {
      obs_index=j;
    //  mu=adromb(0,width,8);
	  f -= -0.5*square(u(1,j))-log(sqrt(2*PI));                         // log of std normal density for epsilon
      f -= -log(sqrt(2*PI))-log(parmat(1,j))-0.5*square(xs(j)/parmat(1,j));        // log of f(x) for half-normal
	//  f-= log(h(xs(j))) - log(mu);
	//  f-= -0.5*square(u(1,j))-0.9189385332046727;
   }  
   cout << "beta = " << beta << endl;
   cout << "sigma = " << sigma << endl;   
   cout << "f = " << f << endl;   
//////////////////////////////   
// Computes reals from betas
//////////////////////////////   
FUNCTION dvar_vector reals(dmatrix& dm, dvar_vector& beta, int ilink)
// dm is the design matrix
// beta is vector of parameters - length matches ncol(dm)
// ilink is type of link function
// u is random values
    dvar_vector tmp;
	if(ilink==1)
        tmp=dm*beta+u(par_index)*sigma(par_index);
    if(ilink==2)
        tmp=exp(dm*beta+u(par_index)*sigma(par_index));
    if(ilink==3)
        tmp=1/(1+exp(-dm*beta-u(par_index)*sigma(par_index)));
    return tmp;
//////////////////////////////////////////////////////////   	
// Computes normalizing constant int fct(x) - 0 to width
//////////////////////////////////////////////////////////   	
FUNCTION dvariable h(const double x)
// x is integration variable
// ifct is index for function read from data
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

FUNCTION dvariable adromb(double a,double b,int ns)
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

FUNCTION dvariable trapzd(double a,double b,int n,const dvariable& s)
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
TOP_OF_MAIN_SECTION
  gradient_structure::set_MAX_NVAR_OFFSET(1001001); 
  gradient_structure::set_NUM_DEPENDENT_VARIABLES(800);
  gradient_structure::set_GRADSTACK_BUFFER_SIZE(900000);
  gradient_structure::set_CMPDIF_BUFFER_SIZE(1000000);
  arrmblsize=5000000;
