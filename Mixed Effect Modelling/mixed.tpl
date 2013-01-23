// Skeleton template for a mixed effects model structure; with logistic regression likelihood
// Jeff Laake; 22 Jan 2013
DATA_SECTION 
    init_int n;                           // number of rows in data
    init_vector y(1,n);                   // vector of responses
    init_int kfixed;                      // number of columns in design matrix(DM) for fixed effects
    init_matrix fixedDM(1,n,1,kfixed);    // fixed effect DM
    init_int nre;                         // number of random effects
    int phase;                            // if no random effects set phase to -1 otherwise 2            
	!! phase=2;
	!! if(nre==0)phase=-1;
    init_int krand;                       // number of columns in random effect DM
    init_matrix randDM(1,n,1,krand);      // random effect DM
    init_imatrix randIndex(1,n,1,krand);  // random effect indices for DM
        
PARAMETER_SECTION
    init_vector Beta(1,kfixed,1);         // parameter vector for fixed effects
    init_vector Tau(1,krand,phase);       // parameter vector for log(sigma)
    objective_function_value g;           // objective function - negative log-likelihood
    random_effects_vector u(1,nre,phase); // random effects vector
 
PROCEDURE_SECTION
    int i;                                                // index variable
    g=0;                                                  // initialize negative log-likelihood
	cout<<"nre = "<< nre <<endl;
	cout<<"kfixed = "<< kfixed <<endl;
	cout<<"krand = "<< krand <<endl;
    for (i=1;i<=nre;i++)                                  // if any random effects compute likelihood contribution for each
        n01_prior(u(i));                                  // u's are N(0,1) distributed
	cout<<"g = "<< g <<endl;

    for(i=1;i<=n;i++)                                     // loop over rows in data
          ll_i(i,Beta,Tau,u(randIndex(i)));               //compute negative log-likelihood for each row with separable function
	cout<<"g = "<< g <<endl;

SEPARABLE_FUNCTION void n01_prior(const prevariable&  u)  // taken from glmmadmb.tpl; uses PI 
 g -= -0.5*log(2.0*PI) - 0.5*square(u);
 
SEPARABLE_FUNCTION void ll_i(const int i, const dvar_vector& Beta,const dvar_vector& Tau,const dvar_vector& u)
    int j;
    dvariable mu;                                         // link value
    mu=fixedDM(i)*Beta;                                   // fixed portion 
    for(j=1;j<=krand;j++)                                 // random portion (if any)
        mu+=randDM(i,j)*u(j)*mfexp(Tau(j));
	
    g-= y(i)*mu - log(1+exp(mu));                         // Bernoulli likelihood with logit link
