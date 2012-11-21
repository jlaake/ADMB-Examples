DATA_SECTION
   init_int n;                        // number of distances
   init_number width;                 // truncation width; if -1 use infinity
   init_vector xs(1,n);               // distances
 PARAMETER_SECTION
   init_number beta;                  // beta parameter for log-sigma;
   init_bounded_number sigeps(0.000001,5);   
                                      // sigma for random effect;             
   number pi=3.14159265358979;
   random_effects_vector u(1,n);      // random effect for scale
   objective_function_value f;        // negative log-likelihood


PROCEDURE_SECTION
   int j;
// loop over each observation computing sum of log-likelihood values
   f=0;
   for (j=1;j<=n;j++)
   {
      ll_j(j,beta,sigeps,u(j));
   }  

SEPARABLE_FUNCTION void ll_j(const int j, const dvariable& beta,const dvariable& sigeps,const dvariable& u)
   dvariable eps=u*sigeps;
   dvariable sigma=exp(beta+eps);
   f -= -0.5*square(u)-log(sqrt(2*pi));  
   if(width==-1)
     f -= -log(sigma*pi/2) - 0.5*square(xs(j)/sigma);    
   else
     f -= -log(sigma*2*pi*(cumd_norm(width)-.5)) - 0.5*square(xs(j)/sigma);  


