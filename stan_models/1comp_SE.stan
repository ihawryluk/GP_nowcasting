functions {
  vector row_sums(int[,] X) {
  int n_rows = dims(X)[1]; // this will give number of rows
  vector[n_rows] s;
  matrix [dims(X)[1], dims(X)[2]] mat_X;
  mat_X = to_matrix(X);
  for (i in 1:n_rows) s[i] = sum(row(mat_X, i));
  return s;
  }
}
data {
    int<lower=0> T;    
    int<lower=0> D;                                                                                                                                                                                                                 
    int n[T,D];
    real x[T]; // index of days (starting at 0)
}
transformed data{        
    real delta0 = 1e-10;
}
parameters {
    real<lower=0> r; // FOR THE NEG-BIN LIKELIHOOD
    
    real<lower=0> delta;
    real<lower=0> rho_gp;
    real<lower=0> alpha_gp;
    matrix[T,D] z;
}
transformed parameters {
    matrix[T,D] lambda;
    matrix[T, T] L_K;
    matrix[T, T] K = cov_exp_quad(x, alpha_gp, rho_gp);
    matrix[T,D] GP;  
    for (i in 1:T){
        K[i, i] = K[i, i] + delta0;
    }
    L_K = cholesky_decompose(K);
    GP =  L_K * z;

    for( t in 1:T )
        {
            for(d in 1:D){
                lambda[t,d] = exp(GP[t,d]);

            }
        }
}
model {
    rho_gp ~ normal(3,1);
    alpha_gp ~ normal(1,1);
    for(d in 1:D)
        z[,d] ~ normal(0,0.1);
    delta ~ normal(0,1e-6);
    // likelihood

    for( t in 1:T )
        {
        for(d in 1:D){
            if (n[t,d] < 1000000)
                {n[t,d] ~ neg_binomial_2(lambda[t,d],r);}}
        }
            
    
    // neg binomial rate
    r ~ gamma(500, 2);
    
}
generated quantities {
    vector [T] sum_n_predict;
    int n_predict[T,D];
    n_predict = n[:,:];
    
    for( t in 1:T )
        {
            for(d in 1:D)
            {
                {
                    n_predict[t,d] = neg_binomial_2_rng(lambda[t,d],r);
                }
            }
        }
    for ( t in 1:T)
    {
        sum_n_predict = row_sums(n_predict);
    }
}