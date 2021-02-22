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
}
parameters {
    real<lower=0> alpha [T];
    simplex[D] beta;
    real<lower=0> r; // FOR THE NEG-BIN LIKELIHOOD
    real<lower=0> tau2alpha; // this is variance so has to be positive
    

}
transformed parameters {
    real<lower=0> lambda[T,D];    
    for( t in 1:T )
        {
            for(d in 1:D){
                lambda[t,d] = exp(alpha[t]) * beta[d];

            }
        }
}
model {
    
    // likelihood

    for( t in 1:T )
        {
        for(d in 1:D){
            if (n[t,d] < 1000000)
                {n[t,d] ~ neg_binomial_2(lambda[t,d],r);}}// NEG BIN
        }
            
    // priors
    
    // neg binomial rate
    r ~ gamma(500, 2);  // FOR THE NEG-BIN LIKELIHOOD
    
    // alpha
    tau2alpha ~ gamma(0.01, 0.01);
    alpha[1] ~ normal(0, sqrt(1/0.001));
    //alpha[2] ~ normal(0, sqrt(1/0.001));

    for( t in 2:T ){
      //alpha[t] ~ normal(sqrt(alpha[t-1]*alpha[t-2]),sqrt(1/tau2alpha));}
    alpha[t] ~ normal(alpha[t-1],sqrt(1/tau2alpha));}


    // beta
    beta ~ dirichlet(rep_vector(0.1, D));
    
}
generated quantities {
    vector [T] sum_n_predict;
    int n_predict[T,D];
    n_predict = n[:,:];
    
    for( t in 1:T )
        {
            for(d in 1:D)
            {
                n_predict[t,d] = neg_binomial_2_rng(lambda[t,d],r);
            }
        }
    for ( t in 1:T)
    {
        sum_n_predict = row_sums(n_predict);
    }
}