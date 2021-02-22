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
    int D1 = D-1;
    int Dmax = D;
}
parameters {
    real<lower=0> delta1;
    real<lower=0> delta2;
    real<lower=0> r;
    vector<lower=0>[1] rho_gp1;
    vector<lower=0>[1] rho_gp2;
    vector<lower=0>[1] rho_gp1b;
    vector<lower=0>[1] rho_gp2b;
    real<lower=0> alpha_gp1;
    real<lower=0> alpha_gp2;
    real<lower=0> alpha_gp1b;
    real<lower=0> alpha_gp2b;
    matrix[T,1] z1;
    matrix[T,D1] z2;
    matrix[T,1] z1b;
    matrix[T,D1] z2b;
}
transformed parameters {
    matrix<lower=0>[T,D] lambda;
    matrix[T,1] GP1;   
    matrix[T,D1] GP2;
    
    {
    matrix[T, T] K1;
    matrix[T, T] K2;
    matrix[T, T] L_K1;
    matrix[T, T] L_K2;
    matrix[T, T] K1b;
    matrix[T, T] K2b;
    matrix[T, T] L_K1b;
    matrix[T, T] L_K2b;
    real sq_alpha1 = square(alpha_gp1);
    real sq_alpha2 = square(alpha_gp2);
    real sq_alpha1b = square(alpha_gp1b);
    real sq_alpha2b = square(alpha_gp2b);
    
    real K1diag = sq_alpha1 + delta1 + delta0;
    real K2diag = sq_alpha2 + delta2 + delta0;
    real K1bdiag = sq_alpha1b + delta1 + delta0;
    real K2bdiag = sq_alpha2b + delta2 + delta0;
    
    for (i in 1:(T-1)) {
      K1[i, i] = K1diag;
      for (j in (i + 1):T) {
        K1[i, j] = sq_alpha1
                * exp(-0.5 * dot_self((x[i] - x[j]) ./ rho_gp1)) ;
        K1[j, i] = K1[i, j];
      }
    }
    K1[T, T] = K1diag;
    
    for (i in 1:(T-1)) {
      K2[i, i] = K2diag;
      for (j in (i + 1):T) { 
        K2[i, j] = sq_alpha2
                * exp(-0.5 * dot_self((x[i] - x[j]) ./ rho_gp2)) ;
        K2[j, i] = K2[i, j];
      }
    }
    K2[T, T] = K2diag;
    
    
    
    
    for (i in 1:(T-1)) {
      K1b[i, i] = K1bdiag;
      for (j in (i + 1):T) {
        K1b[i, j] = sq_alpha1b
                * exp(-0.5 * dot_self((x[i] - x[j]) ./ rho_gp1b)) ;
        K1b[j, i] = K1b[i, j];
      }
    }
    K1b[T, T] = K1bdiag;
    
    for (i in 1:(T-1)) {
      K2b[i, i] = K2bdiag;
      for (j in (i + 1):T) { 
        K2b[i, j] = sq_alpha2b
                * exp(-0.5 * dot_self((x[i] - x[j]) ./ rho_gp2b)) ;
        K2b[j, i] = K2b[i, j];
      }
    }
    K2b[T, T] = K2bdiag;
    

     
    L_K1 = cholesky_decompose(K1);
    L_K2 = cholesky_decompose(K2);
    L_K1b = cholesky_decompose(K1b);
    L_K2b = cholesky_decompose(K2b);
    GP1 =  (L_K1 + L_K1b) * z1;
    GP2 =  (L_K2 + L_K2b) * z2;
    
    
    for( t in 1:T )
        {
            lambda[t,1] = exp(GP1[t,1]);

            for(d in 2:D){
                lambda[t,d] = exp(GP2[t,d-1]);
            }
        }
    }
}
model {
    rho_gp1 ~ normal(T,0.1); // lengthscale
    rho_gp2 ~ normal(D,0.1); // lengthscale
    alpha_gp1 ~ normal(15*3,2); // scale
    alpha_gp2 ~ normal(20*3,2); // scale
    
    rho_gp1b ~ normal(1,0.01); // lengthscale
    rho_gp2b ~ normal(1,0.01); // lengthscale
    alpha_gp1b ~ normal(5*3,1); // scale
    alpha_gp2b ~ normal(1.5*3,1); // scale

    z1[,1] ~ normal(0,0.1);
    z1b[,1] ~ normal(0,0.1);
    for (t in 1:T)
    {
        z2[t,] ~ normal(0,0.1);
        z2b[t,] ~ normal(0,0.1);    
    
    }

        
    delta1 ~ normal(0,1e-6);
    delta2 ~ normal(0,1e-6);
    r ~ gamma(500, 2);
    
    
    for( t in 1:T )
        {
        for(d in 1:Dmax){
           if (n[t,d] < 1000000)
                {n[t,d] ~ neg_binomial_2(lambda[t,d],r);}}// NEG BIN
        }

}
generated quantities {
    vector [T] sum_n_predict;
    int n_predict[T,Dmax];

    for( t in 1:T )
        {
            for(d in 1:Dmax)
            {
                    n_predict[t,d] = neg_binomial_2_rng(lambda[t,d],r);
            }
        }
    sum_n_predict = row_sums(n_predict);
}
