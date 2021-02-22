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
    int Dmax = D;
}
parameters {
    real<lower=0> delta1;
    real<lower=0> r;
    vector<lower=0>[1] rho_gp1;
    vector<lower=0>[1] rho_gp1b;
    real<lower=0> alpha_gp1;
    real<lower=0> alpha_gp1b;
    matrix[T,D] z1;
    matrix[T,D] z1b;
}
transformed parameters {
    matrix<lower=0>[T,D] lambda;
    matrix[T,D] GP1;   
    
    {
    matrix[T, T] K1;
    matrix[T, T] L_K1;
    matrix[T, T] K1b;
    matrix[T, T] L_K1b;
    real sq_alpha1 = square(alpha_gp1);
    real sq_alpha1b = square(alpha_gp1b);
    real tmp = 0;
    
    real K1diag = sq_alpha1 + delta1 + delta0;
    real K1bdiag = sq_alpha1b + delta1 + delta0;
    
    // long part
    for (i in 1:(T-1)) {
      K1[i, i] = K1diag;
      for (j in (i + 1):T) {
        K1[i, j] = sq_alpha1
                * exp(-0.5 * dot_self((x[i] - x[j]) ./ rho_gp1)) ;
        K1[j, i] = K1[i, j];
      }
    }
    K1[T, T] = K1diag;

    // short part
    for (i in 1:(T-1)) {
      K1b[i, i] = K1bdiag;
      for (j in (i + 1):T) {
            tmp = pow(3,0.5) * sqrt(dot_self((x[i] - x[j]) ./ rho_gp1b));
            K1b[i, j] = sq_alpha1b * (1 + tmp) * exp(-1 * tmp) ;
            K1b[j, i] = K1b[i, j];
      }
    }
    K1b[T, T] = K1bdiag;
    

    L_K1 = cholesky_decompose(K1);
    L_K1b = cholesky_decompose(K1b);
    
    GP1 =  (L_K1 + L_K1b) * z1;    
    
    for( t in 1:T )
        {
            for(d in 1:D){
                lambda[t,d] = exp(GP1[t,d]);
            }
        }
    }
}
model {
    rho_gp1 ~ normal(T,0.1); // lengthscale
    alpha_gp1 ~ normal(15,2); // scale
    
    rho_gp1b ~ normal(1,0.01); // lengthscale
    alpha_gp1b ~ normal(5,1); // scale 

for (t in 1:T)
    {
        z1[t,] ~ normal(0,0.1);
        z1b[t,] ~ normal(0,0.1);    
    
    }
        
    delta1 ~ normal(0,1e-6);
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
