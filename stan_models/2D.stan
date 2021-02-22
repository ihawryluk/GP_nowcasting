functions {
  vector row_sums(int[,] X) 
  {
      int n_rows = dims(X)[1]; // this will give number of rows
      vector[n_rows] s;
      matrix [dims(X)[1], dims(X)[2]] mat_X;
      mat_X = to_matrix(X);
      for (i in 1:n_rows) s[i] = sum(row(mat_X, i));
      return s;
  }
  

    matrix kron_mvprod(matrix A, matrix B, matrix V) 
    {
        return transpose(A*transpose(B*V));
    }
  
    matrix calculate_eigenvalues(vector A, vector B, int n1, int n2, real sigma2) 
    {
        matrix[n1,n2] e;
            for(i in 1:n1) 
            {
                for(j in 1:n2) 
                {
                    e[i,j] = (A[i]*B[j] + sigma2);
                }
            }
        return(e);
    }
}


data 
{
    int<lower=0> T;    
    int<lower=0> D;                                                                                                                                                                                                                 
    int n[T,D];
    real time[T]; 
    real delay[D];
}


transformed data
{   
    real delta0 = 1e-9;  
    
    vector[1] rho_gp1_t_dist; 
    vector[1] rho_gp2_t_dist;
    vector[1] rho_gp1_d_dist;
    vector[1] rho_gp2_d_dist;

    rho_gp1_t_dist[1] = T;
    rho_gp2_t_dist[1] = 1;
    rho_gp1_d_dist[1] = D;
    rho_gp2_d_dist[1] = 1;
}


parameters {
    real<lower=0> r; // binomial dispersion
    real<lower=0> alpha_gp1_t;
    real<lower=0> alpha_gp2_t;
    real<lower=0> alpha_gp1_d;
    real<lower=0> alpha_gp2_d;
    real<lower=0> delta1;
    real<lower=0> delta2;
//    vector<lower=0>[1] rho_gp1_t_dist; 
//    vector<lower=0>[1] rho_gp2_t_dist;
//    vector<lower=0>[1] rho_gp1_d_dist;
//    vector<lower=0>[1] rho_gp2_d_dist;
    matrix[T,D] z1;
    matrix[T,D] z2;
}


transformed parameters {
    matrix<lower=0, upper=50000>[T,D] lambda;
    matrix<lower=0>[T,D] lambda_gp1;
    matrix<lower=0>[T,D] lambda_gp2;
    
    matrix[T,D] GP1;//long range
    matrix[T,D] GP2;//short range
    
    matrix[T, T] K1_t;
    matrix[T, T] K2_t;
    matrix[T, T] L_K1_t;
    matrix[T, T] L_K2_t;
    
    matrix[D, D] K1_d;
    matrix[D, D] K2_d;
    matrix[D, D] L_K1_d;
    matrix[D, D] L_K2_d;
    
    real sq_alpha1_t = square(alpha_gp1_t);
    real sq_alpha2_t = square(alpha_gp2_t);
    real sq_alpha1_d = square(alpha_gp1_d);
    real sq_alpha2_d = square(alpha_gp2_d);
    
    real K1diag_t = sq_alpha1_t + delta0 + delta1;
    real K2diag_t = sq_alpha2_t + delta0 ;
    
    real K1diag_d = sq_alpha1_d + delta0;
    real K2diag_d = sq_alpha2_d + delta0 + delta2;
    
 //time   
    for (i in 1:T) 
        {
            K1_t[i, i] = K1diag_t;
            K2_t[i, i] = K2diag_t;
            for (j in (i + 1):T) 
                {
                    K1_t[i, j] = sq_alpha1_t
                            * exp(-0.5 * dot_self((time[i] - time[j]) ./ rho_gp1_t_dist)) ; //long range sq exp                         
                    K2_t[i, j] = sq_alpha2_t
                            * exp( -0.5 * dot_self((time[i] - time[j]) ./ rho_gp2_t_dist));// shortrange  sq exp
                    K1_t[j, i] = K1_t[i, j];
                    K2_t[j, i] = K2_t[i, j];
                }
        }
    K1_t[T, T] = K1diag_t;
    K2_t[T, T] = K2diag_t;
    
 //delay   
    for (i in 1:D) 
        {
            K1_d[i, i] = K1diag_d;
            K2_d[i, i] = K2diag_d;
            for (j in (i + 1):D) 
                {
                    K1_d[i, j] = sq_alpha1_d
                            * exp(-0.5 * dot_self((delay[i] - delay[j]) ./ rho_gp1_d_dist)) ; //long range sq exp
                    K2_d[i, j] = sq_alpha2_d
                            * exp( -0.5 * dot_self((delay[i] - delay[j]) ./ rho_gp2_d_dist));// shortrange  sq exp
                    K1_d[j, i] = K1_d[i, j];
                    K2_d[j, i] = K2_d[i, j];
                }
        }
    K1_d[D, D] = K1diag_d;
    K2_d[D, D] = K2diag_d;

    L_K1_t = cholesky_decompose(K1_t);
    L_K1_d = cholesky_decompose(K1_d);
    L_K2_t = cholesky_decompose(K2_t);
    L_K2_d = cholesky_decompose(K2_d);
    
    GP1 = kron_mvprod(L_K1_d, L_K1_t, z1);
    GP2 = kron_mvprod(L_K2_d, L_K2_t, z2);

    for (t in 1:T)
        {
            for(d in 1:D)
                {
                    lambda_gp1[t,d] = exp(GP1[t,d]);
                    lambda_gp2[t,d] = exp(GP2[t,d]);
                    lambda[t,d] = lambda_gp1[t,d] * lambda_gp2[t,d];
                }
        }
}


model 
{
    delta1 ~ normal(0,1e-7);
    delta2 ~ normal(0,1e-7);
    
    alpha_gp1_t ~ normal(T,5);
    alpha_gp1_d ~ normal(0,1);
    
    alpha_gp2_t ~ normal(D,5);
    alpha_gp2_d ~ normal(0,1);
    

    for(t in 1:T)
    {
        for(d in 1:D)
            {
                {
                   z1[t,d] ~ normal(0,0.1);
                   z2[t,d] ~ normal(0,0.1);
                }
            }
    }
    
    r ~  gamma(400,2);
    for( t in 1:T )
        {
            for(d in 1:D)
                {
                if (n[t,d] < 1000000)
                    {
                       n[t,d] ~ neg_binomial_2(lambda[t,d],r);
                    }
                }
        }
}


generated quantities 
{
    vector [T] sum_n_predict;
    vector [T] sum_n_predict_gp1;
    vector [T] sum_n_predict_gp2;
    int n_predict[T,D];
    int n_predict_alt[T,D];
    int n_predict_gp1[T,D];
    int n_predict_gp2[T,D];


    for(t in 1:T)
        {
            for(d in 1:D)
                {
                   n_predict[t,d] = neg_binomial_2_rng(lambda[t,d],r);
                   n_predict_gp1[t,d] = neg_binomial_2_rng(lambda_gp1[t,d],(r));
                   n_predict_gp2[t,d] = neg_binomial_2_rng(lambda_gp2[t,d],(r));
                }
        }
    sum_n_predict = row_sums(n_predict);
    sum_n_predict_gp1 = row_sums(n_predict_gp1);
    sum_n_predict_gp2 = row_sums(n_predict_gp2);
}
