

void cholesky_forward(const gsl_matrix * cholesky,
                      const gsl_vector * b,
                      gsl_vector * x, const int p){
    gsl_set_error_handler_off();
    int i, j;
    double tmpb, tmpa, tmpsum;
    gsl_vector_set(x, 0, gsl_vector_get(b, 0)/gsl_matrix_get(cholesky, 0, 0));
    for(i=1; i<p; i++){
        tmpb = gsl_vector_get(b, i);
        tmpa = gsl_matrix_get(cholesky, i, i);
        tmpsum = 0.0;
        for(j=0; j<i; j++){
            tmpsum += gsl_vector_get(x, j) * gsl_matrix_get(cholesky, i, j);
        }
        gsl_vector_set(x, i, (tmpb-tmpsum)/tmpa);
    }
    
}


// MVNsamp: function generating random samples from MVN
// with the inverse of Covariance matrix
// and a mean related vector meanM (mean = Cov * meanM)
void MVNsamp(gsl_matrix * theta, const int n,
             const gsl_vector * meanM,
             const gsl_matrix * covInv, const int p,
             const gsl_rng *r){
    gsl_set_error_handler_off();
    int i, j;
    
    // working matrix
    gsl_matrix * work = gsl_matrix_alloc(p, p);
    gsl_matrix_memcpy(work, covInv);
    
    // Cholesky Decomposition
    gsl_vector * y = gsl_vector_alloc(p);
    gsl_vector * x = gsl_vector_alloc(p);
    gsl_vector * ty = gsl_vector_alloc(p);
    
    int status = gsl_linalg_cholesky_decomp(work);
    if(status){
        Rprintf("Failed in Cholesky Decomposition\n");
    }else{

        // Solving for y: covInv * y = meanM
        gsl_linalg_cholesky_solve(work, meanM, y);

        for(i=0; i<n; i++){
            // Sampling from MVN(0, I)
            for(j=0; j<p; j++){
                gsl_vector_set(x, j, gsl_ran_ugaussian(r));
                gsl_vector_set(ty, j, 0.0);
            }
            
            // Forward substitution to solve theta-y
            cholesky_forward(work, x, ty, p);
        
            // Add y to theta
            for(j=0; j<p; j++){
                gsl_matrix_set(theta, i, j, gsl_vector_get(ty, j)+gsl_vector_get(y, j));
            }
        }
        
    }

    // Free Memory
    gsl_matrix_free(work);
    gsl_vector_free(y);
    gsl_vector_free(x);
    gsl_vector_free(ty);
    //gsl_matrix_free(X);
}

// MVN sampling with cholesky decomposition of covInv
void MVNsamp2(gsl_matrix * theta, const int n,
              const gsl_vector * meanM,
              const gsl_matrix * cholesky, const int p,
              const gsl_rng *r){
    gsl_set_error_handler_off();
    int i, j;
    
    // Cholesky Decomposition
    gsl_vector * y = gsl_vector_alloc(p);
    gsl_vector * x = gsl_vector_alloc(p);
    gsl_vector * ty = gsl_vector_alloc(p);
    
    // Solving for y: covInv * y = meanM
    gsl_linalg_cholesky_solve(cholesky, meanM, y);
    
    for(i=0; i<n; i++){
        // Sampling from MVN(0, I)
        for(j=0; j<p; j++){
            gsl_vector_set(x, j, gsl_ran_ugaussian(r));
            gsl_vector_set(ty, j, 0.0);
        }
        
        // Forward substitution to solve theta-y
        cholesky_forward(cholesky, x, ty, p);
        
        // Add y to theta
        for(j=0; j<p; j++){
            gsl_matrix_set(theta, i, j, gsl_vector_get(ty, j)+gsl_vector_get(y, j));
        }
    }
    
    // Free Memory
    gsl_vector_free(y);
    gsl_vector_free(x);
    gsl_vector_free(ty);
    
}
