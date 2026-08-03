
/***************************************************/
/***************** GSL MATRIX INVERSE **************/
/***************************************************/

void gsl_matrix_inv(gsl_matrix *a, gsl_matrix *b)
{
    size_t n=a->size1;
    //size_t m=a->size2;
    
    gsl_matrix *temp1=gsl_matrix_calloc(n,n);
    gsl_matrix_memcpy(temp1,a);
    
    gsl_permutation *p=gsl_permutation_calloc(n);
    int sign=0;
    gsl_linalg_LU_decomp(temp1,p,&sign);
    gsl_matrix *inverse=gsl_matrix_calloc(n,n);
    
    gsl_linalg_LU_invert(temp1, p, inverse);
    gsl_matrix_memcpy(b,inverse);
    
    gsl_permutation_free(p);
    gsl_matrix_free(temp1);
    gsl_matrix_free(inverse);
}


/*****************************************************/
/******** Whether have negative Eigenvalues **********/
/*****************************************************/
int eval_neg(gsl_matrix * A, int n){
    gsl_matrix *work = gsl_matrix_alloc(n, n);
    gsl_matrix_memcpy(work, A);
    gsl_vector *eval = gsl_vector_alloc(n);
    gsl_matrix *evec = gsl_matrix_alloc(n, n);
    gsl_eigen_symmv_workspace * w = gsl_eigen_symmv_alloc(n);
    
    gsl_eigen_symmv(work, eval, evec, w);
    
    int i, check = 0;
    for(i=0; i<n; i++){
        if(gsl_vector_get(eval, i)<0)
            check = -1;
    }
    gsl_matrix_free(work);
    gsl_eigen_symmv_free(w);
    gsl_vector_free(eval);
    gsl_matrix_free(evec);
    return check;
}



/***************************************************/
/******** Sample from Multivaraite Normal **********/
/***************************************************/

void rmvnorm(const gsl_rng *r, const int n, const gsl_vector *mean, const gsl_matrix *var, gsl_vector *result, float * flag){
    gsl_set_error_handler_off();
    int k, i;
    gsl_matrix *work = gsl_matrix_alloc(n,n);
    for(i=0; i<n; i++){
        gsl_vector_set(result, i, 0.0);
    }
    
    gsl_matrix_memcpy(work,var);
    gsl_set_error_handler_off();
    int status = gsl_linalg_cholesky_decomp(work);
    if(status){
        flag[0] = 1.0;
    }else{
        for(k=0; k<n; k++){
            gsl_vector_set( result, k, gsl_ran_ugaussian(r) );
        }
        
        gsl_blas_dtrmv(CblasLower, CblasNoTrans, CblasNonUnit, work, result);
        gsl_vector_add(result, mean);
    }
    
    gsl_matrix_free(work);

}



/***************************************************/
/******** Sample from Wishart Distribution**********/
/***************************************************/

void rwishart(const gsl_rng *r, const unsigned int n, const unsigned int dof, const gsl_matrix *scale, gsl_matrix *result)
{
    gsl_set_error_handler_off();
    unsigned int k,l, i, j;
    gsl_matrix *work = gsl_matrix_calloc(n,n);
    
    for(k=0; k<n; k++){
        gsl_matrix_set( work, k, k, sqrt( gsl_ran_chisq( r, (dof-k) ) ) );
        for(l=0; l<k; l++)
            gsl_matrix_set( work, k, l, gsl_ran_ugaussian(r) );
    }
    gsl_matrix_memcpy(result,scale);
    gsl_linalg_cholesky_decomp(result);
    gsl_blas_dtrmm(CblasLeft,CblasLower,CblasNoTrans,CblasNonUnit,1.0,result,work);
    gsl_blas_dsyrk(CblasUpper,CblasNoTrans,1.0,work,0.0,result);
    
    for(i=0; i<3; i++){
        for(j=0; j<i; j++){
            gsl_matrix_set(result, i, j, gsl_matrix_get(result, j, i));
        }
    }
}



