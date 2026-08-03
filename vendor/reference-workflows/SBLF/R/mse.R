
#' Mean Square Error of SBLF model
#'
#' Compute and print the mean square error (MSE) from the training data
#' and mean square prediction error (MSPE) from the testing data
#'
#' @param fit SBLF fit object, save SBLF() output for this parameter
#'
mse = function(fit) {
  
  if(class(fit) != 'SBLF') {
    stop('Fit must have class SBLF.')
  }
  
  results = fit$posterior_means
  data = fit$data
  
  pred_train = results[['PostMean_Out_train']]
  pred_test = results[['PostMean_Out_test']]
  latent = results[['PostMean_Latent']]
  
  err_train = pred_train - data$ztrain
  mse_train = mean(err_train^2)  
  
  err_test <- pred_test - data$ztest
  mspe_test <- mean(err_test^2)
  
  err <- matrix(c(mse_train, mspe_test), ncol=2)
  colnames(err) <- c("MSE(training)", "MSPE(test)")
  err = as.data.frame(err)
  rownames(err) <- c('')
  
  class(err) <- c('SBLF_mse', 'data.frame')
  return(err)
  
}

print.SBLF_mse = function(x, ...) {
  print.data.frame(x)
}
