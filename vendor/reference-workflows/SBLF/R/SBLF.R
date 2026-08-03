
Rcsblf <- function(outpath, seed, burnin, iter) {
  .Call("csblf", as.character(outpath), as.integer(seed), as.integer(burnin), as.integer(iter))
}

#' Fit Spatial Bayesian Latent Factor Models
#'
#' Use this function to fit spatial Bayesian latent factor models given
#' the prepared input data sets.
#' 
#' The parameters prefixed with "x" are imaging predictor data, and
#' the parameters prefixed with "z" are imaging outcomes.
#' 
#' @param xtrain training data set of imaging predictors. Sample size is `ntrain * (P * image_len)`
#' @param xtest test data set of imaging predictors. sample size is `(ntotal-ntrain) * (P * image_len)`
#' @param ztrain training data set of imaging outcomes. Sample size is `(ntrain * image_len)`
#' @param ztest test data set of imaging outcomes. Sample size is `((ntotal-ntrain) * image_len)`
#' @param voxel_loc matrix of voxel coordinates.
#' @param seed integer random seed value
#' @param burnin integer number of burnin draws
#' @param iter integer total number of draws
#'
#' @details The function expects all data to be entered as matrices with the dimensions
#' outlined in each of the parameter arguments.
#' 
#' 
#' @return A list which contains
#' \itemize{
#'  \item{draws: }{A named list of posterior draws}
#'  \item{posterior_means: }{A named list of posterior means}
#'  \item{data: }{A named list of data sets used for model fitting}
#' }
#'
#' @examples 
#' \donttest{
#'   mod = SBLF(xtrain, xtest, ztrain, ztest, voxel_loc, seed = 1234, burnin = 250, iter = 500)
#'   mse(mod)
#' }
#'
#' @export

SBLF <- function(xtrain, xtest, ztrain, ztest, voxel_loc, seed = NULL, burnin = NULL, iter = NULL) {
  
  if(is.null(seed)) { seed = sample(1:100,1) }
  if(is.null(burnin)) { burnin = 250 }
  if(is.null(iter)) { iter = 500 }
  
  temp_path = paste(tempdir(), '/', sep = '')
  
  data_path = paste(temp_path, 'Data/', sep = '')
  result_path = paste(temp_path, 'Result/', sep = '')
  
  if(file.exists(data_path)) {
    do.call(file.remove, list(list.files(data_path, full.names = TRUE)))
    file.remove(data_path)
  }
  if(file.exists(result_path)) {
    do.call(file.remove, list(list.files(result_path, full.names = TRUE)))
    file.remove(result_path)
  }

  dir.create(data_path)
  dir.create(result_path)

  write.table(x = xtrain, file = paste(data_path, 'ROI_dat_mat.txt', sep = '/'), row.names = FALSE, col.names = FALSE)
  write.table(x = xtest, file = paste(data_path, 'ROI_dat_mat_test.txt', sep = '/'), row.names = FALSE, col.names = FALSE)
  write.table(x = ztrain, file = paste(data_path, 'ROI_task.txt', sep = '/'), row.names = FALSE, col.names = FALSE)
  write.table(x = ztest, file = paste(data_path, 'ROI_task_test.txt', sep = '/'), row.names = FALSE, col.names = FALSE)
  write.table(x = voxel_loc, file = paste(data_path, 'ROI_axes.txt', sep = '/'), row.names = FALSE, col.names = FALSE)

  image_len <- ncol(ztrain) # number of voxels per image
  write.table(x = image_len, file = paste(data_path, 'ROI_count.txt', sep = '/'), row.names = FALSE, col.names = FALSE)

  n_predictors <- ncol(xtrain) / image_len # number of imaging predictors
  n_train <- nrow(xtrain)
  n_test <- nrow(xtest)
  ROI_sizes <- matrix(data = c(1, n_train, image_len, n_predictors, n_test), ncol = 1)
  write.table(x = ROI_sizes, file = paste(data_path, 'ROI_sizes.txt', sep = '/'), row.names = FALSE, col.names = FALSE)
  cat('test')
  #Call c routine, which writes results to tmp/Result/
  creturn = Rcsblf(temp_path, seed, burnin, iter)
  
  #read results into R from tmp/Result/
  results = vector(mode = 'list', length = length(list.files(result_path)))
  names(results) <- list.files(result_path)
  X = list.files(result_path)
  Xnames = gsub('.txt', '', X)
  
  results = 
  lapply(X = X,
    FUN = function(X) {
      as.matrix(read.table(paste0(result_path, X, sep = ''), quote="\"", comment.char=""))
    })
  names(results) <- Xnames
  
  draws = results[-grep('Mean', names(results))]
  posterior_means = results[grep('Mean', names(results))]
  data = list(xtrain = xtrain, xtest = xtest, ztrain = ztrain, 
              ztest = ztest, voxel_loc = voxel_loc)
  
  # return results
  return_list = list(draws = draws, posterior_means = posterior_means, data = data)
  class(return_list) <- 'SBLF'
  return(return_list)
  
}
