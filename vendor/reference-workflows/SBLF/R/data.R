
#' Randomly selected subset of real imaging predictor data for training
#'
#' This matrix contains a randomly selected subset of real imaging predictor 
#' data for training. The data set contains 32 imaging predictors 
#' and 315 voxels per image. `32 * 315 = 10,000` columns, and 50 rows.
#'
"xtrain"

#' Randomly selected subset of real imaging predictor data for testing
#' 
#' This matrix contains a randomly selected subset of real imaging predictor 
#' data for testing. The data set contains 32 imaging predictors 
#' and 315 voxels per image. `32 * 315 = 10,000` columns, and 10 rows.
#'
#'
"xtest"

#' Randomly selected subset of real imaging outcome data for training
#' 
#' This matrix contains a randomly selected subset of real imaging outcome 
#' data for training. The data set contains 1 imaging outcome 
#' and 315 voxels per image. `1 * 315 = 315` columns, and 50 rows.
#' 
#' This outcome is from the Emotion task domain.
#'
"ztrain"

#' Randomly selected subset of real imaging outcome data for testing
#' 
#' This matrix contains a randomly selected subset of real imaging outcome 
#' data for testing. The data set contains 1 imaging outcome 
#' and 315 voxels per image. `1 * 315 = 315` columns, and 10 rows.
#' 
#' This outcome is from the Emotion task domain.
#'
"ztest"

#' Coordinates of voxel locations
#' 
#' Matrix of coordinates of voxel locations. The number of rows should match 
#' the number of voxels. In this example, there are 315 voxels.
#'
"voxel_loc"

