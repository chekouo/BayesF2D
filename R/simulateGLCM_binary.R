#'
#' Simulate GLCMs for binary response
#'
#' Simulate symmetric or non-symmetric GLCMs for binary response data
#' 
#' @param n0 number of samples for y=0 for binary data; defaults to 10
#' @param n1 number of samples for y=1 for binary data; defaults to 10
#' @param m0 list of 4 2x1 means for bivariate normal indicating the highest intensity coordinates in GLCMs for y=0
#' @param S0 list of 4 2x2 covariance matrices of the bivariate normal for y=0
#' @param m1 list of 4 2x1 means for bivariate normal indicating the highest intensity coordinates in GLCMs for y=1
#' @param S1 list of 4 2x2 covariance matrices of the bivariate normal for y=1
#' @param GLCM.type 'sym' or 'nonsym' for symmateric and non-symmetric GLCMs; defaults to symmetric GLCM
#' @param Gam list of 4x4 correlation matrices - 64 non 'nonsym' case and 36 for 'sym' case; defaults to covariance matrices provided as data in the package
#' @param noise noise defaults to 1
#' @param counts a 'n0+n1' vector for the number of pixels assumed in the image; First n0 belong to y=0 and last n1 belong to y=1
#' @export

simulateGLCM_binary = function(n0 = 10, # number of samples for y=0 for binary data; defaults to 10
                               n1 = 10, # number of samples for y=1 for binary data; defaults to 10
                               m0, # list of 4 2x1 means for bivariate normal indicating the highest intensity coordinates in GLCMs for y=0
                               S0, # list of 4 2x2 covariance matrices of the bivariate normal for y=0
                               m1, # list of 4 2x1 means for bivariate normal indicating the highest intensity coordinates in GLCMs for y=1
                               S1, # list of 4 covariance matrices of the bivariate normal for y=1
                               GLCM.type = 'sym', # options: 'sym' or 'nonsym' for symmateric and non-symmetric GLCMs; defaults to symmetric GLCM
                               Gam = NULL, # list of 4x4 correlation matrices - 64 non 'nonsym' case and 36 for 'sym' case; defaults to covariance matrices provided as data in the package
                               noise = 1, # noise defaults to 1
                               counts = NULL # a 'n0+n1' vector for the number of pixels assumed in the image; First n0 belong to y=0 and last n1 belong to y=1
                               ){
  # number of gray levels; defaults to 8
  g_size = 8
  
  # Initialize counts isf not specified
  if(is.null(counts)) counts = sample(4e4:3e6, size = n0+n1)
  
  # If no list of covariance matrices (GAM) supplied use the default in the package
  if(is.null(Gam)) Gam = eval(parse(text=paste('BayesF2D::', GLCM.type, 'Cor', sep=''))) # Replace PACKAGENAME
  
  # response
  y = rep(c(0, 1), c(n0, n1))
  
  # Load the function to simulate one GLCM based on GLCM.type
  # source(paste('./', GLCM.type, 'GLCM.R', sep=''))
  # Assign the name 'Fn' to the function
  # Fn = eval(parse(text=paste(GLCM.type, 'GLCM', sep='')))
  
  # Use the next line instead of the previous two lines when using the package
  Fn = eval(parse(text=paste('BayesF2D::', GLCM.type, 'GLCM', sep=''))) # Replace PACKAGENAME
  
  # simulate GLCMs for y=0
  n0.GLCM = lapply(1:n0, function(i) Fn(m = m0, S = S0, Gam = Gam, noise = noise, g_size = g_size, count = counts[i]))
  # simulate GLCMs for y=1
  n1.GLCM = lapply(1:n1, function(i) Fn(m = m1, S = S1, Gam = Gam, noise = noise, g_size = g_size, count = counts[n0+i]))
  
  # Append GLCMs for y=0 and y=1
  GLCM = n0.GLCM
  for(i in 1:n1) GLCM[[n0+i]] = n1.GLCM[[i]]
  
  # return the GLCMs and the response
  list(GLCM = GLCM,
       y = y)
}
