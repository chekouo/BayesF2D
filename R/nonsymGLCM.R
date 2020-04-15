#'
#' Simulate non-symmetric GLCMs
#'
#' Simulate non-symmetric GLCMs which are correlated across imaging modalities for one subject
#' 
#' @param m list of 4 2x1 means bivariate indicating the highest intensity coordinates in GLCMs
#' @param S list of 4 2x2 covariance matrices of the bivariate normal
#' @param Gam list of 64 4x4 correlation matrices
#' @param noise noise
#' @param g_size number of gray levels
#' @param count number of pixels in the image
#' @export
#' 

nonsymGLCM = function(m, # list of 4 2x1 means bivariate indicating the highest intensity coordinates in GLCMs
                      S, # list of 4 2x2 covariance matrices of the bivariate normal
                      Gam, # list of 64 4x4 correlation matrices
                      noise, # noise
                      g_size, # number of gray levels
                      count # number of pixels in the image
                      ){
  # generate correlated errors for entrywise correlation across imaging modalities
  eps = lapply(1:(g_size^2), function(l) MASS::mvrnorm(1, mu = rep(0, 4), Sigma = noise*Gam[[l]]))
  
  # Generate GLCM for each imaging modality using the algorithm in the manuscript Chekuou et al.
  listGLCM = lapply(1:4,
                    function(r){
                      temp_term = matrix(0,g_size,g_size)
                      l = 0
                      for(k in 1:g_size) for(j in 1:g_size){
                        l = l+1
                        # (k-1,g_size-j), (k, g_size-j), (k, g_size-j+1), (k-1,g_size-j+1)
                        a1 = mvtnorm::pmvnorm(upper = c(k, g_size-j+1),
                                              mean = m[[r]], sigma = S[[r]])
                        a2 = mvtnorm::pmvnorm(upper = c(k, g_size-j),
                                              mean = m[[r]], sigma = S[[r]])
                        a3 = mvtnorm::pmvnorm(upper = c(k-1,g_size-j+1),
                                              mean = m[[r]], sigma = S[[r]])
                        a4 = mvtnorm::pmvnorm(upper = c(k-1,g_size-j),
                                              mean = m[[r]], sigma = S[[r]])
                        p_jk = a1-a2-a3+a4
                        temp_term[j,k] = log(p_jk)+eps[[l]][r]
                      }
                      round(count*exp(temp_term))
                    })
  
  # return the GLCMs as a list
  outGLCM = array(dim=c(4,g_size,g_size))
  for(r in 1:4) outGLCM[r,,] = listGLCM[[r]]
  outGLCM
}