#'Tiger
#'
#' @description Subsampled extended Yale Face Database B. This smaller dataset has 3000 samples for each
#' of the faces. Each 84x96 image patch is stored as a flattened row vector.
#'
#' @docType data
#'
#' @usage data("faces")
#'
#' @format An object of class \code{"spca"}.
#'
#' @keywords face
#'
#' @references Yale Face Database B (Online)
#'
#' @source \href{http://vision.ucsd.edu/~leekc/ExtYaleDatabase/Yale\%20Face\%20Database.htm}{http://vision.ucsd.edu/~leekc/ExtYaleDatabase/Yale\%20Face\%20Database.htm}
#'
#' @examples
#' \dontrun{
#' library("spca")
#' data("faces")
#'
#' #Display first digit
#' face <- matrix(faces[,1], nrow = 84, ncol = 96)
#' image(face[,96:1], col = gray(255:0 / 255))
#' }   
#'
#'
"faces"
