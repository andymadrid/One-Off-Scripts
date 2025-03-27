# Madrid Color Palettes
# a simple function that returns a vector of colors for a specified palette size
# n is a non-zero positive integer (ranging from 1-11) for the number of colors you would like returned

madridPal <- function(n) {

  if (n == 1) { return(c("cadetblue3")) }
  if (n == 2) { return(c("cadetblue3", "darkorange2")) }
  if (n == 3) { return(c("cadetblue3", "darkorange2", "darkgoldenrod2")) }
  if (n == 4) { return(c("cadetblue3", "darkorange2", "darkgoldenrod2", "darkseagreen3")) }
  if (n == 5) { return(c("cadetblue3", "darkorange2", "darkgoldenrod2", "darkseagreen3", "pink3")) }
  if (n == 6) { return(c("cadetblue3", "darkorange2", "darkgoldenrod2", "darkseagreen3", "pink3", "darkgrey")) }
  if (n == 7) { return(c("cadetblue3", "darkorange2", "pink3", "darkseagreen3", "darkgoldenrod2", "darkgrey", "indianred")) }
  if (n == 8) { return(c("cadetblue3", "darkorange2", "pink3", "darkseagreen3", "darkgoldenrod2", "darkgrey", "indianred", "thistle3")) }
  if (n == 9) { return(c("cadetblue3", "darkorange2", "pink3", "darkseagreen3", "darkgoldenrod2", "darkgrey", "indianred", "thistle3", "lightgoldenrod2")) }
  if (n == 10) { return(c("cadetblue3", "darkorange2", "pink3", "darkseagreen3", "darkgoldenrod2", "darkgrey", "indianred", "thistle3", "lightgoldenrod2", "lightcyan2")) }
  if (n == 11) { return(c("cadetblue3", "darkorange2", "pink3", "darkseagreen3", "darkgoldenrod2", "darkgrey", "indianred", "thistle3", "lightgoldenrod2", "lightcyan2", "darkolivegreen3")) }

}
