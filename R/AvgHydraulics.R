### Geomorphic Approach ###

#' A function to execute the Geomorphic Approach hydraulic simulation routine
#'
#' This function executes the hydraulic geometry simulator to evaluate reach-averaged depths and velocities
#' generated at flows less than bankfull conditions.
#' For more information about this model see: McParland et al. (2016) and Gronsdahl et al. (XXXX)

#' @param S channel gradient (m/m)
#' @param wb reach averaged bankfull width (m).
#' @param db reach averaged bankfull depth (m).
#' @param db_max Defaults to NULL. reach averaged maximum bankfull depth (m).  Specifying db_max is preferred to calculate 'b'.
#' @param b User-specified b-value. Defaults to NULL and calculated within model unless specified.
#' @param max_Q maximum discharge (m3/s) to simualted WUA for.  Defaults to 1 m3/s.
#' @param D84 grain size (mm)
#' @param xs_output Defaults to TRUE. An expression specifying whether to produce a .csv and .jpg of the simulated channel cross section.
#' @export
#' @return .csv and .jpeg of channel cross section if specified
#' @return data frame of reach-averaged hydraulics
#' AvgHydraulics()

AvgHydraulics = function(S, wb, db, db_max = NULL, b_value = NULL, max_Q = 1,
                             D84, xs_output = TRUE) {

  # load libraries
  library(dplyr)
  library(zoo)

  ###########################################
  ##### Define Find_U Function #####
  findU = function(Wb, S, D84, depths) {

    deltaX = 0.0001
    Xgrid = Wb * seq(0, 1, deltaX)

    wet.vert = depths[depths >= 0]
    Wi = length(wet.vert) * deltaX * Wb
    Ai = sum(wet.vert * deltaX * Wb)
    di = Ai / Wi
    Pi = sum((diff(wet.vert)^2 + (max(Xgrid) * deltaX) ^2) ^ (1/2))
    Ri = Ai/Pi

    # Ferguson's continuously varying power law
    D.84 = D84 / 1000
    g = 9.81
    a1 = 6.5
    a2 = 2.5
    Res = a1 * a2 * (Ri / D.84) /
      (a1^2 + a2^2 * (Ri / D.84) ^ (5/3)) ^ (1/2)
    Ui = Res * sqrt(g * Ri * S) # Velocity (m/s)

    # Formatting the outputs in a dataframe
    df = data.frame(Ai, Wi, di, Ui)

    return(df)
  }

  #############################################################
  ##### Simulate Hydraulics #####

  # Ferguson model's shape factor (b): define based on specified inputs
  if(is.null(b_value) == FALSE){
    b = b_value
  } else if (is.null(db_max) == FALSE) {
    b = 1 - (db / db_max)
  } else {
    b = (wb / db) / 100
  }

  # stop function execution if error message too high
  try(if(b > 0.7) stop("Error: model will not produce realistic results because b-value unrealistically high"))

  # define grid
  deltaX = 0.0001 # Resolution of grid upon which to calculate Q
  # (as a proportion of wb)
  deltaY = 0.001  # increment by which to change depths when estimating HG

  # Simulate max depth if necessary
  dmax = (1 + b) / (1 - b) * db

  # generate xs_corrdinates
  X = c(0, b * wb, 0.99 * wb, wb)
  Y = 5 * db- c(0, db, dmax, 0)

  # Interpolate the distribution onto an XS raster
  Xgrid = wb * seq(0, 1, deltaX)
  Ygrid = matrix(unlist(approx(X, Y, Xgrid)), ncol = 10001, byrow = TRUE)[2,]

  # Specify the values of the water surface elevation for which to calculate Wi
  Zw = 5 * db - dmax + seq(0.02 * dmax, dmax, deltaY * dmax)

  ######################################################
  # For loop to calculate the width and discharge for each chosen water level
  simulated = data.frame(Q = NA, Ai = NA, Wi = NA, di = NA, Ui = NA)
  results = list()

  for (j in 1:length(Zw)) {
    #j = 20
    depths = Zw[j] - Ygrid # Calculate the depths, for each vertical
    results = findU(wb, S, D84, depths)
    results = c(Q = results[1, 4] * results[1, 1], results)
    simulated[j, ] = results
  }

  ## interpolate outputs to whole numbers
  Q = c(seq(0.001, 0.1, 0.001), seq(0.11, 1, 0.01), seq(1.1, 10, 0.1),
        seq(11, 100, 1), seq(110, 1000, 10), seq(1100, 10000, 100))

  # add modelled hydraulics to output dataframe
  Ai = approx(simulated$Q, simulated$Ai, xout = Q)[2]
  Wi = approx(simulated$Q, simulated$Wi, xout = Q)[2]
  di = approx(simulated$Q, simulated$di, xout = Q)[2]
  Ui = approx(simulated$Q, simulated$Ui, xout = Q)[2]

  mod_hyd = data.frame(Q, Ai = Ai$y, Wi = Wi$y, di = di$y, Ui = Ui$y) %>%
    filter(is.na(Ai) == FALSE) %>% filter(Q <= max_Q)

  #####################################################
  # Prepare graph of cross section
  if(xs_output == TRUE){

    # output coordinates
    if(is.null(db_max) == TRUE){
      plot_y = c(0, (db * -1), (dmax * - 1), 0)
    } else {
      plot_y = c(0, (db * -1), (db_max * - 1), 0)
    }

    # set up x-values to plot
    plot_x = c(0, (b * wb), (0.99 * wb), wb )

    # write channel cross section
    channel_xs = data.frame(x = plot_x, y = plot_y)
    write.csv(channel_xs, "channel_xs.csv", row.names = FALSE)

    # plot simple figure
    jpeg("channel_xs.jpeg", width = 6, height = 4, units = "in", res = 300)
    par(mar = c(4.5, 4.5, 1, 1))
    plot(plot_x, plot_y, type = "l",
         xlab = "Width (m)",
         ylab = "Depth (m)",
         ylim = c((min(plot_y) * 1.2), 0), cex.lab = 0.8, cex.axis = 0.8,
    )
    abline(a = (db * -1), 0, lty = 2, col = "grey")
    legend("bottomleft", col = c("black", "grey"), bty = "n",
           lty = c(1, 2), cex = 0.86,
           legend = c("Channel cross section", "Average depth"))
    dev.off()
  } else {
  }
  return(mod_hyd)
}


