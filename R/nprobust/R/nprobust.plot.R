nprobust.plot <- function(..., alpha=NULL, type=NULL, CItype=NULL,
                          title="", xlabel="", ylabel="",
                          lty=NULL, lwd=NULL, lcol=NULL, pty=NULL, pwd=NULL, pcol=NULL,
                          CIshade=NULL, CIcol=NULL, legendTitle=NULL, legendGroups=NULL) {

  ########################################
  # check how many series are passed in 
  ########################################

  x <- list(...)
  nfig <- length(x)
  if (nfig == 0) stop("Nothing to plot.\n")

  ########################################
  # error handling
  ########################################
  # alpha
  if (length(alpha) == 0) {
    alpha <- rep(0.05, nfig)
  } else if (!all(alpha>0 & alpha<1)) {
    stop("Significance level incorrectly specified.\n")
  } else {
    alpha <- rep(alpha, length.out=nfig)
  }

  # plot type
  if (length(type) == 0) {
    type <- rep("line", nfig)
  } else {
    if (!all(type%in%c("line", "points", "both"))) {
      stop("Plotting type incorrectly specified.\n")
    }
    type <- rep(type, length.out=nfig)
  } 

  # CI type
  if (length(CItype) == 0) {
    CItype <- rep("region", nfig)
  } else {
    if (!all(CItype%in%c("region", "line", "ebar", "all", "none"))) {
      stop("Confidence interval type incorrectly specified.\n")
    }
    CItype <- rep(CItype, length.out=nfig)
  }

  # line style, line width, line color
  if (length(lty) == 0) {
    lty <- rep(1, nfig)
  } else {
    lty <- rep(lty, length.out=nfig)
  }
  if (length(lwd) == 0) {
    lwd <- rep(0.5, nfig)
  } else {
    lwd <- rep(lwd, length.out=nfig)
  }
  if (length(lcol) == 0) {
    lcol <- 1:nfig
  } else {
    lcol <- rep(lcol, length.out=nfig)
  }

  # point style, point width, point color
  if (length(pty) == 0) {
    pty <- rep(1, nfig)
  } else {
    pty <- rep(pty, length.out=nfig)
  }
  if (length(pwd) == 0) {
    pwd <- rep(1, nfig)
  } else {
    pwd <- rep(pwd, length.out=nfig)
  }
  if (length(pcol) == 0) {
    pcol <- lcol
  } else {
    pcol <- rep(pcol, length.out=nfig)
  }

  # CI shade, CI color
  if (length(CIshade) == 0) {
    CIshade <- rep(0.2, nfig)
  } else {
    CIshade <- rep(CIshade, length.out=nfig)
  }
  if (length(CIcol) == 0) {
    CIcol <- lcol
  } else {
    CIcol <- rep(CIcol, length.out=nfig)
  }

  # legend
  # Changes made by Xinwei
  if (length(legendTitle) == 0) {
    legendTitle <- ""
  } else {
    legendTitle <- legendTitle[1]
  }
  if (length(legendGroups) > 0) {
    legendGroups <- rep(legendGroups, length.out=nfig)
    legend_default <- FALSE
  } else {
    legend_default <- TRUE
  }

  ########################################
  # initializing plot
  ########################################
  temp_plot <- ggplot() + theme_bw() #+ theme(legend.position="none")

  CI_l <- CI_r <- tau.us <- eval <- Sname <- NULL

  ########################################
  # looping over input models
  ########################################
  ### all colors
  col_all <- lty_all <- pty_all <- c()
  for (i in 1:nfig) {
    data_x <- data.frame(x[[i]]$Estimate[, c("eval", "tau.us", "tau.bc", "se.us", "se.rb")])
    z_val <- qnorm(1 - alpha[i]/2)
    data_x$CI_l <- data_x$tau.bc - z_val * data_x$se.rb
    data_x$CI_r <- data_x$tau.bc + z_val * data_x$se.rb

    # changes made by Xinwei
    if (legend_default) {
      data_x$Sname <- paste("Series", i, sep=" ")
      legendGroups <- c(legendGroups, data_x$Sname)
    } else {
      data_x$Sname <- legendGroups[i]
    }

    ########################################
    # add CI regions to the plot
    if (CItype[i]%in%c("region", "all"))
      temp_plot <- temp_plot + geom_ribbon(data=data_x, aes(x=eval, ymin=CI_l, ymax=CI_r), alpha=CIshade[i], fill=CIcol[i])

    ########################################
    # add CI lines to the plot
    if (CItype[i]%in%c("line", "all"))
      temp_plot <- temp_plot + geom_line(data=data_x, aes(x=eval, y=CI_l), linetype=2, alpha=CIshade[i], col=CIcol[i]) +
      geom_line(data=data_x, aes(x=eval, y=CI_r), linetype=2, alpha=CIshade[i], col=CIcol[i])

    ########################################
    # add error bars to the plot
    if (CItype[i]%in%c("ebar", "all"))
      temp_plot <- temp_plot + geom_errorbar(data=data_x, aes(x=eval, ymin=CI_l, ymax=CI_r), alpha=CIshade[i], col=CIcol[i], linetype=1)

    ########################################
    # add lines to the plot
    if (type[i]%in%c("line", "both")) {
      temp_plot <- temp_plot + geom_line(data=data_x, aes(x=eval, y=tau.us, colour=Sname, linetype=Sname), size=lwd[i])
    }

    ########################################
    # add points to the plot
    if (type[i]%in%c("points", "both")) {
      temp_plot <- temp_plot + geom_point(data=data_x, aes(x=eval, y=tau.us, colour=Sname, shape=Sname), size=pwd[i])
    }

    if (type[i] == "line") {
      col_all <- c(col_all, lcol[i])
      lty_all <- c(lty_all, lty[i])
      pty_all <- c(pty_all, NA)
    } else if (type[i] == "both") {
      col_all <- c(col_all, lcol[i])
      lty_all <- c(lty_all, lty[i])
      pty_all <- c(pty_all, pty[i])
    } else {
      col_all <- c(col_all, pcol[i])
      lty_all <- c(lty_all, NA)
      pty_all <- c(pty_all, pty[i])
    }
  }

  ########################################
  # change color, line type and point shape back, and customize legend
  ########################################
  
  index <- sort.int(legendGroups, index.return=TRUE)$ix
  temp_plot <- temp_plot + scale_color_manual(values = col_all[index]) +
    scale_linetype_manual(values = lty_all[index]) +
    scale_shape_manual(values = pty_all[index]) +
    guides(colour=guide_legend(title=legendTitle)) +
    guides(linetype=guide_legend(title=legendTitle)) +
    guides(shape=guide_legend(title=legendTitle))

  ########################################
  # add title, x and y labs
  ########################################
  temp_plot <- temp_plot + labs(x=xlabel, y=ylabel) + ggtitle(title)
 
  ########################################
  # return the plot
  ########################################
  return (temp_plot)
}
