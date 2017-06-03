conceptPlot = function(
  bins,
  compoundBoxWidth = 0.2,
  flowPeriod = 3600,
  binGapFactor = 0.05,
  backgroundCol = "black",
  textCol = "white",
  axisFontSize = 12,
  yAxisMax = 0
){

  #!!!!!! reverse the row order of the bins array !!!!!!
  bins = bins[nrow(bins):1,]

  # set the gap height, relative to water storage.
  gapHeight = binGapFactor * sum(bins$waterStorage)

  # calculate the x value (center) for the bins
  bins$rowCenter = cumsum(bins$waterStorage) - bins$waterStorage/2
  bins$rowCenter = bins$rowCenter + 0:(nrow(bins)-1) * gapHeight

  # calculate the max y value (top) of each bin
  bins$rowTop = bins$rowCenter + bins$waterStorage/2

  # set some basics for the plot space
  HZ = ggplot2::ggplot()
  HZ = HZ +
    ggplot2::theme(panel.grid.major = ggplot2::element_blank(), panel.grid.minor = ggplot2::element_blank()) +
    ggplot2::theme(plot.background = ggplot2::element_rect(fill = backgroundCol, colour = backgroundCol), panel.background = ggplot2::element_blank(), axis.text = ggplot2::element_blank())
  HZ = HZ + ggplot2::coord_cartesian(xlim = c(1-compoundBoxWidth/1.9, 1+compoundBoxWidth*3/2))
  if(yAxisMax > 0) {
    HZ = HZ + ggplot2::coord_cartesian(ylim = c(0, yAxisMax))
  }

  # deterning the number of bins (rows)
  nr = nrow(bins)

  ##BINS -- bins are located at x = 1, using the compoundBoxWidth as a width.
  ##Row center, calculated above, is the y-value of the middle of the bin
  HZ = HZ +
    ggplot2::geom_tile(
      data = bins,
      mapping = ggplot2::aes(
        x=1,
        y=rowCenter,
        height = waterStorage,
        width = compoundBoxWidth#,
        #        fill = element
      ),
      #fill = grey(0.7)
      fill = rainbow(nr, end = (nr-6)/nr)
    )

  # calculate flow data
  idx = 1:(nr - 1)
  nIdx = idx + 1
  flowData = data.frame(
    centers = bins$rowTop + gapHeight/2, #c((bins$rowCenter[idx]+bins$rowCenter[nIdx])/2, bins$rowCenter[nr] + (bins$waterStorage[nr]+gapHeight)/2),
    height = gapHeight, #c((bins$waterStorage[idx]+bins$waterStorage[nIdx])/2 , bins$waterStorage[nr]),
    flow = bins$entering*flowPeriod,
    upwelling = bins$returning*flowPeriod
  )
  flowData$returnCenter =  1 + compoundBoxWidth*0.55 + (flowData$flow[nr] - cumsum(flowData$upwelling) + flowData$upwelling/2)*compoundBoxWidth/flowData$height
  #flowData$returnCenter =  1 + compoundBoxWidth*1.05/2 + sum(flowData$upwelling) - cumsum((c(0 , flowData$upwelling[-nr]) + flowData$upwelling)/2)

  ##RETURN FLOW
  HZ = HZ +
    ggplot2::geom_tile(
      data = flowData,
      mapping = ggplot2::aes(
        x = returnCenter,
        y = centers[nr], #(bins$rowCenter + centers[nr] + height[nr]/2)/2,
        height = height[nr], #(centers[nr] - bins$rowCenter) + height[nr]/2,
        width = upwelling*compoundBoxWidth/height
      ),
      #color = "cyan4",
      #fill = "cyan",
      #size = 1
      fill = rainbow(nr, end = (nr-6)/nr)
    )

  ##DOWNWELLING
  HZ = HZ +
    ggplot2::geom_tile(
      data = flowData,
      mapping = ggplot2::aes(
        x = 1, #- compoundBoxWidth*1.05/2 - flow/2,
        y = centers,
        height = height,
        width = flow*compoundBoxWidth/height
      ),
      fill = "cyan",
      color = "cyan4"
      #size = 1
      #fill = c(rainbow(nr, end = (nr-6)/nr)[-1], "blue"),
      #color = "black"
    )

  ## HORIZONTAL LINES
  HZ = HZ +
    ggplot2::geom_segment(
      data = flowData,
      mapping = ggplot2::aes(
        x = 1 + compoundBoxWidth/2,
        xend =  returnCenter,
        y = bins$rowCenter,
        yend = bins$rowCenter
      ),
      color = rainbow(nr, end = (nr-6)/nr),
      #color = "cyan",
      size = 1
    )

  ## VERTICAL LINES
  HZ = HZ +
    ggplot2::geom_segment(
      data = flowData,
      mapping = ggplot2::aes(
        x = returnCenter,
        xend = returnCenter,
        y = bins$rowCenter,
        yend = centers[nr] - height[nr]/2
      ),
      color = rainbow(nr, end = (nr-6)/nr),
      #color = "cyan",
      size = 1
    )

  HZ
  # y.axisDF = plyr::ddply(poolDF, "to", plyr::summarise,  y = median(y))
  # y.axisDF = y.axisDF[with(y.axisDF, order(y)), ]
  #
  # river = river +
  #   ggplot2::scale_y_continuous(breaks = y.axisDF$y, labels = y.axisDF$to, name = "") +
  #   ggplot2::scale_x_continuous(breaks = seq(1.5,9.5,1), labels = as.character(seq(1, 9, 1)), name = "") +
  #   ggplot2::theme(axis.text = ggplot2::element_text(colour = textCol, size = ggplot2::rel(axisFontSize))) +
  #   ggplot2::theme(
  #     legend.background = ggplot2::element_rect(fill = backgroundCol),
  #     legend.text = ggplot2::element_text(size = ggplot2::rel(1.5), colour = backgroundCol),
  #     legend.position = "none"
  #     # legend.position = "top"
  #   )
  #
  # river = river + ggplot2::theme(plot.margin = grid::unit(c(0.5, 1, -0.5, 0.5), "lines"))
  # if(makePDF == TRUE) {
  #   pdf(
  #     fileName,
  #     onefile = T,
  #     paper = "USr",
  #     # paper = "letter"
  #     width = 10.5,
  #     height = 7
  #   )
  #   print(river)
  #   dev.off()
  # } else{
  #   return(river)
  # }
  # river
}
