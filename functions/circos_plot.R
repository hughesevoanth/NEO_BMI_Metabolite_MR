circos_plot =function (track_number, track1_data, track2_data, track3_data, 
    track1_type, track2_type, track3_type, label_column, section_column, 
    order = TRUE, order_column, estimate_column, pvalue_column, 
    pvalue_adjustment, lower_ci, upper_ci, lines_column, lines_type = "o", 
    bar_column, histogram_column, histogram_binsize = 0.01, histogram_densityplot = FALSE, 
    legend = FALSE, track1_label = NA, track2_label = NA, track3_label = NA, 
    pvalue_label = NA, circle_size = 25, track1_height = 0.2, 
    track2_height = 0.2, track3_height = 0.2, equal_axis = FALSE) 
{
    track1 <- 1
    track2 <- 2
    track3 <- 3
    track4 <- 4
    x_axis_index <- 1
    track_axis_reference <- 0
    margins <- c(0.5, 0.5, 0.5, 0.5) * 25
    start_gap <- 17
    start_degree <- 90
    section_track_height <- 0.1
    #discrete_palette <- c("#00378f", "#ffc067", "#894300")
    section_fill_colour <- "snow2"
    section_text_colour <- "black"
    section_line_colour <- "grey"
    section_line_thickness <- 1.5
    section_line_type <- 1
    reference_line_colour <- "grey50"
    reference_line_thickness <- 1.5
    reference_line_type <- 1
    guide_line_colour <- "grey80"
    guide_line_thickness <- 0.8 
    point_pch <- 21
    point_cex <- 1.2
    point_col1 <- discrete_palette[1]
    point_bg1 <- "white"
    #point_col1_sig <- "white"
    point_col1_sig <- discrete_palette[1]
    point_bg1_sig <- discrete_palette[1]
    #point_bg1_sig <- "black"
    point_col2 <- discrete_palette[2]
    point_bg2 <- "white"
    #point_col2_sig <- "white"
    point_col2_sig <- discrete_palette[2]
    point_bg2_sig <- discrete_palette[2]
    #point_bg2_sig <- "black"
    point_col3 <- discrete_palette[3]
    point_bg3 <- "white"
    #point_col3_sig <- "white"
    point_col3_sig <- discrete_palette[3]
    point_bg3_sig <- discrete_palette[3]
    #point_bg3_sig <- "black"
    ci_lwd <- 5
    ci_lty <- 1
    ci_col1 <- discrete_palette[1]
    ci_col2 <- discrete_palette[2]
    ci_col3 <- discrete_palette[3]
    lines_col1 <- discrete_palette[1]
    lines_col2 <- discrete_palette[2]
    lines_col3 <- discrete_palette[3]
    #ci_col1 <- "grey50"
    #ci_col2 <- "grey50"
    #ci_col3 <- "grey50"
    #lines_col1 <- "grey50"
    #lines_col2 <- "grey50"
    #lines_col3 <- "grey50"
    lines_lwd <- 1
    lines_lty <- 1
    y_axis_location <- "left"
    y_axis_tick <- FALSE
    y_axis_tick_length <- 0
    y_axis_label_cex <- 0.75
    label_distance <- 1.5
    label_col <- "black"
    label_cex <- 1
    if (equal_axis == TRUE && track_number == 2) {
        a <- min(track1_data[[lower_ci]])
        b <- min(track1_data[[upper_ci]])
        axis_min <- min(a, b)
        axis_min <- round(axis_min, 3)
        track1_axis_min <- round(axis_min + (axis_min * 0.1), 
            3)
        a <- max(track1_data[[lower_ci]])
        b <- max(track1_data[[upper_ci]])
        axis_max <- max(a, b)
        axis_max <- round(axis_max, 3)
        track1_axis_max <- round(axis_max + (axis_max * 0.01), 
            3)
        a <- min(track2_data[[lower_ci]])
        b <- min(track2_data[[upper_ci]])
        axis_min <- min(a, b)
        axis_min <- round(axis_min, 3)
        track2_axis_min <- round(axis_min + (axis_min * 0.1), 
            3)
        a <- max(track2_data[[lower_ci]])
        b <- max(track2_data[[upper_ci]])
        axis_max <- max(a, b)
        axis_max <- round(axis_max, 3)
        track2_axis_max <- round(axis_max + (axis_max * 0.01), 
            3)
        axis_min <- min(track1_axis_min, track2_axis_min)
        axis_max <- max(track1_axis_max, track2_axis_max)
        track1_axis_min <- axis_min
        track2_axis_min <- axis_min
        track1_axis_max <- axis_max
        track2_axis_max <- axis_max
    }
    else if (equal_axis == TRUE && track_number == 3) {
        a <- min(track1_data[[lower_ci]])
        b <- min(track1_data[[upper_ci]])
        axis_min <- min(a, b)
        axis_min <- round(axis_min, 3)
        track1_axis_min <- round(axis_min + (axis_min * 0.1), 
            3)
        a <- max(track1_data[[lower_ci]])
        b <- max(track1_data[[upper_ci]])
        axis_max <- max(a, b)
        axis_max <- round(axis_max, 3)
        track1_axis_max <- round(axis_max + (axis_max * 0.01), 
            3)
        a <- min(track2_data[[lower_ci]])
        b <- min(track2_data[[upper_ci]])
        axis_min <- min(a, b)
        axis_min <- round(axis_min, 3)
        track2_axis_min <- round(axis_min + (axis_min * 0.1), 
            3)
        a <- max(track2_data[[lower_ci]])
        b <- max(track2_data[[upper_ci]])
        axis_max <- max(a, b)
        axis_max <- round(axis_max, 3)
        track2_axis_max <- round(axis_max + (axis_max * 0.01), 
            3)
        a <- min(track3_data[[lower_ci]])
        b <- min(track3_data[[upper_ci]])
        axis_min <- min(a, b)
        axis_min <- round(axis_min, 3)
        track3_axis_min <- round(axis_min + (axis_min * 0.1), 
            3)
        a <- max(track3_data[[lower_ci]])
        b <- max(track3_data[[upper_ci]])
        axis_max <- max(a, b)
        axis_max <- round(axis_max, 3)
        track3_axis_max <- round(axis_max + (axis_max * 0.01), 
            3)
        axis_min <- min(track1_axis_min, track2_axis_min, track3_axis_min)
        axis_max <- max(track1_axis_max, track2_axis_max, track3_axis_min)
        track1_axis_min <- axis_min
        track2_axis_min <- axis_min
        track3_axis_min <- axis_min
        track1_axis_max <- axis_max
        track2_axis_max <- axis_max
        track3_axis_max <- axis_max
    }
    else if (equal_axis == FALSE && track_number == 2) {
        a <- min(track1_data[[lower_ci]])
        b <- min(track1_data[[upper_ci]])
        axis_min <- min(a, b)
        axis_min <- round(axis_min, 3)
        track1_axis_min <- round(axis_min + (axis_min * 0.1), 
            3)
        a <- max(track1_data[[lower_ci]])
        b <- max(track1_data[[upper_ci]])
        axis_max <- max(a, b)
        axis_max <- round(axis_max, 3)
        track1_axis_max <- round(axis_max + (axis_max * 0.01), 
            3)
        a <- min(track2_data[[lower_ci]])
        b <- min(track2_data[[upper_ci]])
        axis_min <- min(a, b)
        axis_min <- round(axis_min, 3)
        track2_axis_min <- round(axis_min + (axis_min * 0.1), 
            3)
        a <- max(track2_data[[lower_ci]])
        b <- max(track2_data[[upper_ci]])
        axis_max <- max(a, b)
        axis_max <- round(axis_max, 3)
        track2_axis_max <- round(axis_max + (axis_max * 0.01), 
            3)
    }
    else if (equal_axis == FALSE && track_number == 3) {
        a <- min(track1_data[[lower_ci]])
        b <- min(track1_data[[upper_ci]])
        axis_min <- min(a, b)
        axis_min <- round(axis_min, 3)
        track1_axis_min <- round(axis_min + (axis_min * 0.1), 
            3)
        a <- max(track1_data[[lower_ci]])
        b <- max(track1_data[[upper_ci]])
        axis_max <- max(a, b)
        axis_max <- round(axis_max, 3)
        track1_axis_max <- round(axis_max + (axis_max * 0.01), 
            3)
        a <- min(track2_data[[lower_ci]])
        b <- min(track2_data[[upper_ci]])
        axis_min <- min(a, b)
        axis_min <- round(axis_min, 3)
        track2_axis_min <- round(axis_min + (axis_min * 0.1), 
            3)
        a <- max(track2_data[[lower_ci]])
        b <- max(track2_data[[upper_ci]])
        axis_max <- max(a, b)
        axis_max <- round(axis_max, 3)
        track2_axis_max <- round(axis_max + (axis_max * 0.01), 
            3)
        a <- min(track3_data[[lower_ci]])
        b <- min(track3_data[[upper_ci]])
        axis_min <- min(a, b)
        axis_min <- round(axis_min, 3)
        track3_axis_min <- round(axis_min + (axis_min * 0.1), 
            3)
        a <- max(track3_data[[lower_ci]])
        b <- max(track3_data[[upper_ci]])
        axis_max <- max(a, b)
        axis_max <- round(axis_max, 3)
        track3_axis_max <- round(axis_max + (axis_max * 0.01), 
            3)
    }
    else if (equal_axis == FALSE && track_number == 1) {
        a <- min(track1_data[[lower_ci]])
        b <- min(track1_data[[upper_ci]])
        axis_min <- min(a, b)
        axis_min <- round(axis_min, 3)
        track1_axis_min <- round(axis_min + (axis_min * 0.1), 
            3)
        a <- max(track1_data[[lower_ci]])
        b <- max(track1_data[[upper_ci]])
        axis_max <- max(a, b)
        axis_max <- round(axis_max, 3)
        track1_axis_max <- round(axis_max + (axis_max * 0.01), 
            3)
    }
    data <- track1_data
    if (order == TRUE) {
        data <- data[order(data[[section_column]], data[[label_column]]), 
            ]
    }
    else if (order == FALSE) {
        data[[section_column]] <- stats::reorder(data[[section_column]], 
            data[[order_column]])
        data <- data[order(data[[section_column]], data[[label_column]]), 
            ]
    }
    data$x <- with(data, ave(seq_along(data[[section_column]]), 
        data[[section_column]], FUN = seq_along))
    npercat <- as.vector(table(data[[section_column]]))
    getaxis <- function(data) {
        for (i in 1:nrow(data)) {
            data$n[i] <- as.numeric(nrow(subset(data, data[[section_column]] == 
                data[[section_column]][i])))
            data$ncat[i] <- data$x[i]/data$n[i]
        }
        return(data)
    }
    data <- getaxis(data)
    data$section_numbers = factor(data[[section_column]], labels = 1:nlevels(data[[section_column]]))
    gap = c(rep(1, nlevels(data[[section_column]]) - 1), start_gap)
    circlize::circos.clear()
    graphics::par(mar = c(0.6, 0.5, 0.5, 0.5) * circle_size, 
        cex = 0.8, xpd = NA)
    circlize::circos.par(cell.padding = c(0, 0.5, 0, 0.5), start.degree = start_degree, 
        gap.degree = gap, track.margin = c(0.012, 0.012), points.overflow.warning = FALSE, 
        track.height = section_track_height, clock.wise = TRUE)
    circlize::circos.initialize(factors = data$section_numbers, 
        xlim = c(0, 1), sector.width = npercat)
    circlize::circos.trackPlotRegion(factors = data$section_numbers, 
        track.index = track1, x = data$ncat, ylim = c(0, 1), 
        track.height = 0.05, panel.fun = function(x, y) {
            chr = circlize::get.cell.meta.data("sector.index")
            xlim = circlize::get.cell.meta.data("xlim")
            ylim = circlize::get.cell.meta.data("ylim")
            circlize::circos.rect(xlim[1], 0, xlim[2], 1, border = NA, 
                col = section_fill_colour)
            circlize::circos.text(mean(xlim), mean(ylim), chr, 
                cex = 1, facing = "outside", niceFacing = TRUE, 
                col = section_text_colour)
        }, bg.border = NA)
    circlize::circos.trackText(factors = data$section_numbers, 
        track.index = track1, x = data$ncat, y = data$b * 0 + 
            label_distance, labels = data[[label_column]], facing = "reverse.clockwise", 
        niceFacing = TRUE, adj = c(1, 1), col = label_col, cex = label_cex)
    if (track_number >= 1 && track1_type == "points") {
        data <- track1_data
        track.index <- 2
        if (order == TRUE) {
            data <- data[order(data[[section_column]], data[[label_column]]), 
                ]
        }
        else if (order == FALSE) {
            data[[section_column]] <- stats::reorder(data[[section_column]], 
                data[[order_column]])
            data <- data[order(data[[section_column]], data[[label_column]]), 
                ]
        }
        data$x <- with(data, ave(seq_along(data[[section_column]]), 
            data[[section_column]], FUN = seq_along))
        npercat <- as.vector(table(data[[section_column]]))
        getaxis <- function(data) {
            for (i in 1:nrow(data)) {
                data$n[i] <- as.numeric(nrow(subset(data, data[[section_column]] == 
                  data[[section_column]][i])))
                data$ncat[i] <- data$x[i]/data$n[i]
            }
            return(data)
        }
        data <- getaxis(data)
        data$section_numbers = factor(data[[section_column]], 
            labels = 1:nlevels(data[[section_column]]))
        gap = c(rep(1, nlevels(data[[section_column]]) - 1), 
            start_gap)
        for (i in 1:nlevels(data$section_numbers)) {
            data1 = subset(data, section_numbers == i)
            circlize::circos.trackPlotRegion(factors = data1$section_numbers, 
                track.index = track.index, x = data1$ncat, y = data1[[estimate_column]], 
                ylim = c(track1_axis_min, track1_axis_max), track.height = track1_height, 
                bg.border = NA, bg.col = NA, panel.fun = function(x, 
                  y) {
                  circlize::circos.lines(x = x, y = y * 0 + track_axis_reference, 
                    col = reference_line_colour, lwd = reference_line_thickness, 
                    lty = reference_line_type)
                  ###
                  circlize::circos.lines(x = x, y = y * 0 + track_axis_reference + 0.1, 
                    col = guide_line_colour, lwd = guide_line_thickness, 
                    lty = reference_line_type)
                  ###
                  circlize::circos.lines(x = x, y = y * 0 + track_axis_reference+ 0.2, 
                    col = guide_line_colour, lwd = guide_line_thickness, 
                    lty = reference_line_type)
                  ###
                  circlize::circos.lines(x = x, y = y * 0 + track_axis_reference+ 0.3, 
                    col = guide_line_colour, lwd = guide_line_thickness, 
                    lty = reference_line_type)
                  ###
                  circlize::circos.lines(x = x, y = y * 0 + track_axis_reference+ 0.4, 
                   col = guide_line_colour, lwd = guide_line_thickness, 
                    lty = reference_line_type)
                  ##
                  circlize::circos.segments(x0 = data1$ncat, 
                    x1 = data1$ncat, y0 = data1[[estimate_column]] * 
                      0 - -(data1[[lower_ci]]), y1 = data1[[estimate_column]] * 
                      0 + data1[[upper_ci]], col = ci_col1, lwd = ci_lwd, 
                    lty = ci_lty, sector.index = i)
                })
        }
        circlize::circos.trackPoints(factors = subset(data, data[[pvalue_column]] > 
            pvalue_adjustment)$section_numbers, track.index = track.index, 
            x = subset(data, data[[pvalue_column]] > pvalue_adjustment)$ncat, 
            y = subset(data, data[[pvalue_column]] > pvalue_adjustment)[[estimate_column]], 
            cex = point_cex, pch = point_pch, col = point_col1, 
            bg = point_bg1)
        circlize::circos.trackPoints(factors = subset(data, data[[pvalue_column]] < 
            pvalue_adjustment)$section_numbers, track.index = track.index, 
            x = subset(data, data[[pvalue_column]] < pvalue_adjustment)$ncat, 
            y = subset(data, data[[pvalue_column]] < pvalue_adjustment)[[estimate_column]], 
            cex = point_cex, pch = point_pch, col = point_col1_sig, 
            bg = point_bg1_sig)
        circlize::circos.yaxis(side = y_axis_location, sector.index = x_axis_index, 
            track.index = track.index, at = c(track1_axis_min, 
                track_axis_reference, track1_axis_max), tick = y_axis_tick, 
            tick.length = y_axis_tick_length, labels.cex = y_axis_label_cex)
    }
    if (track_number >= 1 && track1_type == "lines") {
        data <- track1_data
        track.index <- 2
        if (order == TRUE) {
            data <- data[order(data[[section_column]], data[[label_column]]), 
                ]
        }
        else if (order == FALSE) {
            data[[section_column]] <- stats::reorder(data[[section_column]], 
                data[[order_column]])
            data <- data[order(data[[section_column]], data[[label_column]]), 
                ]
        }
        data$x <- with(data, ave(seq_along(data[[section_column]]), 
            data[[section_column]], FUN = seq_along))
        npercat <- as.vector(table(data[[section_column]]))
        getaxis <- function(data) {
            for (i in 1:nrow(data)) {
                data$n[i] <- as.numeric(nrow(subset(data, data[[section_column]] == 
                  data[[section_column]][i])))
                data$ncat[i] <- data$x[i]/data$n[i]
            }
            return(data)
        }
        data <- getaxis(data)
        data$section_numbers = factor(data[[section_column]], 
            labels = 1:nlevels(data[[section_column]]))
        gap = c(rep(1, nlevels(data[[section_column]]) - 1), 
            start_gap)
        axis_min <- round(min(data[[lines_column]]), 3)
        axis_min <- round(axis_min + (axis_min * 0.1), 3)
        axis_min_half <- round(axis_min/2, 3)
        axis_max <- round(max(data[[lines_column]]), 3)
        axis_max <- round(axis_max + (axis_max * 0.1), 3)
        axis_max_half <- round(axis_max/2, 3)
        circlize::circos.trackPlotRegion(factors = data$section_numbers, 
            track.index = track.index, x = data$ncat, y = data[[lines_column]], 
            ylim = c(axis_min, axis_max), track.height = track1_height, 
            bg.border = NA, bg.col = NA, panel.fun = function(x, 
                y) {
                circlize::circos.lines(x = x, y = y, col = lines_col1, 
                  lwd = lines_lwd, lty = lines_lty, type = lines_type)
            })
        circlize::circos.yaxis(side = y_axis_location, sector.index = x_axis_index, 
            track.index = track.index, at = c(axis_min, track_axis_reference, 
                axis_max), tick = y_axis_tick, tick.length = y_axis_tick_length, 
            labels.cex = y_axis_label_cex)
    }
    if (track_number >= 1 && track1_type == "bar") {
        data <- track1_data
        track.index <- 2
        if (order == TRUE) {
            data <- data[order(data[[section_column]], data[[label_column]]), 
                ]
        }
        else if (order == FALSE) {
            data[[section_column]] <- stats::reorder(data[[section_column]], 
                data[[order_column]])
            data <- data[order(data[[section_column]], data[[label_column]]), 
                ]
        }
        data$x <- with(data, ave(seq_along(data[[section_column]]), 
            data[[section_column]], FUN = seq_along))
        npercat <- as.vector(table(data[[section_column]]))
        getaxis <- function(data) {
            for (i in 1:nrow(data)) {
                data$n[i] <- as.numeric(nrow(subset(data, data[[section_column]] == 
                  data[[section_column]][i])))
                data$ncat[i] <- data$x[i]/data$n[i]
            }
            return(data)
        }
        data <- getaxis(data)
        data$section_numbers = factor(data[[section_column]], 
            labels = 1:nlevels(data[[section_column]]))
        gap = c(rep(1, nlevels(data[[section_column]]) - 1), 
            start_gap)
        axis_min <- round(min(data[[bar_column]]), 3)
        axis_max <- round(max(data[[bar_column]]), 3)
        circlize::circos.trackPlotRegion(factors = data$section_numbers, 
            track.index = track.index, x = data$ncat, y = data[[bar_column]], 
            ylim = c(axis_min, axis_max), track.height = track1_height, 
            bg.border = NA, bg.col = NA, panel.fun = function(x, 
                y) {
                circlize::circos.lines(x = x, y = y, col = lines_col1, 
                  lwd = 1, lty = lines_lty, type = "s", area = T, 
                  border = "White", baseline = "bottom")
            })
        circlize::circos.yaxis(side = y_axis_location, sector.index = x_axis_index, 
            track.index = track.index, at = c(axis_min, track_axis_reference, 
                axis_max), tick = y_axis_tick, tick.length = y_axis_tick_length, 
            labels.cex = y_axis_label_cex)
    }
    if (track_number >= 1 && track1_type == "histogram") {
        data <- track1_data
        track.index <- 2
        if (order == TRUE) {
            data <- data[order(data[[section_column]], data[[label_column]]), 
                ]
        }
        else if (order == FALSE) {
            data[[section_column]] <- stats::reorder(data[[section_column]], 
                data[[order_column]])
            data <- data[order(data[[section_column]], data[[label_column]]), 
                ]
        }
        data$x <- with(data, ave(seq_along(data[[section_column]]), 
            data[[section_column]], FUN = seq_along))
        npercat <- as.vector(table(data[[section_column]]))
        getaxis <- function(data) {
            for (i in 1:nrow(data)) {
                data$n[i] <- as.numeric(nrow(subset(data, data[[section_column]] == 
                  data[[section_column]][i])))
                data$ncat[i] <- data$x[i]/data$n[i]
            }
            return(data)
        }
        data <- getaxis(data)
        data$section_numbers = factor(data[[section_column]], 
            labels = 1:nlevels(data[[section_column]]))
        gap = c(rep(1, nlevels(data[[section_column]]) - 1), 
            start_gap)
        axis_min <- round(min(data[[histogram_column]]), 3)
        axis_min <- round(axis_min + (axis_min * 0.1), 3)
        axis_min_half <- round(axis_min/2, 3)
        axis_max <- round(max(data[[histogram_column]]), 3)
        axis_max <- round(axis_max + (axis_max * 0.1), 3)
        axis_max_half <- round(axis_max/2, 3)
        circlize::circos.trackHist(factors = data$section_numbers, 
            x = data[[histogram_column]], track.height = track1_height, 
            track.index = NULL, col = lines_col1, border = lines_col1, 
            bg.border = NA, draw.density = histogram_densityplot, 
            bin.size = histogram_binsize)
        circlize::circos.yaxis(side = y_axis_location, sector.index = x_axis_index, 
            track.index = track.index, at = c(axis_min, track_axis_reference, 
                axis_max), tick = y_axis_tick, tick.length = y_axis_tick_length, 
            labels.cex = y_axis_label_cex)
    }
    if (track_number >= 2 && track2_type == "points") {
        data <- track2_data
        track.index <- 3
        if (order == TRUE) {
            data <- data[order(data[[section_column]], data[[label_column]]), 
                ]
        }
        else if (order == FALSE) {
            data[[section_column]] <- stats::reorder(data[[section_column]], 
                data[[order_column]])
            data <- data[order(data[[section_column]], data[[label_column]]), 
                ]
        }
        data$x <- with(data, ave(seq_along(data[[section_column]]), 
            data[[section_column]], FUN = seq_along))
        npercat <- as.vector(table(data[[section_column]]))
        getaxis <- function(data) {
            for (i in 1:nrow(data)) {
                data$n[i] <- as.numeric(nrow(subset(data, data[[section_column]] == 
                  data[[section_column]][i])))
                data$ncat[i] <- data$x[i]/data$n[i]
            }
            return(data)
        }
        data <- getaxis(data)
        data$section_numbers = factor(data[[section_column]], 
            labels = 1:nlevels(data[[section_column]]))
        gap = c(rep(1, nlevels(data[[section_column]]) - 1), 
            start_gap)
        for (i in 1:nlevels(data$section_numbers)) {
            data1 = subset(data, section_numbers == i)
            circlize::circos.trackPlotRegion(factors = data1$section_numbers, 
                track.index = track.index, x = data1$ncat, y = data1[[estimate_column]], 
                ylim = c(track2_axis_min, track2_axis_max), track.height = track2_height, 
                bg.border = NA, bg.col = NA, panel.fun = function(x, 
                  y) {
                  circlize::circos.lines(x = x, y = y * 0 + track_axis_reference, 
                    col = reference_line_colour, lwd = reference_line_thickness, 
                    lty = reference_line_type)
                  ###
                  circlize::circos.lines(x = x, y = y * 0 + track_axis_reference + 0.1, 
                    col = guide_line_colour, lwd = guide_line_thickness, 
                    lty = reference_line_type)
                  ###
                  circlize::circos.lines(x = x, y = y * 0 + track_axis_reference+ 0.2, 
                    col = guide_line_colour, lwd = guide_line_thickness, 
                    lty = reference_line_type)
                  ###
                  circlize::circos.lines(x = x, y = y * 0 + track_axis_reference+ 0.3, 
                    col = guide_line_colour, lwd = guide_line_thickness, 
                    lty = reference_line_type)
                  ###
                  circlize::circos.lines(x = x, y = y * 0 + track_axis_reference+ 0.4, 
                   col = guide_line_colour, lwd = guide_line_thickness, 
                    lty = reference_line_type)
                  ##
                  circlize::circos.segments(x0 = data1$ncat, 
                    x1 = data1$ncat, y0 = data1[[estimate_column]] * 
                      0 - -(data1[[lower_ci]]), y1 = data1[[estimate_column]] * 
                      0 + data1[[upper_ci]], col = ci_col2, lwd = ci_lwd, 
                    lty = ci_lty, sector.index = i)
                })
        }
        circlize::circos.trackPoints(factors = subset(data, data[[pvalue_column]] > 
            pvalue_adjustment)$section_numbers, track.index = track.index, 
            x = subset(data, data[[pvalue_column]] > pvalue_adjustment)$ncat, 
            y = subset(data, data[[pvalue_column]] > pvalue_adjustment)[[estimate_column]], 
            cex = point_cex, pch = point_pch, col = point_col2, 
            bg = point_bg2)
        circlize::circos.trackPoints(factors = subset(data, data[[pvalue_column]] < 
            pvalue_adjustment)$section_numbers, track.index = track.index, 
            x = subset(data, data[[pvalue_column]] < pvalue_adjustment)$ncat, 
            y = subset(data, data[[pvalue_column]] < pvalue_adjustment)[[estimate_column]], 
            cex = point_cex, pch = point_pch, col = point_col2_sig, 
            bg = point_bg2_sig)
        circlize::circos.yaxis(side = y_axis_location, sector.index = x_axis_index, 
            track.index = track.index, at = c(track2_axis_min, 
                track_axis_reference, track2_axis_max), tick = y_axis_tick, 
            tick.length = y_axis_tick_length, labels.cex = y_axis_label_cex)
    }
    if (track_number >= 2 && track2_type == "lines") {
        data <- track2_data
        track.index <- 3
        if (order == TRUE) {
            data <- data[order(data[[section_column]], data[[label_column]]), 
                ]
        }
        else if (order == FALSE) {
            data[[section_column]] <- stats::reorder(data[[section_column]], 
                data[[order_column]])
            data <- data[order(data[[section_column]], data[[label_column]]), 
                ]
        }
        data$x <- with(data, ave(seq_along(data[[section_column]]), 
            data[[section_column]], FUN = seq_along))
        npercat <- as.vector(table(data[[section_column]]))
        getaxis <- function(data) {
            for (i in 1:nrow(data)) {
                data$n[i] <- as.numeric(nrow(subset(data, data[[section_column]] == 
                  data[[section_column]][i])))
                data$ncat[i] <- data$x[i]/data$n[i]
            }
            return(data)
        }
        data <- getaxis(data)
        data$section_numbers = factor(data[[section_column]], 
            labels = 1:nlevels(data[[section_column]]))
        gap = c(rep(1, nlevels(data[[section_column]]) - 1), 
            start_gap)
        axis_min <- round(min(data[[lines_column]]), 3)
        axis_min <- round(axis_min + (axis_min * 0.1), 3)
        axis_min_half <- round(axis_min/2, 3)
        axis_max <- round(max(data[[lines_column]]), 3)
        axis_max <- round(axis_max + (axis_max * 0.1), 3)
        axis_max_half <- round(axis_max/2, 3)
        circlize::circos.trackPlotRegion(factors = data$section_numbers, 
            track.index = track.index, x = data$ncat, y = data[[lines_column]], 
            ylim = c(axis_min, axis_max), track.height = track2_height, 
            bg.border = NA, bg.col = NA, panel.fun = function(x, 
                y) {
                circlize::circos.lines(x = x, y = y, col = lines_col2, 
                  lwd = lines_lwd, lty = lines_lty, type = lines_type)
            })
        circlize::circos.yaxis(side = y_axis_location, sector.index = x_axis_index, 
            track.index = track.index, at = c(axis_min, track_axis_reference, 
                axis_max), tick = y_axis_tick, tick.length = y_axis_tick_length, 
            labels.cex = y_axis_label_cex)
    }
    if (track_number >= 2 && track2_type == "bar") {
        data <- track1_data
        track.index <- 3
        if (order == TRUE) {
            data <- data[order(data[[section_column]], data[[label_column]]), 
                ]
        }
        else if (order == FALSE) {
            data[[section_column]] <- stats::reorder(data[[section_column]], 
                data[[order_column]])
            data <- data[order(data[[section_column]], data[[label_column]]), 
                ]
        }
        data$x <- with(data, ave(seq_along(data[[section_column]]), 
            data[[section_column]], FUN = seq_along))
        npercat <- as.vector(table(data[[section_column]]))
        getaxis <- function(data) {
            for (i in 1:nrow(data)) {
                data$n[i] <- as.numeric(nrow(subset(data, data[[section_column]] == 
                  data[[section_column]][i])))
                data$ncat[i] <- data$x[i]/data$n[i]
            }
            return(data)
        }
        data <- getaxis(data)
        data$section_numbers = factor(data[[section_column]], 
            labels = 1:nlevels(data[[section_column]]))
        gap = c(rep(1, nlevels(data[[section_column]]) - 1), 
            start_gap)
        axis_min <- round(min(data[[bar_column]]), 3)
        axis_max <- round(max(data[[bar_column]]), 3)
        circlize::circos.trackPlotRegion(factors = data$section_numbers, 
            track.index = track.index, x = data$ncat, y = data[[bar_column]], 
            ylim = c(axis_min, axis_max), track.height = track2_height, 
            bg.border = NA, bg.col = NA, panel.fun = function(x, 
                y) {
                circlize::circos.lines(x = x, y = y, col = lines_col2, 
                  lwd = 1, lty = lines_lty, type = "s", area = T, 
                  border = "White", baseline = "bottom")
            })
        circlize::circos.yaxis(side = y_axis_location, sector.index = x_axis_index, 
            track.index = track.index, at = c(axis_min, track_axis_reference, 
                axis_max), tick = y_axis_tick, tick.length = y_axis_tick_length, 
            labels.cex = y_axis_label_cex)
    }
    if (track_number >= 2 && track2_type == "histogram") {
        data <- track2_data
        track.index <- 3
        if (order == TRUE) {
            data <- data[order(data[[section_column]], data[[label_column]]), 
                ]
        }
        else if (order == FALSE) {
            data[[section_column]] <- stats::reorder(data[[section_column]], 
                data[[order_column]])
            data <- data[order(data[[section_column]], data[[label_column]]), 
                ]
        }
        data$x <- with(data, ave(seq_along(data[[section_column]]), 
            data[[section_column]], FUN = seq_along))
        npercat <- as.vector(table(data[[section_column]]))
        getaxis <- function(data) {
            for (i in 1:nrow(data)) {
                data$n[i] <- as.numeric(nrow(subset(data, data[[section_column]] == 
                  data[[section_column]][i])))
                data$ncat[i] <- data$x[i]/data$n[i]
            }
            return(data)
        }
        data <- getaxis(data)
        data$section_numbers = factor(data[[section_column]], 
            labels = 1:nlevels(data[[section_column]]))
        gap = c(rep(1, nlevels(data[[section_column]]) - 1), 
            start_gap)
        axis_min <- round(min(data[[histogram_column]]), 3)
        axis_min <- round(axis_min + (axis_min * 0.1), 3)
        axis_min_half <- round(axis_min/2, 3)
        axis_max <- round(max(data[[histogram_column]]), 3)
        axis_max <- round(axis_max + (axis_max * 0.1), 3)
        axis_max_half <- round(axis_max/2, 3)
        circlize::circos.trackHist(factors = data$section_numbers, 
            x = data[[histogram_column]], track.height = track2_height, 
            track.index = NULL, col = lines_col2, border = lines_col2, 
            bg.border = NA, draw.density = histogram_densityplot, 
            bin.size = histogram_binsize)
        circlize::circos.yaxis(side = y_axis_location, sector.index = x_axis_index, 
            track.index = track.index, at = c(axis_min, track_axis_reference, 
                axis_max), tick = y_axis_tick, tick.length = y_axis_tick_length, 
            labels.cex = y_axis_label_cex)
    }
    if (track_number >= 3 && track3_type == "points") {
        data <- track3_data
        track.index <- 4
        if (order == TRUE) {
            data <- data[order(data[[section_column]], data[[label_column]]), 
                ]
        }
        else if (order == FALSE) {
            data[[section_column]] <- stats::reorder(data[[section_column]], 
                data[[order_column]])
            data <- data[order(data[[section_column]], data[[label_column]]), 
                ]
        }
        data$x <- with(data, ave(seq_along(data[[section_column]]), 
            data[[section_column]], FUN = seq_along))
        npercat <- as.vector(table(data[[section_column]]))
        getaxis <- function(data) {
            for (i in 1:nrow(data)) {
                data$n[i] <- as.numeric(nrow(subset(data, data[[section_column]] == 
                  data[[section_column]][i])))
                data$ncat[i] <- data$x[i]/data$n[i]
            }
            return(data)
        }
        data <- getaxis(data)
        data$section_numbers = factor(data[[section_column]], 
            labels = 1:nlevels(data[[section_column]]))
        gap = c(rep(1, nlevels(data[[section_column]]) - 1), 
            start_gap)
        for (i in 1:nlevels(data$section_numbers)) {
            data1 = subset(data, section_numbers == i)
            circlize::circos.trackPlotRegion(factors = data1$section_numbers, 
                track.index = track.index, x = data1$ncat, y = data1[[estimate_column]], 
                ylim = c(track3_axis_min, track3_axis_max), track.height = track3_height, 
                bg.border = NA, bg.col = NA, panel.fun = function(x, 
                  y) {
                  circlize::circos.lines(x = x, y = y * 0 + track_axis_reference, 
                    col = reference_line_colour, lwd = reference_line_thickness, 
                    lty = reference_line_type)
                  ###
                  circlize::circos.lines(x = x, y = y * 0 + track_axis_reference + 0.1, 
                    col = guide_line_colour, lwd = guide_line_thickness, 
                    lty = reference_line_type)
                  ###
                  circlize::circos.lines(x = x, y = y * 0 + track_axis_reference+ 0.2, 
                    col = guide_line_colour, lwd = guide_line_thickness, 
                    lty = reference_line_type)
                  ###
                  circlize::circos.lines(x = x, y = y * 0 + track_axis_reference+ 0.3, 
                    col = guide_line_colour, lwd = guide_line_thickness, 
                    lty = reference_line_type)
                  ###
                  circlize::circos.lines(x = x, y = y * 0 + track_axis_reference+ 0.4, 
                   col = guide_line_colour, lwd = guide_line_thickness, 
                    lty = reference_line_type)
                  ##
                  circlize::circos.segments(x0 = data1$ncat, 
                    x1 = data1$ncat, y0 = data1[[estimate_column]] * 
                      0 - -(data1[[lower_ci]]), y1 = data1[[estimate_column]] * 
                      0 + data1[[upper_ci]], col = ci_col3, lwd = ci_lwd, 
                    lty = ci_lty, sector.index = i)
                })
        }
        circlize::circos.trackPoints(factors = subset(data, data[[pvalue_column]] > 
            pvalue_adjustment)$section_numbers, track.index = track.index, 
            x = subset(data, data[[pvalue_column]] > pvalue_adjustment)$ncat, 
            y = subset(data, data[[pvalue_column]] > pvalue_adjustment)[[estimate_column]], 
            cex = point_cex, pch = point_pch, col = point_col3, 
            bg = point_bg3)
        circlize::circos.trackPoints(factors = subset(data, data[[pvalue_column]] < 
            pvalue_adjustment)$section_numbers, track.index = track.index, 
            x = subset(data, data[[pvalue_column]] < pvalue_adjustment)$ncat, 
            y = subset(data, data[[pvalue_column]] < pvalue_adjustment)[[estimate_column]], 
            cex = point_cex, pch = point_pch, col = point_col3_sig, 
            bg = point_bg3_sig)
        circlize::circos.yaxis(side = y_axis_location, sector.index = x_axis_index, 
            track.index = track.index, at = c(track3_axis_min, 
                track_axis_reference, track3_axis_max), tick = y_axis_tick, 
            tick.length = y_axis_tick_length, labels.cex = y_axis_label_cex)
    }
    if (track_number >= 3 && track3_type == "lines") {
        data <- track3_data
        track.index <- 4
        if (order == TRUE) {
            data <- data[order(data[[section_column]], data[[label_column]]), 
                ]
        }
        else if (order == FALSE) {
            data[[section_column]] <- stats::reorder(data[[section_column]], 
                data[[order_column]])
            data <- data[order(data[[section_column]], data[[label_column]]), 
                ]
        }
        data$x <- with(data, ave(seq_along(data[[section_column]]), 
            data[[section_column]], FUN = seq_along))
        npercat <- as.vector(table(data[[section_column]]))
        getaxis <- function(data) {
            for (i in 1:nrow(data)) {
                data$n[i] <- as.numeric(nrow(subset(data, data[[section_column]] == 
                  data[[section_column]][i])))
                data$ncat[i] <- data$x[i]/data$n[i]
            }
            return(data)
        }
        data <- getaxis(data)
        data$section_numbers = factor(data[[section_column]], 
            labels = 1:nlevels(data[[section_column]]))
        gap = c(rep(1, nlevels(data[[section_column]]) - 1), 
            start_gap)
        axis_min <- round(min(data[[lines_column]]), 3)
        axis_min <- round(axis_min + (axis_min * 0.1), 3)
        axis_min_half <- round(axis_min/2, 3)
        axis_max <- round(max(data[[lines_column]]), 3)
        axis_max <- round(axis_max + (axis_max * 0.1), 3)
        axis_max_half <- round(axis_max/2, 3)
        circlize::circos.trackPlotRegion(factors = data$section_numbers, 
            track.index = track.index, x = data$ncat, y = data[[lines_column]], 
            ylim = c(axis_min, axis_max), track.height = track3_height, 
            bg.border = NA, bg.col = NA, panel.fun = function(x, 
                y) {
                circlize::circos.lines(x = x, y = y, col = lines_col3, 
                  lwd = lines_lwd, lty = lines_lty, type = lines_type)
            })
        circlize::circos.yaxis(side = y_axis_location, sector.index = x_axis_index, 
            track.index = track.index, at = c(axis_min, track_axis_reference, 
                axis_max), tick = y_axis_tick, tick.length = y_axis_tick_length, 
            labels.cex = y_axis_label_cex)
    }
    if (track_number >= 3 && track3_type == "bar") {
        data <- track1_data
        track.index <- 4
        if (order == TRUE) {
            data <- data[order(data[[section_column]], data[[label_column]]), 
                ]
        }
        else if (order == FALSE) {
            data[[section_column]] <- stats::reorder(data[[section_column]], 
                data[[order_column]])
            data <- data[order(data[[section_column]], data[[label_column]]), 
                ]
        }
        data$x <- with(data, ave(seq_along(data[[section_column]]), 
            data[[section_column]], FUN = seq_along))
        npercat <- as.vector(table(data[[section_column]]))
        getaxis <- function(data) {
            for (i in 1:nrow(data)) {
                data$n[i] <- as.numeric(nrow(subset(data, data[[section_column]] == 
                  data[[section_column]][i])))
                data$ncat[i] <- data$x[i]/data$n[i]
            }
            return(data)
        }
        data <- getaxis(data)
        data$section_numbers = factor(data[[section_column]], 
            labels = 1:nlevels(data[[section_column]]))
        gap = c(rep(1, nlevels(data[[section_column]]) - 1), 
            start_gap)
        axis_min <- round(min(data[[bar_column]]), 3)
        axis_max <- round(max(data[[bar_column]]), 3)
        circlize::circos.trackPlotRegion(factors = data$section_numbers, 
            track.index = track.index, x = data$ncat, y = data[[bar_column]], 
            ylim = c(axis_min, axis_max), track.height = track3_height, 
            bg.border = NA, bg.col = NA, panel.fun = function(x, 
                y) {
                circlize::circos.lines(x = x, y = y, col = lines_col3, 
                  lwd = 1, lty = lines_lty, type = "s", area = T, 
                  border = "White", baseline = "bottom")
            })
        circlize::circos.yaxis(side = y_axis_location, sector.index = x_axis_index, 
            track.index = track.index, at = c(axis_min, track_axis_reference, 
                axis_max), tick = y_axis_tick, tick.length = y_axis_tick_length, 
            labels.cex = y_axis_label_cex)
    }
    if (track_number >= 3 && track3_type == "histogram") {
        data <- track3_data
        track.index <- 4
        if (order == TRUE) {
            data <- data[order(data[[section_column]], data[[label_column]]), 
                ]
        }
        else if (order == FALSE) {
            data[[section_column]] <- stats::reorder(data[[section_column]], 
                data[[order_column]])
            data <- data[order(data[[section_column]], data[[label_column]]), 
                ]
        }
        data$x <- with(data, ave(seq_along(data[[section_column]]), 
            data[[section_column]], FUN = seq_along))
        npercat <- as.vector(table(data[[section_column]]))
        getaxis <- function(data) {
            for (i in 1:nrow(data)) {
                data$n[i] <- as.numeric(nrow(subset(data, data[[section_column]] == 
                  data[[section_column]][i])))
                data$ncat[i] <- data$x[i]/data$n[i]
            }
            return(data)
        }
        data <- getaxis(data)
        data$section_numbers = factor(data[[section_column]], 
            labels = 1:nlevels(data[[section_column]]))
        gap = c(rep(1, nlevels(data[[section_column]]) - 1), 
            start_gap)
        axis_min <- round(min(data[[histogram_column]]), 3)
        axis_min <- round(axis_min + (axis_min * 0.1), 3)
        axis_min_half <- round(axis_min/2, 3)
        axis_max <- round(max(data[[histogram_column]]), 3)
        axis_max <- round(axis_max + (axis_max * 0.1), 3)
        axis_max_half <- round(axis_max/2, 3)
        circlize::circos.trackHist(factors = data$section_numbers, 
            x = data[[histogram_column]], track.height = track3_height, 
            track.index = NULL, col = lines_col3, border = lines_col3, 
            bg.border = NA, draw.density = histogram_densityplot, 
            bin.size = histogram_binsize)
        circlize::circos.yaxis(side = y_axis_location, sector.index = x_axis_index, 
            track.index = track.index, at = c(axis_min, track_axis_reference, 
                axis_max), tick = y_axis_tick, tick.length = y_axis_tick_length, 
            labels.cex = y_axis_label_cex)
    }
    if (legend == TRUE && track_number == 1) {
        legend1 <- ComplexHeatmap::Legend(at = c(track1_label), 
            labels_gp = grid::gpar(fontsize = 15), ncol = 1, 
            border = NA, background = NA, legend_gp = grid::gpar(col = c(discrete_palette[1])), 
            type = "points", pch = 19, size = grid::unit(15, 
                "mm"), grid_height = grid::unit(15, "mm"), grid_width = grid::unit(15, 
                "mm"), direction = "vertical")
    }
    if (legend == TRUE && track_number >= 2) {
        legend1 <- ComplexHeatmap::Legend(at = c(track1_label, 
            track2_label), labels_gp = grid::gpar(fontsize = 15), 
            ncol = 1, border = NA, background = NA, legend_gp = grid::gpar(col = c(discrete_palette[1], 
                discrete_palette[2])), type = "points", pch = 19, 
            size = grid::unit(15, "mm"), grid_height = grid::unit(15, 
                "mm"), grid_width = grid::unit(15, "mm"), direction = "vertical")
    }
    if (legend == TRUE && track_number >= 3) {
        legend1 <- ComplexHeatmap::Legend(at = c(track1_label, 
            track2_label, track3_label), labels_gp = grid::gpar(fontsize = 15), 
            ncol = 1, border = NA, background = NA, legend_gp = grid::gpar(col = c(discrete_palette[1], 
                discrete_palette[2], discrete_palette[3])), type = "points", 
            pch = 19, size = grid::unit(15, "mm"), grid_height = grid::unit(15, 
                "mm"), grid_width = grid::unit(15, "mm"), direction = "vertical")
    }
    if (legend == TRUE) {
        legend2 <- ComplexHeatmap::Legend(at = pvalue_label, 
            labels_gp = grid::gpar(fontsize = 15), ncol = 1, 
            border = NA, background = NA, legend_gp = grid::gpar(col = c("black")), 
            type = "points", pch = 1, size = grid::unit(15, "mm"), 
            grid_height = grid::unit(15, "mm"), grid_width = grid::unit(15, 
                "mm"), direction = "vertical")
        names <- levels(as.factor(data[[section_column]]))
        names <- paste(1:nlevels(data[[section_column]]), names, 
            sep = ". ")
        legend3 <- ComplexHeatmap::Legend(at = names, labels_gp = grid::gpar(fontsize = 15), 
            nrow = 4, ncol = 7, border = NA, background = NA, 
            legend_gp = grid::gpar(col = c("black")), size = grid::unit(15, 
                "mm"), grid_height = grid::unit(15, "mm"), grid_width = grid::unit(10, 
                "mm"), direction = "horizontal")
        legend4 <- ComplexHeatmap::packLegend(legend1, legend2, 
            direction = "vertical", gap = grid::unit(0, "mm"))
        legend <- ComplexHeatmap::packLegend(legend4, legend3, 
            direction = "horizontal", gap = grid::unit(0, "mm"))
        legend_height <- legend@grob[["vp"]][["height"]]
        legend_width <- legend@grob[["vp"]][["width"]]
        grid::pushViewport(grid::viewport(x = grid::unit(0.5, 
            "npc"), y = grid::unit(0.08, "npc"), width = legend_width, 
            height = legend_height, just = c("center", "top")))
        grid::grid.draw(legend)
        grid::upViewport()
    }
}

