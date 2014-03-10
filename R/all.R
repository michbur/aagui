order_boxplot <- function(bp, bp_order) {
  len <- length(bp[[2]])
  
  bp2 <- bp
  bp2[[1]] <- bp[[1]][,bp_order]
  bp2[[3]] <- bp[[3]][,bp_order]
  bp2[[2]] <- bp[[2]][bp_order]
  bp2[[6]] <- bp[[6]][bp_order]
  
  out_ind <- lapply(bp_order, function(x) which(x == bp[[5]]))  
  bp2[[4]] <- bp[[4]][unlist(out_ind)]
  bp2[[5]] <- unlist(lapply(1L:len, function(x) rep(x, length(out_ind[[x]]))))
  
  bp2  
}


order_aggr <- function(central_tendency, group_order, groups) {
  res <- data.frame(t(sapply(split(central_tendency, groups), 
                             function(x) c(mean(x), sd(x), median(x), mad(x))))[group_order,])
  names(res) <- c("mean", "sd", "median", "mad")
  res
}

#group_order - order of declared 'groups'
aggr_data <- function(central_tendency, group_order, groups) {
  if (length(groups) != length(central_tendency))
    stop("'central tendency' and 'groups' must have equal length.", 
         call. = TRUE, domain = NA)    
  
  ord <- order_aggr(central_tendency, group_order, groups)#aggregate and order results  
  
  counts <- data.frame(table(as.vector(groups)))[,2][group_order]
  cbind(ord, counts)
}

test_data <- function(central_tendency, group_order, groups, test, 
                      p_val, alternative) {
  if (p_val >= 1 || p_val <= 0)
    stop("'p_val' must be a value between 0 and 1.", call. = TRUE, domain = NA)
  
  test <- tolower(test)
  if (grepl(test, "wilcoxon")) 
    test <- "wilcoxon"
  
  significance <- switch(test,
                         wilcoxon =  sapply(group_order, function(i) 
                           wilcox.test(central_tendency, central_tendency[groups == i], 
                                       alternative = alternative)$p.value),
                         t = sapply(group_order, function(i) 
                           t.test(central_tendency, central_tendency[groups == i], 
                                  alternative = alternative)$p.value))
  if (is.null(significance))
    stop("Wrong statistical test invoked. Currently only wilcox.test 
         and t.test are implemented.", call. = TRUE, domain = NA)
  significance_pos <- sapply(significance, function(x) x <= p_val)  
  list(significance = significance, significance_pos = significance_pos)
}

calc_freqs <- function(binary, labels) {
  rs <- rowSums(apply(binary, 1, function(x) table(labels[x == 1])))
  s <- sum(rs)
  res <- rs/s
  res
}

setClass("conplot_data", representation(central_tendency = "numeric", 
                                        group_order = "integer", 
                                        groups = "integer", 
                                        aggr = "data.frame", 
                                        freqs = "matrix",
                                        significance = "numeric", 
                                        significance_pos = "logical", 
                                        mid = "integer"))



calc_plots <- function(central_tendency, group_order, groups, binary_data, 
                       binary_labels, mean = TRUE, test = "wilcoxon", p_val = 0.01, 
                       alternative = "less") {

  mid <- 1L #moment index
  if (mean == FALSE) 
    mid <- 3L
  aggr <- aggr_data(central_tendency, group_order, groups)
  test <- test_data(central_tendency, group_order, groups, test, p_val, alternative)
  freqs <- sapply(group_order, function(y) calc_freqs(binary_data[y == groups,], 
                                                      binary_labels))
  dimnames(freqs)[[2]] <- group_order 
  
  conplot <- new("conplot_data")
  slot(conplot, "central_tendency") <- central_tendency
  slot(conplot, "group_order") <- group_order
  slot(conplot, "groups") <- groups
  slot(conplot, "aggr") <- aggr
  slot(conplot, "freqs") <- freqs
  slot(conplot, "significance") <- test[["significance"]]
  slot(conplot, "significance_pos") <- test[["significance_pos"]]
  slot(conplot, "mid") <- mid
  conplot
}


#concrete count axis
conaxis <- function(conplot_data, side = 3, col = "red") {
  len <- 1L:length(conplot_data@aggr[["counts"]])
  axis(side, at=len, labels=conplot_data@aggr[["counts"]], cex.axis = 0.5, 
       tcl = -0.3, mgp = c(3, 0.35, 0))
  axis(side, at=len[conplot_data@significance_pos], 
       labels=conplot_data@aggr[["counts"]][conplot_data@significance_pos], 
       cex.axis = 0.5, tcl = -0.3, lwd.ticks=2, col.ticks = col, 
       col.axis = col, mgp = c(3, 0.35, 0))
  mtext("Counts", side, cex = 0.5, line = 0.8, mgp = c(3, 0.35, 0)) 
}


conbox <- function(conplot_data, plot = TRUE, points = TRUE, line = TRUE, ...) {
  
  NOC <- max(conplot_data@group_order)
  bp <- boxplot(conplot_data@central_tendency ~ conplot_data@groups, plot = F)
  bp <- order_boxplot(bp, conplot_data@group_order)
  pts <- matrix(c(1L:NOC, conplot_data@aggr[["mean"]]), ncol = 2)
  coords <- list(boxplot = bp, points = pts)
  
  if (plot) {
    plot(NA, NA, xlim = c(1L,NOC), ylim = c(min(conplot_data@central_tendency), 
                                            max(conplot_data@central_tendency)), xaxt = "n", ylab = "", xlab = "")
    bxp(bp, boxwex = 0.035*NOC + 0.01, 
        xaxt = "n", boxfill = NA, add = T, ...)
    conaxis(conplot_data)
    
    if (points) #points represent mean
      points(pts, pch = 19, col = 3, cex = 1.5)
    
    if (line) #central tendency of all central tendencies
      if (conplot_data@mid == 1) {
        abline(h = mean(conplot_data@central_tendency), col = "red")
      } else abline(h = median(conplot_data@central_tendency), col = "red")     
    invisible(coords)
  }
  else coords
}

conerrbox <- function(centrum, border, at = 1, width, border_par, line_par) {
  x_coord <- c(at - 0.5*width, at + 0.5*width)
  do.call(rect, c(list(xleft = x_coord[1], ybottom = centrum-border, xright = x_coord[2], 
                       ytop = centrum+border), border_par))
  do.call(lines, c(list(x = x_coord, y = c(centrum, centrum)), line_par))
}

conerrboxes <- function(conplot_data, width = 0.4, errbox_border = list(), 
                        errbox_line = list()) {
  if (length(errbox_border) != 0) 
    border_par = errbox_border else border_par <- list()
  if (!("border" %in% names(border_par))) {
    bord_cols <- rep_len(c("skyblue1", "grey54"), max(conplot_data@group_order))
    if (max(conplot_data@group_order) %% 2 == 1) 
      bord_cols <- c(bord_cols, "skyblue1")
    border_pars <- lapply(bord_cols, function(x) c(list(border = x), border_par))
  }
  if (length(errbox_line) != 0) 
    box_line_par = errbox_line else box_line_par <- list()
  invisible(sapply(conplot_data@group_order, 
                   function(x) conerrbox(conplot_data@aggr[x,][conplot_data@mid], 
                                         conplot_data@aggr[x,][conplot_data@mid + 1], at = x, width, 
                                         border_par = border_pars[[x]], line_par = box_line_par)))
}

constrip <- function(conplot_data, plot = TRUE, errbox = TRUE, line = TRUE, 
                     xlab = NULL, ylab = NULL, ...) {
  
  NOC <- max(conplot_data@group_order)
  d <- split(conplot_data@central_tendency, 
             conplot_data@groups)[conplot_data@group_order]
  
  if (plot) {
    if (is.null(xlab)) 
      xlab <- ""
    if (is.null(ylab)) 
      ylab <- ""    
    
    plot(NA, NA, xlab = xlab, ylab = ylab, xlim = c(1L,NOC), 
         ylim = c(min(conplot_data@central_tendency), 
                  max(conplot_data@central_tendency)), xaxt = "n", ...)
    sapply(1L:NOC, function(at) sapply(d[[at]], function(y) points(at, y)))
    conaxis(conplot_data)
    
    if (errbox == TRUE) { 
      conerrboxes(conplot_data, width = 0.033*NOC + 0.01)      
    }
    
    if (line) #central tendency of all central tendencies
      if (conplot_data@mid == 1) {
        abline(h = mean(conplot_data@central_tendency), col = "red")
      } else abline(h = median(conplot_data@central_tendency), col = "red")  
    invisible(d)
  } else d
}



conbar <- function(conplot_data, plot = TRUE, legend = TRUE, text = TRUE, width = 0.4, 
                   barplot_par = list(), text_par = list(), xlab = "", ylab = "", ...) {
  NOC <- max(conplot_data@group_order)

    if (plot) {
      plot(NA, NA, ylab = ylab, xlab = xlab, xlim = c(1L,NOC), ylim = c(0, 1.3), 
           axes = FALSE, bty="n", ...)
      
      
      labels <- rownames(conplot_data@freqs)
      n <- nrow(conplot_data@freqs)
      w <- width
      
      cols = rainbow(n)
      if (!is.null(names(barplot_par))) {
        if ("col" %in% names(barplot_par)) {
          cols <- barplot_par[["col"]] 
          barplot_par <- barplot_par[-which(names(barplot_par) == "col")]
        } 
      }
      if (length(cols) != n) {
        msg <- paste0("'col' must be a vector of colors with length", n, ".")
        stop(msg, call. = TRUE, domain = NA)
      }
      data <- rbind(0, apply(rbind(conplot_data@freqs), 2, cumsum))
      sapply(1L:ncol(data), function(x) sapply(1L:n, function(i) do.call(rect, 
                                                                         c(list(x - w, 
                                                                                data[i,x], 
                                                                                x + w, 
                                                                                data[i + 1,x]), 
                                                                           col = cols[i], 
                                                                           barplot_par))))
      if (legend) {
        at_ptsx <- seq(from = 1, to = NOC, length.out = 6)
        at_ptsy <- 1 + ceiling(1:n/6)*0.1
        at_txt <- seq(from = 1 + NOC*0.007, to = NOC, length.out = 6)
        sapply(1L:n, function(at) points(at_ptsx[at], at_ptsy[n], col = cols[at], pch = 15, 
                                         cex = 1.3))
        sapply(1L:n, function(at) text(at_txt[at], at_ptsy[n], labels = labels[at], pos = 4,
                                       cex = 0.9))
      }
      
      if (text) conbartext(conplot_data, cols, text_par = text_par)
      
      axis(2, at = 0:5/5, labels = 0:5*20)
      mtext("Frequency [%]", side = 2, at = 0.5, line = 3)
      invisible(conplot_data@freqs)
    } else conplot_data@freqs
}

get_contr_col <-function (color) { 
  rgb <- col2rgb(color)/255 
  L <- c(0.2, 0.6, 0) %*% rgb 
  ifelse(L >= 0.2, "#000060", "#FFFFA0") 
}

round_rounded <-function (vec, n) {
  r_vec = round(vec, n)
  if (sum(vec) - sum(r_vec) != 0) {
    r_vec[which.max(vec - r_vec)] <- r_vec[which.max(vec - r_vec)] + sum(vec) - sum(r_vec)
  }
  r_vec
}


conbartext <- function(conplot_data, contr_color = NULL, at = NULL, text_par = list()) {
  
  n <- nrow(conplot_data@freqs)
  data <- rbind(0, apply(conplot_data@freqs, 2, cumsum))
  pos <- sapply(1L:ncol(data), function(x) sapply(1L:n, function(i) mean(data[i:(i+1),x])))
  labels <- apply(conplot_data@freqs, 2, function(x) round_rounded(x, 2)*100)
  
  if (is.null(at)) 
    at <- 1L:ncol(pos)
  
  if (is.null(contr_color)) 
    contr_color <- gray.colors(n)
  cols <- sapply(contr_color, get_contr_col)
  
  sapply(1L:ncol(pos), function(x) sapply(1L:n, function(y) if (labels[y,x] != 0)
    do.call(text, c(list(at[x], pos[y,x], labels[y,x], col = cols[y]), text_par))))
  list(y = posy, labels = labels, colors = cols)
}


#concrete cluster lines
conclines <- function(conplot_data, ylab = "", conclines_par = list()) {
  conclines_par <- set_default("col", "red", conclines_par)
  
  scaled_centr <- conplot_data@aggr[[conplot_data@mid]]
  NOC <- max(conplot_data@group_order)
  
  par(new = TRUE)
  plot(c(1, NOC), c(0, max(conplot_data@aggr[[conplot_data@mid]])), axes = FALSE, xlab = "", 
       ylab = "", type="n")
  sapply(1L:NOC, function(x) do.call(lines, c(list(x = c(x, x), y = c(0, scaled_centr[x])), 
                                              conclines_par)))
  axis(side = 4, at = axTicks(4), labels=axTicks(4), cex.axis = 0.95, 
       col.axis=conclines_par[["col"]], 
       col = conclines_par[["col"]])
  mtext(ylab, 4, col = conclines_par[["col"]], line = 2)
  axis(side = 1, at = 1L:NOC, labels = rep("", NOC), pos = c(1, NOC), cex.axis = 0.95, 
       lwd.ticks = 0, 
       col = conclines_par[["col"]])
  #par(new = TRUE)
}

set_default <- function(arg_name, default_val, arg_list) {
  if (!arg_name %in% names(arg_list))
    arg_list[[arg_name]] <- default_val
  arg_list
}



concaption  <- function(conplot_data, orientation, col = NULL, text_par = list(),
                        ticks_par = list()) {
  text_par <- set_default("las", 1, text_par)
  text_par <- set_default("cex", 0.75, text_par)
  
  ticks_par <- set_default("lwd.ticks", 1, ticks_par)
  ticks_par <- set_default("tcl", -0.8, ticks_par)
  
  if (!is.null(col)) {
    ticks_par <- set_default("col", col, ticks_par)
    text_par <- set_default("col", col, text_par)
  }
  
  do.call(mtext, c(list(conplot_data@group_order, orientation, 
                        at = 1L:length(conplot_data@group_order), line = 1), text_par))
  do.call(axis, c(list(orientation, at = 1L:length(conplot_data@group_order), 
                       labels = FALSE, las=2, lwd = 0), ticks_par))
}


#concrete quick cluster
conqclus <- function(conclust, conplot_data, orientation = 3, conclines = TRUE,
                     conclines_par = list(), perp_line_par = list(), 
                     paral_lines_par = list()) {
  plot_cluster(conclust, orientation, perp_line_par, paral_lines_par, caption = FALSE)
  if (conclines) {
    if (!"ylab" %in% names(conclines_par)) conclines_par[["ylab"]] <- ""
    conclines(conplot_data, conclines_par[["ylab"]], conclines_par = conclines_par)
  }
  
  #improve it
  invisible(0)
}


concombo <- function(conplot_data, dendrogram, ploty = "1", up_par = list(), 
                     bot_par = list()) {
  
  ploty <- switch(ploty,
                  "1" = c("conbox", "conqclus"),
                  "2" = c("constrip", "conqclus"),
                  "3" = c("conbar", "conqclus"))
  
  old_par <- par("fig") #save old parameter
  par(fig = c(0, 0.95, 0.4, 1))
  do.call(ploty[1], c(list(conplot_data), up_par))
  par(fig = c(0, 0.95, 0, 0.58), new = TRUE)
  do.call(ploty[2], c(list(dendrogram, conplot_data), bot_par))
  par(fig = old_par)
}

draw_lines <- function(coords, vertical = TRUE, perp_line_par, paral_lines_par) {
  x1 <- c(coords[1], coords[1])
  y1 <- c(coords[3], coords[4])
  x2 <- c(coords[2], coords[2])
  y2 <- c(coords[5], coords[6])
  x3 <- coords[1:2]
  y3 <- c(coords[4], coords[4])
  if (vertical == TRUE) {
    do.call(lines, c(list(x1, y1), paral_lines_par))
    do.call(lines, c(list(x2, y2), paral_lines_par)) 
    do.call(lines, c(list(x3, y3), perp_line_par))
  } else {
    do.call(lines, c(list(y1, x1), paral_lines_par))
    do.call(lines, c(list(y2, x2), paral_lines_par))
    do.call(lines, c(list(y3, x3), perp_line_par))  
  }
}



coords <- function(mergeRow, merge, order, n_height, height, coords_table) {
  xc <- c()
  yc <- c()
  test <- sum(mergeRow > 0)
  if (test == 0) {
    xc <- c(xc, sapply(abs(mergeRow), function(x) match(x, order)))
    yc <- c(yc, 0, n_height, 0, n_height)
  }
  if (test == 1) {
    xc <- c(xc, match(abs(mergeRow[1]), order), mean(coords_table[mergeRow[2],1:2]))
    yc <- c(yc, 0, n_height, height[mergeRow[2]], n_height)
  }
  if (test == 2) {
    xc <- c(xc, mean(coords_table[mergeRow[1],1:2]), mean(coords_table[mergeRow[2],1:2]))
    yc <- c(yc, height[mergeRow[1]], n_height, height[mergeRow[2]], n_height)
  }
  c(xc, yc)
}

#nonpretty version


#pretty
plot_conclust <- function(conclust, orientation = 1, perp_line_par = list(), 
                         paral_lines_par = list(), caption = TRUE) {
  merge <- conclust@merge
  order <- conclust@order
  height <- conclust@height
  labels <- conclust@labels
  
  if (is.vector(merge)) merge <- matrix(merge, ncol = 2)
  
  len <- length(order)
  n <- length(height)
  
  coords_table <- matrix(-1L, ncol = 6, nrow = n)
  
  for (i in 1L:n) {
    coords_table[i,] <- coords(merge[i,], merge, order, height[i], height, coords_table)
  }
  
  if (orientation == 3 | orientation == 4) coords_table[,3:6] <- coords_table[,3:6]*-1
  
  if (orientation == 3 | orientation == 1) {
    plot(c(1, len), c(0, coords_table[length(height),4]), cex = 0, ylab='Height',xlab='', axes=F)
    apply(coords_table, 1, function(x) draw_lines(x, vertical = TRUE, perp_line_par, paral_lines_par))
    axis(2, at = axTicks(2), labels = abs(axTicks(2)))
  }
  
  if (orientation == 2 | orientation == 4) {
    plot(c(0, coords_table[length(height),4]), c(1, len), cex = 0, xlab='Height',ylab='', axes=F)
    apply(coords_table, 1, function(x) draw_lines(x, vertical = FALSE, perp_line_par, paral_lines_par))
    axis(1, at = axTicks(1), labels = abs(axTicks(1)))
  }
  
  if (caption == TRUE) {
    mtext(order, orientation, at = 1L:len, las=1, cex = 0.75, line = 1)
    axis(orientation, at = 1L:len, labels = FALSE, las=2, lwd.ticks = 1, lwd = 0, tcl = -0.8)
  }
  
  coords_table
}


find_row <- function(vector_match, vector_matrix) 
  ceiling(sapply(vector_match, function(x) match(x, vector_matrix))/2)

#minimum node height
treetab <- function(clust_data) {
  merge <- clust_data@merge
  height <- clust_data@height
  order <- clust_data@order
  labels <- clust_data@labels
  
  merge_vector <-  as.vector(t(merge))
  
  n_l <- length(labels) # length of leaves  
  n_n <- length(height) # length of nodes  
  
  leaf_id <- 1L:n_l
  node_id <- 1L:n_n + n_l
  node_id_short <- 1L:n_n
  
  leaf_node_names <- c(labels, sapply(1L:n_n, function(x) paste0("node", x)))
  parent_node_id <- c(find_row(leaf_id * -1, merge_vector), find_row(node_id_short, merge_vector)) + n_l
  dendr_table <- data.frame(id = c(leaf_id, node_id), parent_id = parent_node_id, height = c(rep(0, n_l), height),
                            order = c(order, rep(NA, n_n)))
  rownames(dendr_table) <- leaf_node_names
  dendr_table
}

tabtree <- function(x) {
  node_indices <- grepl("node", rownames(x))
  height <- x$height[node_indices]
  labels <- rownames(x)[!node_indices]
  order <- x$order[!node_indices]
  
  n_lab <- length(labels)
  n_node <- length(height)
  
  merge <- t(sapply(1L:n_node, function(i) as.integer(x$id[x$parent_id == x$id[node_indices][i]]))[-3,])
  merge[merge <= n_lab] <- merge[merge <= n_lab]*as.integer(-1)
  merge[merge > 0] <- merge[merge > 0] - n_lab
  storage.mode(new_merge) <- "integer"
  
  data <- list(merge = merge, height = height, order = order, labels = labels, method = "tabletree", 
               call = match.call(), dist.method = "unknown") 
  class(data) <- "hclust"
  data
}

order_cluster <- function(merge) {
  if (is.vector(merge)) abs(merge)
  else {
    n <- nrow(merge)
    negative1 <- sum(merge[,1] < 0)
    negative2 <- sum(merge[,2] < 0)
    
    l_order <- lapply(1L:negative2, function(x) merge[x,])
    
    for (i in 1L:n) {
      test <- sum(merge[i,] > 0)
      if (test == 0) {
        l_order[[i]] <- merge[i, ]
      } 
      if (test == 1) {
        l_order[[i]] <- c(l_order[[merge[i,2]]], merge[i, 1])
      }
      if (test == 2) {
        l_order[[i]] <- c(l_order[[merge[i,1]]], l_order[[merge[i,2]]])
      }
    }
    
    abs(l_order[[n]])
  }
}

trim_dendrogram <- function(conclust, min_height, s_names = FALSE) {  
  
  merge <- conclust@merge
  merge_vector <-  as.vector(t(merge))
  height <- conclust@height
  order <- conclust@order
  
  if (max(height) <= min_height)
    stop("Minimum height cannot equal or exceed maximum height of trimmed dendrogram.", call. = TRUE, domain = NA)
  
  gr <- cutree(conclust, h = min_height)
  vec_gr <- as.vector(gr)
  #change to list of vectors, but keep possibility of naming
  new_names <- as.vector(sapply(split(names(gr), vec_gr), function(x) paste0(x, collapse="+")))
  
  nodes_id <- which(height > min_height)
  new_merge <- merge[nodes_id, ]
  new_leaves_in <- new_merge < min(nodes_id)
  new_merge[new_leaves_in] <- (1:length(new_leaves_in))[order(new_merge[new_leaves_in])]*-1
  new_nodes_in <- new_merge > 0
  new_merge[new_nodes_in] <- (new_merge[new_nodes_in] - min(nodes_id) + 1)
  
  storage.mode(new_merge) <- "integer"
  #   data <- list(merge = new_merge, height = height[nodes_id] - min_height, order = order_cluster(new_merge),
  #                labels = new_names, method = "trim", call = match.call(), dist.method = conclust@dist.method) 
  #   class(data) <- "hclust" 
  #   if (s_names == TRUE) {
  #     data$labels <- as.vector(sapply(data$labels, function(x) substr(x, 0, 8)))
  #     data$long_labels <- new_names
  #   }
  slot(conclust, "merge") <- new_merge
  slot(conclust, "height") <- height[nodes_id] - min_height
  slot(conclust, "order")  <- order_cluster(new_merge)
  slot(conclust, "labels") <- new_names
  slot(conclust, "method") <- "Trimmed dendrogram"
  slot(conclust, "call") <- match.call()
  slot(conclust, "prev.height") <- c(slot(conclust, "prev.height"), max(height))
  
  if (s_names == TRUE) {
    slot(conclust, "labels") <- as.vector(sapply(new_names, function(x) substr(x, 0, 8)))
    slot(conclust, "long.labels") <- new_names
  }
  conclust
}

setOldClass("hclust")

setClass("conclust", contains = "hclust", representation(merge = "matrix", order = "integer", 
                                                         height = "numeric", prev.height = "numeric", 
                                                         labels = "character", long.labels = "character",
                                                         method = "character", dist.method = "character",
                                                         call = "call"))


as.conclust <- function(x){
  if (class(x) != "hclust")
    stop("Currently only 'hclust' class is supported.", call. = TRUE, domain = NA)
  
  if (!exists("x"))
    stop("Argument 'x' is missing.", call. = TRUE, domain = NA)
  
  if (is.null(x$labels))
    x$labels <- as.character(1L:max(x$order))
  
  if (is.na(x$method))
    x$method <- "method unknown"
  conclust <- new("conclust")
  old_slots <- names(x)
  for (slot_name in old_slots) slot(conclust, slot_name) <- x[[slot_name]]
  conclust
}

setGeneric("as.hclust")
setMethod("as.hclust", signature(x = "conclust"), function(x, expanded = FALSE) {
  hclust <- list()
  labels <- c("merge","height","order","labels","method","call","dist.method")
  if (expanded == T) labels <- slotNames(x)
  for (slot_name in labels) hclust[[slot_name]] <- slot(x, slot_name)
  class(hclust) <- "hclust"
  hclust
})

setGeneric("print")
setMethod("print", signature(x = "conclust"), function(x) {
  if (!identical(x@call, character(0)))
    cat("\nCall:\n", deparse(x@call), "\n\n", sep = "")
  if (!identical(x@method, character(0)))
    cat("Cluster method    :", x@method, "\n")
  if (!identical(x@dist.method, character(0)))
    cat("Distance          :", x@dist.method, "\n")
  if (!identical(x@prev.height, numeric(0))) 
    cat("Previous height(s):", x@prev.height, "\n")
  cat("Number of objects :", length(x@height) + 1, "\n")
  cat("\n")
  invisible(x)
})

setGeneric("show")
setMethod("show", signature(object = "conclust"), function(object) {
  if (!identical(object@call, character(0)))
    cat("\nCall:\n", deparse(object@call), "\n\n", sep = "")
  if (!identical(object@method, character(0)))
    cat("Cluster method    :", object@method, "\n")
  if (!identical(object@dist.method, character(0)))
    cat("Distance          :", object@dist.method, "\n")
  if (!identical(object@prev.height, numeric(0))) 
    cat("Previous height(s):", object@prev.height, "\n")
  cat("Number of objects :", length(object@height) + 1, "\n")
  cat("\n")
  invisible(object)
})

setGeneric("treetab")
setMethod("treetab", signature("hclust"), function(clust_data) {  
  treetab(as.conclust(clust_data))
})

setGeneric("cutree")
setMethod("cutree", signature(tree = "conclust"), function(tree, k, h) {  
  tree <- as.hclust(tree)
  cutree(tree, k, h)
})

setGeneric("concreate", function(central_tendency, group_order, groups, dendrogram, 
                                 height, formula, data, binary_data, binary_labels, 
                                 ...) {
  stop("Error text", call. = TRUE, domain = NA)
})

setMethod("concreate", signature(central_tendency = "ANY", group_order = "missing", 
                                 groups = "missing",
                                 dendrogram = "ANY", height = "numeric",
                                 formula = "missing", data = "missing"), 
          function(central_tendency, dendrogram, height, binary_data, 
                   binary_labels, ...) {
            if (!class(dendrogram) %in% c("hclust", "conclust"))
              stop("Invalid type of dendrogram", call. = TRUE, domain = NA)
            trimmed <- contrim(dendrogram, min_height = height, s_names = TRUE)
            calc_plots(central_tendency, group_order = trimmed@order, 
                       groups = cutree(dendrogram, h = height), 
                       binary_data = binary_data, binary_labels = binary_labels, ...)
          })

setMethod("concreate", signature(central_tendency = "ANY", group_order = "integer", 
                                 groups = "integer",
                                 dendrogram = "missing", height = "missing", 
                                 formula = "missing", data = "missing"), 
          function(central_tendency, group_order, groups, binary_data, 
                   binary_labels, ...) {
            calc_plots(central_tendency, group_order = group_order, 
                       groups = groups, binary_data = binary_data, 
                       binary_labels = binary_labels, ...)
          })

setMethod("concreate", signature(central_tendency = "missing", group_order = "integer", 
                                 groups = "integer", dendrogram = "missing", 
                                 height = "missing", formula = "formula"), 
          function(group_order, groups, formula, data, binary_labels, ...) {
            if (length(formula) != 3) 
              stop("'formula' incorrect")
            if (exists("data")) {
              if (is.matrix(data)) 
                data <- as.data.frame(data)
              mf <- model.frame(formula, data)
            } else {
              mf <- model.frame(formula)
            }
            response <- attr(attr(mf, "terms"), "response")
  calc_plots(mf[[response]], group_order = group_order, 
             groups = groups, binary_data = mf[-response], 
             binary_labels = binary_labels, ...)
          })

setMethod("concreate", signature(central_tendency = "missing", group_order = "missing", 
                                 groups = "missing", dendrogram = "ANY", 
                                 height = "numeric", formula = "formula"), 
          function(dendrogram, height, formula, data, binary_labels, ...) {
            if (!class(dendrogram) %in% c("hclust", "conclust"))
              stop("Invalid type of dendrogram", call. = TRUE, domain = NA)
            if (length(formula) != 3) 
              stop("'formula' incorrect")
            if (!missing(data)) {
              if (is.matrix(data)) 
                data <- as.data.frame(data)
              mf <- model.frame(formula, data)
            } else {
              mf <- model.frame(formula)
            }
            response <- attr(attr(mf, "terms"), "response")
            trimmed <- contrim(dendrogram, min_height = height, s_names = TRUE)
            calc_plots(mf[[response]], group_order = trimmed@order, 
                       groups = cutree(dendrogram, h = height), 
                       binary_data = mf[-response], 
                       binary_labels = binary_labels, ...)
          })


setGeneric("contrim", function(conclust, min_height, s_names = FALSE) {
  stop("This function requires an object of classes 'conclust' or 'hclust'.", 
       call. = TRUE, domain = NA)
})

setMethod("contrim", signature(conclust = "hclust"), function(conclust, min_height, s_names = FALSE) {  
  trim_dendrogram(as.conclust(conclust), min_height, s_names)
})
setMethod("contrim", signature(conclust = "conclust"), function(conclust, min_height, 
                                                             s_names = FALSE) {  
  trim_dendrogram(conclust, min_height, s_names)
})

setGeneric("plot_cluster", function(conclust, orientation = 1, perp_line_par = list(), 
                                    paral_lines_par = list(), caption = TRUE) {
  stop("This function requires an object of class 'conclust' or 'hclust'.", 
       call. = TRUE, domain = NA)
})
setMethod("plot_cluster", signature(conclust = "conclust"), 
          function(conclust, orientation, perp_line_par, paral_lines_par) {  
  plot_conclust(conclust, orientation, perp_line_par, paral_lines_par)
})

setMethod("plot_cluster", signature(conclust = "hclust"), 
          function(conclust, orientation, perp_line_par, paral_lines_par) {  
  plot_conclust(as.conclust(conclust), orientation, perp_line_par, paral_lines_par)
})



eagg <- function(id_vars, measure_var) {
  aggregation <- aggregate(measure_var, by = id_vars, function(x) cbind(mean(x), median(x), sd(x), mad(x), length(x)))
  n <- ncol(aggregation)
  result <- data.frame(aggregation[-n], aggregation[[n]][,-5], aggregation[[n]][,3]/aggregation[[n]][,1]*100, aggregation[[n]][,4]/aggregation[[n]][,2]*100, 
                       aggregation[[n]][,5], aggregation[[n]][,1] - aggregation[[n]][,2]) 
  names(result)[(ncol(result)-7):ncol(result)] <- c("mean", "median", "sd", "mad", "CV", "MADCV", "length", "delta")
  result
}





