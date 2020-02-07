.addIsDECol <- function(
  data,
  testCol = "padj",
  alpha,
  lfcCol = "log2FoldChange",
  lfcThreshold = 0L
) {
  # test: P value or S value
  test <- data[[testCol]]
  # lfc: log2 fold change cutoff
  lfc <- data[[lfcCol]]
  isDE <- mapply(
    test = test,
    lfc = lfc,
    FUN = function(test, lfc) {
      if (any(is.na(c(test, lfc)))) {
        # nonsignificant
        0L
      } else if (test < alpha & lfc > lfcThreshold) {
        # upregulated
        1L
      } else if (test < alpha & lfc < -lfcThreshold) {
        # downregulated
        -1L
      } else {
        0L
      }
    },
    SIMPLIFY = TRUE,
    USE.NAMES = FALSE
  )
  isDE <- as.factor(isDE)
  data[["isDE"]] <- isDE
  data
}

plot_niceMA <- function(
  object,
  alpha,
  lfcThreshold = 0L,
  genes = NULL,
  gene2symbol = NULL,
  ntop = 0L,
  direction = c("both", "up", "down"),
  pointColor = "gray50",
  sigPointColor = c(
    upregulated = "purple",
    downregulated = "orange"
  ),
  return = c("ggplot", "data.frame")
) {
  direction <- match.arg(direction)
  
  
  return <- match.arg(return)
  
  # Check to see if we should use `sval` instead of `padj`
  
  testCol <- "padj"
  
  lfcCol <- "log2FoldChange"
  
  data <- res %>%
    as.data.frame() %>%
    rownames_to_column("geneID") %>%
    as_tibble() %>%
    janitor::clean_names(case = "lower_camel") %>%
    # Remove genes with very low expression
    dplyr::filter(!!sym("baseMean") >= 1L) %>%
    mutate(rankScore = abs(!!sym("log2FoldChange"))) %>%
    arrange(desc(!!sym("rankScore"))) %>%
    mutate(rank = row_number()) %>%
    .addIsDECol(
      testCol = testCol,
      alpha = alpha,
      lfcCol = lfcCol,
      lfcThreshold = lfcThreshold
    )
  
  if (direction == "up") {
    data <- data[data[[lfcCol]] > 0L, , drop = FALSE]
  } else if (direction == "down") {
    data <- data[data[[lfcCol]] < 0L, , drop = FALSE]
  }
  
  # Gene-to-symbol mappings
  if (is.data.frame(gene2symbol)) {
    labelCol <- "geneName"
    data <- left_join(data, gene2symbol, by = "geneID")
  } else {
    labelCol <- "geneID"
  }
  
  # Early return data frame, if desired
  if (return == "data.frame") {
    data <- data %>%
      as.data.frame() %>%
      column_to_rownames("geneID")
    return(data)
  }
  
  xFloor <- data[["baseMean"]] %>%
    min() %>%
    log10() %>%
    floor()
  xCeiling <- data[["baseMean"]] %>%
    max() %>%
    log10() %>%
    ceiling()
  xBreaks <- 10L ^ seq(from = xFloor, to = xCeiling, by = 1L)
  
  p <- ggplot(
    data = data,
    mapping = aes(
      x = !!sym("baseMean"),
      y = !!sym(lfcCol),
      color = !!sym("isDE")
    )
  ) +
    geom_hline(
      yintercept = 0L,
      size = 0.5,
      color = pointColor
    ) +
    geom_point(size = 1L) +
    scale_x_continuous(
      breaks = xBreaks,
      limits = c(1L, NA),
      trans = "log10"
    ) +
    scale_y_continuous(breaks = pretty_breaks()) +
    annotation_logticks(sides = "b") +
    guides(color = FALSE) +
    labs(
      x = "mean expression across all samples",
      y = "log2 fold change"
    )
  
  p <- p +
    scale_color_manual(
      values = c(
        # nonsignificant
        "0" = pointColor,
        # upregulated
        "1" = sigPointColor[[1L]],
        # downregulated
        "-1" = sigPointColor[[2L]]
      )
    )
  
  
  
  p
}
