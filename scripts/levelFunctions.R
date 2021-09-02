factor_relevel <- function(x, reference) {
  relev <- relevel(as.factor(gsub("^\\*", "", as.character(x))),
                   ref = reference[cur_column()])
  return(relev)
}

color_levels <- function(df) {
  all_lvls <- lvls_colors <- lapply(df, levels)
  cond_colors <- hue_pal()(length(unlist(all_lvls[["condition"]])))
  rest_colors <- hue_pal()(length(unlist(all_lvls)))
  all_colors <- c(cond_colors, 
                  rest_colors[-c(which(rest_colors %in% cond_colors))])
  i <- 1
  for (n in names(lvls_colors)) {
    x <- lvls_colors[[n]]
    lvls_colors[[n]] <- setNames(all_colors[i:(i - 1 + length(x))], x)
    i <- length(x) + 1
  }
  return(lvls_colors)
}