meta_plot_generate_density_ci <- function(relative_points, point_weight, adjust, enableCI, CI_ResamplingTime = 1000, CI_interval = c(0.025, 0.975)) {
  density_dataframe <- data.frame()
  for (group_name in names(relative_points)) {
    valid <- is.finite(point_weight[[group_name]]) & point_weight[[group_name]] > 0 & is.finite(relative_points[[group_name]])
    point_weight[[group_name]] <- point_weight[[group_name]][valid]
    relative_points[[group_name]] <- relative_points[[group_name]][valid]
    if (length(relative_points[[group_name]]) < 2L) next
    site_weight <- point_weight[[group_name]] / sum(point_weight[[group_name]])
    site_id <- relative_points[[group_name]]
    fit <- suppressWarnings(stats::density(site_id, adjust = adjust, from = 0, to = 1, n = 256, weight = site_weight))
    tmp <- data.frame(x = fit$x, density = fit$y / (sum(diff(fit$x) * (head(fit$y, -1) + tail(fit$y, -1))) / 2), group = rep(group_name, length(fit$y)))
    if (enableCI) {
      point_index <- seq_along(site_id)
      names(point_index) <- names(site_id)
      point_grouped <- split(point_index, names(point_index))
      fit_resampled <- replicate(CI_ResamplingTime, {
        resampled_groups <- sample(names(point_grouped), replace = TRUE)
        resampled_index <- unlist(point_grouped[resampled_groups])
        suppressWarnings(stats::density(site_id[resampled_index], adjust = adjust, from = 0, to = 1, n = 256, weight = site_weight[resampled_index])$y)
      })
      fit_quantile <- apply(fit_resampled, 1, stats::quantile, CI_interval)
      tmp$confidenceDown <- fit_quantile[1, ]
      tmp$confidenceUp <- fit_quantile[2, ]
    }
    density_dataframe <- rbind(density_dataframe, tmp)
  }
  density_dataframe
}
