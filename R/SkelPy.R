#' @import dplyr
#' @import ggplot2
#' @import jsonlite
#' @import magrittr
#' @import purrr
#' @import sf
#' @import sfnetworks
#' @import tibble
#' @import tidygraph
#' @import patchwork
#' @importFrom magrittr %>%
NULL




# Minimal fallback for to_spatial_smooth()
# ------------------------------------------------------------------
if (!exists("to_spatial_smooth")) {
  to_spatial_smooth <- function(graph, summarise_attributes = "first") {
    graph
  }
}

#   Folder Selection Helper functions, allows users to select Calculation folder

#' Select folder containing JSON files
#'
#' Uses RStudio dialog if available, otherwise terminal input.
#' @return A character string: path to folder
#' @export
select_json_folder <- function() {
  if (requireNamespace("rstudioapi", quietly = TRUE) &&
      rstudioapi::isAvailable()) {
    message("Please select the folder containing your JSON files.")
    dir <- rstudioapi::selectDirectory()
    return(dir)
  } else {
    message("Enter full path to the folder containing your JSON files:")
    dir <- readline("> ")
    return(dir)
  }
}

#' Load JSON files from a folder
#'
#' @param base_dir Folder that contains .json files
#' @return Named character vector of JSON file paths
#' @export
load_json_files <- function(base_dir) {

  if (!dir.exists(base_dir)) {
    stop("The folder does not exist: ", base_dir)
  }

  json_files <- list.files(base_dir, pattern = "\\.json$", full.names = TRUE)

  if (length(json_files) == 0) {
    stop("No JSON files found in: ", base_dir,
         "\nMake sure the folder contains files ending in .json")
  }

  names(json_files) <- basename(json_files)

  message("Loaded ", length(json_files), " JSON files from:\n", base_dir)
  return(json_files)
}

#   Core Processing

#' Safe wrapper around eRgot(), prevents crashes from early timepoints with no networks
#' @export
safe_eRgot <- function(path) {
  tryCatch(
    {
      eRgot(path)
    },
    error = function(e) {
      message("Skipping file due to error in ", basename(path), ": ", e$message)
      return(NULL)
    }
  )
}

#' Run eRgot() on all JSON files
#' @export
run_ergot_on_all <- function(json_files) {

  network_results <- lapply(json_files, safe_eRgot)

  keep <- !sapply(network_results, is.null)
  network_results <- network_results[keep]
  json_files <- json_files[keep]
  names(network_results) <- basename(json_files)

  message("Successfully processed ", length(network_results), " valid JSON files.")

  if (length(network_results) == 0) {
    stop("All calls to eRgot() failed. Check the earlier error messages.")
  }

  list(
    network_results = network_results,
    json_files      = json_files
  )
}

#   Summary / Metrics Helpers

#' Cleaning summary for hyphae networks
#' @export
compute_cleaning_summary <- function(network_results) {
  purrr::map_dfr(
    seq_along(network_results),
    function(i) {
      item <- network_results[[i]]$network_hyphae$cleaning_summary

      tibble::tibble(
        file           = names(network_results)[i],
        original_edges = item$original[1],
        cleaned_edges  = item$processed[1],
        original_nodes = item$original[2],
        cleaned_nodes  = item$processed[2]
      )
    }
  )
}

#' Node-level metrics (hyphae)
#' @export
compute_node_stats <- function(network_results) {
  purrr::map_dfr(
    seq_along(network_results),
    function(i) {
      net <- network_results[[i]]$network_hyphae$clean_net

      net %>%
        activate(nodes) %>%
        mutate(
          file       = names(network_results)[i],
          degree     = centrality_degree(),
          betweenness= centrality_betweenness(),
          closeness  = centrality_closeness(),
          component  = group_components()
        ) %>%
        as_tibble()
    }
  )
}

#' Edge-level metrics (hyphae)
#' @export
compute_edge_stats <- function(network_results) {
  purrr::map_dfr(
    seq_along(network_results),
    function(i) {
      net <- network_results[[i]]$network_hyphae$clean_net

      edge_df <- net %>%
        activate(edges) %>%
        mutate(
          file   = names(network_results)[i],
          length = as.numeric(sf::st_length(geometry))
        ) %>%
        as_tibble()

      edge_df$geometry <- NULL
      edge_df$.tidygraph_edge_index <- NULL

      edge_df
    }
  )
}

#' Hyphae-level metrics
#' @export
compute_hyphae_summary <- function(network_results) {
  purrr::map_dfr(
    seq_along(network_results),
    function(i) {

      fname <- names(network_results)[i]
      net   <- network_results[[i]]$network_hyphae$clean_net

      # Edge lengths
      edge_df <- net %>%
        activate(edges) %>%
        mutate(edge_length = as.numeric(sf::st_length(geometry))) %>%
        as_tibble()

      edge_df$geometry <- NULL
      edge_df$.tidygraph_edge_index <- NULL

      n_edges       <- nrow(edge_df)
      total_length  <- sum(edge_df$edge_length, na.rm = TRUE)
      mean_length   <- mean(edge_df$edge_length, na.rm = TRUE)
      median_length <- median(edge_df$edge_length, na.rm = TRUE)
      sd_length     <- sd(edge_df$edge_length, na.rm = TRUE)
      min_length    <- min(edge_df$edge_length, na.rm = TRUE)
      max_length    <- max(edge_df$edge_length, na.rm = TRUE)

      # Nodes
      node_df <- net %>%
        activate(nodes) %>%
        mutate(degree = centrality_degree()) %>%
        as_tibble()

      node_df$geometry <- NULL

      n_nodes    <- nrow(node_df)
      n_tips     <- sum(node_df$degree == 1, na.rm = TRUE)
      n_branches <- sum(node_df$degree >= 3, na.rm = TRUE)

      tibble::tibble(
        file              = fname,
        n_edges           = n_edges,
        total_length      = total_length,
        mean_edge_length  = mean_length,
        median_edge_length= median_length,
        sd_edge_length    = sd_length,
        min_edge_length   = min_length,
        max_edge_length   = max_length,
        n_nodes           = n_nodes,
        n_tips            = n_tips,
        n_branches        = n_branches,
        tip_fraction      = n_tips / n_nodes,
        branch_fraction   = n_branches / n_nodes,
        tip_density       = n_tips / total_length
      )
    }
  )
}

#   Plot Saving

#' Save cleaned hyphae plots
#' @export
save_hyphae_plots <- function(network_results, plots_dir = "plots_json") {

  if (length(network_results) == 0) {
    message("No network results available; skipping plot saving.")
    return(invisible(NULL))
  }

  dir.create(plots_dir, showWarnings = FALSE)

  for (i in seq_along(network_results)) {
    p <- network_results[[i]]$clean_maps$hyphae_plot

    ggsave(
      filename = paste0(names(network_results)[i], "_hyphae_clean.png"),
      plot     = p,
      path     = plots_dir,
      width    = 10,
      height   = 8,
      dpi      = 300
    )
  }

  message("Saved cleaned hyphae plots to: ", plots_dir)
}


#   Main Pipeline Wrapper


#' Run full FySkel pipeline on a folder of JSON files
#'
#' @param base_dir Path to folder with JSON files (optional)
#' @param output_dir Output directory for CSV + PNG files
#' @return A list with results and paths
#' @export
run_fyskel_pipeline <- function(base_dir = NULL, output_dir = getwd()) {

  if (is.null(base_dir)) {
    base_dir <- select_json_folder()
  }

  json_files <- load_json_files(base_dir)
  res        <- run_ergot_on_all(json_files)
  network_results <- res$network_results

  old_wd <- getwd()
  on.exit(setwd(old_wd), add = TRUE)
  setwd(output_dir)

  cleaning_summary_all <- compute_cleaning_summary(network_results)
  node_stats_all       <- compute_node_stats(network_results)
  edge_stats_all       <- compute_edge_stats(network_results)
  hyphae_summary       <- compute_hyphae_summary(network_results)

  write.csv(cleaning_summary_all, "cleaning_summary_hyphae.csv", row.names = FALSE)
  write.csv(node_stats_all,      "node_level_metrics_hyphae.csv", row.names = FALSE)
  write.csv(edge_stats_all,      "edge_level_metrics_hyphae.csv", row.names = FALSE)
  write.csv(hyphae_summary,      "hyphae_network_summary.csv",   row.names = FALSE)

  save_hyphae_plots(
    network_results,
    plots_dir = file.path(output_dir, "plots_json")
  )

  message("FySkel pipeline complete. Outputs written to: ", output_dir)

  invisible(list(
    network_results = network_results,
    json_files      = res$json_files,
    base_dir        = base_dir,
    output_dir      = output_dir
  ))
}
