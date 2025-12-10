#' Create spatial networks use `sfnetworks` for use with `tidygraphs`
#'
#' @description
#' Load, create spatial objects, clean and create spatial networks from HyPhy raw data and
#'
eRgot <- function(x){

  # import the json data.
  raw_data <- jsonlite::fromJSON(x)

  # create a spatial object of the data
  spat_obs <- R_get(raw_data)

  # map the input data set.
  unclean_maps <- unclean_mapR(spat_obs)

  # create a spatial network
  organ_types <- c('sclerotiaPrimordia', 'network')
  network_sclerotia <- networkR(spat_obs, organ_types[1])
  network_hyphae <- networkR(spat_obs, organ_types[2])

  # map the clean data set.
  clean_maps <- clean_mapper(
    network_sclerotia[["clean_net"]],
    network_hyphae[["clean_net"]])
  # return all objects.
  return(
    list(
      raw_data = raw_data, # return raw json information
      sf_unclean = spat_obs, # the raw json to first spatial object

      unclean_maps = unclean_maps, # maps of the raw spatial object
      clean_maps = clean_maps, # maps of the cleaned sf NETWORK
      network_plots = network_plotR(unclean_maps, clean_maps),

      #finally the actual networks.
      network_sclerotia = network_sclerotia,
      network_hyphae = network_hyphae
    )
  )
}




#' Convert a HyPhy skeletonized JSON object into an sf object
#'
#' @description
#' Converts an object in the HyPhy 'SkeletonizationTool\Skeletons\Calculations' dir
#' into an sf object. Allowing for downstream analysis in sf networks and facilitates
#' plotting with ggplot2.
#' @param x a HyPhy object read in from jsonlite::fromJSON
#' @export
R_get <- function(x){

  organ_types <- c('sclerotiaPrimordia', 'network')
  spat_obs <- vector(mode = 'list', length = length(organ_types))
  names(spat_obs) <- organ_types

  for (i in seq_along(organ_types)){

    points_json <- x[[organ_types[i]]][['vectorized']][['points']]
    lines_json <- x[[organ_types[i]]][['vectorized']][['lines']]

    spat_obs[[i]] <- bicycle_day(points_json, lines_json)
  }

  return(spat_obs)
}


#' Construct a spatial node–edge network from HyPhy JSON output
#'
#' @description
#' Converts HyPhy-style JSON lists of points and lines into an sf spatial
#' representation containing nodes and edges for downstream network analysis.
#'
#' @param points_json A data frame or matrix of point coordinates (x, y)
#'   extracted from HyPhy JSON.
#' @param lines_json A list of integer vectors, each containing the indices
#'   of points that form a single edge.
#'
#' @return A list with:
#' \describe{
#'   \item{nodes}{sf POINT object of node coordinates}
#'   \item{edges}{sf LINESTRING object of reconstructed edges}
#' }
#'
#' @export
#'


bicycle_day <- function(points_json, lines_json){

  # load the nodes
  nodes <- setNames(
    data.frame(
      0:(nrow(points_json)-1),
      points_json
    ), c('PointIndx', 'x', 'y')
  ) |>
    sf::st_as_sf(coords = c('x', 'y'), remove = F)

  # read in the edges
  edges <- data.frame(  # note 0 based indexing
    LineIndx = rep((seq_along(lines_json)-1), lengths(lines_json)),
    PointIndx = unlist(lines_json)
  )

  # each edges is composed of points, we will join the point indices back to the
  # edges so we can reconstruct the edges
  edges <- merge(
    edges,
    nodes[,c('PointIndx', 'geometry')],
    by = 'PointIndx',
    all.x = TRUE) # left join in base R syntax...

  # create an explicitly spatial object, group by the LINE ID, to capture all nodes.
  edges <- edges |> # now use the points
    dplyr::group_by(LineIndx) |>
    dplyr::summarize(geometry = sf::st_combine(geometry)) |>
    sf::st_as_sf() |>
    sf::st_cast("LINESTRING")

  return(
    list(
      nodes = nodes,
      edges = edges
    )
  )
}


#' Plot an imported HyPhy object
#'
#' @description
#' Plots HyPhy imported data
#' @param x an sf object, the output from eRgot, containing both sclerotiaPrimordia and hyphal networks
#' @export
unclean_mapR <- function(x){

  sp <- 'sclerotiaPrimordia'
  net <- 'network'

  # Base plot theme
  base_plot <- function() {
    list(
      ggplot2::theme_minimal(),
      ggplot2::xlim(0, 1),
      ggplot2::ylim(0, 1)
    )
  }

  # Common geoms for each type
  sp_geoms <- list(
    ggplot2::geom_sf(data = x[[sp]][['edges']], color = '#360A14'),
    ggplot2::geom_sf(data = x[[sp]][['nodes']],
                     fill = NA, shape = 21, size = 1, alpha = 0.6)
  )

  net_geoms <- list(
    ggplot2::geom_sf(data = x[[net]][['edges']], color = '#0F7173'),
    ggplot2::geom_sf(data = x[[net]][['nodes']],
                     fill = NA, shape = 21, size = 1, alpha = 0.6)
  )

  list(
    sclerotiaPrimordia_plot = ggplot2::ggplot() + sp_geoms +
      ggplot2::labs(title = paste('original', sp)) + base_plot(),

    hyphae_plot = ggplot2::ggplot() + net_geoms +
      ggplot2::labs(title = paste('original hyphae')) + base_plot(),

    fullSlide_plot = ggplot2::ggplot() + sp_geoms + net_geoms +
      ggplot2::labs(title = paste('original', sp, '&', 'hyphae')) + base_plot()
  )

}


#' Convert a spatial HyPhy object into a network, including network cleaning and simplification
#'
#' @description
#' Given HyPhy data which has been converted into a spatial object via `HyPhy` create a spatial network,
#' and simplify it to remove internal nodes along un-branched edges.
#' @param x an sf object, the output from eRgot, containing both sclerotiaPrimordia and hyphal networks
#' @param organ character, the organ type to process. Character, one of 'sclerotiaPrimordia' or 'hyphae'
#' @export
networkR <- function(x, organ){

  # decrease precision of coordinates for faster spatial operations
  nodes_sf <- x[[organ]][['nodes']] |>
    sf::st_set_precision(1e6)

  edges_sf <- x[[organ]][['edges']] |>
    sf::st_set_precision(1e6)

  # Extract start and end points from edges
  edge_starts <- sf::st_line_sample(edges_sf, sample = 0) |> sf::st_cast("POINT")
  edge_ends <- sf::st_line_sample(edges_sf, sample = 1) |> sf::st_cast("POINT")

  # for each edge, find the node which is closest to it's start and end.
  edges_sf$from = sf::st_nearest_feature(edge_starts, nodes_sf)
  edges_sf$to = sf::st_nearest_feature(edge_ends, nodes_sf)
  net <- sfnetworks::sfnetwork(nodes_sf, edges_sf, directed = FALSE)

  # now create a spatial network
  net_cleaned <- net |>
    tidygraph::activate("edges") |>
    dplyr::filter(!sf::st_is_empty(geometry)) |>
    tidygraph::convert(to_spatial_smooth, summarise_attributes = "first")

  # remove internal nodes, i.e. those associated with a 'bend' which is connected
  # by linestrings, which must be staight (e.g. lines) and replace the connections
  # with MULTILINESTRINGS, connections of linestrings which allow for an individual
  # hyphae to bend.
  net_cleaned <- net_cleaned |>
    tidygraph::activate("edges") |>
    sf::st_as_sf() |>
    sf::st_cast("MULTILINESTRING") |>  # Cast to MULTILINESTRING first
    sf::st_line_merge() |>
    sfnetworks::as_sfnetwork(directed = FALSE)

  cleaning_summary <- clean_summary(net, net_cleaned)

  return(
    list(
      original_net = net,
      clean_net = net_cleaned,
      cleaning_summary = cleaning_summary
    )
  )

}

## summary table on difference between the unclean and cleaned data sets.
clean_summary <- function(net, net_cleaned){
  ### difference between the number of uncleaned and cleaned edges.
  ct_raw_edges <- net %>%
    tidygraph::activate("edges") %>%
    tibble::as_tibble() %>%
    nrow()
  ct_clean_edges <- net_cleaned %>%
    tidygraph::activate("edges") %>%
    tibble::as_tibble() %>%
    nrow()

  ### difference between the number of uncleaned and cleaned nodes
  ct_raw_nodes <- net %>%
    tidygraph::activate("nodes") %>%
    tibble::as_tibble() %>%
    nrow()
  ct_clean_nodes <- net_cleaned %>%
    tidygraph::activate("nodes") %>%
    tibble::as_tibble() %>%
    nrow()

  ## want the difference between cleaned-edges:nodes uncleaned-edges:nodes to be minimized.
  data.frame(
    portion = c('Edges', 'Nodes'),
    original = c(ct_raw_edges, ct_raw_nodes),
    processed = c(ct_clean_edges, ct_clean_nodes),
    difference = c(
      ct_raw_edges - ct_clean_edges,
      ct_raw_nodes - ct_clean_nodes
    ),
    node_edge_ratio = c(
      ct_raw_nodes / ct_raw_edges,
      ct_clean_nodes / ct_clean_edges

    )
  )
}

clean_mapper <- function(network_sclerotia, network_hyphae){

  sp <- 'sclerotiaPrimordia'
  net <- 'network'

  sp_nodes = network_sclerotia |> tidygraph::activate("nodes") |> sf::st_as_sf()
  sp_edges = network_sclerotia |> tidygraph::activate("edges") |> sf::st_as_sf()

  hyp_nodes = network_hyphae |> tidygraph::activate("nodes") |> sf::st_as_sf()
  hyp_edges = network_hyphae |> tidygraph::activate("edges") |> sf::st_as_sf()

  # Base plot theme
  base_plot <- function() {
    list(
      ggplot2::theme_minimal(),
      ggplot2::xlim(0, 1),
      ggplot2::ylim(0, 1)
    )
  }
  # Common geoms for each type
  sp_geoms <- list(
    ggplot2::geom_sf(data = sp_edges, color = '#360A14'),
    ggplot2::geom_sf(data = sp_nodes,
                     fill = NA, shape = 21, size = 1, alpha = 0.6)
  )

  net_geoms <- list(
    ggplot2::geom_sf(data = hyp_edges, color = '#0F7173'),
    ggplot2::geom_sf(data = hyp_nodes,
                     fill = NA, shape = 21, size = 1, alpha = 0.6)
  )

  list(
    sclerotiaPrimordia_plot = ggplot2::ggplot() + sp_geoms +
      ggplot2::labs(title = paste('clean', sp)) + base_plot(),

    hyphae_plot = ggplot2::ggplot() + net_geoms +
      ggplot2::labs(title = paste('clean hyphae')) + base_plot(),

    fullSlide_plot = ggplot2::ggplot() + sp_geoms + net_geoms +
      ggplot2::labs(title = paste('clean', sp, '&', 'hyphae')) + base_plot()
  )
}


network_plotR <- function(unclean_maps, clean_maps){

  pls <- list(
    unclean_maps[['sclerotiaPrimordia_plot']],
    unclean_maps[['hyphae_plot']],
    unclean_maps[['fullSlide_plot']],
    clean_maps[['sclerotiaPrimordia_plot']],
    clean_maps[['hyphae_plot']],
    clean_maps[['fullSlide_plot']]
  )

  patchwork::wrap_plots(pls) +
    patchwork::plot_layout(nrow = 2) +
    patchwork::plot_annotation(
      title = 'Networks',
      theme = ggplot2::theme(plot.title = ggplot2::element_text(size = 16, face = 'bold'))
    )
}
