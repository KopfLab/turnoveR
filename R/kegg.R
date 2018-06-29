# pathway functions ====

#' Fetch KEGG pathway maps
#'
#' Downloads basic KEGG pathway maps
#'
#' @param pathways data frame with columns kegg_org_id and kegg_pathway_id
#' @param download_folder which folder to download to (relative to working directory or absolute), will be created if it does not exist
#' @param overwrite whether to overwrite existing files
#' @param quiet parameter passed to \link[utils]{download.file}
#' @param ... additional parameters passed on to \link[utils]{download.file}
#' @return returns the unaltered pathways data frame
#' @export
tor_fetch_kegg_pathway_maps <- function(
  pathways, download_folder = "kegg", overwrite = FALSE,
  kgml = TRUE, species_png = TRUE, general_png = TRUE, quiet = TRUE, ...) {

  # safety checks
  if (missing(pathways) || !is.data.frame(pathways))
    stop("no pathways data frame provided", call. = FALSE)
  if (!"kegg_org_id" %in% names(pathways))
    stop("pathways data frame does not have the required 'kegg_org_id' column", call. = FALSE)
  if (!"kegg_pathway_id" %in% names(pathways))
    stop("pathways data frame does not have the required 'kegg_pathway_id' column", call. = FALSE)
  if (!dir.exists(download_folder)) dir.create(download_folder, recursive = TRUE)

  # base url
  png_pattern <- "http://www.genome.jp/kegg/pathway/%s/%s.png"
  kgml_pattern <- "http://rest.kegg.jp/get/%s/kgml"

  # download maps
  urls <- pathways %>%
    select(kegg_org_id, kegg_pathway_id) %>%
    unique() %>%
    mutate(
      kgml_url = sprintf(kgml_pattern, kegg_pathway_id),
      species_png_url = sprintf(png_pattern, kegg_org_id, kegg_pathway_id),
      general_png_url = sprintf(png_pattern, "map", str_replace(kegg_pathway_id, kegg_org_id, "map"))
    ) %>%
    gather(type, url, ends_with("url")) %>%
    # download everything only once (i.e. maps)
    select(-kegg_org_id) %>%
    unique() %>%
    arrange(kegg_pathway_id) %>%
    # filter depending on what is requested
    filter(type != "kgml_url" | kgml,
           type != "species_png_url" | species_png,
           type != "general_png_url" | general_png) %>%
    mutate(
      path = file.path(download_folder,
                       ifelse(type == "kgml_url",
                              str_c(basename(dirname(url)),".xml"),
                              basename(url))),
      exists = file.exists(path)
    )

  # use info
  if (!quiet) {
    glue("Info: requesting {nrow(urls)} files for {length(unique(urls$kegg_pathway_id))} pathways... ",
         "{if (overwrite) 'overwriting' else 'skipping'} {sum(urls$exists)} existing files...") %>%
      message()
  }

  # do download
  urls <- urls %>% filter(overwrite | !exists)
  if (nrow(urls) > 0) {
    urls %>%
    group_by(url, path) %>%
    do({
      download.file(unique(.$url), destfile = unique(.$path), quiet = quiet, ...)
      data_frame()
    })
  }

  return(pathways)
}


#' Fetch KEGG pathway details
#'
#' Convenience wrapper to accept a pathways data frame and focus on the most relevant columns for pathways (uses \link{tor_fetch_kegg_details} internally).
#'
#' @param pathways data frame with column kegg_pathway_id
#' @export
tor_fetch_kegg_pathway_details <- function(pathways, quiet = FALSE) {

  # safety checks
  if (missing(pathways) || !is.data.frame(pathways))
    stop("no pathways data frame provided", call. = FALSE)
  if (!"kegg_pathway_id" %in% names(pathways))
    stop("pathways data frame does not have the required 'kegg_pathway_id' column", call. = FALSE)
  if (!all(str_detect(pathways$kegg_pathway_id, "^[a-z]{2,}[0-9]{4,}")))
    glue("all pathway KEGG IDs must have the pattern '<org><pathway>' ",
         "where <org> is the lower-case organism code (e.g. 'eco' for E. coli) ",
         "and <pathway> is the pathway id (e.g. '00020' for the TCA cycle)") %>%
    stop(call. = FALSE)

  info <- tor_fetch_kegg_details(pathways$kegg_pathway_id, quiet = quiet)

  # discard unneded (or confusing) columns if they exist
  discard_cols <- c("organism")
  discard_cols <- intersect(discard_cols, names(info))
  if ( length(discard_cols) > 0)  {
    info[discard_cols] <- NULL
  }

  # arrange as makes most sense
  info %>%
    select(kegg_id, ko_pathway, pathway_map, name, class, description, everything())
}

# species functions ====

#' Fetch KEGG species
#'
#' @export
tor_fetch_kegg_species <- function(quiet = FALSE) {
  # user info
  if (!quiet)
    glue("Info: querying KEGG API for species list...") %>% message()
  tbl_df(keggList("organism"))
}

#' Fetch all KEGG pathways for a species
#'
#' @param org_id KEGG organism code (typically 3-4 characters), if not provided, will be searched based on \code{species_name} using \link{tor_fetch_kegg_species}
#' @param species_name species name or regular expression to identify the organism of interest
#' @export
tor_fetch_kegg_pathways_for_species <- function(org_id = find_by_species_name(), species_name = NULL, quiet = FALSE) {

  find_by_species_name <- function() {
    if (is.null(species_name))
      stop("species name is required to find organism code", call. = FALSE)
    org <- tor_fetch_kegg_species() %>% filter(str_detect(species, species_name))
    if (nrow(org) == 0)
      glue("could not identify species based on '{species_name}'") %>% stop(call. = FALSE)
    else if (nrow(org) > 1)
      glue("found more than one ({nrow(org)}) species that match '{species_name}'") %>% stop(call. = FALSE)
    return(org$organism)
  }

  # user info
  if (!quiet)
    glue("Info: querying KEGG API for species' genes and pathways (this may take a few seconds)...") %>% message()

  gene_pathways <-
    keggLink("pathway", org_id) %>%
    enframe("kegg_gene_id", "kegg_pathway_id")

  pathways <- keggList("pathway", org_id) %>%
    enframe("kegg_pathway_id", "pathway")

  uniprot_to_kegg <- keggConv(org_id, "uniprot") %>%
    enframe("uniprot_id", "kegg_gene_id") %>%
    mutate(uniprot_id = str_replace(uniprot_id, "up:", ""))

  species_name <- keggGet(gene_pathways$kegg_gene_id[1])[[1]]$ORGANISM[1]

  glue("Info: found {nrow(pathways)} pathways, ",
       "{length(unique(gene_pathways$kegg_gene_id))} genes, ",
       "and {nrow(gene_pathways)} gene-pathway links ",
       "for {species_name}") %>% message()

  # combine data frames
  pathways %>%
    left_join(gene_pathways, by = "kegg_pathway_id") %>%
    left_join(uniprot_to_kegg, by = "kegg_gene_id") %>%
    select(uniprot_id, kegg_gene_id, kegg_pathway_id, pathway) %>%
    arrange(kegg_gene_id) %>%
    mutate(
      kegg_org_id = org_id,
      kegg_pathway_id = str_replace(kegg_pathway_id, "path:", ""),
      pathway = str_replace(pathway, str_c(" - ", species_name), "")
    ) %>%
    select(kegg_org_id, everything())
}

# general functions =====

#' Fetch KEGG details
#'
#' Helper function to convert output of \link[KEGGREST]{keggGet} to a nested data frame that is easy to work with in the tidyverse.
#' To unnest individual nested columns and preserve others use the \code{.drop=FALSE} parameter in \link[tidyr]{unnest}
#'
#' @param kegg_id one or more KEGG identifiers. Can only do 10 requests at a time (see \link[KEGGREST]{keggGet}) so if there are more than 10, splits up into multiple queries.
#' @param unnest_single_values whether to unnest single values (values that have a none or only a single entry for all retrieve records)
#' @export
tor_fetch_kegg_details <- function(kegg_id, unnest_single_values = TRUE, ..., quiet = FALSE) {
  # safety checks
  if (missing(kegg_id)) stop("no KEGG API entries provided", call. = FALSE)


  # user info
  kegg_id <- unique(kegg_id)
  n_queries <- ceiling(length(kegg_id)/10)

  if (!quiet)
    glue("Info: querying KEGG API for {length(kegg_id)} unique KEGG entries (requires {n_queries} queries) ...") %>%
      message()

  # query KEGG API
  query_results <-
    data_frame(kegg_id = kegg_id) %>%
    mutate(
      n = row_number(),
      request_group = as.integer(ceiling(n/10))
    ) %>%
    group_by(request_group) %>%
    do({
      entries <- .$kegg_id
      info <- list()
      tryCatch(
        info <- keggGet(entries),
        error = function(e) {
          glue("KEGG API could not process the request for '{collapse(entries, sep = \"', '\")}' - skipping... ") %>%
          warning(e$message, immediate. = TRUE, call. = FALSE)
        })
      data_frame(kegg_info = info)
    }) %>%
    ungroup(request_group)

  # safety check
  if (nrow(query_results) == 0)
    stop("KEGG API queries did not return any results, cannot process further", call. = FALSE)

  # unpack list results
  info_df_wide <- unpack_lists_data_frame(
    query_results, column = kegg_info,
    unnest_single_values = unnest_single_values)

  # check for entry
  if ("entry" %in% names(info_df_wide))
    info_df_wide <- select(info_df_wide, kegg_id = entry, everything())

  # return data frame
  return(info_df_wide)
}
