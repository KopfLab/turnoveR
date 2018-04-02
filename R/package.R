#' @keywords internal
"_PACKAGE"

#' @importFrom rlang !! !!! sym syms quo_text quo enquo
#' @importFrom readr read_csv cols_only col_character col_double read_rds write_rds
#' @importFrom readxl read_excel
#' @importFrom tibble enframe
#' @importFrom dplyr as_data_frame tbl_df %>% filter mutate data_frame left_join select rename bind_rows arrange group_by ungroup do everything starts_with ends_with summarize
#' @importFrom glue glue collapse
#' @importFrom stringr str_c str_detect str_replace str_to_lower str_interp
#' @importFrom purrr safely map map2 map_lgl map_dbl map_int map2_int map_chr map2_chr
#' @importFrom tidyr nest unnest spread gather
#' @importFrom broom tidy glance
#' @importFrom KEGGREST keggGet keggConv keggLink keggList
#' @importFrom DiagrammeR mermaid
#' @importFrom utils URLencode
#' @importFrom RCurl getURL
NULL


#' Generate flow chart
#'
#' To render icons correctly in RMarkdown files, make sure to include \code{htmltools::tagList(rmarkdown::html_dependency_font_awesome())}.
#'
#' @export
#' @param height height of the chart
#' @param width width of the chart
#' @param ... passed on to \link[DiagrammeR]{mermaid}
tor_show_package_structure <- function(height = NULL, width = height, ...) {
  tor_prep <-
    list(
      tor_function(
        name = "tor_filter_decoy",
        tooltip = "filters so that the DECOY proteins are > 1% of the total number of proteins identified",
        exists = FALSE
      ),
      tor_function(
        name = "tor_read_svm_data_file",
        tooltip = "read in an SVM data file"
      ),
      tor_function(
        name = "tor_calculate_spectral_fit_quality",
        tooltip = "",
        exists = FALSE
      ),
      tor_function(
        name = "tor_filter_peptides_by_spectral_fit_quality",
        tooltip = ""
      ),
      tor_function(
        name = "tor_add_uniprot_info",
        tooltip = ""
      ),
      tor_function(
        name = "tor_add_abundance_info",
        tooltip = "adds the protein abundance information to the peptides including calculation of relative abundances by counts and by mass",
        exists = FALSE
      ),
      tor_function(
        name = "tor_add_metadata",
        tooltip = ""
      ),
      tor_function(
        name = "tor_add_kegg_pathways_info",
        tooltip = "optional addition of KEGG pathway information",
        exists = FALSE
      )
    )

  tor_db <-
    list(
      tor_function(
        name = "tor_fetch_uniprot_species",
        tooltip = ""
      ),
      tor_function(
        name = "tor_fetch_uniprot_proteins",
        tooltip = ""
      ),
      tor_function(
        name = "tor_fetch_kegg_species",
        tooltip = ""
      ),
      tor_function(
        name = "tor_fetch_kegg_pathways_for_species",
        tooltip = ""
      ),
      tor_function(
        name = "tor_fetch_kegg_pathway_details",
        tooltip = ""
      )
    )

  tor_calc <-
    list(
      tor_function(
        name = "tor_calculate_label_rate",
        tooltip = "calculates the labeling rates from fitting an exponential curve to either individual peptide time courses or across all peptides in a protein"
      ),
      tor_function(
        name = "tor_calculate_labeled_fraction",
        tooltip = "calculates labeled fraction from 15N and 14N areas"
      ),
      tor_function(
        name = "tor_calculate_degradation_dissipation",
        tooltip = ""
      )
    )

  tor_plot <-
    list(
      tor_function(
        name = "tor_plot_label_rate_error",
        tooltip = ""
      ),
      tor_function(
        name = "tor_plot_label_rate_hist",
        tooltip = ""
      ),
      tor_function(
        name = "tor_plot_labeling_curves",
        tooltip = ""
      ),
      tor_function(
        name = "tor_plot_on_kegg_pathway",
        tooltip = "visualizes the data on top of a KEGG pathway",
        exists = FALSE
      ),
      tor_function(
        name = "tor_plot_comparison",
        tooltip = "generates comparison plots for multiple data sets",
        exists = FALSE
      )
    )

  tor_export <-
    list(
      tor_function(
        name = "tor_export_data",
        tooltip = "",
        exists = FALSE
      )
    )

  'graph TB
  %% non-turnoveR steps
  subgraph upstream of turnoveR
  mzXML["fa:fa-file-archive-o raw data files (mzXML)"]
  unknown["fa:fa-file-o unknown output (???)"]
  psms["fa:fa-file-text-o protein abundances (psms.csv)"]
  iso["fa:fa-file-text-o peptide 15N/14N (iso.csv)"]
  svm["fa:fa-file-text-o SVM results (svm_pred_results.csv)"]
  XTST("fa:fa-calculator X!Tandem & SpectraST")
  vlad("fa:fa-bolt Vlad script?")
  massacre("fa:fa-calculator Massacre")
  vladsvm("fa:fa-bolt Vlad SVM script")
  mzXML-->XTST; XTST-->unknown; unknown-->vlad; vlad-->psms
  mzXML-->massacre; massacre-->iso; iso-->vladsvm; vladsvm-->svm
  end

  %% databases
  subgraph turnoveR: database functions
  ${db_functions}
  uniprot{"fa:fa-database uniprot"}
  uniprot--"taxon ID, name"-->tor_fetch_uniprot_species
  tor_fetch_uniprot_species--"uniprot ID, gene, mass, etc."-->tor_fetch_uniprot_proteins
  kegg{"fa:fa-database KEGG"}
  kegg-->tor_fetch_kegg_species
  tor_fetch_kegg_species-->tor_fetch_kegg_pathways_for_species
  end

  %% turnoverR data prep
  subgraph turnoverR: data preparation
  ${prep_functions}
  metadata["fa:fa-file-excel-o metadata file (xlsx)"]
  metadata-->tor_add_metadata
  psms-->tor_filter_decoy
  tor_filter_decoy-->tor_add_abundance_info
  svm-->tor_read_svm_data_file
  tor_read_svm_data_file-->tor_filter_peptides_by_spectral_fit_quality
  iso-->tor_calculate_spectral_fit_quality
  tor_calculate_spectral_fit_quality-->tor_filter_peptides_by_spectral_fit_quality
  tor_fetch_uniprot_proteins-->tor_add_uniprot_info
  tor_filter_peptides_by_spectral_fit_quality-->tor_add_uniprot_info
  tor_add_uniprot_info-->tor_add_abundance_info
  tor_add_abundance_info-->tor_add_metadata
  tor_fetch_kegg_pathway_details-->tor_add_kegg_pathways_info
  tor_add_metadata-->tor_add_kegg_pathways_info
  end

  %% turnoverR calculations
  subgraph turnoveR: calculations
  ${calc_functions}
  tor_add_kegg_pathways_info-->tor_calculate_labeled_fraction
  tor_fetch_kegg_pathways_for_species-->tor_fetch_kegg_pathway_details
  tor_calculate_labeled_fraction-->tor_calculate_label_rate
  tor_calculate_label_rate-->tor_calculate_degradation_dissipation
  growth_rate(("fa:fa-info growth rate"))
  growth_rate-->tor_calculate_degradation_dissipation
  end

  %% turnoveR visualization
  subgraph turnoveR: visualization
  ${plot_functions}
  tor_calculate_degradation_dissipation-->tor_plot_label_rate_error
  tor_calculate_degradation_dissipation-->tor_plot_label_rate_hist
  tor_calculate_degradation_dissipation-->tor_plot_labeling_curves
  tor_calculate_degradation_dissipation-->tor_plot_on_kegg_pathway
  tor_calculate_degradation_dissipation-->tor_plot_comparison
  end

  %% turnoveR export
  subgraph turnoveR: export
  ${export_functions}
  export["fa:fa-file-excel-o turnoveR export (xlsx)"]
  tor_calculate_labeled_fraction-->tor_export_data
  tor_calculate_label_rate-->tor_export_data
  tor_calculate_degradation_dissipation-->tor_export_data
  tor_export_data-->export
  end

  %% tooltips
  click mzXML callback "raw files converted from .wiff and .wiff.scan"
  click XTST callback "identify peptides from spectral information, matched against E.coli proteome database"
  click unknown callback "not sure yet what this looks like"
  click psms callback "summary of protein abundanceas across all samples"
  click massacre callback "estimates the area of 14N and 15N peaks by trying to fit to theoretical isotope substitution spectra"
  click iso callback "Massacre output files with 15N and 14N areas, one file per sample"
  ${tooltips}

  %% node styles
  classDef files fill:#5eb9f9,stroke:#000,stroke-width:2px;
  class mzXML,unknown,psms,iso,svm,protein_map,metadata,export files
  classDef programs fill:#ffa449,stroke:#000,stroke-width:2px;
  class XTST,massacre programs
  classDef functions fill:#b5e835,stroke:#000,stroke-width:2px;
  class ${existing_funcs} functions
  classDef missing fill:#ffaca3,stroke:#000,stroke-width:2px;
  class ${missing_funcs} missing

  %% edge styles
  linkStyle default stroke:black, fill:none, stroke-width: 2px
  ' %>%
    # tooltips
    str_interp(
      list(
        prep_functions = tor_prep %>% map_chr(~.x$node) %>% paste(collapse = "\n"),
        calc_functions = tor_calc %>% map_chr(~.x$node) %>% paste(collapse = "\n"),
        db_functions = tor_db %>% map_chr(~.x$node) %>% paste(collapse = "\n"),
        plot_functions = tor_plot %>% map_chr(~.x$node) %>% paste(collapse = "\n"),
        export_functions = tor_export %>% map_chr(~.x$node) %>% paste(collapse = "\n"),
        existing_funcs = c(tor_prep, tor_calc, tor_db, tor_plot, tor_export) %>%
          map_chr(~if(.x$exists) {.x$name} else {NA_character_} ) %>% na.omit() %>%
          paste(collapse = ","),
        missing_funcs = c(tor_prep, tor_calc, tor_db, tor_plot, tor_export) %>%
          map_chr(~if(!.x$exists) {.x$name} else {NA_character_} ) %>% na.omit() %>%
          paste(collapse = ","),
        tooltips = c(tor_prep, tor_calc, tor_db, tor_plot, tor_export) %>%
          map_chr(~.x$tooltip) %>% paste(collapse = "\n")
      )
    ) %>%
    # render
    mermaid(width = width, height = height)
}

# flow chart helper functions ======

# structuring functions for flow chart
make_tooltip <- function(name, tooltip, url = NULL) {
  if (is.null(url)) url <- "callback"
  if (nchar(tooltip) == 0) tooltip <- name
  if (name != tooltip)
    tooltip <- paste0(name, ":<br>", tooltip)
  sprintf('click %s "%s" "%s"', name, url, tooltip)
}
make_tor_tooltip <- function(name, tooltip = name, exists = TRUE) {
  base_url <- "https://turnover.kopflab.org/reference"
  url <- if(exists) sprintf("%s/%s.html", base_url, name) else NULL
  make_tooltip(name, tooltip, url)
}
make_tor_node <- function(name) {
  sprintf('%s("fa:fa-bolt %s")', name, name)
}

tor_function <- function(name, tooltip = name, exists = TRUE) {
  list(name = name, node = make_tor_node(name), tooltip = make_tor_tooltip(name, tooltip, exists), exists = exists)
}

