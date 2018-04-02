

#' Fetch uniprot species
#'
#' Searches the uniprot data base for bacterial species based on the seacrh term.
#'
#' @param search_term taxa search term
#' @export
tor_fetch_uniprot_species <- function(search_term, quiet = FALSE) {

  if (missing(search_term))
    stop("no species search term supplied", call. = FALSE)

  query <- glue('
PREFIX up:<http://purl.uniprot.org/core/>
PREFIX taxon:<http://purl.uniprot.org/taxonomy/>
PREFIX rdfs:<http://www.w3.org/2000/01/rdf-schema#>
SELECT ?taxon_id ?taxon_name
WHERE
{
	?taxon a up:Taxon .
	?taxon up:scientificName ?taxon_name .
    # bacteria only
	?taxon rdfs:subClassOf taxon:2 .
  	BIND (STRAFTER(STR(?taxon),"taxonomy/") AS ?taxon_id) .
  	FILTER regex(?taxon_name,"[search_term]","i")
}', .open = '[', .close = ']')

  # run query
  taxa <- run_uniprot_csv_query(
    query = query,
    message = glue("taxa with '{search_term}' in the name"),
    read_cache = FALSE,
    save_file_name = glue("uniprot_last_taxa_search.csv"),
    col_types = cols(
      taxon_id = col_integer(),
      taxon_name = col_character()
    ),
    quiet = quiet)

  if (nrow(taxa) == 0) {
    glue("no taxa match the search term, please try a different search") %>%
      warning(immediate. = TRUE, call. = FALSE)
  }

  return(taxa)
}



#' Fetch uniprot information
#'
#' @param taxon the ID of the species
#' @export
tor_fetch_uniprot_proteins <- function(taxon, read_cache = TRUE, cache_dir = "cache", quiet = FALSE) {

  if (missing(taxon))
    stop("no taxon provided to search for", call. = FALSE)

  # fetch from server
  query <- glue('
PREFIX up:<http://purl.uniprot.org/core/>
PREFIX taxon:<http://purl.uniprot.org/taxonomy/>
PREFIX skos:<http://www.w3.org/2004/02/skos/core#>
SELECT ?prot_id, ?uniprot_id, ?prot_name, ?gene, ?mass
WHERE {
  	?protein a up:Protein .
    # organism filter
    ?protein up:organism taxon:[taxon] .
    ##  	?protein up:organism/up:scientificName ?organism .
    # mnemonic (for mapping to massacre output)
    ?protein up:mnemonic ?prot_id .
    # recommended protein name
  	OPTIONAL { ?protein up:recommendedName/up:fullName ?prot_name } .
    ## maybe: recommended EC number - but can have multiple
    ## OPTIONAL { ?protein up:recommendedName/up:ecName ?rec_ec_nr } .
    ## maybe? short name (NOTE: can be multiple, i.e. lead to column duplicates)
    ##   OPTIONAL { ?protein up:recommendedName/up:shortName ?rec_prot_short } .
    # gene
    OPTIONAL { ?protein up:encodedBy/skos:prefLabel ?gene } .
    ## maybe? locus name, can be multipe i.e. lead to column duplicates
    ##   OPTIONAL { ?protein up:encodedBy/up:locusName ?locus } .
    ## maybe? amino acid sequence - can be multiple, i.e. lead to column duplicates
    ##   ?protein up:sequence/rdf:value ?seq .
    # protein mass in Da
    ?protein up:sequence/up:mass ?mass_int .
    BIND(STR(?mass_int) AS ?mass) .
    # protein id
  	BIND (STRAFTER(STR(?protein),"uniprot/") AS ?uniprot_id) .
    ## maybe: metacyc (if available)
    # 	OPTIONAL {
    #    	?protein rdfs:seeAlso ?metacyc .
    # 		?metacyc up:database <http://purl.uniprot.org/database/BioCyc> .
    #     BIND(STRAFTER(STR(?metacyc), "MetaCyc:") AS ?metacyc_id) .
    #   	FILTER regex(?metacyc,"MetaCyc","i")
    #   } .
    ## filter for testing purposes
    ##   #?protein up:mnemonic ?mnemonic .
    ##  	#FILTER regex(?mnemonic,"^cysk","i")
}
', .open = '[', .close = ']')

  # run query
  proteins <- run_uniprot_csv_query(
    query = query,
    message = glue("proteins of taxon {taxon}"),
    read_cache = read_cache, cache_dir = cache_dir,
    save_file_name = glue("uniprot_{taxon}_proteins.csv"),
    col_types = cols(
      prot_id = col_character(),
      uniprot_id = col_character(),
      prot_name = col_character(),
      gene = col_character(),
      mass = col_integer()
    ),
    quiet = quiet)

  # safety checks
  reps <- proteins %>% group_by(uniprot_id) %>% mutate(n = n()) %>% filter(n > 1)
  if (nrow(reps) > 0) {
    glue("{nrow(reps)} of the {nrow(proteins)} records have more than one entry",
         " - care must be taken during merging with protein data") %>%
      warning(immediate. = TRUE, call. = FALSE)
  }
  if (nrow(proteins) == 0) {
    glue("no protein records recovered, please double check the taxon ID ({taxon}) is valid") %>%
      warning(immediate. = TRUE, call. = FALSE)
  }

  return(proteins)
}


# run a uniprot SPARQL query requesting csv data
# @param ... parameters to read_csv
run_uniprot_csv_query <- function(query, message = NULL,
                                  read_cache = TRUE, cache_dir = "cache",
                                  save_file_name = "uniprot.csv", quiet = FALSE, ...) {

  # uniprot SPARQL endpoint
  endpoint <- "http://sparql.uniprot.org/sparql"

  # caching
  if (!dir.exists(cache_dir)) dir.create(cache_dir, recursive = TRUE)
  cache_path <- file.path(cache_dir, save_file_name)

  # run query
  if (!read_cache || !file.exists(cache_path)) {

    if (!quiet) {
      if (!is.null(message))
        glue("Info: querying uniprot database for {message}... ") %>%
        message(appendLF = FALSE)
      else
        message("Info: querying uniprot database... ", appendLF = FALSE)
    }

    # retrieve data
    csv_data <- query %>%
      # URL encode query
      URLencode(reserved = TRUE) %>%
      # + replacement
      str_replace("\\+", "%2B") %>%
      # generate query
      { str_c(endpoint, "?query=", .) } %>%
      # request data (different format options)
      #getURL(httpheader = c(Accept = "application/sparql-results+xml"))
      #getURL(httpheader = c(Accept = "Accept: application/json"))
      getURL(httpheader = c(Accept = "Accept: text/csv"))

    # store in file
    csv_data %>% cat(file = cache_path)

    # user info
    if (!quiet) {
      glue("retrieved {str_count(csv_data, '\\n') - 1L} records"
         #  " - data stored in '{cache_path}'"
         ) %>%
        message()
    }

  } else if (!quiet) {
    # info about cache
    if (!is.null(message))
      glue("Info: reading uniprot {message} from cached file '{cache_path}'... ") %>%
      message()
    else
      glue("Info: reading uniprot data from cached file '{cache_path}'... ") %>%
      message()
  }

  read_csv(cache_path, ...)
}

