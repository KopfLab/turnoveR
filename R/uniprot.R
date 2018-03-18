#' Query Uniprot Database
#' @description access uniprot database and download protein information for organism of choise
#' @param taxon the taxon ID of the organism desired for uniprot query
#' @param read_cache TRUE if ok to read information from prior download if already cached
#' @param cahce_dir directory for the cache to be stored

fetch_uniprot_proteins <- function(taxon = 83333, read_cache = TRUE, cache_dir = "cache", quiet = FALSE) {

  # safety checks
  if (missing(taxon)) stop("need to supply a taxon ID number", call. =FALSE)

  if (!dir.exists(cache_dir)) dir.create(cache_dir, recursive = TRUE)

  cache_path <- file.path(cache_dir, glue("uniprot_{taxon}_proteins.rds"))


  if (!quiet && read_cache && file.exists(cache_path)) {

    # info message that we're reading from cache
      glue("Info: Accessing cache {cache_dir} for uniprot information. To redownload
           from Uniprot server, include 'read_cache = FALSE' in function parameter.") %>%
        message()

    # read from cache
    data <- read_rds(cache_path)

  } else {

    # add some info message that this might take a while
    if (!quiet) {
      glue("Info: Accessing uniprot server for download...this may take a few seconds") %>%
        message()
    }

    # fetch from server
    query <- '
PREFIX up:<http://purl.uniprot.org/core/>
PREFIX taxon:<http://purl.uniprot.org/taxonomy/>
PREFIX skos:<http://www.w3.org/2004/02/skos/core#>
SELECT ?prot_id, ?uniprot_id, ?rec_prot_name, ?rec_ec_nr, ?gene, ?metacyc_id
WHERE {
  	?protein a up:Protein .
    # organism filter
    ?protein up:organism taxon:83333 .
#  	?protein up:organism/up:scientificName ?organism .
    # mnemonic (for mapping to massacre output)
    ?protein up:mnemonic ?prot_id .
    # recommended name
  	OPTIONAL { ?protein up:recommendedName/up:fullName ?rec_prot_name } .
    OPTIONAL { ?protein up:recommendedName/up:ecName ?rec_ec_nr } .
#   # maybe? short name (NOTE: can be multiple for the same enzyme, i.e. lead to column duplicates)
#   OPTIONAL { ?protein up:recommendedName/up:shortName ?rec_prot_short } .
    # gene
    OPTIONAL { ?protein up:encodedBy/skos:prefLabel ?gene } .
#   # maybe? amino acid sequence (NOTE: can also be multiple sequences, i.e. lead to column duplicates)
#   ?protein up:sequence/rdf:value ?seq .

    # protein id
  	BIND (STRAFTER(STR(?protein),"uniprot/") AS ?uniprot_id) .
    # metacyc (if available)
  	OPTIONAL {
     	?protein rdfs:seeAlso ?metacyc .
  		?metacyc up:database <http://purl.uniprot.org/database/BioCyc> .
      BIND(STRAFTER(STR(?metacyc), "MetaCyc:") AS ?metacyc_id) .
    	FILTER regex(?metacyc,"MetaCyc","i")
    } .
  #   #?protein up:mnemonic ?mnemonic .
  #  	#FILTER regex(?mnemonic,"^cysk","i")
}'

    qd <- SPARQL("http://sparql.uniprot.org/sparql", query)
    data <- tbl_df(qd$results)
    write_rds(x = data, path = cache_path)
  }

  return(data)
}
