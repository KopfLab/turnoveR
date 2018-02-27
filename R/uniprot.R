

#' Fetch uniprot information
#'
#' @export
tor_fetch_uniprot_proteins <- function(taxon = 83333, read_cache = TRUE, cache_dir = "cache") {

  if (!dir.exists(cache_dir)) dir.create(cache_dir, recursive = TRUE)

  cache_path <- file.path(cache_dir, glue("uniprot_{taxon}_proteins.rds"))


  if (read_cache && file.exists(cache_path)) {

    # add some info message that we're reading from catch (and could do read_cache=FALSE to force new load from server)

    # read from cache
    data <- read_rds(cache_path)
  } else {

    # add some info message that this might take a while


    # fetch from server
    query <- glue('
PREFIX up:<http://purl.uniprot.org/core/>
PREFIX taxon:<http://purl.uniprot.org/taxonomy/>
PREFIX skos:<http://www.w3.org/2004/02/skos/core#>
SELECT ?prot_id, ?uniprot_id, ?rec_prot_name, ?rec_ec_nr, ?gene, ?locus, ?metacyc_id
WHERE {
  	?protein a up:Protein .
    # organism filter
    ?protein up:organism taxon:[taxon] .
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
    OPTIONAL { ?protein up:encodedBy/up:locusName ?locus } .
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
}', .open = '[', .close = ']')

    qd <- SPARQL("http://sparql.uniprot.org/sparql", query)
    data <- tbl_df(qd$results)
    write_rds(x = data, path = cache_path)
  }

  return(data)
}
