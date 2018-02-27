context("Testing KEGG interaction functions")

test_that("Testing KEGG function errors", {

  expect_error(tor_fetch_kegg_details(), "no.*entries provided")
  expect_message(tryCatch(tor_fetch_kegg_details("DNE"),
                          warning = function(w) {}, error = function(e) {}),
                 "querying KEGG API")
  expect_warning(tryCatch(tor_fetch_kegg_details("DNE"), error = function(e) {}),
                          "KEGG API could not process")
  expect_error(suppressWarnings(tor_fetch_kegg_details("DNE")),
               " queries did not return any results")

  expect_error(tor_fetch_kegg_pathway_details(), "no pathways data frame provided")
  expect_error(tor_fetch_kegg_pathway_details(data_frame()), "does not have the required 'kegg_pathway_id' column")
  expect_error(tor_fetch_kegg_pathway_details(data_frame(kegg_pathway_id = "DNE")), "must have the pattern")

  expect_error(tor_fetch_kegg_pathway_maps(), "no pathways data frame provided")
  expect_error(tor_fetch_kegg_pathway_maps(data_frame()), "does not have the required 'kegg_org_id' column")
  expect_error(tor_fetch_kegg_pathway_maps(data_frame(kegg_org_id=1)), "does not have the required 'kegg_pathway_id' column")
})


test_that("Test KEGG querying", {

  # species
  expect_message(orgs <- tor_fetch_kegg_species(), "querying.*species list")
  expect_true(all(c("T.number", "organism", "species") %in% names(orgs)))

  # pathways
  expect_message(pathways <- tor_fetch_kegg_pathways_for_species(org_id = "eco"),
                 "querying.*genes")
  expect_true(all(c("kegg_org_id", "uniprot_id", "kegg_gene_id", "kegg_pathway_id", "pathway") %in% names(pathways)))

  # pathway details
  expect_message(details <- tor_fetch_kegg_pathway_details(pathways[1:2,]),
                 "querying.*entries")
  expect_true(all(c("kegg_id", "pathway_map", "name", "description") %in% names(details)))

})
