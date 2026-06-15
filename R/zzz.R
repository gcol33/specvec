.onLoad <- function(libname, pkgname) {
  ## cover scales
  register_cover_scale("percent", .scale_percent)
  register_cover_scale("proportion", .scale_proportion)
  register_cover_scale("braun_blanquet", .scale_braun_blanquet)

  ## weightings
  register_weighting("counts", .w_counts,
                     input = "species_species", supports_abundance = FALSE)
  register_weighting("ppmi", .w_ppmi,
                     input = "species_species", supports_abundance = FALSE)
  register_weighting("abundance_pmi", .w_abundance_pmi,
                     input = "species_species", supports_abundance = TRUE)
  register_weighting("chi_square", .w_chi_square,
                     input = "plot_species", supports_abundance = FALSE)
  register_weighting("clr", .w_clr,
                     input = "species_species", supports_abundance = TRUE)

  ## factorizations
  register_factorization("eigen", .f_eigen, kind = "matrix", accepts = "sym")
  register_factorization("svd", .f_svd, kind = "matrix", accepts = "implicit")
  register_factorization("glove", .f_glove, kind = "objective", accepts = "counts")

  ## method presets
  register_method("ca", weighting = "chi_square", factorization = "svd",
                  input = "plot_species", native_output = "species",
                  supports_abundance = FALSE)
  register_method("pmi", weighting = "ppmi", factorization = "eigen",
                  input = "species_species", native_output = "species",
                  supports_abundance = FALSE)
  register_method("abund_pmi", weighting = "abundance_pmi", factorization = "eigen",
                  input = "species_species", native_output = "species",
                  supports_abundance = TRUE)
  register_method("glove", weighting = "counts", factorization = "glove",
                  input = "species_species", native_output = "species",
                  supports_abundance = FALSE)
  register_method("clr", weighting = "clr", factorization = "eigen",
                  input = "species_species", native_output = "species",
                  supports_abundance = TRUE)

  invisible(NULL)
}
