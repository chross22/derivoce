#' Regional indices this package can derive
#'
#' The named region-scale indices, what each measures, what it needs as input,
#' and where the concept comes from. All are derived from gridded fields rather
#' than downloaded, which distinguishes them from the climate indices
#' `datamatch::attach_climate_index()` serves.
#'
#' @section Why three ways of measuring the same thing:
#' Scotian Shelf inflow to the Gulf of Maine appears here three times, because
#' the literature measures it three ways and they answer different questions.
#'
#' A **transport** across a section measures the crossing itself, and is the only
#' one that gives a direction and a flux. A **water-mass fraction** measures how
#' much of the water present came from there, which is what matters for nutrients
#' and is available even where velocities are not. A **box anomaly** measures the
#' most visible consequence, freshening, and is the most robust to poor velocity
#' data but the least specific about cause.
#'
#' They disagree in informative ways. A strong inflow with a normal salinity
#' anomaly means the arriving water was not unusually fresh, which is a real
#' finding about the upstream shelf rather than a contradiction.
#'
#' @param markdown return a markdown table instead of printing a formatted one,
#'   for pasting into documentation
#' @return a data frame of class `derivoce_index_table`, or a character string of
#'   markdown when `markdown = TRUE`
#' @examples
#' derived_indices()
#'
#' cat(derived_indices(markdown = TRUE))
#' @seealso [scotian_shelf_inflow()], [northeast_channel_inflow()],
#'   [water_mass_fraction()], [eastern_gom_salinity()]
#' @export
derived_indices <- function(markdown = FALSE) {
  entries <- list(
    list(
      name = "scotian_shelf_inflow",
      measures = "Scotian Shelf Water crossing into the Gulf of Maine at Cape Sable",
      method = "transport",
      needs = "UO, VO",
      units = "m^2/s",
      sign = "positive into the Gulf of Maine",
      source = "Feng et al. 2016; Wang et al. 2022",
      description = paste(
        "Volume transport per unit depth of the Nova Scotia Current across a",
        "fixed line off Cape Sable, where shelf water either turns into the",
        "Gulf or continues past it. Informally the Scotian Shelf crossover.",
        "Scotian Shelf Water supplies over half the Gulf's freshwater budget,",
        "and is the cold, fresh, nutrient-poor counterpart to slope water",
        "entering through the Northeast Channel. The section is fixed because",
        "an index named for a place is defined by that place.")
    ),
    list(
      name = "northeast_channel_inflow",
      measures = "slope water entering the Gulf of Maine through the Northeast Channel",
      method = "transport",
      needs = "UO, VO",
      units = "m^2/s",
      sign = "positive into the Gulf of Maine",
      source = "Ramp et al. 1985; Du et al. 2022; Silver et al. 2023",
      description = paste(
        "Volume transport per unit depth across the Northeast Channel, between",
        "Georges Bank and Browns Bank. The Channel is the Gulf's deep",
        "connection to the slope and its main source of dissolved inorganic",
        "nutrients. A separate index from the Cape Sable inflow rather than a",
        "variant of it: the two alternate episodically, and the contrast is the",
        "point. Note that the classic estimates are of the deep flow, and a",
        "surface field can run the other way. The inflow is also modulated by",
        "Gulf Stream warm-core rings, whose formation rate nearly doubled after",
        "2000, so a long record spans two regimes rather than one.")
    ),
    list(
      name = "water_mass_fraction",
      measures = "what fraction of each cell is a named water mass",
      method = "T-S endmember mixing",
      needs = "SST, SSS, and two endmembers",
      units = "fraction, 0 to 1",
      sign = "1 is entirely the first endmember",
      source = "Townsend et al. 2015",
      description = paste(
        "Projects each cell's temperature and salinity onto the mixing line",
        "between two endmembers. Used for the Gulf's deep water, a varying",
        "mixture of Labrador Slope Water and Warm Slope Water whose proportion",
        "sets the nutrient supply, and for surface water as a mixture of",
        "Scotian Shelf Water with slope water. The endmembers have no default",
        "because they vary by region, season, and year. Check the residual: a",
        "fraction is returned even for water that is not a mixture of these two",
        "masses at all.")
    ),
    list(
      name = "eastern_gom_salinity",
      measures = "freshening in the eastern Gulf of Maine, where inflow first appears",
      method = "box anomaly",
      needs = "SSS",
      units = "PSU",
      sign = "negative is fresher, indicating stronger inflow",
      source = "Grodsky et al. 2025",
      description = paste(
        "Surface salinity anomaly over the eastern Gulf of Maine, the region",
        "Scotian Shelf Water reaches first after rounding Cape Sable. Follows",
        "the approach of Grodsky et al., who built it from SMAP satellite",
        "salinity, but computes it from whatever gridded salinity is supplied,",
        "so it is not a reproduction of their series. The signal is a winter",
        "one, so an annual mean dilutes it. A box mean shows that conditions",
        "changed, not that water moved.")
    )
  )

  dictionary <- do.call(rbind, lapply(entries, function(e) {
    as.data.frame(e, stringsAsFactors = FALSE)
  }))
  rownames(dictionary) <- NULL

  if (markdown) return(index_markdown(dictionary))

  class(dictionary) <- c("derivoce_index_table", "data.frame")
  dictionary
}

#' @param x a `derivoce_index_table`
#' @param ... ignored
#' @rdname derived_indices
#' @export
print.derivoce_index_table <- function(x, ...) {
  cat("Regional indices derivoce can compute\n")
  cat(strrep("-", 70), "\n", sep = "")

  plain <- as.data.frame(x)
  table <- data.frame(name = plain$name, method = plain$method,
                      needs = plain$needs, units = plain$units)
  print.data.frame(table, row.names = FALSE, right = FALSE)

  for (i in seq_len(nrow(plain))) {
    cat("\n", plain$name[i], "\n", sep = "")
    cat("  measures: ", plain$measures[i], "\n", sep = "")
    cat("  sign:     ", plain$sign[i], "\n", sep = "")
    cat("  source:   ", plain$source[i], "\n", sep = "")
  }

  cat("\nFull descriptions: as.data.frame(derived_indices())$description\n")
  cat("Markdown table:    cat(derived_indices(markdown = TRUE))\n")
  invisible(x)
}

#' Render the index dictionary as a markdown table
#'
#' @param dictionary the plain data frame of entries
#' @return a single character string, newline separated
#' @keywords internal
index_markdown <- function(dictionary) {
  columns <- c(name = "Index", measures = "Measures", method = "Method",
               needs = "Needs", units = "Units", source = "Source")

  escape <- function(z) gsub("|", "\\|", z, fixed = TRUE)
  row <- function(cells) paste0("| ", paste(cells, collapse = " | "), " |")

  lines <- c(
    row(columns),
    row(rep("---", length(columns))),
    vapply(seq_len(nrow(dictionary)), function(i) {
      cells <- vapply(names(columns), function(column) {
        value <- escape(dictionary[[column]][i])
        if (column == "name") paste0("`", value, "()`") else value
      }, character(1))
      row(cells)
    }, character(1))
  )

  paste0(paste(lines, collapse = "\n"), "\n")
}
