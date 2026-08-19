#' Describe derived columns as EML metadata
#'
#' Builds the attribute table that Ecological Metadata Language wants, for the
#' columns this package added. Derived covariates are the hardest part of a
#' dataset to document, because their meaning lives in how they were computed
#' rather than in what was measured: `SST_grad` is degrees per kilometre by
#' central differences on a lon/lat lattice, and nothing about the column name
#' or its numbers says so. That knowledge is already in this package, so it may
#' as well be emitted in the form an archive can read.
#'
#' The result is a data frame in the shape `EML::set_attributes()` consumes.
#' derivoce does not depend on the `EML` package and does not write XML — it
#' hands over the table, and the EML package turns it into a document:
#'
#' ```r
#' attributes <- eml_attributes(env)
#' EML::set_attributes(attributes, col_classes = eml_col_classes(env))
#' ```
#'
#' @section Units it cannot know:
#' A derived unit is usually a function of the source unit — a gradient of
#' temperature is degrees per kilometre, a gradient of chlorophyll is mg/m^3 per
#' kilometre — so the source unit has to be known before the derived one can be.
#'
#' The variables `datamatch` serves are known already: the Copernicus physics
#' and biogeochemistry variables, the seafloor terrain from
#' `attach_bathymetry()`, and the climate indices from
#' `attach_climate_index()`. A workflow built on those needs no `units`
#' argument. Pass one for anything else, or to override a default.
#'
#' Anything still unresolved comes back with `unit = NA` and an
#' `attributeDefinition` that names the gap, rather than a plausible guess: an
#' archived dataset with confidently wrong units is worse than one with an
#' obvious hole in it.
#'
#' Source columns that this package did not create are described as
#' `"not derived by derivoce"` for the same reason. They are yours to document.
#'
#' @section Custom units:
#' EML validates units against a standard dictionary of 195 entries, and several
#' of the quantities here are not in it — per second, per second squared, per
#' day, and metres squared per second squared among them. Those must be declared
#' alongside the attributes or the document will not validate.
#' [eml_custom_units()] returns exactly the declarations the attribute table
#' needs, and returns none if every unit used happened to be standard.
#'
#' @param env_dat an `sf` POINT object that has been through this package
#' @param vars columns to describe, or `NULL` for every covariate column
#' @param units named character vector giving the EML unit of source columns,
#'   for example `c(TEMP_INSITU = "celsius")`. Merged over the defaults for the
#'   variables datamatch serves, which are already known
#' @return a data frame with one row per column, carrying `attributeName`,
#'   `attributeDefinition`, `measurementScale`, `domain`, `unit` and
#'   `numberType`
#' @examples
#' \dontrun{
#' env <- horizontal_gradient(env, "SST")
#' env <- eke(env)
#'
#' attributes <- eml_attributes(env, units = c(SST = "celsius"))
#' custom <- eml_custom_units(attributes)
#' }
#' @seealso [eml_custom_units()], [derived_indices()]
#' @export
eml_attributes <- function(env_dat, vars = NULL, units = NULL) {
  available <- covariate_columns(env_dat)
  vars <- vars %||% available

  missing <- setdiff(vars, available)
  if (length(missing) > 0) {
    stop("Column(s) not present: ", paste(missing, collapse = ", "),
         "\nAvailable: ", paste(available, collapse = ", "), call. = FALSE)
  }

  # datamatch's own variables are known, so a normal workflow needs no `units`
  # at all; anything passed overrides the defaults for the columns it names.
  known <- datamatch_units()
  units <- c(units, known[setdiff(names(known), names(units))])

  rows <- lapply(vars, function(v) describe_column(v, env_dat[[v]], units))
  out <- do.call(rbind, rows)
  rownames(out) <- NULL
  out
}

#' Declarations for the units EML does not know
#'
#' EML validates against a fixed dictionary, and a unit outside it must be
#' declared as a `customUnit` or the document fails validation. Several
#' quantities in this package are outside it, so this returns the declarations
#' for whichever ones an attribute table actually used.
#'
#' Each row carries what EML requires of a custom unit: an `id` matching the
#' `unit` in the attribute table, the `unitType` it belongs to, the SI unit it
#' derives from, the multiplier to reach SI, and a description.
#'
#' @param attributes a table from [eml_attributes()], or a character vector of
#'   unit names
#' @return a data frame of custom unit declarations, with no rows if every unit
#'   used is standard
#' @examples
#' \dontrun{
#' attributes <- eml_attributes(env)
#' EML::set_unitList(eml_custom_units(attributes))
#' }
#' @seealso [eml_attributes()]
#' @export
eml_custom_units <- function(attributes) {
  used <- if (is.data.frame(attributes)) attributes$unit else attributes
  used <- unique(used[!is.na(used)])

  custom <- eml_unit_registry()
  keep <- custom$id %in% used
  out <- custom[keep, , drop = FALSE]
  rownames(out) <- NULL
  out
}

#' Column classes for EML, in the order the attributes are given
#'
#' `EML::set_attributes()` wants a `col_classes` vector alongside the table, one
#' entry per attribute, drawn from "numeric", "character", "factor" and "Date".
#'
#' @param env_dat an `sf` POINT object
#' @param vars columns, or `NULL` for every covariate column
#' @return a character vector
#' @examples
#' \dontrun{
#' EML::set_attributes(eml_attributes(env), col_classes = eml_col_classes(env))
#' }
#' @seealso [eml_attributes()]
#' @export
eml_col_classes <- function(env_dat, vars = NULL) {
  vars <- vars %||% covariate_columns(env_dat)
  vapply(vars, function(v) {
    value <- env_dat[[v]]
    if (is.numeric(value)) "numeric"
    else if (is.logical(value)) "numeric"
    else if (is.factor(value)) "factor"
    else if (inherits(value, "Date")) "Date"
    else "character"
  }, character(1), USE.NAMES = FALSE)
}

#' Describe one column
#'
#' @param name the column name
#' @param value the column
#' @param units named vector of source units
#' @return a one-row data frame
#' @keywords internal
describe_column <- function(name, value, units) {
  spec <- match_derived(name, units)

  numeric_like <- is.numeric(value) || is.logical(value)
  data.frame(
    attributeName = name,
    attributeDefinition = spec$definition,
    measurementScale = if (numeric_like) "ratio" else "nominal",
    domain = if (numeric_like) "numericDomain" else "textDomain",
    unit = if (numeric_like) spec$unit else NA_character_,
    numberType = if (numeric_like) spec$number else NA_character_,
    stringsAsFactors = FALSE
  )
}

#' Match a column name against what this package produces
#'
#' Derived columns are named by convention rather than registered, so the name
#' is what identifies them. Fixed names are matched outright and suffixes are
#' matched against the stem, which also yields the source column whose unit the
#' derived unit is built from.
#'
#' @param name the column name
#' @param units named vector of source units
#' @return a list with `definition`, `unit` and `number`
#' @keywords internal
match_derived <- function(name, units) {
  fixed <- eml_fixed_registry()
  if (name %in% names(fixed)) {
    entry <- fixed[[name]]
    return(list(definition = entry$definition, unit = entry$unit,
                number = entry$number %||% "real"))
  }

  for (entry in eml_suffix_registry()) {
    hit <- regmatches(name, regexpr(entry$pattern, name))
    if (length(hit) == 0 || hit[1] != substring(name, nchar(name) -
                                                nchar(hit[1]) + 1)) next
    source <- substring(name, 1, nchar(name) - nchar(hit[1]))
    if (!nzchar(source)) next
    return(list(
      definition = paste0(entry$definition, " Derived by derivoce from `",
                          source, "`."),
      unit = entry$unit(source, units),
      number = entry$number %||% "real"
    ))
  }

  list(definition = paste0("`", name, "`: not derived by derivoce, so its ",
                           "meaning and units are not known here and should ",
                           "be described by the data author."),
       unit = source_unit(name, units), number = "real")
}

#' Unit of a source column, if the caller supplied one
#'
#' @param name the source column
#' @param units named vector of source units
#' @return an EML unit name, or `NA`
#' @keywords internal
source_unit <- function(name, units) {
  if (is.null(units) || !name %in% names(units)) return(NA_character_)
  unname(units[[name]])
}

#' A derived unit built from a source unit
#'
#' Returns `NA` when the source unit is unknown, rather than guessing. A unit is
#' the one piece of metadata that is worse to state confidently and wrongly than
#' to leave blank.
#'
#' @param source the source column name
#' @param units named vector of source units
#' @param per what the source unit is divided by, as an EML-style suffix
#' @return an EML unit name, or `NA`
#' @keywords internal
derived_unit <- function(source, units, per = NULL) {
  base <- source_unit(source, units)
  if (is.na(base)) return(NA_character_)
  if (is.null(per)) return(base)
  paste0(base, per)
}
