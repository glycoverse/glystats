# Check if required packages are available
.check_pkg_available <- function(pkg, install_hint = NULL) {
  if (!requireNamespace(pkg, quietly = TRUE)) {
    if (is.null(install_hint)) {
      install_hint <- stringr::str_glue("install.packages('{pkg}')")
    }
    cli::cli_abort(c(
      "Package {.pkg {pkg}} is required for this analysis.",
      "i" = "Install it with: {.code {install_hint}}"
    ))
  }
}
