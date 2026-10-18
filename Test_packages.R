# To check that all the R packages required by the simulation scripts are installed please run Test_packages.R
# The test function reports, for each package, whether it is installed and its version
# Note: extrafont is optional (it is loaded only if available)

check_packages <- function(min_aftcoop_version = "0.2.4") {
  required <- c("survival", "Hmisc", "flexsurv", "cvTools", "viridis",
                "tidyverse", "hrbrthemes", "ggpubr", "gridExtra", "survAUC",
                "systemfonts","AFTCoop")
  optional <- c("extrafont", "devtools")

  
  pkgs <- c(required, optional)
  installed <- vapply(pkgs, requireNamespace, logical(1), quietly = TRUE)
  version <- vapply(pkgs, function(p) {
    if (installed[[p]]) as.character(utils::packageVersion(p)) else NA_character_
  }, character(1))

  res <- data.frame(package = pkgs,
                    type = c(rep("required", length(required)), rep("optional", length(optional))),
                    installed = installed,
                    version = version,
                    row.names = NULL)
  print(res)

  missing_req <- required[!installed[required]]
  missing_opt <- optional[!installed[optional]]

  if (length(missing_opt) > 0) {
    message("Optional packages not installed: ", paste(missing_opt, collapse = ", "))
  }

  if (installed[["AFTCoop"]] && utils::packageVersion("AFTCoop") < min_aftcoop_version) {
    warning("AFTCoop version ", utils::packageVersion("AFTCoop"),
            " is installed; results refer to version ", min_aftcoop_version)
  }

  if (length(missing_req) > 0) {
    cran_missing <- setdiff(missing_req, "AFTCoop")
    if (length(cran_missing) > 0) {
      message("Install with: install.packages(c(\"", paste(cran_missing, collapse = "\", \""), "\"))")
    }
    if ("AFTCoop" %in% missing_req) {
      message("Install AFTCoop with: devtools::install_github(\"angeclau/AFTCoop\")")
    }
    stop("Missing required packages: ", paste(missing_req, collapse = ", "))
  }

  message("All required packages are installed.")
  invisible(res)
}

check_packages()
