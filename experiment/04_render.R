#!/usr/bin/env Rscript
#
# 04_render.R -- build the report, with or without pandoc.
#
#   EXP_CONFIG=quick Rscript experiment/04_render.R
#
# rmarkdown::render() needs pandoc, which ships with RStudio but is absent on a
# plain R installation and on most cluster nodes. When it is missing we fall
# back to knitr::knit(), which produces GitHub-readable markdown plus a figure
# directory and needs nothing extra. The Rmd defines a `params` fallback so the
# same source works under both paths.

setwd(file.path(getwd(), "experiment"))
config <- Sys.getenv("EXP_CONFIG", "quick")
Sys.setenv(EXP_CONFIG = config)

if (rmarkdown::pandoc_available("1.12.3")) {
  message("pandoc found; rendering HTML")
  rmarkdown::render("04_report.Rmd", params = list(config = config),
                    output_file = sprintf("04_report_%s.html", config),
                    quiet = TRUE)
  message("wrote experiment/04_report_", config, ".html")
} else {
  message("pandoc not available; knitting to markdown instead")
  out <- sprintf("04_report_%s.md", config)
  knitr::opts_chunk$set(fig.path = sprintf("figures/%s-", config))
  knitr::knit("04_report.Rmd", output = out, quiet = TRUE)
  message("wrote experiment/", out, " (figures in experiment/figures/)")
}
