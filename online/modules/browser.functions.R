browser_module_environment <- environment()
browser_module_files <- c(
  "browser.core.R",
  "browser.references.R",
  "browser.intervals.R",
  "browser.transcript.mapping.R",
  "browser.transcript.table.R",
  "browser.webapp.R"
)

for (browser_module_file in browser_module_files) {
  source(
    app_path("modules", "browser", browser_module_file),
    local = browser_module_environment,
    encoding = "UTF-8"
  )
}

rm(browser_module_environment, browser_module_files, browser_module_file)
