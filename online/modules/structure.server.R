structure_module_environment <- environment()
structure_module_files <- c(
  "structure.config.R",
  "structure.selection.R",
  "structure.transcript.R",
  "structure.sequence.R",
  "structure.job.storage.R",
  "structure.job.runtime.R",
  "structure.module.R"
)

for (structure_module_file in structure_module_files) {
  source(
    app_path("modules", "structure", structure_module_file),
    local = structure_module_environment,
    encoding = "UTF-8"
  )
}

rm(structure_module_environment, structure_module_files, structure_module_file)
