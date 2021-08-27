estimate_eqtm_completion <- function(log, interval = 60) {
  while(!any(grepl("Success!", log_file <- readLines(log)))) {
    x <- as.numeric(
      sub(".*= (.*) %$", "\\1", tail(log_file[grep("Completion", log_file)], 1))
    )

    start <- strptime(
      x = sub("^##-+ (.*) -+##$", "\\1", log_file[grep("eQTM", log_file)[[1]] - 1]),
      format = "%c",
      tz = "UTC"
    )

    chunks <- as.numeric(sub(".*: ", "", log_file[grep("chunks", log_file)]))

    time_chunk <- as.numeric(
      difftime(Sys.time(), start, units = "mins")
    ) / (x / 100 * chunks)

    message(sprintf(
      "Estimated time of completion: %s (total = %s)",
      start + as.difftime(tim = time_chunk * chunks, units = "mins"),
      format(as.difftime(time_chunk * chunks / 60, units = "hours"), digits = 4)
    ))

    Sys.sleep(interval)
  }
}
# estimate_eqtm_completion(log = "logs/16.log", interval = 60)
