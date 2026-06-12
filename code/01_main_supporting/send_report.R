send_report <- function(params) {
  # Optionally email a rendered report as an attachment. Disabled by default
  # (params$send_report = FALSE); all delivery details come from params so this is
  # reusable outside any specific organisation's mail server.
  if (!isTRUE(params$send_report)) return(invisible(NULL))

  library(emayili)

  `%or%` <- function(x, default) if (is.null(x)) default else x

  Sys.sleep(5) # give a just-triggered render a moment to finish writing the file

  email <- emayili::envelope() %>%
    emayili::from(addr = params$report_from) %>%
    emayili::to(params$report_recipients) %>%
    emayili::subject(subject = params$report_subject %or% "Model run complete")

  # attach whatever report files are configured (skip any that are missing)
  for (path in params$report_attachments) {
    if (file.exists(path)) email <- emayili::attachment(email, path = path)
    else warning("send_report: attachment not found, skipping: ", path)
  }

  smtp <- emayili::server(host     = params$smtp_host,
                          port     = params$smtp_port %or% 25,
                          insecure = isTRUE(params$smtp_insecure),
                          reuse    = FALSE)
  smtp(email, verbose = TRUE)

  return(invisible(NULL))
}
