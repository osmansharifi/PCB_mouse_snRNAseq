#' ============================================================================
#' Logging Utilities for PCB snRNA-seq Analysis
#' ============================================================================

#' Log Message with Timestamp
log_message <- function(message, level = "INFO", verbose = TRUE) {
  
  if (!verbose) return(invisible(NULL))
  
  timestamp <- format(Sys.time(), "%Y-%m-%d %H:%M:%S")
  
  icons <- list(
    INFO = "ℹ",
    SUCCESS = "✓",
    WARNING = "⚠",
    ERROR = "✗"
  )
  
  if (!level %in% names(icons)) {
    level <- "INFO"
  }
  
  icon <- icons[[level]]
  formatted_msg <- sprintf("[%s] %s %s", timestamp, icon, message)
  
  cat(formatted_msg, "\n")
  invisible(NULL)
}

#' Log Progress
log_progress <- function(current, total, prefix = "") {
  
  percent <- (current / total) * 100
  bar_length <- 20
  filled <- round((current / total) * bar_length)
  
  bar <- paste0(
    "[",
    paste0(rep("=", filled), collapse = ""),
    paste0(rep("-", bar_length - filled), collapse = ""),
    "]"
  )
  
  msg <- sprintf("%s %s %d/%d (%.1f%%)", prefix, bar, current, total, percent)
  cat("\r", msg, sep = "")
  
  if (current == total) cat("\n")
  invisible(NULL)
}

#' Log DEG Summary
log_deg_summary <- function(deg_results, cell_type = NULL) {
  
  n_deg <- nrow(deg_results)
  n_up <- sum(deg_results$avg_log2FC > 0, na.rm = TRUE)
  n_down <- sum(deg_results$avg_log2FC < 0, na.rm = TRUE)
  
  if (!is.null(cell_type)) {
    msg <- sprintf(
      "%s: %d DEGs (%d up-regulated, %d down-regulated)",
      cell_type, n_deg, n_up, n_down
    )
  } else {
    msg <- sprintf(
      "%d DEGs (%d up-regulated, %d down-regulated)",
      n_deg, n_up, n_down
    )
  }
  
  log_message(msg, level = "SUCCESS")
}

#' Log File Save
log_file_save <- function(file_path, n_rows = NULL, description = NULL) {
  
  if (!is.null(n_rows) & !is.null(description)) {
    msg <- sprintf(
      "Saved %s (%d rows) to: %s",
      description, n_rows, file_path
    )
  } else if (!is.null(n_rows)) {
    msg <- sprintf("Saved %d rows to: %s", n_rows, file_path)
  } else if (!is.null(description)) {
    msg <- sprintf("Saved %s to: %s", description, file_path)
  } else {
    msg <- sprintf("Saved to: %s", file_path)
  }
  
  log_message(msg, level = "SUCCESS")
}

#' Print Session Info
log_session_info <- function() {
  
  session_info <- sessionInfo()
  
  cat("\n")
  cat(paste0(rep("=", 60), collapse = ""), "\n")
  cat("SESSION INFORMATION\n")
  cat(paste0(rep("=", 60), collapse = ""), "\n")
  print(session_info)
  cat(paste0(rep("=", 60), collapse = ""), "\n\n")
}
