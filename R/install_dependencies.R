installDependencies <- function() {
  # Define required packages and their versions
  required_pkgs <- list(
    "kableExtra" = "1.3.4",
    "geosphere" = "1.5-10",
    "sp" = "1.4-5",
    "RandomFieldsUtils" = "1.2.5",
    "RandomFields" = "3.3.14"
  )
  
  for (pkg in names(required_pkgs)) {
    required_version <- required_pkgs[[pkg]]
    
    # Check if package is installed and meets version requirement
    if (!requireNamespace(pkg, quietly = TRUE)) {
      
      # Construct the path to the local archive
      archive_path <- file.path("./inst/archives", paste0(pkg, "_", required_version, ".tar.gz"))
      
      if (file.exists(archive_path)) {
        message(paste("Installing", pkg, "from local archive..."))
        install.packages(archive_path, repos = NULL, type = "source")
      } else {
        message(paste("Installing", pkg, "from CRAN..."))
        tryCatch({
          install.packages(pkg)
        }, error = function(e) {
          message(paste("Failed to install", pkg, ":", e$message))
        })
      }
    } else {
      message(paste(pkg, "is already installed and up-to-date."))
    }
  }
}
