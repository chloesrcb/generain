
Sys.setenv(RETICULATE_PYTHON = Sys.which("python"))
Sys.setenv(RETICULATE_LOG_LEVEL = "ERROR")

stopifnot(py_module_available("tensorflow"))

cat("✅ TensorFlow visible depuis R\n")