# install.packages("reticulate")
library(reticulate)
py_version <- "3.9.18"
path_to_python <- reticulate::install_python(version=py_version)
path_to_python <- "/home/serrec/.pyenv/versions/3.9.18/bin/python3.9"
reticulate::virtualenv_remove("pinnEV_env")
# Create a virtual envionment 'pinnEV_env' with Python 3.9.18. Install tensorflow
# within this environment.
reticulate::virtualenv_create(envname = 'pinnEV_env',
                              python=path_to_python,
                              version=py_version)

path <- paste0(reticulate::virtualenv_root(),"/pinnEV_env/bin/python")
Sys.setenv(RETICULATE_PYTHON = "/home/serrec/.virtualenvs/pinnEV_env/bin/python")
Sys.setenv(RETICULATE_LOG_LEVEL = "DEBUG")
tf_version = "2.13.1"
reticulate::use_virtualenv("pinnEV_env", required = T)
tensorflow::install_tensorflow(method="virtualenv", envname="pinnEV_env",
                               version=tf_version) #Install version of tensorflow in virtual environment
keras::install_keras(method = c("virtualenv"), envname = "pinnEV_env",version=tf_version)

keras::is_keras_available() #Check if keras is available

#Install spektral 1.3.0 - this is for the graph convolutional neural networks. Remove all code hereafter if not necessary.
# reticulate::virtualenv_install("pinnEV_env",
#                                packages = "spektral", version="1.3.0")

intall.packages("caret")