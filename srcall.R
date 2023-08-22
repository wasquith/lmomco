setwd(dirname(rstudioapi::getSourceEditorContext()$path))

# Change to the directory containing the top-level of the package
# and then source this file to source all of the R files.
# This utility is strictly for the core developer to test things
# without ever having to actually install the package!
for(f in list.files("./R", pattern="\\.R")) source(paste0("./R/", f))
