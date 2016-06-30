# env-setup
get.base.path <- function() {
  return('.')
}

get.path <- function(someRelativePath) {
  filePath <- file.path(get.base.path(), someRelativePath)

  if(!file.exists(filePath))
    dir.create(filePath, recursive = T)

  filePath
}
