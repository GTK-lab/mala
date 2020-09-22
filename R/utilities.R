#' Resolve file or URL
#' @keywords internal
#'
#' @param path File path or URL.
#'
#' @details If \code{path} is a local file path, assert that it
#'   exists, then just return it. If \code{path} is a URL, see if it
#'   has been downloaded before by checking if the file
#'   \code{basename(path)} exists. If it exists, just return
#'   \code{basename(path)}. Otherwise, download the file at
#'   \code{path} to \code{basename(path)} and return
#'   \code{basename(path)}.
#'
#' @importFrom assertthat assert_that is.string
#' @importFrom RCurl url.exists
#' @importFrom utils download.file
download_if_url <- function(path) {
  assert_that(is.string(path))

  # If path is a local file path...
  if (!url.exists(path)) {
    assert_that(file.exists(path))
    return(path)
  }

  # Otherwise, first check if it had already been downloaded.
  if (file.exists(basename(path))) return(basename(path))

  # If not already downloaded, download it.
  download.file(path, basename(path), "auto")

  basename(path)
}

#' Warn if current bedtools version is older than recommended
#' @keywords internal
#'
#' @param major Numeric 'x' part of x.y.z.
#' @param minor Numeric 'y' part of x.y.z.
#' @param patch Numeric 'z' part of x.y.z.
#'
#' @return TRUE if current version < major.minor.patch, FALSE otherwise.
#'
#' @details If TRUE, this function produces a warning using the warning function.
#'
#' @importFrom assertthat assert_that is.count
#' @importFrom magrittr %>%
#' @importFrom utils compareVersion
warn_if_outdated_bedtools <- function(major=2, minor=28, patch=0) {
  assert_that(is.count(major))
  assert_that(is.count(minor))
  assert_that(is.count(patch))

  version <- get_bedtools_version()
  # Take advantage of an already implemented version comparison function.
  version_current <- version %>%
    { sprintf("%d.%d-%d", .[1], .[2], .[3]) }
  version_recommended <- sprintf("%d.%d-%d", major, minor, patch)
  is_outdated <- compareVersion(version_current, version_recommended) == -1

  if (is_outdated) {
    paste0(
        "Using an outdated version of bedtools (v%d.%d.%d). ",
        "v%d.%d.%d or higher is recommended.") %>%
      sprintf(version[1], version[2], version[3], major, minor, patch) %>%
      warning()
  }

  is_outdated
}

#' Get the current version of bedtools
#' @keywords internal
#'
#' @param command Command string which produces a bedtools version string.
#'
#' @return Length-three numeric vector of major, minor, and patch versions
#'   respectively.
#'
#' @importFrom assertthat assert_that is.string
#' @importFrom magrittr %>%
#' @importFrom stringr str_extract_all
get_bedtools_version <- function(command="bedtools --version") {
  assert_that(is.string(command))

  command %>%
    system(intern=TRUE) %>%
    str_extract_all("\\d+") %>%
    unlist() %>%  # str_extract_all gives length-1 list of length-3 vector
    as.numeric()
}
