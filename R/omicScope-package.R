#' @importFrom utils packageDescription
.onAttach <- function(libname, pkgname) {
    pkgVersion <- packageDescription(pkgname, fields = "Version")

    packageStartupMessage("Welcome to use omicScope for RNA-seq analysis.")
    packageStartupMessage(paste("The version of omicScope:",
                                pkgVersion,
                                "\nAny advice or suggestions please contact with me: 3219030654@stu.cpu.edu.cn.",
                                sep = " "))
}
