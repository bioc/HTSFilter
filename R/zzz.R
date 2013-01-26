.onAttach <- function(libname, pkgname)
{
	if(.Platform$OS.type=="windows" && 
		.Platform$GUI=="Rgui") {

		winMenuAddItem("Vignettes", "HTSFilter",
		"shell.exec(system.file(\"doc\",\"HTSFilter.pdf\", package=\"HTSFilter\"))")

	}
}