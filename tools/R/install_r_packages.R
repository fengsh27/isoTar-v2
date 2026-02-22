install_if_missing <- function(pkg) {
    if (!requireNamespace(pkg, quietly = TRUE)) {
        tryCatch(
            install.packages(pkg, dependencies = FALSE),
            error = function(e) message("Skipping ", pkg, ": ", e$message)
        )
    }
}

install_if_missing("rsvg")
install_if_missing("V8")
install_if_missing("gsubfn")
install_if_missing("optparse")
install_if_missing("png")

if (!requireNamespace("BiocManager", quietly = TRUE))
    install.packages("BiocManager", dependencies = FALSE)
BiocManager::install("RamiGO", ask = FALSE, update = FALSE)
BiocManager::install("RCytoscape", ask = FALSE, update = FALSE)
BiocManager::install("BiocGenerics", ask = FALSE, update = FALSE)
BiocManager::install("graph", ask = FALSE, update = FALSE)
