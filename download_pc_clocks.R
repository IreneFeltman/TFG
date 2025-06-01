# download_pc_clocks.R

# Instalar y cargar paquetes necesarios
if (!requireNamespace("httr", quietly = TRUE)) {
    install.packages("httr")
}
library(httr)

# URLs de los archivos en GitHub
base_url <- "https://raw.githubusercontent.com/MorganLevineLab/PC-Clocks/main/"
files <- c(
    "run_calcPCClocks.R",
    "run_calcPCClocks_Accel.R",
    "cgHorvathNew.csv"
)

# Descargar cada archivo
for (file in files) {
    url <- paste0(base_url, file)
    response <- GET(url)
    
    if (status_code(response) == 200) {
        writeBin(content(response, "raw"), file)
        cat("Archivo descargado:", file, "\n")
    } else {
        cat("Error al descargar:", file, "\n")
    }
}

cat("Descarga de archivos PC-Clocks completada \n") 