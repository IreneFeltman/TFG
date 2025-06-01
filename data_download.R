# Cargar el instalador/cargador de paquetes
source("install_packages.R")
library(GEOquery)

# Descargar el objeto GEO Series (GSE) con todas las plataformas
GSE_ID <- "GSE224447"
gse_list <- getGEO(GSE_ID, GSEMatrix = TRUE)

# Si hay múltiples plataformas, gse_list será una lista
if (!is.list(gse_list)) gse_list <- list(gse_list)

for (gse in gse_list) {
  platform <- annotation(gse)
  
  # Guardar matrices de expresión
  meth_matrix <- exprs(gse)
  out_file_meth <- paste0(GSE_ID, "-", platform, "_methylation_matrix.csv")
  write.csv(meth_matrix, file = out_file_meth, row.names = TRUE)
  cat("Saved methylation matrix for platform", platform, "as", out_file_meth, "\n")
  
  # Guardar matrices de metadatos (donde se relacionan las muestras con sus características)
  # Será especialmente útil para visualizar los datos agrupados
  # Heterocrónico vs Isocrónico
  pheno_data <- pData(gse)
  out_file_pheno <- paste0(GSE_ID, "-", platform, "_sample_metadata.csv")
  write.csv(pheno_data, file = out_file_pheno, row.names = TRUE)
  cat("Saved sample metadata for platform", platform, "as", out_file_pheno, "\n")
}
