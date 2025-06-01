# clock_script.R

# Aceptar el archivo de matriz de metilación como argumento de línea de comandos
args <- commandArgs(trailingOnly = TRUE)
meth_file <- if (length(args) > 0) args[1] else "GSE224447-GPL28271_methylation_matrix.csv"  

# Cargar el instalador de paquetes
source("install_packages.R")

# Cargar la matriz
meth_matrix <- read.csv(meth_file, check.names = FALSE)
# Renombrar la primera celda a var
if (colnames(meth_matrix)[1] %in% c("X", "", "V1")) {
  colnames(meth_matrix)[1] <- "var"
}
# Para comprobar que el nombre de los CpG corresponde a cg_____
# fue especialmente útil para comprobar si la matriz habia cargado correctamente
# y en algunas pruebas con datos de otros GEO tienen otro formato para nombrar las islas
if (!all(grepl("^cg", as.character(meth_matrix[[1]])[1:5]))) {
  stop("Error, los valores de identificación para las islas tienen otro formato, los valores que aparecen son: ", paste(head(meth_matrix[[1]]), collapse=", "))
}
meth_df <- meth_matrix

# Para no tener en cuenta los valores vacíos o que pongan NA
meth_matrix <- meth_matrix[rowSums(is.na(meth_matrix) | meth_matrix == "") < ncol(meth_matrix), ]

# Código de debbugeo de la línea anterior, ya que dio error en varias ocasiones 
if (nrow(meth_matrix) == 0 || ncol(meth_matrix) == 0) {
  cat("Limpieza de valores vacíos", meth_file, "\n")
  quit(save = "no")
}

# Asegurar que los IDs de CpG no tienen espacios, en ocasiones al descargar otros datos de GEO, los IDs tienen espacios
meth_df[[1]] <- trimws(as.character(meth_df[[1]]))
meth_df[[1]] <- sub("_.*", "", meth_df[[1]])

# Comprobante de los nombres, fue también muy útil para asegurarme de que los datos se cargaron correctamente
print("Nombre primeras columnas:")
print(colnames(meth_df))
print("Nombre primeras filas:")
print(head(meth_df))
print("Primeros valores de la primera columna:")
print(head(meth_df[[1]]))

# Cargado de la tabla de coeficientes (Horvath)
coef_file <- "TableS.csv"
datCoef <- read.csv(coef_file, sep = ";", check.names = FALSE)

# Siguiendo las instrucciones del árticulo, se seleccionan solo las columnas de coeficientes de ratón
datCoef <- datCoef[, c(1, 10:17)]

# Base del análisis, se utiliza la función multivariatePredictorCoef para predecir la edad
data.table::setDTthreads(1) 
multivariatePredictorCoef <- function(dat0, datCOEF, imputeValues = FALSE) {
  datout <- data.frame(matrix(NA, nrow = ncol(dat0) - 1, ncol = ncol(datCOEF) - 1))
  match1 <- match(datCOEF[-1, 1], dat0[, 1])
  if (sum(!is.na(match1)) == 0) stop("Comprobar que los datos que se utilizan son cg ID.")
  dat1 <- dat0[match1, ]
  row.names1 <- as.character(dat1[, 1])
  dat1 <- dat1[, -1]
  if (imputeValues) {
    if (!requireNamespace("impute", quietly = TRUE)) install.packages("impute")
    dat1 <- impute::impute.knn(data = as.matrix(dat1), k = 10)[[1]]
  }
  for (i in 1:ncol(dat1)) {
    for (j in 2:ncol(as.matrix(datCOEF))) {
      datout[i, j - 1] <- sum(dat1[, i] * datCOEF[-1, j], na.rm = TRUE) + datCOEF[1, j]
    }
  }
  colnames(datout) <- colnames(datCOEF)[-1]
  rownames(datout) <- colnames(dat0)[-1]
  datout <- data.frame(SampleID = colnames(dat0)[-1], datout)
  datout
}

# Aplicacion del análisis
datPredictions <- multivariatePredictorCoef(meth_df, datCoef, imputeValues = FALSE)

# Renombrar las columnas
colnames(datPredictions) <- gsub(pattern = "Coef", replacement = "DNAm", x = colnames(datPredictions))
colnames(datPredictions) <- gsub(pattern = "DNAm\\.", replacement = "Edad_Predicha_", x = colnames(datPredictions))

# Transformacion log-lineal inversa para relojes humano-ratas
F.inverse <- Vectorize(function(y, maturity) {
  if (is.na(y) | is.na(maturity)) {return(NA)}
  k <- 1.5
  if (y < 0) {(maturity + k) * exp(y) - k} else {(maturity + k) * y + maturity}
})

# Aplicar la inversa según las instrucciones
maturity <- 0.21917808
if ("Edad_Predicha_FinalHumanRatPanTissueLogLinearAge" %in% colnames(datPredictions)) {
  datPredictions$Edad_Predicha_FinalHumanRatPanTissueLogLinearAge <- F.inverse(datPredictions$Edad_Predicha_FinalHumanRatPanTissueLogLinearAge, maturity = maturity)
}
if ("Edad_Predicha_FinalHumanRatBloodTissueLogLinearAge" %in% colnames(datPredictions)) {
  datPredictions$Edad_Predicha_FinalHumanRatBloodTissueLogLinearAge <- F.inverse(datPredictions$Edad_Predicha_FinalHumanRatBloodTissueLogLinearAge, maturity = maturity)
}

# Guardar predicciones con un nombre único para cada matriz de entrada
out_file <- sub("_methylation_matrix.csv", "_DNAmAge_predictions.csv", meth_file)
write.csv(datPredictions, out_file, row.names = FALSE)

cat("Predicciones de edad DNAm guardadas como", out_file, "\n")
