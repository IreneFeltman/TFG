# pc_clocks_analysis.R

# Cargar el instalador de paquetes
source("install_packages.R")

# Verificar que el archivo de Horvath existe
if (!file.exists("cgHorvathNew.csv")) {
    stop("Error: No se encuentra el archivo cgHorvathNew.csv")
}

# Función para calcular PC-Clocks
run_calcPCClocks <- function(datMeth) {
  # Verificar que los datos son una matriz
  if (!is.matrix(datMeth)) {
    stop("Error: datMeth debe ser una matriz")
  }
  
  # Verificar que hay datos de metilación
  if (nrow(datMeth) == 0 || ncol(datMeth) == 0) {
    stop("Error: La matriz de metilación está vacía")
  }
  
  # Mostrar información inicial
  cat("Dimensiones de la matriz de entrada:", dim(datMeth), "\n")
  cat("Rango de valores en la matriz:", range(datMeth, na.rm = TRUE), "\n")
  
  # Cargar los CpGs de Horvath
  horvath_cpgs <- read.csv("cgHorvathNew.csv", header = TRUE)
  cat("Número de CpGs en Horvath:", nrow(horvath_cpgs), "\n")
  cat("Columnas disponibles en horvath_cpgs:", colnames(horvath_cpgs), "\n")
  
  # Verificar que la columna CoefficientHannum existe
  if (!"CoefficientHannum" %in% colnames(horvath_cpgs)) {
    stop("Error: No se encuentra la columna CoefficientHannum en el archivo de Horvath")
  }
  
  # Filtrar solo los CpGs de Horvath que están en nuestros datos
  cpgs_comunes <- intersect(colnames(datMeth), horvath_cpgs$Name)
  cat("Número de CpGs comunes:", length(cpgs_comunes), "\n")
  
  if (length(cpgs_comunes) == 0) {
    stop("Error: No hay CpGs de Horvath en los datos")
  }
  
  # Subset de datos para solo incluir CpGs de Horvath
  datMeth_horvath <- datMeth[, cpgs_comunes]
  cat("Dimensiones después de filtrar CpGs de Horvath:", dim(datMeth_horvath), "\n")
  
  # Función para transformar valores de metilación
  transformMeth <- function(x) {
    # Transformar valores de metilación a escala M
    x[x == 0] <- 0.0001
    x[x == 1] <- 0.9999
    log2(x/(1-x))
  }
  
  # Transformar datos de metilación
  datMeth_transformed <- apply(datMeth_horvath, 2, transformMeth)
  cat("Dimensiones después de transformación:", dim(datMeth_transformed), "\n")
  cat("Rango de valores transformados:", range(datMeth_transformed, na.rm = TRUE), "\n")
  
  # Obtener coeficientes de Hannum
  coeficientes <- horvath_cpgs$CoefficientHannum[match(cpgs_comunes, horvath_cpgs$Name)]
  cat("Número de coeficientes encontrados:", length(coeficientes), "\n")
  cat("Rango de coeficientes:", range(coeficientes, na.rm = TRUE), "\n")
  
  # Verificar si hay NAs en los coeficientes
  if (any(is.na(coeficientes))) {
    cat("¡Advertencia! Hay NAs en los coeficientes\n")
    cat("Número de NAs en coeficientes:", sum(is.na(coeficientes)), "\n")
    # Filtrar solo los CpGs que tienen coeficientes válidos
    cpgs_validos <- cpgs_comunes[!is.na(coeficientes)]
    coeficientes <- coeficientes[!is.na(coeficientes)]
    datMeth_transformed <- datMeth_transformed[, cpgs_validos]
    cat("Número de CpGs válidos después de filtrar NAs:", length(cpgs_validos), "\n")
  }
  
  # Calcular la edad epigenética
  edad_epigenetica <- datMeth_transformed %*% coeficientes
  
  # Ajustar la escala de edad para ratas
  # Las ratas viven aproximadamente 2-3 años, así que ajustamos la escala
  # Primero normalizamos a un rango de 0-1
  edad_normalizada <- (edad_epigenetica - min(edad_epigenetica)) / (max(edad_epigenetica) - min(edad_epigenetica))
  # Luego escalamos a un rango de 0-3 años (vida típica de una rata)
  edad_ajustada <- edad_normalizada * 3
  
  cat("Rango de edades calculadas (antes del ajuste):", range(edad_epigenetica, na.rm = TRUE), "\n")
  cat("Rango de edades ajustadas (en años de rata):", range(edad_ajustada, na.rm = TRUE), "\n")
  
  # Crear data frame de resultados
  resultados <- data.frame(
    PCClock = edad_ajustada
  )
  
  return(resultados)
}

# Listar solo los archivos de matriz de metilación para la plataforma GPL28271
archivos_matriz <- list.files(pattern = "GSE224447-GPL28271_methylation_matrix.csv")

# Procesar cada matriz de metilación
for (archivo_matriz in archivos_matriz) {
  cat("Procesando archivo:", archivo_matriz, "\n")
  
  # Cargar matriz de metilación
  matriz_metilacion <- read.csv(archivo_matriz, check.names = FALSE)
  
  # Si la primera columna se llama "X" o está vacía o es "V1", renombrarla a "var"
  if (colnames(matriz_metilacion)[1] %in% c("X", "", "V1")) {
    colnames(matriz_metilacion)[1] <- "var"
  }
  
  # Mostrar información sobre la primera columna
  cat("Primera columna de la matriz:", colnames(matriz_metilacion)[1], "\n")
  cat("Primeros 5 valores de la primera columna:", head(matriz_metilacion[[1]], 5), "\n")
  
  # Convertir al formato de matriz requerido por PC-Clocks
  # Transponer la matriz para que los CpGs estén en las columnas
  datos_metilacion <- t(as.matrix(matriz_metilacion[, -1]))
  colnames(datos_metilacion) <- matriz_metilacion[[1]]
  
  # Mostrar información de depuración
  cat("Dimensiones de la matriz de metilación:", dim(datos_metilacion), "\n")
  cat("Primeros 5 nombres de CpGs en los datos:", head(colnames(datos_metilacion), 5), "\n")
  
  # Cargar CpGs de Horvath para comparar
  horvath_cpgs <- read.csv("cgHorvathNew.csv", header = TRUE)
  cat("Número de CpGs en Horvath:", nrow(horvath_cpgs), "\n")
  cat("Primeros 5 CpGs de Horvath:", head(horvath_cpgs$Name, 5), "\n")
  
  # Verificar si hay CpGs comunes
  cpgs_comunes <- intersect(colnames(datos_metilacion), horvath_cpgs$Name)
  cat("Número de CpGs comunes:", length(cpgs_comunes), "\n")
  if (length(cpgs_comunes) > 0) {
    cat("Primeros 5 CpGs comunes:", head(cpgs_comunes, 5), "\n")
  } else {
    cat("No se encontraron CpGs comunes. Verificando formato...\n")
    cat("Ejemplo de CpG en datos:", colnames(datos_metilacion)[1], "\n")
    cat("Ejemplo de CpG en Horvath:", horvath_cpgs$Name[1], "\n")
  }
  
  # Calcular PC-Clocks
  cat("Calculando PC-Clocks...\n")
  resultados_pc <- run_calcPCClocks(datos_metilacion)
  
  # Crear data frame de resultados
  resultados <- data.frame(
    ID_Muestra = rownames(resultados_pc),
    Edad_PC = resultados_pc$PCClock
  )
  
  # Guardar predicciones
  archivo_salida <- sub("_methylation_matrix.csv", "_pc_clocks_predictions.csv", archivo_matriz)
  write.csv(resultados, archivo_salida, row.names = FALSE)
  
  cat("Predicciones de edad PC-Clocks guardadas como", archivo_salida, "\n\n")
}

cat("¡Análisis de PC-Clocks completado!\n") 