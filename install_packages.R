# install_packages.R

# Lista de paquetes requeridos
paquetes <- c(
  "GEOquery",   # Para descargar datos GEO
  "tidyverse",  # Para manipulación y visualización de datos
  "data.table", # Para lectura/escritura rápida de datos
  "ggpubr"      # Para visualizaciones mejoradas de ggplot2 (boxplots, comparaciones estadísticas)
)

# Instalar paquetes faltantes
instalar <- paquetes[!(paquetes %in% installed.packages()[,"Package"])]
if(length(instalar)) install.packages(instalar)

# Cargar todos los paquetes
lapply(paquetes, require, character.only = TRUE) 