# TFG
Código del Análisis epigenético
Realizado en R, para ejecutar cada script se presionan las teclas "Ctrl + Shift + S" una vez se haya abierto en la ventana de RStudio. 

Para realizar el análisis se tienen que seguir los siguientes pasos:

1. Instalar los paquetes ejecutando el script "install_packages.R"
   
2. Descargar los datos de GEO mediante el script de "data_download.R", al guardarse la matriz y metadatos de todas las plataformas puede tomar un poco de tiempo
   
3. En caso de que se quiera hacer el analisis mediante el reloj del modelo de Horvath, se ejecuta "clock_script.R", en caso de que se quiera realizar el analisis
   del reloj de PC-Clocks sera necesario abrir primero "download_pc_clocks.R" y una vez descargados los archivos se ejecutara "pc_clocks_analysis.R". El resultado de
   con las predicciones aparecera en formato csv en la carpeta donde se haya ejecutado el codigo.


Para visualizar exclusivamente los datos y no realizar el analisis, he subido los resultados de los analisis, adjuntados en formato csv con los nombres de 
"GSE22447-GPL28271_DNAmAge_predictions.csv", "GSE22447-GPL28271_pc_clocks_predictions.csv" y "GSE22447-GPL28271_sample_metadata.csv".
Simplemente es necesario ejecutar "install_packages.R" en caso de que no se haya realizado el analisis previamente y, posteriormente, se abre 
bien "visualizacion.R" para ver los resultados del reloj de Horvath o "visalizacion_pc_clocks.R" para ver los resultados del modelo PC-Clock.

Se han adjuntado todas las graficas resultantes del analisis en formato png. 
