# Cargar las librerías necesarias
library(ggplot2)
library(dplyr)
library(readr)
library(rstatix)
library(ggpubr)

# Leer los archivos CSV
predictions <- read_csv("GSE224447-GPL28271_DNAmAge_predictions.csv")
metadata <- read_csv("GSE224447-GPL28271_sample_metadata.csv")

# Renombrar columna para que coincida con merge
colnames(predictions)[1] <- "geo_accession"

# Unir datos
data <- merge(predictions, metadata, by = "geo_accession")

# Crear un campo "group" personalizado con base en el título
data$group <- case_when(
  grepl("YoungHetRec", data$title) ~ "Heterochronic Recovery",
  grepl("YoungHet", data$title) ~ "Heterochronic",
  grepl("YoungIsoRec", data$title) ~ "Isochronic Recovery",
  grepl("YoungIso", data$title) ~ "Isochronic",
  TRUE ~ "Other"
)

# Comprobar la distribución de grupos
print("\nDistribución de grupos:")
print(table(data$group))


#  ANOVA 
perform_anova <- function(df, group1, group2, age_column) {
  # Filtrar datos para los dos grupos
  sub_data <- df %>%
    filter(group %in% c(group1, group2)) %>%
    select(group, !!sym(age_column))
  
  # Realizar ANOVA y prueba t
  anova_result <- aov(as.formula(paste(age_column, "~ group")), data = sub_data)
  t_test <- t.test(as.formula(paste(age_column, "~ group")), data = sub_data)
  
  return(list(
    anova = summary(anova_result),
    t_test = t_test,
    p_value = t_test$p.value
  ))
}

# Función para crear gráfico de barras con error bars
plot_bar_comparison <- function(df, group1, group2, age_column = "DNAm.Age", title) {
  # Filtrar datos
  sub_data <- df %>%
    filter(group %in% c(group1, group2)) %>%
    select(group, !!sym(age_column))
  
  # Calcular estadísticas para las barras de error
  stats_data <- sub_data %>%
    group_by(group) %>%
    summarise(
      mean = mean(!!sym(age_column), na.rm = TRUE),
      sd = sd(!!sym(age_column), na.rm = TRUE),
      n = n(),
      se = sd/sqrt(n)
    )
  
  # Definir colores personalizados
  custom_colors <- c(
    "Heterochronic" = "skyblue",
    "Heterochronic Recovery" = "skyblue",
    "Isochronic" = "indianred",
    "Isochronic Recovery" = "indianred"
  )
  
  # Crear el gráfico de barras
  p <- ggplot(stats_data, aes(x = group, y = mean, fill = group)) +
    geom_bar(stat = "identity", alpha = 0.7, width = 0.5) +
    geom_errorbar(aes(ymin = mean - se, ymax = mean + se), 
                  width = 0.2,
                  position = position_dodge(0.9)) +
    geom_jitter(data = sub_data, 
                aes(x = group, y = .data[[age_column]]),
                width = 0.1, 
                alpha = 0.5, 
                size = 2) +
    labs(title = title,
         x = "",
         y = "Predicted DNAm Age (HumanRatPanTissueLogLinearAge)",
         subtitle = paste("n =", nrow(sub_data))) +
    theme_minimal() +
    theme(legend.position = "none",
          plot.title = element_text(size = 14, face = "bold"),
          plot.subtitle = element_text(size = 10),
          axis.text = element_text(size = 12),
          axis.title = element_text(size = 12)) +
    scale_fill_manual(values = custom_colors)
  
  return(p)
}

# Función para crear boxplot
plot_boxplot <- function(df, group1, group2, age_column = "DNAm.Age", title) {
  # Filtrar datos
  sub_data <- df %>%
    filter(group %in% c(group1, group2)) %>%
    select(group, !!sym(age_column))
  
  # Realizar análisis estadístico
  stats <- perform_anova(sub_data, group1, group2, age_column)
  
  # Definir colores personalizados
  custom_colors <- c(
    "Heterochronic" = "skyblue",
    "Heterochronic Recovery" = "skyblue",
    "Isochronic" = "indianred",
    "Isochronic Recovery" = "indianred"
  )
  
  # Crear el gráfico
  p <- ggplot(sub_data, aes(x = group, y = .data[[age_column]], fill = group)) +
    geom_boxplot(alpha = 0.7, outlier.shape = NA, width = 0.5) +
    stat_boxplot(geom = "errorbar", width = 0.2) +  # Añadir barras horizontales en los bigotes
    labs(title = title, 
         x = "", 
         y = "Predicted DNAm Age (HumanRatPanTissueLogLinearAge)",
         subtitle = paste("n =", nrow(sub_data), "\np-value =", format.pval(stats$p_value, digits = 3))) +
    theme_minimal() +
    theme(legend.position = "none",
          plot.title = element_text(size = 14, face = "bold"),
          plot.subtitle = element_text(size = 10),
          axis.text = element_text(size = 12),
          axis.title = element_text(size = 12)) +
    scale_fill_manual(values = custom_colors)
  
  # Añadir anotación de significancia estadística
  if (stats$p_value < 0.05) {
    p <- p + 
      stat_compare_means(method = "t.test", 
                        comparisons = list(c(group1, group2)),
                        label = "p.signif",
                        size = 5)
  }
  
  return(list(
    plot = p,
    stats = stats,
    is_significant = stats$p_value < 0.05
  ))
}

# Definir las comparaciones a realizar, en este caso será las recuperaciones y 
# edades en el tratamiento de cada grupo (isocrónico vs heterocrónico)
comparisons <- list(
  list(
    group1 = "Isochronic Recovery",
    group2 = "Heterochronic Recovery",
    age_column = "DNAm.FinalHumanRatPanTissueLogLinearAge",
    title = "Heterochronic Recovery vs Isochronic Recovery",
    filename = "recovery_comparison"
  ),
  list(
    group1 = "Heterochronic",
    group2 = "Isochronic",
    age_column = "DNAm.FinalHumanRatPanTissueLogLinearAge",
    title = "Heterochronic vs Isochronic",
    filename = "heterochronic_vs_isochronic"
  )
)

# Realizar las comparaciones y crear gráficos
for (comp in comparisons) {
  # Crear y guardar boxplot
  box_result <- plot_boxplot(data, comp$group1, comp$group2, comp$age_column, comp$title)
  ggsave(paste0(comp$filename, ".png"), box_result$plot, width = 8, height = 6)
  
  # Crear y guardar gráfico de barras
  bar_plot <- plot_bar_comparison(data, comp$group1, comp$group2, comp$age_column, comp$title)
  ggsave(paste0(comp$filename, "_bar.png"), bar_plot, width = 8, height = 6)
  
  # Imprimir resultados estadísticos
  print(paste("\nEstadísticas para", comp$title, ":"))
  print(box_result$stats$anova)
  print(paste("p-value =", format.pval(box_result$stats$p_value, digits = 3)))
}
