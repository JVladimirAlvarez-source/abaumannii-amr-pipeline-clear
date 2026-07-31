#' ==============================================================================
#' @title Pipeline de Extracción y Control de Calidad de Fenotipos AMR (BV-BRC)
#' @description Descarga, filtra y audita genomas y fenotipos de laboratorio.
#' @author Tu Nombre / Grupo de Investigación
#' @date 2026-07-30
#' ==============================================================================

# ------------------------------------------------------------------------------
# 1. GESTIÓN DE DEPENDENCIAS
# ------------------------------------------------------------------------------
if (!require("pacman")) install.packages("pacman")                   
pacman::p_load(
  httr,       # Solicitudes HTTP
  jsonlite,   # Parsing de JSON
  dplyr,      # Manipulación de datos
  purrr,      # Programación funcional
  tidyr,      # Remodelado de datos (pivot_wider)
  readr       # Exportación eficiente de CSV
)

# ------------------------------------------------------------------------------
# 2. CONFIGURACIÓN CENTRALIZADA
# ------------------------------------------------------------------------------
CONFIG <- list(
  TAXON_ID        = 470,                  # Taxón objetivo (470 = A. baumannii)
  CHUNK_SIZE      = 100,                  # Tamaño de lote para la API
  MIN_SUSCEPTIBLE = 100,                  # Umbral mínimo de aislados susceptibles
  MIN_RESISTANT   = 100,                  # Umbral mínimo de aislados resistentes
  QUALITY_FILTER  = c("good"),            # Calidad de ensamblado
  STATUS_FILTER   = c("complete", "wgs"), # Estado del genoma
  DIR_OUTPUT      = "resultados/"         # Carpeta de salida
)

# Crear directorio de salida si no existe
if (!dir.exists(CONFIG$DIR_OUTPUT)) dir.create(CONFIG$DIR_OUTPUT, recursive = TRUE)

# Endpoints de la API
URL_GENOME <- "https://www.bv-brc.org/api/genome/"
URL_AMR    <- "https://www.bv-brc.org/api/genome_amr/"

# ------------------------------------------------------------------------------
# 3. FUNCIONES AUXILIARES
# ------------------------------------------------------------------------------

#' Consultar la API de BV-BRC por lotes de IDs
fetch_bvbrc_in_chunks <- function(endpoint, ids_vector, select_fields = "", extra_rql = "") {
  
  lotes <- split(ids_vector, ceiling(seq_along(ids_vector) / CONFIG$CHUNK_SIZE))
  lista_resultados <- list()
  
  cat(sprintf("🔎 Consultando %d elementos en %d lotes...\n", length(ids_vector), length(lotes)))
  
  for (i in seq_along(lotes)) {
    cadena_ids <- paste(lotes[[i]], collapse = ",")
    
    # Construcción RQL
    rql <- paste0("in(genome_id,(", cadena_ids, "))")
    if (nchar(select_fields) > 0) rql <- paste0(rql, "&select(", select_fields, ")")
    if (nchar(extra_rql) > 0)     rql <- paste0(rql, "&", extra_rql)
    rql <- paste0(rql, "&limit(10000)")
    
    res <- GET(url = paste0(endpoint, "?", rql), add_headers("Accept" = "application/json"))
    
    if (status_code(res) == 200) {
      texto <- content(res, as = "text", encoding = "UTF-8")
      df_chunk <- fromJSON(texto)
      if (is.data.frame(df_chunk) && nrow(df_chunk) > 0) {
        lista_resultados[[i]] <- df_chunk
      }
    } else {
      warning(sprintf("⚠️ Fallo en lote %d. Código HTTP: %d", i, status_code(res)))
    }
    
    if (i %% 20 == 0 || i == length(lotes)) {
      cat(sprintf("⏳ Lote %d de %d procesado (%.1f%%)\n", i, length(lotes), (i / length(lotes)) * 100))
    }
  }
  
  if (length(lista_resultados) > 0) {
    return(bind_rows(lista_resultados))
  } else {
    return(data.frame())
  }
}

# ==============================================================================
# 4. EJECUCIÓN DEL PIPELINE
# ==============================================================================

cat("=========================================================\n")
cat("🚀 PASO 1: Extrayendo TODOS los genome_ids\n")
cat("=========================================================\n")

rql_genomas <- paste0("eq(taxon_lineage_ids,", CONFIG$TAXON_ID, ")&select(genome_id,genome_name,taxon_id)&limit(25000)")

res_genomas <- GET(
  url = paste0(URL_GENOME, "?", rql_genomas),
  add_headers("Accept" = "application/json")
)

if (status_code(res_genomas) == 200) {
  df_genomas <- fromJSON(content(res_genomas, as = "text", encoding = "UTF-8"))
  todos_los_genome_ids <- unique(df_genomas$genome_id)
  
  cat("✅ Total de genomas/cepas encontrados:", length(todos_los_genome_ids), "\n")
  cat("✅ Total de Taxon IDs distintos encontrados:", length(unique(df_genomas$taxon_id)), "\n\n")
} else {
  stop("❌ Error al recuperar la lista de genomas. Código HTTP: ", status_code(res_genomas))
}

cat("=========================================================\n")
cat("📡 PASO 2: Extrayendo fenotipos de laboratorio\n")
cat("=========================================================\n")

df_amr_todos <- fetch_bvbrc_in_chunks(
  endpoint = URL_AMR,
  ids_vector = todos_los_genome_ids,
  extra_rql = "ne(evidence,Computational%20Prediction)"
)

cat("\n=========================================================\n")
cat("📊 PASO 3: Consolidando y guardando datos\n")
cat("=========================================================\n")

if (nrow(df_amr_todos) > 0) {
  df_fenotipos_reales <- df_amr_todos %>%
    filter(!grepl("Computational", evidence, ignore.case = TRUE))
  
  df_fenotipos_clean <- df_fenotipos_reales %>%
    mutate(across(where(is.list), ~ map_chr(.x, ~ paste(.x[!is.na(.x)], collapse = "; "))))
  
  archivo_salida <- file.path(CONFIG$DIR_OUTPUT, "A_baumannii_TODOS_LOS_FENOTIPOS_LABORATORIO.csv")
  write.csv(df_fenotipos_clean, archivo_salida, row.names = FALSE)
  
  cat("🎯 ¡ÉXITO EN EXTRACCIÓN!\n")
  cat("---------------------------------------------------------\n")
  cat("• Total de registros de laboratorio (filas):", nrow(df_fenotipos_clean), "\n")
  cat("• Cepas/genomas únicos con fenotipo real:", length(unique(df_fenotipos_clean$genome_id)), "\n")
  cat("• Archivo guardado en:", archivo_salida, "\n")
  cat("---------------------------------------------------------\n\n")
  
  cat("📋 Top 10 antibióticos con más ensayos de laboratorio:\n")
  print(head(sort(table(df_fenotipos_clean$antibiotic), decreasing = TRUE), 10))
  
  # Filtrar NAs y valores vacíos
  df_filtrado_sin_na <- df_fenotipos_clean %>%
    filter(!is.na(resistant_phenotype) & resistant_phenotype != "" & resistant_phenotype != "NA") %>%
    filter(!is.na(laboratory_typing_method) & laboratory_typing_method != "" & laboratory_typing_method != "NA")
  
  archivo_salida_sin_na <- file.path(CONFIG$DIR_OUTPUT, "A_baumannii_FENOTIPOS_LABORATORIO_SIN_NA.csv")
  write.csv(df_filtrado_sin_na, archivo_salida_sin_na, row.names = FALSE)
  
  cat("\n=========================================================\n")
  cat("📊 REPORTE DE DATOS TRAS ELIMINAR NAs\n")
  cat("=========================================================\n")
  cat("• Total de registros de laboratorio (filas):", nrow(df_filtrado_sin_na), "\n")
  cat("• Cepas/genomas únicos con fenotipo real:", length(unique(df_filtrado_sin_na$genome_id)), "\n")
  cat("• Archivo guardado en:", archivo_salida_sin_na, "\n")
  cat("=========================================================\n")
} else {
  stop("⚠️ No se encontraron fenotipos de laboratorio.")
}

cat("\n=========================================================\n")
cat("🧬 PASO 4: Verificación de Calidad Bioinformática y Estado\n")
cat("=========================================================\n")

ids_unicos_filtrados <- unique(df_filtrado_sin_na$genome_id)

df_metadatos_genomas <- fetch_bvbrc_in_chunks(
  endpoint = URL_GENOME,
  ids_vector = ids_unicos_filtrados,
  select_fields = "genome_id,genome_name,genome_quality,genome_status"
)

if (nrow(df_metadatos_genomas) > 0) {
  df_genomas_alta_calidad <- df_metadatos_genomas %>%
    filter(tolower(genome_quality) %in% CONFIG$QUALITY_FILTER) %>%
    filter(tolower(genome_status) %in% CONFIG$STATUS_FILTER)
  
  ids_alta_calidad <- unique(df_genomas_alta_calidad$genome_id)
  
  df_fenotipos_final_filtrado <- df_filtrado_sin_na %>%
    filter(genome_id %in% ids_alta_calidad)
  
  archivo_genomas_calidad <- file.path(CONFIG$DIR_OUTPUT, "lista_genomas_calidad_GOOD_WGS_COMPLETE.csv")
  archivo_fenotipos_final <- file.path(CONFIG$DIR_OUTPUT, "A_baumannii_FENOTIPOS_FINAL_ALTA_CALIDAD.csv")
  
  write.csv(df_genomas_alta_calidad, archivo_genomas_calidad, row.names = FALSE)
  write.csv(df_fenotipos_final_filtrado, archivo_fenotipos_final, row.names = FALSE)
  
  cat("\n=========================================================\n")
  cat("🎯 REPORTE FINAL DE CONTROL DE CALIDAD\n")
  cat("=========================================================\n")
  cat("• Genomas iniciales con fenotipo sin NA:", length(ids_unicos_filtrados), "\n")
  cat("• Genomas validados con Calidad 'Good' y Ensamblado 'Complete/WGS':", length(ids_alta_calidad), "\n")
  cat("• Filas totales de ensayos de laboratorio retenidos:", nrow(df_fenotipos_final_filtrado), "\n")
  cat("• Lista de genomas filtrados guardada en:", archivo_genomas_calidad, "\n")
  cat("• Dataset fenotípico final guardado en:", archivo_fenotipos_final, "\n")
  cat("=========================================================\n\n")
  
  cat("📊 Desglose de estados de ensamblado en los genomas de alta calidad:\n")
  print(table(df_genomas_alta_calidad$genome_status))
} else {
  stop("❌ No se pudieron recuperar los metadatos de calidad de la API.")
}

cat("\n=========================================================\n")
cat("💊 PASO 5: Filtrado de Antibióticos (>=100 S y >=100 R)\n")
cat("=========================================================\n")

df_conteo_abx <- df_fenotipos_final_filtrado %>%
  mutate(fenotipo_norm = tolower(trimws(resistant_phenotype))) %>%
  filter(fenotipo_norm %in% c("susceptible", "resistant", "susceptible/intermediate", "resistant/intermediate")) %>%
  mutate(categoria = ifelse(grepl("susceptible", fenotipo_norm), "Susceptible", "Resistente")) %>%
  group_by(antibiotic, categoria) %>%
  summarise(n_aislados = n_distinct(genome_id), .groups = "drop") %>%
  tidyr::pivot_wider(
    names_from = categoria, 
    values_from = n_aislados, 
    values_fill = 0
  )

if (!"Susceptible" %in% colnames(df_conteo_abx)) df_conteo_abx$Susceptible <- 0
if (!"Resistente" %in% colnames(df_conteo_abx)) df_conteo_abx$Resistente <- 0

df_abx_retenidos <- df_conteo_abx %>%
  filter(Susceptible >= CONFIG$MIN_SUSCEPTIBLE & Resistente >= CONFIG$MIN_RESISTANT) %>%
  arrange(desc(Susceptible + Resistente))

write.csv(df_abx_retenidos, file.path(CONFIG$DIR_OUTPUT, "tabla_antibioticos_filtrados_ge100.csv"), row.names = FALSE)

cat("📋 TABLA DE ANTIBIÓTICOS RETENIDOS (>= 100 S y >= 100 R):\n")
cat("---------------------------------------------------------\n")
print(as.data.frame(df_abx_retenidos))
cat("---------------------------------------------------------\n")
cat("• Total de antibióticos que superaron el umbral:", nrow(df_abx_retenidos), "\n\n")

antibioticos_validos <- df_abx_retenidos$antibiotic

df_fenotipos_abx_filtrados <- df_fenotipos_final_filtrado %>%
  filter(antibiotic %in% antibioticos_validos)

ids_genomas_finales <- unique(df_fenotipos_abx_filtrados$genome_id)

cat("=========================================================\n")
cat("🌍 PASO 6: Extracción de Geografía (Continente y País)\n")
cat("=========================================================\n")

df_geo_genomas <- fetch_bvbrc_in_chunks(
  endpoint = URL_GENOME,
  ids_vector = ids_genomas_finales,
  select_fields = "genome_id,geographic_group,isolation_country"
)

if (nrow(df_geo_genomas) > 0) {
  if (!"geographic_group" %in% colnames(df_geo_genomas)) df_geo_genomas$geographic_group <- NA
  if (!"isolation_country" %in% colnames(df_geo_genomas)) df_geo_genomas$isolation_country <- NA
  
  df_geo_clean <- df_geo_genomas %>%
    mutate(
      Continente = ifelse(is.na(geographic_group) | geographic_group == "", "No especificado", geographic_group),
      Pais = ifelse(is.na(isolation_country) | isolation_country == "", "No especificado", isolation_country)
    )
  
  tabla_geografica <- df_geo_clean %>%
    group_by(Continente, Pais) %>%
    summarise(Numero_de_Genomas = n_distinct(genome_id), .groups = "drop") %>%
    arrange(Continente, desc(Numero_de_Genomas))
  
  archivo_geo <- file.path(CONFIG$DIR_OUTPUT, "distribucion_genomas_continente_pais.csv")
  write.csv(tabla_geografica, archivo_geo, row.names = FALSE)
  
  cat("\n=========================================================\n")
  cat("🎯 RESUMEN GEOGRÁFICO DE GENOMAS RETENIDOS\n")
  cat("=========================================================\n")
  print(as.data.frame(tabla_geografica))
  cat("---------------------------------------------------------\n")
  cat("• Archivo guardado en:", archivo_geo, "\n")
  cat("=========================================================\n")
} else {
  cat("⚠️ No se pudieron extraer los metadatos geográficos de los genomas.\n")
}

# ------------------------------------------------------------------------------
# AUDITORÍA OPCIONAL: Si la variable ids_223_descartados existe
# ------------------------------------------------------------------------------
if (exists("ids_223_descartados")) {
  cat("\n=========================================================\n")
  cat("📊 AUDITORÍA: Fenotipos de los 223 Genomas Excluidos\n")
  cat("=========================================================\n")
  
  df_223_fenotipos <- df_fenotipos_final_filtrado %>%
    filter(genome_id %in% ids_223_descartados) %>%
    mutate(fenotipo_norm = tolower(trimws(resistant_phenotype))) %>%
    filter(fenotipo_norm %in% c("susceptible", "resistant", "susceptible/intermediate", "resistant/intermediate")) %>%
    mutate(Categoria = ifelse(grepl("susceptible", fenotipo_norm), "Susceptible", "Resistente"))
  
  resumen_global_223 <- df_223_fenotipos %>%
    group_by(Categoria) %>%
    summarise(Genomas_Unicos = n_distinct(genome_id), .groups = "drop")
  
  cat("📈 RESUMEN GLOBAL (N° de genomas únicos con al menos un ensayo S o R):\n")
  print(as.data.frame(resumen_global_223))
  cat("---------------------------------------------------------\n\n")
  
  tabla_abx_223_SR <- df_223_fenotipos %>%
    group_by(antibiotic, Categoria) %>%
    summarise(n_aislados = n_distinct(genome_id), .groups = "drop") %>%
    tidyr::pivot_wider(
      names_from = Categoria, 
      values_from = n_aislados, 
      values_fill = 0
    )
  
  if (!"Susceptible" %in% colnames(tabla_abx_223_SR)) tabla_abx_223_SR$Susceptible <- 0
  if (!"Resistente" %in% colnames(tabla_abx_223_SR)) tabla_abx_223_SR$Resistente <- 0
  
  tabla_abx_223_SR <- tabla_abx_223_SR %>%
    mutate(Total_Genomas = Susceptible + Resistente) %>%
    arrange(desc(Total_Genomas))
  
  archivo_audit <- file.path(CONFIG$DIR_OUTPUT, "audit_223_genomas_fenotipos_SR.csv")
  write.csv(tabla_abx_223_SR, archivo_audit, row.names = FALSE)
  
  cat("💊 DESGLOSE POR ANTIBIÓTICO EXCLUIDO (Susceptibles vs Resistentes):\n")
  print(as.data.frame(tabla_abx_223_SR))
  
  cat("\n=========================================================\n")
  cat("💾 Archivo guardado en:", archivo_audit, "\n")
  cat("=========================================================\n")
}

