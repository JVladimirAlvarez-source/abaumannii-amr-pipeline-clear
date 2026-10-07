#' ==============================================================================
#' @title Pipeline de Extracción y Control de Calidad de Fenotipos AMR (BV-BRC)
#' @description Pipeline bioinformático modular para el procesamiento de genomas y 
#'              fenotipos de Acinetobacter baumannii desde la API de BV-BRC.
#' @author Jhonny Vladimir Alvarez Poma / IIFB
#' @date 2026-10-30
#' ==============================================================================

# ------------------------------------------------------------------------------
# 1. GESTIÓN DE DEPENDENCIAS
# ------------------------------------------------------------------------------
if (!require("pacman")) install.packages("pacman")                    
pacman::p_load(
  httr,       # Solicitudes HTTP robustas
  jsonlite,   # Parsing y estructuración de JSON
  dplyr,      # Manipulación y transformación de datos (tidyverse)
  purrr,      # Programación funcional y manejo de listas
  tidyr,      # Remodelado de datos (pivot_wider)
  readr       # Exportación eficiente de archivos planos
)

# ------------------------------------------------------------------------------
# 2. CONFIGURACIÓN CENTRALIZADA
# ------------------------------------------------------------------------------
CONFIG <- list(
  TAXON_ID        = 470,                  # Taxón objetivo (470 = Acinetobacter baumannii)
  CHUNK_SIZE      = 100,                  # Tamaño de lote optimizado para consultas API
  MIN_SUSCEPTIBLE = 100,                  # Umbral analítico mínimo: aislados susceptibles
  MIN_RESISTANT   = 100,                  # Umbral analítico mínimo: aislados resistentes
  QUALITY_FILTER  = c("good"),            # Criterio de calidad de ensamblado
  STATUS_FILTER   = c("complete", "wgs"), # Criterio de estado del genoma
  DIR_OUTPUT      = "resultados_analiticos/" # Directorio centralizado de salida
)

# Crear directorio de salida si no existe de forma recursiva
if (!dir.exists(CONFIG$DIR_OUTPUT)) {
  dir.create(CONFIG$DIR_OUTPUT, recursive = TRUE)
}

# Endpoints oficiales de la API de BV-BRC
URL_GENOME <- "https://www.bv-brc.org/api/genome/"
URL_AMR    <- "https://www.bv-brc.org/api/genome_amr/"

# ------------------------------------------------------------------------------
# 3. FUNCIONES AUXILIARES (CAPA DE CONEXIÓN)
# ------------------------------------------------------------------------------

#' Consulta modular a la API de BV-BRC mediante lotes (chunks) de IDs
#' @param endpoint URL del servicio API
#' @param ids_vector Vector con los identificadores (genome_id) a consultar
#' @param select_fields Campos específicos a recuperar (opcional)
#' @param extra_rql Filtros RQL adicionales (opcional)
#' @return Data frame consolidado con los resultados
fetch_bvbrc_in_chunks <- function(endpoint, ids_vector, select_fields = "", extra_rql = "") {
  
  lotes <- split(ids_vector, ceiling(seq_along(ids_vector) / CONFIG$CHUNK_SIZE))
  lista_resultados <- list()
  
  cat(sprintf("🔎 Consultando %d elementos en %d lotes...\n", length(ids_vector), length(lotes)))
  
  for (i in seq_along(lotes)) {
    cadena_ids <- paste(lotes[[i]], collapse = ",")
    
    # Construcción de la sentencia RQL
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
      warning(sprintf("⚠️ Alerta en lote %d. Código HTTP: %d", i, status_code(res)))
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
# 4. EJECUCIÓN DEL PIPELINE ANALÍTICO
# ==============================================================================

cat("=========================================================\n")
cat("🚀 PASO 1: Extracción del Universo Inicial de Genomas\n")
cat("=========================================================\n")

rql_genomas <- paste0("eq(taxon_lineage_ids,", CONFIG$TAXON_ID, ")&select(genome_id,genome_name,taxon_id)&limit(25000)")

res_genomas <- GET(
  url = paste0(URL_GENOME, "?", rql_genomas),
  add_headers("Accept" = "application/json")
)

if (status_code(res_genomas) == 200) {
  df_genomas <- fromJSON(content(res_genomas, as = "text", encoding = "UTF-8"))
  todos_los_genome_ids <- unique(df_genomas$genome_id)
  cat("✅ Total de genomas detectados en taxón objetivo:", length(todos_los_genome_ids), "\n\n")
} else {
  stop("❌ Error crítico al recuperar genomas. Código HTTP: ", status_code(res_genomas))
}


cat("=========================================================\n")
cat("📡 PASO 2: Extracción y Limpieza de Fenotipos de Laboratorio\n")
cat("=========================================================\n")

df_amr_todos <- fetch_bvbrc_in_chunks(
  endpoint = URL_AMR,
  ids_vector = todos_los_genome_ids,
  extra_rql = "ne(evidence,Computational%20Prediction)"
)

if (nrow(df_amr_todos) > 0) {
  df_fenotipos_reales <- df_amr_todos %>%
    filter(!grepl("Computational", evidence, ignore.case = TRUE))
  
  df_fenotipos_clean <- df_fenotipos_reales %>%
    mutate(across(where(is.list), ~ map_chr(.x, ~ paste(.x[!is.na(.x)], collapse = "; "))))
  
  # Persistencia de datos crudos iniciales
  write.csv(df_fenotipos_clean, file.path(CONFIG$DIR_OUTPUT, "abaumannii_phenotypes_raw.csv"), row.names = FALSE)
  
  # Filtrado estructural eliminando registros sin fenotipo o método de tipificación
  df_fenotipos_sin_na <- df_fenotipos_clean %>%
    filter(!is.na(resistant_phenotype) & resistant_phenotype != "" & resistant_phenotype != "NA") %>%
    filter(!is.na(laboratory_typing_method) & laboratory_typing_method != "" & laboratory_typing_method != "NA")
  
  write.csv(df_fenotipos_sin_na, file.path(CONFIG$DIR_OUTPUT, "abaumannii_phenotypes_clean.csv"), row.names = FALSE)
} else {
  stop("⚠️ No se encontraron registros de fenotipos experimentales.")
}


cat("\n=========================================================\n")
cat("🧬 PASO 3: Validación de Calidad Bioinformática y Ensamblado\n")
cat("=========================================================\n")

ids_unicos_fenotipo <- unique(df_fenotipos_sin_na$genome_id)

df_metadatos_genomas <- fetch_bvbrc_in_chunks(
  endpoint = URL_GENOME,
  ids_vector = ids_unicos_fenotipo,
  select_fields = "genome_id,genome_name,genome_quality,genome_status"
)

df_genomas_alta_calidad <- df_metadatos_genomas %>%
  filter(tolower(genome_quality) %in% CONFIG$QUALITY_FILTER) %>%
  filter(tolower(genome_status) %in% CONFIG$STATUS_FILTER)

ids_alta_calidad <- unique(df_genomas_alta_calidad$genome_id)

# Persistencia de la lista base de genomas de calidad aprobada
write.csv(df_genomas_alta_calidad, file.path(CONFIG$DIR_OUTPUT, "abaumannii_genomes_quality_pass.csv"), row.names = FALSE)

df_fenotipos_filtrados_calidad <- df_fenotipos_sin_na %>%
  filter(genome_id %in% ids_alta_calidad)

cat("• Genomas validados con calidad 'Good' y estado 'Complete/WGS':", length(ids_alta_calidad), "\n")


cat("\n=========================================================\n")
cat("🌍 PASO 4: Pre-filtro Geográfico por Criterios de Exclusión\n")
cat("=========================================================\n")

df_geo_genomas <- fetch_bvbrc_in_chunks(
  endpoint = URL_GENOME,
  ids_vector = ids_alta_calidad,
  select_fields = "genome_id,geographic_group,isolation_country"
)

if (!"geographic_group" %in% colnames(df_geo_genomas)) df_geo_genomas$geographic_group <- NA
if (!"isolation_country" %in% colnames(df_geo_genomas)) df_geo_genomas$isolation_country <- NA

df_geo_evaluado <- df_geo_genomas %>%
  mutate(
    Continente = ifelse(is.na(geographic_group) | geographic_group == "", "No especificado", geographic_group),
    Pais = ifelse(is.na(isolation_country) | isolation_country == "", "No especificado", isolation_country),
    # Criterio de exclusión: Perú, Honduras O (No especificado en ambas dimensiones simultáneamente)
    excluido_geo = (Pais %in% c("Peru", "Honduras", "Perú")) | (Continente == "No especificado" & Pais == "No especificado")
  )

ids_geo_excluidos <- df_geo_evaluado %>% filter(excluido_geo) %>% pull(genome_id)
ids_geo_aprobados <- df_geo_evaluado %>% filter(!excluido_geo) %>% pull(genome_id)

# Exportar archivo independiente de registros excluidos por geografía
df_fenotipos_geo_excluidos <- df_fenotipos_filtrados_calidad %>%
  filter(genome_id %in% ids_geo_excluidos)
write.csv(df_fenotipos_geo_excluidos, file.path(CONFIG$DIR_OUTPUT, "audit_excluded_geographic.csv"), row.names = FALSE)

df_fenotipos_post_geo <- df_fenotipos_filtrados_calidad %>%
  filter(genome_id %in% ids_geo_aprobados)

cat("• Genomas excluidos por criterio geográfico:", length(ids_geo_excluidos), "\n")
cat("• Genomas aprobados tras filtro geográfico:", length(ids_geo_aprobados), "\n")


cat("\n=========================================================\n")
cat("💊 PASO 5: Filtrado de Antibióticos Robustos (>= 100 S y >= 100 R)\n")
cat("=========================================================\n")

df_conteo_antibioticos <- df_fenotipos_post_geo %>%
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

if (!"Susceptible" %in% colnames(df_conteo_antibioticos)) df_conteo_antibioticos$Susceptible <- 0
if (!"Resistente" %in% colnames(df_conteo_antibioticos)) df_conteo_antibioticos$Resistente <- 0

df_antibioticos_robustos <- df_conteo_antibioticos %>%
  filter(Susceptible >= CONFIG$MIN_SUSCEPTIBLE & Resistente >= CONFIG$MIN_RESISTANT) %>%
  arrange(desc(Susceptible + Resistente))

write.csv(df_antibioticos_robustos, file.path(CONFIG$DIR_OUTPUT, "antibiotics_summary_robust.csv"), row.names = FALSE)

antibioticos_validos <- df_antibioticos_robustos$antibiotic
df_fenotipos_analiticos <- df_fenotipos_post_geo %>%
  filter(antibiotic %in% antibioticos_validos)

ids_genomas_analiticos <- unique(df_fenotipos_analiticos$genome_id)
ids_abx_excluidos <- setdiff(ids_geo_aprobados, ids_genomas_analiticos)

cat("• Antibióticos que superaron el umbral de robustez:", nrow(df_antibioticos_robustos), "\n")
cat("• Genomas finales retenidos para el set analítico:", length(ids_genomas_analiticos), "\n")


cat("\n=========================================================\n")
cat("📋 PASO 6: Generación de Auditoría y Datasets Analíticos Finales\n")
cat("=========================================================\n")

# 1. Auditoría Integral del Ciclo de Vida de los Genomas
df_auditoria_ciclo_vida <- data.frame(
  genome_id = c(ids_genomas_analiticos, ids_geo_excluidos, ids_abx_excluidos),
  motivo_exclusion = c(
    rep("APROBADO (Cumplió filtros de calidad, geografía y umbrales analíticos)", length(ids_genomas_analiticos)),
    rep("EXCLUIDO_GEOGRAFIA (Perú, Honduras o 'No especificado' dual en región/país)", length(ids_geo_excluidos)),
    rep("EXCLUIDO_ANTIBIOTICOS (Insuficientes ensayos en antibióticos robustos)", length(ids_abx_excluidos))
  )
)
write.csv(df_auditoria_ciclo_vida, file.path(CONFIG$DIR_OUTPUT, "audit_genomes_lifecycle.csv"), row.names = FALSE)

# 2. Datasets Analíticos Principales (Estándar Profesional sin sufijos hardcodeados)
write.csv(df_fenotipos_analiticos, file.path(CONFIG$DIR_OUTPUT, "abaumannii_phenotypes_analytical.csv"), row.names = FALSE)

df_genomas_analiticos <- df_genomas_alta_calidad %>%
  filter(genome_id %in% ids_genomas_analiticos)
write.csv(df_genomas_analiticos, file.path(CONFIG$DIR_OUTPUT, "abaumannii_genomes_analytical.csv"), row.names = FALSE)

# 3. Resumen y Distribución Geográfica Limpia Final
df_geografia_analitica <- df_geo_evaluado %>%
  filter(genome_id %in% ids_genomas_analiticos) %>%
  mutate(
    Continente = ifelse(is.na(geographic_group) | geographic_group == "", "No especificado", geographic_group),
    Pais = ifelse(is.na(isolation_country) | isolation_country == "", "No especificado", isolation_country)
  )

tabla_distribucion_geo <- df_geografia_analitica %>%
  group_by(Continente, Pais) %>%
  summarise(Numero_de_Genomas = n_distinct(genome_id), .groups = "drop") %>%
  arrange(Continente, desc(Numero_de_Genomas))

write.csv(tabla_distribucion_geo, file.path(CONFIG$DIR_OUTPUT, "geographic_distribution_summary.csv"), row.names = FALSE)

cat("🎯 ¡PIPELINE EJECUTADO Y ARQUITECTURA DE DATOS APLICADA EXITOSAMENTE!\n")
cat("• Todos los reportes estandarizados se encuentran listos en la ruta:", CONFIG$DIR_OUTPUT, "\n")
cat("=========================================================\n")