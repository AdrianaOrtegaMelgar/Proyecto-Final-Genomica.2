# Universidad del Valle de Guatemala
# Genómica Aplicada
# Integrantes:
# Daniella Navarro - 231179
# Alexa Castro - 231194
# María Fernanda Silva - 221035
# Adriana Ortega - 231001

# IMPACTO DE LA HETEROGENEIDAD MOLECULAR EN LA IDENTIFICACIÓN DE GENES
# DIFERENCIALMENTE EXPRESADOS EN CÁNCER DE MAMA:
# ANÁLISIS GLOBAL VERSUS ESTRATIFICADO POR SUBTIPO (TCGA-BRCA)

# DESCRIPCIÓN
# -----------------------------------------------------------------------------
# En este proyecto se compara dos estrategias de análisis de expresión diferencial en
# cáncer de mama utilizando datos públicos de RNA-seq del proyecto TCGA-BRCA, por medio de 2 enfoques:

#   Enfoque 1 - global:       todas las muestras tumorales vs. tejido normal
#   Enfoque 2 - estratificado: cada subtipo molecular vs. tejido normal 
#                              (Luminal A, Luminal B, HER2-enriquecido, Basal)

# El objetivo es cuantificar en qué medida el agrupamiento de subtipos
# biológicamente distintos modifica el conjunto de genes identificados como
# diferencialmente expresados
# -----------------------------------------------------------------------------


# ESTRUCTURA DEL SCRIPT
# -----------------------------------------------------------------------------
#   Instalación y carga de paquetes
#   1. Verificación de subtipos moleculares
#   2. Consulta y descarga de datos (GDC)
#   3. Control de calidad y filtrado
#   4. Análisis de expresión diferencial global
#   5. Análisis de expresión diferencial estratificado
#   6. Comparación entre enfoques
#   7. Enriquecimiento funcional (GO / KEGG) y comparación de rutas
#   8. Contraste con biomarcadores reportados
# -----------------------------------------------------------------------------


# FUENTE DE DATOS
# -----------------------------------------------------------------------------
# GDC Data Portal — https://portal.gdc.cancer.gov
# Proyecto: TCGA-BRCA | Workflow: STAR - Counts
# -----------------------------------------------------------------------------


# INSTALACIÓN Y CARGA DE PAQUETES
# -----------------------------------------------------------------------------
# IMPORTANTE:
# 1) La instalación se hace una sola vez, en una sesión limpia de R.
# 2) Si aparece un error indicando que S4Vectors no puede descargarse,
#    reinicia RStudio (Session > Restart R) y ejecuta la instalación antes
#    de cargar cualquier paquete con library().
# 3) No ejecutes BiocManager::install() en cada corrida del análisis.

# INSTALACIÓN (descomentar y ejecutar por separado solo si faltan paquetes):
# if (!requireNamespace("BiocManager", quietly = TRUE))
#   install.packages("BiocManager")
# BiocManager::install(c(
#   "TCGAbiolinks", "DESeq2", "SummarizedExperiment",
#   "clusterProfiler", "org.Hs.eg.db", "AnnotationDbi", "ashr"
# ), ask = FALSE, update = FALSE)
# install.packages(c("dplyr", "ggplot2", "pheatmap", "UpSetR", "ggrepel"))

# CARGA DE LIBRERÍAS
suppressPackageStartupMessages({
  library(TCGAbiolinks)
  library(SummarizedExperiment)
  library(DESeq2)
  library(clusterProfiler)
  library(org.Hs.eg.db)
  library(AnnotationDbi)
  library(dplyr)
  library(ggplot2)
})

# Semilla para reproducibilidad
set.seed(123)

# Carpetas de salida
dir.create("resultados", showWarnings = FALSE)
dir.create("figuras",    showWarnings = FALSE)
# -----------------------------------------------------------------------------


# 1. VERIFICACIÓN DE SUBTIPOS MOLECULARES
# -----------------------------------------------------------------------------
# Este paso se ejecuta ANTES de la descarga de datos de expresión para
# confirmar la disponibilidad y completitud de la anotación de subtipos.

# Descargar tabla de subtipos

subtipos <- TCGAquery_subtype(tumor = "BRCA")

# Inspección general

cat("Dimensiones de la tabla de subtipos:", dim(subtipos), "\n")
cat("Número de muestras anotadas:", nrow(subtipos), "\n\n")

# Distribución por subtipo PAM50

cat("Distribución de subtipos moleculares (PAM50):\n")
print(table(subtipos$BRCA_Subtype_PAM50, useNA = "ifany"))

# Cuantificar datos faltantes

n_na <- sum(is.na(subtipos$BRCA_Subtype_PAM50))
cat("\nMuestras sin clasificación PAM50:", n_na,
    sprintf("(%.1f%%)", 100 * n_na / nrow(subtipos)), "\n")

# Guardar para referencia

write.csv(subtipos, "resultados/01_tabla_subtipos.csv", row.names = FALSE)
# -----------------------------------------------------------------------------


# 2. CONSULTA Y DESCARGA DE DATOS
# -----------------------------------------------------------------------------

# 2.1 Construir la consulta al GDC (no descarga aún) 

query <- GDCquery(
  project       = "TCGA-BRCA",
  data.category = "Transcriptome Profiling",
  data.type     = "Gene Expression Quantification",
  workflow.type = "STAR - Counts"
)

# 2.2 Revisar qué se descargaría

resultados_query <- getResults(query)

cat("Total de archivos disponibles:", nrow(resultados_query), "\n\n")
cat("Distribución por tipo de muestra:\n")
print(table(resultados_query$sample_type))


# 2.3 Descarga 
# los archivos quedan guardados en la carpeta GDCdata/

# Comprobar que la consulta sigue disponible
nrow(getResults(query))

# Revisar cuánto espacio hay en el disco
system("df -h ~")

# descargar GDCdata
GDCdownload(
  query,
  method = "api",
  directory = "GDCdata",
  files.per.chunk = 5
)

# Verificar el objeto de consulta
getResults(query) |> dim()

# Revisar la carpeta de trabajo
getwd()

# Verificar la carpeta de descarga
dir.exists("GDCdata")

# Contar los archivos descargados
length(list.files(
  "GDCdata",
  recursive = TRUE,
  full.names = TRUE
))

# 2.4 Construir el objeto SummarizedExperiment

datos_brca <- GDCprepare(query)

# 2.5 Guardar el objeto  (así no se repite la descarga)

saveRDS(datos_brca, "resultados/02_datos_brca_crudos.rds")

# Para volver a cargarlo sin volver a descargar:
# datos_brca <- readRDS("resultados/02_datos_brca_crudos.rds")

# 2.6 Inspección inicial

cat("Dimensiones del objeto:", dim(datos_brca), "\n")
cat("  Genes:  ", nrow(datos_brca), "\n")
cat("  Muestras:", ncol(datos_brca), "\n\n")

cat("Tipos de muestra disponibles:\n")
print(table(datos_brca$sample_type))
# -----------------------------------------------------------------------------


# 3. PREPARACIÓN, CONTROL DE CALIDAD Y FILTRADO
# -----------------------------------------------------------------------------
conteos_raw <- assay(datos_brca, "unstranded")
metadata_raw <- as.data.frame(colData(datos_brca))

# Garantizar que las filas de metadata coincidan con las columnas de conteos
if (!"barcode" %in% names(metadata_raw)) {
  metadata_raw$barcode <- colnames(conteos_raw)
}
metadata_raw <- metadata_raw[match(colnames(conteos_raw), rownames(metadata_raw)), , drop = FALSE]
stopifnot(nrow(metadata_raw) == ncol(conteos_raw))

# Solo tumor primario y tejido normal sólido
validas <- metadata_raw$sample_type %in% c("Primary Tumor", "Solid Tissue Normal")
conteos_raw <- conteos_raw[, validas, drop = FALSE]
metadata_raw <- metadata_raw[validas, , drop = FALSE]

# Cruce de subtipo por identificador de paciente TCGA (primeros 12 caracteres)
metadata_raw$patient_id <- substr(as.character(metadata_raw$barcode), 1, 12)
columna_paciente <- if ("patient" %in% names(subtipos)) "patient" else
  if ("submitter_id" %in% names(subtipos)) "submitter_id" else NA_character_
if (is.na(columna_paciente) || !"BRCA_Subtype_PAM50" %in% names(subtipos)) {
  stop("No se encontraron las columnas esperadas en subtipos. Revise names(subtipos).")
}
subtipos$patient_id <- substr(as.character(subtipos[[columna_paciente]]), 1, 12)
metadata_raw$subtipo <- as.character(subtipos$BRCA_Subtype_PAM50[
  match(metadata_raw$patient_id, subtipos$patient_id)
])
metadata_raw$subtipo[metadata_raw$sample_type == "Solid Tissue Normal"] <- "Normal_Tissue"

cat("\\nDistribución de muestras tras anotar subtipo:\\n")
print(table(metadata_raw$sample_type, metadata_raw$subtipo, useNA = "ifany"))

# ENFOQUE GLOBAL: todos los tumores primarios con o sin clasificación PAM50
metadata_global <- metadata_raw
metadata_global$sample_type <- factor(
  metadata_global$sample_type,
  levels = c("Solid Tissue Normal", "Primary Tumor")
)
conteos_global <- conteos_raw

# ENFOQUE ESTRATIFICADO: cuatro subtipos PAM50 + normales
grupos_analisis <- c("LumA", "LumB", "Her2", "Basal", "Normal_Tissue")
conservar_sub <- !is.na(metadata_raw$subtipo) & metadata_raw$subtipo %in% grupos_analisis
metadata_sub <- metadata_raw[conservar_sub, , drop = FALSE]
conteos_sub <- conteos_raw[, conservar_sub, drop = FALSE]
metadata_sub$subtipo <- factor(metadata_sub$subtipo, levels = grupos_analisis)
metadata_sub$subtipo <- relevel(metadata_sub$subtipo, ref = "Normal_Tissue")

stopifnot(all(colnames(conteos_global) == rownames(metadata_global)))
stopifnot(all(colnames(conteos_sub) == rownames(metadata_sub)))

# Filtrado por baja expresión separado para cada conjunto analítico
filtrar_genes <- function(mat, meta, grupo) {
  tab <- table(grupo)
  n_min <- min(tab[tab > 0])
  keep <- rowSums(mat >= 10) >= n_min
  cat("Genes antes:", nrow(mat), "| después:", sum(keep), "\\n")
  mat[keep, , drop = FALSE]
}
conteos_global <- filtrar_genes(conteos_global, metadata_global, metadata_global$sample_type)
conteos_sub <- filtrar_genes(conteos_sub, metadata_sub, metadata_sub$subtipo)

saveRDS(conteos_global, "resultados/03_conteos_global.rds")
saveRDS(metadata_global, "resultados/03_metadata_global.rds")
saveRDS(conteos_sub, "resultados/03_conteos_subtipos.rds")
saveRDS(metadata_sub, "resultados/03_metadata_subtipos.rds")
# -----------------------------------------------------------------------------

# 4. ANÁLISIS DE EXPRESIÓN DIFERENCIAL — ENFOQUE GLOBAL
# -----------------------------------------------------------------------------
# Todas las muestras tumorales (sin distinguir subtipo) vs. tejido normal

# 4.1 Construir el objeto DESeqDataSet

dds_global <- DESeqDataSetFromMatrix(
  countData = conteos_global,
  colData   = metadata_global,
  design    = ~ sample_type
)

cat("Objeto creado:", dim(dds_global), "\n")
cat("Nivel de referencia:", levels(dds_global$sample_type)[1], "\n")


# 4.2 Ejecutar DESeq2

dds_global <- DESeq(dds_global)

saveRDS(dds_global, "resultados/04_dds_global.rds")


# 4.3 Extraer resultados

res_global <- results(
  dds_global,
  contrast = c("sample_type", "Primary Tumor", "Solid Tissue Normal"),
  alpha    = 0.05
)

res_global <- lfcShrink(
  dds_global,
  contrast = c("sample_type", "Primary Tumor", "Solid Tissue Normal"),
  res      = res_global,
  type     = "ashr"
)

summary(res_global)


# 4.4 Convertir a data frame y aplicar criterios de significancia

df_global <- as.data.frame(res_global)
df_global$gene_id <- rownames(df_global)

# Anexar símbolo del gen
df_global$gene_name <- rowData(datos_brca)$gene_name[
  match(df_global$gene_id, rowData(datos_brca)$gene_id)
]

# Criterios: padj < 0.05 y |log2FC| > 1
deg_global <- df_global %>%
  filter(!is.na(padj), padj < 0.05, abs(log2FoldChange) > 1) %>%
  arrange(padj)

cat("\n--- RESULTADOS ENFOQUE GLOBAL ---\n")
cat("Genes analizados:            ", nrow(df_global), "\n")
cat("Genes diferencialmente expr.:", nrow(deg_global), "\n")
cat("  Sobreexpresados:  ", sum(deg_global$log2FoldChange > 0), "\n")
cat("  Subexpresados:    ", sum(deg_global$log2FoldChange < 0), "\n")

cat("\nTop 10 genes por significancia:\n")
print(head(deg_global[, c("gene_name", "log2FoldChange", "padj")], 10))


# 4.5 Guardar resultados

write.csv(df_global,  "resultados/04_global_completo.csv", row.names = FALSE)
write.csv(deg_global, "resultados/04_global_significativos.csv", row.names = FALSE)
saveRDS(df_global,  "resultados/04_df_global.rds")
saveRDS(deg_global, "resultados/04_deg_global.rds")
# -----------------------------------------------------------------------------


# 5. ANÁLISIS DE EXPRESIÓN DIFERENCIAL — ENFOQUE ESTRATIFICADO
# -----------------------------------------------------------------------------
# Cada subtipo molecular vs. tejido normal, de forma independiente

# 5.1 Construir el objeto con diseño por subtipo

dds_sub <- DESeqDataSetFromMatrix(
  countData = conteos_sub,
  colData   = metadata_sub,
  design    = ~ subtipo
)

cat("Niveles del factor subtipo:\n")
print(levels(dds_sub$subtipo))


# 5.2 Ejecutar DESeq2 

dds_sub <- DESeq(dds_sub)

saveRDS(dds_sub, "resultados/05_dds_estratificado.rds")


# 5.3 Extraer resultados para cada subtipo

subtipos_analizar <- c("LumA", "LumB", "Her2", "Basal")

resultados_sub <- list()
deg_sub        <- list()

for (st in subtipos_analizar) {
  
  cat("\n=== Procesando:", st, "===\n")
  
  res <- results(
    dds_sub,
    contrast = c("subtipo", st, "Normal_Tissue"),
    alpha    = 0.05
  )
  
  res <- lfcShrink(
    dds_sub,
    contrast = c("subtipo", st, "Normal_Tissue"),
    res      = res,
    type     = "ashr"
  )
  
  df <- as.data.frame(res)
  df$gene_id   <- rownames(df)
  df$gene_name <- rowData(datos_brca)$gene_name[
    match(df$gene_id, rowData(datos_brca)$gene_id)
  ]
  df$subtipo <- st
  
  deg <- df %>%
    filter(!is.na(padj), padj < 0.05, abs(log2FoldChange) > 1) %>%
    arrange(padj)
  
  resultados_sub[[st]] <- df
  deg_sub[[st]]        <- deg
  
  cat("  Genes diferencialmente expresados:", nrow(deg), "\n")
  cat("    Sobreexpresados:", sum(deg$log2FoldChange > 0), "\n")
  cat("    Subexpresados:  ", sum(deg$log2FoldChange < 0), "\n")
  
  write.csv(df,  paste0("resultados/05_", st, "_completo.csv"), row.names = FALSE)
  write.csv(deg, paste0("resultados/05_", st, "_significativos.csv"), row.names = FALSE)
}


# 5.4 Tabla resumen

resumen <- data.frame(
  Enfoque = c("Global", subtipos_analizar),
  n_muestras_tumor = c(
    sum(metadata_global$sample_type == "Primary Tumor"),
    sapply(subtipos_analizar, function(x) sum(metadata_sub$subtipo == x))
  ),
  n_DEG = c(
    nrow(deg_global),
    sapply(subtipos_analizar, function(x) nrow(deg_sub[[x]]))
  ),
  n_sobreexpresados = c(
    sum(deg_global$log2FoldChange > 0),
    sapply(subtipos_analizar, function(x) sum(deg_sub[[x]]$log2FoldChange > 0))
  ),
  n_subexpresados = c(
    sum(deg_global$log2FoldChange < 0),
    sapply(subtipos_analizar, function(x) sum(deg_sub[[x]]$log2FoldChange < 0))
  )
)

cat("\n--- TABLA RESUMEN ---\n")
print(resumen)

write.csv(resumen, "resultados/05_resumen_DEG.csv", row.names = FALSE)


# 5.5 Guardar listas completas

saveRDS(resultados_sub, "resultados/05_resultados_subtipos.rds")
saveRDS(deg_sub,        "resultados/05_deg_subtipos.rds")
# -----------------------------------------------------------------------------


# 6. COMPARACIÓN ENTRE ENFOQUES
# -----------------------------------------------------------------------------
# Cuantificar cuánto se solapan los genes identificados por cada estrategia

# 6.1 Extraer los conjuntos de genes

genes_global <- deg_global$gene_id

genes_por_subtipo <- lapply(deg_sub, function(x) x$gene_id)
names(genes_por_subtipo) <- names(deg_sub)

# Unión de todos los genes detectados en el enfoque estratificado
genes_estratificado <- unique(unlist(genes_por_subtipo))

cat("Genes detectados:\n")
cat("  Enfoque global:       ", length(genes_global), "\n")
cat("  Enfoque estratificado:", length(genes_estratificado), "(unión de los 4 subtipos)\n")


# 6.2 Solapamiento global vs. estratificado

compartidos    <- intersect(genes_global, genes_estratificado)
solo_global    <- setdiff(genes_global, genes_estratificado)
solo_estratif  <- setdiff(genes_estratificado, genes_global)

cat("\n--- SOLAPAMIENTO ---\n")
cat("Compartidos:                   ", length(compartidos),
    sprintf("(%.1f%%)", 100 * length(compartidos) / length(union(genes_global, genes_estratificado))), "\n")
cat("Exclusivos del global:         ", length(solo_global), "\n")
cat("Exclusivos del estratificado:  ", length(solo_estratif), "\n")


# 6.3 Genes detectados en TODOS los subtipos (señal compartida)

genes_en_todos <- Reduce(intersect, genes_por_subtipo)

cat("\nGenes detectados en los 4 subtipos:", length(genes_en_todos), "\n")


# 6.4 Genes exclusivos de cada subtipo
# Estos son el hallazgo central (señal que el análisis global no captura)

cat("\n--- GENES EXCLUSIVOS POR SUBTIPO ---\n")

exclusivos_subtipo <- list()

for (st in names(genes_por_subtipo)) {
  otros <- unlist(genes_por_subtipo[names(genes_por_subtipo) != st])
  excl  <- setdiff(genes_por_subtipo[[st]], otros)
  
  exclusivos_subtipo[[st]] <- excl
  
  # ¿Cuántos de esos NO aparecen en el análisis global?
  perdidos <- setdiff(excl, genes_global)
  
  cat(sprintf("%-6s: %4d exclusivos | %4d ausentes del análisis global\n",
              st, length(excl), length(perdidos)))
}

# 6.5 Tabla maestra de comparación
# Para cada gen (log2FC en cada análisis y si fue significativo)

tabla_comparacion <- data.frame(
  gene_id   = df_global$gene_id,
  gene_name = df_global$gene_name,
  log2FC_global = df_global$log2FoldChange,
  padj_global   = df_global$padj,
  sig_global    = df_global$gene_id %in% genes_global
)

for (st in subtipos_analizar) {
  df_st <- resultados_sub[[st]]
  idx   <- match(tabla_comparacion$gene_id, df_st$gene_id)
  
  tabla_comparacion[[paste0("log2FC_", st)]] <- df_st$log2FoldChange[idx]
  tabla_comparacion[[paste0("padj_", st)]]   <- df_st$padj[idx]
  tabla_comparacion[[paste0("sig_", st)]]    <- tabla_comparacion$gene_id %in% genes_por_subtipo[[st]]
}

# Número de subtipos en que cada gen fue significativo
tabla_comparacion$n_subtipos_sig <- rowSums(
  tabla_comparacion[, paste0("sig_", subtipos_analizar)]
)

# Clasificación
tabla_comparacion$categoria <- with(tabla_comparacion, ifelse(
  sig_global & n_subtipos_sig == 4, "Compartido (global + todos)",
  ifelse(sig_global & n_subtipos_sig > 0, "Global + algunos subtipos",
         ifelse(sig_global & n_subtipos_sig == 0, "Solo global",
                ifelse(!sig_global & n_subtipos_sig > 0, "Solo estratificado (señal diluida)",
                       "No significativo")))))

cat("\n--- CLASIFICACIÓN DE GENES ---\n")
print(table(tabla_comparacion$categoria))

write.csv(tabla_comparacion, "resultados/06_tabla_comparacion.csv", row.names = FALSE)


# 6.6 Genes con señal diluida (este es el hallazgo central)

senal_diluida <- tabla_comparacion %>%
  filter(categoria == "Solo estratificado (señal diluida)") %>%
  arrange(desc(n_subtipos_sig))

cat("\nGenes con señal diluida:", nrow(senal_diluida), "\n")
cat("\nTop 15 por magnitud de cambio en algún subtipo:\n")

columnas_lfc <- paste0("log2FC_", subtipos_analizar)

matriz_lfc <- abs(
  as.matrix(senal_diluida[, columnas_lfc, drop = FALSE])
)

senal_diluida$max_lfc <- Reduce(
  pmax,
  lapply(seq_len(ncol(matriz_lfc)), function(i) matriz_lfc[, i]),
  na.rm = TRUE
)

senal_diluida$max_lfc[
  !is.finite(senal_diluida$max_lfc)
] <- NA_real_

print(head(
  senal_diluida[order(-senal_diluida$max_lfc),
                c("gene_name", "log2FC_global", "log2FC_Basal", "log2FC_Her2", "n_subtipos_sig")],
  15
))

write.csv(senal_diluida, "resultados/06_senal_diluida.csv", row.names = FALSE)


#  6.7 Guardar conjuntos para las visualizaciones

saveRDS(list(
  global        = genes_global,
  por_subtipo   = genes_por_subtipo,
  compartidos   = compartidos,
  solo_global   = solo_global,
  solo_estratif = solo_estratif,
  en_todos      = genes_en_todos,
  exclusivos    = exclusivos_subtipo
), "resultados/06_conjuntos_genes.rds")


# 7. ENRIQUECIMIENTO FUNCIONAL GO Y KEGG

# 7.1 Definir genes de fondo por contraste. Genes efectivamente evaluados, no todos los genes del genoma. 
contrastes <- c(list(Global = list(completo = df_global, deg = deg_global)),
               lapply(subtipos_analizar, function(st) {
                 list(completo = resultados_sub[[st]], deg = deg_sub[[st]])
               }))
names(contrastes) <- c("Global", subtipos_analizar)

todos_ensembl <- unique(sub("\\..*$", "", unlist(lapply(contrastes, function(x) {
  x$completo$gene_id[!is.na(x$completo$padj)]
}))))
todos_ensembl <- todos_ensembl[!is.na(todos_ensembl) & nzchar(todos_ensembl)]

mapa_entrez <- AnnotationDbi::mapIds(
  org.Hs.eg.db, keys = todos_ensembl, keytype = "ENSEMBL",
  column = "ENTREZID", multiVals = "first"
)
cat("\nEnsembl evaluados:", length(todos_ensembl),
    "| con Entrez:", sum(!is.na(mapa_entrez)), "\n")
write.csv(data.frame(ensembl = names(mapa_entrez), entrez = unname(mapa_entrez)),
          "resultados/07_mapeo_ensembl_entrez.csv", row.names = FALSE)

a_entrez <- function(ids) {
  ids <- unique(sub("\\..*$", "", as.character(ids)))
  ids <- ids[!is.na(ids) & nzchar(ids)]
  unique(unname(na.omit(mapa_entrez[ids])))
}

# 7.2 Analizar por separado genes sobreexpresados y subexpresados.

tablas_enriquecimiento <- list()
resumen_conjuntos <- list()

for (enfoque in names(contrastes)) {
  completo <- contrastes[[enfoque]]$completo
  deg <- contrastes[[enfoque]]$deg
  fondo <- a_entrez(completo$gene_id[!is.na(completo$padj)])

  for (direccion in c("Sobreexpresados", "Subexpresados")) {
    ids <- if (direccion == "Sobreexpresados") {
      deg$gene_id[!is.na(deg$log2FoldChange) & deg$log2FoldChange > 0]
    } else {
      deg$gene_id[!is.na(deg$log2FoldChange) & deg$log2FoldChange < 0]
    }
    genes <- intersect(a_entrez(ids), fondo)
    resumen_conjuntos[[length(resumen_conjuntos) + 1L]] <- data.frame(
      Enfoque = enfoque, Direccion = direccion,
      n_DEG = length(unique(ids)), n_Entrez = length(genes),
      n_fondo_Entrez = length(fondo)
    )
    cat("\n", enfoque, direccion, ":", length(genes),
        "Entrez DEG /", length(fondo), "Entrez evaluados\n")

    if (length(genes) < 10L || length(fondo) < 10L) next

    analisis <- list(
      GO_BP = function() clusterProfiler::enrichGO(
        gene = genes, universe = fondo, OrgDb = org.Hs.eg.db,
        keyType = "ENTREZID", ont = "BP", pAdjustMethod = "BH",
        pvalueCutoff = 0.05, qvalueCutoff = 1, readable = FALSE
      ),
      KEGG = function() clusterProfiler::enrichKEGG(
        gene = genes, universe = fondo, organism = "hsa",
        keyType = "ncbi-geneid", pAdjustMethod = "BH",
        pvalueCutoff = 0.05, qvalueCutoff = 1
      )
    )

    for (base in names(analisis)) {
      resultado <- tryCatch(analisis[[base]](), error = function(e) {
        warning(sprintf("Enriquecimiento %s (%s, %s): %s",
                        base, enfoque, direccion, conditionMessage(e)))
        NULL
      })
      tabla <- if (is.null(resultado)) data.frame() else as.data.frame(resultado)
      if (nrow(tabla) == 0L) next
      tabla <- tabla[!is.na(tabla$p.adjust) & tabla$p.adjust < 0.05, , drop = FALSE]
      if (nrow(tabla) == 0L) next

      tabla$Enfoque <- enfoque
      tabla$Direccion <- direccion
      tabla$Base <- base
      tabla <- tabla[order(tabla$p.adjust), , drop = FALSE]
      clave <- paste(base, enfoque, direccion, sep = "_")
      tablas_enriquecimiento[[clave]] <- tabla
      write.csv(tabla, paste0("resultados/07_", clave, ".csv"), row.names = FALSE)
      cat("  ", base, ":", nrow(tabla), "términos (FDR < 0.05)\n")
    }
  }
}

write.csv(dplyr::bind_rows(resumen_conjuntos),
          "resultados/07_resumen_genes_mapeados.csv", row.names = FALSE)

# 7.3 Tabla conjunta para comparar funciones entre global y subtipos.

if (length(tablas_enriquecimiento) > 0L) {
  enriquecimiento <- dplyr::bind_rows(tablas_enriquecimiento)
  write.csv(enriquecimiento, "resultados/07_enriquecimiento_completo.csv",
            row.names = FALSE)

  top_terminos <- enriquecimiento %>%
    dplyr::group_by(Base, Direccion, Enfoque) %>%
    dplyr::slice_min(order_by = p.adjust, n = 10, with_ties = FALSE) %>%
    dplyr::ungroup()
  write.csv(top_terminos, "resultados/07_top10_por_contraste.csv",
            row.names = FALSE)

  resumen_rutas <- enriquecimiento %>%
    dplyr::count(Base, Direccion, Enfoque, name = "n_terminos_FDR05")
  write.csv(resumen_rutas, "resultados/07_resumen_rutas.csv",
            row.names = FALSE)

  # Dos mapas de calor: hasta 12 términos de mayor evidencia por dirección.
  for (base in c("GO_BP", "KEGG")) {
    datos_base <- enriquecimiento[enriquecimiento$Base == base, , drop = FALSE]
    if (nrow(datos_base) == 0L) next
    mejores <- datos_base %>%
      dplyr::group_by(Direccion, ID, Description) %>%
      dplyr::summarise(min_fdr = min(p.adjust), .groups = "drop") %>%
      dplyr::group_by(Direccion) %>%
      dplyr::slice_min(min_fdr, n = 12, with_ties = FALSE) %>%
      dplyr::ungroup()
    grafico <- dplyr::inner_join(datos_base, mejores[, c("Direccion", "ID")],
                                by = c("Direccion", "ID"))
    grafico$etiqueta <- paste0(grafico$Description, " [", grafico$ID, "]")
    grafico$Enfoque <- factor(grafico$Enfoque,
                              levels = c("Global", subtipos_analizar))
    p <- ggplot(grafico, aes(x = Enfoque, y = etiqueta,
                             fill = -log10(pmax(p.adjust, .Machine$double.xmin)))) +
      geom_tile(color = "white") +
      facet_wrap(~ Direccion, ncol = 1, scales = "free_y") +
      scale_fill_gradient(low = "#CD6E94", high = "#F4A652",
                          name = expression(-log[10](FDR))) +
      labs(title = paste("Rutas enriquecidas:", base), x = NULL, y = NULL,
           caption = "Solo términos con FDR < 0.05; celda vacía = no enriquecido") +
      theme_minimal(base_size = 10) +
      theme(axis.text.x = element_text(angle = 35, hjust = 1),
            panel.grid = element_blank())
    ggsave(paste0("figuras/07_mapa_rutas_", base, ".png"), p,
           width = 12, height = 10, dpi = 300, limitsize = FALSE)
  }
} else {
  message("No hubo términos enriquecidos con FDR < 0.05; revisar mapeo y resultados DEG.")
}



# 8. CONTRASTE CON BIOMARCADORES REPORTADOS EN LA LITERATURA
# -----------------------------------------------------------------------------
# No se codifica una lista de biomarcadores inventada. Complete el archivo
# resultados/biomarcadores_referencia.csv con genes documentados en artículos
# revisados por el equipo. Debe contener una columna llamada gene_name y otra
# llamada referencia (DOI o cita APA breve).
#
# Ejemplo de estructura:
# gene_name,referencia
# GEN1,Autor et al. (Año), DOI...
#
# El archivo se puede crear manualmente y luego ejecutar este bloque.

archivo_biomarcadores <- "resultados/biomarcadores_referencia.csv"
if (file.exists(archivo_biomarcadores)) {
  biomarcadores <- read.csv(archivo_biomarcadores, stringsAsFactors = FALSE)
  if (!"gene_name" %in% names(biomarcadores)) {
    warning("El archivo de biomarcadores debe incluir la columna 'gene_name'.")
  } else {
    biomarcadores$gene_name <- toupper(trimws(biomarcadores$gene_name))
    tabla_comparacion$gene_name_upper <- toupper(trimws(tabla_comparacion$gene_name))
    contraste_biomarcadores <- tabla_comparacion %>%
      filter(gene_name_upper %in% biomarcadores$gene_name) %>%
      left_join(biomarcadores, by = c("gene_name_upper" = "gene_name"))
    write.csv(contraste_biomarcadores,
              "resultados/08_contraste_biomarcadores.csv", row.names = FALSE)
    cat("Biomarcadores de referencia encontrados en resultados:",
        nrow(contraste_biomarcadores), "\\n")
  }
} else {
  message("Se omitió el contraste de biomarcadores: cree ",
          archivo_biomarcadores,
          " con genes respaldados por la literatura y su referencia.")
}

# 9. GRÁFICOS BÁSICOS PARA PRESENTAR LOS RESULTADOS
# -----------------------------------------------------------------------------
# Gráfico resumen de genes significativos por enfoque
resumen_plot <- resumen
resumen_plot$Enfoque <- factor(resumen_plot$Enfoque,
                               levels = c("Global", subtipos_analizar))
p_resumen <- ggplot(resumen_plot, aes(x = Enfoque, y = n_DEG)) +
  geom_col() +
  labs(title = "Genes diferencialmente expresados por enfoque",
       x = "Enfoque de análisis", y = "Número de genes (FDR < 0.05 y |log2FC| > 1)") +
  theme_minimal()
ggsave("figuras/09_resumen_genes_DEG.png", p_resumen,
       width = 8, height = 5, dpi = 300)

# Exportar sesión y versiones para facilitar la reproducibilidad
writeLines(capture.output(sessionInfo()), "resultados/sessionInfo.txt")
cat("\\nAnálisis terminado. Revise las carpetas 'resultados' y 'figuras'.\\n")
