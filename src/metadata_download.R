# Import libraries
library(recount3)

# Obtener la lista de proyectos disponibles en recount3
human_projects <- available_projects()

# Mostrar los valores únicos de la columna "project_home" para ver de dónde provienen los datos
unique(human_projects$project_home)

# Filtrar solo los proyectos que pertenecen a TCGA (The Cancer Genome Atlas)
human_projects[human_projects$project_home == "data_sources/tcga", ]

# Seleccionar la información del proyecto específico "LAML" (Leucemia Mieloide Aguda)
# y asegurarse de que proviene de la fuente de datos correcta
project_info <- subset(
    human_projects,
    project == "LAML" & project_type == "data_sources"
)

# Crear un objeto RangedSummarizedExperiment (RSE) con los datos de LAML
rse_LAML <- create_rse(project_info)

# Calcular los conteos de lectura (read counts) y almacenarlos en el assay "counts" 
# dentro del objeto RSE
assay(rse_LAML, "counts") <- compute_read_counts(rse_LAML)
