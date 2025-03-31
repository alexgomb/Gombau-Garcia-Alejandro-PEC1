Tiene 63 features (los metabolitos) y 77 muestras (los pacientes). La matriz de datos se almacena en el assay "counts", los nombres de las filas corresponden a cada metabolito (por ejemplo, X1.6.Anhydro.beta.D.glucose, X1.Methylnicotinamide, etc.) y en el rowData se guarda información adicional, en este caso solo el nombre del metabolito bajo la columna "Metabolito". Por otro lado, colData contiene dos variables: "Patient.ID" y "Muscle.loss", que representan los metadatos de cada paciente.


class: SummarizedExperiment 
dim: 63 77 
metadata(0):
assays(1): counts
rownames(63): X1.6.Anhydro.beta.D.glucose X1.Methylnicotinamide ... pi.Methylhistidine tau.Methylhistidine
rowData names(1): Metabolito
colnames(77): PIF_178 PIF_087 ... NETL_003_V1 NETL_003_V2
colData names(2): Patient.ID Muscle.loss
