'''
NAME
    data_download

VERSION
    1.1


AUTHOR
    Gabriel Ramirez-Vilchis

DESCRIPTION
    Gets data from the SRA database using the BioPython library

CATEGORY
    Data download

USAGE
    python data_download.py

ARGUMENTS
    None

SEE ALSO

'''

# Importar librerías
from Bio import Entrez

# Registrar un correo electrónico
Entrez.email = "gramirez@lcg.unam.mx"

# Extraer una entrada de GenBank desde la base de datos SRA
sra_id = "SRR8615998"
with Entrez.efetch(db="sra", 
                   id=sra_id,
                   rettype="xml", 
                   retmode="text") as handle:
    data = handle.read()  # Leer los datos

# Decodificar si es necesario
if isinstance(data, bytes):  
    data = data.decode("utf-8")  # Convertir de bytes a string

# Guardar los datos en un archivo
with open(f"../results/{sra_id}.xml", "w", encoding="utf-8") as output:
    output.write(data)
