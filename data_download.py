'''
NAME
    data_download

VERSION
    1.0


AUTHOR
    Gabriel Ramirez-Vilchis

DESCRIPTION
    Gets data from the SRA database using the BioPython library

CATEGORY
    Data download

USAGE
    python src/data_download.py

ARGUMENTS
    None

SEE ALSO

'''

# Importar librerías
from Bio import Entrez

# Registrar un correo electrónico
Entrez.email = "gramirez@lcg.unam.mx"

# Extraer una entrada de GenBank desde la base de datos SRA
with Entrez.efetch(db="sra", 
                   id="SRR8615998",
                   rettype="xml", 
                   retmode="text") as handle:
    data = handle.read()  # Leer los datos

# Decodificar si es necesario
if isinstance(data, bytes):  
    data = data.decode("utf-8")  # Convertir de bytes a string

# Guardar los datos en un archivo
with open("SRR8615998.xml", "w", encoding="utf-8") as output:
    output.write(data)