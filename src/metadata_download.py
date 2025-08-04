import requests
import json

def obtener_metadatos_ena(srr_id):
    """
    Obtiene los metadatos de ENA para un identificador SRR dado.

    Args:
        srr_id (str): El identificador SRR.

    Returns:
        list: Una lista de diccionarios con los metadatos, o None si ocurre un error.
    """
    url = f"https://www.ebi.ac.uk/ena/portal/api/filereport?accession={srr_id}&result=read_run&fields=study_accession,sample_accession,experiment_accession,run_accession,tax_id,scientific_name,instrument_platform,instrument_model,library_layout,library_source,library_selection,read_count,base_count,experiment_title,sample_title&format=json"
    try:
        response = requests.get(url)
        response.raise_for_status()
        return response.json()
    except requests.exceptions.RequestException as e:
        print(f"Error al obtener metadatos para {srr_id}: {e}")
        return None

def guardar_metadatos(metadatos, archivo_salida):
    """
    Guarda los metadatos en un archivo JSON.

    Args:
        metadatos (list): Una lista de diccionarios con los metadatos.
        archivo_salida (str): La ruta del archivo de salida.
    """
    with open(archivo_salida, "w") as f:
        json.dump(metadatos, f, indent=4)

if __name__ == "__main__":
    srr_ids = ["SRR8615901", "SRR8615876", "SRR8615924"]
    metadatos_totales = []
    for srr_id in srr_ids:
        metadatos = obtener_metadatos_ena(srr_id)
        if metadatos:
            metadatos_totales.extend(metadatos)
    guardar_metadatos(metadatos_totales, "metadatos_ena.json")