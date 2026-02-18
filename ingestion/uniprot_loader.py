import requests
import json
import os
from tqdm import tqdm
import time

# ==============================
# CONFIG
# ==============================
GENES = ["BRCA1", "BRCA2", "TP53", "PIK3CA", "ERBB2"]
OUTPUT_DIR = "data/uniprot"
BASE_URL = "https://rest.uniprot.org/uniprotkb/search"

os.makedirs(OUTPUT_DIR, exist_ok=True)

HEADERS = {
    "User-Agent": "genai-breast-cancer-genomics/1.0"
}

# ==============================
# UNIPROT FETCH
# ==============================
def fetch_uniprot(gene: str):
    """
    Fetch reviewed (Swiss-Prot) UniProt entry for a gene
    """
    query = f"gene_exact:{gene} AND reviewed:true"
    params = {
        "query": query,
        "format": "json",
        "size": 1
    }

    try:
        response = requests.get(BASE_URL, params=params, headers=HEADERS, timeout=20)
        response.raise_for_status()
        data = response.json()
    except Exception as e:
        print(f"[ERROR] UniProt request failed for {gene}: {e}")
        return None

    results = data.get("results", [])
    if not results:
        print(f"[WARN] No reviewed UniProt entry found for {gene}")
        return None

    entry = results[0]

    # ==============================
    # PARSING
    # ==============================
    protein_name = ""
    functions = []
    subcellular = []
    diseases = []

    # Protein name
    protein_desc = entry.get("proteinDescription", {})
    recommended = protein_desc.get("recommendedName", {})
    protein_name = recommended.get("fullName", {}).get("value", "")

    # Comments
    for comment in entry.get("comments", []):
        ctype = comment.get("commentType")

        # Function
        if ctype == "FUNCTION":
            for txt in comment.get("texts", []):
                functions.append(txt.get("value", ""))

        # Subcellular location
        if ctype == "SUBCELLULAR LOCATION":
            for loc in comment.get("subcellularLocations", []):
                location = loc.get("location", {}).get("value", "")
                if location:
                    subcellular.append(location)

        # Disease association
        if ctype == "DISEASE":
            disease = comment.get("disease", {})
            desc = disease.get("description", "")
            if desc:
                diseases.append(desc)

    # Fallback if function missing
    if not functions:
        for comment in entry.get("comments", []):
            if comment.get("commentType") == "SIMILARITY":
                functions.append("Function inferred from sequence similarity")

    record = {
        "gene": gene,
        "uniprot_id": entry.get("primaryAccession", ""),
        "protein_name": protein_name,
        "function": " ".join(functions),
        "subcellular_location": sorted(list(set(subcellular))),
        "disease_association": sorted(list(set(diseases)))
    }

    return record


# ==============================
# MAIN
# ==============================
if __name__ == "__main__":
    for gene in tqdm(GENES, desc="Fetching UniProt annotations"):
        record = fetch_uniprot(gene)
        if record:
            out_path = os.path.join(OUTPUT_DIR, f"{gene}.json")
            with open(out_path, "w", encoding="utf-8") as f:
                json.dump(record, f, indent=2, ensure_ascii=False)
        time.sleep(1)  # polite rate limiting
