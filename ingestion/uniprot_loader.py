import argparse
import json
import os
import time
from pathlib import Path

import requests
from requests.adapters import HTTPAdapter
from tqdm import tqdm
from urllib3.util.retry import Retry


GENES = ["BRCA1", "BRCA2", "TP53", "PIK3CA", "ERBB2"]
OUTPUT_DIR = Path("data/uniprot")
BASE_URL = "https://rest.uniprot.org/uniprotkb/search"
REQUEST_DELAY = float(os.getenv("UNIPROT_REQUEST_DELAY", "0.5"))

OUTPUT_DIR.mkdir(parents=True, exist_ok=True)

HEADERS = {"User-Agent": "genai-breast-cancer-genomics/1.0"}


def build_session() -> requests.Session:
    retry = Retry(
        total=3,
        backoff_factor=0.7,
        status_forcelist=[429, 500, 502, 503, 504],
        allowed_methods=["GET"],
    )
    adapter = HTTPAdapter(max_retries=retry)

    session = requests.Session()
    session.headers.update(HEADERS)
    session.mount("https://", adapter)
    session.mount("http://", adapter)
    return session


def fetch_uniprot(gene: str, session: requests.Session) -> dict | None:
    query = f"gene_exact:{gene} AND reviewed:true"
    params = {"query": query, "format": "json", "size": 1}

    try:
        response = session.get(BASE_URL, params=params, timeout=20)
        response.raise_for_status()
        data = response.json()
    except Exception as exc:
        print(f"[ERROR] UniProt request failed for {gene}: {exc}")
        return None

    results = data.get("results", [])
    if not results:
        print(f"[WARN] No reviewed UniProt entry found for {gene}")
        return None

    entry = results[0]
    protein_desc = entry.get("proteinDescription", {})
    recommended = protein_desc.get("recommendedName", {})
    protein_name = recommended.get("fullName", {}).get("value", "")

    functions: list[str] = []
    subcellular: list[str] = []
    diseases: list[str] = []

    for comment in entry.get("comments", []):
        comment_type = comment.get("commentType")

        if comment_type == "FUNCTION":
            for text in comment.get("texts", []):
                functions.append(text.get("value", ""))

        if comment_type == "SUBCELLULAR LOCATION":
            for location_data in comment.get("subcellularLocations", []):
                location = location_data.get("location", {}).get("value", "")
                if location:
                    subcellular.append(location)

        if comment_type == "DISEASE":
            description = comment.get("disease", {}).get("description", "")
            if description:
                diseases.append(description)

    if not functions:
        for comment in entry.get("comments", []):
            if comment.get("commentType") == "SIMILARITY":
                functions.append("Function inferred from sequence similarity")

    return {
        "gene": gene,
        "uniprot_id": entry.get("primaryAccession", ""),
        "protein_name": protein_name,
        "function": " ".join(functions),
        "subcellular_location": sorted(set(subcellular)),
        "disease_association": sorted(set(diseases)),
    }


def parse_args() -> argparse.Namespace:
    parser = argparse.ArgumentParser(description="Fetch UniProt gene annotations.")
    parser.add_argument("--genes", nargs="+", default=GENES)
    return parser.parse_args()


if __name__ == "__main__":
    args = parse_args()
    session = build_session()

    for gene in tqdm(args.genes, desc="Fetching UniProt annotations"):
        record = fetch_uniprot(gene, session)
        if record:
            output_path = OUTPUT_DIR / f"{gene}.json"
            with open(output_path, "w", encoding="utf-8") as file:
                json.dump(record, file, indent=2, ensure_ascii=False)
        time.sleep(REQUEST_DELAY)
