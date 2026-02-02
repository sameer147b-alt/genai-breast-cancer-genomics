from Bio import Entrez
import json
import os
import time
from tqdm import tqdm
from urllib.error import HTTPError
import http.client

# =========================
# REQUIRED BY NCBI
# =========================
Entrez.email = "sameer.147b@gmail.com"  # use your real email

# =========================
# CONFIG
# =========================
OUTPUT_DIR = "data/pubmed"
os.makedirs(OUTPUT_DIR, exist_ok=True)

GENES = ["BRCA1", "BRCA2", "TP53", "PIK3CA", "ERBB2"]

REQUEST_DELAY = 0.4   # seconds (NCBI-safe)
MAX_RETRIES = 3       # retry on connection drop
RETMAX = 50           # keep quality high

# =========================
# PUBMED SEARCH
# =========================
def search_pubmed(query, retmax=RETMAX):
    handle = Entrez.esearch(
        db="pubmed",
        term=query,
        retmax=retmax
    )
    record = Entrez.read(handle)
    handle.close()
    return record["IdList"]

# =========================
# FETCH ABSTRACT (ROBUST)
# =========================
def fetch_abstract(pmid):
    for attempt in range(MAX_RETRIES):
        try:
            time.sleep(REQUEST_DELAY)

            handle = Entrez.efetch(
                db="pubmed",
                id=pmid,
                rettype="abstract",
                retmode="xml"
            )
            record = Entrez.read(handle)
            handle.close()

            article = record["PubmedArticle"][0]["MedlineCitation"]["Article"]

            title = article.get("ArticleTitle", "")
            abstract = " ".join(
                article.get("Abstract", {}).get("AbstractText", [])
            )
            journal = article["Journal"]["Title"]
            year = article["Journal"]["JournalIssue"]["PubDate"].get("Year", "")

            return {
                "pmid": pmid,
                "title": title,
                "abstract": abstract,
                "journal": journal,
                "year": year
            }

        except (HTTPError, http.client.RemoteDisconnected):
            if attempt == MAX_RETRIES - 1:
                return None
            time.sleep(2)  # wait before retry

        except Exception:
            return None

# =========================
# INGEST PER GENE
# =========================
def ingest_gene(gene):
    query = f"({gene}[Title/Abstract]) AND (breast cancer[Title/Abstract])"
    pmids = search_pubmed(query)

    records = []

    for pmid in tqdm(pmids, desc=f"Fetching {gene}"):
        data = fetch_abstract(pmid)
        if data and data["abstract"]:
            data["gene"] = gene
            records.append(data)

    output_path = os.path.join(OUTPUT_DIR, f"{gene}.json")
    with open(output_path, "w", encoding="utf-8") as f:
        json.dump(records, f, indent=2)

# =========================
# MAIN
# =========================
if __name__ == "__main__":
    for gene in GENES:
        ingest_gene(gene)
