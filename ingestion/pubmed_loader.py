import argparse
import http.client
import json
import os
import re
import time
from pathlib import Path
from urllib.error import HTTPError

from Bio import Entrez
from tqdm import tqdm


OUTPUT_DIR = Path("data/pubmed")
OUTPUT_DIR.mkdir(parents=True, exist_ok=True)

GENES = ["BRCA1", "BRCA2", "TP53", "PIK3CA", "ERBB2"]

RETMAX = int(os.getenv("PUBMED_RETMAX", "50"))
BATCH_SIZE = int(os.getenv("PUBMED_BATCH_SIZE", "10"))
MAX_RETRIES = int(os.getenv("PUBMED_MAX_RETRIES", "3"))

DEFAULT_EMAIL = "noreply@example.com"
Entrez.email = os.getenv("NCBI_EMAIL", DEFAULT_EMAIL)
NCBI_API_KEY = os.getenv("NCBI_API_KEY")
if NCBI_API_KEY:
    Entrez.api_key = NCBI_API_KEY
REQUEST_DELAY = float(
    os.getenv("PUBMED_REQUEST_DELAY", "0.12" if NCBI_API_KEY else "0.34")
)


def search_pubmed(gene: str, retmax: int) -> list[str]:
    query = f"({gene}[Title/Abstract]) AND (breast cancer[Title/Abstract])"
    handle = Entrez.esearch(db="pubmed", term=query, retmax=retmax)
    record = Entrez.read(handle)
    handle.close()
    return record.get("IdList", [])


def chunked(items: list[str], size: int) -> list[list[str]]:
    return [items[i : i + size] for i in range(0, len(items), size)]


def extract_year(article: dict) -> str:
    pub_date = article.get("Journal", {}).get("JournalIssue", {}).get("PubDate", {})
    if pub_date.get("Year"):
        return str(pub_date["Year"])

    medline_date = str(pub_date.get("MedlineDate", ""))
    match = re.search(r"\b(19|20)\d{2}\b", medline_date)
    return match.group(0) if match else ""


def parse_article(pubmed_article: dict, gene: str) -> dict | None:
    article = pubmed_article.get("MedlineCitation", {}).get("Article", {})
    if not article:
        return None

    pmid = str(pubmed_article.get("MedlineCitation", {}).get("PMID", ""))
    title = str(article.get("ArticleTitle", ""))

    abstract_parts = article.get("Abstract", {}).get("AbstractText", [])
    abstract = " ".join(str(part) for part in abstract_parts).strip()
    if not abstract:
        return None

    journal = str(article.get("Journal", {}).get("Title", ""))
    year = extract_year(article)

    return {
        "pmid": pmid,
        "title": title,
        "abstract": abstract,
        "journal": journal,
        "year": year,
        "gene": gene,
    }


def fetch_abstract_batch(pmids: list[str], gene: str) -> list[dict]:
    if not pmids:
        return []

    for attempt in range(MAX_RETRIES):
        try:
            time.sleep(REQUEST_DELAY)
            handle = Entrez.efetch(
                db="pubmed",
                id=",".join(pmids),
                rettype="abstract",
                retmode="xml",
            )
            record = Entrez.read(handle)
            handle.close()

            parsed_records = []
            for pubmed_article in record.get("PubmedArticle", []):
                parsed = parse_article(pubmed_article, gene)
                if parsed:
                    parsed_records.append(parsed)
            return parsed_records

        except (HTTPError, http.client.RemoteDisconnected):
            if attempt == MAX_RETRIES - 1:
                return []
            time.sleep(2)
        except Exception:
            return []

    return []


def ingest_gene(gene: str, retmax: int, batch_size: int) -> None:
    pmids = search_pubmed(gene, retmax=retmax)
    batches = chunked(pmids, max(1, batch_size))

    records = []
    for batch in tqdm(batches, desc=f"Fetching {gene}", unit="batch"):
        records.extend(fetch_abstract_batch(batch, gene))

    output_path = OUTPUT_DIR / f"{gene}.json"
    with open(output_path, "w", encoding="utf-8") as file:
        json.dump(records, file, indent=2, ensure_ascii=False)

    print(f"[OK] {gene}: saved {len(records)} records to {output_path}")


def parse_args() -> argparse.Namespace:
    parser = argparse.ArgumentParser(description="Fetch PubMed abstracts by gene.")
    parser.add_argument("--genes", nargs="+", default=GENES)
    parser.add_argument("--retmax", type=int, default=RETMAX)
    parser.add_argument("--batch-size", type=int, default=BATCH_SIZE)
    return parser.parse_args()


if __name__ == "__main__":
    args = parse_args()

    if Entrez.email == DEFAULT_EMAIL:
        print(
            "[WARN] Using fallback NCBI email. Set NCBI_EMAIL for compliant usage."
        )

    for gene in args.genes:
        ingest_gene(gene, retmax=args.retmax, batch_size=args.batch_size)
