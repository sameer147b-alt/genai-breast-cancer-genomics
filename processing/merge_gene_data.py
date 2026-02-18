import json
from pathlib import Path

# =========================
# PATH CONFIG
# =========================
BASE_DIR = Path(__file__).resolve().parents[1]

PUBMED_DIR = BASE_DIR / "data" / "pubmed"
UNIPROT_DIR = BASE_DIR / "data" / "uniprot"
OUTPUT_DIR = BASE_DIR / "data" / "processed"

OUTPUT_DIR.mkdir(parents=True, exist_ok=True)

GENES = ["BRCA1", "BRCA2", "TP53", "PIK3CA", "ERBB2"]


def load_json(path):
    with open(path, "r", encoding="utf-8") as f:
        return json.load(f)


def build_llm_text(gene, uniprot, pubmed_records):
    sections = []

    sections.append(f"Gene: {gene}")
    sections.append(f"Protein name: {uniprot.get('protein_name', 'N/A')}")
    sections.append(f"UniProt ID: {uniprot.get('uniprot_id', 'N/A')}")

    if uniprot.get("function"):
        sections.append("Biological function:")
        sections.append(uniprot["function"])

    if uniprot.get("subcellular_location"):
        sections.append(
            "Subcellular location: " +
            ", ".join(uniprot["subcellular_location"])
        )

    if uniprot.get("disease_association"):
        sections.append("Disease relevance:")
        for d in uniprot["disease_association"]:
            sections.append(f"- {d}")

    sections.append("Relevant breast cancer literature:")

    for paper in pubmed_records[:10]:
        sections.append(f"Title: {paper.get('title', '')}")
        sections.append(f"Journal: {paper.get('journal', '')} ({paper.get('year', '')})")
        sections.append(f"Abstract: {paper.get('abstract', '')}")

    return "\n".join(sections)


def merge_gene(gene):
    pubmed_path = PUBMED_DIR / f"{gene}.json"
    uniprot_path = UNIPROT_DIR / f"{gene}.json"

    if not pubmed_path.exists() or not uniprot_path.exists():
        print(f"[WARN] Missing data for {gene}")
        return

    pubmed_data = load_json(pubmed_path)
    uniprot_data = load_json(uniprot_path)

    llm_text = build_llm_text(gene, uniprot_data, pubmed_data)

    merged = {
        "gene": gene,
        "metadata": {
            "gene": gene,
            "uniprot_id": uniprot_data.get("uniprot_id"),
            "protein_name": uniprot_data.get("protein_name"),
            "sources": ["UniProt", "PubMed"]
        },
        "llm_context": llm_text
    }

    output_file = OUTPUT_DIR / f"{gene}.json"
    with open(output_file, "w", encoding="utf-8") as f:
        json.dump(merged, f, indent=2)

    print(f"[OK] Created {output_file.name}")


if __name__ == "__main__":
    print("Starting gene data merge...")
    for gene in GENES:
        merge_gene(gene)
    print("Merge completed.")
