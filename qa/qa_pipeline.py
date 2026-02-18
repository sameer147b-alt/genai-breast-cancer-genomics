import argparse
import json
from pathlib import Path

import faiss
import numpy as np
from sentence_transformers import SentenceTransformer


BASE_DIR = Path(__file__).resolve().parents[1]
VECTORSTORE_DIR = BASE_DIR / "vectorstore" / "faiss_index"

INDEX_PATH = VECTORSTORE_DIR / "index.faiss"
DOCUMENTS_PATH = VECTORSTORE_DIR / "documents.json"
MANIFEST_PATH = VECTORSTORE_DIR / "manifest.json"

DEFAULT_MODEL = "sentence-transformers/all-MiniLM-L6-v2"


def parse_args() -> argparse.Namespace:
    parser = argparse.ArgumentParser(description="Offline QA over breast cancer genomics docs.")
    parser.add_argument("--query", type=str, help="Run a single query and exit.")
    parser.add_argument("--top-k", type=int, default=2)
    parser.add_argument("--model-name", default=None)
    return parser.parse_args()


def load_assets(model_name_override: str | None):
    if not INDEX_PATH.exists() or not DOCUMENTS_PATH.exists():
        raise FileNotFoundError(
            "Vectorstore missing. Run processing/merge_gene_data.py and "
            "embeddings/build_embeddings.py first."
        )

    manifest = {}
    if MANIFEST_PATH.exists():
        manifest = json.loads(MANIFEST_PATH.read_text(encoding="utf-8"))

    model_name = model_name_override or manifest.get("model_name", DEFAULT_MODEL)
    model = SentenceTransformer(model_name)

    index = faiss.read_index(str(INDEX_PATH))
    documents = json.loads(DOCUMENTS_PATH.read_text(encoding="utf-8"))

    if index.ntotal != len(documents):
        raise ValueError(
            f"Index/document mismatch: index={index.ntotal}, docs={len(documents)}"
        )

    return model, index, documents


def retrieve(query: str, model, index, documents: list[dict], top_k: int):
    top_k = max(1, min(top_k, len(documents)))
    query_embedding = model.encode(
        [query], convert_to_numpy=True, normalize_embeddings=True
    )
    query_embedding = np.asarray(query_embedding, dtype=np.float32)

    scores, indices = index.search(query_embedding, top_k)
    results = []
    for score, idx in zip(scores[0], indices[0]):
        if idx < 0:
            continue
        doc = documents[idx]
        results.append({"score": float(score), "doc": doc})
    return results


def format_result(result: dict) -> str:
    doc = result["doc"]
    metadata = doc.get("metadata", {})
    gene = metadata.get("gene", "N/A")
    protein = metadata.get("protein_name", "N/A")
    sources = metadata.get("sources", [])
    sources_text = ", ".join(sources) if isinstance(sources, list) else str(sources)
    snippet = doc.get("page_content", "").replace("\n", " ").strip()[:500]

    return (
        f"Gene: {gene} | Protein: {protein} | Sources: {sources_text} | "
        f"Similarity: {result['score']:.3f}\n{snippet}..."
    )


def run_query(query: str, top_k: int, model, index, documents) -> None:
    results = retrieve(query, model, index, documents, top_k)
    print("\nTop matches:")
    for rank, result in enumerate(results, start=1):
        print(f"{rank}. {format_result(result)}")
    print()


def main() -> None:
    args = parse_args()
    model, index, documents = load_assets(args.model_name)

    if args.query:
        run_query(args.query, args.top_k, model, index, documents)
        return

    print("\nBreast Cancer Genomics QA (vectorstore)")
    print("Type 'exit' to quit\n")

    while True:
        query = input("Ask a question: ").strip()
        if query.lower() == "exit":
            break
        run_query(query, args.top_k, model, index, documents)


if __name__ == "__main__":
    main()
