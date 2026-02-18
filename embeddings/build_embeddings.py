import argparse
import hashlib
import json
from pathlib import Path

import faiss


BASE_DIR = Path(__file__).resolve().parents[1]
PROCESSED_DIR = BASE_DIR / "data" / "processed"
VECTORSTORE_DIR = BASE_DIR / "vectorstore" / "faiss_index"

INDEX_PATH = VECTORSTORE_DIR / "index.faiss"
DOCUMENTS_PATH = VECTORSTORE_DIR / "documents.json"
MANIFEST_PATH = VECTORSTORE_DIR / "manifest.json"
LEGACY_PATH = VECTORSTORE_DIR / "index.pkl"

DEFAULT_MODEL = "sentence-transformers/all-MiniLM-L6-v2"


def parse_args() -> argparse.Namespace:
    parser = argparse.ArgumentParser(description="Build FAISS index from processed data.")
    parser.add_argument("--model-name", default=DEFAULT_MODEL)
    parser.add_argument("--batch-size", type=int, default=32)
    parser.add_argument("--force", action="store_true")
    return parser.parse_args()


def load_documents() -> tuple[list[dict], list[Path]]:
    files = sorted(PROCESSED_DIR.glob("*.json"))
    if not files:
        raise ValueError(f"No processed JSON files found in {PROCESSED_DIR}")

    documents = []
    for file_path in files:
        data = json.loads(file_path.read_text(encoding="utf-8"))
        documents.append(
            {
                "source_file": file_path.name,
                "page_content": data["llm_context"],
                "metadata": data.get("metadata", {}),
            }
        )
    return documents, files


def compute_source_hash(files: list[Path]) -> str:
    hasher = hashlib.sha256()
    for file_path in files:
        hasher.update(file_path.name.encode("utf-8"))
        hasher.update(file_path.read_bytes())
    return hasher.hexdigest()


def load_manifest() -> dict:
    if not MANIFEST_PATH.exists():
        return {}
    return json.loads(MANIFEST_PATH.read_text(encoding="utf-8"))


def is_up_to_date(manifest: dict, source_hash: str, model_name: str, doc_count: int) -> bool:
    if not manifest:
        return False
    return (
        manifest.get("source_hash") == source_hash
        and manifest.get("model_name") == model_name
        and manifest.get("doc_count") == doc_count
        and INDEX_PATH.exists()
        and DOCUMENTS_PATH.exists()
    )


def save_manifest(source_hash: str, model_name: str, doc_count: int) -> None:
    manifest = {
        "source_hash": source_hash,
        "model_name": model_name,
        "doc_count": doc_count,
    }
    MANIFEST_PATH.write_text(json.dumps(manifest, indent=2), encoding="utf-8")


def build_index(documents: list[dict], model_name: str, batch_size: int) -> faiss.Index:
    import numpy as np
    from sentence_transformers import SentenceTransformer

    model = SentenceTransformer(model_name)
    texts = [doc["page_content"] for doc in documents]

    embeddings = model.encode(
        texts,
        batch_size=batch_size,
        show_progress_bar=True,
        convert_to_numpy=True,
        normalize_embeddings=True,
    )
    embeddings = np.asarray(embeddings, dtype=np.float32)

    index = faiss.IndexFlatIP(embeddings.shape[1])
    index.add(embeddings)
    return index


def main() -> None:
    args = parse_args()
    VECTORSTORE_DIR.mkdir(parents=True, exist_ok=True)

    print("Embedding build started")
    print(f"Processed dir: {PROCESSED_DIR}")
    print(f"Vectorstore dir: {VECTORSTORE_DIR}")

    documents, files = load_documents()
    source_hash = compute_source_hash(files)
    manifest = load_manifest()

    if not args.force and is_up_to_date(
        manifest, source_hash, args.model_name, len(documents)
    ):
        print("Embeddings are up to date. Skipping rebuild.")
        return

    print(f"Building index for {len(documents)} documents")
    index = build_index(documents, model_name=args.model_name, batch_size=args.batch_size)

    faiss.write_index(index, str(INDEX_PATH))
    DOCUMENTS_PATH.write_text(json.dumps(documents, indent=2), encoding="utf-8")
    save_manifest(source_hash, args.model_name, len(documents))

    if LEGACY_PATH.exists():
        LEGACY_PATH.unlink()

    print("Vector index built successfully")
    print(f"Index: {INDEX_PATH}")
    print(f"Documents: {DOCUMENTS_PATH}")
    print(f"Manifest: {MANIFEST_PATH}")


if __name__ == "__main__":
    main()
