import json
from pathlib import Path

from langchain_community.vectorstores import FAISS
from langchain_core.documents import Document
from langchain_huggingface import HuggingFaceEmbeddings


# -----------------------
# PATHS
# -----------------------
BASE_DIR = Path(__file__).resolve().parent.parent
PROCESSED_DIR = BASE_DIR / "data" / "processed"
VECTORSTORE_DIR = BASE_DIR / "vectorstore" / "faiss_index"

VECTORSTORE_DIR.mkdir(parents=True, exist_ok=True)

print("Embedding script started")
print(f"BASE_DIR: {BASE_DIR}")
print(f"PROCESSED_DIR EXISTS: {PROCESSED_DIR.exists()}")

# -----------------------
# LOAD DOCUMENTS
# -----------------------
documents = []

for file in sorted(PROCESSED_DIR.glob("*.json")):
    print(f"Loading {file.name}")
    data = json.loads(file.read_text(encoding="utf-8"))
    documents.append(
        Document(
            page_content=data["llm_context"],
            metadata=data["metadata"],
        )
    )

if not documents:
    raise ValueError(f"No processed JSON files found in: {PROCESSED_DIR}")

print(f"Loaded {len(documents)} documents")

# -----------------------
# EMBEDDINGS
# -----------------------
embeddings = HuggingFaceEmbeddings(
    model_name="sentence-transformers/all-MiniLM-L6-v2"
)

# -----------------------
# VECTOR STORE
# -----------------------
vectorstore = FAISS.from_documents(documents, embeddings)
vectorstore.save_local(str(VECTORSTORE_DIR))

print("Vectorstore successfully created")
print(f"Saved to: {VECTORSTORE_DIR}")
