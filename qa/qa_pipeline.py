from pathlib import Path

from langchain_community.vectorstores import FAISS
from langchain_huggingface import HuggingFaceEmbeddings


BASE_DIR = Path(__file__).resolve().parents[1]
VECTORSTORE_DIR = BASE_DIR / "vectorstore" / "faiss_index"


def load_vectorstore():
    index_file = VECTORSTORE_DIR / "index.faiss"
    store_file = VECTORSTORE_DIR / "index.pkl"

    if not index_file.exists() or not store_file.exists():
        raise FileNotFoundError(
            "Vectorstore not found. Run processing/merge_gene_data.py and "
            "embeddings/build_embeddings.py first."
        )

    embeddings = HuggingFaceEmbeddings(
        model_name="sentence-transformers/all-MiniLM-L6-v2"
    )
    return FAISS.load_local(
        str(VECTORSTORE_DIR),
        embeddings,
        allow_dangerous_deserialization=True,
    )


def retrieve_answer(query, vectorstore, top_k=2):
    return vectorstore.similarity_search(query, k=top_k)


def format_answer(doc):
    metadata = doc.metadata or {}
    gene = metadata.get("gene", "N/A")
    protein = metadata.get("protein_name", "N/A")
    source = metadata.get("sources", [])
    source_text = ", ".join(source) if isinstance(source, list) else str(source)
    snippet = doc.page_content[:500].replace("\n", " ").strip()

    return (
        f"Gene: {gene} | Protein: {protein} | Sources: {source_text}\n"
        f"{snippet}..."
    )


def main():
    print("\nBreast Cancer Genomics QA (vectorstore)")
    print("Type 'exit' to quit\n")

    vectorstore = load_vectorstore()

    while True:
        query = input("Ask a question: ").strip()
        if query.lower() == "exit":
            break

        docs = retrieve_answer(query, vectorstore)
        print("\nTop matches:")
        for idx, doc in enumerate(docs, start=1):
            print(f"{idx}. {format_answer(doc)}")
        print()


if __name__ == "__main__":
    main()
