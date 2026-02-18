# genai-breast-cancer-genomics

GenAI-powered, retrieval-first research assistant for breast cancer genomics.
It ingests curated gene annotations from UniProt and abstracts from PubMed,
builds a local FAISS index, and answers questions with evidence-grounded context.

## What This Project Does

- Fetches gene-level protein/function/disease annotations from UniProt.
- Fetches breast-cancer-focused PubMed abstracts for a fixed gene set.
- Merges both sources into normalized per-gene context documents.
- Builds a fast local vector index for semantic retrieval.
- Runs an offline QA loop over the indexed corpus.

## Current Gene Set

`BRCA1`, `BRCA2`, `TP53`, `PIK3CA`, `ERBB2`

## Architecture

1. `ingestion/uniprot_loader.py`
2. `ingestion/pubmed_loader.py`
3. `processing/merge_gene_data.py`
4. `embeddings/build_embeddings.py`
5. `qa/qa_pipeline.py`

Data artifacts:

- `data/uniprot/*.json`
- `data/pubmed/*.json`
- `data/processed/*.json`
- `vectorstore/faiss_index/index.faiss`
- `vectorstore/faiss_index/documents.json`
- `vectorstore/faiss_index/manifest.json`

## Performance and Deployment Optimizations

- Batched PubMed `efetch` requests (fewer round trips, faster ingestion).
- Connection-reuse + retry strategy for UniProt requests.
- Lightweight runtime: pure `faiss + sentence-transformers` retrieval path.
- Incremental embedding rebuilds via content-hash manifest checks.
- CLI flags for query mode, top-k control, model override, and force rebuilds.

## Quickstart

### 1) Create environment

```bash
python -m venv .venv
.venv\Scripts\activate
pip install -r requirements.txt
```

### 2) (Recommended) Set NCBI identity

```bash
set NCBI_EMAIL=you@example.com
set NCBI_API_KEY=your_ncbi_api_key
```

`NCBI_API_KEY` is optional, but increases allowed request throughput.

### 3) Run full pipeline

```bash
python ingestion/uniprot_loader.py
python ingestion/pubmed_loader.py
python processing/merge_gene_data.py
python embeddings/build_embeddings.py
python qa/qa_pipeline.py
```

### 4) One-shot QA query

```bash
python qa/qa_pipeline.py --query "What is the role of BRCA1 in breast cancer?" --top-k 3
```

## Useful CLI Options

### PubMed ingestion

```bash
python ingestion/pubmed_loader.py --genes BRCA1 TP53 --retmax 30 --batch-size 10
```

Env controls:

- `PUBMED_RETMAX`
- `PUBMED_BATCH_SIZE`
- `PUBMED_MAX_RETRIES`
- `PUBMED_REQUEST_DELAY`

### Embedding/index build

```bash
python embeddings/build_embeddings.py --batch-size 32
python embeddings/build_embeddings.py --force
```

### QA

```bash
python qa/qa_pipeline.py --top-k 5
python qa/qa_pipeline.py --query "How does HER2 affect breast cancer?"
```

## Deployment Notes

- First run will download the embedding model; subsequent runs use local cache.
- Keep `vectorstore/faiss_index/*` as a build artifact for fast startup in demos.
- For API/serving integration, wrap `qa/qa_pipeline.py` retrieval functions in a web layer (FastAPI/Flask) and preload model + FAISS index once at process start.

## Limitations

- Research support only. Not for diagnosis, prognosis, or treatment decisions.
- Retrieval quality depends on source completeness and selected genes.
- Output reflects indexed data snapshot, not real-time literature updates.

## License / Usage

Add your preferred license before public release if this repository is intended for external use.
