## Scope Definition

This project focuses on **breast cancer genomics**, with an emphasis on
interpreting the biological role of well-characterized cancer-associated genes
using GenAI-based retrieval and reasoning.

### Initial Gene Set
BRCA1, BRCA2, TP53, PIK3CA, ERBB2, PTEN, CDH1, GATA3, MAP3K1, AKT1,
ESR1, FOXA1, NF1, RB1, CCND1, FGFR1, MYC, ARID1A, TBX3, KMT2C

### Use Case
Research and interpretation support — **not diagnosis or prediction**.

# GenAI-Powered Breast Cancer Genomics Research Assistant

## Overview
This project implements a GenAI-driven research support system for interpreting
breast cancer–associated genes by integrating cancer genomics data and biomedical
literature using a Retrieval-Augmented Generation (RAG) architecture.

The system is designed to assist researchers by providing evidence-grounded,
citation-backed explanations of gene function and cancer relevance.

---

## Problem Statement
Interpreting the role of cancer-associated genes requires synthesizing information
from genomics datasets, functional annotations, and a rapidly growing body of
biomedical literature. This process is time-consuming and prone to fragmented
understanding.

---

## Why GenAI
Large Language Models combined with retrieval-based pipelines enable:
- Context-aware synthesis of biological knowledge
- Evidence-grounded reasoning with citations
- Scalable exploration of biomedical literature

However, LLMs are constrained by hallucination risks, which this system explicitly
addresses through controlled retrieval and evaluation.

---

## Data Sources
- TCGA (Breast Cancer mutation summaries)
- UniProt (Gene and protein annotations)
- PubMed (Breast cancer–specific abstracts)

---

## System Architecture
[Architecture diagram to be added]

---

## Example Queries
- What is the role of BRCA1 mutations in breast cancer?
- How does PIK3CA contribute to tumor progression?

---

## Evaluation
The system is evaluated using manual relevance assessment, citation accuracy,
and hallucination checks on curated gene-level queries.

---

## Limitations & Ethics
- Not intended for clinical diagnosis or treatment decisions
- Dependent on available literature and curated datasets
- LLM-generated interpretations may omit emerging evidence

---

## Future Work
- Variant-level interpretation
- Pathway integration
- Biomedical knowledge graph expansion
