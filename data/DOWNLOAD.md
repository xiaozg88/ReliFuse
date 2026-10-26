# Large artifacts not stored in this repository

Two inputs are too large for plain Git hosting. Both are needed only if you want
to re-run the structure path (`data_list.pt`) or the sequence path
(`esm_meanpool_*.npz`) from scratch; every figure and table in the manuscript can
be regenerated from the files shipped in `data/` alone.

## 1. AlphaFold contact-graph dataset — `data_list.pt` (≈1.05 GB)

PyTorch Geometric dataset of 8,069 residue contact graphs
(nodes = residues with 32-dim features; edges = Cα–Cα < 8.0 Å).

* Download: `<Zenodo DOI — to be inserted on acceptance>`
* Place at: `data/data_list.pt`
* Verify: 8,069 graphs, `num_node_features = 32`, `num_classes = 2`

## 2. Protein-language-model mean-pool caches (≈90 MB total)

| File | Encoder | Dim |
|---|---|---|
| `esm_meanpool_8069_esm2_150m.npz` | ESM-2 150M | 640 |
| `esm_meanpool_8069_esm2_650m_faesm.npz` | ESM-2 650M | 1280 |
| `esm_meanpool_8069_esmc300m.npz` | ESM-C 300M | 1152 |

Each archive contains `names`, `y`, `splits`, `X` (mean-pooled embeddings for
8,067 proteins; two proteins have no cached embedding).

* Download: `<Zenodo DOI — to be inserted on acceptance>`
* Or regenerate: run the encoder over the UniProt FASTA of the 8,069 proteins,
  mean-pool the final-layer representations, and save under the same keys.

## 3. Suggested upload layout

When creating the Zenodo record, upload the two artifacts as separate files with
the exact names above and put the DOI in this file, in `README.md`, and in the benchmark release
data card. A `.zenodo.json` with the manuscript metadata is
recommended so the DOI resolves to the paper.
