# DOC-2-009 R3: "The Reliability Map of Protein Embeddings" - can per-assay ESM-2 reliability be predicted from cheap features? (frozen protocol, lock-1)

Written 2026-10-09 IST and committed BEFORE any model of per-assay rho on features has been fit or opened.

Question: per-assay Spearman rho of ESM-2 650M zero-shot variant scoring varies a lot across assays. Do cheap, pre-scoring features (sequence length, wild-type log-likelihood, taxon, assay class, protein annotation depth) predict which assays the model is reliable on, beyond a trivial taxon + assay-class baseline, on held-out proteins and on an external-year cohort?

## Data (frozen, from DOC-2-073, same run)
- per_assay.csv (204 ProteinGym v1 assays; rho, wtll, seq_len, n_scored), exposures.tsv (UniProt_ID, taxon, year, annotation_score, n_pubmed). md5 in INPUT_HASHES.md5.
- Disclosure: the rho values exist from DOC-2-073 and were examined there only against annotation exposure, and in DOC-2-077 against sequence length (HONEST NEGATIVE; raw Spearman log(seq_len) vs rho = -0.232). This protocol does not hide that: seq_len was already looked at.
- Drop assays with n_scored < 50. Discovery = year <= 2021; External = year >= 2022 (as in 073/077). 2 proteins overlap the two sets (DYR_ECOLI, SRC_HUMAN); they are kept in discovery and removed from external for this analysis.

## Models (all fixed in advance)
- y = rho. Ridge regression (lambda = 1, standardized features), pure numpy.
- B0: intercept only. B1: taxon dummies + Tsuboyama indicator (DMS_id contains "Tsuboyama_2023").
- M: B1 + log(seq_len) + wtll + annotation_score + log(1+n_pubmed) (NaN annotation -> 0 with a missing flag).
- Evaluation 1 (discovery): leave-one-protein-out (group = UniProt_ID) out-of-fold R2, reference = training-fold mean.
- Evaluation 2 (external): fit on all discovery, predict external once; R2 reference = discovery mean.
- CI: protein-cluster bootstrap of the external delta R2 (M minus B1), 10,000 resamples, seed 12345, percentile 95%. Permutation: shuffle the M-only features across proteins (2,000 permutations, seed 12345) for the delta.

## Gates
- G1: discovery LOPO R2(M) - R2(B1) >= +0.05.
- G2: external R2(M) > 0 and delta R2 (M - B1) >= +0.05 with bootstrap CI lower bound > 0.
- G3: external permutation p < 0.05 for the delta.
Labels: POSITIVE if G1, G2, G3 all pass; HONEST NEGATIVE otherwise. No other label. A "partly" result is reported as the table of numbers with the label HONEST NEGATIVE.

## Rules
No simulated data, no stubs. Smoke test only on synthetic rho. analysis.py is run once. Protein-level clusters only; wild-type-marginal scoring; CPU fp32. Deviation budget: one tolerance line; any amendment is a new dated AMENDMENT-N.md committed before outcomes exist.
Out of scope / not run: any claim about why reliability differs; second model; MSA-depth features (not acquired).
