# DOC-2-009 R3 results (frozen analysis, PROTOCOL.md lock-1 = fdfbfc3; no post-hoc changes)

Label (set mechanically by analysis.py): HONEST NEGATIVE. G1, G2 and G3 all fail. Not evidence of absence: the external CI is wide, see below. Independent gate review is pending; no claim is made until it clears.

## Results (ESM-2 650M wt-marginal rho for 204 assays; ridge lambda 1)
- Sets: discovery n=98 assays (year <= 2021). External n=103 assays (year >= 2022, after removing the 2 overlap proteins DYR_ECOLI and SRC_HUMAN, 3 external assays).
- G1 (discovery leave-one-protein-out R2): B0 (intercept) 0.000, B1 (taxon + Tsuboyama) 0.155, M (B1 + log length, WTLL, annotation depth, n_pubmed) 0.155. Delta M-B1 = -0.0002, needed >= +0.05. Fail.
- G2 (external, fit on discovery once): R2(M) = +0.021, R2(B1) = -0.002, delta = +0.023, bootstrap 95% CI [-0.147, +0.184], needed delta >= +0.05 with CI lower bound > 0. Fail. The bootstrap resamples external proteins only, with the discovery fit held fixed, so this CI omits discovery-training variance and is narrower than a full-uncertainty CI; "wide" understates the uncertainty. It contains 0 and +0.05, so a modest gain is not excluded. The +0.021 external R2 for M is a gate failure with a CI containing 0, not a positive finding.
- G3 (permutation of the M-only features across proteins, 2000 permutations): p = 0.253, needed < 0.05. Fail. The permutation shuffles the M-only features across discovery (training-side) proteins and measures the external delta, i.e. it tests the feature-to-rho link in training, as the protocol specifies.
- Descriptive, not a gate: in discovery, the baseline B1 has out-of-fold R2 = 0.155, but on the external cohort it had R2 -0.002, so the baseline did not transfer either. **Amended wording (see Amendment note below):** B1's discovery R2 came from taxon alone, not from taxon + assay class.

## Blind spots and limits (kept as written in the protocol)
- Length is not blind: the log(seq_len) vs rho relation was already examined in DOC-2-077 (raw Spearman -0.232 over 204 assays; HONEST NEGATIVE there). WTLL and annotation exposure were used in DOC-2-073.
- The 2 discovery/external overlap proteins were removed from external only, as pre-stated.
- Protein-level groups only (no family clustering); wild-type-marginal scoring (a different estimand from masked-marginal); CPU fp32.
- Not run: MSA-depth features (not acquired), second model, any causal claim about why reliability differs.
- n is about 200 assays, so the test is low power; the label is a gate outcome, not a finding that reliability cannot be predicted.

## Receipts
- run_log.txt: command, UTC start and end (2026-10-09T14:04:43Z to 14:06:37Z), library versions (python 3.10.12, numpy 2.2.6, pandas 2.3.3), input md5 check, the RESULT_JSON line, exit code 0.
- Input check in run_log.txt shows per_assay.csv and exposures.tsv OK. The check listed exclusions.tsv as FAILED only because the file was not copied into the run directory. analysis.py does not read exclusions.tsv. input_check.txt is a post-run re-check (committed after the run): it shows the files match now, not at run time.
- "analysis.py was run once on real data after a smoke test on synthetic random rho" is a builder statement, not repo evidence (the synthetic test left no record). The lock-before-run order is shown by the lock-1 tag at fdfbfc3 and the repo commit history; the independent gate re-ran the analysis bit-identically.

## Amendment note (wording only, 2026-10-09 IST, after review)
Triggered by the DOC-2-009 R4 finding (source: https://github.com/uditakankananonononono/doc-2-009-r4-assay-descriptors-reliability, main HEAD fc92b6f622fc9f37e6ae6f4fb7194593f9a33f14, RESULTS.md). All 63 Tsuboyama 2023 assays sit in the external split (0 in discovery). The Tsuboyama indicator in B1 and M was therefore constant 0 in discovery and could not be learned. B1's discovery leave-one-protein-out R2 of 0.155 came from taxon alone. The earlier descriptions "taxon + Tsuboyama" (PROTOCOL.md B1 definition, and the "taxon + assay class explain about 15%" line) overstated what B1 learned. The label (HONEST NEGATIVE) and all gate numbers are unchanged. PROTOCOL.md is left as committed at lock-1; this note is the correction.
