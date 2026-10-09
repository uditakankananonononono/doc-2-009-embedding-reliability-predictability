# DOC-2-009 R3 results (frozen analysis, PROTOCOL.md lock-1 = fdfbfc3; no post-hoc changes)

Label (set mechanically by analysis.py): HONEST NEGATIVE. G1, G2 and G3 all fail. Not evidence of absence: the external CI is wide, see below. Independent gate review is pending; no claim is made until it clears.

## Results (ESM-2 650M wt-marginal rho for 204 assays; ridge lambda 1)
- Sets: discovery n=98 assays (year <= 2021). External n=103 assays (year >= 2022, after removing the 2 overlap proteins DYR_ECOLI and SRC_HUMAN, 3 external assays).
- G1 (discovery leave-one-protein-out R2): B0 (intercept) 0.000, B1 (taxon + Tsuboyama) 0.155, M (B1 + log length, WTLL, annotation depth, n_pubmed) 0.155. Delta M-B1 = -0.0002, needed >= +0.05. Fail.
- G2 (external, fit on discovery once): R2(M) = +0.021, R2(B1) = -0.002, delta = +0.023, bootstrap 95% CI [-0.147, +0.184], needed delta >= +0.05 with CI lower bound > 0. Fail. The CI is wide and contains +0.05, so a modest gain is not excluded.
- G3 (permutation of the M-only features across proteins, 2000 permutations): p = 0.253, needed < 0.05. Fail.
- Descriptive, not a gate: in discovery, taxon + assay class alone explain about 15% of out-of-fold variance in per-assay rho (B1 R2 = 0.155), but on the external cohort that baseline had R2 -0.002, so the baseline did not transfer either.

## Blind spots and limits (kept as written in the protocol)
- Length is not blind: the log(seq_len) vs rho relation was already examined in DOC-2-077 (raw Spearman -0.232 over 204 assays; HONEST NEGATIVE there). WTLL and annotation exposure were used in DOC-2-073.
- The 2 discovery/external overlap proteins were removed from external only, as pre-stated.
- Protein-level groups only (no family clustering); wild-type-marginal scoring (a different estimand from masked-marginal); CPU fp32.
- Not run: MSA-depth features (not acquired), second model, any causal claim about why reliability differs.
- n is about 200 assays, so the test is low power; the label is a gate outcome, not a finding that reliability cannot be predicted.

## Receipts
- run_log.txt: command, UTC start and end (2026-10-09T14:04:43Z to 14:06:37Z), library versions (python 3.10.12, numpy 2.2.6, pandas 2.3.3), input md5 check, the RESULT_JSON line, exit code 0.
- Input check in run_log.txt shows per_assay.csv and exposures.tsv OK. The check listed exclusions.tsv as FAILED only because the file was not copied into the run directory. analysis.py does not read exclusions.tsv. input_check.txt is the check re-run afterwards with all three files present.
- analysis.py was run once on real data after the smoke test on synthetic random rho. These run statements are builder statements; the lock-before-run order is shown by the lock-1 tag at fdfbfc3 and the repo commit history.
