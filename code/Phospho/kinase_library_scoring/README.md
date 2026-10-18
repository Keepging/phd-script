# Kinase Library scoring of phosphosites

Scores every site in `input/site_meta.tsv` against all Kinase Library PSSMs
(Johnson et al. 2023 S/T atlas; Yaron-Barir et al. 2024 Y atlas) with the PyPI
package `kinase-library` 1.8.0. No enrichment, no statistics: scores and percentiles only.

## Steps

1. `fetch_uniprot.py` – UniProt REST (`/uniprotkb/search`, TSV) canonical sequences for all
   3,768 accessions in `site_meta.tsv`. Writes `uniprot_sequences.tsv` (entry/sequence version,
   modification dates, sequence) and `uniprot_release.json` (`x-uniprot-release`,
   `x-uniprot-release-date` response headers, fetch time, secondary→primary accession map).
2. `score_sites.py`
   * parses `site_id` = `ACC_[STY]pos`, checks the residue at `pos` (1-based) in the canonical
     sequence. `residue_check`: `match` / `mismatch_STY` (another phosphoacceptor at that position,
     still scored, with the UniProt residue) / `mismatch_non_STY` / `position_out_of_range` /
     `accession_not_found` (the last three are not scored).
   * `window_-5_+4`: residues −5…+4 around the site, `_` outside the protein (positions used by
     the S/T matrices). `substrate_15mer`: −7…+7 window in kinase-library's required format
     (odd length, phosphoacceptor at the centre, lower-case); this is what is fed to the package
     because tyrosine matrices also read position +5.
   * S/T sites → 311 ser_thr matrices; Y sites → 93 tyrosine matrices (78 canonical Y kinases +
     15 dual-specificity `*_TYR` matrices, `non_canonical=True`).
   * package defaults otherwise: `st_fav=True`, `pp=False` (no phospho-priming), percentiles
     relative to the bundled Ochoa et al. scored phosphoproteome. `score` is log2.

## Outputs (`output/`)

* `kinase_scores_long.parquet` – `site_id, kinase, kinase_family, score, percentile`
  (`kinase_family` from the package's `kinome_information.tsv` FAMILY column).
* `site_sequences.tsv` – `site_id, protein, position, residue_expected, residue_in_uniprot,
  residue_check, window_-5_+4, substrate_15mer, uniprot_entry, sequence_version, sequence_length`.
* `uniprot_sequences.tsv`, `uniprot_release.json` – provenance.

## Run

```bash
bash run.sh            # creates .venv with Python 3.12 via uv, fetches, scores -> output/
```
