# SIFT: Search Conjoint

## Measuring Preferences and Search Costs in a Single Survey Instrument

Dan Yavorsky (GBK Collective; Rady School of Management, UC San Diego) and Eric T. Bradlow (The Wharton School).

A survey instrument in which the respondent faces a list-screen of alternatives,
clicks through to inspect individual listings, and then chooses. The resulting
click-and-choose data are fit with a sequential search model, recovering both
preference parameters — as conjoint does — and search costs, for which conjoint
has no analogue. Conjoint is the zero-search-cost special case.

## Status

Skeleton manuscript. The introduction is intentionally blank; the literature
review holds annotated bullets pending a citation-verification pass; the
remaining sections are section and subsection headers with scoping notes in
HTML comments.

**All references in `references.bib` and `sections/litreview.qmd` are
unverified** — reconstructed from memory, flagged [H]/[M]/[L] for confidence.
Verify before circulating.

## Build

```
quarto render
```

Outputs to `_manuscript/`.

## Layout

- `index.qmd` — front matter, abstract, section includes
- `sections/` — one file per section
- `latex_commands.qmd` — shared macros
- `references.bib` — bibliography (unverified)
- `lit/` — PDFs of cited work
- `notebooks/` — computational notebooks for the manuscript project type
