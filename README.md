# SIFT: Search Conjoint

## Measuring Preferences and Search Costs in a Single Survey Instrument

Dan Yavorsky (GBK Collective; Rady School of Management, UC San Diego) and Eric T. Bradlow (The Wharton School).

A survey instrument in which the respondent faces a list-screen of alternatives,
clicks through to inspect individual listings, and then chooses. The resulting
click-and-choose data are fit with a sequential search model, recovering both
preference parameters — as conjoint does — and search costs, for which conjoint
has no analogue. Conjoint is the zero-search-cost special case.

## The name

**SIFT — Search Intensity and Feature Tradeoffs.**

The expansion names both deliverables, which is the whole pitch: preferences
*and* search costs from one instrument. "Feature tradeoffs" deliberately borrows
conjoint's own vocabulary, which makes the nesting argument quietly.

The bare word carries the pitch — it is one syllable, it is verbable ("we ran a
SIFT", "the SIFT data"), and it describes what the respondent actually does. The
expansion is meant to appear once, on first use, and then be retired. The test
applied when choosing it: TURF works as an acronym precisely because nobody ever
expands it. So the criterion is whether an expansion can survive a single
footnote and never be spoken again — which argues for the version most
defensible in front of a referee, not the one most memorable in a pitch.

Alternates considered and set aside:

| candidate | why not |
|---|---|
| Sequential Inspection and Feature Tradeoffs | leads with mechanism, truer to the Weitzman protocol, less client-friendly |
| Shopping Information Foraging Task | nods to information foraging theory (Pirolli & Card), but "Task" names an activity rather than a method |
| Search-Inferred Feature Tradeoffs | cleanest grammar, but hides the search cost, which is the actual differentiator |
| Search in Field Testing | **rejected** — in this literature "field" means a real market experiment with real money, so this promises the opposite of what a survey instrument delivers |

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
