# JOSS paper

A JOSS-length (≤ 3 pages) paper for Jexpresso.

| file | contents |
|---|---|
| `paper.tex` | the paper |
| `headers.tex` | preamble — **local build only** (see below) |
| `paper.bib` | BibTeX for the same references, for the Markdown workflow |

## Build

```bash
cd paper && pdflatex paper.tex
```

Produces `paper.pdf` at **3 pages**, which is the length this text is tuned to.
Two things control that and are easy to disturb:

- The bibliography `\itemsep` in `headers.tex` (currently `0.18em`). The
  reference list is ~13 entries, so a change here moves the last page boundary.
- The `geometry` margins. The JOSS sidebar lives in a wide right margin
  (`right=2.05in`, `marginparwidth=1.5in`); widening the text block reflows
  everything.

If you add a paragraph, expect to re-check the page count rather than assume it.

## Note on `headers.tex`

When a paper is submitted through the normal JOSS pipeline, the Open Journals
bot compiles `paper.md` + `paper.bib` with **its own** template and the
`headers.tex` here is ignored. This file exists only so the `.tex` compiles
locally for length checking; it is a reconstruction of the environments the JOSS
template provides (`CSLReferences`, `\ExternalLink`, `linky`, the listing style),
not the official template.

If you prefer the standard workflow, port `paper.tex` to `paper.md` with YAML
front matter for the author/affiliation block and cite with `[@key]` against
`paper.bib`.
