# JOSS paper

The JOSS paper for Jexpresso.

| file | contents |
|---|---|
| `paper.md` | the paper, JOSS Markdown source (the submission path) |
| `paper.tex` | the same paper in LaTeX |
| `headers.tex` | preamble — **local build only** (see below) |
| `paper.bib` | BibTeX for the references, shared by both |

## Figure to supply

`results_collage.png` is **not in the repository** — it is a placeholder for a
composite figure of representative results (atmospheric flows with cloud
microphysics, mountain and gravity waves, turbulent boundary layers, high-speed
compressible flow, MHD, radiative transfer). Drop the file into this directory
and both sources pick it up.

Until it exists, `paper.tex` draws a framed grey box in its place via
`\IfFileExists`, so the document still compiles; once the real image is in,
that guard can be reduced to the bare `\includegraphics` line. `paper.md`
references the file directly, so **pandoc will fail to build a PDF until the
image is present** — that is deliberate, so a missing figure cannot be
submitted unnoticed.

The JOSS logo is handled the same way: `joss-logo.png` (or `logo.png`) is
expected in this directory and is gitignored rather than committed.

## Build

```bash
cd paper && pdflatex paper.tex
```

## Length

The binding JOSS constraint is the **word count, 1750 maximum** — not the page
count. Current state, counting prose only (no code listings, display math or
references):

| source | words |
|---|---|
| `paper.md` | ~1230 |
| `paper.tex` | ~1220 |

Check it after any edit with:

```bash
python3 - <<'EOF'
import re
s = open("paper.md").read().split('---', 2)[2]
s = re.sub(r"```.*?```|\$\$.*?\$\$", " ", s, flags=re.S)
s = re.sub(r"\$[^$]*\$|\[@[^\]]*\]|[#*`\[\]]", " ", s).replace("# References", "")
print(len([w for w in s.split() if re.search(r"[A-Za-z]", w)]), "words / 1750")
EOF
```

The local `pdflatex` build currently runs to 5 pages. That number is **not**
meaningful for the submission: the Open Journals bot compiles `paper.md` with
its own template, so its margins, fonts and bibliography sizing — not the ones
in `headers.tex` — decide the real page count. If you nonetheless want the local
build shorter, the cheapest cuts are the Key Features list (it partly restates
the narrative) and the Euler display equation (the MHD listing already carries
the "scaling up" argument).

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
