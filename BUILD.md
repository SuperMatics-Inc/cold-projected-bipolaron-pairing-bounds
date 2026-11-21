# BUILD: how the released PDF was built, and what CI checks

The manuscript is `paper/main.tex` with `paper/refs.bib`; the released PDF is
`paper/main.pdf`. The PDF is released byte for byte as reviewed: this repository
never rebuilds it on the way to a release, and a rebuild produces a different
file whenever the TeX distribution differs.

## 1. Build recipe

```
cd paper
SOURCE_DATE_EPOCH=1795132800 tectonic -X compile --keep-logs main.tex
```

- Tool: Tectonic 0.17.0. Packages resolve from the Tectonic bundle cache; the
  document needs `article`, `amsmath`, `amsthm`, `mathtools`, `geometry`,
  `hyperref`, `microtype`, `array` and `etoolbox`.
- Bibliography: `refs.bib` with `\bibliographystyle{plain}`, processed by
  Tectonic's internal BibTeX run.
- Result of this recipe: exit code 0 and a 20-page PDF, with no overfull box, no
  underfull box, no undefined reference or citation and no other warning line in
  the log.
- Determinism: the document is reproducible from the sources with the fixed
  `SOURCE_DATE_EPOCH`; the compressed byte stream can differ between TeX
  distributions, so the released hash below is the hash of the reviewed build.

## 2. Identifier policy

`paper/main.tex` leaves `\SELFDOI` empty (line 13) and prints the fallback phrase
"the record reserved for it" where the record's own DOI belongs. The DOI is
recorded in metadata only: in `README.md`, in `CITATION.cff` and in the deposit
metadata. The reason is that the reviewed PDF must be released unchanged, so the
DOI is not filled into the source and the PDF is not rebuilt. The companion
preprint is cited through `\STRONGDOI` (`10.5281/zenodo.24414616`).

## 3. Replay

```
pip install sympy mpmath numpy scipy
python replay/replay.py                 # exit 0
python replay/replay.py --all-mutants   # exit 0
python replay/replay.py --numerical     # exit 0
```

Observed final lines:

- default: `CERTIFIED REPLAY PASSED (5 exact-rational & outward-interval check groups).`
- `--all-mutants`: `ALL 16 GUARDED MUTANTS REJECTED (each exited 1).`
- `--numerical`: the two-site Table A (thermal crossover at `t = 0.5`, `U = 20`,
  `Omega = 1`) and Table B (six parameter points of its own), after the default
  checks.

`replay/README.md` is the guide pinned since review, and `replay/SHA256SUMS`
pins the program and that guide.

## 4. Policy checks on the released source

| Check | Command | Result |
|---|---|---|
| exactly one peer-review sentence | `grep -c 'has not been peer reviewed' paper/main.tex` | `1` |
| no em dash or en dash characters | `grep -cP '\x{2014}\|\x{2013}\|\x{2015}' paper/main.tex` | `0` |
| no lane, agent, model or process names | `grep -ni 'stage\|astra\|lane\|claude\|gpt\|anthropic\|opus\|sonnet\|agent\|zeta' paper/main.tex` | no match |
| no draft or review banners | `grep -ni 'draft' paper/main.tex` | one match, inside the AI disclosure ("were drafted and edited with AI assistance"), no banner |
| no review tags | `grep -ni 'review' paper/main.tex` | three matches, all inside the AI disclosure and the peer-review sentence ("Independent proof reviews", "direction and review", "has not been peer reviewed") |
| byline | `\author{Davis Rens\\SuperMatics Inc.\\ORCID 0000-0001-5839-5197}` | present once |
| own DOI placeholder | `\newcommand{\SELFDOI}{}` with `\record{\SELFDOI}{the record reserved for it}` | empty on purpose, see Section 2 |
| AI disclosure | `\section*{AI assistance}` | present |

## 5. Released files and hashes (SHA-256)

| File | Bytes | SHA-256 |
|---|---|---|
| `paper/main.pdf` | 225252 | `b92ee523c4947cce05e9f1b7906cf6abe47685568982759c4ec185c98df7633c` |
| `paper/main.tex` | 79365 | `54577845049ab85c2c3cb1e5a0acc9acbaae88526510a695ce155a07445ac84a` |
| `paper/refs.bib` | 2656 | `97580d0cd5975b22d6816c541c4835533fbfe4d8e5dee216a946e00225a9f51b` |
| `replay/replay.py` | 30284 | `4de1a9451dab5715a74001bad155b33da4a8b11f88b474dc9ed43640036d54fb` |
| `replay/README.md` | 3822 | `712c891ab3716dec163650ebd0e0dec893c87d6ebc8b7e3ec4b11dfc1bca5314` |
| `replay/SHA256SUMS` | 152 | `9f1fa739d0b3e245090bbdb64f523ee843010a03009191a5216313f9da03ea07` |

`MANIFEST.sha256` holds the digest of every tracked file except itself; verify
with `sha256sum -c MANIFEST.sha256` on Linux, or with Python's `hashlib` on
Windows, where MSYS `sha256sum` reads the PDF in text mode and reports a spurious
mismatch.

## 6. Continuous integration

Two jobs run on every push and on the tag:

1. *Replay the certificate program on the release archive*: verifies
   `MANIFEST.sha256` and `replay/SHA256SUMS` on the checked-out tree, checks
   that the manifest covers every tracked file, builds the release archive with
   `git archive`, unpacks it and verifies both pin files again inside it,
   installs `sympy`, `mpmath`, `numpy` and `scipy`, and runs
   `replay/replay.py`, `replay/replay.py --all-mutants` and
   `replay/replay.py --numerical` on the unpacked tree, requiring each to exit 0
   and the default, mutant and numerical outputs to carry their expected lines.
2. *Verify the released PDF and its text*: extracts the text with `pdftotext`,
   requires the title, the byline, the ORCID, the AI disclosure, the peer-review
   sentence, the DOI placeholder phrase and the cited identifiers to be present,
   rejects an em dash in the text, and records the SHA-256 and MD5 of the
   released PDF. The job never rebuilds the PDF.
