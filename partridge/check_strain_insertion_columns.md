# Reading `check_strain_insertion.py`'s output tsv

Each row is one candidate insertion breakpoint from the evidence `.fa`/`.fa.gz`
file. The script tests whether that insertion already exists in the query
genome by aligning its clipped (non-reference) sequence against a window
extracted from both genomes.

## The breakpoint itself

| Column | What it is |
|---|---|
| `seqname`, `pos`, `strand` | The breakpoint's location, taken straight from the evidence file's header. `strand` is the LTR orientation isa_hmmer2 recorded — not related to `query_strand` below. |
| `right_clip_len`, `left_clip_len` | Length (bp) of the clip sequence used for every alignment in this row — the longest `RIGHT_CLIP`/`LEFT_CLIP` read at this breakpoint. |

## Where the query genome was searched

| Column | What it is |
|---|---|
| `query_scaffold`, `query_start`, `query_end` | The region in the query genome that was compared against. If `--same-scaffold` was used, this is just the same coordinates on the same scaffold name. Otherwise, it's wherever minimap2 found the best match for the reference window. |
| `query_strand` | The orientation of that match. Unrelated to the breakpoint's own `strand` column above. |
| `orthology_hit_span_diff` | Only filled in when minimap2 was used to find the query window. It's (matched span in query) − (span searched in reference) from that one alignment. A value close to the insertion's expected size (a few kb for an IAP) is, on its own, a hint the query has extra sequence there — i.e. the insertion. A value near zero means no such gap was found. |

## A coarse sanity check (not used for the verdict)

| Column | What it is |
|---|---|
| `window_identity_pct`, `window_aligned_len` | How well the *whole* reference window matches the *whole* query window, as one direct alignment. This is background context only — it is not used to decide `note`. A low value across most of the window usually just means the neighborhood is repeat-dense, not that anything is wrong. |

## The actual evidence

The clip sequence is aligned separately against each window. This is what
`note` is based on.

| Column | What it is |
|---|---|
| `right_clip_vs_ref_identity_pct`, `left_clip_vs_ref_identity_pct` | How well the clip matches the **reference** window. This is a negative control — the clip is insertion sequence, so it shouldn't match the reference well. A strong match here is a red flag, not good news. |
| `right_clip_vs_ref_aligned_len`, `left_clip_vs_ref_aligned_len` | How many bp of the clip that alignment covered. |
| `right_clip_vs_query_identity_pct`, `left_clip_vs_query_identity_pct` | How well the clip matches the **query** window. This is the real test — a strong, specific match here (high identity **and** covering most of the clip) is what "the insertion is already in the query genome" looks like. |
| `right_clip_vs_query_aligned_len`, `left_clip_vs_query_aligned_len` | How many bp of the clip that alignment covered. |
| `right_clip_vs_query_pos`, `left_clip_vs_query_pos` | Where inside the query window (local coordinates, not genome coordinates) the best match landed. |

## The verdict

| Column | What it is |
|---|---|
| `note` | The only column that makes a call — everything else is the evidence behind it. It compares each `_vs_query` identity/coverage pair against `--min-identity`/`--min-coverage` (default 85% / 60%). |

`note` will read one of:
- **"likely present in query genome"** — both clips clear the threshold against the query window, and neither clears it against the reference. Consistent with a strain-polymorphic insertion, not a real novel/somatic one.
- **"likely absent from query genome (consistent with novel/somatic)"** — neither clip clears the threshold anywhere. Consistent with a genuine novel insertion.
- **"ambiguous — only one clip end matches query genome"** — only `RIGHT_CLIP` or only `LEFT_CLIP` hit, not both. Worth a manual look.
- **"inconclusive (clip also matches reference window — check window/repeat content)"** — a clip matched the *reference* window too well. Usually means the region is repetitive enough that the clip is picking up unrelated copies, not a specific hit; the result can't be trusted as-is.
- **"no orthologous window found"** — minimap2 couldn't locate this locus in the query genome at all (only possible without `--same-scaffold`).

## Rule of thumb

Ignore everything except the four `_vs_query`/`_vs_ref` identity+coverage
pairs and `note` on a first pass. The rest exists so you can check *why* a
row got the verdict it did, or investigate an "ambiguous"/"inconclusive"
row by hand — e.g. by pulling the matching sequences out of the companion
`<output>.windows.fa` file and looking at them directly.
