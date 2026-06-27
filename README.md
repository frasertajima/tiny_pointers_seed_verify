# seedverify — amortising an expensive verifier with a cheap token gate

> A study of **when** it pays to replace a brute-force `O(N)` scan with a cheap discrete
> shortlist plus an exact check on the survivors — and, just as importantly, **when it
> does not.** This document is meant to be read start-to-finish: it builds the idea from
> first principles, derives the one inequality that governs the whole thing, measures it,
> and then spends as much effort on the failure modes as on the wins.

---

## Table of contents

1. [The problem: when comparison is the bottleneck](#1-the-problem-when-comparison-is-the-bottleneck)
2. [The move: seed and verify](#2-the-move-seed-and-verify)
3. [The economic model (the thesis, derived)](#3-the-economic-model-the-thesis-derived)
4. [What the lab measures](#4-what-the-lab-measures)
5. [Where the simple model breaks — three caveats that keep it honest](#5-where-the-simple-model-breaks--three-caveats-that-keep-it-honest)
6. [The fit test](#6-the-fit-test)
7. [Where this points (natural fits)](#7-where-this-points-natural-fits)
8. [Usage](#8-usage)
9. [Files](#9-files)

---

## 1. The problem: when comparison is the bottleneck

A surprising number of workloads have the same shape. You have a large set of `N` items.
For a query (often *each* item, giving `N²`), you must run an **expensive, exact pairwise
operation** and keep the best matches:

- compare a DNA read against a genome database (sequence alignment),
- decide whether two customer records are the same person (edit distance, a learned matcher),
- find the near-duplicate of a document in a web-scale crawl,
- rank passages for a question with a cross-encoder, or have an **LLM judge** each candidate,
- classify a handwritten digit by its distance to every training example (`k`-NN).

The naïve solution computes the expensive operation `N` times per query. When `N` is large
and the operation is costly, that is simply intractable — you cannot run an LLM against ten
million documents to answer one question.

The key observation: **most of those `N` comparisons are wasted.** The overwhelming
majority of items are obviously irrelevant; the expensive operation spends its time
confirming non-matches. If we could *cheaply* discard the obvious non-matches first, we
would only pay the expensive price on a handful of plausible candidates.

## 2. The move: seed and verify

That is the entire idea, and it is old and well-proven:

> **Seed** — a cheap, approximate, discrete filter that proposes a small candidate set.
> **Verify** — the exact, expensive operation, run *only* on those candidates.

You see it everywhere good systems live: **BLAST** seeds on exact `k`-mer matches and only
then extends with expensive alignment; **record linkage** "blocks" records by a cheap key
and only then compares within a block; **image / vector retrieval** shortlists with an
approximate index and then re-ranks exactly. Our own sibling labs are the same engine in
different clothes:

| lab | item | shingle (the window) | token |
|---|---|---|---|
| text search (COBOLMM) | a document | 3 characters | a folded trigram bucket |
| `kmerstash` | a DNA read | `k` bases | a 2-bit-packed `u64` `k`-mer |
| `digitstash` | an MNIST image | a 3×3 pixel patch | 9 bits → a pattern in `0..512` |
| **`seedverify`** | anything | any shingling | any small integer token |

The mechanism is always: **shingle** each item into a set of small integer **tokens**;
build an **inverted index** (token → list of item ids containing it); weight each token by
its rarity (**IDF**, `ln(N/df)`) so that ubiquitous tokens count for nothing and rare,
discriminative ones dominate; for a query, walk the posting lists of its tokens,
accumulate IDF mass onto the items they hit, and return the **top-`budget`** items. Those
are the seed. The verify then runs the real, exact operation on just those `budget` items.

Crucially, **the verify is unchanged** — it computes exactly what the brute-force scan
would. The index never approximates the *answer*; it only approximates *which items are
worth the expensive operation*. That is the honest difference from a pure approximate-NN
system, and it means every error is a single, measurable thing: *did the shortlist keep
the true best item?* (its **recall**).

## 3. The economic model (the thesis, derived)

This lab exists to answer one quantitative question: **how much of the algorithmic saving
actually shows up as wall-clock time?** The two quantities are not the same, and the gap
between them is the whole story.

Let:

- `N` = number of items, `B` = candidate budget (`B ≪ N`),
- $C_v$ = wall-clock cost of **one** exact verify operation,
- $W_{\text{seed}}$ = wall-clock cost of generating the shortlist for one query.

**Brute force** runs `N` verifies per query:

$$T_{\text{exhaustive}} = N \cdot C_v$$

**Seed-and-verify** pays the seed once, then `B` verifies:

$$T_{\text{stash}} = W_{\text{seed}} + B \cdot C_v$$

The **operation-count speedup** — the clean algorithmic number everyone quotes — is just

$$S_{\text{op}} = \frac{N}{B}$$

But the **wall-clock speedup** is what you actually feel:

$$S_{\text{wall}} = \frac{T_{\text{exhaustive}}}{T_{\text{stash}}}
= \frac{N \cdot C_v}{W_{\text{seed}} + B \cdot C_v}$$

Divide top and bottom by $B \cdot C_v$ and write the **per-candidate seed cost**
$w = W_{\text{seed}} / B$:

$$\boxed{\,S_{\text{wall}} = \frac{S_{\text{op}}}{1 + \dfrac{w}{C_v}}\,}
\qquad\Longleftrightarrow\qquad
\frac{S_{\text{wall}}}{S_{\text{op}}} = \frac{C_v}{C_v + w}$$

That ratio on the right — **the fraction of the algorithmic saving you keep** — is the
thesis in one expression. It is a saturating function of $C_v / w$, the ratio of
**verification cost to seeding cost**:

- when $C_v \gg w$ (verify is expensive relative to a token lookup) → fraction → **1**: you
  keep almost the entire `N/B` speedup;
- when $C_v \approx w$ → fraction ≈ **½**;
- when $C_v \ll w$ (verify is nearly free) → fraction → **0**: the seed overhead dominates
  and the algorithmic win evaporates in practice.

**This is why MNIST was "disappointing".** There, `digitstash` reported ~300× fewer
operations but only ~13× wall-clock — a kept fraction of `13/300 ≈ 0.043`. Plug that in:
$w / C_v \approx 22$. The per-candidate seed cost was about **twenty times** the verify
cost. That is not the method failing; it is a direct, predictable consequence of running
it with the cheapest verify in existence (an `L₂` distance). MNIST is the **left edge** of
the curve. The whole point of this lab is to traverse the curve and watch the fraction
climb toward 1 as $C_v$ grows.

## 4. What the lab measures

Two experiments, both pure-Rust, std-only, reproducible from a seed.

### E1 — the cost-sweep: dialling $C_v$ directly

A synthetic token corpus; each query is a noisy copy of one target document (so a small
shortlist *should* keep it). The verify is exact Jaccard similarity **plus a tunable
busy-work knob** (`--verify-work`) that stands in for everything from an almost-free
distance to a seconds-scale model call. We sweep the knob and watch $S_{\text{wall}}$ rise
toward the flat $S_{\text{op}} = N/B$ ceiling.

```
  N=20000, budget=100  ⇒  op-count ceiling N/B = 200×
  verify_work    exhaustive     stash      wall speedup     fraction kept     accuracy
  ----------------------------------------------------------------------------------------
       0           240 ms       2.7 ms        90×              45%            1.00 / 1.00
      50           271 ms       2.2 ms        126×             63%            1.00 / 1.00
     200           465 ms       3.4 ms        137×             68%            1.00 / 1.00
     800          1226 ms       7.5 ms        163×             81%            1.00 / 1.00
    3200          4444 ms      22.7 ms        196×             98%            1.00 / 1.00
   12800         16780 ms      85   ms        198×             99%  ◄──       1.00 / 1.00
```

The curve is exactly $S_{\text{op}} / (1 + w/C_v)$: a kept fraction that saturates at 1 as
the verify gets expensive. (Note the *floor* here is 90×, not MNIST's 13×, because this
scalar Jaccard verify is already costly relative to our cheap seed — see §5.1. The
*shape*, not the floor, is the universal part.)

### E2 — record linkage: a *real* expensive verify

The textbook non-ML fit, with nothing stubbed. Synthetic person records
(`"Mary Tajima 4821 Lakeside"`), typo'd query copies. Seed = **char-trigram blocking**
(IDF, df-capped). Verify = **Levenshtein edit distance**, genuinely `O(L²)` — the classic
blocking-then-compare comparator.

```
  20000 records, 50 typo'd queries (2 edits each), budget = 50
  exhaustive : recall@1 100%   20000 edit-distance calls/query   17330 µs/query
  stash      : recall@1 100%      50 edit-distance calls/query      80 µs/query
  ⇒ 400× fewer comparisons, 215× wall-clock, candidate recall 100%, zero accuracy lost
```

Edit distance is expensive enough that $C_v \gg w$, so the wall-clock speedup (215×) lands
close to the comparison-count ratio (400×) — the right edge of the E1 curve, reached with
a real algorithm instead of a knob.

### E3 — corpus-wide LLM alert scan: the expensive verify is an *inference*

The payoff case, on real data and a real model. The corpus is real SEC 10-K risk-factor
passages (COBOLMM's `SEC10K.DAT`, 40 tickers). The verify is a **local LLM** (`gemma4:e2b`
via ollama) judging *"does this passage discuss `<risk>`?"*. The seed is char-trigram
**blocking** on the alert phrase. An exhaustive scan asks the model about **every** passage;
the stash asks it only about the shortlist. Measured over 200 passages, alert = *"material
weakness in internal control over financial reporting"*, ≈620 ms per LLM call:

```
                  LLM calls   matches   wall time
  exhaustive         200        10        124 s     (the intractable O(N) baseline)
  stash  (budget 10)  10         7          6 s      20× fewer calls,  recall 70%
  stash  (budget 40)  40         9         25 s       5× fewer calls,  recall 90%
  stash  (budget 80)  86        10         53 s     2.3× fewer calls,  recall 100%
  ⇒ projected to a 10,000-passage corpus: exhaustive ≈ 103 min, stash ≈ 12–53 s
```

Two lessons, both honest, both predicted by the theory above:

1. **Wall-clock tracks calls exactly.** At ≈620 ms each, the inference utterly dwarfs the
   seed, so $C_v \gg w$ and the op-count saving *is* the wall-clock saving — the far right
   of the E1 curve, now with a real model. This is the "LLM too slow corpus-wide" problem
   (alert scanning, filing scoring) turned from intractable into an overnight-or-less batch.

2. **The recall ceiling (§5.3) is real, visible, and tunable.** A *lexical* trigram seed
   under a *semantic* verify first found only **4 of 10** — the LLM accepts paraphrases
   ("deficiencies in our financial reporting controls") that share no rare trigrams with
   the alert phrase. The built-in recall-vs-budget sweep (a free, no-extra-inference
   computation once the matches are known) exposes the whole speed↔recall trade, and
   **seed expansion** — adding synonym terms to the *blocker*, never to the LLM prompt —
   lifted the ceiling from **40% → 100%**. Where recall still falls short at a large budget,
   the matches are *genuinely* semantic and you have found exactly the point at which to
   swap in an embedding seed. The engine tells you when you have hit that wall.

```bash
seedverify bench --kind alert --backend stub    # reproducible: keyword oracle, no model
seedverify bench --kind alert --backend ollama --model gemma4:e2b   # the real verify
```

### E4 — a dense seed + hybrid fusion: recovering the semantic matches

E3 left a clean problem: a *lexical* trigram seed under a *semantic* LLM verify is blind to
paraphrases (§5.3). E4 adds the other arm of a RAG stack — a **dense embedding seed**
(`embeddinggemma` via ollama, cosine over cached 768-d vectors) — and fuses it with the
trigram seed by **weighted reciprocal-rank fusion**. Verify and ground truth are unchanged:
the LLM labels every passage once (cached), so the per-seed recall sweep is then free.

Two alerts expose two regimes (200 SEC passages, `gemma4:e2b` verify):

```
  "material weakness…"  (filings use the phrase near-verbatim — lexically anchored)
    budget   lexical  embedding  RRF(equal)  RRF(emb×3)
      10       78%       89%        100%        100%
      20       89%      100%        100%        100%

  "depends on a small number of large customers"  (filings say "customer concentration",
                                                   "a limited number of customers", …)
    budget   lexical  embedding  RRF(equal)  RRF(emb×3)
      10        0%       33%          0%         33%
      20        0%       67%         33%         67%
    ceiling     33%     100%        100%        100%   ← lexical reaches 1 of 3; dense, all 3
```

Three findings, each one a reason the lab was worth building:

1. **The dense arm is essential for semantic alerts.** Where the lexical seed caps at 33%
   (blind to paraphrases), the embedding seed reaches 100% — exactly the boundary §5.3
   predicted, now crossed with the right tool.
2. **Naive equal-weight fusion can *hurt*.** On the semantic alert, blending the weak lexical
   arm into the strong dense one **dilutes** it (budget 20: equal-RRF 33% vs embedding 67%).
   "Always blend" is not unconditionally right.
3. **Weight the dense arm.** RRF with the embedding arm weighted ×3 tracks the *better* arm
   in **both** regimes (lexical-anchored: matches equal-RRF's 100%; semantic: recovers the
   embedding's 67%). That is the robust default.

This is the full local RAG stack — **sparse retriever ∪ dense retriever → LLM reader** —
each tool doing the bit it is best at, with the fusion weighted so the dominant arm wins.
The embedding model is a *different, small* local model (a decoder LLM cannot serve good
embeddings); the corpus vectors are embedded once and cached, so the dense arm's `O(N)` cost
is paid a single time.

```bash
seedverify bench --kind hybrid --backend ollama --model gemma4:e2b --embed-model embeddinggemma
```

### E5 — `scan`: the production tool, and the COBOLMM menu item

E1–E4 are *benchmarks* (they run an exhaustive pass to measure recall). The `scan`
command is the **production** tool: it does **no** exhaustive pass — it builds the
dense-weighted hybrid shortlist, has the local LLM judge only those, and prints the
passages it flags. That is the whole proposition, delivered: spend the slow model on a
handful of candidates, not the corpus.

```bash
seedverify scan --alert "material weakness in internal control over financial reporting" \
                --sec /…/SEC10K.DAT --budget 40 --backend ollama \
                --model gemma4:e2b --embed-model embeddinggemma --emb-weight 3
# → prints the matching passages (ticker + snippet) found via 40 LLM calls, not thousands
```

It is wired into **COBOLMM** as menu item **AL** ("Scan filings for an alert"): the
launcher (`main_menu/alert_scan/`) resolves the configured local model (`LLM_MODEL_SMALL`),
prompts for a corpus (10-K / 10-Q), an alert phrase, and a budget, then runs the engine.
Measured through the live menu: a material-weakness alert over 600 SEC passages returned 12
relevant matches via 20 LLM calls — **30× fewer** than judging every passage. The corpus is
embedded once and cached, so subsequent scans are seed-fast. An overnight corpus-wide LLM
sweep becomes a coffee-break query, with the budget as the analyst's coverage dial.

## 5. Where the simple model breaks — three caveats that keep it honest

The boxed formula is a model, and a model's value is knowing where it lies to you. Three
system-level effects complicate the clean curve; a careful lab keeps an eye on all three.

### 5.1 Hardware amortisation & SIMD — "operations" are not all priced equally

The model treats $C_v$ and $w$ as wall-clock costs, and that hides a brutal asymmetry on
real CPUs. An `L₂` distance over a dense float vector is **embarrassingly parallel**: a
modern compiler vectorises the loop with AVX-512/AMX, crunching many lanes per cycle, with
perfectly predictable, contiguous memory access and prefetching. The verify's *effective*
$C_v$ is therefore far smaller than its raw operation count suggests.

The seed does the opposite. Traversing sparse posting lists is **pointer-chasing**: random
memory access, cache misses, branch mispredictions on the df-cap test, and a final sort of
the touched set — none of which vectorises. The seed's *effective* $w$ is therefore
**larger** than its operation count suggests.

So the very thing that makes a dense verify "cheap" (it maps onto vector hardware) is also
what makes the algorithmic win hard to realise: the hardware is already executing the
brute-force scan with ruthless efficiency, while it executes the clever index poorly. **A
300× operation reduction yielding only a 13× wall-clock speedup is the hardware
compensating for brute force.** Practical consequence: this engine pays off most against
verifies that are *intrinsically* sequential and un-vectorisable — alignment, edit
distance, tree/graph walks, a model forward pass — and least against ones that are a tight,
vectorisable numeric kernel. Always quote *wall-clock* costs in the ratio, never op counts.

### 5.2 The index-overhead threshold — $W_{\text{seed}}$ is not constant in $N$

The derivation quietly treated $W_{\text{seed}}$ as fixed. It is not. The seed walks the
posting lists of the query's tokens, and a posting list's length is the token's document
frequency `df`. We cap traversal at `df_cap = frac · N` — but that cap **scales with `N`**.
As the corpus grows, the surviving "moderately common" tokens have longer and longer
posting lists, the touched-candidate set grows, and the final top-`B` sort grows with it.
In the worst case (queries full of mid-frequency tokens) the per-query seed cost trends
**linear in `N`**, the very thing we were trying to escape.

For small corpora (MNIST's 60k, our 20k) this is negligible — $w \ll C_v$ at any real
verify cost, so it never shows. At **corporate / web scale** (`N` in the millions) the
seed can develop its own bottleneck and the realised speedup *plateaus* well below `N/B`.
This lab deliberately reports `seed_walked` (posting entries traversed per query) so the
term is visible, not hidden. Mitigations are well-known and are the next engineering
frontier if you push `N` up: tighter `df_cap`, skip-pointers / block-max (WAND) traversal,
early-termination once the top-`B` cannot change, capping or sampling pathological lists,
and sharding. The thesis still holds — but $w$ in the formula must be measured at *your*
scale, not extrapolated from a toy one.

### 5.3 The recall ceiling — the third axis you cannot wish away

E1 and E2 both show accuracy pinned at 1.0, which makes it tempting to treat the gate as
free. It is not free; it was *paid for by the data*. The clean formula assumes a **perfect
gate**: the true best item is always inside the shortlist (recall = 1), so the verify on
the shortlist returns the same answer as the full scan. The moment that assumption slips,
shrinking `B` to go faster starts **dropping true matches before the verify ever sees
them** — silent false negatives.

So the real object is not a curve but a **surface with three axes**:

> **Wall-clock speed ⟷ Verifier cost ⟷ Recall ceiling.**

You can buy speed by shrinking `B` (or tightening `df_cap`), but only down to the recall
the *data* allows. That recall ceiling is high exactly when the precondition of §6 (condition 3) holds
— true matches share rare tokens — and it collapses when it doesn't. Our sibling result is
the cautionary tale: `digitstash` on **Fashion-MNIST** found the shortlist could *not* keep
the true neighbour on textured images, and accuracy fell — same engine, recall-limited
data. `seedverify` holds recall at 1.0 only because its synthetic and record data satisfy
the precondition strongly. **Read E1/E2 as the cost axis with recall held fixed by
construction; read Fashion-MNIST as the recall axis. The honest tool lives on the whole
surface.**

## 6. The fit test

Putting §3–§5 together, this engine is the right tool when **all four** hold:

1. **An expensive exact operation runs pairwise over a large set.** (Otherwise there is
   nothing to amortise.)
2. **The verify is costly *relative to a token lookup*** — ideally intrinsically
   sequential / un-vectorisable (§5.1), so $C_v \gg w$ and you keep most of `N/B`.
3. **Items shingle into discrete tokens, and true matches share *rare* ones** — so a small
   shortlist keeps them and the recall ceiling (§5.3) is high.
4. **You need the *exact* answer on the survivors** — so you genuinely verify, which is the
   honest edge over a pure approximate-NN index.

**Anti-fit:** dense embeddings under a trustworthy metric. There the metric is already
good (high recall needs no structural prior), the distance is vectorisable (low $C_v$,
§5.1), and an ANN index (HNSW/FAISS) wins outright with no verify needed. Forcing a
token-gate there buys nothing — it is the left edge of every one of the three axes at once.

## 7. Where this points (natural fits)

Ranked by payoff, leading with the case that maximises $C_v / w$:

- **Corpus-wide LLM operations (the big one).** Verify = an LLM call (seconds). Each pruned
  candidate is an entire inference *not* run, and an inference is maximally
  un-vectorisable-from-your-side, so you keep nearly the whole `N/B` ratio. This turns an
  intractable `O(N)` "ask the model about every document" into an `O(budget)` overnight
  batch — directly the "LLM too slow corpus-wide" problem (alert scanning across a feed,
  scoring every filing). *(This is E3, **measured above** on real SEC filings: 5–20× fewer
  LLM calls, recall tunable to 100% via budget + seed expansion.)*
- **Record linkage / fuzzy join / dedup / entity resolution** — blocking + an exact
  comparator. E2 *is* this.
- **Near-duplicate / clone / plagiarism detection** — shingle → candidates → exact diff.
- **Bioinformatics** — BLAST is seed-and-verify; shown in `kmerstash`.
- **Cross-encoder reranking, threat-intel / log triage** — cheap retrieve, expensive verify.

## 8. Usage

```bash
make            # build (pure Rust, std-only, no crates)
make demo       # run both experiments
make bench      # the experiment CSVs
make test       # unit tests
make run-notebook   # build + execute seedverify_demo.ipynb end-to-end
```

E1 knobs: `--corpus N` (`N`), `--vocab`, `--doclen`, `--queries`, `--keep`/`--noise`
(query overlap with its target), `--budget B`, `--dfcap-frac` (the `frac` in §5.2),
`--works` (comma-separated verify-cost knob values), `--seed`.
E2 knobs: `--records N`, `--queries`, `--edits` (typos per query), `--budget`,
`--dfcap-frac`, `--seed`. The binary prints CSV to stdout and diagnostics to stderr.

## 9. Files

| File | Role |
|---|---|
| `src/index.rs` | the generic IDF token inverted index; `candidates()` returns the shortlist **and** `seed_walked` (the §5.2 term) |
| `src/verify.rs` | the exact ops: `spin` (tunable $C_v$), `jaccard_sorted`, `levenshtein` |
| `src/data.rs` | std-only xorshift RNG + synthetic corpus / queries / typo'd records / char-trigram tokenizer |
| `src/main.rs` | CLI: `bench --kind cost-sweep \| linkage`, `demo` |
| `build_notebook.py` | generates `seedverify_demo.ipynb` (the plotted version of §4) |
| `PLAN.md` | the thesis, the fit test, and the E1/E2/E3 roadmap |

**Honest scope.** The data is synthetic on purpose: the thesis (§3) is about the
*economics* of seed-and-verify — how an operation-count saving converts to wall-clock as
verify cost grows — which is dataset-independent. The three caveats (§5) are where real
deployments differ, and the lab surfaces each term (`seed_walked`, the cost knob, the
recall columns) so you can re-measure them at your own scale rather than trust a toy
extrapolation. Real corpora — records, filings, LLM judging — are instances of the same
curve, read off the same three axes.
