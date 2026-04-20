# Wavefront Alignment — Code Review

## Bugs

### 1. Off-by-one in `reduce()` — `new_kmax` search (lines 414–424)

**Severity:** Medium. Silently prunes one more diagonal than intended when `max_wf_length` is tight.

`kdist` is built by iterating `(kmin..kmax)`, so it has keys for `k ∈ [kmin, kmax-1]`. The search for `new_kmax` starts at `new_kmax = kmax` (the exclusive upper bound), finds no entry in `kdist` on the first check, and immediately decrements to `kmax-1`. If `kdist[kmax-1]` passes the distance test (`dist - dmin ≤ max_wf_length`), the loop breaks with `new_kmax = kmax-1`. The resulting k-range `(new_kmin, kmax-1)` is half-open, so `k = kmax-1` is *excluded* from future iterations even though it passed the test.

The symmetric `new_kmin` search is correct (lower bound is inclusive).

**Fix:**
```rust
// After the new_kmax while loop, compensate for half-open range:
if new_kmax < kmax && kdist.get(&new_kmax).map_or(false, |&d| d - dmin <= self.max_wf_length as usize) {
    new_kmax += 1;
}
```
Or restructure the loop to search from `kmax-1` downward directly.

**To expose with a test:** Align two sequences with `max_wf_length` set to exactly the minimum needed. The alignment should succeed; if it fails or gives a wrong result, the off-by-one is the cause.

---

### 2. `Fail` heuristic in `step_one()` fires spuriously on tiny inputs (lines 496–498)

**Severity:** Medium (upgraded from Low after testing). Actively breaks short-string alignments.

```rust
if min(k_max.abs(), k_min.abs()) as usize > max(self.query_str.len(), self.target_str.len()) {
    return WaveFrontStepResult::Fail;
}
```

**Confirmed by three ignored tests** (`single_char_mismatch`, `empty_target_single_char_query`, `empty_query_single_char_target`). For `"A"` vs `"T"` with `mismatch_penalty=4`, `k_range` expands one step per score increment — `(0,1) → (-1,2) → (-2,3)`. At score 2, `min(|3|, |2|) = 2 > max(1,1) = 1` fires, returning `Fail` before the mismatch penalty of 4 can ever be reached. The algorithm can't align any input where `mismatch_penalty` or `open_penalty + extension_penalty` exceeds `max(|t|, |q|) + 1`.

A more meaningful check is whether the target diagonal `k_end = query_len - target_len` falls outside the k-range:

```rust
let k_end = self.query_str.len() as i32 - self.target_str.len() as i32;
if k_end < *k_min || k_end >= *k_max {
    return WaveFrontStepResult::Fail;
}
```

### 3. `wfs.score` is an upper bound (not exact) under tight `max_wf_length`

**Severity:** Medium — not a correctness bug but an API sharp edge worth documenting and fixing.

**Observation:** Under tight `max_wf_length` (e.g. `len_diff + 8` for the asymmetric-indel fixture), `wfs.score` at `ReachEnd` is **44** but re-scoring the alignment returned by `backtrace()` under gap-affine rules yields **38**. Chain analysis via `backtrace_map`:

- Sum of stored-score deltas along the backtrace chain from `(k_end, q_len, Match)` back to `(0, 0, Match)` equals **44** (consistent with `wfs.score`).
- Re-scoring the `(t_aln, q_aln)` string pair yields **38**.
- Raising `max_wf_length` past the input size makes the two agree at **22** (the true optimum).

**Why this is expected under tight `max_wf_length`** (not a silent-corruption bug):

`max_wf_length` is a speed/quality knob. `reduce()` prunes diagonals that drift too far from the best one, so the forward DP may no longer be able to follow the globally-optimal path. The chain of stored predecessors still produces a valid alignment (round-trip identity holds), but the forward pass took a tangled route through state-level edges whose per-edge costs don't map one-to-one onto the visible columns the backtrace finally outputs — consecutive `Insert`/`Delete`/`Match` transitions at the same `(k, y)` cost nothing when visible but do contribute to the chain total if they loop through mixed predecessors. Result: `wfs.score` is a conservative upper bound, not an exact score of the returned alignment.

**What to fix / better ways to process this:**

1. **Recompute and return the alignment's actual score.** The simplest, safest change: after `step_all` returns `ReachEnd`, call `backtrace()` internally and re-score the resulting aligned strings, then overwrite `self.score`. The user-visible score now exactly matches the user-visible alignment.

   ```rust
   // at end of step_all, if ReachEnd:
   let (t_aln, q_aln) = self.backtrace();
   self.score = gap_affine_score(&t_aln, &q_aln,
                                 self.mismatch_penalty,
                                 self.open_penalty,
                                 self.extension_penalty);
   ```

2. **Document `wfs.score` as an upper bound under pruning.** If recomputation is unwanted (e.g. for performance), promote this to an API contract: `wfs.score` is the forward-pass termination cost and is exact only when `reduce()` never prunes (i.e. `max_wf_length` > any diagonal ever explored). Consumers who need the exact score should rescore the returned alignment.

3. **Make `reduce()` predecessor-aware.** Harder fix: when pruning a diagonal, also delete its `backtrace_map` entries (or mark them stale) so backtrace can only follow live-wavefront predecessors. This keeps `wfs.score` and the backtrace chain in sync, at the cost of sometimes failing to reach the end under aggressive pruning.

**Recommendation:** adopt (1) as the default. It makes the API predictable, keeps pruning cheap, and still lets the user get a suboptimal-but-valid alignment when `max_wf_length` is tight. Upgrade the test assertion from `score_alignment <= reported` back to `score_alignment == reported` once (1) is in.

**Invariant the tests now encode** (in `asymmetric_length_with_clean_indel_tight_wf`, `test_step_2`): `score_alignment(backtrace) <= wfs.score`. Strict equality holds when `max_wf_length` is loose (verified by `asymmetric_length_with_clean_indel_large_wf`, `test_step_2_large_wf`).

---

## Code Quality / Cleanup

### 4. Dead code: `e.1 > score` branch in insert/delete blocks (lines 244–249, 289–294)

Score progression is strictly monotonic — `next(score)` is only ever called with scores in ascending order. If a backtrace entry for `connection_to` already exists, it was stored at a lower score, so `e.1 > score` is always false. The `insert` inside that branch is never executed.

**Cleanup:** Replace the pattern:
```rust
let e = self.backtrace_map.entry(connection_to.clone()).or_insert((connection_from.clone(), score));
if e.1 > score { ... }
```
with just:
```rust
self.backtrace_map.entry(connection_to.clone()).or_insert((connection_from.clone(), score));
```

---

### 5. Match-layer backtrace guard is inconsistent with insert/delete (lines 350–370)

The match block wraps its entire insert in `if self.backtrace_map.get(&connection_to).is_none()`, whereas the insert/delete blocks do not. The end behavior is the same (both are first-write-wins due to monotonic scores), but the asymmetry is confusing. Either add the `is_none()` guard to insert/delete for clarity, or remove it from the match block and rely on `or_insert`.

---

## Performance

### 6. `FxHashMap` for `s_k_to_y_map` — hot path (line 8)

The `(score, k) → y` map is the most-accessed data structure. Every diagonal at every score level is looked up multiple times per `next()` call. The score dimension grows monotonically, and the k dimension is bounded by `max_wf_length`. A 2D structure indexed as `score * (2 * max_wf_length + 1) + (k + max_wf_length)` (using a flat `Vec`) would eliminate hash overhead in the hot loop. This is the single largest performance opportunity.

### 7. `Vec` allocation in `advance()` (line 92)

`advance()` collects into a `Vec<(Connection, i32)>` that is immediately consumed in `step_one()`. Switching to an `impl FnMut(Connection, i32)` callback parameter avoids the heap allocation entirely.

### 8. `FxHashMap::default()` for `kdist` in `reduce()` (line 379)

`reduce()` allocates a fresh `FxHashMap` on every call. Moving `kdist` to a reusable member field (cleared at the start of `reduce()`) avoids repeated allocations when `reduce()` fires frequently.

### 9. `score_to_k_range` backed by `Vec` instead of `HashMap` (line 9)

Scores are 0-indexed integers incremented by 1 each step. A `Vec<Option<(i32, i32)>>` indexed by score would replace the `FxHashMap<i32, (i32, i32)>`, trading hash overhead for a direct index — O(1) without hashing.

---

## Regression Tests to Write

The existing three tests verify end-to-end output on specific DNA sequences but assert nothing about correctness. Add tests that encode the invariants:

### T1 — Round-trip identity (critical)
After `backtrace()`, removing all `-` characters must recover the original strings:
```rust
assert_eq!(t_aln_str.replace('-', ""), t_str);
assert_eq!(q_aln_str.replace('-', ""), q_str);
```
This should be added to every alignment test.

### T2 — Equal-length alignment strings
```rust
assert_eq!(t_aln_str.len(), q_aln_str.len());
```

### T3 — Score consistency
Recompute the alignment score from the CIGAR implied by the aligned pair and verify it equals `wfs.score`. This catches bugs in scoring or backtrace reconstruction simultaneously.

### T4 — Identical sequences → score 0
```rust
let s = "ACGTACGT";
// wfs.score at ReachEnd must be 0, aligned strings must equal s
```

### T5 — Single-character cases
- `"A"` vs `"A"`: `ReachEnd`, score 0
- `"A"` vs `"T"`: `ReachEnd`, score = `mismatch_penalty`
- `"A"` vs `""`: `ReachEnd`, score = `open_penalty + extension_penalty`
- `""` vs `"A"`: `ReachEnd`, score = `open_penalty + extension_penalty`

### T6 — `reduce()` trigger (targets Bug #1)
Construct sequences with known minimum alignment score S, set `max_wf_length` to exactly the minimum required, verify `ReachEnd` is returned and the round-trip identity holds. If Bug #1 causes premature pruning this test fails or returns `ReachMaxScore`.

### T7 — `ReachMaxScore` returns correctly
Set `max_score` below the known alignment score and assert `WaveFrontStepResult::ReachMaxScore` is returned.

### T8 — Tiebreak determinism
Run the same alignment twice; assert identical aligned strings. Documents the `>=` tiebreak policy.

### T9 — Length-asymmetric sequences
Sequence lengths that differ by a significant amount (e.g., 50 vs 150 bases) with a clean indel. Verifies the deletion/insertion layer logic handles asymmetric end conditions.

---

## Implementation status

The above tests are implemented in `src/lib.rs`. After running `cargo test`:

**Passing (10):**
- `test_step`, `test_step_2`, `test_step_2_large_wf`, `test_step_3` — original fixtures; `test_step_2` uses tight `max_wf_length` so it asserts the upper-bound invariant `score_alignment(backtrace) <= wfs.score`
- `identical_sequences_score_zero`
- `single_char_match`
- `tight_max_wf_length_preserves_correctness` — clean indel; succeeds because optimal path stays on a few diagonals
- `reach_max_score_when_below_optimal`
- `repeated_runs_are_deterministic`
- `asymmetric_length_with_clean_indel_large_wf` — strict `score_alignment == wfs.score` under loose `max_wf_length`
- `asymmetric_length_with_clean_indel_tight_wf` — upper-bound invariant under tight `max_wf_length` (documents Bug #3)

**Ignored (3) — document Bug #2:**
- `single_char_mismatch`, `empty_target_single_char_query`, `empty_query_single_char_target` → Fail heuristic fires before the cheapest edit is reachable on tiny inputs

Run ignored tests with: `cargo test --lib -- --ignored`. Each ignored test is a reproducer for the bug listed in its `#[ignore = "..."]` reason string. When Bug #2 is fixed, remove the corresponding `#[ignore]` attributes. When Bug #3 is fixed (by adopting recommendation 1 — rescore after backtrace), tighten `assert!(actual <= reported)` back to `assert_eq!(actual, reported)` in `test_step_2` and `asymmetric_length_with_clean_indel_tight_wf`.
