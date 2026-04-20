# Changelog

## 0.2.1

### Performance

- **`advance()` is now callback-based.** The old signature returned `Vec<(Connection, i32)>` that `step_one()` consumed immediately. The new signature takes `impl FnMut(Connection, i32)` and streams each edge into the caller's sink, eliminating the per-step heap allocation.
- **`backtrace_map` densified.** The `FxHashMap<(i32, u32, AlnLayer), ...>` is replaced by three dense `Vec<Vec<PredEntry>>` (one per layer, outer `[y]`, inner `[k + k_offset]`, with `layer == PRED_NONE` marking absent). Each cell is 16 bytes, so memory is proportional to `q_len × k_width × 3`; hash lookups become direct Vec indexing on the hot path. Every backtrace hit inside `next()` and `step_one()` is now a pointer dereference.
- **`reduce()` scratch buffer reused.** `kdist` is now a `Vec<u32>` field on `WaveFronts` (size `k_width`, sentinel `u32::MAX`). `reduce()` writes into it, reads for the `new_kmin`/`new_kmax` search, then resets only the slots it touched — no per-call allocation.
- **Preallocate `s_k_to_y` and `score_to_k_range` rows upfront.** `WaveFront::new_with_capacity` now eagerly fills `capacity` rows instead of lazily pushing one per score. Removes per-step malloc jitter when the caller provides a meaningful capacity hint.
- **Dependencies pruned.** `rustc-hash` and `simple_logger` removed from `Cargo.toml` — both were only used by the HashMap-era code and the old fixture-only tests.

### Breaking change

- `WaveFronts::backtrace_map` (public field) is removed. The new dense store is private and accessed through `backtrace()`. If downstream code read the field directly, call `backtrace()` instead.

## 0.2.0

### Performance

- **Flat-`Vec` storage for the two hot data structures.** `s_k_to_y` is now a `Vec<Vec<u32>>` (outer indexed by `score`, inner a fixed-width window of size `t_len + q_len + 3` indexed by `k + k_offset`, with sentinel `u32::MAX` for absent). `score_to_k_range` is a `Vec<Option<(i32, i32)>>` indexed by `score`. Both replace `FxHashMap`, removing hash overhead from the inner loop of `next()`, `advance()`, and `reduce()` — the hottest paths in the algorithm. Out-of-window reads (e.g. `k±1` at the boundary) return `None` as before; writes are guaranteed in-window by the existing `x ≤ t_len && y ≤ q_len` guard.
- **`advance()` takes `&[u8]` directly.** `step_one()` passes `target_str.as_bytes()` / `query_str.as_bytes()` once per step instead of shadowing inside the inner loop.
- **Dead code removed.** The `if e.1 > score` branch in the insert/delete blocks (unreachable given monotonic `score`) and the asymmetric `is_none()` guard on the match block collapsed into single `entry(...).or_insert(...)` / `contains_key(...)` calls.

### Bug fixes

- **`reduce()` off-by-one on the upper-bound search.** The farthest-right passing diagonal was silently excluded from the reduced k-range because the half-open upper bound wasn't made exclusive after a passing break. Fixed by `new_kmax += 1` on the passing-break path.
- **`reduce()` collapsed the range when no match-layer entries existed at the current score.** With an empty `kdist` the new-bounds search would drift to a single-cell range and effectively kill the wavefront. `reduce()` now early-returns when `kdist` is empty and leaves the range untouched.
- **`reduce()` could prune the target diagonal `k_end = query_len - target_len`** out of the active range, leaving the algorithm unable to close out the alignment. `reduce()` now always preserves `k_end` in `[new_kmin, new_kmax)`.
- **`step_one()`'s `Fail` heuristic (`min(|k_max|, |k_min|) > max(|t|, |q|)`) fired spuriously on tiny inputs**, causing `"A"` vs `"T"` and single-char empty-vs-nonempty cases to return `Fail` before the cheapest edit was reachable. The heuristic was unsound in general and has been removed; `max_score` remains as the termination safety net.
- **`wfs.score` at `ReachEnd` could overstate the cost of the alignment returned by `backtrace()`** under tight `max_wf_length`. Internal state-layer transitions in the forward DP contributed to the chain cost without producing visible alignment columns. `step_all()` now rescores the returned alignment under gap-affine rules on `ReachEnd` and overwrites `self.score`, so the reported score always matches the cost of the alignment the user sees.

### Tests

Invariant-based test suite added in `src/lib.rs` covering round-trip identity (`t_aln.replace("-","") == t_str`), equal-length alignment strings, score consistency (`score_alignment(backtrace) == wfs.score`), single-character and empty-string edge cases, tight `max_wf_length` behavior, `ReachMaxScore` on unreachable scores, determinism across repeated runs, and asymmetric-length indels.

### Notes

See `IMPROVEMENT.md` for the original bug analysis, remaining cleanup items, and performance opportunities that have not been implemented.
