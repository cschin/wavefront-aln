# Changelog

## 0.2.0

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
