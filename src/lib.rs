use log::debug;
use rustc_hash::FxHashMap;

use std::cmp::{max, min};

/// Sentinel value used in the flat `s_k_to_y` / `kdist` Vecs to mean "no
/// entry". Chosen so it can never collide with a real `y` (bounded by `q_len`)
/// or a real distance (bounded by `max(t_len, q_len)`).
const SENTINEL: u32 = u32::MAX;

fn gap_affine_score(t_aln: &str, q_aln: &str, mismatch: i32, open: i32, ext: i32) -> i32 {
    let tb = t_aln.as_bytes();
    let qb = q_aln.as_bytes();
    let mut score = 0i32;
    let mut in_gap_t = false;
    let mut in_gap_q = false;
    for i in 0..tb.len() {
        let t_gap = tb[i] == b'-';
        let q_gap = qb[i] == b'-';
        if t_gap {
            score += if in_gap_t { ext } else { open + ext };
            in_gap_t = true;
            in_gap_q = false;
        } else if q_gap {
            score += if in_gap_q { ext } else { open + ext };
            in_gap_q = true;
            in_gap_t = false;
        } else {
            if tb[i] != qb[i] {
                score += mismatch;
            }
            in_gap_t = false;
            in_gap_q = false;
        }
    }
    score
}

/// Dense, score-major storage for one wavefront layer. Both maps are flat Vecs
/// indexed by `score` (outer) and `k + k_offset` (inner); missing entries use
/// the `SENTINEL` value or `None`, avoiding the hash overhead of FxHashMap on
/// the hot path.
struct WaveFront {
    // s_k_to_y[score as usize][(k + k_offset) as usize] = y or SENTINEL
    s_k_to_y: Vec<Vec<u32>>,
    score_to_k_range: Vec<Option<(i32, i32)>>,
    k_width: usize,
    k_offset: i32,
}
#[derive(PartialEq, Eq, Hash, Clone, Debug)]
pub enum AlnLayer {
    Insert,
    Delete,
    Match,
}

#[derive(Eq, PartialEq, Debug)]
pub enum WaveFrontStepResult {
    ReachMaxScore,
    ReachEnd,
    Fail,
    Continue,
}

pub struct WaveFronts<'a> {
    target_str: &'a str,
    query_str: &'a str,
    insertion_layer: WaveFront,
    deletion_layer: WaveFront,
    match_layer: WaveFront,
    pub max_wf_length: u32,
    pub mismatch_penalty: i32,
    pub open_penalty: i32,
    pub extension_penalty: i32,
    pub score: i32,
    /// Sparse backtrace store keyed by the destination anchor. Value holds the
    /// predecessor anchor and the edge score at which it was reached. Only
    /// cells actually visited by the (possibly pruned) wavefront are
    /// materialised, so memory scales with explored cells — not with
    /// `t_len × q_len` as a dense store would.
    backtrace_map: FxHashMap<Anchor, (Anchor, i32)>,
    /// Reusable scratch buffer for `reduce()`; indexed by `k + k_offset`, with
    /// `SENTINEL` marking "no entry". Cleared in-place after each `reduce()`
    /// call so we don't allocate a fresh FxHashMap on every invocation.
    kdist: Vec<u32>,
}
type Anchor = (i32, u32, AlnLayer);
type Connection = (Anchor, Anchor); // Connection.0 is the best score at the end point

impl WaveFront {
    fn new_with_capacity(t_len: usize, q_len: usize, capacity: usize) -> Self {
        // k = y - x with y ∈ [0, q_len] and x ∈ [0, t_len] gives
        // k ∈ [-t_len, q_len]. Add a 1-slot margin on each side so reads at
        // k-1 / k+1 from boundary diagonals stay in-bounds instead of needing
        // extra conditionals in the hot loop.
        let k_offset = t_len as i32 + 1;
        let k_width = t_len + q_len + 3;
        // Reserve spine capacity only; rows are allocated lazily by
        // `ensure_score` as scores are actually reached. Upfront row alloc
        // would cost `capacity × k_width × 4 B` per layer × 3 layers, which
        // is tens of MB per alignment even for modest N.
        let s_k_to_y: Vec<Vec<u32>> = Vec::with_capacity(capacity);
        let score_to_k_range: Vec<Option<(i32, i32)>> = Vec::with_capacity(capacity);
        let mut wf = WaveFront {
            s_k_to_y,
            score_to_k_range,
            k_width,
            k_offset,
        };
        wf.set_range(0, (0, 1));
        wf
    }

    #[inline]
    fn k_idx(&self, k: i32) -> Option<usize> {
        let idx = k + self.k_offset;
        if idx < 0 || (idx as usize) >= self.k_width {
            None
        } else {
            Some(idx as usize)
        }
    }

    fn ensure_score(&mut self, score: i32) {
        if score < 0 {
            return;
        }
        let s = score as usize;
        while self.s_k_to_y.len() <= s {
            self.s_k_to_y.push(vec![SENTINEL; self.k_width]);
        }
        while self.score_to_k_range.len() <= s {
            self.score_to_k_range.push(None);
        }
    }

    #[inline]
    fn get_y(&self, score: i32, k: i32) -> Option<u32> {
        if score < 0 {
            return None;
        }
        let row = self.s_k_to_y.get(score as usize)?;
        let i = self.k_idx(k)?;
        let y = row[i];
        if y == SENTINEL {
            None
        } else {
            Some(y)
        }
    }

    fn set_y(&mut self, score: i32, k: i32, y: u32) {
        self.ensure_score(score);
        let s = score as usize;
        let i = self.k_idx(k).expect("k out of bounds in set_y");
        self.s_k_to_y[s][i] = y;
    }

    #[inline]
    fn get_range(&self, score: i32) -> Option<(i32, i32)> {
        if score < 0 {
            return None;
        }
        self.score_to_k_range
            .get(score as usize)
            .copied()
            .flatten()
    }

    fn set_range(&mut self, score: i32, range: (i32, i32)) {
        self.ensure_score(score);
        self.score_to_k_range[score as usize] = Some(range);
    }

    /// Extend matches along every active diagonal at `score`, streaming each
    /// `(to_anchor, from_anchor)` pair into `sink` as we go. The callback-based
    /// API avoids the intermediate `Vec<(Connection, i32)>` allocation that
    /// the old signature produced on every call.
    fn advance<F>(&mut self, t: &[u8], q: &[u8], score: i32, mut sink: F)
    where
        F: FnMut(Connection, i32),
    {
        debug!("advance score: {}", score);
        let (k_min, k_max) = self.get_range(score).expect("no range");
        for k in k_min..k_max {
            if let Some(y) = self.get_y(score, k) {
                let mut xs = (y as i32 - k) as usize; // k = y - x
                let mut ys = y as usize;
                while xs < t.len() && ys < q.len() && t[xs] == q[ys] {
                    xs += 1;
                    ys += 1;
                }
                debug!("advance: k:{} y: {} -> {}", k, y, ys);
                if ys > y as usize {
                    self.set_y(score, k, ys as u32);
                    sink(
                        ((k, ys as u32, AlnLayer::Match), (k, y, AlnLayer::Match)),
                        score,
                    );
                }
            }
        }
    }
}

impl<'a> WaveFronts<'a> {
    pub fn new(
        t_str: &'a str,
        q_str: &'a str,
        max_wf_length: u32,
        mismatch_penalty: i32,
        open_penalty: i32,
        extension_penalty: i32,
    ) -> Self {
        Self::new_with_capacity(
            t_str,
            q_str,
            max_wf_length,
            mismatch_penalty,
            open_penalty,
            extension_penalty,
            0,
        )
    }

    pub fn new_with_capacity(
        t_str: &'a str,
        q_str: &'a str,
        max_wf_length: u32,
        mismatch_penalty: i32,
        open_penalty: i32,
        extension_penalty: i32,
        capacity: usize,
    ) -> Self {
        let t_len = t_str.len();
        let q_len = q_str.len();
        let insertion_layer = WaveFront::new_with_capacity(t_len, q_len, capacity);
        let deletion_layer = WaveFront::new_with_capacity(t_len, q_len, capacity);
        let mut match_layer = WaveFront::new_with_capacity(t_len, q_len, capacity);
        match_layer.set_y(0, 0, 0);

        let k_width = match_layer.k_width;

        WaveFronts {
            target_str: t_str,
            query_str: q_str,
            insertion_layer,
            deletion_layer,
            match_layer,
            max_wf_length,
            mismatch_penalty,
            open_penalty,
            extension_penalty,
            score: 0,
            backtrace_map: FxHashMap::default(),
            kdist: vec![SENTINEL; k_width],
        }
    }

    // --- Sparse backtrace accessors -------------------------------------------

    #[inline]
    fn bt_get(&self, k: i32, y: u32, layer: AlnLayer) -> Option<(Anchor, i32)> {
        self.backtrace_map.get(&(k, y, layer)).cloned()
    }

    #[inline]
    fn bt_contains(&self, k: i32, y: u32, layer: AlnLayer) -> bool {
        self.backtrace_map.contains_key(&(k, y, layer))
    }

    /// First-write-wins insert used by `next()`. Matches `Entry::or_insert`:
    /// if the destination anchor already has a predecessor we keep the first
    /// one, preserving the optimal-score ordering the algorithm relies on.
    fn bt_insert_if_absent(&mut self, to: &Anchor, from: &Anchor, score: i32) {
        self.backtrace_map
            .entry(to.clone())
            .or_insert_with(|| (from.clone(), score));
    }


    pub fn next(&mut self, score: i32) {
        let (match_k_min, match_k_max) = self
            .match_layer
            .get_range(score - self.mismatch_penalty)
            .unwrap_or((0, 1));

        let (open_k_min, open_k_max) = self
            .match_layer
            .get_range(score - self.open_penalty - self.extension_penalty)
            .unwrap_or((0, 1));

        let (i_k_min, i_k_max) = self
            .insertion_layer
            .get_range(score - self.extension_penalty)
            .unwrap_or((0, 1));

        let (d_k_min, d_k_max) = self
            .deletion_layer
            .get_range(score - self.extension_penalty)
            .unwrap_or((0, 1));

        let k_max = max(match_k_max, max(open_k_max, max(i_k_max, d_k_max))) + 1;
        let k_min = min(match_k_min, min(open_k_min, min(i_k_min, d_k_min))) - 1;
        debug!("score:{} k_min:{}, k_max:{}", score, k_min, k_max);
        self.insertion_layer.set_range(score, (k_min, k_max));
        self.deletion_layer.set_range(score, (k_min, k_max));
        self.match_layer.set_range(score, (k_min, k_max));

        let t_len = self.target_str.len();
        let q_len = self.query_str.len();

        (k_min..k_max).for_each(|k| {
            debug!("k: {:?}", k);

            // Insert: extend query (y advances, x stays).
            let e1 = self
                .match_layer
                .get_y(score - self.open_penalty - self.extension_penalty, k - 1);
            let e2 = self
                .insertion_layer
                .get_y(score - self.extension_penalty, k - 1);

            let connection = match (e1, e2) {
                (Some(y1), None) => Some((
                    (k, y1 + 1, AlnLayer::Insert),
                    (k - 1, y1, AlnLayer::Match),
                )),
                (None, Some(y2)) => Some((
                    (k, y2 + 1, AlnLayer::Insert),
                    (k - 1, y2, AlnLayer::Insert),
                )),
                (Some(y1), Some(y2)) => {
                    if y1 >= y2 {
                        Some((
                            (k, y1 + 1, AlnLayer::Insert),
                            (k - 1, y1, AlnLayer::Match),
                        ))
                    } else {
                        Some((
                            (k, y2 + 1, AlnLayer::Insert),
                            (k - 1, y2, AlnLayer::Insert),
                        ))
                    }
                }
                (None, None) => None,
            };

            if let Some((connection_to, connection_from)) = connection {
                let x = (connection_to.1 as i32 - connection_to.0) as usize;
                let y = connection_to.1 as usize;
                if x <= t_len && y <= q_len {
                    self.insertion_layer.set_y(score, k, connection_to.1);
                    self.bt_insert_if_absent(&connection_to, &connection_from, score);
                }
            }

            // Delete: extend target (x advances, y stays).
            let e1 = self
                .match_layer
                .get_y(score - self.open_penalty - self.extension_penalty, k + 1);
            let e2 = self
                .deletion_layer
                .get_y(score - self.extension_penalty, k + 1);
            let connection = match (e1, e2) {
                (Some(y1), None) => {
                    Some(((k, y1, AlnLayer::Delete), (k + 1, y1, AlnLayer::Match)))
                }
                (None, Some(y2)) => {
                    Some(((k, y2, AlnLayer::Delete), (k + 1, y2, AlnLayer::Delete)))
                }
                (Some(y1), Some(y2)) => {
                    if y1 >= y2 {
                        Some(((k, y1, AlnLayer::Delete), (k + 1, y1, AlnLayer::Match)))
                    } else {
                        Some(((k, y2, AlnLayer::Delete), (k + 1, y2, AlnLayer::Delete)))
                    }
                }
                (None, None) => None,
            };

            if let Some((connection_to, connection_from)) = connection {
                let x = (connection_to.1 as i32 - connection_to.0) as usize;
                let y = connection_to.1 as usize;
                if x <= t_len && y <= q_len {
                    self.deletion_layer.set_y(score, k, connection_to.1);
                    self.bt_insert_if_absent(&connection_to, &connection_from, score);
                }
            }

            // Match: mismatch on diagonal, or free close from Insert/Delete.
            let e1 = self.match_layer.get_y(score - self.mismatch_penalty, k);
            let e2 = self.insertion_layer.get_y(score, k);
            let e3 = self.deletion_layer.get_y(score, k);
            let connection = match (e1, e2, e3) {
                (Some(y1), None, None) => {
                    Some(((k, y1 + 1, AlnLayer::Match), (k, y1, AlnLayer::Match)))
                }
                (None, Some(y2), None) => {
                    Some(((k, y2, AlnLayer::Match), (k, y2, AlnLayer::Insert)))
                }
                (None, None, Some(y3)) => {
                    Some(((k, y3, AlnLayer::Match), (k, y3, AlnLayer::Delete)))
                }
                (Some(y1), Some(y2), None) => {
                    if y1 + 1 >= y2 {
                        Some(((k, y1 + 1, AlnLayer::Match), (k, y1, AlnLayer::Match)))
                    } else {
                        Some(((k, y2, AlnLayer::Match), (k, y2, AlnLayer::Insert)))
                    }
                }
                (Some(y1), None, Some(y3)) => {
                    if y1 + 1 >= y3 {
                        Some(((k, y1 + 1, AlnLayer::Match), (k, y1, AlnLayer::Match)))
                    } else {
                        Some(((k, y3, AlnLayer::Match), (k, y3, AlnLayer::Delete)))
                    }
                }
                (None, Some(y2), Some(y3)) => {
                    if y2 >= y3 {
                        Some(((k, y2, AlnLayer::Match), (k, y2, AlnLayer::Insert)))
                    } else {
                        Some(((k, y3, AlnLayer::Match), (k, y3, AlnLayer::Delete)))
                    }
                }
                (Some(y1), Some(y2), Some(y3)) => {
                    if y1 + 1 >= y2 && y1 + 1 >= y3 {
                        Some(((k, y1 + 1, AlnLayer::Match), (k, y1, AlnLayer::Match)))
                    } else if y2 >= y3 {
                        Some(((k, y2, AlnLayer::Match), (k, y2, AlnLayer::Insert)))
                    } else {
                        Some(((k, y3, AlnLayer::Match), (k, y3, AlnLayer::Delete)))
                    }
                }
                (None, None, None) => None,
            };

            if let Some((connection_to, connection_from)) = connection {
                if !self.bt_contains(connection_to.0, connection_to.1, connection_to.2.clone()) {
                    let x = (connection_to.1 as i32 - connection_to.0) as usize;
                    let y = connection_to.1 as usize;
                    if x <= t_len && y <= q_len {
                        self.match_layer.set_y(score, k, connection_to.1);
                        self.bt_insert_if_absent(&connection_to, &connection_from, score);
                    }
                }
            }
        })
    }

    fn reduce(&mut self) {
        let (kmin, kmax) = self
            .match_layer
            .get_range(self.score)
            .expect("match-layer range missing at self.score");
        let k_end = self.query_str.len() as i32 - self.target_str.len() as i32;
        debug!("reduce, kmin-kmax: ({}):({}), {}", kmin, kmax, kmax - kmin);
        if (kmax - kmin) as u32 <= self.max_wf_length {
            return;
        }

        let t_len = self.target_str.len();
        let q_len = self.query_str.len();
        let k_offset = self.match_layer.k_offset;
        let k_width = self.match_layer.k_width;
        let score = self.score;
        // Convert a k to an in-bounds index into `kdist`, or None if k falls
        // outside the pre-sized window. `k_offset`/`k_width` were chosen so
        // every valid wavefront position fits, but the wavefront range can
        // temporarily include a 1-slot margin on either side.
        let kdist_idx = |k: i32| -> Option<usize> {
            let idx = k + k_offset;
            if idx < 0 || (idx as usize) >= k_width {
                None
            } else {
                Some(idx as usize)
            }
        };

        let mut dmin = usize::MAX;
        let mut any_entry = false;

        for k in kmin..kmax {
            let Some(y) = self.match_layer.get_y(score, k) else {
                continue;
            };
            let y = y as usize;
            if y > q_len {
                continue;
            }
            let dy = q_len - y;
            let x = (y as i32 - k) as usize;
            if x > t_len {
                continue;
            }
            let dx = t_len - x;
            let max_d = max(dx, dy);
            if let Some(k_idx) = kdist_idx(k) {
                self.kdist[k_idx] = max_d as u32;
                any_entry = true;
            }
            dmin = min(dmin, max_d);
        }

        if !any_entry {
            return;
        }

        let max_wf = self.max_wf_length as usize;

        let mut new_kmin = kmin;
        while new_kmin < kmax - 1 {
            let d = kdist_idx(new_kmin).map(|i| self.kdist[i]).unwrap_or(SENTINEL);
            if d == SENTINEL {
                new_kmin += 1;
                continue;
            }
            if (d as usize) - dmin <= max_wf {
                break;
            }
            new_kmin += 1;
        }

        let mut new_kmax = kmax;
        while new_kmax > new_kmin + 1 {
            let d = kdist_idx(new_kmax).map(|i| self.kdist[i]).unwrap_or(SENTINEL);
            if d == SENTINEL {
                new_kmax -= 1;
                continue;
            }
            if (d as usize) - dmin <= max_wf {
                new_kmax += 1;
                break;
            }
            new_kmax -= 1;
        }

        // Reset the scratch slots we wrote so the next reduce() sees SENTINEL
        // everywhere. We only touch (kmin..kmax), not the whole k_width.
        for k in kmin..kmax {
            if let Some(k_idx) = kdist_idx(k) {
                self.kdist[k_idx] = SENTINEL;
            }
        }

        // Never prune the target diagonal — the optimal path must end there.
        // Without this, reduce() can drift the range far from k_end and leave
        // the algorithm unable to close out the alignment.
        let new_kmin = min(new_kmin, k_end);
        let new_kmax = max(new_kmax, k_end + 1);

        debug!(
            "score: {}, kmin:{}, kmax:{}, new_kmin:{}, new_kmax:{}",
            score, kmin, kmax, new_kmin, new_kmax
        );

        self.match_layer.set_range(score, (new_kmin, new_kmax));

        let (i_kmin, i_kmax) = self.insertion_layer.get_range(score).unwrap();
        self.insertion_layer.set_range(
            score,
            (max(new_kmin, i_kmin), min(new_kmax, i_kmax)),
        );

        let (d_kmin, d_kmax) = self.deletion_layer.get_range(score).unwrap();
        self.deletion_layer.set_range(
            score,
            (max(new_kmin, d_kmin), min(new_kmax, d_kmax)),
        );
    }

    fn step_one(&mut self, max_score: Option<i32>) -> WaveFrontStepResult {
        let t_bytes = self.target_str.as_bytes();
        let q_bytes = self.query_str.as_bytes();
        let score = self.score;

        // Split-borrow: advance() needs &mut match_layer, the sink needs
        // &mut backtrace_map (disjoint fields).
        let match_layer = &mut self.match_layer;
        let backtrace_map = &mut self.backtrace_map;

        match_layer.advance(t_bytes, q_bytes, score, |(to, from), edge_score| {
            // Update-if-better: step_one's match extension should keep the
            // edge with the lowest score (matches the pre-0.2.1 semantics
            // that picked the earliest score at which a cell was reached).
            let entry = backtrace_map.entry(to).or_insert((from.clone(), edge_score));
            if edge_score < entry.1 {
                *entry = (from, edge_score);
            }
        });

        let k_end = self.query_str.len() as i32 - self.target_str.len() as i32;
        if let Some(y) = self.match_layer.get_y(self.score, k_end) {
            if y as usize >= self.query_str.len() {
                return WaveFrontStepResult::ReachEnd;
            }
        }

        self.score += 1;
        self.next(self.score);
        self.reduce();
        debug!("score: {}", self.score);
        debug!("---");
        if let Some(max_score) = max_score {
            if self.score >= max_score {
                WaveFrontStepResult::ReachMaxScore
            } else {
                WaveFrontStepResult::Continue
            }
        } else {
            WaveFrontStepResult::Continue
        }
    }

    pub fn step_all(&mut self, max_score: Option<i32>) -> WaveFrontStepResult {
        let result = loop {
            match self.step_one(max_score) {
                WaveFrontStepResult::Continue => {
                    continue;
                }
                WaveFrontStepResult::ReachMaxScore => break WaveFrontStepResult::ReachMaxScore,
                WaveFrontStepResult::Fail => break WaveFrontStepResult::Fail,
                WaveFrontStepResult::ReachEnd => break WaveFrontStepResult::ReachEnd,
            }
        };
        if result == WaveFrontStepResult::ReachEnd {
            let (t_aln, q_aln) = self.backtrace();
            self.score = gap_affine_score(
                &t_aln,
                &q_aln,
                self.mismatch_penalty,
                self.open_penalty,
                self.extension_penalty,
            );
        }
        result
    }

    pub fn backtrace(&mut self) -> (String, String) {
        let mut anchors = Vec::<Anchor>::new();

        let mut connect_to = (
            self.query_str.len() as i32 - self.target_str.len() as i32,
            self.query_str.len() as u32,
            AlnLayer::Match,
        );
        debug!(
            "seed: {:?} {} {}",
            connect_to,
            self.target_str.len(),
            self.query_str.len()
        );
        anchors.push(connect_to.clone());
        while let Some((pred_anchor, _score)) =
            self.bt_get(connect_to.0, connect_to.1, connect_to.2.clone())
        {
            assert!(pred_anchor != connect_to);
            connect_to = pred_anchor;
            anchors.push(connect_to.clone());
        }
        anchors.reverse();
        debug!("anchors length: {:?}", anchors.len());
        let t_str = self.target_str.to_string();
        let q_str = self.query_str.to_string();
        let mut t_aln_str = Vec::<&str>::new();
        let mut q_aln_str = Vec::<&str>::new();
        let mut previous_anchor: Option<Anchor> = None;
        anchors.into_iter().for_each(|(k, y, layer)| {
            if let Some((pk, py, pl)) = previous_anchor.clone() {
                let px = (py as i32 - pk) as usize;
                let py = py as usize;
                let x = (y as i32 - k) as usize;
                let y = y as usize;
                debug!(
                    "connection: {} {} {:?} -> {} {} {:?} ( {}, {} -> {}, {} )",
                    px, py, pl, x, y, layer, pk, py, k, y
                );
                if x > px || y > py {
                    if x > px {
                        t_aln_str.push(&t_str[px..x]);
                    } else {
                        t_aln_str.push("-");
                    }
                    if y > py {
                        q_aln_str.push(&q_str[py..y]);
                    } else {
                        q_aln_str.push("-");
                    }
                }
            };
            previous_anchor = Some((k, y, layer))
        });
        let t_aln_str = t_aln_str.join("");
        let q_aln_str = q_aln_str.join("");
        (t_aln_str, q_aln_str)
    }
}

#[cfg(test)]
mod tests {
    use super::*;

    fn score_alignment(t_aln: &str, q_aln: &str, mismatch: i32, open: i32, ext: i32) -> i32 {
        assert_eq!(
            t_aln.len(),
            q_aln.len(),
            "alignment strings must be equal length"
        );
        let tb = t_aln.as_bytes();
        let qb = q_aln.as_bytes();
        let mut score = 0i32;
        let mut in_gap_t = false;
        let mut in_gap_q = false;
        for i in 0..tb.len() {
            let t_gap = tb[i] == b'-';
            let q_gap = qb[i] == b'-';
            assert!(!(t_gap && q_gap), "all-gap column is invalid");
            if t_gap {
                score += if in_gap_t { ext } else { open + ext };
                in_gap_t = true;
                in_gap_q = false;
            } else if q_gap {
                score += if in_gap_q { ext } else { open + ext };
                in_gap_q = true;
                in_gap_t = false;
            } else {
                if tb[i] != qb[i] {
                    score += mismatch;
                }
                in_gap_t = false;
                in_gap_q = false;
            }
        }
        score
    }

    fn assert_valid_alignment(t_str: &str, q_str: &str, t_aln: &str, q_aln: &str) {
        assert_eq!(
            t_aln.replace('-', ""),
            t_str,
            "target round-trip identity failed"
        );
        assert_eq!(
            q_aln.replace('-', ""),
            q_str,
            "query round-trip identity failed"
        );
        assert_eq!(
            t_aln.len(),
            q_aln.len(),
            "alignment strings differ in length"
        );
    }

    fn run_alignment(
        t_str: &str,
        q_str: &str,
        max_wf_length: u32,
        mismatch: i32,
        open: i32,
        ext: i32,
        max_score: Option<i32>,
    ) -> (WaveFrontStepResult, Option<(String, String, i32)>) {
        let cap = std::cmp::max(16, std::cmp::max(t_str.len(), q_str.len()) >> 5);
        let mut wfs =
            WaveFronts::new_with_capacity(t_str, q_str, max_wf_length, mismatch, open, ext, cap);
        let result = wfs.step_all(max_score);
        let aln = if result == WaveFrontStepResult::ReachEnd {
            let score = wfs.score;
            let (t_aln, q_aln) = wfs.backtrace();
            Some((t_aln, q_aln, score))
        } else {
            None
        };
        (result, aln)
    }

    // Original regression fixtures, strengthened with invariant checks.

    #[test]
    fn test_step() {
        let t_str = "ACATACATGAAAAAAGTTGCATGAAACCCCAAAAGTTGCATGAAACATACATGAAAATACATGAAAGTTGCATGAAACATACATGAAAAAAGTTGCATGAAACCCCATACATGAAAGTTGCATGAA";
        let q_str = "ACATACATGAAAAAAGTTGCATGAAAAAACATACATGAAAGTTGCATGAAACATACATGAAAAAAGTTGCAAAAGTTGCATGAAACATACATGAAAATGAAAAAACATACATGAAAGTTGCATGAA";
        let (result, aln) = run_alignment(t_str, q_str, 40, 2, 2, 1, Some(1024));
        assert_eq!(result, WaveFrontStepResult::ReachEnd);
        let (t_aln, q_aln, score) = aln.unwrap();
        assert_valid_alignment(t_str, q_str, &t_aln, &q_aln);
        assert_eq!(score_alignment(&t_aln, &q_aln, 2, 2, 1), score);
    }

    #[test]
    fn test_step_2() {
        let t_str =
            "ACATACATGAAAAAAGTTGCATGAAACCCCAAAAGTTGCATGAAACATACATGAAAAAAAATGAAAGTTGCATGAAAATTTT";
        let q_str = "ACATACATGAAAAAAGTTGCATGAAACCCCAAAAGTTGCATGAAACATACATGAAAAATGAAAGTAAAATGAAAGTTGCATGAATGAAATGGTACATACATGAAAGTTGCAGGGG";
        let len_diff = (t_str.len() as i32 - q_str.len() as i32).unsigned_abs();
        let (result, aln) = run_alignment(t_str, q_str, len_diff, 9, 2, 1, Some(1024));
        assert_eq!(result, WaveFrontStepResult::ReachEnd);
        let (t_aln, q_aln, score) = aln.unwrap();
        assert_valid_alignment(t_str, q_str, &t_aln, &q_aln);
        assert_eq!(score_alignment(&t_aln, &q_aln, 9, 2, 1), score);
    }

    #[test]
    fn test_step_2_large_wf() {
        let t_str =
            "ACATACATGAAAAAAGTTGCATGAAACCCCAAAAGTTGCATGAAACATACATGAAAAAAAATGAAAGTTGCATGAAAATTTT";
        let q_str = "ACATACATGAAAAAAGTTGCATGAAACCCCAAAAGTTGCATGAAACATACATGAAAAATGAAAGTAAAATGAAAGTTGCATGAATGAAATGGTACATACATGAAAGTTGCAGGGG";
        let (result, aln) = run_alignment(t_str, q_str, 256, 9, 2, 1, Some(1024));
        assert_eq!(result, WaveFrontStepResult::ReachEnd);
        let (t_aln, q_aln, score) = aln.unwrap();
        assert_valid_alignment(t_str, q_str, &t_aln, &q_aln);
        assert_eq!(score_alignment(&t_aln, &q_aln, 9, 2, 1), score);
    }

    #[test]
    fn test_step_3() {
        let t_str = "GAAAGACCTGAAAGATCACGGTGCCTTCATTTCAACTGTGAGACATGAAGTAATTTTCCCAAATCTACAACATTAAGATATGGTGCAATAAGGACCAGAT";
        let q_str = "CTCCAACACGAGATTACCCAACCCAGGAGCAAGGAAATCAGTAACTTCCTCCCTATAACTTGGAATGTGGGTGGAGGGGTTCATAGTTCTCCCTGAGTGA";
        let len_diff = (t_str.len() as i32 - q_str.len() as i32).unsigned_abs();
        let max_wf_length = len_diff.max(1) * 2;
        let (result, aln) = run_alignment(t_str, q_str, max_wf_length, 4, 2, 1, Some(1024));
        assert_eq!(result, WaveFrontStepResult::ReachEnd);
        let (t_aln, q_aln, score) = aln.unwrap();
        assert_valid_alignment(t_str, q_str, &t_aln, &q_aln);
        assert_eq!(score_alignment(&t_aln, &q_aln, 4, 2, 1), score);
    }

    #[test]
    fn identical_sequences_score_zero() {
        let s = "ACGTACGTACGTACGTACGTACGT";
        let (result, aln) = run_alignment(s, s, 32, 4, 2, 1, Some(128));
        assert_eq!(result, WaveFrontStepResult::ReachEnd);
        let (t_aln, q_aln, score) = aln.unwrap();
        assert_eq!(score, 0);
        assert_eq!(t_aln, s);
        assert_eq!(q_aln, s);
    }

    #[test]
    fn single_char_match() {
        let (result, aln) = run_alignment("A", "A", 4, 4, 2, 1, Some(32));
        assert_eq!(result, WaveFrontStepResult::ReachEnd);
        let (t_aln, q_aln, score) = aln.unwrap();
        assert_eq!(score, 0);
        assert_eq!(t_aln, "A");
        assert_eq!(q_aln, "A");
    }

    #[test]
    fn single_char_mismatch() {
        let mismatch = 4;
        let (result, aln) = run_alignment("A", "T", 4, mismatch, 2, 1, Some(32));
        assert_eq!(result, WaveFrontStepResult::ReachEnd);
        let (t_aln, q_aln, score) = aln.unwrap();
        assert_valid_alignment("A", "T", &t_aln, &q_aln);
        assert_eq!(score, mismatch);
        assert_eq!(score_alignment(&t_aln, &q_aln, mismatch, 2, 1), score);
    }

    #[test]
    fn empty_target_single_char_query() {
        let open = 2;
        let ext = 1;
        let (result, aln) = run_alignment("", "A", 4, 4, open, ext, Some(32));
        assert_eq!(result, WaveFrontStepResult::ReachEnd);
        let (t_aln, q_aln, score) = aln.unwrap();
        assert_valid_alignment("", "A", &t_aln, &q_aln);
        assert_eq!(score, open + ext);
        assert_eq!(score_alignment(&t_aln, &q_aln, 4, open, ext), score);
    }

    #[test]
    fn empty_query_single_char_target() {
        let open = 2;
        let ext = 1;
        let (result, aln) = run_alignment("A", "", 4, 4, open, ext, Some(32));
        assert_eq!(result, WaveFrontStepResult::ReachEnd);
        let (t_aln, q_aln, score) = aln.unwrap();
        assert_valid_alignment("A", "", &t_aln, &q_aln);
        assert_eq!(score, open + ext);
        assert_eq!(score_alignment(&t_aln, &q_aln, 4, open, ext), score);
    }

    // Targets the reduce() off-by-one: a tight max_wf_length forces pruning to
    // fire. If the far edge diagonal is silently discarded, the run fails or
    // produces a wrong score. A single 4-base insertion should score open+4*ext.
    #[test]
    fn tight_max_wf_length_preserves_correctness() {
        let t_str = "ACGTACGTACGTACGTACGTACGT";
        let q_str = "ACGTACGTACGTACGTACGTACGTACGT";
        let len_diff = (q_str.len() - t_str.len()) as u32;
        let (result, aln) = run_alignment(t_str, q_str, len_diff, 4, 2, 1, Some(64));
        assert_eq!(result, WaveFrontStepResult::ReachEnd);
        let (t_aln, q_aln, score) = aln.unwrap();
        assert_valid_alignment(t_str, q_str, &t_aln, &q_aln);
        assert_eq!(score, 2 + 4);
        assert_eq!(score_alignment(&t_aln, &q_aln, 4, 2, 1), score);
    }

    #[test]
    fn reach_max_score_when_below_optimal() {
        let t_str = "AAAAAAAAAA";
        let q_str = "TTTTTTTTTT";
        let (result, aln) = run_alignment(t_str, q_str, 8, 4, 2, 1, Some(5));
        assert_eq!(result, WaveFrontStepResult::ReachMaxScore);
        assert!(aln.is_none());
    }

    #[test]
    fn repeated_runs_are_deterministic() {
        let t_str = "ACGTACGTACGTACGTACGT";
        let q_str = "ACGTACGTAAAACGTACGTACGT";
        let (_, a1) = run_alignment(t_str, q_str, 16, 4, 2, 1, Some(64));
        let (_, a2) = run_alignment(t_str, q_str, 16, 4, 2, 1, Some(64));
        let (t1, q1, s1) = a1.unwrap();
        let (t2, q2, s2) = a2.unwrap();
        assert_eq!(s1, s2);
        assert_eq!(t1, t2);
        assert_eq!(q1, q2);
    }

    #[test]
    fn asymmetric_length_with_clean_indel_large_wf() {
        let t_str = "ACGTACGTACGTACGTACGTACGTACGTACGT";
        let q_str = "ACGTACGTACGTACGTGGGGGGGGGGGGGGGGGGGGACGTACGTACGTACGT";
        let max_wf_length = 256;
        let (result, aln) = run_alignment(t_str, q_str, max_wf_length, 4, 2, 1, Some(256));
        assert_eq!(result, WaveFrontStepResult::ReachEnd);
        let (t_aln, q_aln, score) = aln.unwrap();
        assert_valid_alignment(t_str, q_str, &t_aln, &q_aln);
        assert_eq!(score_alignment(&t_aln, &q_aln, 4, 2, 1), score);
        assert_eq!(score, 2 + 20);
    }

    #[test]
    fn asymmetric_length_with_clean_indel_tight_wf() {
        let t_str = "ACGTACGTACGTACGTACGTACGTACGTACGT";
        let q_str = "ACGTACGTACGTACGTGGGGGGGGGGGGGGGGGGGGACGTACGTACGTACGT";
        let len_diff = (q_str.len() - t_str.len()) as u32;
        let max_wf_length = len_diff + 8;
        let (result, aln) = run_alignment(t_str, q_str, max_wf_length, 4, 2, 1, Some(256));
        assert_eq!(result, WaveFrontStepResult::ReachEnd);
        let (t_aln, q_aln, score) = aln.unwrap();
        assert_valid_alignment(t_str, q_str, &t_aln, &q_aln);
        assert_eq!(score_alignment(&t_aln, &q_aln, 4, 2, 1), score);
    }

}
