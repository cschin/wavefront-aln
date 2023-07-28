use log::debug;
use rustc_hash::FxHashMap;

use std::cmp::{max, min};

#[derive(Default)]
struct WaveFront {
    s_k_to_y_map: FxHashMap<(i32, i32), u32>,
    score_to_k_range: FxHashMap<i32, (i32, i32)>,
}
#[derive(PartialEq, Eq, Hash, Clone, Debug)]
pub enum AlnLayer {
    Insert,
    Delete,
    Match,
}

pub struct WaveFronts<'a> {
    target_str: &'a str,
    query_str: &'a str,
    insertion_layer: WaveFront,
    deletion_layer: WaveFront,
    match_layer: WaveFront,
    pub min_wf_length: u32,
    pub mismatch_penalty: i32,
    pub open_penalty: i32,
    pub extension_penalty: i32,
    pub score: i32,
    pub backtrace_map: FxHashMap<(i32, u32, AlnLayer), ((i32, u32, AlnLayer), i32)>,
}
type Anchor = (i32, u32, AlnLayer);
type Connection = (Anchor, Anchor); // Connection.0 is the best score at the end point

impl WaveFront {
    fn new() -> Self {
        let mut wf = WaveFront::default();
        wf.score_to_k_range.insert(0, (0, 1));
        wf
    }

    fn new_with_capacity(capacity: usize) -> Self {
        let s_k_to_y_map = FxHashMap::with_capacity_and_hasher(capacity * 64, Default::default());
        let score_to_k_range = FxHashMap::with_capacity_and_hasher(capacity, Default::default());
        let mut wf = WaveFront {
            s_k_to_y_map,
            score_to_k_range
        };
        wf.score_to_k_range.insert(0, (0, 1));
        wf
    }

    fn advance<'a>(&mut self, t: &'a str, q: &'a str, score: i32) -> Vec<(Connection, i32)> {
        debug!("advance score: {}", score);
        let q = q.as_bytes();
        let t = t.as_bytes();
        let &(k_min, k_max) = self.score_to_k_range.get(&score).expect("no range");
        (k_min..k_max)
            .filter_map(|k| {
                if let Some(&y) = self.s_k_to_y_map.get(&(score, k)) {
                    let mut xs = (y as i32 - k) as usize; // k = y - x
                    let mut ys = y as usize;
                    loop {
                        if xs >= t.len() || ys >= q.len() || t[xs] != q[ys] {
                            break;
                        }
                        xs += 1;
                        ys += 1;
                    }
                    debug!("advance: k:{} y: {} -> {}", k, y, ys);
                    if ys > y as usize {
                        self.s_k_to_y_map.insert((score, k), ys as u32);

                        Some((
                            ((k, ys as u32, AlnLayer::Match), (k, y, AlnLayer::Match)),
                            score,
                        ))
                    } else {
                        None
                    }
                } else {
                    None
                }
            })
            .collect::<Vec<_>>()
    }
}

impl<'a> WaveFronts<'a> {
    pub fn new(
        t_str: &'a str,
        q_str: &'a str,
        min_wf_length: u32,
        mismatch_penalty: i32,
        open_penalty: i32,
        extension_penalty: i32,
    ) -> Self {
        let insertion_layer = WaveFront::new();
        let deletion_layer = WaveFront::new();
        let mut match_layer = WaveFront::new();
        match_layer.s_k_to_y_map.insert((0, 0), 0);

        WaveFronts {
            target_str: t_str,
            query_str: q_str,
            insertion_layer,
            deletion_layer,
            match_layer,
            min_wf_length,
            mismatch_penalty,
            open_penalty,
            extension_penalty,
            score: 0,
            backtrace_map: FxHashMap::default(),
        }
    }

    pub fn new_with_capacity(
        t_str: &'a str,
        q_str: &'a str,
        min_wf_length: u32,
        mismatch_penalty: i32,
        open_penalty: i32,
        extension_penalty: i32,
        capacity: usize
    ) -> Self {
        let insertion_layer = WaveFront::new_with_capacity(capacity);
        let deletion_layer = WaveFront::new_with_capacity(capacity);
        let mut match_layer = WaveFront::new_with_capacity(capacity);
        match_layer.s_k_to_y_map.insert((0, 0), 0);

        WaveFronts {
            target_str: t_str,
            query_str: q_str,
            insertion_layer,
            deletion_layer,
            match_layer,
            min_wf_length,
            mismatch_penalty,
            open_penalty,
            extension_penalty,
            score: 0,
            backtrace_map: FxHashMap::default(),
        }
    }

    pub fn next(&mut self, score: i32) {
        let (match_k_min, match_k_max) = *self
            .match_layer
            .score_to_k_range
            .get(&(score - self.mismatch_penalty))
            .unwrap_or(&(0, 1));

        let (open_k_min, open_k_max) = *self
            .match_layer
            .score_to_k_range
            .get(&(score - self.open_penalty - self.extension_penalty))
            .unwrap_or(&(0, 1));

        let (i_k_min, i_k_max) = *self
            .insertion_layer
            .score_to_k_range
            .get(&(score - self.extension_penalty))
            .unwrap_or(&(0, 1));

        let (d_k_min, d_k_max) = *self
            .deletion_layer
            .score_to_k_range
            .get(&(score - self.extension_penalty))
            .unwrap_or(&(0, 1));

        let k_max = max(match_k_max, max(open_k_max, max(i_k_max, d_k_max))) + 1;
        let k_min = min(match_k_min, min(open_k_min, min(i_k_min, d_k_min))) - 1;
        debug!("score:{} k_min:{}, k_max:{}", score, k_min, k_max);
        self.insertion_layer
            .score_to_k_range
            .insert(score, (k_min, k_max));
        self.deletion_layer
            .score_to_k_range
            .insert(score, (k_min, k_max));
        self.match_layer
            .score_to_k_range
            .insert(score, (k_min, k_max));

        (k_min..k_max).for_each(|k| {
            debug!("k: {:?}", k);

            let e1 = self
                .match_layer
                .s_k_to_y_map
                .get(&(score - self.open_penalty - self.extension_penalty, k - 1));

            let e2 = self
                .insertion_layer
                .s_k_to_y_map
                .get(&(score - self.extension_penalty, k - 1));
            debug!(
                "I {:?} {:?} {:?}",
                e1, e2, &self.insertion_layer.s_k_to_y_map
            );

            let connection = match (e1, e2) {
                (Some(y1), None) => Some((
                    score,
                    (k, *y1 + 1, AlnLayer::Insert),
                    (k - 1, *y1, AlnLayer::Match),
                )),
                (None, Some(y2)) => Some((
                    score,
                    (k, *y2 + 1, AlnLayer::Insert),
                    (k - 1, *y2, AlnLayer::Insert),
                )),
                (Some(y1), Some(y2)) => {
                    if *y1 >= *y2 {
                        Some((
                            score,
                            (k, *y1 + 1, AlnLayer::Insert),
                            (k - 1, *y1, AlnLayer::Match),
                        ))
                    } else {
                        Some((
                            score,
                            (k, *y2 + 1, AlnLayer::Insert),
                            (k - 1, *y2, AlnLayer::Insert),
                        ))
                    }
                }
                (None, None) => None,
            };

            if let Some((score, connection_to, connection_from)) = connection {
                let x = (connection_to.1 as i32 - connection_to.0) as usize;
                let y = connection_to.1 as usize;
                if x < self.target_str.len() && y < self.query_str.len() {
                    self.insertion_layer
                        .s_k_to_y_map
                        .insert((score, k), connection_to.1);
                    let e = self
                        .backtrace_map
                        .entry(connection_to.clone())
                        .or_insert((connection_from.clone(), score));
                    if e.1 > score {
                        assert!(connection_to != connection_from);
                        self.backtrace_map
                            .insert(connection_to, (connection_from, score));
                    }
                }
            };

            let e1 = self
                .match_layer
                .s_k_to_y_map
                .get(&(score - self.open_penalty - self.extension_penalty, k + 1));

            let e2 = self
                .deletion_layer
                .s_k_to_y_map
                .get(&(score - self.extension_penalty, k + 1));
            debug!(
                "D {:?} {:?} {:?}",
                e1, e2, &self.deletion_layer.s_k_to_y_map
            );
            let connection = match (e1, e2) {
                (Some(y1), None) => {
                    Some(((k, *y1, AlnLayer::Delete), (k + 1, *y1, AlnLayer::Match)))
                }
                (None, Some(y2)) => {
                    Some(((k, *y2, AlnLayer::Delete), (k + 1, *y2, AlnLayer::Delete)))
                }
                (Some(y1), Some(y2)) => {
                    if *y1 >= *y2 {
                        Some(((k, *y1, AlnLayer::Delete), (k + 1, *y1, AlnLayer::Match)))
                    } else {
                        Some(((k, *y2, AlnLayer::Delete), (k + 1, *y2, AlnLayer::Delete)))
                    }
                }
                (None, None) => None,
            };

            if let Some((connection_to, connection_from)) = connection {
                let x = (connection_to.1 as i32 - connection_to.0) as usize;
                let y = connection_to.1 as usize;
                if x < self.target_str.len() && y < self.query_str.len() {
                    self.deletion_layer
                        .s_k_to_y_map
                        .insert((score, k), connection_to.1);
                    let e = self
                        .backtrace_map
                        .entry(connection_to.clone())
                        .or_insert((connection_from.clone(), score));
                    if e.1 > score {
                        assert!(connection_to != connection_from);
                        self.backtrace_map
                            .insert(connection_to, (connection_from, score));
                    }
                }
            };

            let e1 = self
                .match_layer
                .s_k_to_y_map
                .get(&(score - self.mismatch_penalty, k));
            let e2 = self.insertion_layer.s_k_to_y_map.get(&(score, k));
            let e3 = self.deletion_layer.s_k_to_y_map.get(&(score, k));
            debug!(
                "M {:?} {:?} {:?} {:?}",
                e1, e2, e3, &self.match_layer.s_k_to_y_map
            );
            let connection = match (e1, e2, e3) {
                (Some(y1), None, None) => {
                    Some(((k, *y1 + 1, AlnLayer::Match), (k, *y1, AlnLayer::Match)))
                }
                (None, Some(y2), None) => {
                    Some(((k, *y2, AlnLayer::Match), (k, *y2, AlnLayer::Insert)))
                }
                (None, None, Some(y3)) => {
                    Some(((k, *y3, AlnLayer::Match), (k, *y3, AlnLayer::Delete)))
                }

                (Some(y1), Some(y2), None) => {
                    if *y1 + 1 >= *y2 {
                        Some(((k, *y1 + 1, AlnLayer::Match), (k, *y1, AlnLayer::Match)))
                    } else {
                        Some(((k, *y2, AlnLayer::Match), (k, *y2, AlnLayer::Insert)))
                    }
                }
                (Some(y1), None, Some(y3)) => {
                    if *y1 + 1 >= *y3 {
                        Some(((k, *y1 + 1, AlnLayer::Match), (k, *y1, AlnLayer::Match)))
                    } else {
                        Some(((k, *y3, AlnLayer::Match), (k, *y3, AlnLayer::Delete)))
                    }
                }
                (None, Some(y2), Some(y3)) => {
                    if *y2 >= *y3 {
                        Some(((k, *y2, AlnLayer::Match), (k, *y2, AlnLayer::Insert)))
                    } else {
                        Some(((k, *y3, AlnLayer::Match), (k, *y3, AlnLayer::Delete)))
                    }
                }

                (Some(y1), Some(y2), Some(y3)) => {
                    if *y1 + 1 >= *y2 && *y1 + 1 >= *y3 {
                        Some(((k, *y1 + 1, AlnLayer::Match), (k, *y1, AlnLayer::Match)))
                    } else if *y2 >= *y3 {
                        Some(((k, *y2, AlnLayer::Match), (k, *y2, AlnLayer::Insert)))
                    } else {
                        Some(((k, *y3, AlnLayer::Match), (k, *y3, AlnLayer::Delete)))
                    }
                }

                (None, None, None) => None,
            };

            if let Some((connection_to, connection_from)) = connection {
                if self.backtrace_map.get(&connection_to).is_none() {
                    let x = (connection_to.1 as i32 - connection_to.0) as usize;
                    let y = connection_to.1 as usize;
                    if x < self.target_str.len() && y < self.query_str.len() {
                        self.match_layer
                            .s_k_to_y_map
                            .insert((score, k), connection_to.1);
                        let e = self
                            .backtrace_map
                            .entry(connection_to.clone())
                            .or_insert((connection_from.clone(), score));

                        if e.1 > score {
                            assert!(connection_to != connection_from);
                            self.backtrace_map
                                .insert(connection_to, (connection_from, score));
                        }
                    }
                }
            };
        })
    }

    fn reduce(&mut self) {
        let (kmin, kmax) = *self.match_layer.score_to_k_range.get(&self.score).unwrap();
        debug!("reduce, kmin-kmax: ({}):({}), {}", kmin, kmax, kmax - kmin);
        if (kmax - kmin) as u32 > self.min_wf_length {
            let mut dmin = usize::MAX;
            let mut kdist = FxHashMap::<i32, usize>::default();

            (kmin..kmax).for_each(|k| {
                debug!("score: {}, k: {}", self.score, k);
                if !self.match_layer.s_k_to_y_map.contains_key(&(self.score, k)) {
                    return;
                }
                let y = *self.match_layer.s_k_to_y_map.get(&(self.score, k)).unwrap() as usize;
                if y > self.query_str.len() {
                    return;
                }
                debug!("query: {} y: {}\n", self.query_str.len(), y);
                let dy = self.query_str.len() - y;
                let x = (y as i32 - k) as usize;
                if x > self.target_str.len() {
                    return;
                }
                let dx = self.target_str.len() - x;
                let max_d = max(dx, dy);
                kdist.insert(k, max_d);
                dmin = min(dmin, max_d);
            });

            let mut new_kmin = kmin;
            while new_kmin < kmax {
                if !kdist.contains_key(&new_kmin) {
                    new_kmin += 1;
                    continue;
                }
                if *kdist.get(&new_kmin).unwrap() - dmin <= self.min_wf_length as usize {
                    break;
                }
                new_kmin += 1;
            }

            let mut new_kmax = kmax;
            while new_kmax > new_kmin {
                if !kdist.contains_key(&new_kmax) {
                    new_kmax -= 1;
                    continue;
                }
                if *kdist.get(&new_kmax).unwrap() - dmin <= self.min_wf_length as usize {
                    break;
                }
                new_kmax -= 1;
            }
            debug!(
                "score: {}, kmin:{}, kmax:{}, new_kmin:{}, new_kmax:{}",
                self.score, kmin, kmax, new_kmin, new_kmax
            );

            *self
                .match_layer
                .score_to_k_range
                .get_mut(&self.score)
                .unwrap() = (new_kmin, new_kmax);

            let (i_kmin, i_kmax) = *self
                .insertion_layer
                .score_to_k_range
                .get(&self.score)
                .unwrap();

            *self
                .insertion_layer
                .score_to_k_range
                .get_mut(&self.score)
                .unwrap() = (max(new_kmin, i_kmin), min(new_kmax, i_kmax));

            let (d_kmin, d_kmax) = *self
                .deletion_layer
                .score_to_k_range
                .get(&self.score)
                .unwrap();

            *self
                .deletion_layer
                .score_to_k_range
                .get_mut(&self.score)
                .unwrap() = (max(new_kmin, d_kmin), min(new_kmax, d_kmax));
        }
    }

    fn step_one(&mut self) -> Option<i32> {
        self.match_layer
            .advance(self.target_str, self.query_str, self.score)
            .into_iter()
            .for_each(|((connection_to, connection_from), score)| {
                if let Some(e) = self.backtrace_map.get(&connection_to) {
                    if score < e.1 {
                        self.backtrace_map
                            .insert(connection_to, (connection_from, score));
                    }
                } else {
                    self.backtrace_map
                        .insert(connection_to, (connection_from, score));
                }
            });
        debug!("backtrack_map: {:?}", self.backtrace_map.clone());
        debug!("match wf 0: {:?}", self.match_layer.s_k_to_y_map.clone());
        if let Some(y) = self.match_layer.s_k_to_y_map.get(&(
            self.score,
            self.query_str.len() as i32 - self.target_str.len() as i32,
        )) {
            if *y as usize >= self.query_str.len() {
                return None;
            }
        }

        self.score += 1;
        self.next(self.score);
        self.reduce();
        debug!("score: {}", self.score);
        debug!("match wf 1: {:?}", self.match_layer.s_k_to_y_map.clone());
        debug!("---");
        Some(self.score)
    }

    pub fn step_all(&mut self) {
        loop {
            if self.step_one().is_none() {
                break;
            }
        }
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
        while let Some(connect_from) = self.backtrace_map.get(&connect_to) {
            assert!(connect_from.0 != connect_to);
            connect_to = connect_from.clone().0;
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

    #[test]
    fn test_step() {
        use simple_logger::SimpleLogger;
        SimpleLogger::new().init().unwrap();
        let t_str = "ACATACATGAAAAAAGTTGCATGAAACCCCAAAAGTTGCATGAAACATACATGAAAATACATGAAAGTTGCATGAAACATACATGAAAAAAGTTGCATGAAACCCCATACATGAAAGTTGCATGAA";
        let q_str = "ACATACATGAAAAAAGTTGCATGAAAAAACATACATGAAAGTTGCATGAAACATACATGAAAAAAGTTGCAAAAGTTGCATGAAACATACATGAAAATGAAAAAACATACATGAAAGTTGCATGAA";
        let mut wfs = WaveFronts::new_with_capacity(t_str, q_str, 40, 2, 2, 1, t_str.len() >> 4);
        wfs.step_all();
        let (t_aln_str, q_aln_str) = wfs.backtrace();
        println!("{}", t_aln_str);
        println!("{}", q_aln_str);
    }

    #[test]
    fn test_step_2() {
        use simple_logger::SimpleLogger;
        SimpleLogger::new().init().unwrap();
        let t_str = "ACATACATGAAAAAAGTTGCATGAAACCCCAAAAGTTGCATGAAACATACATGAAAAAAAATGAAAGTTGCATGAA";
        let q_str = "ACATACATGAAAAAAGTTGCATGAAACCCCAAAAGTTGCATGAAACATACATGAAAAATGAAAGTAAAATGAAAGTTGCATGAATGCATACATGAAAGTTGCA";
        let mut wfs = WaveFronts::new_with_capacity(t_str, q_str, 40, 2, 2, 1, t_str.len() >> 4);
        wfs.step_all();
        let (t_aln_str, q_aln_str) = wfs.backtrace();
        println!("{}", t_aln_str);
        println!("{}", q_aln_str);
    }
}
