use crate::core::error::CcmError;
use rand::prelude::*;
use rand_chacha::ChaCha8Rng;
use std::collections::BTreeSet;

#[derive(Debug, Clone)]
pub struct CcmOutput {
    pub aest: Vec<f64>,
    pub rho: Vec<f64>,
    pub sdevrho: Vec<f64>,
    pub lobs: Vec<usize>,
    pub fullinfo: Vec<Vec<f64>>,
}

fn get_acceptable_lib_ccm(a: &[f64], e: usize, tau: usize, plengtht: usize) -> (Vec<usize>, Vec<usize>) {
    let gapdist = tau * e.saturating_sub(1);
    let mut acceptable = vec![1.0_f64; a.len()];

    for i in 0..a.len() {
        if !a[i].is_finite() {
            acceptable[i] = 0.0;
        }
    }

    for i in 1..=gapdist {
        for idx in (0..a.len()).rev() {
            let shifted_finite = if idx >= i { a[idx - i].is_finite() } else { false };
            if !shifted_finite {
                acceptable[idx] = 0.0;
            }
        }
    }

    let acceptablelib: Vec<usize> = acceptable
        .iter()
        .enumerate()
        .filter_map(|(idx, v)| if *v > 0.0 && idx <= (plengtht.saturating_sub(1)) { Some(idx) } else { None })
        .collect();

    let acceptablelib2: Vec<usize> = acceptablelib
        .iter()
        .copied()
        .filter(|x| *x < (plengtht.saturating_sub(1)).saturating_sub(tau))
        .collect();

    (acceptablelib, acceptablelib2)
}

fn get_rho(a: &[f64], aest: &[f64], acceptablelib: &[usize]) -> f64 {
    if acceptablelib.is_empty() {
        return 0.0;
    }

    let xbar = acceptablelib.iter().map(|i| a[*i]).sum::<f64>() / acceptablelib.len() as f64;
    let ybar = acceptablelib.iter().map(|i| aest[*i]).sum::<f64>() / acceptablelib.len() as f64;

    let mut xyxybar = 0.0;
    let mut xxbarsq = 0.0;
    let mut yybarsq = 0.0;

    for idx in acceptablelib {
        xyxybar += (a[*idx] - xbar) * (aest[*idx] - ybar);
        xxbarsq += (a[*idx] - xbar).powi(2);
        yybarsq += (aest[*idx] - ybar).powi(2);
    }

    let denom = xxbarsq.sqrt() * yybarsq.sqrt();
    if denom == 0.0 {
        return 0.0;
    }

    let rhocalc = xyxybar / denom;
    if (-1.0..=1.0).contains(&rhocalc) {
        rhocalc
    } else {
        0.0
    }
}

fn nearest_value(sorted: &[usize], target: usize) -> Option<usize> {
    if sorted.is_empty() {
        return None;
    }
    let mut best = sorted[0];
    let mut best_diff = best.abs_diff(target);
    for v in sorted.iter().copied() {
        let diff = v.abs_diff(target);
        if diff < best_diff {
            best_diff = diff;
            best = v;
        }
    }
    Some(best)
}

pub fn ccm_boot(
    a_in: &[f64],
    b_in: &[f64],
    e: usize,
    tau: usize,
    desired_l_in: Option<&[usize]>,
    iterations: usize,
) -> Result<CcmOutput, CcmError> {
    if a_in.is_empty() || b_in.is_empty() {
        return Err(CcmError::InvalidInput("A and B must not be empty"));
    }
    if a_in.len() != b_in.len() {
        return Err(CcmError::InvalidInput("A and B must have same length"));
    }

    let mut a = a_in.to_vec();
    let mut b = b_in.to_vec();
    let length_a = a.len();
    let plengtht = a.iter().filter(|x| x.is_finite()).count().min(a.len());

    let (acceptablelib, acceptablelib2) = get_acceptable_lib_ccm(&a, e, tau, plengtht);
    let lengthacceptablelib = acceptablelib.len();
    let from_idx = tau * e.saturating_sub(1);

    let mut desired_l: Vec<usize> = if let Some(values) = desired_l_in {
        values.iter().map(|x| x + e.saturating_sub(2)).collect()
    } else {
        let start = from_idx + e + 1;
        let end_excl = length_a.saturating_sub(e) + 2;
        (start..end_excl).collect()
    };

    let mut valid_l = Vec::new();
    for dl in desired_l.drain(..) {
        if let Some(v) = nearest_value(&acceptablelib2, dl) {
            valid_l.push(v);
        }
    }
    valid_l.sort_unstable();
    valid_l.dedup();
    desired_l = valid_l;

    if tau.saturating_mul(e + 1) > lengthacceptablelib {
        return Ok(CcmOutput {
            aest: vec![f64::NAN; a.len()],
            rho: vec![f64::NAN],
            sdevrho: vec![f64::NAN],
            lobs: vec![usize::MAX],
            fullinfo: vec![vec![f64::NAN]],
        });
    }

    for v in &mut a {
        if !v.is_finite() {
            *v = 0.0;
        }
    }
    for v in &mut b {
        if !v.is_finite() {
            *v = 0.0;
        }
    }

    let n_neighbors = e + 1;
    let min_weight = 0.000001_f64;
    let mut rng = ChaCha8Rng::seed_from_u64(42);

    let mut lpos_set = BTreeSet::new();
    let mut rho_results: Vec<Vec<f64>> = Vec::new();
    let mut aest_results: Vec<Vec<f64>> = Vec::new();

    for _ in 0..iterations {
        let mut rho_iter = vec![0.0_f64; desired_l.len()];
        let mut aest = vec![0.0_f64; a.len()];

        for (lidx, l0) in desired_l.iter().copied().enumerate() {
            let mut l = l0;
            if l < (from_idx + e + 1) {
                l = from_idx + e + 1;
            }
            if l >= plengtht {
                l = plengtht.saturating_sub(1);
            }
            let to = l;
            if to < from_idx {
                continue;
            }

            let lib_size = to - from_idx + 1;
            let mut lib_indices = vec![0_usize; lib_size];
            for v in &mut lib_indices {
                let idx = rng.random_range(0..acceptablelib.len());
                *v = acceptablelib[idx];
            }

            aest.fill(0.0);

            for i in acceptablelib.iter().copied() {
                if i < from_idx || i >= b.len() {
                    continue;
                }

                let mut distances = vec![0.0_f64; lib_size];
                for (row, lib_idx) in lib_indices.iter().copied().enumerate() {
                    let mut dist = 0.0;
                    for k in 0..e {
                        let p1 = i.saturating_sub(tau * k);
                        let p2 = lib_idx.saturating_sub(tau * k);
                        dist += (b[p1] - b[p2]).powi(2);
                    }
                    distances[row] = dist.sqrt();
                    if lib_idx == i {
                        distances[row] = f64::INFINITY;
                    }
                }

                let mut ordered: Vec<(usize, f64)> = distances.iter().copied().enumerate().collect();
                ordered.sort_by(|a, b| a.1.partial_cmp(&b.1).unwrap_or(std::cmp::Ordering::Equal));
                let take_n = n_neighbors.min(ordered.len());
                if take_n == 0 {
                    continue;
                }

                let neighbor_rows: Vec<usize> = ordered.iter().take(take_n).map(|(idx, _)| *idx).collect();
                let neighbor_dists: Vec<f64> = neighbor_rows.iter().map(|idx| distances[*idx]).collect();
                let neighbors: Vec<usize> = neighbor_rows.iter().map(|idx| lib_indices[*idx]).collect();

                let distsv = neighbor_dists[0];
                if distsv != 0.0 && distsv.is_finite() {
                    let mut u = vec![0.0_f64; take_n];
                    for j in 0..take_n {
                        u[j] = (-neighbor_dists[j] / distsv).exp();
                    }
                    let sumu: f64 = u.iter().sum();
                    let mut w = vec![0.0_f64; take_n];
                    for j in 0..take_n {
                        w[j] = (u[j] / sumu).max(min_weight);
                    }
                    let sumw: f64 = w.iter().sum();
                    for j in 0..take_n {
                        w[j] /= sumw;
                        aest[i] += a[neighbors[j]] * w[j];
                    }
                } else {
                    let mut w = vec![min_weight; take_n];
                    for j in 0..take_n {
                        if neighbor_dists[j] == 0.0 {
                            w[j] = 1.0;
                        }
                    }
                    let sumw: f64 = w.iter().sum();
                    for j in 0..take_n {
                        w[j] /= sumw;
                        aest[i] += a[neighbors[j]] * w[j];
                    }
                }
            }

            rho_iter[lidx] = get_rho(&a, &aest, &acceptablelib);
        }

        for lob in desired_l.iter().map(|v| v.saturating_sub(e).saturating_add(1)) {
            lpos_set.insert(lob);
        }

        rho_results.push(rho_iter);
        aest_results.push(aest);
    }

    let lobs: Vec<usize> = lpos_set.into_iter().collect();

    let rows = desired_l.len();
    let cols = iterations;
    let mut rho_mat = vec![vec![f64::NAN; cols]; rows];
    for c in 0..cols {
        for r in 0..rows {
            rho_mat[r][c] = rho_results[c][r];
        }
    }

    let mut rho_means = vec![0.0_f64; rows];
    let mut rho_sdev = vec![0.0_f64; rows];
    for r in 0..rows {
        let row = &rho_mat[r];
        let mean = row.iter().sum::<f64>() / row.len() as f64;
        rho_means[r] = mean;
        let var = row.iter().map(|x| (x - mean).powi(2)).sum::<f64>() / row.len() as f64;
        rho_sdev[r] = var.sqrt();
    }

    let mut aest_avg = vec![0.0_f64; a.len()];
    for v in &aest_results {
        for (idx, x) in v.iter().enumerate() {
            aest_avg[idx] += *x;
        }
    }
    if !aest_results.is_empty() {
        for x in &mut aest_avg {
            *x /= aest_results.len() as f64;
            if *x == 0.0 {
                *x = f64::NAN;
            }
        }
    }

    Ok(CcmOutput {
        aest: aest_avg,
        rho: rho_means,
        sdevrho: rho_sdev,
        lobs,
        fullinfo: rho_mat,
    })
}
