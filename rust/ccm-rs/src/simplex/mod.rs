use crate::core::error::CcmError;

#[derive(Debug, Clone)]
pub struct SsrOutput {
    pub aest: Vec<f64>,
    pub rho: f64,
    pub acceptablelib: Vec<usize>,
    pub plengthacceptablelib: usize,
}

fn get_acceptable_lib(a: &[f64], e: usize, tau: usize, predstep: usize) -> Vec<usize> {
    let gapdist = tau * (e.saturating_sub(1)) + predstep;
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

    acceptable
        .iter()
        .enumerate()
        .filter_map(|(idx, v)| if *v > 0.0 { Some(idx) } else { None })
        .collect()
}

fn getorder_ssr(
    distances: &[f64],
    e: usize,
    acceptablelib: &[usize],
    i: usize,
    predstep: usize,
) -> (Vec<usize>, usize) {
    let mut nneigh = 1_usize;
    let mut n = 0_usize;
    let length = acceptablelib.len();

    if length == 0 {
        return (vec![0; e + 1], 0);
    }

    if acceptablelib[0] == i {
        n = 1;
    }
    if n >= length {
        n = length.saturating_sub(1);
    }

    let mut neighbors = vec![0_usize; e + 1];
    neighbors[0] = acceptablelib[n];

    if length <= predstep {
        return (neighbors, nneigh);
    }

    for iii in n..(length - predstep) {
        let ii = acceptablelib[iii];
        let mut trip = false;

        for j in 0..nneigh {
            if distances[ii] < distances[neighbors[j]] && ii != i && j > 0 {
                let mut k = nneigh;
                while k > j {
                    if k < e + 1 {
                        neighbors[k] = neighbors[k - 1];
                    }
                    k -= 1;
                }
                neighbors[j] = ii;
                trip = true;
                break;
            }
        }

        if !trip && nneigh < (e + 1) && ii != i && neighbors[nneigh - 1] != ii {
            neighbors[nneigh] = ii;
            if nneigh < e + 1 {
                nneigh += 1;
            }
        }
    }

    (neighbors, nneigh)
}

pub fn ssr_pred_boot(
    a_in: &[f64],
    b_in: Option<&[f64]>,
    e: usize,
    tau: usize,
    predstep: usize,
) -> Result<SsrOutput, CcmError> {
    if a_in.is_empty() {
        return Err(CcmError::InvalidInput("A must not be empty"));
    }

    let mut a = a_in.to_vec();
    let mut b = b_in.unwrap_or(a_in).to_vec();

    if b.len() != a.len() {
        return Err(CcmError::InvalidInput("A and B must have the same length"));
    }

    let finite_a_count = a.iter().filter(|x| x.is_finite()).count();
    let equal_finite_count = a
        .iter()
        .zip(b.iter())
        .filter(|(x, y)| x.is_finite() && y.is_finite() && (*x == *y))
        .count();

    let repvec = if equal_finite_count == finite_a_count && a.len() == b.len() {
        1
    } else {
        0
    };

    let acceptablelib = get_acceptable_lib(&a, e, tau, predstep);
    let lengthacceptablelib = acceptablelib.len();

    if tau.saturating_mul(e + 1).saturating_add(predstep) >= lengthacceptablelib {
        return Ok(SsrOutput {
            aest: vec![f64::NAN; a.len()],
            rho: f64::NAN,
            acceptablelib,
            plengthacceptablelib: lengthacceptablelib,
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

    let nneigh = e + 1;
    let mut aest = vec![0.0_f64; a.len()];
    let mut maxdist = 0.0_f64;

    for i in acceptablelib.iter().copied() {
        let mut distances = vec![f64::INFINITY; b.len()];

        if repvec == 1 {
            if lengthacceptablelib > predstep {
                for j in acceptablelib.iter().copied().take(lengthacceptablelib - predstep) {
                    if (j > i + predstep) || (j <= i.saturating_sub(e)) {
                        let mut dist = 0.0;
                        for k in 0..e {
                            let ai = i.saturating_sub(tau * k);
                            let bj = j.saturating_sub(tau * k);
                            dist += (a[ai] - b[bj]).powi(2);
                        }
                        distances[j] = dist.sqrt();
                        if distances[j] > maxdist {
                            maxdist = 999_999_999.0 * distances[j];
                        }
                    } else {
                        distances[j] = maxdist;
                    }
                }
            }
        } else {
            for j in acceptablelib.iter().copied() {
                let mut dist = 0.0;
                for k in 0..e {
                    let ai = i.saturating_sub(tau * k);
                    let bj = j.saturating_sub(tau * k);
                    dist += (a[ai] - b[bj]).powi(2);
                }
                distances[j] = dist.sqrt();
            }
        }

        let (neighbors, found_nneigh) = getorder_ssr(&distances, e, &acceptablelib, i, predstep);
        if found_nneigh < nneigh {
            aest[i] = 0.0;
            continue;
        }

        let distsv = distances[neighbors[0]];
        let mut sumaest = 0.0;

        if distsv != 0.0 {
            let mut u = vec![0.0_f64; nneigh];
            for j in 0..nneigh {
                u[j] = (-distances[neighbors[j]] / distsv).exp();
            }
            let sumu: f64 = u.iter().sum();

            let mut w = vec![0.0_f64; nneigh];
            for j in 0..nneigh {
                w[j] = (u[j] / sumu).max(0.000001);
            }
            let sumw: f64 = w.iter().sum();

            for j in 0..nneigh {
                w[j] /= sumw;
                let idx = (neighbors[j] + predstep).min(b.len() - 1);
                sumaest += b[idx] * w[j];
            }
        } else {
            let mut w = vec![0.0_f64; nneigh];
            let mut sumw = 0.0;
            for j in 0..nneigh {
                w[j] = if distances[neighbors[j]] == 0.0 { 1.0 } else { 0.000001 };
                sumw += w[j];
            }
            for j in 0..nneigh {
                w[j] /= sumw;
                sumaest += a[neighbors[j]] * w[j];
            }
        }

        aest[i] = sumaest;
    }

    let mut aest_out = vec![0.0_f64; a.len()];
    if predstep < a.len() {
        aest_out[predstep..a.len()].copy_from_slice(&aest[..(a.len() - predstep)]);
    }
    for v in &mut aest_out {
        if *v == 0.0 {
            *v = f64::NAN;
        }
    }

    let mut xv = Vec::new();
    let mut yv = Vec::new();
    for idx in 0..a.len() {
        if a[idx].is_finite() && aest_out[idx].is_finite() {
            xv.push(a[idx]);
            yv.push(aest_out[idx]);
        }
    }

    let rho = if xv.len() > 1 {
        let xbar = xv.iter().sum::<f64>() / xv.len() as f64;
        let ybar = yv.iter().sum::<f64>() / yv.len() as f64;
        let mut xyxybar = 0.0;
        let mut xxbarsq = 0.0;
        let mut yybarsq = 0.0;
        for i in 0..xv.len() {
            xyxybar += (xv[i] - xbar) * (yv[i] - ybar);
            xxbarsq += (xv[i] - xbar).powi(2);
            yybarsq += (yv[i] - ybar).powi(2);
        }
        let denom = xxbarsq.sqrt() * yybarsq.sqrt();
        if denom == 0.0 {
            0.0
        } else {
            xyxybar / denom
        }
    } else {
        f64::NAN
    };

    Ok(SsrOutput {
        aest: aest_out,
        rho,
        acceptablelib,
        plengthacceptablelib: lengthacceptablelib,
    })
}
