use ccm_rs::{ccm, simplex};

#[test]
fn ssr_kernel_returns_finite_rho_or_nan() {
    let a: Vec<f64> = (0..80).map(|i| (i as f64 / 10.0).sin()).collect();
    let out = simplex::ssr_pred_boot(&a, None, 2, 1, 1).expect("ssr should run");

    assert_eq!(out.aest.len(), a.len());
    assert!(!out.acceptablelib.is_empty());
    assert!(out.rho.is_finite() || out.rho.is_nan());
}

#[test]
fn ccm_kernel_returns_expected_shapes() {
    let a: Vec<f64> = (0..120).map(|i| (i as f64 / 8.0).sin()).collect();
    let b: Vec<f64> = (0..120)
        .map(|i| (i as f64 / 8.0).sin() + 0.15 * (i as f64 / 5.0).cos())
        .collect();

    let out = ccm::ccm_boot(&a, &b, 2, 1, None, 8).expect("ccm should run");

    assert_eq!(out.aest.len(), a.len());
    assert_eq!(out.rho.len(), out.sdevrho.len());
    assert!(!out.fullinfo.is_empty());
    assert_eq!(out.fullinfo.len(), out.rho.len());
    assert_eq!(out.fullinfo[0].len(), 8);
}
