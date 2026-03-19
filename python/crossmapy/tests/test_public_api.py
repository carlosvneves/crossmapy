from crossmapy import CCM_boot, SSR_pred_boot, ccmtest, make_ccm_data


def test_imports_are_available():
    assert callable(make_ccm_data)
    assert callable(SSR_pred_boot)
    assert callable(CCM_boot)
    assert callable(ccmtest)
