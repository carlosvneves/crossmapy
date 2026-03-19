"""Compatibility shim for CCM module path."""

from .ccm.service import CCM_boot, ccmtest

__all__ = ["CCM_boot", "ccmtest"]
