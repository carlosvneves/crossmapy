# crossmapy

Monorepo with a two-layer architecture:

- `python/crossmapy`: public Python API and wrapper
- `rust/ccm-rs`: computational core for CCM kernels

The current baseline copies the existing multispatialCCM computational logic into the new structure and keeps a safe Python fallback while Rust kernels are wired.
