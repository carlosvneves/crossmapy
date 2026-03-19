# crossmapy

Monorepo with a two-layer architecture:

- `python/crossmapy`: public Python API and wrapper
- `rust/ccm-rs`: computational core for CCM kernels

Rust kernels are now implemented in `ccm-rs` and consumed through the Python adapter, with Python fallback retained for resilience.
