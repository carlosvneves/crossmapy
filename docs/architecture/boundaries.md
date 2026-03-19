# Architectural Boundaries

- `python/crossmapy` is the only public API surface.
- `rust/ccm-rs` is the computational core and must remain language-agnostic.
- `crates/ccm-rs-py` is a thin binding layer and must avoid domain logic.
- Algorithm changes happen in `rust/ccm-rs/src/*` and are exposed through adapters.
