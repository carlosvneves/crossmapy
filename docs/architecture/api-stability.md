# API Stability Policy

- Public Python imports under `crossmapy` are compatibility targets.
- Rust core APIs are semver-managed and can evolve with explicit migration notes.
- Adapter-level error mapping must preserve deterministic Python exceptions.
