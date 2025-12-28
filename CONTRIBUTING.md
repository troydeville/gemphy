# Contributing to gemphy

Thank you for your interest in contributing to the Geometric Encoded Medium (GEM) Physics library! This project aims to derive fundamental constants from geometric impedance in a horn torus vacuum. Contributions are welcome from physicists, Rust developers, and enthusiasts—whether it's code, docs, tests, or theory refinements.

## Code of Conduct

Please follow the [Contributor Covenant Code of Conduct](https://www.contributor-covenant.org/version/2/1/code_of_conduct.html) to ensure a positive environment for everyone.

## How to Contribute

### Reporting Issues

- Use GitHub Issues for bugs, feature requests, or questions.
- Search existing issues first to avoid duplicates.
- Provide details: steps to reproduce, expected vs. actual behavior, Rust version, and relevant logs.

### Submitting Pull Requests (PRs)

1. **Fork the Repo**: Click "Fork" on GitHub and clone your fork: `git clone https://github.com/YOUR_USERNAME/gemphy.git`
2. **Create a Branch**: `git checkout -b feature/your-feature` or `fix/your-bugfix`
3. **Make Changes**:
   - Follow Rust style: Run `cargo fmt` and `cargo clippy` before committing.
   - Add tests for new features (use `#[test]` and approx crate for floats).
   - Update docs: Add Rustdoc comments (///) to new publics; build with `cargo doc`.
   - For physics additions (e.g., new derivations), include validations against CODATA.
4. **Commit**: Use clear messages, e.g., "feat: add black hole phase simulation" (follow Conventional Commits).
5. **Push and PR**: `git push origin your-branch`, then create a PR on the main repo. Describe changes, link issues, and note breaking changes.
6. **Review**: Address feedback; CI will run tests automatically.

### Development Setup

- Rust: Stable toolchain (install via rustup.rs).
- Dependencies: Listed in Cargo.toml; run `cargo build` to install.
- Tests: `cargo test -- --nocapture` for verbose output.
- Docs: `cargo doc --open` to view locally.
- Binaries: Run examples like `cargo run --bin electron_proton_action`.

### Areas for Contribution

- **Core Physics**: Improve derivations (e.g., better proton mass precision) or add predictions (dark matter as impedance artifacts).
- **Simulations**: Expand GemSystem for multi-particle dynamics or visualizations.
- **Bindings**: Add Python bindings via pyo3 for broader use.
- **Docs**: Expand mdBook guide with tutorials on horn torus math or complex phases.
- **Tests/Benchmarks**: Add more CODATA validations or perf tests with criterion.rs.

All contributions are under GPL-3.0. Questions? Open an issue or reach out on X (@troy_deville).

Thanks for helping build GEM!
