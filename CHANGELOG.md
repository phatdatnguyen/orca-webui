# Changelog

## 0.1.11 — Unreleased

### Added

- Show **Running calculation...** in the shared status area while an ORCA calculation
  runs, then replace it with the completion, error, or stop message.

### Fixed

- Run ORCA from the selected working directory using an absolute executable path.
  Validate input and output paths before rotating logs, restore the previous log if
  launch fails, and reject concurrent runs within one app instance.
- Preserve calculation logs on startup and restrict transient cleanup to `static/`.
  Use unique structure-viewer filenames and honor an explicitly configured run directory.
- Reject filenames and symlinks that point outside the selected working directory.
  Prevent uploads and conformer generation from overwriting existing files.
- Save editor contents atomically to the file actually opened in the editor. Clear
  stale selections, editor contents, and loaded results when switching directories.
- Start successfully without the optional `gradio_molecule2d` editor, using SMILES
  entry. Validate conformer settings and skip failed or unconverged optimizations.
- Parse standard and legacy XYZ headers without dropping atoms or mistaking the
  comment for coordinates; reject malformed or truncated structures.
- Convert molecular-orbital energies to hartree and display both spin channels for
  unrestricted calculations. Handle missing optional energy, frequency,
  thermochemistry, and transition metadata; validate transition arrays before plotting.
- Align transition orbital numbering with the MO table, reject empty data exports,
  and clear cached results after loading errors.
- Render equivalent-partner NMR multiplets using Pascal weights directly, fall back
  to defaults for nonfinite reference shieldings, and keep the chemical-shift axis descending.
- Escape error text inserted into HTML status messages and add regression coverage
  for execution, workspace handling, structure parsing, and results.

### Release preparation

- Align package metadata and the runtime version at 0.1.11.
- Update user and contributor documentation, including version-specific build,
  validation, and publication commands. Publication is a separate maintainer step.
