# AdVentSeq — project rules

## Versioning: bump the version before every commit

This project uses semantic versioning (`MAJOR.MINOR.PATCH`). Before creating any
commit that changes the package code, bump the version **in the same commit**:

- **PATCH** (`x.y.Z+1`) — bug fixes and other backward-compatible corrections.
- **MINOR** (`x.Y+1.0`) — new, backward-compatible features.
- **MAJOR** (`X+1.0.0`) — backward-incompatible / breaking changes.

Update the version in **both** places, keeping them identical:

- `pyproject.toml` → `version = "..."`
- `AdVentSeq/__init__.py` → `__version__ = "..."`

Changes that do not touch the package (e.g. docs-only edits, `examples/`, CI, or
this file) do not require a version bump.

## Git: never create branches automatically

Do **not** create a new git branch on your own — commit on the current branch
(including `main`) instead. Only create a branch when I explicitly ask for one.
