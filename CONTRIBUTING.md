# Contributing

Open an issue with a reproducible numerical example and a description of the model assumptions. Contributions should include unit tests for any new physical model or numerical algorithm. Submit pull requests against the default branch; the GitHub Actions test matrix must pass. Describe normalization, boundary conditions, stability restrictions, and references in documentation. Do not claim general or experimental validation without evidence.

## Developer install

```bash
python -m pip install -e '.[test]'
python -m pytest -q
```
