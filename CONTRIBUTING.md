# Contributing

1. Branch from `main`: `git checkout -b <name>/<short-topic>`.
2. Keep changes focused. Add or update tests in `tests/`.
3. Before pushing, run:
   ```bash
   ruff check . && ruff format --check .
   python -m pytest
   ```
4. Open a pull request into `main` and request a review from at least one other team member.

Conventions:

- Follow the paper's notation: δ for distance-to-distress, (λ<sub>c</sub>, γ) for the common factor,
  β<sub>j</sub>(q) for polynomial roots. Firm loadings are **ω** (`omega` in code). Do not call them
  β, because β already means the roots.
- Notebooks go in `notebooks/` with outputs cleared before committing.
- Never commit WRDS / S&P Global data.
