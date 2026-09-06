# Contributing

Contributions are welcome, including bug reports, documentation fixes and new
functionality.

## Reporting a problem

Open an issue at
<https://github.com/mickwh2764/par2discovery/issues>. For a bug, please include
the package version (`python -c "import par2; print(par2.__version__)"`), your
Python version, and a minimal series or CSV that reproduces the behaviour. A
statistical objection — for example that an estimate is biased in a regime the
benchmark does not cover — is a valid bug report and does not need code.

## Seeking support

Usage questions belong in an issue too, so that the answer is findable. For
anything relating to the patent boundary described in [NOTICE](NOTICE), email
mickwh@msn.com.

## Contributing code

1. Fork the repository and create a branch.
2. Install in editable mode with the test dependencies:
   `pip install -e . && pip install pytest`.
3. Add or update tests in `tests/`; new numerical behaviour needs a test whose
   expected value is derived analytically or from a known-answer simulation
   rather than from the current output.
4. Run `python -m pytest` — the suite must pass.
5. Open a pull request describing what changed and why.

Contributions are accepted under the Apache-2.0 licence of the project.
