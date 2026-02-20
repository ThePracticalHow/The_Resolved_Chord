# Contributing to LOTUS

Thank you for your interest in contributing to LOTUS (Lens Orbifold Theory of the Unified Spectrum)! This document provides guidelines for contributors.

## Code of Conduct

This project follows a code of conduct to ensure a welcoming environment for all contributors. By participating, you agree to:
- Be respectful and inclusive
- Focus on constructive feedback
- Accept responsibility for mistakes
- Show empathy towards other contributors

## Ways to Contribute

### 🐛 Bug Reports
- Use the [issue tracker](https://github.com/ThePracticalHow/The_Resolved_Chord/issues) on GitHub
- Include detailed steps to reproduce the bug
- Specify your Python version, operating system, and numpy version
- Include error messages and stack traces

### 💡 Feature Requests
- Open an [issue](https://github.com/ThePracticalHow/The_Resolved_Chord/issues) with the "enhancement" label
- Clearly describe the proposed feature and its benefits
- Consider how it fits with the project's goals

### 🛠️ Code Contributions
- Fork the repository
- Create a feature branch from `main`
- Make your changes
- Add tests for new functionality
- Ensure all tests pass
- Submit a pull request

### 📚 Documentation
- Improve existing documentation
- Add examples or tutorials
- Translate documentation to other languages
- Fix typos or clarify confusing sections

### 🧪 Testing
- Add test cases for existing functionality
- Improve test coverage
- Test on different platforms/Python versions

## Development Setup

1. **Clone the repository:**
   ```bash
   git clone https://github.com/ThePracticalHow/The_Resolved_Chord.git
   cd The_Resolved_Chord/05_Project_LENG/public-release
   ```

2. **Set up a virtual environment:**
   ```bash
   python -m venv venv
   source venv/bin/activate  # On Windows: venv\Scripts\activate
   ```

3. **Install dependencies:**
   ```bash
   pip install -e ".[dev]"
   ```

4. **Run tests:**
   ```bash
   python -m pytest tools/falsification/ -v
   ```

## Code Style

- Follow [PEP 8](https://pep8.org/) style guidelines
- Use meaningful variable and function names
- Add docstrings to all public functions and classes
- Keep functions focused on single responsibilities
- Use type hints where appropriate

## Testing Guidelines

- Write tests for all new functionality
- Ensure tests are independent and reproducible
- Use descriptive test names that explain what they're testing
- Test both success and failure cases
- Aim for high test coverage

## Commit Guidelines

- Use clear, descriptive commit messages
- Start with a verb in imperative mood (e.g., "Add", "Fix", "Update")
- Reference issue numbers when applicable
- Keep commits focused on single changes

Example:
```
Fix neutrino mass calculation in cosmology.py

- Corrected the formula for m_ν₃ calculation
- Added test case for the fix
- Updated documentation

Closes #123
```

## Pull Request Process

1. **Create a Pull Request** from your feature branch to `main`
2. **Fill out the PR template** with:
   - Clear description of changes
   - Link to related issues
   - Screenshots/demos if applicable
3. **Ensure CI passes** - all tests must pass
4. **Request review** from maintainers
5. **Address feedback** and make requested changes
6. **Merge** once approved

## Scientific Contributions

Since LOTUS is a scientific project, contributions involving physics or mathematics should:

- Include references to relevant literature
- Provide mathematical justification for changes
- Consider implications for existing predictions
- Maintain consistency with the overall framework

### Prediction Count Messaging Policy

To keep release messaging consistent across docs and scripts:

- Use **87 predictions** for current framework descriptions and user-facing summaries.
- Keep older counts (e.g., **48**) only for historical releases or archived artifacts.
- Label older counts explicitly as **historical** whenever they appear.

### AI Collaboration Guardrails

For AI-assisted edits (Copilot, Claude, etc.), keep handoffs deterministic:

- Multi-agent protocol reference (Claude, Composer, Copilot, Gemini, Grok, Perplexity, ChatGPT): `../MULTI_AGENT_HYGIENE_PROTOCOL.md`

- Update canonical docs first: `README.md` (public claims), `MASTER_CODE_INDEX.md` (ops map), `V11_INVENTORY.md` (numbered registry).
- Preserve append-only timelines/history docs; add new dated entries instead of rewriting prior logs.
- Run robustness checks before PR or release:
   ```bash
   python tools/release_gate.py
   python -m pytest tools/falsification/ -m "not slow" -q
   ```
- Prefer `tools/auto_clean.py` for cleanup over ad-hoc shell deletes to keep behavior cross-platform.

## Recognition

All contributors will be:
- Listed in the contributors file
- Acknowledged in release notes
- Recognized for their scientific contributions

## Questions?

If you have questions about contributing, please:
- Check existing [issues](https://github.com/ThePracticalHow/The_Resolved_Chord/issues) and [discussions](https://github.com/ThePracticalHow/The_Resolved_Chord/discussions)
- Open a new discussion for general questions
- Contact the maintainers directly for sensitive matters

Thank you for contributing to LOTUS! 🚀