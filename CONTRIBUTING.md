# Contributing Guidelines

## Welcome! 👋

Thank you for your interest in contributing to this PCB snRNA-seq analysis!

## Getting Started

### For Small Changes (Documentation, typos)
1. Fork the repository
2. Make your edits
3. Submit a pull request with a clear description

### For Code Changes
1. Fork the repository
2. Create a feature branch: `git checkout -b feature/your-feature-name`
3. Make your changes following code standards
4. Test your changes
5. Submit a pull request

## Code Standards

Please follow professional R standards:

- Use `snake_case` for variables and functions
- Add comprehensive headers to scripts with purpose, author, date
- Include error handling with try-catch
- Use logging statements (from scripts/utils/logging.R)
- Use relative paths (not absolute paths)
- Set random seeds for reproducibility

## Workflow for Code Changes

1. Create feature branch
```bash
   git checkout -b feature/your-feature
```

2. Make changes and test
```bash
   Rscript scripts/path/to/your_script.R
```

3. Follow code standards
   - Comprehensive header
   - Error handling
   - Logging statements
   - Reproducibility parameters

4. Commit with clear message
```bash
   git commit -m "feat: brief description

   - Detail what changed
   - Why it changed
   - Any implications"
```

5. Push to your fork
```bash
   git push origin feature/your-feature
```

6. Submit pull request on GitHub

## Pull Request Guidelines

- **Title**: Start with type (feat:, fix:, docs:, refactor:)
- **Description**: Explain what changed and why
- **Testing**: Note what testing was done
- **Issues**: Reference related issues (#123)

## Questions?

- Check `README.md` for usage
- See `docs/` for documentation
- Search existing issues/PRs

---

Thank you for contributing! 🙌
