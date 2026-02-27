# Registering Jexpresso in the Julia General Registry

This document explains how to register Jexpresso in the Julia General Registry, which will allow users to install it using `add Jexpresso` instead of cloning the repository manually.

## Prerequisites

Before registering, ensure that:

1. ✅ The `Project.toml` has a valid UUID
2. ✅ The `Project.toml` includes Julia version compatibility (e.g., `julia = "1.11"`)
3. ✅ All dependencies have appropriate compat bounds in the `[compat]` section
4. ✅ The package is properly structured with `src/Jexpresso.jl` as the main module
5. ✅ All changes are committed and pushed to the main branch
6. ✅ The repository is hosted on GitHub

All of these prerequisites are now satisfied.

## Registration Process

### Step 1: Create a GitHub Release

Before registering with JuliaRegistrator, you should create a release on GitHub:

1. Go to your repository on GitHub: https://github.com/smarras79/Jexpresso
2. Click on "Releases" → "Create a new release"
3. Create a new tag (e.g., `v0.1.0`) - the version should match what's in `Project.toml`
4. Add a release title and description
5. Click "Publish release"

### Step 2: Register with JuliaRegistrator

There are two ways to register the package:

#### Option A: Using the JuliaRegistrator GitHub App (Recommended)

1. Install the [JuliaRegistrator GitHub App](https://github.com/apps/julia-registrator) on your repository
2. Comment on a commit or pull request (on the master/main branch) with:
   ```
   @JuliaRegistrator register
   ```
3. The bot will automatically create a pull request to the Julia General registry
4. Wait for automated checks to pass and for a registry maintainer to merge

#### Option B: Using the Registrator Web Interface

1. Go to https://juliahub.com/ui/Registrator
2. Sign in with your GitHub account
3. Enter your repository URL: `https://github.com/smarras79/Jexpresso`
4. Select the branch to register from (usually `master` or `main`)
5. Click "Submit"

### Step 3: Wait for Registry PR Review

After triggering registration, JuliaRegistrator will:

1. Create a pull request in the [General registry](https://github.com/JuliaRegistries/General)
2. Run automated checks (AutoMerge) to verify the package meets requirements
3. If checks pass, a registry maintainer will review and merge the PR
4. Once merged, the package becomes available via `add Jexpresso`

This process typically takes a few hours to a few days, depending on whether any issues are found.

## Common Issues and Solutions

### Issue 1: Missing Compat Entries

**Error**: "Package has dependencies without compat bounds"

**Solution**: Ensure all dependencies in `[deps]` have corresponding entries in `[compat]`. This has been addressed in the latest `Project.toml`.

### Issue 2: Invalid Version Bounds

**Error**: "Compat bounds are invalid"

**Solution**: Use proper semantic versioning in compat bounds:
- `"1"` means `≥1.0.0, <2.0.0`
- `"1.2"` means `≥1.2.0, <2.0.0`
- `"^1.2.3"` means `≥1.2.3, <2.0.0`
- `"=1.2.3"` means exactly version 1.2.3 (use sparingly)

### Issue 3: Package Name Already Registered

**Error**: "Package name already exists in registry"

**Solution**: Choose a different package name or contact the existing package maintainer.

### Issue 4: Repository Not on GitHub

**Error**: "Repository must be hosted on GitHub"

**Solution**: The Julia General registry currently only supports GitHub-hosted repositories.

## After Registration

Once Jexpresso is registered:

1. **Update documentation**: Update `README.md` and `INSTALLATION.md` to reflect that users can now use:
   ```julia
   using Pkg
   Pkg.add("Jexpresso")
   ```

2. **Announce the release**: Notify users through:
   - GitHub release notes
   - Julia Discourse forum
   - Social media or project website

3. **Future updates**: For subsequent versions:
   - Update the version in `Project.toml`
   - Create a new GitHub release
   - Comment `@JuliaRegistrator register` on a commit/PR
   - JuliaRegistrator will automatically create a new registry PR

## Version Numbering

Follow [Semantic Versioning](https://semver.org/):

- **Major version** (x.0.0): Breaking changes that affect backward compatibility
- **Minor version** (0.x.0): New features that are backward compatible
- **Patch version** (0.0.x): Bug fixes that are backward compatible

Current version in `Project.toml`: `0.1.0`

## Useful Resources

- [JuliaRegistrator Documentation](https://github.com/JuliaRegistries/Registrator.jl)
- [Julia Package Guidelines](https://julialang.github.io/Pkg.jl/v1/creating-packages/)
- [General Registry](https://github.com/JuliaRegistries/General)
- [Julia Discourse - Package Announcements](https://discourse.julialang.org/c/package-announcements/10)

## Contacts

For questions about registration:
- Simone Marras: smarras@njit.edu
- Yassine Tissaoui: tissaoui@wisc.edu
- Hang Wang: hang.wang@njit.edu

For questions about the Julia registry system:
- [Julia Discourse](https://discourse.julialang.org/)
- [Julia Slack](https://julialang.org/slack/)
