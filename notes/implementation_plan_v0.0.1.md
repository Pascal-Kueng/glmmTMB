# v0.0.1 implementation plan: separable AR1 member-time covariance

## Goal

Implement a fully usable first version of separable member x time covariance structures in `glmmTMB` for:

- `unxar1(membertime(member, time) + 0 | group)`
- `homcsxar1(membertime(member, time) + 0 | group)`

The target covariance for one group is:

```text
Corr((member_i, time_t), (member_j, time_s))
  = Corr_member(i, j) * phi ^ abs(t - s)
```

For Gaussian models, `dispformula = ~0` gives the pure separable residual covariance. A non-zero Gaussian dispersion model adds independent nugget variance. For non-Gaussian models, the separable term is a latent Gaussian structured random effect on the link scale.

## Scope and positioning

`homcsxar1()` is the clean v0.0.1 structure for exchangeable/indistinguishable dyads. It estimates one member SD, one same-time member correlation, and one AR(1) time parameter.

For `homcsxar1()`, the first `membertime()` coordinate is an exchangeable member-position index. Labels may be arbitrary, but they must be stable within group across time so the AR(1) margin follows the same individual. If the mean model uses `Idiff`, use the same stable member assignment for `membertime()`. v0.0.1 should warn about this.

`unxar1()` is a useful but restricted distinguishable-dyad structure. It estimates role/member-specific SDs and an unstructured same-time member correlation matrix, but it still imposes one common AR(1) parameter across roles:

```text
phi_A = phi_B
```

For two distinguishable members, it also implies symmetric lagged cross-role covariance:

```text
Corr(A_t, B_{t+k}) = Corr(B_t, A_{t+k}) = rho_AB * phi^k
```

Therefore v0.0.1 should not claim to fully solve distinguishable dyadic ILD covariance. It should also not claim that `homcsxar1()` is the final exchangeable-dyad dynamic model: it is the clean separable exchangeable AR(1) target, but it has no cross-lag transition parameter.

The priority follow-up is a small VAR(1) family:

- `var1diag()` for distinguishable dyads, with role-specific AR parameters and correlated innovations;
- full `var1()` for distinguishable dyads, adding cross-lagged transition parameters;
- a symmetric/exchangeable VAR(1), tentatively `var1sym()` or `var1exch()`, that is invariant to swapping arbitrary dyad-member labels.

## Non-goals for v0.0.1

- No OU/continuous-time margin.
- No general `sep(member = ..., time = ...)` formula syntax.
- No nonseparable VAR(1) structures, including distinguishable `var1diag()`/`var1()` or symmetric exchangeable VAR(1).
- No mixed distinguishable/exchangeable groups in one covariance term.
- No Kronecker-optimized likelihood; use dense per-group MVN likelihood first.
- No large vignette set. Keep docs minimal and add fuller documentation after validation.

## User-facing syntax

Use `membertime()` coordinates with member first and discrete time second:

```r
fit_un <- glmmTMB(
  y ~ x + unxar1(membertime(member, time) + 0 | group),
  data = dat,
  family = gaussian(),
  dispformula = ~0
)

fit_homcs <- glmmTMB(
  y ~ x + homcsxar1(membertime(member, time) + 0 | group),
  data = dat,
  family = gaussian(),
  dispformula = ~0
)
```

The alias `mt(member, time)` is equivalent to `membertime(member, time)`.

`member` and `time` may be numeric or factor values accepted by `membertime()`, but after parsing the coordinate grid must be complete, sorted, and representable as `M * T` random-effect columns.

For incomplete longitudinal data, define the coordinate factor against the intended full member x time grid, not only against observed rows if some member-time positions are absent globally. A safe pattern is:

```r
members <- sort(unique(dat$member))
times <- sort(unique(dat$time))
grid <- expand.grid(member = members, time = times)

full_member_time <- glmmTMB::numFactor(grid$member, grid$time)
obs_member_time <- glmmTMB::numFactor(dat$member, dat$time)
dat$member_time <- factor(as.character(obs_member_time),
                          levels = levels(full_member_time))
```

Missing observations within groups are allowed: a group does not need observed rows for every member-time cell. The global coordinate levels must still describe a valid full member x time grid so the covariance structure is well defined.

## Parameterization

### `unxar1`

For `M` member levels and `T` time levels:

- block size: `M * T`
- parameters: `M` member log-SDs + `M * (M - 1) / 2` member correlation parameters + `1` AR1 parameter
- theta length:

```text
M + M * (M - 1) / 2 + 1
```

Member covariance:

- log-SDs use `exp(theta)`
- member correlations use TMB's `density::UNSTRUCTURED_CORR_t`

Time covariance:

```text
phi = theta_phi / sqrt(1 + theta_phi^2)
R_time[t, s] = phi ^ abs(t - s)
```

Full covariance:

```text
Sigma = diag(sd_member_by_coordinate) * Corr_sep * diag(sd_member_by_coordinate)
Corr_sep[a, b] = Corr_member[m_a, m_b] * phi ^ abs(time_a - time_b)
```

where `sd_member_by_coordinate[a] = sd_member[m_a]`.

### `homcsxar1`

For `M` member levels and `T` time levels:

- block size: `M * T`
- parameters: common log-SD + common member correlation + AR1 parameter
- theta length: `3`

Member correlation should use the existing homogeneous compound-symmetry transform:

```text
a = 1 / (M - 1)
rho = plogis(theta_member) * (1 + a) - a
```

This gives the valid range `[-1 / (M - 1), 1]`.

Time covariance uses the same AR1 transform as above.

## Required file changes

### `glmmTMB/src/glmmTMB.cpp`

Add covariance enum entries after the existing values to avoid changing old serialized model meanings:

```cpp
unxar1_covstruct = 16,
homcsxar1_covstruct = 17
```

Extend `per_term_info` with enough metadata for separable structures:

- `int sepNumMembers`
- `int sepNumTimes`
- `vector<Type> sepMembers`
- `vector<Type> sepTimes`

Add parsing of these optional fields in `terms_t`.

Add a dense MVN branch in `termwise_nll()` for both structures:

1. Validate `blockSize == sepNumMembers * sepNumTimes`.
2. Build member correlation.
3. Build full `blockSize x blockSize` separable correlation from parsed member/time coordinates.
4. Build coordinate-level SD vector.
5. Evaluate one MVN density per group.
6. Implement simulation for `random_simcode`; support `fix_simcode` and `zero_simcode` by matching behavior of existing structures where practical.
7. Populate `term.corr` and `term.sd` for `VarCorr()`.

Implementation can initially use `density::MVNORM_t<Type>` plus `density::VECSCALE_t`.

### `glmmTMB/R/glmmTMB.R`

Update `getReStruc()`:

- Add `parFun()` entries:

```r
"unxar1" = M + M * (M - 1) / 2 + 1
"homcsxar1" = 3
```

where `M` must be derived from parsed coordinate levels, not from `blockSize` alone.

- For `unxar1` and `homcsxar1`, parse `reTrms$cnms[[i]]` with `parseNumLevels()`.
- Require exactly two coordinate columns.
- Treat column 1 as member and column 2 as time.
- Require no intercept.
- Require a complete Cartesian grid in the coordinate factor levels. This is a global grid-level validation, not a requirement that every group has every observation.
- Require one random-effect column per coordinate.
- Require time values to be sorted with unit-spaced integer positions for v0.0.1.
- Store in the term list:

```r
tmp$sepNumMembers <- length(unique(members))
tmp$sepNumTimes <- length(unique(times))
tmp$sepMembers <- match(members, sort(unique(members)))
tmp$sepTimes <- match(times, sort(unique(times)))
```

Add a dedicated validation helper, e.g. `checkSepAr1Coord()`, so all errors are specific and actionable. Error messages should include both what failed and the expected construction pattern.

Required validation errors:

- If levels cannot be parsed:

```text
unxar1()/homcsxar1() require two-dimensional member-time coordinates. Use e.g. unxar1(membertime(member, time) + 0 | group), or precompute dat$member_time <- glmmTMB::numFactor(dat$member, dat$time) and specify unxar1(member_time + 0 | group).
```

- If coordinates are not two-dimensional:

```text
unxar1()/homcsxar1() expect membertime(member, time) or numFactor(member, time) with exactly two coordinates: member first, discrete AR1 time second. The parsed coordinate factor has <k> coordinate column(s).
```

- If an intercept is included:

```text
unxar1()/homcsxar1() must be specified without an intercept, e.g. unxar1(membertime(member, time) + 0 | group). The random-effect columns must be one column per member-time coordinate.
```

- If the global coordinate levels are not a complete member x time grid:

```text
unxar1()/homcsxar1() require coordinate factor levels to form a complete member x time grid. Missing observed rows within groups are allowed, but the factor levels must include all intended member-time combinations. Use membertime(member, time) for complete observed grids. For globally missing cells, build levels from a full grid, e.g. full <- glmmTMB::numFactor(grid$member, grid$time); dat$member_time <- factor(as.character(glmmTMB::numFactor(dat$member, dat$time)), levels = levels(full)).
```

- If time values are not valid for AR1 v0.0.1:

```text
unxar1()/homcsxar1() v0.0.1 require the second numFactor coordinate to be discrete unit-spaced time positions, e.g. 1, 2, ..., T. For irregular continuous time, use is not supported until the OU implementation.
```

- If member/time order appears reversed, where this can be detected heuristically:

```text
unxar1()/homcsxar1() interpret numFactor coordinate 1 as member and coordinate 2 as time. The parsed grid looks like time may have been supplied first. Use glmmTMB::numFactor(member, time).
```

Do not rely on the heuristic for correctness; the definitive checks are dimensionality, Cartesian completeness, and unit-spaced time values.

Also add warning orchestration after conditional `condReStruc` is constructed:

- detect whether any conditional term has block code `unxar1` or `homcsxar1`
- warn based on `family$family` and original `dispformula`

Warning behavior:

- Gaussian + `dispformula = ~0`: no dispersion warning. `homcsxar1()` may still warn about stable exchangeable member-position assignment.
- Gaussian + `dispformula != ~0`: warn that the model estimates separable latent/residual covariance plus independent nugget variance; suggest `dispformula = ~0` for pure separable residual covariance.
- Gaussian + structured dispersion formula such as `~ 0 + role`: use a more specific warning that this is role-/covariate-specific nugget variance.
- Families without dispersion (`binomial`, `poisson`, `truncated_poisson`, `bell`): warn only if the user supplies a non-default/non-trivial `dispformula`; note that `dispformula` is ignored and the separable term remains a latent Gaussian random effect.
- Non-Gaussian families with dispersion: warn if `dispformula` is non-default or `~0`; do not recommend `~0` as a standard fix. Explain that `dispformula` controls family-specific dispersion, not Gaussian residual variance.

Important implementation detail: `mkTMBStruc()` currently converts `dispformula = ~0` for Gaussian into an internal mapped `~1`, and resets no-dispersion families to `~0`. Warning decisions must use `dispformula.orig` before these internal rewrites.

### `glmmTMB/R/enum.R`

Regenerate via:

```sh
make enum-update
```

Do not edit by hand.

### Tests

Add focused tests, preferably in `glmmTMB/tests/testthat/test-varstruc.R` or a new `test-xar1-covstruct.R`.

Required v0.0.1 tests:

1. Formula parsing accepts `unxar1(member_time + 0 | group)`.
2. Formula parsing accepts `homcsxar1(member_time + 0 | group)`.
3. Invalid coordinate input errors:
   - non-`numFactor`/unparseable levels
   - one-dimensional coordinate
   - three-dimensional coordinate
   - missing member-time grid cells
   - non-unit/non-sorted time positions
   - intercept included
- actionable error messages mention `membertime(member, time)`, `numFactor(member, time)`, `+ 0`, and full-grid factor levels
4. Matrix construction through `VarCorr()`:
   - for `unxar1`, same-member lag `k` correlation equals `phi^k`
   - cross-member lag `k` correlation equals `rho_member * phi^k`
   - role-specific SDs repeat across times
   - document/test that `unxar1` uses one common `phi` across distinguishable roles
   - for `homcsxar1`, member SDs are equal and cross-member correlations are exchangeable
5. Simulation/fitting smoke tests:
   - small balanced Gaussian dyad data
   - small unbalanced Gaussian dyad data with member-time observations missing within some groups
   - `dispformula = ~0`
   - convergence object is a `glmmTMB`
   - `VarCorr()` has expected dimensions and finite correlations
6. Warning tests:
   - Gaussian + `dispformula = ~0`: no separable-dispersion warning
   - Gaussian + default `~1`: warning
   - Gaussian + structured dispersion formula: warning
   - Poisson/binomial with explicit non-default `dispformula`: ignored-dispersion warning
   - Negative-binomial with `dispformula = ~0` or structured dispersion: family-specific warning

Use deterministic parameter checks where possible via fixed `start`/`map` or `doFit = FALSE`; avoid fragile full recovery tests in v0.0.1.

### Minimal documentation

Update the covariance-structure list in the `glmmTMB()` roxygen block.

Add a compact example to existing covariance documentation or a short developer note:

- `membertime(member, time)` or `mt(member, time)` should be used for complete observed grids.
- precomputed `numFactor(member, time)` factors are the explicit fallback for globally incomplete grids.
- first coordinate is member, second is discrete AR1 time.
- use `+ 0`.
- use `dispformula = ~0` for pure Gaussian separable residual covariance.
- non-Gaussian interpretation is latent Gaussian on the link scale.
- `homcsxar1()` is the clean exchangeable/indistinguishable dyad structure.
- `unxar1()` is a restricted distinguishable structure with a common AR parameter, not a full distinguishable VAR model.
- role-specific AR parameters and asymmetric lagged cross-role covariance are planned for `var1diag()`, not v0.0.1.

Full dyadic APIM/DIM vignettes are out of scope for v0.0.1.

### Publication-readiness additions

For v0.0.1 to be usable in a real analysis manuscript, include a small but complete validation and reporting layer, not just the covariance implementation.

Required additions:

- Parameter extraction guidance:
  - document how to recover `phi`, member SDs, and member correlations from `VarCorr()` and/or `getME(fit, "theta")`
  - include helper examples for dyadic `unxar1` and `homcsxar1`
  - make parameter ordering explicit
  - label `unxar1` as a common-AR distinguishable model in all examples
- Confidence intervals:
  - verify that `confint()` works for the new theta parameters at least with Wald intervals
  - document profile-CI limitations if profiling is slow or unstable
  - ensure parameter names in output are interpretable enough for users to report
- Simulation validation:
  - balanced Gaussian dyads with known `rho_member` and `phi`
  - unbalanced/missing-observation Gaussian dyads
  - exchangeable dyads with arbitrary member-label swaps for `homcsxar1`
  - restricted distinguishable dyads showing that both roles share the same recovered AR parameter under `unxar1`
  - enough replications to show bias and convergence behavior for manuscript-supporting examples
- Cross-check against reference models:
  - compare `unxar1` to `us(member_time + 0 | group)` for small `T`, where the estimated unstructured covariance should show the same separable pattern
  - compare fixed-parameter likelihood/covariance construction against an explicit `kronecker()` matrix in R
  - optionally compare one or two small Gaussian cases against a standalone `mvtnorm` likelihood or Stan reference model
- Missing-data behavior:
  - demonstrate that missing outcome rows within groups are handled correctly when the global coordinate factor levels are complete
  - state that missing predictor/group/member/time values are handled by the existing model-frame/`na.action` workflow and may drop rows before model construction
  - include a test where one group is missing one partner-day and another group is complete
- Diagnostics and failure modes:
  - document boundary estimates for `phi` and member correlations
  - document what non-positive-definite Hessian warnings mean in this context
  - recommend checking convergence, `fit$sdr$pdHess`, and sensitivity to starting values for publication analyses
  - state that evidence of role-specific inertia should motivate `var1diag()`/VAR-style modeling rather than `unxar1()`
- Reproducible example:
  - include a compact, runnable simulated dyad example that estimates and reports the target quantities
  - keep it small enough for tests/docs, but realistic enough to show the publication workflow

Nice-to-have additions:

- A convenience extractor, e.g. an internal or exported helper that returns a named list with `sd_member`, `cor_member`, `phi`, and the full separable correlation matrix. This is not strictly required if `VarCorr()` is clear, but it would reduce user error.
- A short comparison note explaining why `homcs(member + 0 | group:time) + ar1(time + 0 | person)` is not equivalent to `homcsxar1(member_time + 0 | group)`.
- A warning when `M * T` is large enough that dense MVN evaluation will be slow or memory-heavy.

### Priority VAR(1) extensions after v0.0.1

The next implementation priority is a nonseparable VAR(1) family for dyadic ILD. There should be both distinguishable and exchangeable/symmetric targets.

#### Distinguishable diagonal VAR(1)

The first distinguishable-dyad target is:

```r
var1diag(membertime(role, time) + 0 | group)
```

Target process:

```text
A_t = phi_A * A_{t-1} + eps_A,t
B_t = phi_B * B_{t-1} + eps_B,t
Cov(eps_A,t, eps_B,t) != 0
```

This should support:

- role-specific AR parameters;
- role-specific stationary variances;
- correlated same-time innovations;
- asymmetric lagged cross-role covariance when `phi_A != phi_B`;
- the same Gaussian/non-Gaussian latent-process interpretation and dispersion warnings used for v0.0.1.

This is the main path beyond `unxar1()` when roles are meaningful but cross-lagged transition parameters are not yet needed.

#### Distinguishable full VAR(1)

The later full distinguishable target is:

```r
var1(membertime(role, time) + 0 | group)
```

This should allow a full 2 x 2 transition matrix:

```text
A_t = phi_AA * A_{t-1} + phi_AB * B_{t-1} + eps_A,t
B_t = phi_BA * A_{t-1} + phi_BB * B_{t-1} + eps_B,t
```

The transition matrix must be parameterized so its eigenvalues remain inside the unit circle.

#### Exchangeable symmetric VAR(1)

Exchangeable/indistinguishable dyads also need a nonseparable VAR(1) target, because `homcsxar1()` is separable and does not model symmetric cross-lagged influence. Tentative names:

```r
var1sym(membertime(member, time) + 0 | group)
var1exch(membertime(member, time) + 0 | group)
```

Target process in member coordinates:

```text
A_t = alpha * A_{t-1} + gamma * B_{t-1} + eps_A,t
B_t = gamma * A_{t-1} + alpha * B_{t-1} + eps_B,t
```

with exchangeable innovations:

```text
Var(eps_A,t) = Var(eps_B,t)
Cov(eps_A,t, eps_B,t) unrestricted within the exchangeable bounds
```

This is equivalent to fitting independent AR(1) dynamics in sum/difference coordinates:

```text
S_t = phi_S * S_{t-1} + eps_S,t
D_t = phi_D * D_{t-1} + eps_D,t
```

where `phi_S = alpha + gamma` and `phi_D = alpha - gamma`. The stationarity constraints are therefore `abs(phi_S) < 1` and `abs(phi_D) < 1`. This parameterization is attractive because it preserves exchangeability and gives a simple stable transformation.

Initial implementation constraints should probably be:

- exactly two members/roles;
- unit-spaced discrete time;
- stationary covariance construction;
- dense MVN likelihood first;
- `var1diag()` has no cross-lagged transition parameters;
- symmetric exchangeable VAR(1) has one label-swap-invariant cross-lag parameter;
- unrestricted cross-lagged transition parameters are deferred to full distinguishable `var1()`.

A later full distinguishable bivariate VAR(1) should allow unrestricted cross-lagged transition parameters and will require a stable parameterization of the 2 x 2 transition matrix.

## Development sequence

1. Add C++ enum values and regenerate `R/enum.R`.
2. Add R-side coordinate parsing and parameter counts.
3. Add minimal C++ dense covariance branch for `homcsxar1` as the exchangeable dyad target.
4. Compile and fit a tiny balanced Gaussian exchangeable dyad model.
5. Add restricted/common-AR `unxar1` by reusing TMB unstructured member correlation.
6. Add `VarCorr()` reporting checks.
7. Add dispersion-warning helper and warning tests.
8. Add invalid-input tests.
9. Run focused tests:

```sh
Rscript -e 'devtools::test("glmmTMB", filter = "varstruc|sep")'
```

10. Run broader covariance tests:

```sh
Rscript -e 'devtools::test("glmmTMB", filter = "varstruc|VarCorr|simulate_new")'
```

11. If stable, run full package tests:

```sh
make test
```

## Reliability criteria for v0.0.1

v0.0.1 is considered workable when:

- both new structures compile and fit balanced Gaussian dyad examples
- both new structures fit unbalanced dyad examples where some member-time observations are missing within groups
- expected covariance identities are verified from `VarCorr()`
- invalid coordinate specifications fail early on the R side
- validation errors give enough guidance to fix malformed `numFactor()` inputs
- warning behavior matches the planned Gaussian/non-Gaussian interpretation
- publication-facing examples show parameter recovery and missing-observation handling
- examples and documentation clearly state that `unxar1` is a restricted/common-AR distinguishable model
- confidence intervals and parameter reporting are either verified or their limitations are explicitly documented
- existing covariance tests still pass

## Known risks

- Dense MVN scales poorly for large `M * T`; acceptable for dyads and moderate diary length in v0.0.1.
- `VarCorr()` will report the full member-time covariance matrix, which can be large but is useful for validation.
- Missing observations within groups are supported. What is unsupported is an ill-defined coordinate factor whose levels do not describe a complete global member x time grid.
- `unxar1` should not be oversold as a full distinguishable dyadic ILD covariance model. It does not allow role-specific AR parameters or cross-lagged transition effects.
- The name `unxar1` uses `un` for the member margin. If upstream maintainers prefer glmmTMB's existing `us` terminology, rename to `usxar1` before PR discussion.

## Pre-PR cleanup

This plan is a fork-local implementation aid. Before opening a `glmmTMB` pull request:

- remove `notes/implementation_plan.md` and `notes/implementation_plan_v0.0.1.md` from the final patch unless maintainers explicitly ask for them;
- keep only package-relevant code, tests, and generated documentation;
- convert any useful planning text into the GitHub issue, PR description, help text, or a focused vignette;
- decide separately whether exploratory scripts in `notes/` should remain local, be removed, or be converted into formal tests.
