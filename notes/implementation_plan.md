# Project plan: member × time covariance structures for dyadic and small-group intensive longitudinal models in `glmmTMB`

## 1. Working motivation

Standard dyadic intensive longitudinal models need to account for several forms of dependence at once:

1. stable dyad-level non-independence,
2. same-occasion partner/member dependence,
3. temporal autocorrelation within persons,
4. cross-partner/member dependence over time.

Current high-level R workflows can approximate these components, for example by combining same-day member covariance and person-level AR(1), but this is additive and does not reproduce a separable partner/member × time covariance structure such as

\[
\Sigma_{\text{within group}} = \Sigma_{\text{member}} \otimes \Sigma_{\text{time}}.
\]

For dyads, the key residual/latent covariance target is:

\[
\operatorname{Cov}(e_{m,t}, e_{m',s}) = \Sigma_{\text{member}}[m,m'] \times \Sigma_{\text{time}}[t,s].
\]

For an AR(1) time margin:

\[
\Sigma_{\text{time}}[t,s] = \phi^{|t-s|}.
\]

For an OU / continuous-time exponential margin:

\[
\Sigma_{\text{time}}[t,s] = \exp(-\lambda |t-s|).
\]

This directly gives cross-partner lagged dependence, for example:

\[
\operatorname{Corr}(A_t, B_{t+k}) = \rho_{AB}\phi^k.
\]

The immediate goal is to implement this as a new covariance structure in a fork of `glmmTMB`, initially for dyads and then with a path toward triads and arbitrary small groups.

Important scope clarification: separable member × time covariance is a clean target for exchangeable/indistinguishable dyads and a useful restricted model for distinguishable dyads, but it is not the full distinguishable dyadic ILD covariance model. For distinguishable dyads, a fully clean model should allow member-specific autoregression and eventually cross-lagged dependence. That requires nonseparable bivariate VAR-style covariance structures, discussed below as a priority extension.

## 2. Why start with `glmmTMB`

`glmmTMB` is the preferred first target because:

- it already has a covariance-structure architecture for latent random effects;
- it already supports covariance structures such as `ar1`, `ou`, `cs`, `homcs`, `us`, `toep`, `propto`, and `equalto`;
- covariance structures in `glmmTMB` work as latent Gaussian random effects on the linear predictor scale, so they can be used with non-Gaussian families;
- it is much faster than full Bayesian Stan models for iteration, simulation, and repeated testing;
- the project has documented developer guidance for adding covariance structures.

We may still build a small standalone Stan implementation later as a validation/reference model, but the first implementation target should be a `glmmTMB` fork.

## 3. Repository and implementation target

Clone the `glmmTMB` repository:

```bash
git clone https://github.com/glmmTMB/glmmTMB.git
```

Most changes should be inside the `glmmTMB` package code, not in the separate `TMB` package. `glmmTMB` uses TMB as the computational backend. A new covariance structure will likely require:

- adding a new covariance structure enum on the C++ side;
- modifying the C++ negative log-likelihood contribution for structured random effects;
- adding R-side parsing/plumbing;
- adding documentation;
- adding tests and simulation examples.

Relevant implementation areas to inspect:

- `glmmTMB/src/glmmTMB.cpp`
- covariance-structure enum / registry
- `termwise_nll`
- R parsing of covariance structures
- documentation for covariance structures
- tests involving `ar1`, `ou`, `us`, `cs`, `homcs`, `propto`, and `equalto`

## 4. Conceptual covariance design

### 4.1 General separable structure

The long-term covariance target is:

\[
\Sigma_{\text{group}} = \Sigma_{\text{member}} \otimes \Sigma_{\text{time}}.
\]

This should eventually generalize to:

```text
<memberstructure>x<timestructure>(membertime(member_coord, time_coord) + 0 | groupID)
```

where:

- `member_coord` indexes the member/role position within the group;
- `time_coord` indexes the time point or continuous time value;
- `groupID` is the dyad/triad/group identifier.

### 4.2 Restricted distinguishable separable margin: `un`

For clarity and future extension, a restricted distinguishable-role separable implementation can use an unstructured member margin:

```r
unxar1(...)
unxou(...)
```

For two-member distinguishable dyads, `un` estimates:

\[
\Sigma_{\text{member}} =
\begin{pmatrix}
\sigma_1^2 & \rho_{12}\sigma_1\sigma_2 \\
\rho_{12}\sigma_1\sigma_2 & \sigma_2^2
\end{pmatrix}.
\]

For three-member distinguishable triads, `un` estimates:

\[
\Sigma_{\text{member}} =
\begin{pmatrix}
\sigma_1^2 & \rho_{12}\sigma_1\sigma_2 & \rho_{13}\sigma_1\sigma_3 \\
\rho_{12}\sigma_1\sigma_2 & \sigma_2^2 & \rho_{23}\sigma_2\sigma_3 \\
\rho_{13}\sigma_1\sigma_3 & \rho_{23}\sigma_2\sigma_3 & \sigma_3^2
\end{pmatrix}.
\]

This is flexible for the member covariance at the same occasion and easiest to generalize from dyads to triads and larger small groups. However, with an AR(1) time margin it still imposes a common temporal autocorrelation across members/roles:

\[
\phi_1 = \phi_2 = \cdots = \phi_M.
\]

For two distinguishable members, it also implies symmetric lagged cross-member correlations:

\[
\operatorname{Corr}(A_t, B_{t+k}) =
\operatorname{Corr}(B_t, A_{t+k}) =
\rho_{AB}\phi^k.
\]

This is a useful parsimonious model, but it should not be described as the clean endpoint for distinguishable dyadic ILD.

### 4.3 Exchangeable / indistinguishable groups

For genuinely indistinguishable dyads or triads, an unconstrained `un` member margin is **not** the exchangeable model if the member labels are arbitrary. It would estimate label-specific variances/correlations for arbitrary labels.

Therefore, the broad design should include an exchangeability-constrained version, implemented either as a separate covariance structure or as a constrained special case:

```r
homcsxar1(...)
homcsxou(...)
```

or, if maintainers prefer consistency with existing terminology:

```r
homcsxar1(...)
homcsxou(...)
```

For exchangeable dyads:

\[
\Sigma_{\text{member}} =
\sigma^2
\begin{pmatrix}
1 & \rho \\
\rho & 1
\end{pmatrix}.
\]

For exchangeable triads:

\[
\Sigma_{\text{member}} =
\sigma^2
\begin{pmatrix}
1 & \rho & \rho \\
\rho & 1 & \rho \\
\rho & \rho & 1
\end{pmatrix}.
\]

The exchangeable/homogeneous-CS member margin is the correct model when group members are statistically indistinguishable and member labels are arbitrary.

### 4.4 Distinguishable and indistinguishable dyads in the same model

A future goal is to allow both distinguishable and indistinguishable dyads/groups in the same analysis. This is important for datasets that contain, for example, both mixed-gender distinguishable couples and same-gender exchangeable couples, or both role-defined and role-ambiguous dyads.

Possible future approaches:

1. **Multiple separable terms applied to subsets**
   - One `unx*` term for role-distinguishable groups.
   - One `homcsx*` term for exchangeable groups.
   - Requires a safe way to apply covariance terms to subsets, likely through generated grouping factors or design columns.

2. **Single separable term with group-level pattern metadata**
   - The covariance structure receives a per-group indicator such as `member_pattern = "exchangeable"` vs `"role_defined"`.
   - Internally, the member covariance parameterization is mapped differently by group type.
   - More elegant, but likely much harder to integrate into existing `glmmTMB` covariance parsing.

3. **Formula helper layer**
   - A higher-level helper constructs valid `glmmTMB` formulas and internal factors from variables such as `role_index`, `member_in_dyad`, and `dyad_type`.
   - This may be more realistic than trying to support all mixed distinguishability patterns directly in the covariance syntax.

This should be considered a future extension. The prototype should focus first on exchangeable dyads with `homcsxar1()` and a restricted/common-AR distinguishable model with `unxar1()`.

### 4.5 Priority nonseparable extension for distinguishable dyads: diagonal VAR(1)

For distinguishable dyads, the next priority after the separable v0.0.1 structures should be a nonseparable diagonal VAR(1)-style covariance, tentatively:

```r
var1diag(membertime(role_index, diaryday) + 0 | coupleID)
```

or, if the naming should emphasize covariance rather than process equations:

```r
var1_diag(numFactor(role_index, diaryday) + 0 | coupleID)
```

For two roles:

\[
\begin{aligned}
A_t &= \phi_A A_{t-1} + \epsilon_{A,t}, \\
B_t &= \phi_B B_{t-1} + \epsilon_{B,t},
\end{aligned}
\]

with correlated innovations:

\[
\operatorname{Cov}(\epsilon_{A,t}, \epsilon_{B,t}) \neq 0.
\]

This allows:

- role-specific stationary variances;
- role-specific autoregressive parameters;
- same-time innovation correlation;
- asymmetric lagged cross-role covariance when \(\phi_A \neq \phi_B\).

This is a cleaner distinguishable-dyad covariance model than `unxar1()` because it does not force the same temporal persistence for both roles.

For a stationary diagonal VAR(1), the covariance identities are:

\[
\operatorname{Var}(A_t) = \frac{\sigma_{\epsilon,A}^2}{1-\phi_A^2},
\quad
\operatorname{Var}(B_t) = \frac{\sigma_{\epsilon,B}^2}{1-\phi_B^2},
\]

\[
\operatorname{Cov}(A_t, B_t) =
\frac{\sigma_{\epsilon,AB}}{1-\phi_A\phi_B}.
\]

For positive lags:

\[
\operatorname{Cov}(A_t, B_{t+k}) =
\phi_B^k \operatorname{Cov}(A_t, B_t),
\]

\[
\operatorname{Cov}(B_t, A_{t+k}) =
\phi_A^k \operatorname{Cov}(A_t, B_t).
\]

This structure is still simpler than a full bivariate VAR(1), but it captures the key distinguishable-dyad requirement that role A and role B can have different inertia.

### 4.6 Longer-term full bivariate VAR(1)

A full bivariate VAR(1) covariance would allow cross-lagged dynamics:

\[
\begin{pmatrix}
A_t \\
B_t
\end{pmatrix}
=
\begin{pmatrix}
\phi_{AA} & \phi_{AB} \\
\phi_{BA} & \phi_{BB}
\end{pmatrix}
\begin{pmatrix}
A_{t-1} \\
B_{t-1}
\end{pmatrix}
+
\begin{pmatrix}
\epsilon_{A,t} \\
\epsilon_{B,t}
\end{pmatrix}.
\]

This is the most complete distinguishable dyadic ILD covariance target, but it requires a stable parameterization of the 2 × 2 transition matrix. The eigenvalues of the transition matrix must lie inside the unit circle. This makes it a larger implementation and validation project than `var1diag()`.

## 5. User-facing syntax strategy

### 5.1 Prototype syntax

For the first prototype, use an explicit member-time coordinate helper:

```r
unxar1(membertime(role_index, diaryday) + 0 | coupleID)
```

or, for exchangeable dyads:

```r
homcsxar1(membertime(member_in_dyad, diaryday) + 0 | coupleID)
```

`membertime(member, time)` is a convenience wrapper around a two-dimensional coordinate factor. The covariance structure needs a Cartesian member × time coordinate grid, not additive member and time columns. The shorter alias `mt(member, time)` is equivalent.

### 5.2 Why not `(member + time + 0 | coupleID)`

This is **not** the desired structure:

```r
(member + time + 0 | coupleID)
```

because it creates additive random-effect columns:

```r
member_1, member_2, time_1, time_2, ...
```

The separable covariance needs latent positions indexed by member × time:

```r
member_1_time_1, member_2_time_1, member_1_time_2, member_2_time_2, ...
```

The implementation should therefore use a single coordinate factor that encodes both dimensions.

### 5.3 Longer-term syntax

Manual construction of `numFactor(member, time)` is awkward and error-prone. The user-facing helper hides this:

```r
unxar1(mt(role_index, diaryday) + 0 | coupleID)
```

where `mt()` constructs and validates the member-time coordinate factor.

A more ambitious future syntax could be:

```r
sep(member = un(role_index), time = ar1(diaryday) | coupleID)
```

or:

```r
kron_cov(un(role_index), ar1(diaryday) | coupleID)
```

However, this requires more formula-parser work and should not be the first implementation target.

## 6. Prototype and extension implementation scope

### 6.1 Phase 1: v0.0.1 separable AR1 structures

Implement:

```r
homcsxar1(membertime(member_in_dyad, diaryday) + 0 | coupleID)
unxar1(membertime(role_index, diaryday) + 0 | coupleID)
```

Initial assumptions:

- the first coordinate has two or more member/role levels;
- `diaryday` is ordered, discrete, and unit-spaced;
- The time margin is AR(1): \(\phi^{|t-s|}\).
- `homcsxar1()` uses a homogeneous compound-symmetric member margin and is the clean v0.0.1 target for exchangeable/indistinguishable dyads.
- `unxar1()` uses an unstructured member margin and is a restricted/common-AR distinguishable model, not the full distinguishable endpoint.
- The covariance structure is a latent Gaussian random effect in `glmmTMB`.
- For Gaussian outcomes, the exact residual-covariance version should use `dispformula = ~ 0`.
- For non-Gaussian outcomes, the separable term is interpreted as a latent Gaussian member × time process on the link scale.

For exchangeable dyads, `member_in_dyad` can be derived from `Idiff`, which is coded arbitrarily but stably within dyad:

```r
df$member_in_dyad <- ifelse(df$Idiff < 0, 1L, 2L)
df$member_day <- glmmTMB::numFactor(df$member_in_dyad, df$diaryday)
```

Important: keep `Idiff` and `member_in_dyad` conceptually separate.

- `Idiff` is a signed contrast used in sum/difference random effects.
- `member_in_dyad` is a coordinate index used in the covariance structure.

### 6.2 Phase 2: priority distinguishable dyad extension with diagonal VAR(1)

Implement a nonseparable diagonal VAR(1)-style covariance for distinguishable dyads:

```r
var1diag(membertime(role_index, diaryday) + 0 | coupleID)
```

Initial assumptions:

- exactly two roles for the first implementation;
- unit-spaced discrete time;
- role-specific AR parameters \(\phi_A\), \(\phi_B\);
- unstructured 2 × 2 innovation covariance;
- stationary covariance construction, not conditional likelihood from an observed initial state;
- Gaussian `dispformula` warnings should follow the same residual/latent/nugget logic as the separable structures.

This should be treated as the main path to a clean distinguishable-dyad covariance model.

### 6.3 Phase 3: full bivariate VAR(1)

Implement a full bivariate VAR(1) covariance after `var1diag()` is validated:

```r
var1(membertime(role_index, diaryday) + 0 | coupleID)
```

This adds cross-lagged parameters \(\phi_{AB}\) and \(\phi_{BA}\). The central challenge is a stable parameterization of the 2 × 2 transition matrix. This should not be folded into v0.0.1.

### 6.4 Phase 4: OU / continuous-time exponential margin

Implement:

```r
unxou(membertime(role_index, time_value) + 0 | coupleID)
homcsxou(membertime(member_in_dyad, time_value) + 0 | coupleID)
```

The time margin becomes:

\[
R_{ts} = \exp(-\lambda |t-s|).
\]

This supports irregular time spacing, but it is still a separable covariance model, not a full continuous-time state-space/DSEM model.

### 6.5 Phase 5: generalization to triads and arbitrary small groups

Extend member index from two levels to \(M\) levels.

For distinguishable groups:

```r
unxar1(membertime(role_index, diaryday) + 0 | groupID)
```

For exchangeable groups:

```r
homcsxar1(membertime(member_in_group, diaryday) + 0 | groupID)
```

For triads, `un` estimates three variances and three pairwise correlations in the member margin. `exch` / homogeneous CS estimates one variance and one common member correlation.

### 6.6 Phase 6: general separable covariance algebra

Long-term goal:

Use explicit names that follow the coordinate order:

```text
<memberstructure>x<timestructure>()
```

Potential combinations:

- `unxar1()`
- `unxou()`
- `unxma()`
- `unxarma()`
- `homcsxar1()`
- `homcsxou()`
- `diagxar1()`
- `csxar1()`
- `toepxar1()`

A fully general syntax such as `sep(un(role), ar1(time) | group)` should only be considered after the simpler explicit structures are validated and maintainers are consulted.

## 7. Model formulas we are working toward

### 7.1 Exchangeable dyads: APIM/DIM-style mean model with sum-difference random effects

For exchangeable / indistinguishable dyads, the mean model is APIM-style but can be interpreted equivalently in DIM terms. The random-effects structure includes a dyad-level sum block and a person-within-dyad difference block.

```r
df$member_in_dyad <- ifelse(df$Idiff < 0, 1L, 2L)
df$member_day <- glmmTMB::numFactor(df$member_in_dyad, df$diaryday)

fit_exchangeable <- glmmTMB(
  closeness ~ 1 + 
    diaryday_c +

    # Within-person APIM
    provided_support_actor_cwp + 
    provided_support_partner_cwp +

    # Between-person APIM
    provided_support_actor_cbp +
    provided_support_partner_cbp +

    # Dyad-level intercept and slopes for time-varying predictors
    (1 + diaryday_c +
       provided_support_actor_cwp +
       provided_support_partner_cwp | coupleID) +

    # Both partners' deviations from these dyad-level means and slopes.
    # Separate random-effects block to keep these uncorrelated with the dyad-level block.
    (0 + Idiff +
       I(Idiff * diaryday_c) +
       I(Idiff * provided_support_actor_cwp) +
       I(Idiff * provided_support_partner_cwp) | coupleID) +

    # Hypothetical exchangeable separable residual/latent structure
    homcsxar1(member_day + 0 | coupleID),

  data = df,
  family = gaussian(),

  # For exact Gaussian residual CS(member) x AR1(time):
  dispformula = ~ 0
)
```

If maintainers prefer existing terminology, the term may instead be called:

```r
homcsxar1(member_day + 0 | coupleID)
```

### 7.2 Distinguishable dyads: role-specific APIM with restricted `unxar1()`

For distinguishable dyads, define a role coordinate. Example: male/female dyads.

```r
df$role_index <- ifelse(df$is_male == 1, 1L, 2L)
df$role_day <- glmmTMB::numFactor(df$role_index, df$diaryday)

fit_distinguishable <- glmmTMB(
  closeness ~ 

    # Role-specific intercepts
    0 + is_male + is_female +

    # Role-specific time slopes
    is_male:diaryday_c +
    is_female:diaryday_c +

    # Between-person APIM
    is_male:provided_support_actor_cbp + 
    is_male:provided_support_partner_cbp +
    is_female:provided_support_actor_cbp + 
    is_female:provided_support_partner_cbp +

    # Within-person APIM
    is_male:provided_support_actor_cwp + 
    is_male:provided_support_partner_cwp +
    is_female:provided_support_actor_cwp + 
    is_female:provided_support_partner_cwp +

    # Dyad-level role-specific random effects
    (0 + is_male + is_female +  
       is_male:diaryday_c +
       is_female:diaryday_c +
       is_male:provided_support_actor_cwp + 
       is_male:provided_support_partner_cwp +
       is_female:provided_support_actor_cwp + 
       is_female:provided_support_partner_cwp | coupleID) +

    # Restricted/common-AR distinguishable separable residual/latent structure
    unxar1(role_day + 0 | coupleID),

  data = df,
  family = gaussian(),

  # For exact Gaussian residual UN(member) x AR1(time):
  dispformula = ~ 0
)
```

This model allows role-specific latent variances and same-time role covariance, but it assumes one common AR(1) parameter for both roles. It should be described as a restricted distinguishable model:

\[
\phi_{\text{male}} = \phi_{\text{female}}.
\]

It also imposes symmetric lagged cross-role covariance:

\[
\operatorname{Cov}(\text{male}_t, \text{female}_{t+k})
\propto \phi^k
\]

and

\[
\operatorname{Cov}(\text{female}_t, \text{male}_{t+k})
\propto \phi^k.
\]

For a cleaner distinguishable-dyad covariance model, the next target should be `var1diag()`.

### 7.3 Distinguishable dyads: priority `var1diag()` target

The priority extension for distinguishable dyads is:

```r
df$role_index <- ifelse(df$is_male == 1, 1L, 2L)
df$role_day <- glmmTMB::numFactor(df$role_index, df$diaryday)

fit_distinguishable_var1diag <- glmmTMB(
  closeness ~
    0 + is_male + is_female +
    is_male:diaryday_c +
    is_female:diaryday_c +
    is_male:provided_support_actor_cwp +
    is_male:provided_support_partner_cwp +
    is_female:provided_support_actor_cwp +
    is_female:provided_support_partner_cwp +
    var1diag(role_day + 0 | coupleID),
  data = df,
  family = gaussian(),
  dispformula = ~ 0
)
```

This would allow role-specific temporal persistence:

\[
\phi_{\text{male}} \neq \phi_{\text{female}}.
\]

It is the natural next implementation target before attempting a full bivariate VAR(1) with cross-lagged dynamics.

For Gaussian outcomes, if `dispformula` is not set to `~ 0`, a structured covariance term becomes a latent process plus independent residual/nugget variance. For separable terms:

\[
\Sigma_{\text{total}} = \Sigma_{\text{member}} \otimes \Sigma_{\text{time}} + \sigma_\epsilon^2 I.
\]

This may be appropriate if independent measurement error is intended, but it is not the pure residual Kronecker covariance model.

### 7.4 Gaussian distinguishable dyads with role-specific nugget variance

This can be meaningful:

```r
fit_distinguishable_nugget <- glmmTMB(
  closeness ~ ... +
    unxar1(role_day + 0 | coupleID),
  data = df,
  family = gaussian(),
  dispformula = ~ 0 + is_male + is_female
)
```

Interpretation:

\[
\Sigma_{\text{total}} = \Sigma_{\text{member}} \otimes \Sigma_{\text{time}} + \Sigma_{\text{nugget}}.
\]

where `dispformula = ~ 0 + is_male + is_female` models role-specific independent Gaussian residual/measurement-error variance.

This is a valid model if intentionally chosen. The package should warn that it is no longer a pure separable residual covariance model, not reject it.

## 8. Dispersion-warning design

The package should include interpretation warnings for `dispformula` when separable covariance terms are used. These should be warnings, not errors.

### 8.1 Gaussian family with separable covariance term

If `family = gaussian()` and a member-time product covariance term is present:

- If `dispformula = ~ 0`: no warning.
- If `dispformula` is missing/default `~ 1`: warn.
- If `dispformula = ~ role` or any other non-zero dispersion model: warn, but acknowledge this can be appropriate if the user intentionally wants independent nugget/measurement-error variance.

Suggested warning text for default/non-zero Gaussian dispersion:

```text
A Gaussian model with a separable member-by-time covariance term and a non-zero dispersion model estimates a separable latent covariance process plus independent residual/nugget variance: Sigma_total = Sigma_sep + Sigma_nugget. This may be appropriate if independent measurement error is intended. If your goal is a pure residual CS/UN(member) x AR1/OU(time) covariance structure, set dispformula = ~0.
```

Suggested warning text for role-specific Gaussian dispersion:

```text
You specified a Gaussian separable member-by-time covariance term together with a structured dispersion model. This estimates role-/covariate-specific independent nugget variance in addition to the separable covariance structure. This is valid if intentional, but it is not a pure separable residual covariance model. Use dispformula = ~0 for the pure CS/UN(member) x AR1/OU(time) structure.
```

### 8.2 Families that ignore `dispformula`

For families without a dispersion parameter, such as common binomial/Poisson-type families where `dispformula` is ignored by `glmmTMB`, warn if the user explicitly supplies a non-default `dispformula`.

Suggested warning:

```text
The selected family has no estimated dispersion parameter, so dispformula is ignored by glmmTMB. The separable member-by-time covariance term is still included and remains interpretable as a latent Gaussian structured random effect on the link scale.
```

Do not warn if `dispformula` is not supplied and the family ignores dispersion internally.

### 8.3 Non-Gaussian families with dispersion parameters

For families such as negative binomial, beta, Gamma, Tweedie, etc., `dispformula` controls family-specific dispersion, not an independent Gaussian residual nugget.

If `dispformula = ~ 0` or another unusual dispersion model is specified, warn carefully:

```text
You specified a separable member-by-time covariance term with a non-Gaussian family that has a family-specific dispersion parameter. Unlike Gaussian models, dispformula does not remove a Gaussian residual nugget; it changes the outcome distribution's dispersion model. Check that this is intended. The separable covariance term remains a latent Gaussian structured random effect on the link scale.
```

For negative-binomial models, do **not** recommend automatically setting `dispformula = ~ 0`. In non-Gaussian models, the separable covariance term and the family-specific dispersion parameter have different roles.

### 8.4 No success messages

When the user specifies the recommended pure Gaussian setup correctly, for example:

```r
family = gaussian(), dispformula = ~ 0
```

no informational message should be printed. The package should warn only when a specification is likely to be misunderstood.

## 9. Testing and validation plan

### 9.1 Matrix-construction tests

For each structure, verify that the implied covariance matrix equals the expected Kronecker product.

For `unxar1()`:

\[
\Sigma = \Sigma_{\text{member,UN}} \otimes R_{\text{AR1}}.
\]

Check:

```r
Corr(role1_t, role1_t+k) = phi^k
Corr(role2_t, role2_t+k) = phi^k
Corr(role1_t, role2_t+k) = rho_12 * phi^k
```

For `homcsxar1()`:

```r
Corr(member_i_t, member_j_t+k) = rho_member * phi^k, i != j
```

For OU:

```r
Corr(member_i_t, member_j_s) = rho_ij * exp(-lambda * abs(t - s))
```

### 9.2 Special-case tests

- `rho_member = 0`: independent member-specific time processes.
- `phi = 0`: same-occasion member covariance only.
- `rho_member = 0` and `phi = 0`: independent latent effects.
- exchangeable dyads are invariant to swapping arbitrary member labels.
- distinguishable dyads are not invariant to swapping role labels unless role labels and fixed effects are swapped consistently.

For `var1diag()` once implemented:

- \(\phi_A = \phi_B\) with compatible innovation covariance should reduce to the same lag pattern as a separable common-AR distinguishable model.
- \(\phi_A \neq \phi_B\) should produce asymmetric lagged cross-role covariance.
- \(\phi_A = 0\) or \(\phi_B = 0\) should remove persistence for that role only.
- innovation correlation \(= 0\) should remove same-time innovation covariance while preserving role-specific AR persistence.
- estimates near \(|\phi| = 1\) should produce clear convergence/boundary diagnostics.

### 9.3 Simulation recovery tests

Simulate data under known parameters and fit the model.

Vary:

- number of dyads/groups: 50, 100, 300;
- number of days: 7, 14, 28, 56;
- partner/member correlation: -0.4, 0, 0.3, 0.6;
- temporal correlation: 0, 0.3, 0.7;
- missingness: 0%, 20%, 40%;
- family: Gaussian first, then Poisson and negative binomial.

Evaluate:

- convergence status;
- positive-definite Hessian;
- bias;
- RMSE;
- interval coverage;
- boundary estimates;
- runtime.

### 9.4 Comparison models

Compare the new separable structure against current approximations:

```r
# Additive approximation
homcs(member + 0 | coupleID:diaryday) + ar1(diaryday + 0 | personID)

# Fully unstructured member-time covariance, small T only
us(member_day + 0 | coupleID)

# Fixed known covariance matrix
propto(member_day + 0 | coupleID, V_kron)
```

The new separable structure should recover cross-partner lagged covariance parsimoniously, whereas additive same-day + person-AR terms do not impose:

\[
\operatorname{Corr}(A_t, B_{t+k}) = \rho_{AB}\phi^k.
\]

For distinguishable dyads, add comparisons once `var1diag()` exists:

```r
# Restricted common-AR distinguishable model
unxar1(role_day + 0 | coupleID)

# Priority nonseparable distinguishable model
var1diag(role_day + 0 | coupleID)

# Full unstructured member-time covariance, small T only
us(role_day + 0 | coupleID)
```

The expected improvement of `var1diag()` over `unxar1()` is recovery of role-specific inertia and asymmetric lagged cross-role covariance.

### 9.5 Warning tests

Add tests for warning behavior:

- Gaussian + `unxar1()`/`homcsxar1()` + `dispformula = ~ 0`: no warning.
- Gaussian + `unxar1()`/`homcsxar1()` + missing/default dispersion: warning about nugget variance.
- Gaussian + `unxar1()`/`homcsxar1()` + `dispformula = ~ role`: warning that this is a role-specific nugget and may be intentional.
- Poisson/binomial + explicit `dispformula`: warning that `dispformula` is ignored but the separable term remains interpretable.
- Negative binomial + `dispformula = ~ 0`: warning that dispersion is family-specific and should be checked.

## 10. Documentation plan

Documentation should include:

1. Conceptual vignette: why additive same-day covariance + person-level AR is not equivalent to separable member × time covariance.
2. Dyadic APIM/DIM vignette:
   - exchangeable dyads with `homcsxar1()`;
   - restricted distinguishable dyads with `unxar1()`;
   - priority distinguishable extension with `var1diag()` once implemented.
3. Residual-vs-latent interpretation vignette:
   - Gaussian pure residual covariance with `dispformula = ~ 0`.
   - Gaussian latent process plus nugget with `dispformula != ~ 0`.
   - Non-Gaussian latent Gaussian structured random effect.
4. Simulation vignette: recovery of partner-time covariance.
5. Distinguishable-dyad limitation note:
   - `unxar1()` assumes a common AR parameter across roles.
   - `var1diag()` is the planned clean distinguishable covariance target.
   - full bivariate VAR(1) is the later target for cross-lagged covariance.
6. Developer notes: covariance construction, parameter mapping, row ordering, and missingness handling.

## 11. Open design questions

1. Should exchangeable product structures use `homcsx*` names, or should maintainers prefer a different spelling?
2. Should the first implementation expose only explicit names like `unxar1()` and `homcsxar1()`, or attempt a more general `sep()` syntax?
3. How should mixed distinguishability within the same dataset be represented in formula syntax?
4. Should `var1diag()` be named as a covariance structure (`var1diag`) or as a separable-family extension name (`ar1_diag`, `diag_var1`, etc.)?
5. Should `var1diag()` be parameterized by innovation SDs/correlation plus AR parameters, or by stationary marginal SDs/correlation plus AR parameters?
6. What is the best stable parameterization for the full 2 × 2 VAR(1) transition matrix?
7. Should OU be implemented immediately after AR1, or after the distinguishable `var1diag()` implementation?
8. Can the implementation exploit Kronecker algebra for balanced data, while falling back to observed submatrices for missing/unbalanced data?
9. How should non-Gaussian dispersion warnings be calibrated for families where `dispformula = ~ 0` is technically allowed but rarely intended?
10. Should a small Stan reference implementation be built in parallel for validation, especially for `var1diag()` and full VAR(1)?

## 12. External references to check while implementing

- `glmmTMB` covariance structures vignette: `ar1`, `ou`, `cs`, `homcs`, `us`, `propto`, `equalto`, and examples using `dispformula = ~ 0`.
- `glmmTMB` help for `dispformula`: ignored for families without dispersion; `dispformula = ~ 0` fixes Gaussian residual variance near zero.
- `glmmTMB` hacking/developer vignette: adding a new covariance structure via C++ enum / `termwise_nll` and R-side plumbing.
- Existing `glmmTMB` tests for covariance structures.
- TMB documentation on multivariate Gaussian random effects and parameter transformations.

## 13. Immediate next actions

1. Fork and clone `glmmTMB`.
2. Build and run current `glmmTMB` tests locally.
3. Inspect implementation of `ar1`, `ou`, `homcs`, `cs`, `us`, and `propto`.
4. Create a minimal C++ prototype for `homcsxar1()` and `unxar1()`.
5. Use `membertime(member_or_role_index, diaryday)` as the first prototype coordinate encoding.
6. Simulate small balanced exchangeable and restricted distinguishable dyad datasets and verify covariance recovery.
7. Add warning checks for Gaussian and non-Gaussian `dispformula` behavior.
8. Add explicit documentation that `unxar1()` is not a full distinguishable VAR model.
9. Begin design notes for `var1diag()` as the next priority extension.
10. Open a `glmmTMB` GitHub issue before submitting a PR, explaining:
   - dyadic ILD motivation;
   - why current additive structures are not equivalent;
   - proposed minimal syntax;
   - simulation validation;
   - v0.0.1 scope as exchangeable plus restricted distinguishable separable covariance;
   - priority future extension to `var1diag()` and full bivariate VAR(1);
   - future generalization to triads and arbitrary small groups.

## 14. Pre-PR cleanup

These planning notes are useful while developing the fork, but they should not be included in a `glmmTMB` pull request unless maintainers explicitly ask for a design note.

Before opening a PR:

1. Remove fork-local planning files from the final patch, including:
   - `notes/implementation_plan.md`
   - `notes/implementation_plan_v0.0.1.md`
2. Decide whether exploratory scripts in `notes/` should be removed, kept only locally, or converted into official tests/vignettes.
3. Move any user-facing material into the appropriate package files:
   - concise help text in `glmmTMB/R/glmmTMB.R`
   - generated Rd files in `glmmTMB/man/`
   - tests in `glmmTMB/tests/testthat/`
   - a vignette only if the example is stable and maintainer-facing scope allows it.
4. Use the plans as source material for the GitHub issue or PR description rather than adding them as package files.
