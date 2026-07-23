# Complete Code Walkthrough: `us() x ar1()` Separable Term

This file is a review helper for one concrete model term:

```r
separable(us(0 + member) %x% ar1(0 + time) | group)
```

It follows that term through the relevant glmmTMB runtime path, from the first
place this branch touches the formula to the last place the separable term is
handled for likelihood/reporting.

Markers:

```text
[NEW/MODIFIED]  Code added or changed by this branch.
[EXISTING]      Pre-existing glmmTMB/reformulas/TMB path that this branch uses.
```

This walkthrough is intentionally concrete.  It is not a generic separable
covariance design document.  The generic review guide is:

```text
notes/separable-kron-work/separable-kron-code-review-guide.md
```

## How To Use This File

If you are not familiar with glmmTMB internals, read this file in this order:

```text
1. Read the "One-Screen Mental Model" below.
2. Read "End-To-End Data Shape Summary" near the end.
3. Then read sections 1 through 18 in order.
4. For each section, open the named source file only if you want to verify the
   quoted code against the actual branch.
5. Use Appendix A to check small support hooks that are not part of the main fit
   path.
```

You do not need to understand all of lme4, reformulas, or TMB before reading
this.  The important thing is to track one object as it changes shape:

```text
public formula
  -> rewritten formula
  -> model-frame sepgrid factor
  -> random-effect design metadata
  -> R covariance metadata
  -> C++ per-term metadata
  -> theta parsing
  -> likelihood evaluation
  -> reporting data
```

## One-Screen Mental Model

The user writes:

```r
separable(us(0 + member) %x% ar1(0 + time) | group)
```

glmmTMB cannot directly build a random-effect design matrix from that product
syntax.  So this branch compiles it to:

```r
separable(sepgrid(member, time) + 0 | group, list(...metadata...))
```

The two pieces have different jobs:

```text
sepgrid(member, time) + 0 | group
  Builds the actual latent random-effect columns:
  one column per member-time cell per group.

list(...metadata...)
  Tells getReStruc()/C++ what covariance model those columns should have:
  us for member, ar1 for time, scale on the member margin.
```

For a dataset with:

```text
2 members
3 time points
2 groups
```

each group has this latent vector:

```text
m1_t1, m2_t1, m1_t2, m2_t2, m1_t3, m2_t3
```

and C++ evaluates:

```text
member covariance  x  AR(1) time correlation
```

without building the full dense covariance matrix during likelihood evaluation.

## Toy Example Carried Through The Shapes

If you only understand one section before reading the code, make it this one.
It shows the same object at each layer with small numbers.

This is a small concrete example:

```text
M = 2 members: m1, m2
T = 3 time points: t1, t2, t3
G = 2 groups: g1, g2
```

Translation:

```text
There are six latent random effects per group:

2 members * 3 times = 6 latent cells

With two groups, the separable random-effect term has two independent
six-cell blocks, one for g1 and one for g2.
```

The user writes:

```r
y ~ 1 + separable(us(0 + member) %x% ar1(0 + time) | group)
```

### Public Formula To Internal Formula

The public separable term:

```r
separable(us(0 + member) %x% ar1(0 + time) | group)
```

is compiled to the internal shape:

```r
separable(
  sepgrid(member, time) + 0 | group,
  list(
    grid = c("member", "time"),
    margins = data.frame(
      struc = c("us", "ar1"),
      var = c("member", "time")
    )
  )
)
```

What this means:

```text
sepgrid(member, time) + 0 | group
  tells existing random-effect machinery which latent effects exist.

list(...)
  tells separable-specific code what covariance those latent effects have.
```

No statistics happens in this rewrite.  It is just a compilation step:

```text
user-friendly syntax -> internal design syntax plus metadata
```

### `sepgrid()` Levels

With member first and time second, the first coordinate runs fastest:

```text
level 1: (1,1) = m1_t1
level 2: (2,1) = m2_t1
level 3: (1,2) = m1_t2
level 4: (2,2) = m2_t2
level 5: (1,3) = m1_t3
level 6: (2,3) = m2_t3
```

The same ordering as a table:

| flat row `k` | member index | time index | cell label |
|---:|---:|---:|---|
| 1 | 1 | 1 | `m1_t1` |
| 2 | 2 | 1 | `m2_t1` |
| 3 | 1 | 2 | `m1_t2` |
| 4 | 2 | 2 | `m2_t2` |
| 5 | 1 | 3 | `m1_t3` |
| 6 | 2 | 3 | `m2_t3` |

So each within-group random-effect vector is ordered as:

```text
u = [u_m1_t1, u_m2_t1, u_m1_t2, u_m2_t2, u_m1_t3, u_m2_t3]
```

Read this as:

```text
first all members at time 1,
then all members at time 2,
then all members at time 3.
```

The flat index rule is:

```text
k = member_index + M * time_index
```

using zero-based C++ indices:

```text
k = i0 + 2 * i1
```

### Metadata From `getReStruc()`

For this toy example:

```text
sepDims         = c(2, 3)
sepCodes        = c(us, ar1)
sepDensityKinds = c(dense_corr, ar1)
sepDispatch     = dense_ar1
sepScaleMode    = margin
sepScaleSpec    = 0
blockSize       = 6
blockReps       = 2
```

`sepScaleSpec = 0` means the first coordinate, `member`, supplies the SDs.

Plain English:

```text
Each group-specific random-effect vector has 6 entries.
The first grid direction is member.
The second grid direction is time.
The member direction gets an unstructured correlation.
The time direction gets AR(1).
The member direction also supplies the standard deviations.
```

### Theta Layout

For `M = 2`, `us(member)` has:

```text
2 log-SDs
1 member-correlation parameter
```

`ar1(time)` has:

```text
1 AR(1) parameter
```

So the separable theta vector is:

```text
theta = [
  log_sd_m1,
  log_sd_m2,
  theta_member_cor_m1_m2,
  theta_time_ar1
]
```

As a review table:

| theta slot | consumed by | transformed object | role |
|---:|---|---|---|
| `theta[0]` | `us(member)` | `sd_m1 = exp(theta[0])` | member 1 SD |
| `theta[1]` | `us(member)` | `sd_m2 = exp(theta[1])` | member 2 SD |
| `theta[2]` | `us(member)` | `R_member[1,2]` through `UNSTRUCTURED_CORR_t` | member correlation |
| `theta[3]` | `ar1(time)` | `phi = theta[3] / sqrt(1 + theta[3]^2)` | AR(1) time correlation |

Do not confuse `theta` with the random effects `u`.

```text
theta:
  covariance parameters estimated by the optimizer.

u:
  latent random-effect values for the groups.
```

C++ transforms this to:

```text
sd_m1 = exp(log_sd_m1)
sd_m2 = exp(log_sd_m2)
rho_member = interpreted by TMB UNSTRUCTURED_CORR_t
phi_time = theta_time_ar1 / sqrt(1 + theta_time_ar1^2)
```

For a two-level `us()` margin, the unstructured correlation has only one
correlation parameter.  Conceptually it is the member correlation:

```text
R_member =
  [ 1    rho ]
  [ rho  1   ]
```

The time margin is:

```text
R_time =
  [ 1       phi     phi^2 ]
  [ phi     1       phi   ]
  [ phi^2   phi     1     ]
```

### C++ `U` Shape

`allterms_nll()` passes C++ an array `U` with:

```text
rows    = blockSize = 6
columns = blockReps = 2
```

For group `g1`, the first column is:

```text
U[, 1] =
  row 1: u_g1_m1_t1
  row 2: u_g1_m2_t1
  row 3: u_g1_m1_t2
  row 4: u_g1_m2_t2
  row 5: u_g1_m1_t3
  row 6: u_g1_m2_t3
```

For group `g2`, the second column has the same within-group ordering:

```text
U[, 2] =
  row 1: u_g2_m1_t1
  row 2: u_g2_m2_t1
  row 3: u_g2_m1_t2
  row 4: u_g2_m2_t2
  row 5: u_g2_m1_t3
  row 6: u_g2_m2_t3
```

The `U` object as a matrix:

| row | cell within group | `U[, 1]`: group `g1` | `U[, 2]`: group `g2` |
|---:|---|---|---|
| 1 | `m1_t1` | `u_g1_m1_t1` | `u_g2_m1_t1` |
| 2 | `m2_t1` | `u_g1_m2_t1` | `u_g2_m2_t1` |
| 3 | `m1_t2` | `u_g1_m1_t2` | `u_g2_m1_t2` |
| 4 | `m2_t2` | `u_g1_m2_t2` | `u_g2_m2_t2` |
| 5 | `m1_t3` | `u_g1_m1_t3` | `u_g2_m1_t3` |
| 6 | `m2_t3` | `u_g1_m2_t3` | `u_g2_m2_t3` |

Why two columns?

```text
G = 2 in this toy example.
If G = 10, U would have 6 rows and 10 columns.
Each column is one group's member-time random-effect vector.
```

### C++ Reshape To `z`

`separable_dense_ar1_nll()` reshapes that flat vector to:

```text
z[member, time]
```

and divides by the member SD:

```text
z[m1, t1] = u_m1_t1 / sd_m1
z[m2, t1] = u_m2_t1 / sd_m2
z[m1, t2] = u_m1_t2 / sd_m1
z[m2, t2] = u_m2_t2 / sd_m2
z[m1, t3] = u_m1_t3 / sd_m1
z[m2, t3] = u_m2_t3 / sd_m2
```

For one group column, the reshaped `z` array looks like:

| `z[member, time]` | `t1` | `t2` | `t3` |
|---|---|---|---|
| `m1` | `u_m1_t1 / sd_m1` | `u_m1_t2 / sd_m1` | `u_m1_t3 / sd_m1` |
| `m2` | `u_m2_t1 / sd_m2` | `u_m2_t2 / sd_m2` | `u_m2_t3 / sd_m2` |

That table is the same data as the flat vector, just viewed as a
member-by-time array.  The C++ code constructs this array separately for each
group column in `U`.

This standardization means:

```text
U is on the original random-effect scale.
z is on the correlation scale.
```

After dividing by SDs, the remaining covariance structure is just correlation:

```text
member correlation x time correlation
```

The scale adjustment added to the negative log-likelihood is:

```text
log(sd_m1) + log(sd_m2) +
log(sd_m1) + log(sd_m2) +
log(sd_m1) + log(sd_m2)
```

or:

```text
3 * log(sd_m1) + 3 * log(sd_m2)
```

### C++ Separable Density Call

Because the array dimensions are:

```text
dimension 0 = member
dimension 1 = time
```

TMB receives the density arguments in reverse order:

```cpp
SEPARABLE(AR1(phi_time), UNSTRUCTURED_CORR(member_cor))(z)
```

This evaluates the standardized covariance:

```text
R_time %x% R_member
```

without constructing that full matrix in the likelihood path.

With `M = 2` and `T = 3`, the standardized covariance has the block form:

| time block | covariance block |
|---|---|
| `t1` with `t1` | `1 * R_member` |
| `t1` with `t2` | `phi * R_member` |
| `t1` with `t3` | `phi^2 * R_member` |
| `t2` with `t2` | `1 * R_member` |
| `t2` with `t3` | `phi * R_member` |
| `t3` with `t3` | `1 * R_member` |

This is what `R_time %x% R_member` means for the current storage order: each
time-time entry multiplies the full member correlation matrix.

The reversal is important:

```text
array is indexed as z[member, time]
SEPARABLE gets densities as time first, member second
```

That is a TMB convention.  The code has tests because this is easy to get wrong.

### Full Correlation Only For Reporting

For this toy example, the reported full correlation matrix is:

|  | `m1_t1` | `m2_t1` | `m1_t2` | `m2_t2` | `m1_t3` | `m2_t3` |
|---|---|---|---|---|---|---|
| `m1_t1` | `1` | `rho` | `phi` | `rho*phi` | `phi^2` | `rho*phi^2` |
| `m2_t1` | `rho` | `1` | `rho*phi` | `phi` | `rho*phi^2` | `phi^2` |
| `m1_t2` | `phi` | `rho*phi` | `1` | `rho` | `phi` | `rho*phi` |
| `m2_t2` | `rho*phi` | `phi` | `rho` | `1` | `rho*phi` | `phi` |
| `m1_t3` | `phi^2` | `rho*phi^2` | `phi` | `rho*phi` | `1` | `rho` |
| `m2_t3` | `rho*phi^2` | `phi^2` | `rho*phi` | `phi` | `rho` | `1` |

This full matrix is useful for tests and `VarCorr()` output, but the likelihood
uses `SEPARABLE()` instead.

Noob check:

```text
Same member, adjacent time:
  corr(m1_t1, m1_t2) = phi

Different member, same time:
  corr(m1_t1, m2_t1) = rho

Different member, adjacent time:
  corr(m1_t1, m2_t2) = rho * phi

Different member, two time steps apart:
  corr(m1_t1, m2_t3) = rho * phi^2
```

## Toy Trace Through Every Runtime Step

The detailed sections below explain the actual code.  This trace keeps the toy
example attached to every major runtime step so the reader does not have to
reconstruct it.

```text
Toy setup:
  M = 2 members
  T = 3 times
  G = 2 groups

Public term:
  separable(us(0 + member) %x% ar1(0 + time) | group)
```

Step 1, `glmmTMB()` formula rewrite:

```text
Before:
  separable(us(0 + member) %x% ar1(0 + time) | group)

After:
  separable(sepgrid(member, time) + 0 | group, list(...metadata...))

Stored call:
  still shows the user-facing product syntax.
```

Step 2, product parser:

```text
us(0 + member)
  -> struc = "us", var = "member"

ar1(0 + time)
  -> struc = "ar1", var = "time"

grid variables:
  c("member", "time")
```

Step 3, `sepgrid()` model-frame column:

```text
sepgrid(member, time) has 6 levels:
  (1,1), (2,1), (1,2), (2,2), (1,3), (2,3)

These correspond to:
  m1_t1, m2_t1, m1_t2, m2_t2, m1_t3, m2_t3
```

Step 4, model-frame level preservation:

```text
If t2 were globally unobserved but still a factor level, sepgrid() would still
keep (1,2) and (2,2).

That matters because t1 to t3 should be a two-step AR(1) lag, not a one-step
lag.
```

Step 5, `mkReTrms()` design construction:

```text
The random-effect design has 6 columns for this separable term:
  one column per member-time cell.

For G = 2, there are two independent group blocks.
For G = 10, the same 6-column block is repeated across 10 groups.
```

Step 6, `splitForm()` metadata alignment:

```text
reTrmClasses[i]:
  "separable"

reTrmFormulas[[i]]:
  sepgrid(member, time) + 0 | group

reTrmAddArgs[[i]]:
  list(grid = c("member", "time"),
       margins = data.frame(struc = c("us", "ar1"),
                            var = c("member", "time")))

The important review point:
  the metadata stays attached to this exact separable term.
```

Step 7, `getReStruc()` grid parsing:

```text
cnms:
  "(1,1)", "(2,1)", "(1,2)", "(2,2)", "(1,3)", "(2,3)"

parseNumLevels(cnms):
  matrix with rows:
    1 1
    2 1
    1 2
    2 2
    1 3
    2 3

sepDims:
  c(2, 3)

blockSize:
  6
```

Step 8, `getReStruc()` theta count:

```text
us(member):
  2 log-SDs
  1 correlation parameter

ar1(time):
  1 correlation parameter

blockNumTheta:
  2 + 1 + 1 = 4
```

Step 9, R metadata sent to C++:

```text
sepDims         = c(2, 3)
sepCodes        = c(us_covstruct, ar1_covstruct)
sepDensityKinds = c(dense_corr_sep, ar1_sep)
sepDispatch     = dense_ar1_dispatch
sepScaleMode    = margin_sep_scale
sepScaleSpec    = 0
```

Step 10, C++ `per_term_info`:

```text
blockCode     = separable_covstruct
blockSize     = 6
blockReps     = 2
blockNumTheta = 4

sep* fields:
  copied from R into C++ vectors.
```

Step 11, `allterms_nll()` slicing:

```text
U has shape:
  6 rows x 2 columns

theta segment has length:
  4
```

Step 12, C++ metadata check:

```text
sepDims length:
  2

sepDims product:
  2 * 3 = 6 = blockSize

scale mode:
  margin

scale margin:
  coordinate 0 = member
```

Step 13, C++ theta parsing:

```text
theta[0] = log_sd_m1
theta[1] = log_sd_m2
theta[2] = member correlation parameter
theta[3] = AR(1) parameter

sd = c(exp(theta[0]), exp(theta[1]))
phi = theta[3] / sqrt(1 + theta[3]^2)
```

Step 14, dense member density:

```text
sep.dense_code = us_covstruct

C++ builds:
  density::UNSTRUCTURED_CORR_t<Type> dense_density(sep.us_corr_params)
```

Step 15, C++ reshape and scaling:

```text
Flat U column for each group:
  u_m1_t1, u_m2_t1, u_m1_t2, u_m2_t2, u_m1_t3, u_m2_t3

z array:
  z[m1,t1] = u_m1_t1 / sd_m1
  z[m2,t1] = u_m2_t1 / sd_m2
  z[m1,t2] = u_m1_t2 / sd_m1
  z[m2,t2] = u_m2_t2 / sd_m2
  z[m1,t3] = u_m1_t3 / sd_m1
  z[m2,t3] = u_m2_t3 / sd_m2

logscale:
  3*log(sd_m1) + 3*log(sd_m2)

The C++ loop repeats this reshape/scaling/density evaluation once per group
column.  With G = 2, the same per-group logscale contribution is added twice,
once for g1 and once for g2.
```

Step 16, C++ separable likelihood:

```text
Because dimensions are [member, time], TMB gets densities as:

  SEPARABLE(AR1(phi), dense_density)(z)

This evaluates:
  R_time %x% R_member

without building that full matrix in the likelihood path.
```

Step 17, reporting:

```text
term.sd:
  sd_m1, sd_m2, sd_m1, sd_m2, sd_m1, sd_m2

term.corr if fullCor == 1:
  6 x 6 Kronecker correlation matrix in member-fastest order.
```

Step 18, current limitations:

```text
simulate(fit):
  errors for separable terms.

predict(fit, newdata = ...):
  errors for separable terms.

Reason:
  those paths need extra grid encoding/simulation contracts that are not
  implemented yet.
```

## Common Names In This Walkthrough

```text
M
  number of member levels.

T
  number of time levels.

G
  number of group levels.

blockSize
  number of random effects for one group = M * T.

blockReps
  number of repeated blocks = G.

cnms
  random-effect column names created by mkReTrms(); for sepgrid these are
  parseable labels like "(1,1)", "(2,1)", "(1,2)", ...

theta
  unconstrained covariance parameters optimized by TMB.

U
  C++ view of random effects for one term:
  rows = blockSize, columns = blockReps.

z
  standardized two-dimensional array used by TMB SEPARABLE().
```

## What A Reviewer Should Be Able To Check

After reading this file, a reviewer should be able to answer:

```text
Does the public syntax become the intended complete grid?
Does the complete grid preserve missing intended time levels?
Does each separable term keep its own metadata?
Does getReStruc() count theta correctly for us x ar1?
Does C++ consume theta in the same order R counted it?
Does C++ reshape the random effects in the same order sepgrid created them?
Does SEPARABLE() get the dimensions in the right order?
Is scale applied exactly once per latent cell?
Is full Kronecker construction avoided in the likelihood path?
```

## Source Comments Mirrored In This Guide

The source files now contain short comments at the highest-risk handoff points.
This guide expands those comments rather than copying them literally.

| source comment topic | where it appears in source | where this guide expands it |
|---|---|---|
| `aa[[i]]` carries separable metadata, not ordinary rank metadata | `glmmTMB/R/glmmTMB.R`, `getReStruc()` | Sections 5, 8, and "Why The Current Implementation Is Hacky But Reviewable" |
| separable theta count is margin-based, not flattened-block-based | `glmmTMB/R/glmmTMB.R`, `parFun()` | "Theta Layout", Step 8, and Section 8 |
| `cnms` are parseable `sepgrid()` cell labels | `glmmTMB/R/utils_covstruct.R`, `.sep_restruc_info()` | "`sepgrid()` Levels", Step 7, and Section 8 |
| theta is counted and consumed in product/margin order | `glmmTMB/R/utils_covstruct.R`, `.sep_restruc_info()` and `glmmTMB/src/glmmTMB.cpp`, `parse_separable_dense_ar1()` | "Theta Layout", Step 13, and Sections 13-14 |
| `U` is `blockSize x blockReps`: cells by groups | `glmmTMB/src/glmmTMB.cpp`, `separable_dense_ar1_nll()` | "C++ `U` Shape", Step 11, and Section 15 |
| each group column is reshaped to `z[member, time]` | `glmmTMB/src/glmmTMB.cpp`, `separable_dense_ar1_nll()` | "C++ Reshape To `z`" and Section 15 |
| reported correlation is margin-0 correlation times margin-1 correlation | `glmmTMB/src/glmmTMB.cpp`, `report_separable_dense_ar1()` | "Full Correlation Only For Reporting" and Section 16 |
| `kronecker()` order follows `sepgrid()` storage order | `glmmTMB/tests/testthat/test-separable.R`, `make_sep_case()` | "C++ Separable Density Call", "Full Correlation Only For Reporting", and Appendix B |
| dense likelihood comparison uses two group blocks | `glmmTMB/tests/testthat/test-separable.R`, `expect_separable_dense_nll()` | "C++ `U` Shape", Step 11, and Appendix B |

If the code comments feel too compressed, use the table above as an index: the
source says the rule, and this file shows the concrete two-member, three-time,
two-group example behind the rule.

Notation in later code blocks:

```text
Comments copied from source:
  These are ordinary source comments shown inside the code excerpt.

Guide-only notes:
  These appear outside the source excerpts, usually as tables or text blocks.
  They are not meant to be copied into R/C++ files; they are here only to make
  the walkthrough easier to read.
```

## Rationale Map

This table explains why each major step exists and why the branch implements it
in the current way.  The detailed sections below show the code.

```text
Step: Rewrite public product syntax early.
Why:
  model.frame() and mkReTrms() must see a concrete random-effect design before
  getReStruc() runs.
Why this way:
  Rewriting to sepgrid(...) lets existing glmmTMB/lme4-style machinery build
  the complete member-time design without a new random-effect design compiler.
Tradeoff:
  The formula is temporarily used as both design syntax and metadata transport.
Cleaner future:
  A parser API could return structured nested covariance metadata directly while
  still giving mkReTrms() a design expression.

Step: Use sepgrid(member, time).
Why:
  The separable latent vector needs one effect for every intended member-time
  cell, including globally unobserved intended cells.
Why this way:
  A factor with complete Cartesian-product levels fits the existing model-matrix
  and random-effect machinery.
Tradeoff:
  It only supports complete indicator-grid margins, not arbitrary marginal
  random-slope designs.
Cleaner future:
  Build marginal design matrices explicitly and combine them row-wise.

Step: Preserve sepgrid levels in model.frame().
Why:
  Dropping globally unobserved time levels would compress AR(1) distances.
Why this way:
  Temporarily disable level dropping, then drop ordinary factor levels manually.
Tradeoff:
  This is a narrow special case in glmmTMB() model-frame construction.
Cleaner future:
  Store intended grid metadata separately and make grid preservation an explicit
  part of separable term construction.

Step: Carry metadata through reTrmAddArgs as list(...).
Why:
  getReStruc() needs to know which margin is us, which is ar1, which margin
  carries scale, and how many theta parameters to allocate.
Why this way:
  reformulas::splitForm() already keeps extra arguments attached to the exact
  random-effect term, so multiple separable terms remain aligned.
Tradeoff:
  The generated list(...) call is not a beautiful public representation.
Cleaner future:
  reformulas/glmmTMB could return first-class structured covariance metadata.

Step: Validate and count theta in .sep_restruc_info().
Why:
  C++ should receive compact, already-validated metadata rather than re-parsing
  R formulas.
Why this way:
  R has the formula context, grid column names, and registry of supported
  margins; C++ should focus on numeric likelihood evaluation.
Tradeoff:
  The R and C++ enum/code tables must stay synchronized.
Cleaner future:
  Generate or centralize shared codes so R/C++ mappings cannot drift.

Step: Use dense_ar1 dispatch.
Why:
  The current implementation supports only one dense correlation margin crossed
  with one AR(1) margin.
Why this way:
  Dispatch separates "what evaluator is needed" from coordinate order, so
  us x ar1 and ar1 x us can share one evaluator.
Tradeoff:
  Adding dense x dense or ar1 x ar1 still requires new C++ branches.
Cleaner future:
  Add more dispatch entries and margin density builders behind the same
  metadata contract.

Step: Put scale on the dense margin.
Why:
  If both margins had arbitrary scale, the Kronecker product would be
  non-identifiable.
Why this way:
  For the motivating member-time model, member-specific SDs plus time
  correlation-only AR(1) are identifiable and useful.
Tradeoff:
  No global/product/cell scale modes yet.
Cleaner future:
  Add a general cell_sd builder supporting global(), product(), and cell().

Step: Parse theta in C++ in margin order.
Why:
  The optimizer gives C++ one flat theta vector.
Why this way:
  R counts parameters in margin order, and C++ consumes them in the same order,
  including reversed margin syntax.
Tradeoff:
  R and C++ must agree exactly on parameter counts and ordering.
Cleaner future:
  More shared helpers/tests around theta layout for each dispatch.

Step: Use TMB SEPARABLE().
Why:
  It evaluates the Kronecker/separable Gaussian density without building and
  factorizing the full covariance matrix.
Why this way:
  The dense member margin and AR(1) time margin are already available as TMB
  density objects.
Tradeoff:
  TMB expects SEPARABLE arguments in reverse array-dimension order, which is an
  easy place to make mistakes.
Cleaner future:
  Keep explicit ordering tests and possibly wrap the order reversal in a helper.

Step: Scale manually in C++.
Why:
  SEPARABLE() is evaluated on standardized correlation-scale random effects.
Why this way:
  Manual division by SD and addition of sum(log SD) is explicit and works for
  the current margin-scale mode.
Tradeoff:
  The loop knows too much about the current scale mode.
Cleaner future:
  Build a cell_sd vector once, then use either the same manual scaling or
  TMB's VECSCALE() after ordering is well-tested.

Step: Build full correlation only for reporting/testing.
Why:
  VarCorr() and tests need reported SD/correlation objects.
Why this way:
  Constructing the full matrix outside the likelihood path makes correctness
  easy to check without slowing optimization.
Tradeoff:
  Full printed matrices are not suitable for long time series.
Cleaner future:
  Report compact separable components by default, with optional full expansion.

Step: Disable simulation and predict(newdata).
Why:
  Those paths need additional contracts for drawing/writing separable arrays and
  encoding new data against the fitted grid.
Why this way:
  Explicit errors are safer than silently using a wrong grid or layout.
Tradeoff:
  Some downstream workflows are unavailable for now.
Cleaner future:
  Implement simulation via standardized separable arrays and implement newdata
  encoding against stored fitted grid metadata.
```

## Example Term

Assume the user fits something like:

```r
glmmTMB(
  y ~ 1 + separable(us(0 + member) %x% ar1(0 + time) | group),
  data = dd
)
```

The intended latent random-effect block for each `group` is a complete
`member x time` grid:

```text
member1_time1, member2_time1, ..., memberM_time1,
member1_time2, member2_time2, ..., memberM_time2,
...
```

The covariance is:

```text
D_member * (R_time %x% R_member) * D_member
```

where `us(member)` supplies member SDs plus member correlation parameters, and
`ar1(time)` supplies only the AR(1) correlation parameter.

## 1. `glmmTMB()` First Rewrites The Public Formula

Location:

```text
glmmTMB/R/glmmTMB.R
```

### Existing Context

`glmmTMB()` sets the formula environment and stores the call before constructing
the model frame and random-effect structures.

### New/Modified Code

Local explanation:

```text
What this code does:
  Saves the user-facing formula, then rewrites the internal formula.

Why it exists:
  Downstream model-frame and random-effect construction need an ordinary
  grid-shaped random-effect term, not the public product-margin syntax.

Why this implementation:
  Keeping user_formula in the stored call preserves readable user-facing output,
  while formula gets the internal representation needed by glmmTMB internals.
```

```r
environment(formula) <- parent.frame()
user_formula <- formula
formula <- rewrite_separable_formula(formula)
call$formula <- mc$formula <- user_formula
```

The same pattern is applied to `ziformula` and `dispformula`:

```r
environment(ziformula) <- environment(formula)
user_ziformula <- ziformula
ziformula <- rewrite_separable_formula(ziformula)
call$ziformula <- user_ziformula

environment(dispformula) <- environment(formula)
user_dispformula <- dispformula
dispformula <- rewrite_separable_formula(dispformula)
call$dispformula <- user_dispformula
```

### What Happens

The formula used internally is rewritten, but the stored call remains the
user-facing product syntax.  This avoids exposing generated metadata calls in
`fit$call`.

For the example term, the internal formula becomes conceptually:

```r
separable(
  sepgrid(member, time) + 0 | group,
  list(
    grid = c("member", "time"),
    margins = data.frame(struc = c("us", "ar1"),
                         var = c("member", "time"),
                         stringsAsFactors = FALSE)
  )
)
```

The first argument is the random-effect design expression that existing lme4-ish
machinery can turn into a complete grid design.  The second argument is
term-local metadata carried through `reformulas::splitForm()` as
`reTrmAddArgs`.

### Why This Is Done Here

This must happen before `model.frame()` and `mkReTrms()` because those existing
functions need to see `sepgrid(member, time) + 0 | group` in order to build the
right factor levels, model matrix columns, `cnms`, and block size.

Doing this only inside `getReStruc()` would be too late: the random-effect design
matrix would already have been built.

## 2. R Formula Rewrite Helpers

Location:

```text
glmmTMB/R/utils_covstruct.R
```

### 2.1 `sepgrid()` Builds The Complete Latent Grid

Marker: `[NEW/MODIFIED]`

Local explanation:

```text
What this code does:
  Converts member/time coordinates into one factor with complete Cartesian
  product levels.

Why it exists:
  The latent separable random-effect block needs every intended member-time
  cell, even if some cells are not observed in the response data.

Why this implementation:
  A factor with encoded numeric levels works with existing model.matrix(),
  mkReTrms(), and parseNumLevels() machinery.
```

```r
sepgrid <- function(x, ...) {
    y <- data.frame(x, ...)

    ok <- vapply(y, function(z) is.numeric(z) || is.factor(z) ||
                   is.character(z) || is.integer(z), logical(1))
    if (!all(ok))
        stop("All arguments to 'sepgrid' must be numeric, factor, integer, or character.")

    levs <- lapply(y, function(z) {
        if (is.factor(z)) levels(z) else sort(unique(z[!is.na(z)]))
    })
    vals <- Map(function(z, lev) match(if (is.factor(z)) as.character(z) else z, lev),
                y, levs)
    vals <- as.data.frame(vals)

    asChar <- function(y) {
        is_na <- !stats::complete.cases(y)
        y <- lapply(y, as.character)
        ans <- do.call("paste", c(y, list(sep=",")))
        ans <- paste0("(", ans, ")")
        ans[is_na] <- NA_character_
        ans
    }

    grid <- do.call(expand.grid, c(lapply(levs, seq_along),
                                   list(KEEP.OUT.ATTRS = FALSE)))
    factor(asChar(vals), levels = asChar(grid))
}
```

Explanation:

```text
Input:
  member, time

Output:
  a factor whose levels are the complete Cartesian product of member levels and
  time levels.

Important:
  factor inputs preserve unused levels.  This is what keeps globally unobserved
  intended time points in the AR(1) grid.

Ordering:
  expand.grid() varies the first coordinate fastest, so the C++ flat index is
  k = member_index + n_member * time_index.
```

This function is evaluated by the ordinary model-frame/model-matrix path after
the formula rewrite.

### 2.2 Small Call-Tree Helpers

Marker: `[NEW/MODIFIED]`

Local explanation:

```text
What this code does:
  Gives the parser stable one-line call names and deparsed fragments.

Why it exists:
  The separable parser walks unevaluated formula calls, not evaluated objects.

Why this implementation:
  Small local helpers avoid introducing a broader formula parsing dependency for
  the prototype.
```

```r
.sep_deparse <- function(x) deparse1(x, collapse = "", width.cutoff = 500L)

.sep_call_name <- function(x) {
    if (!is.call(x)) return(NULL)
    .sep_deparse(x[[1]])
}
```

Explanation:

```text
.sep_deparse()
  Gives a stable one-line text representation of a formula/call fragment.

.sep_call_name()
  Returns the call head, e.g. "separable", "%x%", "us", or "ar1".
```

These are deliberately small because the parser is local to this prototype.

### 2.3 Finding Generated `sepgrid()` Columns

Marker: `[NEW/MODIFIED]`

Local explanation:

```text
What this code does:
  Finds every sepgrid(...) call in rewritten formulas and records the generated
  model-frame column names.

Why it exists:
  glmmTMB must preserve unused levels for sepgrid columns while preserving
  ordinary factor-dropping behavior elsewhere.

Why this implementation:
  Searching the rewritten formula is enough because sepgrid(...) is explicit by
  this point and has the same deparsed name that model.frame() uses.
```

```r
.sep_find_call <- function(x, name) {
    if (!is.call(x)) return(NULL)
    if (identical(.sep_call_name(x), name)) return(x)
    for (i in seq_along(x)[-1]) {
        ans <- .sep_find_call(x[[i]], name)
        if (!is.null(ans)) return(ans)
    }
    NULL
}

.sep_find_calls <- function(x, name) {
    if (!is.call(x)) return(list())
    ans <- if (identical(.sep_call_name(x), name)) list(x) else list()
    for (i in seq_along(x)[-1]) {
        ans <- c(ans, .sep_find_calls(x[[i]], name))
    }
    ans
}

.sepgrid_colnames <- function(...) {
    forms <- list(...)
    calls <- unlist(lapply(forms, function(f) {
        if (!inherits(f, "formula")) return(list())
        .sep_find_calls(f[[length(f)]], "sepgrid")
    }), recursive = FALSE)
    unique(vapply(calls, .sep_deparse, character(1)))
}
```

Explanation:

```text
After rewrite, the formula contains sepgrid(member, time).  glmmTMB needs to
know which model-frame columns came from sepgrid() so it can preserve their full
level sets even when ordinary unused factor levels are dropped.
```

For the example, `.sepgrid_colnames()` returns:

```text
"sepgrid(member, time)"
```

### 2.4 Product-Margin Parser

Marker: `[NEW/MODIFIED]`

Local explanation:

```text
What this code does:
  Parses us(0 + member) and ar1(0 + time) into margin structure/variable pairs.

Why it exists:
  The public syntax is margin-oriented, but the backend needs explicit metadata:
  margin type, margin variable, and coordinate order.

Why this implementation:
  It accepts only structure(0 + one_variable) so the product syntax cannot be
  mistaken for random-slope support that the backend does not implement.
```

```r
.sep_is_zero <- function(x) {
    is.numeric(x) && length(x) == 1L && isTRUE(unname(x) == 0)
}

.sep_product_margin_var <- function(x) {
    if (is.call(x) && identical(.sep_call_name(x), "+") && length(x) == 3L) {
        if (.sep_is_zero(x[[2]]) && is.name(x[[3]])) return(.sep_deparse(x[[3]]))
        if (.sep_is_zero(x[[3]]) && is.name(x[[2]])) return(.sep_deparse(x[[2]]))
    }
    stop("separable() product margins currently require exactly one ",
         "no-intercept variable, e.g. us(0 + role) %x% ar1(0 + day).")
}

.sep_product_margin_spec <- function(x) {
    if (!is.call(x) || length(x) != 2L)
        stop("separable() product margins must look like us(0 + role) ",
             "or ar1(0 + day).")
    structure(.sep_product_margin_var(x[[2]]), names = .sep_deparse(x[[1]]))
}
```

For the example:

```text
us(0 + member)  -> named vector c(us = "member")
ar1(0 + time)   -> named vector c(ar1 = "time")
```

This is intentionally restrictive.  It rejects random-slope-like margins such
as:

```r
us(0 + member + x)
```

because those would require a separate marginal design compiler rather than the
current complete-grid indicator design.

### 2.5 Internal Grid Formula Builder

Marker: `[NEW/MODIFIED]`

Local explanation:

```text
What this code does:
  Builds sepgrid(member, time) + 0 | group as an unevaluated formula call.

Why it exists:
  mkReTrms() needs an ordinary random-effect design expression.

Why this implementation:
  Constructing a call keeps the rewrite inside R's normal formula machinery
  instead of creating strings that would need to be parsed later.
```

```r
.sepgrid_bar_call <- function(vars, group) {
    grid_call <- as.call(c(list(as.name("sepgrid")), lapply(vars, as.name)))
    as.call(list(as.name("|"),
                 as.call(list(as.name("+"), grid_call, 0)),
                 group))
}
```

For the example:

```text
vars  = c("member", "time")
group = group
```

returns:

```r
sepgrid(member, time) + 0 | group
```

This is the design expression `mkReTrms()` can understand.

### 2.6 Metadata Normalization And Registry

Marker: `[NEW/MODIFIED]`

Local explanation:

```text
What this code does:
  Normalizes margin metadata and declares what each supported margin contributes.

Why it exists:
  getReStruc() must count theta and C++ must know which density builder to use.

Why this implementation:
  A small registry keeps formula parsing separate from parameter-count and
  density-kind knowledge.
```

```r
.sep_spec_df <- function(x, what = "margin") {
    if (is.null(x)) return(NULL)
    if (is.matrix(x) || is.data.frame(x)) {
        ans <- as.data.frame(x, stringsAsFactors = FALSE)
    } else {
        ans <- data.frame(struc = unname(names(x)),
                          var = unname(x),
                          stringsAsFactors = FALSE)
    }
    if (!all(ans$struc %in% names(.sep_margin_registry))) {
        bad <- unique(ans$struc[!ans$struc %in% names(.sep_margin_registry)])
        stop("Unsupported separable() ", what, ": ", paste(bad, collapse = ", "))
    }
    ans
}

.sep_parse_spec <- function(x) {
    if (is.list(x) && !is.null(x$grid) && !is.null(x$margins)) {
        x$margins <- as.matrix(x$margins)
        return(x)
    }
    stop("Internal separable() metadata is missing or malformed.")
}

.sep_margin_registry <- list(
    cs = list(
        code = "cs",
        density_kind = "dense_corr",
        can_scale = TRUE,
        n_scale = function(n) as.integer(n),
        n_corr = function(n) 1L
    ),
    homcs = list(
        code = "homcs",
        density_kind = "dense_corr",
        can_scale = TRUE,
        n_scale = function(n) 1L,
        n_corr = function(n) 1L
    ),
    us = list(
        code = "us",
        density_kind = "dense_corr",
        can_scale = TRUE,
        n_scale = function(n) as.integer(n),
        n_corr = function(n) as.integer(n * (n - 1L) / 2L)
    ),
    ar1 = list(
        code = "ar1",
        density_kind = "ar1",
        can_scale = FALSE,
        n_scale = function(n) 0L,
        n_corr = function(n) 1L
    )
)
```

For `us x ar1`:

```text
us:
  dense_corr margin
  can carry scale
  n scale parameters = n_member log-SDs
  n correlation parameters = n_member * (n_member - 1) / 2

ar1:
  AR(1) margin
  cannot carry scale inside separable()
  n scale parameters = 0
  n correlation parameters = 1
```

For `cs x ar1`, the dense margin is like `homcs` in that it has one
compound-symmetry correlation parameter, but like `us` in that it has one
scale parameter per member level.

### 2.7 R/C++ Integer Codes

Marker: `[NEW/MODIFIED]`

Local explanation:

```text
What this code does:
  Defines compact integer codes that R passes to C++.

Why it exists:
  TMB data structures use numeric codes more naturally than R strings.

Why this implementation:
  The codes separate exact covariance names (us/ar1) from broader evaluator
  kinds (dense correlation / AR1).
```

```r
.sep_density_kind_code <- c(dense_corr = 1L, ar1 = 2L)

.sep_dispatch_code <- c(dense_ar1 = 1L)

.sep_scale_mode_code <- c(
    margin = 1L
)
```

These must match the C++ enums.

For the example:

```text
sepDensityKinds = c(1, 2)
sepDispatch     = 1
sepScaleMode    = 1
sepScaleSpec    = 0
```

### 2.8 Supported Pair Lookup

Marker: `[NEW/MODIFIED]`

Local explanation:

```text
What this code does:
  Checks whether the two margins match an implemented separable evaluator.

Why it exists:
  The prototype supports dense_corr x ar1, not every possible product.

Why this implementation:
  The lookup is order-insensitive for dispatch but leaves coordinate order in
  sepCodes/sepDensityKinds for C++ reshaping.
```

```r
.sep_supported_pairs <- list(
    list(
        dispatch = "dense_ar1",
        density_kinds = c("dense_corr", "ar1"),
        allowed_codes = list(dense_corr = c("cs", "homcs", "us"))
    )
)

.sep_pair_dispatch <- function(regs) {
    kinds <- vapply(regs, `[[`, character(1), "density_kind")
    codes <- vapply(regs, `[[`, character(1), "code")
    for (pair in .sep_supported_pairs) {
        if (!identical(unname(sort(kinds)),
                       unname(sort(pair$density_kinds)))) {
            next
        }
        allowed <- pair$allowed_codes
        allowed_ok <- TRUE
        for (kind in names(allowed)) {
            codes_for_kind <- codes[kinds == kind]
            allowed_ok <- allowed_ok &&
                length(codes_for_kind) > 0L &&
                all(codes_for_kind %in% allowed[[kind]])
        }
        if (allowed_ok) return(pair$dispatch)
    }
    NA_character_
}
```

For the example:

```text
kinds = c("dense_corr", "ar1")
codes = c("us", "ar1")
dispatch = "dense_ar1"
```

The lookup is order-insensitive, so `ar1(...) %x% us(...)` would select the same
dispatch while preserving coordinate order in other metadata.

### 2.9 Scale Resolution

Marker: `[NEW/MODIFIED]`

Local explanation:

```text
What this code does:
  Decides which margin supplies absolute standard deviations.

Why it exists:
  A Kronecker product with arbitrary scale on both margins is not identifiable.

Why this implementation:
  The current member-time use case has scale on the dense member margin and
  correlation-only AR(1) over time.
```

```r
.sep_margin_label <- function(x) {
    x <- as.data.frame(x, stringsAsFactors = FALSE)
    paste0(x$struc, "(", x$var, ")", collapse = " x ")
}

.sep_scale_info <- function(margins, regs, scale = NULL) {
    can_scale <- vapply(regs, `[[`, logical(1), "can_scale")
    scale_candidates <- which(can_scale)

    if (is.null(scale)) {
        scale_margin <- scale_candidates
    } else {
        scale_spec <- .sep_spec_df(scale, "scale margin")
        if (nrow(scale_spec) != 1L) {
            stop("separable() scale must be a single margin call such as ",
                 "scale = us(0 + member).")
        }
        scale_margin <- which(margins$struc == scale_spec$struc &
                              margins$var == scale_spec$var)
        if (length(scale_margin) != 1L) {
            stop("separable() scale must match one of the specified margins, ",
                 "for example scale = us(0 + member) when us(0 + member) ",
                 "is a margin.")
        }
        if (!regs[[scale_margin]]$can_scale) {
            stop("separable() scale = ", scale_spec$struc, "(",
                 scale_spec$var, ") selects a correlation-only margin. ",
                 "Use a scale-capable margin such as cs(), homcs(), or us().")
        }
    }

    if (length(scale_margin) == 0L) {
        stop("separable() margins ", .sep_margin_label(margins),
             " define only a correlation product. This needs an explicit ",
             "overall scale, but scale = global() is not implemented yet.")
    }
    if (length(scale_margin) > 1L) {
        stop("separable() requires exactly one scale-carrying margin. ",
             "Specify it explicitly with scale = us(0 + member) or use one ",
             "scale-capable margin with one correlation-only margin.")
    }

    list(
        mode = "margin",
        mode_code = as.integer(.sep_scale_mode_code[["margin"]]),
        spec = as.integer(scale_margin - 1L),
        margin = scale_margin
    )
}
```

For `us x ar1` with no explicit `scale`:

```text
can_scale        = c(TRUE, FALSE)
scale_margin     = 1
sepScaleMode     = 1
sepScaleSpec     = 0
```

`us(member)` supplies absolute SDs.  `ar1(time)` supplies only correlation.

### 2.10 Unsupported Pair Errors

Marker: `[NEW/MODIFIED]`

Local explanation:

```text
What this code does:
  Turns unsupported separable combinations into targeted R errors.

Why it exists:
  Users should fail before C++ when they request ar1 x ar1, dense x dense, or
  invalid scale ownership.

Why this implementation:
  The error path distinguishes scale-identifiability problems from missing C++
  evaluator support.
```

```r
.sep_stop_unsupported_pair <- function(margins, regs, scale = NULL) {
    can_scale <- vapply(regs, `[[`, logical(1), "can_scale")

    if (sum(can_scale) == 0L) {
        stop("separable() margins ", .sep_margin_label(margins),
             " currently define only a correlation product. A future global ",
             "scale mode is planned, for example scale = global(), but it is ",
             "not implemented yet.")
    }

    if (is.null(scale) && sum(can_scale) > 1L) {
        stop("More than one separable() margin can carry scale in ",
             .sep_margin_label(margins), ". Please specify the scale margin ",
             "explicitly, for example scale = ", margins$struc[which(can_scale)[1]],
             "(", margins$var[which(can_scale)[1]], "). This covariance pair ",
             "is also outside the current dense x ar1 prototype.")
    }

    if (!is.null(scale)) {
        scale <- .sep_spec_df(scale, "scale margin")
        scale_match <- which(margins$struc == scale$struc &
                             margins$var == scale$var)
        if (length(scale_match) == 1L && !regs[[scale_match]]$can_scale) {
            stop("separable() scale = ", scale$struc, "(",
                 scale$var, ") selects a correlation-only margin. ",
                 "Use a scale-capable margin such as cs(), homcs(), or us(), ",
                 "or wait for a future global scale mode.")
        }
        if (length(scale_match) == 1L) {
            stop("separable() does not yet implement the covariance pair ",
                 .sep_margin_label(margins), ". The explicit scale selector ",
                 "chooses the scale margin, but it does not enable unsupported ",
                 "density combinations.")
        }
    }

    stop("separable() currently supports cs(0 + member), homcs(0 + member), ",
         "or us(0 + member) crossed with ar1(0 + time).")
}
```

This is not used for valid `us x ar1`, but it is part of the parser path: if the
supported-pair lookup fails, errors are reported before C++ sees the term.

### 2.11 `getReStruc()` Metadata Builder

Marker: `[NEW/MODIFIED]`

Local explanation:

```text
What this code does:
  Converts rewritten formula metadata plus mkReTrms column names into the final
  R-side separable term metadata.

Why it exists:
  C++ should receive dimensions, codes, theta count, and scale info without
  needing to understand R formulas.

Why this implementation:
  getReStruc() is the existing boundary where glmmTMB converts random-effect
  parser output into TMB term metadata.
```

```r
.sep_restruc_info <- function(spec, cnms, blksize) {
    ## Centralize the R-side contract for currently supported separable terms.
    ## This keeps `getReStruc()` from knowing how many parameters each margin
    ## contributes, and gives future margins one obvious place to declare their
    ## parameter counts and density kind.
    spec <- .sep_parse_spec(spec)

    ## `cnms` are the model-matrix column names created by
    ## `sepgrid(member, time) + 0`.  Their numeric labels define the latent cell
    ## grid that C++ will later reshape.  For two members and three time points:
    ##
    ##   cnms    = "(1,1)", "(2,1)", "(1,2)", "(2,2)", "(1,3)", "(2,3)"
    ##   dims    = c(2, 3)
    ##   blksize = 6
    ##
    ## `prod(dims) == blksize` is what proves that the random-effect block is a
    ## complete rectangle rather than a ragged collection of observed cells.
    coords <- parseNumLevels(cnms)
    if (ncol(coords) != 2L)
        stop("separable() currently requires a two-dimensional sepgrid().")
    dims <- as.integer(apply(coords, 2, function(z) length(unique(z))))
    if (prod(dims) != blksize)
        stop("separable() requires a complete rectangular sepgrid().")

    margins <- .sep_spec_df(spec$margins)
    if (nrow(margins) != 2L)
        stop("separable() currently requires exactly two margin structures.")

    pair <- margins$struc
    regs <- .sep_margin_registry[pair]
    dispatch <- .sep_pair_dispatch(regs)
    if (is.na(dispatch)) {
        .sep_stop_unsupported_pair(margins, regs, spec$scale)
    }
    if (!identical(margins$var, spec$grid)) {
        stop("The separable() margin variables must match sepgrid() variables. ",
             "Use, for example, ",
             "separable(homcs(0 + member) %x% ar1(0 + time) | group).")
    }

    scale_info <- .sep_scale_info(margins, regs, spec$scale)

    ## Count theta in margin order, matching both the user's product syntax and
    ## the C++ parser.  Only the scale-carrying margin contributes log-SDs; every
    ## supported margin contributes its own correlation parameters.
    ##
    ## Example:
    ##   separable(us(0 + member) %x% ar1(0 + time) | group)
    ##   theta = member_logsd..., member_corr..., time_ar1_phi
    ntheta <- sum(vapply(seq_along(regs), function(i) {
        nscale <- if (i == scale_info$margin) regs[[i]]$n_scale(dims[[i]]) else 0L
        nscale + regs[[i]]$n_corr(dims[[i]])
    }, integer(1)))
    density_kind <- vapply(regs, `[[`, character(1), "density_kind")

    list(
        dims = dims,
        codes = as.integer(vapply(pair, function(z) .valid_covstruct[[z]], numeric(1))),
        density_kinds = as.integer(.sep_density_kind_code[density_kind]),
        dispatch = as.integer(.sep_dispatch_code[dispatch]),
        scale_mode = scale_info$mode_code,
        scale_spec = scale_info$spec,
        ntheta = as.integer(ntheta),
        density_kind = density_kind,
        margins = margins,
        scale = spec$scale
    )
}
```

For `n_member = M`, `n_time = T`, `us x ar1`:

```text
dims           = c(M, T)
codes          = c(us = 1, ar1 = 3)
density_kinds  = c(1, 2)
dispatch       = 1
scale_mode     = 1
scale_spec     = 0
ntheta         = M + M*(M - 1)/2 + 1
```

`cnms` are the random-effect column names created by existing `mkReTrms()`.
They are parseable because the rewritten design uses `sepgrid()`.

Guide-only reading notes:

| code symbol | toy value | what to check mentally |
|---|---|---|
| `spec` | list with `grid`, `margins`, optional `scale` | came from the rewritten formula metadata, not from model data |
| `cnms` | `"(1,1)", "(2,1)", ...` | defines flat cell order |
| `coords` | 6-row numeric matrix | parsed version of `cnms` |
| `dims` | `c(2, 3)` | first dimension has 2 member levels, second has 3 time levels |
| `blksize` | `6` | number of random-effect cells inside one group |
| `regs` | registry entries for `us` and `ar1` | tells R how many parameters each margin contributes |
| `dispatch` | `dense_ar1` | tells C++ which evaluator branch to use |
| `scale_info$margin` | `1` in R indexing | member margin carries SDs |
| `scale_info$spec` | `0` in C++ indexing | same margin, converted to zero-based indexing |

The most important invariant in this function is:

```text
prod(dims) == blksize
```

For the toy example:

```text
prod(c(2, 3)) = 6
```

If that check failed, the code would not know how to reshape one group's flat
random-effect vector into a rectangular `z[member, time]` array.

### 2.12 Product Spec Compiler

Marker: `[NEW/MODIFIED]`

Local explanation:

```text
What this code does:
  Compiles the public separable(us(...) %x% ar1(...) | group) call into a grid
  expression plus structured metadata.

Why it exists:
  The public syntax is clear for users, but the backend needs a complete-grid
  design and explicit margin metadata.

Why this implementation:
  It keeps all product-syntax restrictions and validation in one place.
```

```r
.sep_make_product_spec <- function(bar_expr, scale = NULL) {
    ## Compile the public product syntax
    ##
    ##   separable(us(0 + role) %x% ar1(0 + day) | group, scale = us(0 + role))
    ##
    ## into the same complete-grid representation used by the current backend.
    ## Only simple no-intercept one-variable margins are accepted for now; this
    ## keeps the door open for a later marginal-design compiler without implying
    ## that random-slope margins are already supported.
    if (!is.call(bar_expr) || !identical(.sep_call_name(bar_expr), "|") ||
        length(bar_expr) != 3L) {
        stop("separable() product syntax must look like ",
             "separable(us(0 + role) %x% ar1(0 + day) | group, ...).")
    }
    prod_expr <- bar_expr[[2]]
    if (!is.call(prod_expr) || !identical(.sep_call_name(prod_expr), "%x%") ||
        length(prod_expr) != 3L) {
        stop("separable() product syntax requires exactly two margins joined ",
             "by %x%, e.g. us(0 + role) %x% ar1(0 + day).")
    }

    margin_calls <- list(prod_expr[[2]], prod_expr[[3]])
    margins <- c(.sep_product_margin_spec(margin_calls[[1]]),
                 .sep_product_margin_spec(margin_calls[[2]]))
    if (anyDuplicated(unname(margins))) {
        stop("separable() product margins must use distinct variables.")
    }
    if (!all(names(margins) %in% names(.sep_margin_registry))) {
        bad <- unique(names(margins)[!names(margins) %in% names(.sep_margin_registry)])
        stop("Unsupported separable() margin: ", paste(bad, collapse = ", "))
    }

    scale_spec <- NULL
    if (!is.null(scale)) {
        scale_vec <- .sep_product_margin_spec(scale)
        scale_spec <- .sep_spec_df(scale_vec, "scale margin")
    }

    m <- .sep_spec_df(margins)

    ## `grid_expr` is the formula fragment that existing random-effect machinery
    ## can already process:
    ##
    ##   sepgrid(member, time) + 0 | group
    ##
    ## The separate `margins`/`scale` metadata is carried beside that fragment so
    ## the columns still get the separable covariance rather than an ordinary
    ## unstructured covariance on the flattened grid.
    list(grid = unname(margins),
         margins = m,
         scale = scale_spec,
         grid_expr = .sepgrid_bar_call(unname(margins), bar_expr[[3]]))
}
```

For the example:

```text
bar_expr = us(0 + member) %x% ar1(0 + time) | group
margins  = c(us = "member", ar1 = "time")
grid     = c("member", "time")
scale    = NULL
grid_expr = sepgrid(member, time) + 0 | group
```

### 2.13 Metadata Call Builder

Marker: `[NEW/MODIFIED]`

Local explanation:

```text
What this code does:
  Turns structured metadata into an unevaluated list(...) call.

Why it exists:
  reformulas::splitForm() already has a term-local reTrmAddArgs channel.

Why this implementation:
  It is less fragile than a positional sidecar because metadata stays attached
  to the exact rewritten separable term.
```

```r
.sep_spec_call <- function(spec) {
    chr_vec_call <- function(x) as.call(c(list(as.name("c")), as.list(x)))
    margins <- as.data.frame(spec$margins, stringsAsFactors = FALSE)
    call_args <- list(
        as.name("list"),
        grid = chr_vec_call(spec$grid),
        margins = as.call(list(
            as.name("data.frame"),
            struc = chr_vec_call(margins$struc),
            var = chr_vec_call(margins$var),
            stringsAsFactors = FALSE
        ))
    )
    if (!is.null(spec$scale)) {
        scale <- as.data.frame(spec$scale, stringsAsFactors = FALSE)
        call_args$scale <- as.call(list(
            as.name("data.frame"),
            struc = chr_vec_call(scale$struc),
            var = chr_vec_call(scale$var),
            stringsAsFactors = FALSE
        ))
    }
    as.call(call_args)
}
```

This turns structured metadata into an unevaluated `list(...)` call.  The reason
is pragmatic: `reformulas::splitForm()` already keeps additional arguments
attached to the correct random-effect term in `reTrmAddArgs`.

This is ugly but robust for:

```r
separable(...) + separable(...) + separable(...)
```

because each rewritten term carries its own metadata.

### 2.14 Formula Tree Rewrite

Marker: `[NEW/MODIFIED]`

Local explanation:

```text
What this code does:
  Walks the formula tree and replaces public separable product calls with the
  internal grid-plus-metadata form.

Why it exists:
  The rewrite must happen before model.frame(), splitForm(), and mkReTrms().

Why this implementation:
  It changes only separable(...) calls and recurses through larger formulas, so
  ordinary fixed/random effects keep using the existing path.
```

```r
.rewrite_separable_expr <- function(x) {
    if (!is.call(x)) return(x)
    if (identical(.sep_call_name(x), "separable")) {
        args <- as.list(x[-1])
        nms <- names(args)
        if (is.null(nms)) nms <- rep("", length(args))
        scale_i <- which(nms == "scale")
        if (length(scale_i) > 1L)
            stop("separable() accepts at most one scale argument.")
        scale_arg <- if (length(scale_i)) args[[scale_i]] else NULL
        named_i <- which(nzchar(nms) & nms != "scale")
        if (length(named_i) > 0L) {
            stop("separable() only accepts named argument scale.")
        }

        unnamed_args <- args[!nzchar(nms)]
        if (length(unnamed_args) == 1L) {
            spec <- .sep_make_product_spec(unnamed_args[[1]], scale = scale_arg)
            return(as.call(list(as.name("separable"),
                                spec$grid_expr,
                                .sep_spec_call(spec))))
        }

        if (length(unnamed_args) == 2L && is.call(unnamed_args[[2]]) &&
            identical(.sep_call_name(unnamed_args[[2]]), "list") &&
            is.null(scale_arg)) {
            return(x)
        }

        stop("separable() requires separable(",
             "margin1(0 + variable) %x% margin2(0 + variable) | group).")
    }
    for (i in seq_along(x)[-1]) x[[i]] <- .rewrite_separable_expr(x[[i]])
    x
}

rewrite_separable_formula <- function(f) {
    if (!inherits(f, "formula")) return(f)
    f[[length(f)]] <- .rewrite_separable_expr(f[[length(f)]])
    f
}
```

This is where the public product syntax is rewritten.  It also preserves the
internal `separable(grid, list(...))` form if one is already present.

## 3. Model Frame: Preserve Complete `sepgrid()` Levels

Location:

```text
glmmTMB/R/glmmTMB.R
```

### Existing Context

`glmmTMB()` constructs a combined formula for `model.frame()` from conditional,
zero-inflation, and dispersion formulas.  Existing code then evaluates the
model frame and may drop unused factor levels.

### New/Modified Code

Local explanation:

```text
What this code does:
  Evaluates the model frame without dropping sepgrid levels, then drops ordinary
  unused factor levels by hand.

Why it exists:
  Dropping sepgrid levels would remove intended latent cells and change AR(1)
  distances.

Why this implementation:
  It limits special handling to generated sepgrid columns and leaves ordinary
  factor behavior as close as possible to existing glmmTMB behavior.
```

```r
formList <- list(formula, ziformula, dispformula)
sepgrid_cols <- .sepgrid_colnames(formula, ziformula, dispformula)
for (i in seq_along(formList)) {
    f <- formList[[i]]
    f <- noSpecials(sub_specials(f), delete=FALSE, specials = c(names(.valid_covstruct), "s"))
    formList[[i]] <- f
}
combForm <- do.call(addForm,formList)
```

Then:

```r
mf$formula <- combForm
if (length(sepgrid_cols) > 0L && control$drop_unused_levels) {
    mf$drop.unused.levels <- FALSE
}
fr <- eval(mf,envir=environment(formula),enclos=parent.frame())
if (length(sepgrid_cols) > 0L && control$drop_unused_levels) {
    for (nm in names(fr)) {
        if (is.factor(fr[[nm]]) && !(nm %in% sepgrid_cols)) {
            fr[[nm]] <- droplevels(fr[[nm]])
        }
    }
}
```

### What Happens

For the example, `sepgrid(member, time)` becomes a model-frame column named:

```text
sepgrid(member, time)
```

If `time` is a factor with levels `1:30` but day 17 has no observed rows,
`sepgrid()` still has the day-17 cells in its levels.  The code above prevents
`model.frame()` from dropping those levels.  It then restores ordinary
level-dropping behavior for non-`sepgrid()` factors.

## 4. Existing Random-Effect Machinery Builds The Design

Location:

```text
glmmTMB/R/glmmTMB.R, getXReTrms()
```

Marker: `[EXISTING]`

Local explanation:

```text
What this code does:
  Builds the sparse random-effect design objects from the rewritten formula.

Why it matters here:
  This is where sepgrid(member, time) becomes actual random-effect columns.

Why this existing path is reused:
  Reusing mkReTrms avoids writing a new random-effect model-matrix builder for
  the prototype.
```

```r
RHSForm(ranform) <- subbars(RHSForm(reOnly(formula)))

if (has_re) {
    mf$formula <- ranform
    reTrms <- mkReTrms(no_specials(
        findbars_x(formula)),
        fr, reorder.terms=FALSE, calc.lambdat=FALSE, sparse = TRUE)
} else {
    reTrms <- list(Ztlist = list(), flist = list(), cnms = list(),
                   theta = list())
}
```

Explanation:

```text
findbars_x(formula)
  Finds random-effect terms after the rewrite.

no_specials(...)
  Removes the covariance-structure wrapper so mkReTrms sees an ordinary
  random-effect term.

mkReTrms(...)
  Builds Z, Ztlist, flist, cnms, Gp, and block sizes.
```

For the example, `mkReTrms()` sees the rewritten design:

```r
sepgrid(member, time) + 0 | group
```

and therefore builds one random-effect column per complete member-time cell.

## 5. Existing `splitForm()` Carries The Metadata

Location:

```text
glmmTMB/R/glmmTMB.R, getXReTrms()
reformulas::splitForm()
```

Marker: `[EXISTING]`

Local explanation:

```text
What this code does:
  Splits the rewritten formula into random-effect formulas, covariance classes,
  and extra arguments.

Why it matters here:
  The extra argument is the generated metadata list for this exact separable
  term.

Why this existing path is reused:
  reTrmAddArgs already preserves term-local alignment across multiple random
  effects.
```

```r
ss <- splitForm(formula, specials = c(names(.valid_covstruct), "s"))
```

Existing `reformulas::splitForm()` returns:

```text
fixedFormula
reTrmFormulas
reTrmAddArgs
reTrmClasses
```

For the rewritten example:

```text
reTrmClasses[[i]]  = "separable"
reTrmFormulas[[i]] = sepgrid(member, time) + 0 | group
reTrmAddArgs[[i]]  = separable(list(...metadata...)) shape internally
```

More precisely, glmmTMB later evaluates the second element of `reTrmAddArgs`.

Marker: `[EXISTING]`

```r
get_arg <- function(v) {
  if (length(v) == 1) return(NA_real_)
  payload <- v[[2]]
  res <- tryCatch(eval(payload, envir = fr,
                       enclos = environment(formula)),
                  error = function(e)
                    stop("can't evaluate argument ",
                         sQuote(deparse(payload)),
                         call. = FALSE))
  return(res)
}
aa <- lapply(ss$reTrmAddArgs, get_arg)
```

For the example, `aa[[i]]` becomes an R list:

```r
list(
  grid = c("member", "time"),
  margins = data.frame(struc = c("us", "ar1"),
                       var = c("member", "time"),
                       stringsAsFactors = FALSE)
)
```

Important naming point:

```text
aa is not newly invented for separable().

In existing glmmTMB code, the same `aa` pathway is also used for extra
covariance-structure arguments such as reduced-rank rank.  For separable terms,
`aa[[i]]` is instead the evaluated metadata payload:

  grid
  margins
  optional scale selector

That is why getReStruc() has to branch before ordinary theta counting.  It must
turn this list into compact numeric fields before C++ sees the term.
```

For this toy model, the handoff is:

| object | value | role |
|---|---|---|
| `ss[[i]]` | `"separable"` | says this random-effect term uses the new covariance code |
| `aa[[i]]$grid` | `c("member", "time")` | says which variables built the `sepgrid()` factor |
| `aa[[i]]$margins$struc` | `c("us", "ar1")` | says which covariance structure belongs to each coordinate |
| `aa[[i]]$margins$var` | `c("member", "time")` | cross-checks that margins match the grid variables |
| `reTrms$cnms[[i]]` | `"(1,1)", "(2,1)", ...` | gives the realized cell order and dimensions |

## 6. Existing `reXterms` Construction

Location:

```text
glmmTMB/R/glmmTMB.R, getXReTrms()
```

Marker: `[EXISTING]`

Local explanation:

```text
What this code does:
  Builds terms objects for the random-effect left-hand-side expressions.

Why it matters here:
  The separable term's expression is already sepgrid(member, time) + 0, so the
  existing terms path can evaluate it.

Why no new separable code is needed here:
  The covariance-specific logic is carried in metadata; the design expression is
  ordinary after rewriting.
```

```r
termsfun <- function(x) {
    ff <- eval(substitute( ~ foo, list(foo = x[[2]])))
    tt <- try(terms(ff, data=fr), silent=TRUE)
    if (inherits(tt,"try-error")) {
        stop(
            sprintf("can't evaluate RE term %s: simplify?",
                    sQuote(deparse(ff)))
        )
    }
    tt
}

drop_s <- function(f, a) {
    if (identical(a[[1]], as.symbol('s'))) NA else termsfun(f)
}
reXterms <- Map(drop_s, ss$reTrmFormulas, ss$reTrmAddArgs)
```

For the example, `reXterms[[i]]` is based on:

```r
~ sepgrid(member, time) + 0
```

This is pre-existing glmmTMB machinery.  The separable branch mostly needs the
`cnms` and block sizes from `reTrms`, not special handling here.

## 7. `mkTMBStruc()` Calls `getReStruc()`

Location:

```text
glmmTMB/R/glmmTMB.R, mkTMBStruc()
```

Marker: `[EXISTING]`

Local explanation:

```text
What this code does:
  Runs getXReTrms() for cond/zi/disp components and then asks getReStruc() to
  summarize each random-effect term for TMB.

Why it matters here:
  This is the handoff from formula/design objects to covariance-structure
  metadata.

Why this existing path is reused:
  Separable terms should behave like other glmmTMB covariance structures at the
  component level.
```

```r
condList  <- getXReTrms(formula, mf, fr, type="conditional", contrasts=contrasts, sparse=sparseX[["cond"]],
                        old_smooths = old_smooths$cond)
ziList    <- getXReTrms(ziformula, mf, fr, type="zero-inflation", contrasts=contrasts, sparse=sparseX[["zi"]],
                        old_smooths = old_smooths$zi)
dispList  <- getXReTrms(dispformula, mf, fr, type="dispersion", contrasts=contrasts, sparse=sparseX[["disp"]],
                        old_smooths = old_smooths$disp)
```

Then:

```r
condReStruc <- with(condList, getReStruc(reTrms, ss, aa, reXterms, fr, fc_list[[1]]))
ziReStruc <- with(ziList, getReStruc(reTrms, ss, aa, reXterms, fr, fc_list[[2]]))
dispReStruc <- with(dispList, getReStruc(reTrms, ss, aa, reXterms, fr, fc_list[[3]]))
```

For a conditional `us x ar1` term, the relevant call is the first one.

## 8. `getReStruc()` Adds Separable Metadata

Location:

```text
glmmTMB/R/glmmTMB.R
```

### Existing Context

`getReStruc()` already computes:

```text
blockReps
blockSize
blockNumTheta
blockCode
simCode
fullCor
```

for each random-effect term.

### New/Modified Code: Precompute `sepInfo`

Local explanation:

```text
What this code does:
  Computes separable-specific metadata before general parameter counting.

Why it exists:
  The separable theta count depends on margin types and grid dimensions, not
  just block size.

Why this implementation:
  It lets the later generic getReStruc loop stay mostly unchanged.
```

```r
## For ordinary covariance structures, `aa` is unused except by `rr`.
## For separable terms, `aa[[i]]` is the metadata list that was injected by
## `rewrite_separable_formula()` as the second argument to the internal
## `separable(sepgrid(...) | group, list(...))` call.  Convert it here into
## compact numeric fields before the main per-term metadata list is built.
sepInfo <- vector("list", length(ss))
for (i in which(ss == "separable")) {
    sepInfo[[i]] <- .sep_restruc_info(aa[[i]], reTrms$cnms[[i]], blksize[i])
}
```

For the example:

```text
aa[[i]]
  metadata list from the rewritten formula

reTrms$cnms[[i]]
  sepgrid column names, e.g. "(1,1)", "(2,1)", "(1,2)", ...

blksize[i]
  n_member * n_time
```

### New/Modified Code: Parameter Count

Local explanation:

```text
What this code does:
  Uses the separable metadata's ntheta for separable terms.

Why it exists:
  A full us covariance over M*T cells would have far too many parameters; the
  separable us x ar1 model has only member SD/correlation plus AR(1) phi.

Why this implementation:
  The existing parFun switch already centralizes theta counts for covariance
  structures, so separable gets one special branch there.
```

```r
parFun <- function(struc, blksize, blkrank, sep_info) {
    if (struc == "separable") {
        ## Separable theta counts are not a function of the flattened block
        ## size alone.  For example, `us(member) %x% ar1(time)` with
        ## 2 members and 3 times has 4 theta values, not the 21 that an
        ## unstructured 6-cell covariance would need.
        return(sep_info$ntheta)
    }
    switch(as.character(struc),
           "diag" = blksize,
           "us" = blksize * (blksize+1) / 2,
           "cs" = blksize + 1,
           "ar1" = 2,
           "hetar1" = blksize + 1,
           "ou" = 2,
           "exp" = 2,
           "gau" = 2,
           "mat" = 3,
           "toep" = 2 * blksize - 1,
           "rr" = blksize * blkrank - (blkrank - 1) * blkrank / 2,
           "homdiag" = 1,
           "propto" = blksize * (blksize+1) / 2 + 1,
           "homcs" = 2,
           "homtoep" = blksize,
           "equalto" = blksize * (blksize+1) / 2,
           stop(sprintf("undefined number of parameters for covstruct '%s'", struc))
           )
}
blockNumTheta <- mapply(parFun, ss, blksize, blkrank, sepInfo,
                        SIMPLIFY=FALSE)
```

For `us x ar1`:

```text
blockNumTheta = n_member + n_member*(n_member - 1)/2 + 1
```

The `+1` is AR(1) phi's unconstrained parameter.

### New/Modified Code: Attach Metadata For C++

Local explanation:

```text
What this code does:
  Adds separable dimensions, margin codes, dispatch, and scale metadata to the
  term list passed to C++.

Why it exists:
  C++ needs numeric metadata to reshape U, parse theta, select the evaluator,
  and apply scale.

Why this implementation:
  It extends the existing per-term list rather than introducing a separate C++
  data structure for separable terms.
```

```r
} else if(ss[i] == "separable") {
    tmp$sepDims <- sepInfo[[i]]$dims
    tmp$sepCodes <- sepInfo[[i]]$codes
    tmp$sepDensityKinds <- sepInfo[[i]]$density_kinds
    tmp$sepDispatch <- sepInfo[[i]]$dispatch
    tmp$sepScaleMode <- sepInfo[[i]]$scale_mode
    tmp$sepScaleSpec <- sepInfo[[i]]$scale_spec
}
```

For the example:

```text
tmp$blockCode       = separable code
tmp$blockSize       = M * T
tmp$blockReps       = number of group levels
tmp$blockNumTheta   = M + M*(M - 1)/2 + 1
tmp$sepDims         = c(M, T)
tmp$sepCodes        = c(us_covstruct, ar1_covstruct)
tmp$sepDensityKinds = c(dense_corr_sep, ar1_sep)
tmp$sepDispatch     = dense_ar1_dispatch
tmp$sepScaleMode    = margin_sep_scale
tmp$sepScaleSpec    = 0
```

This `tmp` list becomes part of the `terms` data structure passed to C++.

## 9. C++ Enum Codes

Location:

```text
glmmTMB/src/glmmTMB.cpp
```

Marker: `[NEW/MODIFIED]`

Local explanation:

```text
What this code does:
  Adds numeric codes for separable covariance structures and separable dispatch.

Why it exists:
  R passes compact numeric codes to the TMB template.

Why this implementation:
  It keeps exact covariance codes, broad density-kind codes, dispatch codes, and
  scale-mode codes separate so future combinations can be added without
  overloading one integer.
```

```cpp
enum valid_covStruct {
  diag_covstruct = 0,
  us_covstruct   = 1,
  cs_covstruct   = 2,
  ar1_covstruct  = 3,
  ou_covstruct   = 4,
  exp_covstruct = 5,
  gau_covstruct = 6,
  mat_covstruct = 7,
  toep_covstruct = 8,
  rr_covstruct = 9,
  homdiag_covstruct = 10,
  propto_covstruct = 11,
  hetar1_covstruct = 12,
  homcs_covstruct = 13,
  homtoep_covstruct = 14,
  equalto_covstruct = 15,
  separable_covstruct = 16
};

enum separable_density_kind {
  dense_corr_sep = 1,
  ar1_sep = 2
};

enum separable_dispatch {
  dense_ar1_dispatch = 1
};

enum separable_scale_mode {
  margin_sep_scale = 1
};
```

The `us x ar1` term reaches C++ as:

```text
blockCode = separable_covstruct
sepCodes = c(us_covstruct, ar1_covstruct)
sepDensityKinds = c(dense_corr_sep, ar1_sep)
sepDispatch = dense_ar1_dispatch
sepScaleMode = margin_sep_scale
```

## 10. C++ Per-Term Metadata

Location:

```text
glmmTMB/src/glmmTMB.cpp
```

### Existing Context

Each random-effect term is represented in C++ by `per_term_info`.

### New/Modified Code

Local explanation:

```text
What this code does:
  Adds optional sep* metadata fields to each C++ random-effect term.

Why it exists:
  The likelihood evaluator needs separable dimensions and margin metadata in
  addition to the usual block size/repetition information.

Why this implementation:
  Optional fields let non-separable covariance structures keep using the same
  per_term_info object unchanged.
```

```cpp
template <class Type>
struct per_term_info {
  int blockCode;
  int blockSize;
  int blockReps;
  int blockNumTheta;
  int simCode;
  int fullCor;
  matrix<Type> dist;
  vector<Type> times;
  vector<int> sepDims;
  vector<int> sepCodes;
  vector<int> sepDensityKinds;
  vector<int> sepDispatch;
  vector<int> sepScaleMode;
  vector<int> sepScaleSpec;
  matrix<Type> corr;
  vector<Type> sd;
  matrix<Type> fact_load;
};
```

The fields before `sepDims` existed already.  The `sep*` vectors are new.  They
are optional and empty for non-separable covariance structures.

### New/Modified Code: Read Metadata From R

Local explanation:

```text
What this code does:
  Copies sep* fields from the R term list into C++ vectors.

Why it exists:
  getReStruc() prepared these fields as TMB data; the template must read them
  before termwise_nll() can dispatch.

Why this implementation:
  Each field is optional so ordinary covariance structures do not need dummy
  separable metadata.
```

```cpp
SEXP sdims = getListElement(y, "sepDims");
if(!Rf_isNull(sdims)){
  RObjectTestExpectedType(sdims, &Rf_isNumeric, "sepDims");
  (*this)(i).sepDims = asVector<int>(sdims);
}
SEXP scodes = getListElement(y, "sepCodes");
if(!Rf_isNull(scodes)){
  RObjectTestExpectedType(scodes, &Rf_isNumeric, "sepCodes");
  (*this)(i).sepCodes = asVector<int>(scodes);
}
SEXP skinds = getListElement(y, "sepDensityKinds");
if(!Rf_isNull(skinds)){
  RObjectTestExpectedType(skinds, &Rf_isNumeric, "sepDensityKinds");
  (*this)(i).sepDensityKinds = asVector<int>(skinds);
}
SEXP sdispatch = getListElement(y, "sepDispatch");
if(!Rf_isNull(sdispatch)){
  RObjectTestExpectedType(sdispatch, &Rf_isNumeric, "sepDispatch");
  (*this)(i).sepDispatch = asVector<int>(sdispatch);
}
SEXP sscalemode = getListElement(y, "sepScaleMode");
if(!Rf_isNull(sscalemode)){
  RObjectTestExpectedType(sscalemode, &Rf_isNumeric, "sepScaleMode");
  (*this)(i).sepScaleMode = asVector<int>(sscalemode);
}
SEXP sscalespec = getListElement(y, "sepScaleSpec");
if(!Rf_isNull(sscalespec)){
  RObjectTestExpectedType(sscalespec, &Rf_isNumeric, "sepScaleSpec");
  (*this)(i).sepScaleSpec = asVector<int>(sscalespec);
}
```

For the example, all six fields are present.

## 11. Existing `allterms_nll()` Splits Random Effects And Theta

Location:

```text
glmmTMB/src/glmmTMB.cpp
```

Marker: `[EXISTING]`

Local explanation:

```text
What this code does:
  Slices the global random-effect vector and theta vector term by term, then
  calls termwise_nll().

Why it matters here:
  The separable term receives U as a blockSize x blockReps array and theta as
  only that term's covariance parameters.

Why this existing path is reused:
  Once getReStruc() has set blockSize/blockReps/blockNumTheta correctly,
  separable terms fit the existing termwise likelihood framework.
```

```cpp
template <class Type>
Type allterms_nll(vector<Type> &u, vector<Type> theta,
                  vector<per_term_info<Type> >& terms,
                  bool do_simulate = false) {
  Type ans = 0;
  int upointer = 0;
  int tpointer = 0;
  int nr, np = 0, offset;
  for(int i=0; i < terms.size(); i++){
    nr = terms(i).blockSize * terms(i).blockReps;
    bool emptyTheta = ( terms(i).blockNumTheta == 0 );
    offset = ( emptyTheta ? -np : 0 );
    np     = ( emptyTheta ?  np : terms(i).blockNumTheta );
    vector<int> dim(2);
    dim << terms(i).blockSize, terms(i).blockReps;
    array<Type> useg( &u(upointer), dim);
    vector<Type> tseg = theta.segment(tpointer + offset, np);
    ans += termwise_nll(useg, tseg, terms(i), do_simulate);
    upointer += nr;
    tpointer += terms(i).blockNumTheta;
  }
  return ans;
}
```

For the separable term:

```text
U rows    = M * T
U columns = number of group levels
theta     = member log-SDs, member us correlations, AR(1) parameter
```

This existing function calls the new separable branch in `termwise_nll()`.

## 12. C++ Metadata Check

Marker: `[NEW/MODIFIED]`

Local explanation:

```text
What this code does:
  Verifies that R supplied a complete separable metadata contract.

Why it exists:
  Bad metadata would otherwise cause wrong reshaping or out-of-bounds parameter
  parsing in C++.

Why this implementation:
  It checks only the C++ contract; higher-level formula validity is already
  checked in R.
```

```cpp
template <class Type>
void check_separable_metadata(per_term_info<Type>& term) {
  if (term.sepDims.size() != 2 || term.sepCodes.size() != 2 ||
      term.sepDensityKinds.size() != 2 || term.sepDispatch.size() != 1 ||
      term.sepScaleMode.size() != 1 || term.sepScaleSpec.size() != 1)
    error("separable covariance structure is missing margin metadata");
  if (term.sepDims(0) * term.sepDims(1) != term.blockSize)
    error("separable dimensions do not match block size");
  if (term.sepScaleMode(0) != margin_sep_scale)
    error("separable covariance currently supports only margin scale");
  if (term.sepScaleSpec(0) < 0 || term.sepScaleSpec(0) > 1)
    error("separable margin scale index is out of range");
}
```

For the example, this confirms:

```text
two dimensions
two covariance codes
two density kinds
one dispatch
one scale mode
one scale spec
M * T == blockSize
scale mode is margin
scale margin is coordinate 0 or 1
```

## 13. C++ Parsed Parameter Struct

Marker: `[NEW/MODIFIED]`

Local explanation:

```text
What this code does:
  Names the parsed pieces needed by the dense x AR(1) evaluator.

Why it exists:
  termwise_nll() needs a clean handoff between theta parsing, density
  construction, likelihood evaluation, and reporting.

Why this implementation:
  A small struct keeps the separable branch readable without creating a full
  covariance-class hierarchy.
```

```cpp
template <class Type>
struct sep_dense_ar1_pars {
  int dense_margin;
  int ar1_margin;
  int dense_code;
  Type phi;
  vector<Type> sd;
  matrix<Type> dense_corr;
  vector<Type> us_corr_params;
};
```

For `us x ar1`:

```text
dense_margin = 0
ar1_margin   = 1
dense_code   = us_covstruct
phi          = transformed AR(1) parameter
sd           = member SD vector
us_corr_params = TMB unstructured-correlation parameters
dense_corr   = filled later only for reporting via dense_density.cov()
```

## 14. C++ Theta Parser For `us x ar1`

Marker: `[NEW/MODIFIED]`

Local explanation:

```text
What this code does:
  Converts the flat theta vector into member SDs, member correlation parameters,
  and AR(1) phi.

Why it exists:
  TMB optimizes unconstrained theta values, but density objects need transformed
  SD/correlation parameters.

Why this implementation:
  C++ consumes theta in the same margin order that R used to count it, making
  reversed margin syntax possible without hidden reordering.
```

```cpp
bool is_dense_corr_margin(int code) {
  // This is the small generalization added after review: cs, homcs, and us are
  // all dense-correlation margins from the separable evaluator's point of view.
  return code == cs_covstruct || code == homcs_covstruct || code == us_covstruct;
}

template <class Type>
matrix<Type> compound_symmetry_corr(int n, Type corr_transf) {
  // Shared correlation builder for cs and homcs.  They differ in their scale
  // parameters, not in their correlation matrix.
  Type a = Type(1) / (Type(n) - Type(1));
  Type rho = invlogit(corr_transf) * (Type(1) + a) - a;
  matrix<Type> corr(n, n);
  for (int i = 0; i < n; i++)
    for (int j = 0; j < n; j++)
      corr(i, j) = (i == j ? Type(1) : rho);
  return corr;
}

template <class Type>
void parse_separable_dense_margin(int code, int n, vector<Type> theta,
                                  int& theta_pos, bool scale_here,
                                  vector<Type>& sd, matrix<Type>& corr,
                                  vector<Type>& us_corr_params) {
  // Consume exactly one dense margin's theta segment.  The caller has already
  // decided which margin carries scale, so this helper only reads SDs when
  // scale_here is true.
  if (!is_dense_corr_margin(code))
    error("unsupported dense margin for separable covariance structure");

  if (scale_here) {
    if (code == homcs_covstruct) {
      sd.resize(n);
      sd.fill(exp(theta(theta_pos++)));
    } else {
      vector<Type> logsd = theta.segment(theta_pos, n);
      theta_pos += n;
      sd = exp(logsd);
    }
  }

  if (code == cs_covstruct || code == homcs_covstruct) {
    corr = compound_symmetry_corr(n, theta(theta_pos++));
  } else {
    int n_corr = n * (n - 1) / 2;
    us_corr_params = theta.segment(theta_pos, n_corr);
    theta_pos += n_corr;
  }
}

template <class Type>
sep_dense_ar1_pars<Type> parse_separable_dense_ar1(vector<Type> theta,
                                                   per_term_info<Type>& term) {
  // Convert the flat theta vector for a separable dense x AR(1) term into
  // margin-specific objects.  Theta is consumed in sepgrid/margin order, which
  // mirrors the R-side parameter counter:
  //
  //   cs    + ar1:  cs_log_sd..., cs_corr, ar1_phi
  //   homcs + ar1:  log_sd, homcs_corr, ar1_phi
  //   ar1 + homcs:  ar1_phi, log_sd, homcs_corr
  //   us    + ar1:  us_log_sd..., us_corr..., ar1_phi
  //   ar1  + us:    ar1_phi, us_log_sd..., us_corr...
  //
  sep_dense_ar1_pars<Type> out;
  out.dense_margin = -1;
  out.ar1_margin = -1;

  for (int m = 0; m < 2; m++) {
    if (term.sepDensityKinds(m) == dense_corr_sep) out.dense_margin = m;
    if (term.sepDensityKinds(m) == ar1_sep) out.ar1_margin = m;
  }
  // The current evaluator is deliberately narrow: one dense margin supplies the
  // flexible cross-member correlation and one AR(1) margin supplies serial
  // correlation.  R validates that unsupported pairs never reach this point, but
  // keep the check here so malformed metadata fails at the C++ boundary.
  if (out.dense_margin < 0 || out.ar1_margin < 0 ||
      out.dense_margin == out.ar1_margin)
    error("separable covariance currently requires one dense margin and one AR1 margin");
  if (term.sepScaleSpec(0) != out.dense_margin)
    error("separable covariance currently requires scale on the dense margin");

  int n_dense = term.sepDims(out.dense_margin);
  out.dense_code = term.sepCodes(out.dense_margin);
  out.phi = Type(0);
  out.sd.resize(n_dense);
  out.us_corr_params.resize(0);

  int theta_pos = 0;
  for (int m = 0; m < 2; m++) {
    int code = term.sepCodes(m);
    int n = term.sepDims(m);
    bool scale_here = (term.sepScaleSpec(0) == m);

    if (term.sepDensityKinds(m) == ar1_sep) {
      Type corr_transf = theta(theta_pos++);
      // Same unconstrained-to-correlation transform used by existing glmmTMB
      // AR1 and by `get_cor()` for the single-correlation case.
      out.phi = corr_transf / sqrt(Type(1) + pow(corr_transf, 2));
    } else if (term.sepDensityKinds(m) == dense_corr_sep) {
      // cs, homcs, and us all enter the same dense-margin parser.  This is why
      // adding cs did not require a new separable likelihood evaluator.
      parse_separable_dense_margin(code, n, theta, theta_pos, scale_here,
                                   out.sd, out.dense_corr, out.us_corr_params);
    } else {
      error("unsupported dense margin for separable covariance structure");
    }
  }
  if (theta_pos != theta.size())
    error("separable covariance theta parsing mismatch");

  return out;
}
```

For the canonical example:

```text
theta[1:M]                         -> member log-SDs
theta[(M+1):(M + M*(M-1)/2)]       -> us correlation parameters
theta[last]                        -> AR(1) transformed correlation
```

The AR(1) transform is:

```text
phi = theta / sqrt(1 + theta^2)
```

For `us`, TMB's `UNSTRUCTURED_CORR_t` later interprets `us_corr_params`.

Guide-only reading notes:

| parser variable | toy `us(member) %x% ar1(time)` value | meaning |
|---|---|---|
| `out.dense_margin` | `0` | coordinate 0, `member`, is the dense-correlation margin |
| `out.ar1_margin` | `1` | coordinate 1, `time`, is the AR(1) margin |
| `n_dense` | `2` | the dense margin has 2 member levels |
| `theta_pos` before `us` | `0` | parser starts at first theta slot |
| `theta_pos` after `us` log-SDs | `2` | two member SD parameters consumed |
| `theta_pos` after `us` correlation | `3` | one member correlation parameter consumed |
| `theta_pos` after `ar1` | `4` | one AR(1) parameter consumed; parser should be done |

For the default toy order, the parser consumes theta like this:

| parser step | theta slice | stored as |
|---:|---|---|
| 1 | `theta[0:1]` | `out.sd = c(sd_m1, sd_m2)` |
| 2 | `theta[2]` | `out.us_corr_params` |
| 3 | `theta[3]` | `out.phi` |

For reversed syntax, the same parser loop consumes the AR(1) parameter first
because it walks margin order:

```text
separable(ar1(0 + time) %x% us(0 + member) | group)
theta = ar1_phi, member_logsd..., member_corr...
```

This is why R and C++ must agree on margin order.  R counts theta in margin
order, then C++ consumes theta in the same margin order.

## 15. C++ Separable Likelihood

Marker: `[NEW/MODIFIED]`

Local explanation:

```text
What this code does:
  Reshapes each group's flat random effects into a 2D array, standardizes by
  member SDs, and evaluates a TMB separable density.

Why it exists:
  This is the efficient likelihood path for the implied Kronecker covariance.

Why this implementation:
  TMB SEPARABLE() avoids full covariance factorization; manual scaling keeps
  the current margin-scale mode explicit.
```

```cpp
template <class Type, class DenseDensity>
Type separable_dense_ar1_nll(array<Type> &U, vector<Type> sd, Type phi,
                             DenseDensity dense_density,
                             per_term_info<Type>& term,
                             int dense_margin, int ar1_margin) {
  // Evaluate one separable dense-correlation x AR(1) random-effect block.
  //
  // Inputs:
  //   U:
  //     matrix with rows = block entries and columns = grouping levels.
  //     This is exactly how `allterms_nll()` passes every random-effect term.
  //     With sepDims = c(2, 3) and two groups, the shape is:
  //
  //       rows:    m1_t1, m2_t1, m1_t2, m2_t2, m1_t3, m2_t3
  //       columns: g1, g2
  //
  //   sd:
  //     standard deviations along the scale-carrying margin.  For homcs this
  //     vector has one common value repeated for every level; for us it has one
  //     value per level.  The AR(1) margin deliberately has no scale inside
  //     separable(); otherwise the product scale would be non-identifiable.
  //
  //   phi:
  //     AR(1) correlation on the AR(1) margin.
  //
  //   dense_density:
  //     a standardized zero-mean Gaussian density for the dense correlation
  //     margin.  The AR(1) margin is standardized too, so cell scales are
  //     applied manually below.
  //
  // The density being evaluated is:
  //
  //   u_g ~ N(0, D_cell (R_margin2 %x% R_margin1) D_cell)
  //
  // for every grouping level g.  We avoid constructing the full Kronecker matrix
  // by reshaping the flat vector to an array and using TMB's `SEPARABLE`.
  //
  // The two supported orderings are both evaluated here:
  //
  //   sepgrid(dense, ar1):  SEPARABLE(AR1(phi), dense_density)(z)
  //   sepgrid(ar1, dense):  SEPARABLE(dense_density, AR1(phi))(z)
  //
  // TMB's SEPARABLE arguments are supplied in reverse array-dimension order.
  int n0 = term.sepDims(0);
  int n1 = term.sepDims(1);
  vector<int> dim(2);
  dim << n0, n1;
  Type ans = 0;
  int scale_margin = term.sepScaleSpec(0);

  for (int g = 0; g < term.blockReps; g++) {
    array<Type> z(dim);
    Type logscale = 0;
    // Convert the flat block for group g into an array in sepgrid() order.
    // For sepgrid(member, time), `z(i0, i1)` is a member-by-time table:
    //
    //              time 1   time 2   time 3
    //     member 1  z(0,0)  z(0,1)  z(0,2)
    //     member 2  z(1,0)  z(1,1)  z(1,2)
    //
    // `sepgrid()` and the R metadata use first-coordinate-fastest order:
    //   k = i0 + n0 * i1
    //
    // We divide by the scale-margin SD here so that `z` is on the standardized
    // correlation scale expected by `dense_density` and `AR1(phi)`.
    //
    // The log-Jacobian is the sum of log SDs for all cells.  This is the same
    // adjustment performed by TMB's VECSCALE/SCALE helpers, but manual scaling
    // keeps the multidimensional separable case explicit and easy to inspect.
    for (int i1 = 0; i1 < n1; i1++) {
      for (int i0 = 0; i0 < n0; i0++) {
        int k = i0 + n0 * i1;
        int si = (scale_margin == 0 ? i0 : i1);
        z(i0, i1) = U(k, g) / sd(si);
        logscale += log(sd(si));
      }
    }
    if (dense_margin == 0 && ar1_margin == 1) {
      ans += density::SEPARABLE(density::AR1(phi), dense_density)(z) + logscale;
    } else if (dense_margin == 1 && ar1_margin == 0) {
      ans += density::SEPARABLE(dense_density, density::AR1(phi))(z) + logscale;
    } else {
      error("invalid separable dense/AR1 margin order");
    }
  }
  return ans;
}
```

For canonical `us(member) %x% ar1(time)`:

```text
n0 = M
n1 = T
dense_margin = 0
ar1_margin   = 1
scale_margin = 0
```

The flat random-effect block is reshaped by:

```text
k = member_index + M * time_index
```

Each cell is standardized:

```text
z(member, time) = U[k, group] / sd_member[member]
```

The log-density adjustment:

```text
sum(log(sd_member[member])) over all cells
```

is added manually.  This is the Gaussian scale Jacobian.

TMB's `SEPARABLE()` arguments are supplied in reverse array-dimension order:

```cpp
density::SEPARABLE(density::AR1(phi), dense_density)(z)
```

This evaluates:

```text
z ~ N(0, R_time %x% R_member)
```

without building the full Kronecker covariance in the likelihood path.

Guide-only reading notes for the inner loop:

| loop variable | toy values | meaning |
|---|---|---|
| `g` | `0`, then `1` | group column: first `g1`, then `g2` |
| `i0` | `0`, `1` | coordinate 0: member index |
| `i1` | `0`, `1`, `2` | coordinate 1: time index |
| `k = i0 + n0 * i1` | `0..5` | flat row inside the group column |
| `si` | `i0` for default order | which member SD to divide by |
| `z(i0, i1)` | standardized random effect | original `U(k,g)` divided by member SD |

For group `g1`, the loop reads `U` like this:

| `i1` time | `i0` member | `k` | read from `U(k, 0)` | write to `z(i0, i1)` | divide by |
|---:|---:|---:|---|---|---|
| 0 | 0 | 0 | `u_g1_m1_t1` | `z(m1,t1)` | `sd_m1` |
| 0 | 1 | 1 | `u_g1_m2_t1` | `z(m2,t1)` | `sd_m2` |
| 1 | 0 | 2 | `u_g1_m1_t2` | `z(m1,t2)` | `sd_m1` |
| 1 | 1 | 3 | `u_g1_m2_t2` | `z(m2,t2)` | `sd_m2` |
| 2 | 0 | 4 | `u_g1_m1_t3` | `z(m1,t3)` | `sd_m1` |
| 2 | 1 | 5 | `u_g1_m2_t3` | `z(m2,t3)` | `sd_m2` |

Then the same table repeats for group `g2`, except the code reads from
`U(k, 1)`.

## 16. C++ Reporting Objects

Marker: `[NEW/MODIFIED]`

Local explanation:

```text
What this code does:
  Builds the standard deviation vector and optional full correlation matrix for
  VarCorr()/tests.

Why it exists:
  Users and tests need interpretable reported covariance components.

Why this implementation:
  Full matrix construction is acceptable for reporting/testing and is kept out
  of the likelihood path.
```

```cpp
template <class Type>
void report_separable_dense_ar1(vector<Type> sd, matrix<Type> dense_corr,
                                Type phi, per_term_info<Type>& term,
                                int dense_margin, int ar1_margin) {
  // Build report objects for VarCorr().
  //
  // This is not the final desired user-facing display.  For the prototype we
  // report the full standard deviation vector and, when requested, the full
  // Kronecker correlation matrix.  That makes correctness tests straightforward:
  //
  //   reported corr == kronecker(R_slowest_margin, R_fastest_margin)
  //
  // A later cleanup should replace the default printed output with a compact
  // component display to avoid printing very large separable matrices.
  int n0 = term.sepDims(0);
  int n1 = term.sepDims(1);
  int n = n0 * n1;
  int scale_margin = term.sepScaleSpec(0);
  term.sd.resize(n);
  // Repeat each margin SD across the other coordinate.  Which coordinate this
  // is depends on the user's sepgrid() order.
  for (int i1 = 0; i1 < n1; i1++) {
    for (int i0 = 0; i0 < n0; i0++) {
      int k = i0 + n0 * i1;
      int si = (scale_margin == 0 ? i0 : i1);
      term.sd(k) = sd(si);
    }
  }
  if (term.fullCor == 1) {
    // Construct the full correlation only for reporting/testing.  This should
    // not be used in the likelihood path for long time series.
    term.corr.resize(n, n);
    for (int b1 = 0; b1 < n1; b1++) {
      for (int a1 = 0; a1 < n0; a1++) {
        int k1 = a1 + n0 * b1;
        for (int b2 = 0; b2 < n1; b2++) {
          for (int a2 = 0; a2 < n0; a2++) {
            int k2 = a2 + n0 * b2;
            // Correlation between two flattened cells is the product of the
            // correlation along coordinate 0 and the correlation along coordinate 1.
            // For sepgrid(member, time), that is:
            //
            //   corr((member a1, time b1), (member a2, time b2))
            //     = corr_member(a1, a2) * corr_time(b1, b2)
            Type c0 = (dense_margin == 0 ?
                       dense_corr(a1, a2) :
                       pow(phi, abs(a1 - a2)));
            Type c1 = (dense_margin == 1 ?
                       dense_corr(b1, b2) :
                       pow(phi, abs(b1 - b2)));
            term.corr(k1, k2) = c0 * c1;
          }
        }
      }
    }
  } else {
    // Follow the existing compact-report convention for structured covariance
    // terms that do not store the whole correlation matrix.
    term.corr.resize(1,1);
    term.corr(0,0) = NAN;
  }
}
```

For canonical `us x ar1`, this reports:

```text
term.sd:
  member SDs repeated across time

term.corr:
  kronecker(R_time, R_member), if fullCor == 1
```

This full matrix is for reporting/testing only.  It is not used for likelihood
evaluation.

Guide-only reading notes for `term.corr(k1, k2) = c0 * c1`:

For default `sepgrid(member, time)`:

```text
coordinate 0 = member
coordinate 1 = time
dense_margin = 0
ar1_margin   = 1
```

So:

| code expression | default meaning |
|---|---|
| `a1`, `a2` | two member indices |
| `b1`, `b2` | two time indices |
| `k1 = a1 + n0*b1` | flat row for first member-time cell |
| `k2 = a2 + n0*b2` | flat row for second member-time cell |
| `c0 = dense_corr(a1, a2)` | member correlation |
| `c1 = pow(phi, abs(b1 - b2))` | time correlation |
| `term.corr(k1, k2) = c0 * c1` | product correlation for the two cells |

Concrete cells:

| first cell | second cell | `c0` member part | `c1` time part | reported correlation |
|---|---|---|---|---|
| `m1_t1` | `m2_t1` | `rho` | `1` | `rho` |
| `m1_t1` | `m1_t2` | `1` | `phi` | `phi` |
| `m1_t1` | `m2_t2` | `rho` | `phi` | `rho * phi` |
| `m1_t1` | `m2_t3` | `rho` | `phi^2` | `rho * phi^2` |

## 17. C++ `termwise_nll()` Separable Branch

Marker: `[NEW/MODIFIED]`

Local explanation:

```text
What this code does:
  Adds the separable covariance branch to termwise_nll().

Why it exists:
  Existing covariance branches do not know how to combine two marginal density
  objects with TMB SEPARABLE().

Why this implementation:
  R has already validated/encoded the term, so C++ only checks metadata, parses
  theta, constructs density objects, dispatches, and reports.
```

```cpp
else if (term.blockCode == separable_covstruct) {
  if (do_simulate)
    error("simulation is not yet implemented for separable covariance structures");

  check_separable_metadata(term);
  switch (term.sepDispatch(0)) {
  case dense_ar1_dispatch: {
    sep_dense_ar1_pars<Type> sep = parse_separable_dense_ar1(theta, term);
    if (sep.dense_code == cs_covstruct || sep.dense_code == homcs_covstruct) {
      density::MVNORM_t<Type> dense_density(sep.dense_corr);
      ans += separable_dense_ar1_nll(U, sep.sd, sep.phi, dense_density, term,
                                     sep.dense_margin, sep.ar1_margin);
      DISABLE_AD {
        report_separable_dense_ar1(sep.sd, sep.dense_corr, sep.phi, term,
                                   sep.dense_margin, sep.ar1_margin);
      }
    } else if (sep.dense_code == us_covstruct) {
      density::UNSTRUCTURED_CORR_t<Type> dense_density(sep.us_corr_params);
      ans += separable_dense_ar1_nll(U, sep.sd, sep.phi, dense_density, term,
                                     sep.dense_margin, sep.ar1_margin);
      DISABLE_AD {
        report_separable_dense_ar1(sep.sd, dense_density.cov(), sep.phi, term,
                                   sep.dense_margin, sep.ar1_margin);
      }
    } else {
      error("unsupported dense margin for separable covariance structure");
    }
    break;
  }
  default:
    error("unsupported separable covariance dispatch");
  }
}
```

For the example:

```text
term.blockCode       = separable_covstruct
term.sepDispatch(0)  = dense_ar1_dispatch
sep.dense_code       = us_covstruct
```

So C++ takes the `us_covstruct` branch:

```cpp
density::UNSTRUCTURED_CORR_t<Type> dense_density(sep.us_corr_params);
```

That dense member correlation density is combined with `AR1(phi)` by
`separable_dense_ar1_nll()`.

`DISABLE_AD` is used for reporting because `VarCorr()` does not need automatic
differentiation through the reported matrices.

## 18. Last Touch: Existing Reporting Pipeline Reads `term.sd` And `term.corr`

Marker: `[EXISTING]`

Local explanation:

```text
What this code does:
  Existing report/VarCorr machinery consumes term.sd and term.corr.

Why it matters here:
  The separable branch only has to fill the same report fields other covariance
  structures already use.

Why this existing path is reused:
  It avoids adding a separate public reporting pipeline for the prototype.
```

After `termwise_nll()` fills:

```text
term.sd
term.corr
```

the existing glmmTMB report/`VarCorr()` pipeline uses those fields just like it
does for other covariance structures.

This branch's last separable-specific runtime touch is therefore:

```cpp
report_separable_dense_ar1(...)
```

and the `termwise_nll()` separable branch that calls it.

## End-To-End Data Shape Summary

For:

```r
separable(us(0 + member) %x% ar1(0 + time) | group)
```

with `M` members, `T` time points, and `G` groups:

```text
R formula rewrite:
  product syntax -> separable(sepgrid(member, time) + 0 | group, list(...))

model frame:
  sepgrid(member, time) factor with M*T levels

mkReTrms:
  blockSize = M*T
  blockReps = G
  cnms = "(1,1)", "(2,1)", ..., "(M,T)"

getReStruc:
  blockCode = separable
  blockNumTheta = M + M*(M - 1)/2 + 1
  sepDims = c(M, T)
  sepCodes = c(us, ar1)
  sepDensityKinds = c(dense_corr, ar1)
  sepDispatch = dense_ar1
  sepScaleMode = margin
  sepScaleSpec = 0

C++ theta parse:
  theta[1:M] = member log-SDs
  next M*(M-1)/2 = unstructured member correlation parameters
  last = AR(1) transformed correlation

C++ likelihood:
  reshape U[, g] to z[M, T]
  divide by member SDs
  evaluate SEPARABLE(AR1(phi), UNSTRUCTURED_CORR(us_corr))(z)
  add log-scale Jacobian

C++ reporting:
  report member SDs repeated over time
  report full kronecker correlation if fullCor == 1
```

## Why The Current Implementation Is Hacky But Reviewable

The hackiest part is the generated metadata call:

```r
list(grid = ..., margins = ..., scale = ...)
```

inside the rewritten `separable(...)` call.

Why it exists:

```text
reformulas::splitForm() already keeps additional arguments attached to the
specific random-effect term in reTrmAddArgs.
```

Why that is useful:

```text
separable(...) + separable(...) + separable(...)
```

does not need a separate matching system.  Each term carries its own metadata.

Why it is not the final ideal:

```text
The formula is being used as a metadata transport.  A future reformulas/glmmTMB
API could return structured nested covariance metadata directly.
```

What a cleaner future implementation would likely change:

```text
- keep building the random-effect design before getReStruc();
- represent product covariance metadata as a first-class parser result;
- keep term-local alignment guaranteed by the parser;
- later add marginal-design compilation for richer margins.
```

What should not change casually:

```text
The transformation to a complete grid must happen before mkReTrms(), unless
glmmTMB gets a new random-effect design builder for separable product terms.
```

## Appendix A: Small Support Code Written By This Branch

The main walkthrough above follows the fit/likelihood/reporting path.  This
appendix covers the remaining small runtime/support hooks written by the branch
that are relevant to the same `us x ar1` feature.

### A.1 R Covariance Enum

Location:

```text
glmmTMB/R/enum.R
```

Marker: `[NEW/MODIFIED]`

```r
.valid_covstruct <- c(
  diag = 0,
  us   = 1,
  cs   = 2,
  ar1  = 3,
  ou   = 4,
  exp = 5,
  gau = 6,
  mat = 7,
  toep = 8,
  rr = 9,
  homdiag = 10,
  propto = 11,
  hetar1 = 12,
  homcs = 13,
  homtoep = 14,
  equalto = 15,
  separable = 16
)
```

Explanation:

```text
This makes "separable" a recognized covariance-structure class on the R side.
The value 16 must match separable_covstruct in C++.
```

For the example, `splitForm()` classifies the term as `"separable"`, and
`getReStruc()` later converts that to `blockCode = 16`.

### A.2 Export `sepgrid()`

Location:

```text
glmmTMB/NAMESPACE
```

Marker: `[NEW/MODIFIED]`

```text
export(sepgrid)
```

Explanation:

```text
The public product syntax is preferred, but sepgrid() is still a real helper
used in formulas after rewriting and useful for explicit/debugging workflows.
```

### A.3 Roxygen Detail For `sepgrid()`

Location:

```text
glmmTMB/R/utils_covstruct.R
```

Marker: `[NEW/MODIFIED]`

```r
##' \code{sepgrid} is similar to \code{numFactor}, but creates levels for
##' the complete Cartesian product of the supplied coordinate levels.  Factor
##' inputs preserve unused levels; non-factor inputs use sorted observed
##' values.  Use factors with explicit levels when globally unobserved cells
##' are part of the intended separable grid, for example an unobserved day in
##' an AR(1) time series.
```

Explanation:

```text
This documents the most important user-facing footgun: if globally unobserved
intended cells matter, the coordinate must be a factor with explicit levels.
```

### A.4 `predict(newdata = ...)` Guard

Location:

```text
glmmTMB/R/predict.R
```

Marker: `[NEW/MODIFIED]`

```r
has_separable <- any(vapply(object$modelInfo$reStruc, function(x) {
  any(vapply(x, function(y) {
    identical(unname(y$blockCode), unname(.valid_covstruct[["separable"]]))
  }, logical(1)))
}, logical(1)))
if (!is.null(newdata) && has_separable) {
  stop("predict() with newdata is not yet implemented for separable covariance structures. ",
       "The separable margin coordinates must be encoded against the fitted grid; ",
       "use prediction on the original data for now.")
}
```

Explanation:

```text
Prediction on the original data can use the existing fitted object.
Prediction with newdata would need to encode member/time coordinates against
the fitted sepgrid levels.  That contract is not implemented yet, so the branch
fails explicitly instead of silently rebuilding a different latent grid.
```

This is not part of the fitting path, but it is a runtime touchpoint for fitted
models containing the example separable term.

## Appendix B: Code Written But Not Duplicated Here

This walkthrough intentionally does not duplicate every changed documentation,
test, or repository-support line.  Those are still part of the branch, but they
are not runtime code that the `us x ar1` term passes through while
fitting/evaluating.

Review them separately:

```text
glmmTMB/tests/testthat/test-separable.R
glmmTMB/man/glmmTMB.Rd
glmmTMB/man/numFactor.Rd
glmmTMB/man/predict.glmmTMB.Rd
.gitignore
```

The tests are especially important because they verify:

```text
- product syntax rewrite;
- full sepgrid level preservation;
- us x ar1 likelihood equality against dense us();
- reversed margin order;
- unsupported margin/scale errors;
- predict(newdata) and simulation limitations.
```

The two test snippets most useful for understanding the code shape are included
below.

### B.1 Test Code: Kronecker Order Follows `sepgrid()` Order

This excerpt is from `make_sep_case()` in
`glmmTMB/tests/testthat/test-separable.R`.

```r
if (reversed) {
    ## Reversed order checks that the implementation follows the user's
    ## sepgrid/product order rather than assuming member is always first.
    form <- switch(struc,
        cs = y ~ 1 + separable(ar1(0 + time) %x% cs(0 + member) | group),
        homcs = y ~ 1 + separable(ar1(0 + time) %x% homcs(0 + member) | group),
        us = y ~ 1 + separable(ar1(0 + time) %x% us(0 + member) | group)
    )
    dense_form <- y ~ 1 + us(sepgrid(time, member) + 0 | group)
    theta <- c(ar1_to_theta(phi), member_theta)
    ## `sepgrid(time, member)` makes time the fastest coordinate.  The dense
    ## covariance in that flattened order is therefore member blocks, each
    ## containing the full time correlation.
    R_full <- kronecker(R_member, R_time)
    sd_full <- if (length(sd) == 1) rep(sd, n_member * n_time) else rep(sd, each = n_time)
    codes <- unname(c(.valid_covstruct[["ar1"]], .valid_covstruct[[struc]]))
    kinds <- c(2L, 1L)
    scale_spec <- 1L
} else {
    ## Default order used in the walkthrough:
    ## rows inside each group are member-fastest within time:
    ## m1_t1, m2_t1, m1_t2, m2_t2, ...
    form <- switch(struc,
        cs = y ~ 1 + separable(cs(0 + member) %x% ar1(0 + time) | group),
        homcs = y ~ 1 + separable(homcs(0 + member) %x% ar1(0 + time) | group),
        us = y ~ 1 + separable(us(0 + member) %x% ar1(0 + time) | group)
    )
    dense_form <- y ~ 1 + us(sepgrid(member, time) + 0 | group)
    theta <- c(member_theta, ar1_to_theta(phi))
    ## `sepgrid(member, time)` makes member the fastest coordinate.  The dense
    ## covariance in that flattened order is time blocks, each containing the
    ## full member correlation.
    R_full <- kronecker(R_time, R_member)
    sd_full <- if (length(sd) == 1) rep(sd, n_member * n_time) else rep(sd, n_time)
    codes <- unname(c(.valid_covstruct[[struc]], .valid_covstruct[["ar1"]]))
    kinds <- c(1L, 2L)
    scale_spec <- 0L
}
```

Plain interpretation:

| product syntax | flat cell order | dense comparison matrix |
|---|---|---|
| `us(member) %x% ar1(time)` | all members within each time | `kronecker(R_time, R_member)` |
| `ar1(time) %x% us(member)` | all times within each member | `kronecker(R_member, R_time)` |

This is the same order rule used by the C++ likelihood and reporting code:

```text
flat index = coordinate0 + n_coordinate0 * coordinate1
```

### B.2 Test Code: Dense Likelihood Comparison Uses Two Groups

This excerpt is from `expect_separable_dense_nll()` in
`glmmTMB/tests/testthat/test-separable.R`.

```r
expect_separable_dense_nll <- function(case) {
    ## Use two grouping levels so the dense comparison covers the normal
    ## glmmTMB layout: one flat member-time block per group.
    dd <- make_sep_dat(n_member = case$n_member, n_time = case$n_time, n_group = 2)
    dd$y <- 0
    ## `b` concatenates the group-specific random-effect blocks.  With the default
    ## 2 member x 3 time x 2 group case this is 12 values: 6 cells for group 1,
    ## followed by the same 6 cells for group 2.
    b <- seq(-0.4, 0.5, length.out = length(case$sd_full) * nlevels(dd$group))

    sep_nll <- joint_nll_at(case$form, dd, case$theta, b)
    dense_nll <- joint_nll_at(case$dense_form, dd, case$theta_dense, b)

    expect_equal(unname(sep_nll), unname(dense_nll), tolerance = 1e-6)
}
```

For the toy walkthrough dimensions:

```text
length(case$sd_full) = 6 cells per group
nlevels(dd$group)   = 2 groups
length(b)           = 12 random-effect values
```

The expected `b` layout is:

| position range | group | cells |
|---:|---|---|
| 1-6 | `g1` | `m1_t1, m2_t1, m1_t2, m2_t2, m1_t3, m2_t3` |
| 7-12 | `g2` | `m1_t1, m2_t1, m1_t2, m2_t2, m1_t3, m2_t3` |
