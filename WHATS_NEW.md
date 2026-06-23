# What's New — HaloMaker Python

A human-friendly summary of recent improvements to the Python + Fortran
HaloMaker / AdaptaHOP pipeline (RAMSES halo & galaxy finder), and of **GasMaker**,
the new gas post-processor built on top of it.

> Scope tags: **[released]** = committed to the main release line ·
> **[prototype]** = validated but staged in the dev workspace
> (`codex_eval/tmp/gasmaker/`), not yet folded into the release tree ·
> **[inherited]** = from the earlier port, listed for context (provenance
> confirmed with the original author — see `codex_eval/PL_TO_PRIOR_CODEX_QA.md`).

---

## 🛰 GasMaker — new gas post-processor  **[prototype]**

A brand-new tool that reads an existing halo/galaxy catalog plus the RAMSES
AMR/hydro data and computes **gas properties per halo and galaxy** that
HaloMaker alone does not produce.

- **Tier 1 (galaxy/halo gas):** total / cold (T<10⁴ K) / dense (cold & n_H>5/cc)
  gas masses, gas metallicity, gas kinematics, gas angular momentum, and
  per-element chemistry (H, O, Fe, Mg, C, N, Si, S, D) — measured inside the
  stellar radii (r\*, r50, r90) and the virial radius.
- **Tier 2 (halo overdensity):** spherical-overdensity radii and masses
  (r200/m200, r500/m500) and the enclosed dark-matter / star / gas masses at
  r_vir / r200 / r500.
- **Restartable:** a long run can be interrupted and resumed; finished halos are
  skipped and the output stays consistent.
- **Accurate cell–sphere overlap, made practical.** Measuring gas in a sphere
  needs the fraction of each boundary AMR cell that lies inside. The exact
  reference (subdivide every boundary cell to octree depth 8) is correct but
  explosively slow — on the hardest multi-halo case it took **≈1 h 41 m / ≈8.9 GB**.
  GasMaker's **adaptive** method reaches the *same* accuracy (within
  **0.04–0.32%** of the exact reference) in **≈2–3 min / ≈4 GB** by only
  subdividing a boundary cell while it is large relative to the aperture.
  *(This compares two GasMaker overlap strategies — the brute-force reference
  vs the adaptive design that made the tool usable — not a speedup over any
  earlier product.)*
- **Scientific validation in progress:** GasMaker outputs are being checked
  column-by-column against the independent **RUR** reference (calling RUR's own
  `_extend_*` functions with identical apertures) across a large stratified
  sample. Class A mass fields (gas/cold/dense at r_vir and r50/r90, plus DM and
  star masses) currently agree at **machine precision**; r200/r500 differ by a
  known, intended radius-definition choice and are reported separately.

---

## 🌟 HaloMaker — major changes

### Consistent output units  **[released, breaking — `halomaker_units_v2`]**
HaloMaker and GasMaker HDF5 outputs now share one unit system, removing a
1e11× / physical-vs-code-unit mismatch that could silently corrupt joins
(e.g. a gas fraction across the two files):
- masses → `Msun` (was `10^11 Msun`); positions / radii / shape axes → RAMSES
  code units `[0,1)`; angular momentum → `Msun Mpc km/s`; energies →
  `Msun (km/s)^2`; `rho_0` → `Msun/kpc^3`; velocities/SFR/age/metal unchanged.
- Conversion happens **only at the HDF5 write boundary** — internal physics is
  untouched, so values differ from the old output by an exact unit factor and
  nothing else (verified field-by-field; GasMaker results provably unchanged).
- Files are tagged `units_version="halomaker_units_v2"` with `box_*_mpc` and
  per-field `field_units` attrs; GasMaker reads the tag and stays correct on
  both old and new catalogs. The 07206 / NH2 goldens were re-frozen to v2.
- **Breaking:** analysis code assuming the old units must branch on
  `units_version`. Full table in `CATALOG_FORMAT.md`.

### Per-element stellar chemistry  **[released]**
The catalog now carries per-element **stellar** chemistry — `H_star, O_star,
Fe_star, Mg_star, C_star, N_star, Si_star, S_star, D_star` (mass fraction),
mass-weighted over the same stellar members as `metal` (catalog is now
81 fields). The Ra4 reader locates the `chem_*` star columns from the particle
descriptor (no positional assumption); snapshots without stellar chemistry
(e.g. Ra3) get `NaN` and run normally. Validated against the RUR reference
(reader exact per-particle; aggregation matches `metal`). The 9-element list is
shared with GasMaker's gas chemistry via a single `chem_species.py`.

### `dump_DMs` → `dump_members`  **[released, config rename]**
The member-dump flag is renamed to `dump_members` because it writes pos/vel/mass
for **all** members (DM + stars), not just DM. The dormant GalaxyMaker-era
`dump_stars` placeholder is retired. Legacy keys `dump_DMs` and `dump_stars` are
still accepted as aliases, so existing input files keep working.

### Full-box is now the single, clean mode  **[released]**
The legacy zoom-in code path has been fully removed (Python paths **and** the
Fortran interface), so the pipeline runs one well-tested full-box mode.
- Verified the removal changes **no scientific output** — the catalog is
  **bit-for-bit identical** before vs after (2,000 halos, zero regressions).
- Old `zoomin` / `zoombox` config keys are silently ignored for backward
  compatibility.

### Inherited core (the earlier port)  **[inherited]**
For context, these came from the original Python+Fortran port and underpin
everything above:
- **Combined DM + star analysis** (better density estimate, natural galaxy–halo
  connection) with **HDF5** catalog output instead of Fortran binary. This
  Python/f2py *Reference* line was itself a workflow modernization of the
  *Original*, which ran HaloMaker and GalaxyMaker as **separate executables**
  with Fortran-binary output and linked DM halos to stellar galaxies only after
  the fact — the Reference unified them into one structure-finding pass.
- **Extended halo/galaxy properties computed in-line** during the run (the
  `extend` feature): stellar center (px\*,py\*,pz\*), stellar radii (r\*, r50,
  r90), stellar angular momentum, mass-weighted age/metallicity, star-formation
  rate, V_max/R_max, NFW concentration + error, and the DM inner-density slope
  (`inslope`/`inslopeerr`). This is the same *family* as the RUR post-processing
  scripts, re-implemented as an inline subset (not an exact port; gas/hydro
  quantities are GasMaker's job).
  - **Stellar kinematics** are also computed in-line: `vrot`, `sig3d`, `sigcyl`,
    each for all member stars and within r50 / r90. `sig3d` is the
    bulk-subtracted Cartesian dispersion; `sigcyl` is measured in the stellar
    angular-momentum frame after removing the mean cylindrical motions.
  - **Star formation** is stored for both 100 Myr and 10 Myr windows
    (`SFR*` and `SFR10*`), each with total, r50, and r90 apertures.
- **Member storage** in HDF5: per-row `index`/`count` plus a flat `pids` array
  (and `pos`/`vel`/`mass` when `dump_DMs` is on).
- **SSP photometry** in three stellar-population models — **CB07** (Charlot &
  Bruzual 2007), **BC03** (Bruzual & Charlot 2003), and **FSPS** (via
  python-fsps) — all assuming a Chabrier IMF and SSP-normalized to 1 M⊙ formed.
  Magnitudes: SDSS *ugriz* (AB) and Johnson *UBV* + K (Vega). Luminosities use
  `initial_mass` when available (e.g. Ra4) and fall back to current stellar mass
  otherwise; the choice is recorded in the HDF5 `mass_source` metadata. Each SSP
  group also stores r-band luminosity-weighted `age_r`, `metal_r`, `r50_r`, and
  `r90_r`, row-aligned with `/catalog/halo`.

---

## 🔬 Correctness & scientific validation
- **Bit-exact regression gate.** Every code change is checked by re-running the
  pipeline and comparing the catalog to a frozen reference at **rtol=0, atol=0**
  (nan-aware). A change is accepted only if scientific values are identical
  (timing/memory/log order may differ).
- **Frozen golden references.** Two trusted catalogs are pinned as baselines:
  NH2 (Ra4; ~1.4×10⁸ particles, 15,386 halos/subhalos) and 07206 (Ra3).
- **Halo-center fix.** Corrected a center-determination offset so the reported
  center matches the member center-of-mass.
- **Adversarial cross-review.** The two implementers audit each other's work
  (prosecutor / defense / judge); this caught real unit bugs and report-vs-
  reality mismatches before they could reach the release.

---

## ⚡ Performance & memory (HaloMaker core)
- **Memory-lean AdaptaHOP core** **[inherited]** — the neighbor-search core
  avoids a large persistent neighbor table (it recomputes neighbors on demand
  with the same tree traversal) and parallelizes the heavy passes with OpenMP.
  Representative dev3 benchmarks: **07206 full-box 9.27 → 2.21 GiB (−76%)** and
  **39990 full-box 101.75 → 23.29 GiB (−77%)**, with walltime changing only a
  few percent. (A legacy zoom-in run showed a similar drop; the release now uses
  the unified full-box path.) The science algorithm is unchanged; the catalog is
  kept **bit-reproducible** by design (deterministic serial group-ID assignment;
  nodes normalized before export — increasing level first, then decreasing true
  mass within a level).
  - *Implementation detail (for maintainers):* earlier versions kept
    `iparneigh(nhop, npart)` alive across AdaptaHOP stages — with `nhop=100`
    that table dominated peak RSS. dev3 removes it: the density pass still uses
    the same octree neighbor search, but later **local-maxima and saddle stages
    recompute neighbors** instead of reading a stored table. A naive serial
    recomputation was too slow, so those heavy passes are **OpenMP-parallelized**
    while group-ID assignment, symmetry cleanup, compaction, and node export stay
    **deterministic**. In short, dev3 trades a huge long-lived neighbor array for
    parallel recomputation — which is why peak RSS drops ~76–77% while walltime
    stays within a few percent. *(Pitfalls to avoid: do not re-introduce the
    persistent neighbor table, and do not move group-ID assignment inside the
    parallel loop.)*
  - *Development lineage:* the Reference carried the original
    `iparneigh(nhop, npart)` table and stayed memory-heavy. **dev1** proved that
    storing neighbors for *active* particles only preserves the catalog (but
    needed a slow second pass); **dev2** stored them in one pass (faster, still a
    large `nhop × nactive` table); **dev3** removed the persistent table entirely
    and pays for lookups by OpenMP recomputation — the decisive step. (07206
    full-box across the ladder: 9.27 → 6.56 → 6.54 → 2.21 GiB.) So dev3 is *not*
    merely "the OpenMP-fast version" — the real change is dropping the stored
    neighbor table after dev1/dev2 narrowed it.
- **In-place periodicity correction** — `correct_for_periodicity_1d` now edits
  in place (masked) instead of allocating, reducing transient memory.
- **No import-time side effects** — importing the neighbor-tree module no longer
  spawns a multiprocessing manager process, fixing a portability hazard in
  sandboxed environments (and a small startup cost).

---

## 🧹 Micro-optimizations, cleanups & hardening  **[released]**
All verified output-neutral (07206 catalog bit-identical):
- Removed the legacy zoom-in Fortran module (~2,640 lines) and assorted dead
  code: commented-out legacy blocks, three unused `jacobi` implementations,
  superseded `rf()` helper, and dead imports — **~2,950 lines removed total**.
- Replaced an old hand-rolled input parser with a pointer to the live
  `PARAMS`-based parser.
- Small idioms: squared-distance comparison (drop a per-halo `sqrt`),
  `np.count_nonzero` / `np.full` / boolean-mask sums, dict-membership without
  `.keys()`, deduplicated imports.

Inherited release-hardening (context) **[inherited]**:
- The f2py build uses **explicit `.pyf` interfaces** so only public routines are
  exposed to Python — this avoids f2py failures from accidentally wrapping
  internal Fortran derived types, and makes signature changes (e.g.
  `sync_others`) auditable. (The current Fortran clean-up followed this same
  discipline: change the Fortran source, the `.pyf`, and the Python call
  together.)
- Small correctness/usability fixes: robust `#` / `!` comment handling in input
  files, a `sigma_dm` zero-division guard, `megaverbose`-only memory-tracking
  logs, and an unconditional AdaptaHOP particle relink that must not depend on
  `verbose`.

---

## 🛠️ How the work is organized
Development runs as a small "AI research team": a planner/reviewer lead plus two
**competing implementers** whose submissions are independently verified and
cross-reviewed, with all coordination in human-readable markdown. The bit-exact
regression gate and the cross-review are what let fast iteration coexist with
scientific reproducibility.

---

## ⚠️ Known caveats & things to watch
Honest notes (confirmed with the original author) for anyone running or
extending the pipeline:
- **Per-element chemistry (gas + stellar) is validated; aperture-expanded
  variants are not.** GasMaker **gas** chemistry on `r*` (TASK-11) and the
  catalog **stellar** chemistry `*_star` (TASK-13) were both validated against
  the RUR reference at machine precision, with the missing-chemistry fallback
  confirmed (`NaN`, schema preserved). The tracked 9 elements sum to ~0.752
  (remainder ~0.248 = He, consistent with RAMSES). Still open: aperture-resolved
  stellar chemistry (r50/r90 variants are not produced), and NH2 output-60
  format edge cases beyond the validated path.
- **`initial_mass` (Ra4) reader assumes a particle descriptor order.** It locates
  `initial_mass` relative to `metallicity` in `part_file_descriptor.txt`; a
  simulation with a different descriptor order must be verified (the code errors
  out on an unsupported order rather than guessing).
- **Cosmology from the RAMSES header is overwritten, but `H_f` stays an input
  value.** Mixing simulations with different Hubble parameters in one input list
  (or reusing the same output number across simulations) can produce wrong
  cosmology state or output-filename collisions.
- **Exact ties.** Pathological snapshots with many equal-density / equal-distance
  ties can make the deterministic ordering assumption fragile (the code logs a
  warning); compare membership first in that case.
- **Runtime cleanup.** A force-killed run can leave `forkserver` /
  `resource_tracker` / `/dev/shm` residue; `clean_runtime.sh` clears it.

---

## 📋 In progress / next
- All-root RUR scientific validation of GasMaker (Tier 1 + Tier 2). *Done at the
  stratified-sample level (machine-precision PASS); full all-root is the
  remaining optional cap.*
- ~~Folding the validated GasMaker prototype into the release tree.~~ **Done** —
  GasMaker is now a first-class release tool (`GasMaker.py` + `gasmaker/`).
- ~~Stabilize the DM-profile inner-slope fit (`inslope`) against platform/compiler
  floating-point drift.~~ **Done** — roundoff-degenerate shells are now dropped
  from the fit; the residual is within the field-policy tolerance, so the frozen
  goldens were not re-frozen.
- Streaming/chunked verification harness for very large runs (npart > 10⁹).
