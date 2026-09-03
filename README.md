# SKILL-DRIVEN FLEXIBILITY ASSESSMENT — Reference Implementation

This software quantifies how well a production resource functionally aligns with a set of product-process requirements. Resource functionalities and requirements are both represented as bounded regions in a shared, skill-specific property space, containing continuous, discrete and categorical properties as well as cross-property constraints. The software computes the hypervolume of each region, of their intersection, and of their difference, and returns three scalar metrics: compatibility (Λ), utilisation (Γ), and a geometric expansion effort ratio (Δ_E). Together these express how much of the requirement space a resource covers, how much of the resource's own functionality the requirements occupy, and how large the uncovered remainder is relative to what the resource already offers. The metrics are intended as technical screening information for manufacturing system (re)design: ranking brownfield retention or retrofit candidates, and comparing greenfield alternatives. This code is the executable form of the computational procedure specified in the appendix of the reference publication (section 9).

**Contact name**  [Lauren Van De Ginste](mailto:Lauren.VanDeGinste@UGent.be)  
**Authors**  [Lauren Van De Ginste](mailto:Lauren.VanDeGinste@UGent.be), [El-Houssaine Aghezzaf](mailto:ElHoussaine.Aghezzaf@UGent.be), [Johannes Cottyn](mailto:Johannes.Cottyn@UGent.be)  
**ZAP name**  [Johannes Cottyn](mailto:Johannes.Cottyn@UGent.be)  
**Website**  <http://www.ugent.be/m-f/en>, <https://www.isye.be/>

**Download**  <https://github.com/ugent-isye/Skill-DrivenFlexibilityAssessment> — archived version: Zenodo DOI `[TODO]`

The code was developed at **Flanders Make @ UGent (Belgium)** in the context of various research projects on flexibility assessment.

---

## 1. ABSTRACT

The transition towards sustainable and resilient production systems requires manufacturing systems that can accommodate evolving product-process requirements across multiple product generations. Current system (re)design practices nevertheless remain rigid and project-oriented: driven by short-term payback expectations, decision-makers frequently neglect the long-term value of flexibility, favouring equipment replacement over reconfiguration and reuse. The reference publication introduces a skill-driven assessment methodology that formally connects product-process requirements to the functional (re)use potential of production resources. It evaluates the alignment of available resource skill sets with evolving requirements through three scalar metrics: compatibility, utilisation, and a geometric expansion effort ratio. In contrast to the binary feasibility checks prevalent in skill-based matching, these metrics quantify the degree and direction of geometric resource-requirement alignment within skill-specific property spaces, and trace the underlying gaps to their origin.

This repository is the computational artefact of that publication. The paper's appendix, *Case Study Data and Computational Procedure*, specifies the computation as Algorithm 1; `reference_implementation.py` is written as a step-for-step image of that pseudocode. It carries one worked example from the industrial case study — resource NS_600 assessed against Scenario 4 of skill S1, cylindrical grinding — and asserts its own output against the values reported in the paper, so a correct installation is self-evident on the first run.

Reference publication:

> Van De Ginste, L., Aghezzaf, E.-H., & Cottyn, J. (2026). Hypervolume-based scalar metrics for skill-driven flexibility assessment in manufacturing system (re)design. *International Journal of Production Research* (UNDER REVIEW).

## 2. SOFTWARE DEPENDENCIES

**No commercial software is required.** Every dependency is open source and available through conda-forge.

| Component | Version used | Role |
| --- | --- | --- |
| Python | 3.11.15 | — |
| NumPy | 2.4.6 | array handling |
| polytope | 0.2.5 | polyhedral set operations (union, set difference) |
| Normaliz | 3.11.0 | exact volume of a polytope from its half-space representation |

These are the versions the reported output was developed and validated against. The script prints them alongside whatever is actually installed on every run, so a version mismatch is visible immediately as a possible cause of any discrepancy.

**Normaliz is an external executable**, not a Python package. The script locates it, in order, via the `NORMALIZ_EXE` environment variable, the system `PATH`, or the standard conda install location inside the active environment. Installing the supplied conda environment covers all three routes.

`environment.yml` pins the complete environment in which the reported output was produced. It is a Windows (win-64) export: it includes the `ucrt`, `vc14_runtime` and `vs2015_runtime` packages, which exist only for that platform, so `conda env create -f environment.yml` will not solve on Linux or macOS. Use it on Windows for an exact reproduction; elsewhere, create a minimal environment instead — the script itself needs only NumPy, polytope and the Normaliz executable:

```
conda create -n skill-based-flexibility-assessment -c conda-forge python=3.11 numpy polytope normaliz
```

## 3. THE FILES YOU SHOULD GET

| File | Purpose |
| --- | --- |
| `reference_implementation.py` | The computation, plus the worked example |
| `README.md` | This file |
| `environment.yml` | The conda environment the reported output was produced in |
| `LICENSE.txt` | GNU General Public License v3 |

The script is **self-contained**. It carries the worked example's inputs as literals, copied from Appendix Tables A1 and A3 of the publication, so no data files are needed to run it. There is no package to install and nothing to configure beyond the environment.

The full case study covers further resource-scenario pairs. Their inputs are given in the same appendix tables, which are the inputs of record, and the same `compute_metrics` call evaluates any of them.

## 4. GETTING STARTED

Create the environment once, then run the script:

```
conda env create -f environment.yml
conda activate skill-based-flexibility-assessment
cd <the folder containing reference_implementation.py>
python reference_implementation.py
```

On Linux or macOS, replace the first line with the minimal `conda create` command from section 2 — `environment.yml` is a Windows export and will not solve there.

On the classic Anaconda Prompt (cmd.exe) rather than PowerShell, use `cd /d "..."` instead of a plain `cd`: cmd.exe does not switch drives with a plain `cd`, so changing to a folder on another drive silently does nothing.

Expected output: a version check, the four hypervolumes, the three metrics, and a confirmation line:

```
Software versions (tested / installed):
  python     tested=3.11.15    installed=3.11.15    [ok]
  numpy      tested=2.4.6      installed=2.4.6      [ok]
  polytope   tested=0.2.5      installed=0.2.5      [ok]
  normaliz   tested=3.11.0     installed=3.11.0     [ok]

NS_600 vs. Scenario 4 (skill S1, cylindrical grinding)
  mu_R            = 4432331934923.373
  mu_Q            = 5454000.0
  mu_I            = 2430000.0
  mu_D            = 3024000.0
  Lambda          = 0.44554455445544555
  Gamma           = 5.482441377762043e-07
  Delta_E         = 6.82259371454832e-07

Values match the expected values.
```

**Assessing a different resource-scenario pair.** Edit the `NS_600_CR` (resource functionality subsets) and `SCENARIO4_QP` (requirement subsets) literals near the bottom of the script with the bounds and constraint coefficients of the pair you want, taken from Appendix Tables A1 and A3, and adjust `PROPERTIES_S1` if you are working in the S2 property space. `compute_metrics` is unchanged; it takes any such pair. Remove or update the `EXPECTED` assertion, which is specific to the worked example.

The procedure is **deterministic rather than sampled**. There is no sample size, sampling tolerance or convergence criterion to set: each polytope volume is computed by Normaliz in exact rational arithmetic on the supplied half-space representation, and the input preparation and the aggregation of the resulting volumes are done in double precision.

## 5. KNOWN LIMITATIONS/BUGS

**Scope.** The script implements Algorithm 1 and nothing else: Λ, Γ and Δ_E for one resource against one requirement scenario within one skill-specific property space.

**One worked example is carried.** Inputs are literals for NS_600 against Scenario 4 of skill S1. Any other resource-scenario pair has to be entered by hand from Appendix Tables A1 and A3.

**The counting-measure path is not exercised here.** A property space with no continuous properties is measured by counting admitted combinations rather than by volume. Both S1 and S2 have continuous properties, so this path is never taken in the case study and is not covered by the worked example's assertion.

**Cross-property constraints must be affine.** Each constraint becomes a single row of the half-space representation, so it has to be expressible as a weighted sum of properties against a constant, `c·x ≤ b`, with a `≥` relation written by negating the coefficients. The paper's own g1 and g2 are of this form. A nonlinear coupling — a product, ratio or power of two properties — cannot be represented.

**Cross-property constraints must be between continuous properties.** `subset_polytope` skips any constraint that names a discrete or categorical property, rather than folding that property's fixed value at the combination into the constant term. This is consistent with appendix Eq. (A1), which fixes the discrete and categorical values before forming the polytope, but it is a restriction on how such a dependency has to be written. Express it by splitting the resource or requirement into one subset per discrete or categorical combination, each with its own bounds and its own purely continuous constraints — which is how the case study models the external and internal subsets of each grinder.

**Normaliz must be locatable.** There is no bundled fallback. If `NORMALIZ_EXE`, the system `PATH` and the active conda prefix all miss it, every volume call fails.

**Normaliz failures are reported thinly.** The volume routine does not inspect the subprocess return code or its standard error, so any failure surfaces as `normaliz did not produce output for ...` regardless of the underlying cause.

**A 120 second ceiling per volume call.** Normaliz is invoked with a 120 second timeout. This is ample for the case study's seven- and nine-dimensional property spaces, but a substantially higher-dimensional constrained polytope can exceed it.

**Windows PATH/OpenMP hazard.** On Windows, other installed software (CAD/CAM tools, for example) can place its own numeric or OpenMP DLLs earlier on the system `PATH` than the conda environment's own DLLs. This can silently corrupt NumPy's internal state, and because Normaliz is OpenMP-based it then fails sporadically. Symptoms are a crash or hang with no output at all, or occasional `normaliz did not produce output` errors. A normally activated conda session usually prevents this, because activation puts the environment's own directories first on `PATH`. If it still happens, sanitise `PATH` by hand before running:

```powershell
$E = $env:CONDA_PREFIX
$env:PATH = "$E;$E\Library\mingw-w64\bin;$E\Library\usr\bin;$E\Library\bin;$E\Scripts;$E\bin;" +
            "C:\Windows\System32;C:\Windows;C:\Windows\System32\Wbem;" +
            "C:\Windows\System32\WindowsPowerShell\v1.0"
python reference_implementation.py
```

## 6. CHANGE LOG

- **v1.0 (2026)** — initial release, accompanying the reference publication.

## 7. LICENSE

    Copyright (c) 2026, FlandersMake@UGent - ISyE - Ghent University, Belgium

    This program is free software: you can redistribute it and/or modify
    it under the terms of the GNU General Public License as published by
    the Free Software Foundation, either version 3 of the License, or
    (at your option) any later version.

    This program is distributed in the hope that it will be useful,
    but WITHOUT ANY WARRANTY; without even the implied warranty of
    MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
    GNU General Public License for more details.

    You should have received a copy of the GNU General Public License
    along with this program.  If not, see <https://www.gnu.org/licenses/>.

The full licence text is in `LICENSE.txt`.

## 8. CONTACTING THE AUTHOR(S)

We would very much appreciate hearing from you if you use this software and find problems, or if you can think of ways it could be improved — and even (or is that 'especially'?) if you just think it's great. Even if the facility you would like to see appears to be of interest only to you, tell us about it; you'd be surprised how many ideas in that class have a much wider appeal.

See above for further contact information.

We read and consider all mail we receive, even though we may not have time to reply.

## 9. ACKNOWLEDGEMENTS

Reference publication:

> Van De Ginste, L., Aghezzaf, E.-H., & Cottyn, J. (2026). Hypervolume-based scalar metrics for skill-driven flexibility assessment in manufacturing system (re)design. *International Journal of Production Research* (UNDER REVIEW).

Developed at Flanders Make @ UGent (Belgium), within the Industrial Systems Engineering (ISyE) research group of Ghent University, in the context of various research projects on flexibility assessment.

This software relies on [Normaliz](https://www.normaliz.uni-osnabrueck.de/) for its exact volume computations, and on the [polytope](https://github.com/tulip-control/polytope) package for polyhedral set operations.
