# Moment/SOS relaxations for power-system MINLPs (AC-OPF, AC-OTS, AC-UC)

This is research code for studying moment/sum-of-squares (Lasserre) relaxations of AC optimal transmission switching and AC unit commitment. It covers two threads:
1. Selecting relaxation strengthening: per-clique order, bound tightening, QC envelopes.
2. Correlated randomized rounding of binaries using the relaxation's pseudo-moments.

## Layout

| path | contents |
|---|---|
| `src/MomentMINLP.jl` | package module |
| `src/polynomials.jl`, `src/pop.jl` | sparse polynomials, POP container, Ipopt NLP solve of a POP |
| `src/sparsity.jl` | chordal extension / maximal cliques |
| `src/relaxation.jl` | sparse moment relaxation with binary reduction, per-clique orders, PMIs (MOSEK) |
| `src/power.jl` | AC-OPF / AC-OTS / AC-UC POPs from PowerModels data; configuration evaluation with PowerModels + Ipopt |
| `src/rounding.jl` | threshold, independent, Gaussian-copula, pseudo-moment conditional, and diving rounding |
| `scripts/instances.jl` | test instances (all data modifications are documented here) |
| `scripts/enumerate_instances.jl` | ground truth by enumerating all binary configurations |
| `scripts/experiment1.jl` | order 1 vs order 2 × rounding schemes |
| `scripts/analyze_experiment1.jl` | builds `results/experiment1_tables.md` from the JSON outputs |
| `scripts/experiment1b_cliques.jl` | clique augmentation and its effect on joint binary moments / rounding |
| `scripts/experiment2_scaling.jl` | case14/24/30 scaling: mixed-order and augmented sparse relaxations vs NLP relax-and-round (`VARIANTS=base` or `adjacent`) |
| `scripts/analyze_experiment2.jl` | builds `results/experiment2_tables.md` |
| `scripts/validate_formulations.jl` | checks the POP against PowerModels' AC-OPF and reports relaxation bounds |
| `data/` | MATPOWER test cases (from the PowerModels.jl test suite) |
| `docs/formulations.md` | mathematical formulations, as implemented |
| `docs/literature_notes.md` | literature review |
| `results/` | enumeration tables, experiment outputs, and write-ups |

## Setup

Requires Julia ≥ 1.10 and a MOSEK license.

```bash
julia --project=. -e 'using Pkg; Pkg.instantiate()'
julia --project=. scripts/check_mosek.jl
```

MOSEK looks for `~/mosek/mosek.lic`. The package also picks up `~/mosek/<version>/mosek.lic` automatically. Alternatively, set `MOSEKLM_LICENSE_FILE`.

## Typical workflow

```bash
julia --project=. scripts/validate_formulations.jl
julia --project=. scripts/enumerate_instances.jl       # writes results/enum_*.json
julia --project=. scripts/experiment1.jl case5_uc_sym  # writes results/experiment1_*.json
```

Significant findings, positive and negative, are summarized in `results/*.md` and in the shared Overleaf document.
