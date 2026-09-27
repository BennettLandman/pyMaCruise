# TGDM synthetic runtime check

The reproducible benchmark is `scripts/benchmark_tgdm.py`. It runs the full
default MaCRUISE TGDM stage sequence on concentric synthetic tissue fields,
with an ACE skeleton and a zero-valued GVF field. It uses no subject images,
does not write image outputs, and reports the host process's peak resident set
size (including Python and imported libraries).

Run it from the repository root:

```sh
.venv/bin/python scripts/benchmark_tgdm.py --size 25
.venv/bin/python scripts/benchmark_tgdm.py --size 33
.venv/bin/python scripts/benchmark_tgdm.py --size 49
```

Measurements on 2026-09-27 with Python 3.13.2 on an arm64 Mac:

| Cubic volume | Voxels | Full sequence time | Process peak RSS | Finite outputs |
| --- | ---: | ---: | ---: | --- |
| 25³ | 15,625 | 11.0 s | 82.2 MiB | Yes |
| 33³ | 35,937 | 20.3 s | 101.8 MiB | Yes |
| 49³ | 117,649 | 48.5 s | 145.5 MiB | Yes |

Each run completes the selected wrapper's full 790 PDE updates: initial
smoothing 180, initial inner 100, central 190, central shrink 10, final inner
100, and outer 210. The separately reported completed PDE steps are available
in the script's JSON output.

These synthetic checks show the current NumPy implementation runs the full
schedule, produces finite nested outputs, and provide a repeatable runtime /
memory baseline. They do not include GVF calculation, NIfTI serialization, or
mesh extraction, and do not support extrapolating runtime or peak memory to
subject-sized volumes. Stage-wise numerical parity against the Java reference
is reported separately in [tgdm-reference-comparison.md](tgdm-reference-comparison.md).
