# WAVECAR Compressor

Small utility for reducing VASP `WAVECAR` size by lowering the stored plane-wave cutoff and rewriting the coefficient records. It can also analyze how much coefficient norm would be retained at different cutoff fractions before you choose a compression setting.

## Build

```bash
cmake -S . -B build
cmake --build build
```

The executable is written to:

```bash
build/wavecar_compressor
```

## Compress A WAVECAR

```bash
./build/wavecar_compressor compress <input WAVECAR> <output WAVECAR> <cutoff_fraction> [zero_small]
```

`cutoff_fraction` must be between `0` and `1`. Plane waves with kinetic energy above:

```text
cutoff_fraction * ENCUT
```

are removed. The output `WAVECAR` header is rewritten with the reduced `ENCUT`, new record length, and new per-k-point coefficient count.

Example:

```bash
./build/wavecar_compressor compress WAVECAR WAVECAR.small 0.75 true
```

The optional `zero_small` argument accepts `true/false`, `1/0`, `yes/no`, or `on/off`. When enabled, real or imaginary coefficient components with absolute value below `std::numeric_limits<double>::epsilon()` are written as exact zero.

## Analyze Cutoff Fractions

```bash
./build/wavecar_compressor analyze <input WAVECAR> [cutoff_fraction ...]
```

If no cutoff fractions are given, a default grid from `0` to `1` is used. The output is CSV with retained coefficient counts, retained `sum |C_G|^2`, worst per-band retained norm, estimated record length, and estimated output size.

Example:

```bash
./build/wavecar_compressor analyze WAVECAR 0.5 0.75 0.9 0.95 1.0
```

Key columns:

- `retained_norm_fraction`: total retained coefficient norm over all bands, spins, and k-points.
- `min_band_retained_norm_fraction`: worst retained norm for any single band.
- `estimated_size_ratio`: estimated output size divided by input size.

Analysis reads the coefficient payload once, so expect runtime comparable to streaming the whole `WAVECAR`.

## Python Wrapper

The `wavecar_compressor.py` wrapper shells out to the compiled executable.

```python
from wavecar_compressor import analyze, compress

rows = analyze("WAVECAR", [0.75, 0.9, 1.0])
compress("WAVECAR", "WAVECAR.small", 0.75, zero_small=True, create_bz2=True)
```

CLI form:

```bash
./wavecar_compressor.py analyze WAVECAR 0.75 0.9 1.0
./wavecar_compressor.py compress WAVECAR WAVECAR.small 0.75 --zero-small --bz2
```

With `--bz2`, the wrapper writes `WAVECAR.small.bz2` after the raw compressed `WAVECAR.small` is created. Use `--bz2-output <path>` to choose another archive path. Set `WAVECAR_COMPRESSOR_BIN` or pass `--binary` if the executable is not at `build/wavecar_compressor`.

## Notes And Caveats

The produced file is intended as a smaller VASP restart guess, not as an equivalent final wavefunction. Keep `NBANDS`, spin setup, k-points, and structure compatible with the original calculation.

The tool supports VASP coefficient tags `45200` and `45210`. Gamma-half storage is detected from the plane-wave count, but gamma-half ordering can be ambiguous; validate compressed gamma-only files with a small VASP restart before relying on them in production.
