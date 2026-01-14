# joint_BP_plus_PP

C++ implementation of a joint BP decoder with ETS/flip post-processing.

This repository includes:
- `jointbp_ets.cpp` source
- `jointbp_ets` prebuilt binary (macOS arm64)
- Example code parameters and ETS files under `H_P768_J3_L12_dmax3_nc0-3_1-2_seed11579811919164041`

## Build

```sh
c++ -O2 -std=c++17 -o jointbp_ets jointbp_ets.cpp
```

## Quick start

```sh
./jointbp_ets \
  --params H_P768_J3_L12_dmax3_nc0-3_1-2_seed11579811919164041.txt \
  --simulate --p 0.05 --trials 1000 --max-iter 50
```

## Input modes

- `--simulate`: random trials.
- `--sx FILE --sz FILE`: decode a single syndrome.
- `--err FILE`: decode a single known error (mutually exclusive with `--simulate` and `--sx/--sz`).

Input file formats:
- `--sx/--sz`: one `0`/`1` per line.
- `--err`: one `0`/`1`/`2`/`3` per line (`0=I`, `1=X`, `2=Z`, `3=Y`).

## ETS files

By default, ETS files are loaded from a directory named after the params file basename.
For example, `H_P...txt` uses:

- `H_P.../ets_6_2.txt`
- `H_P.../ets_8_2.txt`
- `H_P.../ets_8_2_path4.txt`
- `H_P.../ets_10_2.txt` ... `ets_22_2.txt`

Override with `--ets6-x`, `--ets6-z`, `--ets12-x`, `--ets12-z`, or `--ets8-path4`.

## Help

```sh
./jointbp_ets --help
```

## Notes

- `--no-pp` disables post-processing.
- `--verbose` and `--verbose-all` control per-iteration output.
