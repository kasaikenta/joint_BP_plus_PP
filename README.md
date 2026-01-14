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

## Examples

Decode from syndrome files:

```sh
./jointbp_ets \
  --params H_P768_J3_L12_dmax3_nc0-3_1-2_seed11579811919164041.txt \
  --sx sx.txt --sz sz.txt
```

Decode from a known error:

```sh
./jointbp_ets \
  --params H_P768_J3_L12_dmax3_nc0-3_1-2_seed11579811919164041.txt \
  --err err.txt
```

## Example output (abridged)

Captured with:

```sh
./jointbp_ets \
  --params H_P768_J3_L12_dmax3_nc0-3_1-2_seed11579811919164041.txt \
  --simulate --p 0.05 --trials 2 --max-iter 10 --report-every 2
```

```text
JJJJJ  OOO  III  N   N TTTTT BBBB  PPPP      EEEE TTTTT SSSS
  J   O   O  I   NN  N   T   B   B P   P     E      T  S
  J   O   O  I   N N N   T   BBBB  PPPP      EEEE   T   SSS
J J   O   O  I   N  NN   T   B   B P         E      T     S
JJJ   OOO  III  N   N   T   BBBB  P         EEEE   T  SSSS

AI Report: JointBP ETS
Tips for getting started:
1. Use --help to see all options.
2. Use --simulate for randomized trials.
3. Use --report-every to print progress periodically.


+------------------------------------------------------------------------+
| Input                                                                  |
+------------------------------------------------------------------------+
params_path=H_P768_J3_L12_dmax3_nc0-3_1-2_seed11579811919164041.txt
fi=[(763,435) (679,69) (397,330) (61,18) (697,612) (373,246)]
gi=[(289,496) (257,640) (625,200) (41,524) (193,672) (449,672)]

+------------------------------------------------------------------------+
| Code Summary                                                           |
+------------------------------------------------------------------------+
P=768  J=3  L=12  L2=6  p=0.050000  n=9216  mx=2304  mz=2304
rank_x=2302  rank_z=2302  stab_rank=4604  k=4612  rate=0.500434
X: edges=27648  var_deg=3/3.00/3  check_deg=12/12.00/12
Z: edges=27648  var_deg=3/3.00/3  check_deg=12/12.00/12

...

+------------------------------------------------------------------------+
| PROGRESS                                                               |
+------------------------------------------------------------------------+
|trials      =        2 failures    =        2 bp_fail     =        2    |
|FER         = 1.000000 elapsed_s   =    0.02s latency     =   11.1ms    |
|fps         =    89.47 qbps        =     413k avg_iter    =    10.00    |
|iters       =       20 stab_success=        0 pp_success  =        0    |
|pp_rate     =   0.0000 pp_ets      =        0 pp_flip     =        0    |
|pp_osd      =        0 ets_used    =        0 ets_rate    =   0.0000    |
|ets6        =        0 ets8        =        0 ets8p4      =        0    |
|ets10       =        0 ets12       =        0 ets14       =        0    |
|ets16       =        0 ets18       =        0 ets20       =        0    |
|ets22       =        0                                                  |
+------------------------------------------------------------------------+

+------------------------------------------------------------------------+
| Run Summary                                                            |
+------------------------------------------------------------------------+
trials=2 failures=2 bp_fail=2 pp_success=0 flip_pp_success=0 stab_success=0 ets_used=0 ets_saved=0 ets6_used=0 ets6_saved=0 ets12_used=0 ets12_saved=0 fer=1.00000000 ci95=[0.34237195,1.00000000] iters=20 avg_iter=10.00 avg_latency=11.1ms exact_rate=0.000000 elapsed_s=0.02s
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
