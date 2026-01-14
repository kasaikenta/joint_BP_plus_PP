# joint_BP_plus_PP

C++ implementation of a joint BP decoder with ETS/flip post-processing.

## Paper demo

This repository is a demo of the code construction and decoder proposed in:
https://arxiv.org/abs/2601.08824

Code construction (paper summary):
- Uses permutation matrices with controlled commutativity.
- Enforces orthogonality only where needed while keeping regular check-matrix structure.
- Enables large girth and avoids the usual distance upper bounds seen in more constrained designs.
- The included parameter file corresponds to the paper's example: a girth-8, (3,12)-regular
  [[9216,4612, <=48]] quantum LDPC code.

Decoder overview:
- Belief-propagation (BP) decoding with low-complexity post-processing.
- This implementation provides ETS-based and flip-based post-processing, and can be disabled
  with `--no-pp`.

This repository includes:
- `jointbp_ets.cpp` source
- Example code parameters and ETS files under `H_P768_J3_L12_dmax3_nc0-3_1-2_seed11579811919164041`

## Install and run (git)

```sh
git clone https://github.com/kasaikenta/joint_BP_plus_PP.git
cd joint_BP_plus_PP
c++ -O2 -std=c++17 -o jointbp_ets jointbp_ets.cpp
./jointbp_ets \
  --params H_P768_J3_L12_dmax3_nc0-3_1-2_seed11579811919164041.txt \
  --simulate --p 0.04000000 --trials 1000000 --max-iter 1000 \
  --flip-hist 5 --damping 0.000000 --report-every 20
```

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

Run command:

```sh
./jointbp_ets \
  --params H_P768_J3_L12_dmax3_nc0-3_1-2_seed11579811919164041.txt \
  --simulate --p 0.04000000 --trials 1000000 --max-iter 1000 \
  --flip-hist 5 --damping 0.000000 --report-every 20
```

Output excerpt (startup only; progress and summary omitted):

```text
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
P=768  J=3  L=12  L2=6  p=0.040000  n=9216  mx=2304  mz=2304
rank_x=2302  rank_z=2302  stab_rank=4604  k=4612  rate=0.500434
X: edges=27648  var_deg=3/3.00/3  check_deg=12/12.00/12
Z: edges=27648  var_deg=3/3.00/3  check_deg=12/12.00/12

+------------------------------------------------------------------------+
| Graph Checks                                                           |
+------------------------------------------------------------------------+
hx_hz_orthogonal=true  hx_hz_odd_pairs=0  girth_x=8  girth_z=8

+------------------------------------------------------------------------+
| ETS Load                                                               |
+------------------------------------------------------------------------+
ETS6:   X=48 Z=16 status=ok   
ETS8:   X= 0 Z=64 status=ok   
ETS8p4: X=48 Z= 0 status=ok   
ETS10:  X= 0 Z= 0 status=empty
ETS12:  X=23 Z= 0 status=ok   
ETS14:  X= 0 Z= 0 status=empty
ETS16:  X= 3 Z= 0 status=ok   
ETS18:  X=48 Z= 0 status=ok   
ETS20:  X= 0 Z= 7 status=ok   
ETS22:  X= 0 Z= 0 status=empty
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
