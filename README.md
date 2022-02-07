![tajimas-d](https://github.com/not-a-feature/tajimas_d/raw/main/tajimas-d.png)

Compute the Tajima's-D, Pi-Estimator or Watterson-Estimator for multiple sequences.

![Test Badge](https://github.com/not-a-feature/tajimas_d/actions/workflows/tests.yml/badge.svg)
![Python Version Badge](https://img.shields.io/pypi/pyversions/tajimas_d)
![Download Badge](https://img.shields.io/pypi/dm/tajimas_d.svg)

Tajima's D is a population genetic test statistic created which computes the difference between the mean number of pairwise differences and the number of segregating sites. It is used to determine if a population is expanding or shrinking.

## Tajima's D
Tajima's D is defined as follows:

![Tajima](https://render.githubusercontent.com/render/math?math=\theta_\text{Tajima}=\frac{\theta_{\pi}%20-%20\theta_{W}}{\sqrt{\text{Var}(\theta_{\pi}-\theta_{W})}})

## Pi-Estimator
The π estimator is the average number of pairwise differences between any two sequences:

![Pi](https://render.githubusercontent.com/render/math?math=\theta_{\pi}=\frac{\text{Nr.%20of%20pairwise%20differences}}{\binom{n}{2}})

## Watterson-Estimator
The Watterson estimator is the expected number of segregating sites.

![Watterson](https://render.githubusercontent.com/render/math?math=\theta_{\W}=\frac{\text{Nr.%20of%20segregating%20sites}}{\sum^{n-1}_{i=1}\frac{1}{i}})

## Installation
Using pip  / pip3:
```bash
pip install tajimas_d
```
Or by source:
```bash
git clone git@github.com:not-a-feature/tajimas_d.git
cd tajimas_d
pip install .
```

## How to use

```python
from tajimas_d import tajimas_d, watterson_estimator, pi_estimator

sequences = ["AAAA", "AAAT", "AAGT", "AAGT"]

theta_tajima = tajimas_d(sequences)
theta_pi = pi_estimator(sequences)
theta_w = watterson_estimator(sequences)
```


## Standalone version

Clone the repository and cd into it. Intall the requirements (miniFasta > 2.2, for fasta reading).

```
usage: tajimas-d.py [-h] -f PATH [-p] [-t] [-w]

tajimas-d: Compute the Tajima s-D, Pi-Estimator or Watterson-Estimator for multiple
sequences.

optional arguments:
  -h, --help            show this help message and exit
  -f PATH, --file PATH  Path to fasta file with all sequences.
  -p, --pi              Compute the Pi-Estimator score. (default)
  -t, --tajima          Compute the Pi-Estimator score.
  -w, --watterson       Compute the Watterson-Estimator score.

```

## License
```
Copyright (C) 2021 by Jules Kreuer - @not_a_feature
This piece of software is published unter the GNU General Public License v3.0
TLDR:

| Permissions      | Conditions                   | Limitations |
| ---------------- | ---------------------------- | ----------- |
| ✓ Commercial use | Disclose source              | ✕ Liability |
| ✓ Distribution   | License and copyright notice | ✕ Warranty  |
| ✓ Modification   | Same license                 |             |
| ✓ Patent use     | State changes                |             |
| ✓ Private use    |                              |             |
```
Go to [LICENSE.md](https://github.com/not-a-feature/tajimas_d/blob/main/LICENSE) to see the full version.

## Dependencies
In addition to packages included in Python 3, this piece of software uses 3rd-party software packages for development purposes that are not required in the published version.
Go to [DEPENDENCIES.md](https://github.com/not-a-feature/tajimas_d/blob/main/DEPENDENCIES.md) to see all dependencies and licenses.