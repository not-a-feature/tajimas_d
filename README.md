<img src="https://github.com/not-a-feature/tajimas_d/raw/main/tajimas-d.png" width=300px alt="tajimas_d logo"></img>

Compute the Tajima's-D, Pi-Estimator or Watterson-Estimator for multiple sequences.

This module is now part of the bfx suite. See https://py-bfx.readthedocs.io for more information.

![Test Badge](https://github.com/not-a-feature/tajimas_d/actions/workflows/tests.yml/badge.svg)
![Python Version Badge](https://img.shields.io/pypi/pyversions/tajimas_d)
![Download Badge](https://img.shields.io/pypi/dm/tajimas_d.svg)
[![Code style: black](https://img.shields.io/badge/code%20style-black-000000.svg)](https://github.com/psf/black)

Tajima's D is a population genetic test statistic that computes the difference between the mean number of pairwise differences and the number of segregating sites. It is used to determine whether a population is expanding or shrinking.

## Tajima's D
Tajima's D is defined as follows:
$\theta_\text{Tajima}=\frac{\theta_{\pi}%20-%20\theta_{W}}{\sqrt{\text{Var}(\theta_{\pi}-\theta_{W})}}$

If $\theta_\text{Tajima}<0$, there are many rare variants, indicating an **expanding** population.

Whereas $0<\theta_\text{Tajima}$, indicates an **declining** population as there are many intermediate variants.

A result is consideres significant if $\theta_\text{Tajima}<-2$ or $2<\theta_\text{Tajima}$.

## Pi-Estimator
The π estimator is the average number of pairwise differences between any two sequences:

$\theta_{\pi}=\frac{\text{Nr. of pairwise differences}}{\binom{n}{2}}$

## Watterson-Estimator
The Watterson estimator is the expected number of segregating sites.

$\theta_{W}=\frac{\text{Nr. of segregating sites}}{\Sigma_{i=1}^{n-1}\frac{1}{i}}$

## Installation
Using pip  / pip3:
```bash
pip install tajimas_d
```

Using conda:
```bash
conda install -c bioconda tajimas_d
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
The standalone version requires `miniFasta>=2.2` to be installed.

```
usage: tajimas_d [-h] -f PATH [-p] [-t] [-w]

tajimas_d: Compute Tajima's D, the Pi- or Watterson-Estimator for multiple
sequences.

optional arguments:
  -h, --help            show this help message and exit
  -f PATH, --file PATH  Path to fasta file with all sequences.
  -p, --pi              Compute the Pi-Estimator score.
  -t, --tajima          Compute the Pi-Estimator score. (default)
  -w, --watterson       Compute the Watterson-Estimator score.

```

## License

Copyright (C) 2024 by Jules Kreuer - @not_a_feature

This piece of software is published unter the GNU General Public License v3.0
TLDR:

| Permissions      | Conditions                   | Limitations |
| ---------------- | ---------------------------- | ----------- |
| ✓ Commercial use | Disclose source              | ✕ Liability |
| ✓ Distribution   | License and copyright notice | ✕ Warranty  |
| ✓ Modification   | Same license                 |             |
| ✓ Patent use     | State changes                |             |
| ✓ Private use    |                              |             |

Go to [LICENSE.md](https://github.com/not-a-feature/tajimas_d/blob/main/LICENSE) to see the full version.
