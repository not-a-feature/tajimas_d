"""
tajimas-d: A small Tajima's-D module
See: https://github.com/not-a-feature/tajimas-d
Or:  https://pypi.org/project/tajimas-d/

@author: Jules Kreuer / not_a_feature
License: GPL-3.0
"""

from itertools import combinations
from typing import List
import argparse
import miniFasta
from sys import argv


def __check(sequences: List[str]) -> None:
    """
    Check for valid sequence proportions (1 < n, samle length).

    Input:
        sequences: list of str, list of sequences.
    Returns:
        Raises ValueError
    """
    if len(sequences) < 2:
        raise ValueError("At leat 2 sequences required")

    ref_len = len(sequences[0])
    if not all([len(s) == ref_len for s in sequences]):
        raise ValueError("Sequences are required to have the same length.")


def pi_estimator(sequences: List[str], safe=True) -> float:
    """
    Computes Θ_π, the Pi estimator.
    Θ_π = Number of pairwise differences / binomial(n, 2)

    Input:
        sequences: list of str, list of sequences.
        safe: bool, check if sequences have the same length.
    Returns:
        Θ_π: float: Pi estimator.
    """
    # Check for valid sequence proportions.
    if safe:
        __check(sequences)

    pairwise = combinations(sequences, 2)
    # Pairwise differences
    cs = [sum([not charA == charB for charA, charB in zip(seqA, seqB)]) for seqA, seqB in pairwise]
    n = len(sequences)
    binomial = ((n - 1) * n) / 2  # Binomial(n, 2)

    return sum(cs) / binomial


def __segregating(sequences: List[str]) -> int:
    """
    Counts the number of segregating sites.

    Input:
        sequences: list of str, list of sequences.
    Returns:
        seg_sites: int, number of segregating sites.
    """
    seg_sites = 0
    # For each position in sequence
    for i in range(0, len(sequences[0])):
        s = sequences[0][i]  # Take the char of first sequence as reference
        for seq in sequences:  # For each other sequence
            if not seq[i] == s:
                seg_sites += 1  # Add 1 if not equal
                break  # And stop comparing this position
    return seg_sites


def __harmonic(n: int) -> float:
    """
    Computes the n-1th harmonic number.
    Input:
        n: int
    Returns:
        h: float, n-1th harmonic number
    """
    return sum([1 / i for i in range(1, n)])


def watterson_estimator(sequences: List[str], safe=True) -> float:
    """
    Computes Θ_W, the Watterson estimator.
    Θ_W = Number of segregating sites / (n th harmonic number)
    Input:
        sequences: list of str, list of sequences.
        safe: bool, check for valid sequence proportions.
    Returns:
        Θ_W: float: Watterson estimator.
    """
    # Check for valid sequence proportions.
    if safe:
        __check(sequences)

    seg_sites = __segregating(sequences)
    harmonic = __harmonic(len(sequences))

    return seg_sites / harmonic


def tajimas_d(sequences: List[str]) -> float:
    """
    Computes Tajima's D for a list of sequences.
    See: https://www.ncbi.nlm.nih.gov/pmc/articles/PMC1203831/
    and https://ocw.mit.edu/courses/health-sciences-and-technology/hst-508-quantitative-
        genomics-fall-2005/study-materials/tajimad1.pdf

    Input:
        sequences: list of str, list of sequences.
    Returns:
        Θ_D: float: Tajima's D.
    """
    # Check for valid sequence proportions.
    __check(sequences)

    seg_sites = __segregating(sequences)  # Number of segregatig sites

    # Prevent devision by 0
    if seg_sites == 0:
        return 0

    theta_pi = pi_estimator(sequences, safe=False)  # Pi Estimator

    num_seq = len(sequences)  # Number of sequences
    harmonic = __harmonic(num_seq)  # N-1th harmonic number

    harmonic = sum([1 / i for i in range(1, num_seq)])  # Ref 3, aka. a1
    a2 = sum([1 / (i**2) for i in range(1, num_seq)])  # Ref 4

    b1 = (num_seq + 1) / (3 * (num_seq - 1))  # Ref 8
    b2 = (2 * (num_seq**2 + num_seq + 3)) / (9 * num_seq * (num_seq - 1))  # Ref 9

    c1 = b1 - 1 / harmonic
    c2 = b2 - ((num_seq + 2) / (harmonic * num_seq)) + (a2 / (harmonic**2))

    e1 = c1 / harmonic
    e2 = c2 / ((harmonic**2) + a2)

    delta_Theta = theta_pi - (seg_sites / harmonic)
    tD = delta_Theta / (((e1 * seg_sites) + (e2 * seg_sites * (seg_sites - 1))) ** 0.5)  # Ref 27
    return float(tD)


def parse_args(args):
    """Parse command line parameters.

    Args:
      args (List[str]): command line parameters as list of strings
          (for example  ``["--help"]``).

    Returns:
      :obj:`argparse.Namespace`: command line parameters namespace
    """
    # Parse arguments
    parser = argparse.ArgumentParser(
        description="""tajimas_d: Compute Tajima's D, the Pi- or Watterson-Estimator for multiple sequences."""
    )

    parser.add_argument(
        "-f",
        "--file",
        action="store",
        dest="path",
        help="Path to fasta file with all sequences.",
        required=True,
    )

    parser.add_argument(
        "-p",
        "--pi",
        action="store_true",
        dest="pi",
        help="Compute the Pi-Estimator score. (default)",
        required=False,
    )

    parser.add_argument(
        "-t",
        "--tajima",
        action="store_true",
        dest="tajima",
        help="Compute Tajima's D",
        required=False,
    )

    parser.add_argument(
        "-w",
        "--watterson",
        action="store_true",
        dest="watterson",
        help="Compute the Watterson-Estimator score.",
        required=False,
    )

    return parser.parse_args()


def _main_cli(args):
    args = parse_args(args)
    # Load sequences
    sequences = [mf.body for mf in miniFasta.read(args.path)]

    # Compute
    if args.tajima or not (args.pi or args.watterson):
        print(f"Tajima's D score:\t\t{tajimas_d(sequences)}")

    if args.pi:
        print(f"Pi-Estimator score:\t\t{pi_estimator(sequences)}")

    if args.watterson:
        print(f"Watterson-Estimator score:\t{watterson_estimator(sequences)}")


def run_cli():
    """Calls :func:`_main_cli` passing the CLI arguments extracted from :obj:`sys.argv`.

    This function is used as entry point to create a console script with setuptools.
    """
    _main_cli(argv[1:])


if __name__ == "__main__":
    # ^  This is a guard statement that will prevent the following code from
    #    being executed in the case someone imports this file instead of
    #    executing it as a script.
    #    https://docs.python.org/3/library/__main__.html

    # After installing your project with pip, users can also run your Python
    # modules as scripts via the ``-m`` flag, as defined in PEP 338::
    #
    #     python -m tajimas_d.tajimas_d [ARGS]
    #
    run_cli()
