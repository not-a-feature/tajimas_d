"""
tajimas-d: A small Tajima's-D module
This file is only used for the standalone version.

See: https://github.com/not-a-feature/tajimas_d
Or:  https://pypi.org/project/tajimas_d/

@author: Jules Kreuer / not_a_feature
License: GPL-3.0
"""

import argparse
import miniFasta
from src.tajimas_d._tajimas_d import tajimas_d, watterson_estimator, pi_estimator


if __name__ == "__main__":
    # Parse arguments
    parser = argparse.ArgumentParser(
        description="""
tajimas-d: Compute the Tajima's-D, Pi-Estimator or Watterson-Estimator for multiple sequences.
"""
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
        help="Compute the Pi-Estimator score.",
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

    args = parser.parse_args()

    # Load sequences
    sequences = [mf.body for mf in miniFasta.read(args.path)]

    # Compute
    if args.tajima or not (args.pi or args.watterson):
        print(f"Tajima's D score:\t\t{tajimas_d(sequences)}")

    if args.pi:
        print(f"Pi-Estimator score:\t\t{pi_estimator(sequences)}")

    if args.watterson:
        print(f"Watterson-Estimator score:\t{watterson_estimator(sequences)}")
