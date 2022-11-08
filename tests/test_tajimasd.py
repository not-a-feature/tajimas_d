# Testing Tajima's D Handling
from tajimas_d import tajimas_d, watterson_estimator, pi_estimator
import pytest
from os import path
import miniFasta


def getSeqs(path):
    return [mf.body for mf in miniFasta.read(path)]


allA = getSeqs(path.join(path.dirname(__file__), "allA.fasta"))

# Example from https://ocw.mit.edu/courses/health-sciences-and-technology/hst-508-quantitative-
# genomics-fall-2005/study-materials/tajimad1.pdf
MIT_example = getSeqs(path.join(path.dirname(__file__), "MIT_example.fasta"))


@pytest.mark.parametrize("sequences,expected", [(allA, 0), (MIT_example, -1.446172)])
def test_tajima(sequences, expected):
    assert round(tajimas_d(sequences), 6) == expected


@pytest.mark.parametrize("sequences,expected", [(allA, 0), (MIT_example, 3.888889)])
def test_pi(sequences, expected):
    assert round(pi_estimator(sequences), 6) == expected


@pytest.mark.parametrize(
    "sequences,expected", [(allA, 0), (MIT_example, 5.655772)]
)  # Check if thats true
def test_watterson(sequences, expected):
    assert round(watterson_estimator(sequences), 6) == expected


@pytest.mark.xfail
def test_tajima_diff_len():
    sequences = ["AAAA", "AAAB", "AAA"]
    tajimas_d(sequences)


@pytest.mark.xfail
def test_watterson_diff_len():
    sequences = ["AAAA", "AAAB", "AAA"]
    watterson_estimator(sequences)


@pytest.mark.xfail
def test_pi_diff_len():
    sequences = ["AAAA", "AAAB", "AAA"]
    pi_estimator(sequences)
