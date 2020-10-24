"""
Acid equation:
x[Ha Xb Oc]^m + yH+ + ne -> wXd + zH2O
[a, b, c, d, m] -> (x, y, n, w, z)

Basic equation:
x[Ha Xb Oc]^m + rH2O + ne -> wXd + zH2O + sOH-
[a, b, c, d, m] -> (x, r, n, w, z, s)
"""

from helpers.core import coef_ac, coef_bas


def test_dichromate_to_chromium_acid():
    """ Cr2O7^{2-} -> Cr """
    assert coef_ac([0, 2, 7, 1, -2]) == (1, 14, 12, 2, 7)


def test_chromate_to_chromium_basic():
    """ CrO4^{2-} -> Cr """
    assert coef_bas([0, 1, 4, 1, -2]) == (1, 8, 6, 1, 4, 8)


def test_BiH3_to_Bi_acid():
    """ BiH3 -> Bi """
    assert coef_ac([3, 1, 0, 1, 0]) == (1, -3, -3, 1, 0)


def test_BiH3_to_Bi_basic():
    """ BiH3 -> Bi """
    assert coef_bas([3, 1, 0, 1, 0]) == (1, -3, -3, 1, 0, -3)


def test_water_acid():
    """ H2O -> O2 """
    assert coef_ac([2, 1, 0, 2, 0]) == (1, -2, -2, 1/2, 0)


def test_water_basic():
    """ OH- -> O2 """
    assert coef_bas([1, 1, 0, 2, -1]) == (1, -1, -2, 1/2, 0, -1)
