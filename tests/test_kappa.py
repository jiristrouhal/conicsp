import unittest
from math import pi, tan

import conicsp


class Test_Kappa_For_Lambda_One(unittest.TestCase):

    def setUp(self) -> None:
        self.proj = conicsp.Projector(lambda_=1)

    def test_nonzero_radius_is_kappa_always_1(self):
        for r in (0.0001, 1, 1000):
            for x in (-1, 0, 1, 2, 3.5):
                self.assertAlmostEqual(self.proj.kappa(x, r), 1, 5)


class Test_Kappa_For_Lambda_Greater_Than_One(unittest.TestCase):

    def setUp(self) -> None:
        self.proj = conicsp.Projector(lambda_=3.6)

    def _test_kappa_for_nonzero_radius_is_greater_than_one_for_any_finite_x(self):
        for r in (0.0001, 1, 1000):
            for x in (-100, -1, -0.01, 0, 0.01, 1,  100):
                with self.subTest("", r=r, x=x):
                    self.assertGreater(self.proj.kappa(x, r), 1.00)

    def test_kappa_for_zero_radius_and_positive_x_is_equal_to_one(self):
        for x in (0.01, 1,  100):
            with self.subTest(x=x):
                self.assertEqual(self.proj.kappa(x, 0), 1.00)

    def test_kappa_for_zero_radius_and_nonpositive_x_raises_value_error(self):
        for x in (0, -0.01, -1,  -100):
            with self.subTest(x=x):
                with self.assertRaises(ValueError):
                    self.proj.kappa(x, 0)

    def test_kappa_for_zero_x_is_independent_of_radius_if_radius_nonzero(self):
        expected_val = self.proj.kappa(0, 1)
        for r in (0.0001, 2, 12, 1000):
            with self.subTest(x=0.0):
                self.assertAlmostEqual(self.proj.kappa(0, r), expected_val, 5)


class Test_Kappa_For_Lambda_Between_Half_and_One(unittest.TestCase):

    def setUp(self) -> None:
        self.proj = conicsp.Projector(lambda_=0.7)

    def test_kappa_for_nonzero_radius_is_smaller_than_one_for_any_nonegative_finite_x(self):
        for r in (0.01, 1, 1000):
            for x in (0, 0.01, 1,  100):
                with self.subTest("", r=r, x=x):
                    kappa = self.proj.kappa(x, r)
                    self.assertLess(kappa, 1.00)
                    self.assertGreater(kappa, 0.0)

    def test_kappa_for_zero_radius_and_positive_x_is_equal_to_one(self):
        for x in (0.01, 1,  100):
            with self.subTest(x=x):
                self.assertEqual(self.proj.kappa(x, 0), 1.00)

    def test_kappa_for_zero_radius_and_nonpositive_x_raises_value_error(self):
        for x in (0, -0.01, -1,  -100):
            with self.subTest(x=x):
                with self.assertRaises(ValueError):
                    self.proj.kappa(x, 0)

    def test_kappa_for_zero_x_is_independent_of_radius_if_radius_nonzero(self):
        expected_val = self.proj.kappa(0, 1)
        for r in (0.0001, 2, 12, 1000):
            with self.subTest(x=0.0):
                self.assertAlmostEqual(self.proj.kappa(0, r), expected_val, 5)



class Test_Kappa_For_Lambda_Less_Than_Half(unittest.TestCase):

    def setUp(self) -> None:
        self.proj = conicsp.Projector(lambda_=0.3)

    def test_kappa_for_nonzero_radius_is_smaller_than_one_for_any_finite_nonegative_x(self):
        for r in (0.01, 1, 1000):
            for x in (0, 0.01, 1,  100):
                with self.subTest(r=r, x=x):
                    kappa = self.proj.kappa(x, r)
                    self.assertLess(kappa, 1.00)
                    self.assertGreater(kappa, -1.00)

    def test_for_zero_radius_and_positive_x_is_equal_to_one(self):
        for x in (0.01, 1,  100):
            with self.subTest(x=x):
                self.assertEqual(self.proj.kappa(x, 0), 1.00)

    def test_for_zero_radius_and_nonpositive_x_raises_value_error(self):
        for x in (0, -0.01, -1,  -100):
            with self.subTest(x=x):
                with self.assertRaises(ValueError):
                    self.proj.kappa(x, 0)

    def test_for_zero_x_is_independent_of_radius_if_radius_nonzero(self):
        expected_val = self.proj.kappa(0, 1)
        for r in (0.0001, 2, 12, 1000):
            with self.subTest(x=0.0):
                self.assertAlmostEqual(self.proj.kappa(0, r), expected_val, 5)

    def test_for_some_maximum_positive_x_and_nonzero_radius_is_equal_to_zero(self):
        r = 0.2
        x = r/tan(self.proj.lambda_*pi)
        self.assertAlmostEqual(self.proj.kappa(x, r), 0.0, 5)


if __name__=="__main__":
    unittest.main()