import unittest
from math import pi

from conicsp import Point


class Test_Omega(unittest.TestCase):

    def test_omega_with_zero_revolutions(self):
        self.assertAlmostEqual(Point(0,0,0).omega, 0)
        self.assertAlmostEqual(Point(1,0,0).omega, 0)
        self.assertAlmostEqual(Point(-1,0,0).omega, 0)
        self.assertAlmostEqual(Point(0,1,0).omega, 0)
        self.assertAlmostEqual(Point(0,-1,0).omega, -pi)
        self.assertAlmostEqual(Point(0,0,1).omega, pi/2)
        self.assertAlmostEqual(Point(0,0,-1).omega, -pi/2)
        self.assertAlmostEqual(Point(0,1,1).omega, pi/4)
        self.assertAlmostEqual(Point(0,-1,1).omega, 3*pi/4)
        self.assertAlmostEqual(Point(0,1,-1).omega, -pi/4)
        self.assertAlmostEqual(Point(0,-1,-1).omega, -3*pi/4)
        self.assertAlmostEqual(Point(1,0,1).omega, pi/2)
        self.assertAlmostEqual(Point(1,0,-1).omega, -pi/2)
        self.assertAlmostEqual(Point(1,1,0).omega, 0)
        self.assertAlmostEqual(Point(1,-1,0).omega, -pi)

    def test_omega_with_revolutions_equal_to_one(self):
        self.assertAlmostEqual(Point(0,0,0,1).omega, 2*pi)
        self.assertAlmostEqual(Point(1,0,0,1).omega, 2*pi)
        self.assertAlmostEqual(Point(-1,0,0,1).omega, 2*pi)
        self.assertAlmostEqual(Point(0,1,0,1).omega, 2*pi)
        self.assertAlmostEqual(Point(0,-1,0,1).omega, pi)
        self.assertAlmostEqual(Point(0,0,1,1).omega, 5*pi/2)
        self.assertAlmostEqual(Point(0,0,-1,1).omega, 3*pi/2)
        self.assertAlmostEqual(Point(0,1,1,1).omega, 9*pi/4)
        self.assertAlmostEqual(Point(0,-1,1,1).omega, 11*pi/4)
        self.assertAlmostEqual(Point(0,1,-1,1).omega, 7*pi/4, 5)
        self.assertAlmostEqual(Point(0,-1,-1,1).omega, 5*pi/4)
        self.assertAlmostEqual(Point(1,0,1,1).omega, 5*pi/2)
        self.assertAlmostEqual(Point(1,0,-1,1).omega, 3*pi/2)
        self.assertAlmostEqual(Point(1,1,0,1).omega, 2*pi)
        self.assertAlmostEqual(Point(1,-1,0,1).omega, pi)

    def test_omega_with_revolutions_equal_to_minus_one(self):
        self.assertAlmostEqual(Point(0,0,0,-1).omega, -2*pi)
        self.assertAlmostEqual(Point(1,0,0,-1).omega, -2*pi)
        self.assertAlmostEqual(Point(-1,0,0,-1).omega, -2*pi)
        self.assertAlmostEqual(Point(0,1,0,-1).omega, -2*pi)
        self.assertAlmostEqual(Point(0,-1,0,-1).omega, -3*pi)
        self.assertAlmostEqual(Point(0,0,1,-1).omega, -3*pi/2)
        self.assertAlmostEqual(Point(0,0,-1,-1).omega, -5*pi/2)
        self.assertAlmostEqual(Point(0,1,1,-1).omega, -7*pi/4)
        self.assertAlmostEqual(Point(0,-1,1,-1).omega, -5*pi/4)
        self.assertAlmostEqual(Point(0,1,-1,-1).omega, -9*pi/4, 5)
        self.assertAlmostEqual(Point(0,-1,-1,-1).omega, -11*pi/4)
        self.assertAlmostEqual(Point(1,0,1,-1).omega, -3*pi/2)
        self.assertAlmostEqual(Point(1,0,-1,-1).omega, -5*pi/2)
        self.assertAlmostEqual(Point(1,1,0,-1).omega, -2*pi)
        self.assertAlmostEqual(Point(1,-1,0,-1).omega, -3*pi)


if __name__=="__main__":
    unittest.main()