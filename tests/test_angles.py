import unittest
from math import pi

from conicsp import Point


class Test_Angles(unittest.TestCase):

    def test_zero_revolutions(self):
        self.assertAlmostEqual(Point(0,0,0).omega, 0)
        self.assertAlmostEqual(Point(0,0,0).phi, 0)

        self.assertAlmostEqual(Point(1,0,0).omega, 0)
        self.assertAlmostEqual(Point(1,0,0).phi, 0)

        self.assertAlmostEqual(Point(-1,0,0).omega, 0)
        self.assertAlmostEqual(Point(-1,0,0).phi, -pi)

        self.assertAlmostEqual(Point(0,1,0).omega, 0)
        self.assertAlmostEqual(Point(0,1,0).phi, pi/2)

        self.assertAlmostEqual(Point(0,-1,0).omega, 0)
        self.assertAlmostEqual(Point(0,-1,0).phi, -pi/2)

        self.assertAlmostEqual(Point(0,0,1).omega, pi/2)
        self.assertAlmostEqual(Point(0,0,1).phi, pi/2)

        self.assertAlmostEqual(Point(0,0,-1).omega, -pi/2)
        self.assertAlmostEqual(Point(0,0,-1).phi, pi/2)

        self.assertAlmostEqual(Point(0,1,1).omega, pi/4)
        self.assertAlmostEqual(Point(0,1,1).phi, pi/2)

        self.assertAlmostEqual(Point(0,-1,1).omega, pi/4)
        self.assertAlmostEqual(Point(0,-1,1).phi, -pi/2)

        self.assertAlmostEqual(Point(0,1,-1).omega, -pi/4)
        self.assertAlmostEqual(Point(0,1,-1).phi, pi/2)

        self.assertAlmostEqual(Point(0,-1,-1).omega, -pi/4)
        self.assertAlmostEqual(Point(0,-1,-1).phi, -pi/2)

        self.assertAlmostEqual(Point(1,0,1).omega, pi/2)
        self.assertAlmostEqual(Point(1,0,1).phi, pi/4)

        self.assertAlmostEqual(Point(1,0,-1).omega, -pi/2)
        self.assertAlmostEqual(Point(1,0,-1).phi, pi/4)

        self.assertAlmostEqual(Point(1,1,0).omega, 0)
        self.assertAlmostEqual(Point(1,1,0).phi, pi/4)

        self.assertAlmostEqual(Point(1,-1,0).omega, 0)
        self.assertAlmostEqual(Point(1,-1,0).phi, -pi/4)


class Test_Phi_With_Nonzero_Revolutions(unittest.TestCase):

    def test_phi_with_nonzero_revolutions(self):
        self.assertAlmostEqual(Point(0,1,0,-2).phi, -7*pi/2)
        self.assertAlmostEqual(Point(0,1,0,-1).phi, -3*pi/2)
        self.assertAlmostEqual(Point(0,1,0,1).phi, 5*pi/2)
        self.assertAlmostEqual(Point(0,1,0,2).phi, 9*pi/2)

        self.assertAlmostEqual(Point(1,1,0,-1).phi, -7*pi/4)
        self.assertAlmostEqual(Point(1,-1,0,-1).phi, -9*pi/4)
        self.assertAlmostEqual(Point(1,1,0,1).phi, 9*pi/4)
        self.assertAlmostEqual(Point(1,-1,0,1).phi, 7*pi/4)

        self.assertAlmostEqual(Point(0,1,1,1).phi, 5*pi/2)
        self.assertAlmostEqual(Point(0,1,-1,1).phi, 5*pi/2)
        self.assertAlmostEqual(Point(0,-1,1,1).phi, 3*pi/2)
        self.assertAlmostEqual(Point(0,-1,-1,1).phi, 3*pi/2)

        self.assertAlmostEqual(Point(0,0,1,1).phi, 5*pi/2)
        self.assertAlmostEqual(Point(0,0,-1,1).phi, 5*pi/2)
        self.assertAlmostEqual(Point(0,0,1,-1).phi, -3*pi/2)
        self.assertAlmostEqual(Point(0,0,-1,-1).phi, -3*pi/2)


if __name__=="__main__":
    unittest.main()