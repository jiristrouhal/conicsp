import unittest
from math import pi, sqrt

import conicsp
from conicsp import Point, rot_yz, rot_xy


def xyz_almost_equal(
    case: unittest.TestCase,
    first: tuple[float,float,float],
    second: tuple[float,float,float]
) -> None:

    for x, x_ in zip(first, second):
        case.assertAlmostEqual(x, x_)


class Test_Point_Rotation_In_YZ(unittest.TestCase):

    def test_rotating_origin_has_no_effect(self):
        x = Point(0,0,0)
        self.assertEqual(rot_yz(x, 0).xyz, (0,0,0))
        self.assertEqual(rot_yz(x, pi/2).xyz, (0,0,0))
        self.assertEqual(rot_yz(x, -pi/2).xyz, (0,0,0))

    def test_rotating_point_with_zero_y_and_z_has_no_effect(self):
        x = Point(1,0,0)
        xyz_almost_equal(self, rot_yz(x, 0).xyz, (1,0,0))
        xyz_almost_equal(self, rot_yz(x, pi/2).xyz, (1,0,0))
        xyz_almost_equal(self, rot_yz(x, -pi/2).xyz, (1,0,0))

    def test_rotating_point(self):
        self.assertEqual(rot_yz(Point(0,0,1), pi/4), Point(0,-1/sqrt(2),1/sqrt(2)))
        self.assertEqual(rot_yz(Point(0,1,0), pi/2), Point(0,0,1))
        self.assertEqual(rot_yz(Point(0,1,1), pi/2), Point(0,-1,1))
        self.assertEqual(rot_yz(Point(0,1,1), -pi), Point(0,-1,-1))
        self.assertEqual(rot_yz(Point(0,1,0), -pi), Point(0,-1,0))


class Test_Point_Rotation_In_XY(unittest.TestCase):

    def test_rotating_origin_has_no_effect(self):
        x = Point(0,0,0)
        self.assertEqual(rot_xy(x, 0).xyz, (0,0,0))
        self.assertEqual(rot_xy(x, pi/2).xyz, (0,0,0))
        self.assertEqual(rot_xy(x, -pi/2).xyz, (0,0,0))

    def test_rotating_point_with_zero_x_and_y_has_no_effect(self):
        p = Point(0,0,1)
        self.assertEqual(rot_xy(p, 0).xyz, (0,0,1))
        self.assertEqual(rot_xy(p, pi/2).xyz, (0,0,1))
        self.assertEqual(rot_xy(p, -pi/2).xyz, (0,0,1))

    def test_rotating_point_without_changing_revolutions(self):
        self.assertEqual(rot_xy(Point(0,1,0), pi/4), Point(-1/sqrt(2),1/sqrt(2),0))
        self.assertEqual(rot_xy(Point(1,0,0), pi/2), Point(0,1,0))
        self.assertEqual(rot_xy(Point(1,1,0), pi/2), Point(-1,1,0))
        self.assertEqual(rot_xy(Point(1,1,0), -pi), Point(-1,-1,0))
        self.assertEqual(rot_xy(Point(1,0,0), -pi), Point(-1,0,0))

    def test_rotating_point_with_changing_revolutions(self):
        self.assertEqual(rot_xy(Point(0,1,0), 2*pi), Point(0,1,0,1))
        self.assertEqual(rot_xy(Point(0,1,0), -2*pi), Point(0,1,0,-1))
        self.assertAlmostEqual(rot_xy(Point(0,1,0), 3*pi), Point(0,-1,0,1))
        self.assertAlmostEqual(rot_xy(Point(1,0,0), 3*pi), Point(-1,0,0,1))
        self.assertAlmostEqual(rot_xy(Point(1,0,0, 2), -5*pi), Point(-1,0,0,0))
        self.assertAlmostEqual(rot_xy(Point(1,0,0, 2), -5/2*pi), Point(0,-1,0,1))


class Test_Lambda_One(unittest.TestCase):

    def setUp(self):
        self.proj = conicsp.Projector(lambda_ = 1.0)

    def test_single_image_is_projected_to_the_original_point(self):
        coord_vals = {-1, 0 ,1}
        test_points = list()
        for x in coord_vals:
            for y in coord_vals:
                for z in coord_vals:
                    test_points.append(Point(x,y,z))
        for p in test_points:
            with self.subTest(point = p):
                i = self.proj.point(p)
                self.assertAlmostEqual(i.x, p.x, 5)
                self.assertAlmostEqual(i.y, p.y, 5)
                self.assertAlmostEqual(i.z, p.z, 5)


if __name__=="__main__":  # pragma: no cover
    unittest.main()