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


class Test_Lambda_Two(unittest.TestCase):

    def setUp(self) -> None:
        self.proj = conicsp.Projector(lambda_=2.0)

    def test_points_on_x_axis_are_projected_on_themselves(self):
        points = Point.from_xyz((2,0,0), (1,0,0), (0,0,0), (-1,0,0))
        for p in points:
            with self.subTest(point = p):
                self.assertEqual(self.proj.point(p), p)

    def test_points_on_y_axis_with_base_plane_angle_set_to_0_are_projected_onto_themselves(self):
        points = Point.from_xyz((0,1,0), (0,5,0), (0,2,0))
        for p in points:
            with self.subTest(point = p):
                self.assertEqual(self.proj.point(p, base=0), p)

    def test_points_on_positive_z_axis_with_base_plane_angle_set_to_pi_half_are_projected_onto_themselves(self):
        points = Point.from_xyz((0,0,1), (0,0,5), (0,0,2))
        for p in points:
            with self.subTest(point = p):
                i = self.proj.point(p, base=pi/2)
                self.assertAlmostEqual(i.x, p.x, 5)
                self.assertAlmostEqual(i.y, p.y, 5)
                self.assertAlmostEqual(i.z, p.z, 5)

    def test_point_on_positive_y_with_with_base_minus_pi_is_projected_onto_negative_z(self):
        point = Point(0,2,0)
        image = self.proj.point(point, base=-pi)
        self.assertAlmostEqual(image.x, point.x, 5)
        self.assertAlmostEqual(image.y, 0, 5)
        self.assertAlmostEqual(image.z, -2, 5)

    def test_point_on_positive_y_with_with_base_pi_is_projected_onto_positive_z(self):
        point = Point(0,2,0)
        image = self.proj.point(point, base=-pi)
        self.assertAlmostEqual(image.x, point.x, 5)
        self.assertAlmostEqual(image.y, 0, 5)
        self.assertAlmostEqual(image.z, -2, 5)

    def test_point_on_negative_z_with_with_base_angle_pi_half_is_projected_onto_y_anticlockwise(self):
        point = Point(0,0,-2)
        i = self.proj.point(point, base=pi/2)
        self.assertAlmostEqual(i.y, 2, 5)
        self.assertAlmostEqual(i.z, 0, 5)

    def test_point_on_positive_z_with_with_base_minus_pi_half_is_projected_onto_positive_y(self):
        point = Point(0,0,2,0)
        i = self.proj.point(point, base=-pi/2)
        self.assertAlmostEqual(i.y, 2, 5)
        self.assertAlmostEqual(i.z, 0, 5)


if __name__=="__main__":
    unittest.main()