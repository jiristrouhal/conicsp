import unittest
from math import pi, sqrt

import conicsp
from conicsp import Point


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
        point = Point(-2/sqrt(3),2,0)
        image = self.proj.point(point, base=-pi)
        self.assertAlmostEqual(image.x, point.x, 5)
        self.assertAlmostEqual(image.y, 0, 5)
        self.assertAlmostEqual(image.z, -2, 5)

    def test_point_on_positive_y_with_with_base_pi_is_projected_onto_positive_z(self):
        point = Point(-2/sqrt(3),2,0)
        image = self.proj.point(point, base=-pi)
        self.assertAlmostEqual(image.x, point.x, 5)
        self.assertAlmostEqual(image.y, 0, 5)
        self.assertAlmostEqual(image.z, -2, 5)

    def test_point_on_negative_z_with_with_base_angle_pi_half_is_projected_onto_y_anticlockwise(self):
        point = Point(-2/sqrt(3),0,-2)
        i = self.proj.point(point, base=pi/2)
        self.assertAlmostEqual(i.y, 2, 5)
        self.assertAlmostEqual(i.z, 0, 5)

    def test_point_on_positive_z_with_with_base_minus_pi_half_is_projected_onto_positive_y(self):
        point = Point(-2/sqrt(3),0,2,0)
        i = self.proj.point(point, base=-pi/2)
        self.assertAlmostEqual(i.y, 2, 5)
        self.assertAlmostEqual(i.z, 0, 5)


if __name__=="__main__":  # pragma: no cover
    unittest.main()