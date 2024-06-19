import unittest
from math import pi, sqrt

import conicsp
from conicsp import Point, limit_azimuth


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


class Test_Limiting_Azimuth(unittest.TestCase):

    def test_limiting_azimuth_between_minus_pi_and_pi_has_no_effect_for_any_positive_gamma(self):
        # negativa gamma is not allowed
        for gamma in (pi, 2*pi, 3*pi, 8*pi):
            for az in (0, pi/4, -pi/4, pi/2, -pi/2, 0.75*pi, -0.75*pi, pi, -pi):
                self.assertEqual(limit_azimuth(az, gamma), az)

    def test_azimuth_falling_between_gamma_minus_and_plus_pi_is_unchanged(self):
        gamma = 5*pi
        for az in (5.1*pi, 5*pi, 4.5*pi, 4*pi):
            with self.subTest(az=f"{az/pi} pi", gamma=f"{gamma/pi} pi"):
                self.assertEqual(limit_azimuth(az, gamma), az)

    def test_azimuth_for_gamma_greater_than_one_that_is_not_one_pi_from_the_mulitple_of_gamma_is_set_to_sign_of_azimuth_times_pi(self):
        gamma = 4*pi
        for az in (1.1*pi, 2*pi, 6*pi, 10*pi):
            self.assertEqual(limit_azimuth(az, gamma), pi)
        for az in (-1.1*pi, -2*pi, -6*pi, -10*pi):
            self.assertEqual(limit_azimuth(az, gamma), -pi)

    def test_azimuth_equal_to_multiple_of_gamma_minus_or_plus_pi_is_unchanged(self):
        gamma = 4*pi
        for az in (3*pi, 7*pi, pi, 5*pi):
            with self.subTest(az=f"{az/pi} pi", gamma=f"{gamma/pi} pi"):
                self.assertEqual(limit_azimuth(az, gamma), az)
        for az in (-3*pi, -7*pi, -pi, -5*pi):
            with self.subTest(az=f"{az/pi} pi", gamma=f"{gamma/pi} pi"):
                self.assertEqual(limit_azimuth(az, gamma), az)


# class Test_Limiting_Azimuth_In_Projecting_Point(unittest.TestCase):

#     def setUp(self) -> None:
#         self.proj = conicsp.Projector(lambda_=2.0)

#     def test_point_falling_into_unobservable_section_of_the_space_is_projected_behind_origin(self):
#         point = Point(-sqrt(2)/2, sqrt(2)/2, 0, 1)
#         image = self.proj.point(point, base=pi)
#         self.assertAlmostEqual(image.y, 0)


if __name__=="__main__":  # pragma: no cover
    # runner = unittest.TextTestRunner()
    # runner.run(Test_Limiting_Azimuth("test_point_falling_into_unobservable_section_of_the_space_is_projected_behind_origin"))
    unittest.main()