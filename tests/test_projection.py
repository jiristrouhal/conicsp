import unittest
from math import pi, tan, sqrt

import conicsp
from conicsp import Point, rot

class Test_Point_Rotation(unittest.TestCase):

    def test_rotating_origin_has_no_effect(self):
        x = Point(0,0,0)
        self.assertEqual(rot(x, 0).xyz, (0,0,0))
        self.assertEqual(rot(x, pi/2).xyz, (0,0,0))
        self.assertEqual(rot(x, -pi/2).xyz, (0,0,0))

    def test_rotating_point_with_zero_y_and_z_has_no_effect(self):
        x = Point(1,0,0)
        self.assertEqual(rot(x, 0).xyz, (1,0,0))
        self.assertEqual(rot(x, pi/2).xyz, (1,0,0))
        self.assertEqual(rot(x, -pi/2).xyz, (1,0,0))

    def test_rotating_point_without_increasing_revolutions(self):
        self.assertEqual(rot(Point(0,0,1), pi/4), Point(0,-1/sqrt(2),1/sqrt(2)))
        self.assertEqual(rot(Point(0,1,0), pi/2), Point(0,0,1))
        self.assertEqual(rot(Point(0,1,1), pi/2), Point(0,-1,1))
        self.assertEqual(rot(Point(0,1,1), -pi), Point(0,-1,-1))
        self.assertEqual(rot(Point(0,1,0), -pi), Point(0,-1,0))

    def test_rotating_point_with_changing_revolutions(self):
        self.assertEqual(rot(Point(0,1,0), pi), Point(0,-1,0,1))
        self.assertEqual(rot(Point(0,1,0), 3*pi), Point(0,-1,0,2))
        self.assertEqual(rot(Point(0,1,0), -3*pi), Point(0,-1,0,-1))
        self.assertEqual(rot(Point(0,0,1), pi), Point(0,0,-1,1))
        self.assertEqual(rot(Point(0,0,-1), 3*pi), Point(0,0,1,1))
        self.assertEqual(rot(Point(0,0,-1), -pi), Point(0,0,1,-1))



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


    def test_point_on_positive_y_with_with_base_pi_is_projected_onto_negative_z(self):
        point = Point(0,2,0)
        image = self.proj.point(point, base=pi)
        self.assertAlmostEqual(image.x, point.x, 5)
        self.assertAlmostEqual(image.y, 0, 5)
        self.assertAlmostEqual(image.z, -2, 5)

    def test_point_on_positive_y_with_with_base_minus_pi_is_projected_onto_negative_z(self):
        point = Point(0,2,0)
        image = self.proj.point(point, base=-pi)
        self.assertAlmostEqual(image.x, point.x, 5)
        self.assertAlmostEqual(image.y, 0, 5)
        self.assertAlmostEqual(image.z, -point.y, 5)

    def test_point_on_negative_z_with_with_base_angle_pi_half_is_projected_onto_y_clockwise(self):
        point = Point(0,0,-2)
        i = self.proj.point(point, base=pi/2)
        self.assertAlmostEqual(i.y, -2, 5)
        self.assertAlmostEqual(i.z, 0, 5)

    def test_point_on_positive_z_with_with_base_minus_pi_half_is_projected_onto_positive_y(self):
        point = Point(0,0,2,0)
        i = self.proj.point(point, base=-pi/2)
        self.assertAlmostEqual(i.y, 2, 5)
        self.assertAlmostEqual(i.z, 0, 5)


if __name__=="__main__":
    unittest.main()