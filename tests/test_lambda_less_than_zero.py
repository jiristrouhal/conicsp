import unittest
from math import pi, sqrt

import conicsp
from conicsp import Point


def xyz_almost_equal(
    case: unittest.TestCase,
    first: tuple[float,float,float],
    second: tuple[float,float,float]
) -> None:

    for x, x_ in zip(first, second):
        case.assertAlmostEqual(x, x_)


class Test_Lambda_Two(unittest.TestCase):

    def setUp(self) -> None:
        self.proj = conicsp.Projector(lambda_=0.5)

    def test_origin_is_projected_on_itself(self):
        point = Point(0,0,0)
        self.assertEqual(self.proj.point(point), point)

    def test_point_on_x_axis_is_projected_on_itself(self):
        point = Point(5,0,0,2)
        self.assertEqual(self.proj.point(point), point)

    def test_point_on_positive_half_of_base_plane_is_projected_on_itself(self):
        point = Point(0,2,0)
        self.assertEqual(self.proj.point(point, base=0), point)
        point = Point(0,0,1)
        self.assertEqual(self.proj.point(point, base=pi/2), point)
        point = Point(0,-1,0)
        self.assertEqual(self.proj.point(point, base=pi), point)
        point = Point(0,0,-1)
        self.assertEqual(self.proj.point(point, base=-pi/2), point)

    def test_point_on_negative_half_of_base_plane_is_projected_onto_the_positive(self):
        point = Point(2/sqrt(3),-2,0)
        xyz_almost_equal(self,self.proj.point(point, base=0).xyz, (2/sqrt(3),2,0))
        point = Point(2/sqrt(3),-2,0)
        xyz_almost_equal(self,self.proj.point(point, base=0).xyz, (2/sqrt(3),2,0))

        point = Point(3/sqrt(3),3,0)
        xyz_almost_equal(self,self.proj.point(point, base=pi).xyz, (3/sqrt(3),-3,0))
        point = Point(7/sqrt(3),0,7)
        xyz_almost_equal(self,self.proj.point(point, base=-pi/2).xyz, (7/sqrt(3),0,-7))
        point = Point(7/sqrt(3),0,-7)
        xyz_almost_equal(self,self.proj.point(point, base=pi/2).xyz, (7/sqrt(3),0,7))


if __name__=="__main__":  # pragma: no cover
    unittest.main()