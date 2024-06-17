from __future__ import annotations
from math import sqrt, acos, sin, cos, pi, floor, atan
import dataclasses


Xyz = tuple[float,float,float]


TOLERANCE = 1e-07


@dataclasses.dataclass(frozen=True)
class Point:
    x: float
    y: float
    z: float
    rev: int = 0

    @property
    def xyz(self) -> tuple[float,float,float]:
        """Cartesian coordinates with respect to a base plane."""
        return self.x, self.y, self.z

    @property
    def r(self) -> float:
        """Distance from the origin"""
        return sqrt(self.x**2 + self.y**2 + self.z**2)

    @property
    def omega(self) -> float:
        """Angle between base plane and the plane defined by this point and the x axis."""
        if self.z==0:
            return 0 if self.y>=0 else pi
        else:
            return  sgn(self.z)*acos(self.y/sqrt(self.y**2+self.z**2))

    @property
    def phi(self) -> float:
        """Azimuth measured in the plane defined by this point and the x axis.

        Angle 0 is set to be aligned with the positive x axis.

        For points with nonzero y coordinate the pi/2 angle corresponds to the positive y direction,
        for zero y coordinate the point has always nonnegative azimuth.
        """
        phi = 0.0
        if self.y==0 and self.z==0:
            phi = 0 if self.x>=0 else -pi
        else:
            phi = sgn_0plus(self.y) * acos(self.x/self.r)
        return phi + 2*pi*self.rev


    def __eq__(self, other: object) -> bool:
        if not isinstance(other, Point):
            raise TypeError(f"Equality between point and {type(other)} is not defined.")
        dist_squared = (other.x-self.x)**2 + (other.y-self.y)**2 + (other.z-self.z)**2
        return dist_squared <= TOLERANCE**2 and self.rev==other.rev

    def __sub__(self, other: object) -> float:
        if not isinstance(other, Point):
            raise TypeError(f"Cannot subtract {type(other)} from Point.")
        dist_squared = (other.x-self.x)**2 + (other.y-self.y)**2 + (other.z-self.z)**2
        return sqrt(dist_squared)

    @staticmethod
    def from_xyz(*xyz: Xyz) -> list[Point]:
        """Get point from cartesian coordinates.

        Revolutions are set to zero.
        """
        return [Point(*coords) for coords in xyz]


def sgn(x: float) -> int:
    """Signum function."""
    return 1 if x>0 else (-1 if x<0 else 0)


def sgn_0plus(x: float) -> int:
    return 1 if x>=0 else -1


class Projector:

    def __init__(self, lambda_: float = 1.0) -> None:
        self._lambda = lambda_

    @property
    def lambda_(self) -> float:
        return self._lambda

    def point(self, point: Point, base: float = 0) -> Point:
        point1 = rot_yz(point, -base)
        return rot_yz(point1, base + point1.omega/self._lambda - point1.omega)

    def kappa(self, x: float, radius: float) -> float:
        if radius == 0:
            if x > 0:
                return 1
            else:
                raise ValueError("Cannot compute scaling of circle of radius zero behing singularity and non-unit lambda.")
        elif x==0:
            return self._lambda * sin(pi/(2*self._lambda))
        else:
            c = sqrt(radius**2 + x**2) / radius
            s = x / sqrt(radius**2 + x**2)
            return self._lambda*c*sin(acos(s)/self._lambda)


def rot_xy(point: Point, angle: float) -> Point:
    assert isinstance(point, Point)
    c, s = cos(angle), sin(angle)
    curr_angle = point.omega
    n, _ = rev(curr_angle + angle)
    return Point((c*point.x - s*point.y),  (s*point.x + c*point.y), point.z, rev=n)


def rot_yz(point: Point, angle: float) -> Point:
    c, s = cos(angle), sin(angle)
    return Point(point.x, (c*point.y - s*point.z), (s*point.y + c*point.z), rev=point.rev)


def rev(angle: float) -> tuple[int, float]:
    n = floor((angle+pi)/(2*pi))
    return n, angle-n*2*pi