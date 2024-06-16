from __future__ import annotations
from math import sqrt, acos, sin, cos, pi, floor
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
        return self.x, self.y, self.z

    @property
    def omega(self) -> float:
        omega_0 = 0.0
        if self.z != 0:
            omega_0 = sgn(self.z)*acos(self.y/sqrt(self.y**2+self.z**2))
        else:
            omega_0 = -pi if self.y < 0 else 0
        return omega_0 + 2*pi*self.rev

    def __eq__(self, other: object) -> bool:
        if not isinstance(other, Point):
            raise TypeError
        dist_squared = (other.x-self.x)**2 + (other.y-self.y)**2 + (other.z-self.z)**2
        return dist_squared <= TOLERANCE**2 and self.rev==other.rev

    @staticmethod
    def from_xyz(*xyz: Xyz) -> list[Point]:
        return [Point(*coords) for coords in xyz]


def sgn(x: float) -> int:
    return 1 if x>0 else (-1 if x<0 else 0)


class Projector:

    def __init__(self, lambda_: float = 1.0) -> None:
        self._lambda = lambda_

    @property
    def lambda_(self) -> float:
        return self._lambda

    def point(self, point: Point, base: float = 0) -> list[Point]:
        point1 = rot(point, -base, rev_=False)
        omega = point1.omega
        point_1a = rot(point1, -omega)

        point1 = rot(point_1a, omega/self._lambda)
        point = rot(point1, base, rev_=False)
        return [point]

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

def rot(point: Point, angle: float, rev_: bool = True) -> Point:
    assert isinstance(point, Point)
    c, s = cos(angle), sin(angle)
    if rev_:
        curr_angle = point.omega
        n, _ = rev(curr_angle + angle)
    else:
        n = point.rev
    return Point(point.x, (c*point.y - s*point.z), (s*point.y + c*point.z), rev=n)


def rev(angle: float) -> tuple[int, float]:
    n = floor((angle+pi)/(2*pi))
    return n, angle-n*2*pi