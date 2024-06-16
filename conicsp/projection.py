from math import sqrt, acos, sin, pi


Point = tuple[float, float, float]


class Projector:

    def __init__(self, lambda_: float = 1.0) -> None:
        self._lambda = lambda_

    @property
    def lambda_(self) -> float:
        return self._lambda

    def point(self, point: Point, base: float = 0) -> list[Point]:
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

