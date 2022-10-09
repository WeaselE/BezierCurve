from __future__ import annotations

from typing import List, Tuple, Optional
import sys

from matplotlib import pyplot as plt
from numpy import linspace


class PointError(Exception):
    """
    PointError

    Args:
        Exception (_type_): _description_
    """

    def __repr__(self) -> str:
        exc_type, exc_obj, exc_tb = sys.exc_info()
        return f"Invalid point given.\n{str(exc_type)}\n{str(exc_obj)}"


class BezierCurve:
    """
    Bezier Curve class
    """

    def __init__(
        self,
        point1: Optional[Tuple[float, float]] = (),
        point2: Optional[Tuple[float, float]] = (),
        point3: Optional[Tuple[float, float]] = (),
        point4: Optional[Tuple[float, float]] = (),
        t_point: float = 0.5,
    ) -> None:

        self.point1 = point1
        self.point2 = point2
        self.point3 = point3
        self.point4 = point4
        if not self.point4:
            self.type = "Quadratic"
        elif not self.point3:
            self.type = "Linear"
        else:
            self.type = "Cubic"
        self.curve_points = []
        self.t_point = t_point

    def lerp(
        self,
        l_point1: Tuple[float, float],
        l_point2: Tuple[float, float],
        t_point: float,
    ) -> Tuple[float, float]:

        """
        Calculates Linear Interpolation (lerp) from given 2 points and t position.
        t must be in range 0-1.
        """

        if self.type == "linear" or not self.type:
            l_point1, l_point2 = self._parse(
                t_point=t_point, point1=l_point1, point2=l_point2
            )
        return (1 - self.t_point) * l_point1[0] + self.t_point * l_point2[0], (
            1 - self.t_point
        ) * l_point1[1] + t_point * l_point2[1]

    def quadratic_lerp(
        self,
        q_point1: Tuple[float, float],
        q_point2: Tuple[float, float],
        q_point3: Tuple[float, float],
        t_point: float,
    ) -> Tuple[float, float]:

        """
        Calculates Linear Interpolation (lerp) from given 3 points and t position.
        t must be in range 0-1.
        """

        if self.type == "Quadratic":
            q_point1, q_point2, q_point3 = self._parse(
                t_point=t_point, point1=q_point1, point2=q_point2, point3=q_point3
            )
        return self.lerp(
            l_point1=self.lerp(
                l_point1=q_point1, l_point2=q_point2, t_point=self.t_point
            ),
            l_point2=self.lerp(
                l_point1=q_point2, l_point2=q_point3, t_point=self.t_point
            ),
            t_point=self.t_point,
        )

    def cubic_lerp(
        self,
        c_point1: Tuple[float, float],
        c_point2: Tuple[float, float],
        c_point3: Tuple[float, float],
        c_point4: Tuple[float, float],
        t_point: float,
    ) -> Tuple[float, float]:

        """
        Calculates Linear Interpolation (lerp) from 4 points and t position.
        t must be value between 0 and 1.
        """
        if self.type == "Cubic":
            self._parse(
                t_point=t_point,
                point1=c_point1,
                point2=c_point2,
                point3=c_point3,
                point4=c_point4,
            )
        return self.lerp(
            l_point1=self.quadratic_lerp(
                q_point1=self.point1,
                q_point2=self.point2,
                q_point3=self.point3,
                t_point=self.t_point,
            ),
            l_point2=self.quadratic_lerp(
                q_point1=self.point2,
                q_point2=self.point3,
                q_point3=self.point4,
                t_point=self.t_point,
            ),
            t_point=self.t_point,
        )

    def calculate_points(self, points: int = 20) -> List[Tuple[float, float]]:
        """
        Calculates Bezier Curve points from given position points.

        _summary_

        Args:
            points (int, optional): _description_. Defaults to 20.

        Returns:
            List[Tuple[float, float]]: _description_
        """
        t_range = linspace(0, 1, points)
        self.curve_points = []
        if self.type == "Cubic":
            for t in t_range:
                self.curve_points.append(
                    self.cubic_lerp(
                        c_point1=self.point1,
                        c_point2=self.point2,
                        c_point3=self.point3,
                        c_point4=self.point4,
                        t_point=t,
                    )
                )
            return self.curve_points

        for t in t_range:
            self.curve_points.append(
                self.quadratic_lerp(
                    q_point1=self.point1,
                    q_point2=self.point2,
                    q_point3=self.point3,
                    t_point=t,
                )
            )
        return self.curve_points

    def draw_curve(self) -> None:
        """
        Draws Bezier Curve from points calculated.
        """
        if self.curve_points is None:
            raise AttributeError("Curve points have not been calculated.")
        plt.scatter(*zip(*self.curve_points))
        plt.show()

    def _parse(
        self,
        t_point: float,
        point1: Tuple[float, float],
        point2: Tuple[float, float],
        point3: Optional[Tuple[float, float]] = None,
        point4: Optional[Tuple[float, float]] = None,
    ) -> Tuple:

        # TODO: Use protective statements rather than checking for 4th and 3rd point first.
        # Should be more concise code.
        """
        _summary_

        """

        # Check that t is type float
        try:
            self.t_point = float(t_point)
        except ValueError as e:
            print(e)

        # Check that t is in range 0-1
        if t_point > 1 or t_point < 0:
            raise ValueError("t must be in range 0-1")

        # Check that first 2 points required are tuples
        if not isinstance(point1, Tuple) or not isinstance(point2, Tuple):
            raise TypeError("position point must be a tuple with given x and y.")

        # Check that all first 2 points required contain floats
        try:
            self.point1 = (float(point1[0]), float(point1[1]))
            self.point2 = (float(point2[0]), float(point2[1]))
        except PointError as exc:
            print(exc)

        # Parses if 4 points were given
        if point4:

            # Checking all point values are correct type
            if not isinstance(point3, Tuple) or not isinstance(point4, Tuple):
                raise TypeError("position point must be a tuple with given x and y.")

            # Checks that all tuples contain floats
            try:
                self.point1 = (float(point1[0]), float(point1[1]))
                self.point2 = (float(point2[0]), float(point2[1]))
                self.point3 = (float(point3[0]), float(point3[1]))
                self.point4 = (float(point4[0]), float(point4[1]))
            except PointError as exc:
                print(exc)
            pass

        if point3:

            # Checking all point values are correct type
            if not isinstance(point3, Tuple):
                raise TypeError("position point must be a tuple with given x and y.")

            # Checks that all tuples contain floats
            try:
                point1 = (float(point1[0]), float(point1[1]))
                point2 = (float(point2[0]), float(point2[1]))
                point3 = (float(point3[0]), float(point3[1]))
            except PointError as exc:
                print(exc)
            return (point1, point2, point3)

        return (point1, point2)


if __name__ == "__main__":
    point1 = (0, 0)
    point2 = (2, 5)
    point3 = (5, 2)
    point4 = (7, 10)

    # Quadratic Bezier Curve
    # quad_curve = BezierCurve(point1, point2, point3)
    # quad_curve.calculate_points()
    # quad_curve.draw_curve()

    # Cubic Bezier Curve
    cubic_curve = BezierCurve(point1, point2, point3, point4)
    cubic_curve.calculate_points()
    cubic_curve.draw_curve()
