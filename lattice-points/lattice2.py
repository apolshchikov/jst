import numpy as np
import matplotlib.pyplot as plt
from matplotlib.path import Path
import matplotlib.patches as patches
from matplotlib.collections import PatchCollection
import sympy.geometry as g

import itertools as itr

TOLERANCE = float(1e-6)


def quandrant(p: g.Point):
    x = float(p.x)
    y = float(p.y)
    if x < 0:
        if y < 0:
            # Third quadrant
            return (-1, -1)
        else:
            # Second quadrant
            return (-1, 1)
    else:
        if y < 0:
            # Fourth quadrant
            return (1, -1)
        else:
            # First quadrant
            return (1, 1)


class Triangle(g.Triangle):
    def ordering(self):
        angles = self.angles
        sides = self.sides

        ordering = []
        for side in sides:
            ordering.append(float(side.length))
            ordering.append(float(angles[side.p1]))

        return ordering

    @staticmethod
    def compare_ordering(o1, o2):
        return np.array_equal(o1, o2)

    def congruent(self, other):
        pass

    def standardize(self):
        # Fix p1 to the origin and translate triangle
        points = self.vertices
        fixed = points[0]
        points = [p - fixed for p in points]
        t_fix = g.Triangle(*points)

        # Rotate triangle into Q1 and fix side to x-axis
        origin = g.Point(0, 0)
        origin_sides = []
        for s in t_fix.sides:
            if s.p1 == origin or s.p2 == origin:
                origin_sides.append(s)

        xvector = g.Line(p1=origin, p2=g.Point(1, 0))
        angle_segments = []
        for os in origin_sides:

            if os.p1 == origin:
                other = os.p2
            else:
                other = os.p1
            sg_orient = g.Segment(p1=origin, p2=other)
            angle = xvector.angle_between(sg_orient)
            angle_segments.append((float(angle), sg_orient))

        print(angle_segments)
        rotate_angle, rotate_segment = max(angle_segments, key=lambda t: float(t[0]))
        print(rotate_angle, rotate_segment)

        # Check which quadrant the evaluation segment is in. If Q3 or Q4, do nothing to angle; if Q1 or Q2, negate angle

        if float(rotate_segment.p2.y) > 0:
            rotate_angle = -rotate_angle

        rotate_points = t_fix.vertices

        t_rot = g.Triangle(*[i.rotate(rotate_angle) for i in rotate_points])

        # Reflect the triangle about its base midpoint if the angle from the origin is less than or equal to 45
        origin_angle = float(t_rot.angles[origin])

        base_points = []
        alterior_vertex = None
        for v in t_rot.vertices:
            if float(v.y) == float(0):
                base_points.append(v)
            else:
                alterior_vertex = v

        base = g.Segment(*base_points)

        if origin_angle <= (np.pi / 4):
            midpoint = base.midpoint
            x_dist = alterior_vertex.x - midpoint.x
            new_av = g.Point(x=(midpoint.x - x_dist), y=alterior_vertex.y)
            vertices = base_points + [new_av]
            t = g.Triangle(*vertices)
        else:
            t = t_rot

        return t


def generate_triangles(n: int):
    min_x = -np.power(2, float(n+1))
    max_x = np.power(2, float(n+1))

    min_y = min_x
    max_y = max_x

    #points = []
    #for x in np.arange(min_x, max_x, step=1):
    #    for y in np.arange(min_y, max_y, step=1):
    #        points.append(g.Point(x, y))

    x_dist = int(max_x) - int(min_x) + 1
    y_dist = int(max_y) - int(min_y) + 1
    xs = np.linspace(min_x, max_x, num=x_dist)
    ys = np.linspace(min_y, max_y, num=y_dist)

    grid = np.meshgrid(xs, ys, indexing='xy')
    grid = np.array([grid[0], grid[1]])
    n, rows, cols = grid.shape
    points = [g.Point(*tuple(grid[:, i, j])) for i in range(0, rows) for j in range(0, cols)]

    print(len(points))

    combinations = list(itr.combinations(points, r=3))
    row = lambda c_: [float(1), float(c_.x), float(c_.y)]
    matrix = lambda c: np.array([row(c[0]), row(c[1]), row(c[2])])
    valid_combinations = []

    for c in combinations:
        comp1 = (float(c[1].x) - float(c[0].x)) * (float(c[2].y) - float(c[0].y))
        comp2 = (float(c[1].y) - float(c[0].y)) * (float(c[2].x) - float(c[0].x))
        area = float(0.5) * abs(comp1 - comp2)

        if area == float(np.power(2, n)):
            valid_combinations.append(c)


    #for c in combinations:
    #    mat = matrix(c).T
    #    area = float(0.5) * np.linalg.det(mat)
    #    if float(area) == float(np.power(2, n)):
    #        valid_combinations.append(c)


    # triangles = [Triangle(*c) for c in valid_combinations]
    print(len(combinations))
    return valid_combinations


if __name__ == "__main__":
    print(generate_triangles(2))
