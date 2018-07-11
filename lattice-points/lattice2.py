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

        print(t_fix)

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

    def to_array(self):
        vertices = [[float(v.x), float(v.y)] for v in self.vertices]
        return np.array(vertices)


def generate_lattice(min_x, min_y, max_x, max_y):
    g = np.mgrid[min_x:max_x + 1, min_y:max_y + 1]
    pos = np.vstack(map(np.ravel, g)).T

    return pos


def check_linear(arr):
    test = np.all(arr == arr[0, :], axis=0)
    return any(test)


def calculate_area(arr):
    ones = np.ones((3, 1))
    t = np.hstack((ones, arr))
    area = float(0.5) * abs(np.linalg.det(t.T))

    return area


def triangle(arr):
    v = np.roll(arr, 1) - arr

    print(v)


def generate_combinations(points, n):
    combinations = np.array(list(itr.combinations(points, r=3)))
    valid_combinations = []

    for c in combinations:
        check_nonlinear = not check_linear(c)

        area = calculate_area(c)
        check_area = np.isclose([area - np.power(2, n)], [0], atol=TOLERANCE)

        checks = [check_nonlinear, check_area]
        if all(checks):
            valid_combinations.append(c)

    return np.array(valid_combinations)


def valid_triangles(combinations):
    vt = []
    for c in combinations:
        points = [g.Point(*xy) for xy in c]
        t = Triangle(*points)
        print(c)
        print(t.standardize())


def generate_triangles(n: int):
    min_x = -np.power(2, float(n))
    max_x = np.power(2, float(n))

    min_y = min_x
    max_y = max_x

    points = generate_lattice(min_x, min_y, max_x, max_y)
    combinations = generate_combinations(points, n)

    return len(combinations)


if __name__ == "__main__":
    k = 1
    min_x = -np.power(2, k)
    max_x = np.power(2, k)
    min_y = min_x
    max_y = max_x

    ps = generate_lattice(min_x, min_y, max_x, max_y)
    combs = generate_combinations(ps, k)
    t1 = combs[1:6]
    # valid_triangles(t1)
    triangle(t1[0])
    # [print(generate_triangles(i)) for i in range(1, 3)]
