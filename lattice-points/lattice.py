import numpy as np
import matplotlib.pyplot as plt
from matplotlib.path import Path
import matplotlib.patches as patches
from matplotlib.collections import PatchCollection
import sympy.geometry as g

import itertools as itr

TOLERANCE = float(1e-6)


def _generate_base_region(n: int):
    # Generate the region where all bases can lie in the first quadrant. The region is a quarter circle.
    # Find all lattice points within the region
    min_x = 1
    max_x = np.power(2, float(n + 2) / float(2))
    close = np.isclose([abs(max_x - float(np.round(max_x)))], [0], atol=TOLERANCE)

    if close:
        max_x = max_x + 1
    else:
        max_x = np.floor(max_x)

    min_y = 0

    points = []
    vectors = []
    for x in np.arange(min_x, max_x, step=1):
        max_y = np.sqrt(np.power(max_x, 2) - np.power(x, 2))
        std_max_y = np.floor(max_y) + 1
        for y in np.arange(min_y, std_max_y, step=1):
            origin = g.Point(0, 0)
            p = g.Point(x, y)
            v = g.Segment(p1=origin, p2=p)
            points.append(p)
            vectors.append(v)

    return vectors


def generate_base_region(n: int):
    min_x = -np.power(2, float(n + 1))
    max_x = np.power(2, float(n + 1))

    min_y = min_x
    max_y = max_x

    points = []
    for x in np.arange(min_x, max_x, step=1):
        for y in np.arange(min_y, max_y, step=1):
            points.append(g.Point(x, y))

    combinations = itr.combinations(points, r=3)




def generate_parallelogram(p1, p2, h):
    l1 = g.Segment(p1=p1, p2=p2)
    l2 = l1.perpendicular_line(p1)
    # l2 = g.Segment(l2.p1, l2.p2).unit()
    l4 = l1.perpendicular_line(p2)

    # Get x3 point
    if float(p2.y) == float(0):
        x3 = 0
        y3 = h
    elif float(p2.x) == float(0):
        x3 = -h
        y3 = 0
    else:
        numerator = np.power(h, 2)
        denominator = 1 + np.power((float(p2.x) / float(p2.y)), 2)
        x3 = np.sqrt((numerator / denominator))
        y3 = l2.slope * (x3)
    p3 = g.Point(x3, y3)

    l3 = l2.perpendicular_line(p3)

    p4 = l3.intersection(l4)[0]

    return p1, p2, p3, p4


def ordering(t: g.Triangle):
    angles = t.angles
    sides = t.sides

    ordering = []

    for side in sides:
        ordering.append(side.length)
        ordering.append(angles[side.p1])

    return ordering


def compare_ordering(o1, o2):
    o1 = [float(i) for i in o1]
    o2 = [float(i) for i in o2]
    return np.array_equal(o1, o2)


def congruent(t1: g.Triangle, t2: g.Triangle):
    t1_ordering = ordering(t1)
    t2_ordering = ordering(t2)
    # rot = lambda x: [x[4], x[5], x[0], x[1], x[2], x[3]]
    rot = lambda x: np.roll(x, shift=2)
    ref_x = lambda x: [x[0], x[5], x[4], x[3], x[2], x[1]]

    s1 = t2_ordering
    s2 = rot(s1)
    s3 = rot(s2)
    r1 = ref_x(s1)
    r2 = ref_x(s2)
    r3 = ref_x(s3)

    convert = lambda x: [float(i) for i in x]

    t2_orderings = list(map(convert, [s1, s2, s3, r1, r2, r3]))
    t1_orderings = list(map(convert, [t1_ordering] * len(t2_orderings)))

    return np.array_equal(t1_orderings, t2_orderings)


def generate_triangles(n: int, bases: [g.Segment]):
    area = np.power(2, n)

    triangles = []
    for base in bases:
        height = (2 * area) / float(base.length)
        try:
            p1, p2, p3, p4 = generate_parallelogram(base.p1, base.p2, height)
            al = g.Segment(p3, p4)

            y = lambda x: al.slope * (x - float(p3.x)) + float(p3.y)
            min_x = np.ceil(float(p3.x))
            max_x = np.floor(float(p4.x))

            for x in np.arange(min_x, max_x, step=1):
                x_test = float(x) % 1 == 0
                y_p = y(x)
                y_test = float(y_p) % 1 == 0
                if x_test and y_test:
                    p3 = g.Point(x, y_p)
                    t = g.Triangle(base.p1, base.p2, p3)
                    triangles.append(t)
        except ValueError:
            pass

    acute_output = []
    for t in triangles:
        append = True
        for vertex, angle in t.angles.items():
            if float(angle) > np.pi / 2 or angle <= 0:
                append = False

        if append:
            acute_output.append(t)

    try:
        output = [acute_output[0]]
        for ao in acute_output[1:]:
            additional = []
            for o in output:
                if not congruent(ao, o):
                    additional.append(ao)
            new_output = list(set(output + additional))
            output = new_output

        output = list(set(output))
    except IndexError:
        output = []

    return output


if __name__ == "__main__":
    for i in range(1, 4):
        vectors = generate_base_region(i)
        triangles = generate_triangles(i, vectors)

        print(i, len(triangles))

    # print(_generate_base_region(2))
