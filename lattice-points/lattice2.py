import numpy as np
import matplotlib.pyplot as plt
from matplotlib.path import Path
import matplotlib.patches as patches
from matplotlib.collections import PatchCollection
import sympy.geometry as g
import time

import timeit

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


def generate_lattice(min_x: int, min_y: int, max_x: int, max_y: int) -> np.ndarray:
    g = np.mgrid[min_x:max_x + 1, min_y:max_y + 1]
    pos = np.vstack(map(np.ravel, g)).T

    return pos


def check_linear_arr(arrs):
    verification = lambda x: x == np.repeat(np.array([x[:, 0]]), 3, axis=0).T
    reduction = lambda v: np.all(v, axis=1)
    test = lambda t: reduction(verification(t))

    xs = test(arrs[:, :, 0])
    ys = test(arrs[:, :, 1])
    valid = ~np.any(np.array([xs, ys]).T, axis=1)

    return arrs[valid]


def calculate_area_arr(arrs, area_requirement):
    test = np.append(np.ones((arrs.shape[0], arrs.shape[1], 1)), arrs, axis=2)
    dets = float(0.5) * np.absolute(np.linalg.det(test))

    index_pos = np.isclose(dets, area_requirement, atol=TOLERANCE)

    return arrs[index_pos]


def generate_triangle_info(v):
    crs = np.cross(v[:2], v[2:])
    s = np.linalg.norm(crs)
    if crs > 0:
        c = -np.dot(v[:2], v[2:])
    else:
        c = np.dot(v[:2], v[2:])
    return np.round(np.arctan2(s, c), 30)


def triangle_data(arrs):
    shifted = np.roll(arrs, 1, axis=1)
    vertices = shifted - arrs
    vertices[:, 2, :] = -vertices[:, 2, :]

    sizes = np.apply_along_axis(np.linalg.norm, 2, vertices)
    sizes = np.round(sizes, 30).reshape(arrs.shape[0], 3, 1)
    sizes = np.sort(sizes, axis=1)
    s_vertices = np.roll(vertices, 1, axis=1)

    m_vertices = np.concatenate([vertices, s_vertices], axis=2)

    reverse1 = m_vertices[:, 2, :][:, :2]
    reverse2 = m_vertices[:, 2, :][:, 2:]

    reversed_c = np.concatenate([reverse2, reverse1], axis=1)

    m_vertices[:, 2, :] = reversed_c

    angles = np.apply_along_axis(generate_triangle_info, 2, m_vertices).reshape(arrs.shape[0], 3, 1)
    angles = np.sort(angles, axis=1)

    obtuse = np.where(angles >= np.pi / 2)
    ob_arr = np.zeros(angles.shape)
    ob_arr[obtuse] = True
    mask = np.any(ob_arr, axis=1)

    print(mask)


def originate_triangles(arrs):
    transform = np.array(arrs[:, 0, :]).repeat(3, axis=0).reshape(arrs.shape)
    originated = arrs - transform

    print(arrs, originated)


def triangle(arr: np.ndarray) -> np.ndarray:
    shifted = np.roll(arr, -1, axis=0)
    v = shifted - arr
    v = np.array([v[0], v[1], -v[2]])
    norms = np.apply_along_axis(np.linalg.norm, 1, v)
    shifted_v = np.roll(v, -1, axis=0)
    corresponding_norms = np.roll(norms, 1, axis=0)

    angles = []
    for k, sv in zip(v, shifted_v):
        cos_th = np.dot(k, sv) / (np.linalg.norm(k) * np.linalg.norm(sv))
        angles.append(np.arccos(cos_th))

    angles = np.array(angles)

    angle_test = len(np.where(angles >= np.pi / 2)[0])

    err = 64
    output = np.array([np.round(angles, err), np.round(corresponding_norms, err)]).T
    output = np.sort(output, axis=0)

    alternate_output = np.empty(output.shape)
    alternate_output[:, :] = -1

    if angle_test != 0:
        return alternate_output

    return output


def triangle_arr(arrs):
    result = np.empty(arrs.shape)
    for i in range(len(arrs)):
        result[i] = triangle(arrs[i])

    return result


def generate_combinations(points: np.ndarray, n: int, verbose=False) -> np.ndarray:
    combinations = np.array(list(itr.combinations(points, r=3)))
    v1 = time.time()
    valid_combinations = check_linear_arr(combinations)
    v2 = time.time()
    valid_combinations = calculate_area_arr(valid_combinations, np.power(2, n))
    v3 = time.time()

    if verbose:
        print("Checking Linear: {0}s\nChecking Area: {1}s".format(np.round(v2 - v1, 2), np.round(v3 - v2, 2)))

    return np.array(valid_combinations)


def valid_triangles(combinations: np.ndarray) -> np.ndarray:
    triangles = triangle_arr(combinations)
    triangles = np.sort(triangles, axis=0)
    un = np.unique(triangles, axis=0)

    return un


def generate_triangles(n: int, verbose=False) -> [int, int, int]:
    min_x = -np.power(2, float(n + 1))
    max_x = np.power(2, float(n + 1))

    min_y = min_x
    max_y = max_x

    t1 = time.time()
    points = generate_lattice(min_x, min_y, max_x, max_y)
    t2 = time.time()
    combinations = generate_combinations(points, n, verbose=verbose)
    t3 = time.time()
    triangles = valid_triangles(combinations)
    t4 = time.time()

    if verbose:
        print("Generating Lattices: {0}s\nGenerating Combinations: {1}s\nGenerating Triangles: {2}s".format(
            np.round(t2 - t1, 2), np.round(t3 - t2, 2), np.round(t4 - t3, 2)))

    return n, len(combinations), len(triangles)


if __name__ == "__main__":
    k = range(1, 4)

    # for z in k:
    #     print(generate_triangles(z, verbose=True))

    g = generate_lattice(-2, -2, 2, 2)
    combs = generate_combinations(g, 1)
    z = combs[1:3]

    print(z)
    triangle_data(z)

