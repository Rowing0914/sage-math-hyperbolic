""" this version includes the computation of Dirichlet domain, which is super-expensive...."""

import numpy as np
import multiprocessing

from sage.all import *
from sage.geometry.hyperbolic_space.hyperbolic_interface import HyperbolicPlane
from sage.geometry.hyperbolic_space.hyperbolic_model import moebius_transform
from sage.plot.hyperbolic_regular_polygon import HyperbolicRegularPolygon

PD = HyperbolicPlane().PD()
UHP = HyperbolicPlane().UHP()

""" Note
Float-precision is crucial  to make computation of Dirichlet Domain work...
But, the following precision causes the runtime error of computation of others when sides > 9...
"""
# CC = ComplexField(150)  # don't use this!
# CC = ComplexField(20)  # don't use this!
RR = RealField(5)


def process_data(num_sides, i_angle, base_pt_x, base_pt_y):
    n = int(num_sides)

    # base point
    p_base = PD.get_point(CC(base_pt_x + base_pt_y * I))

    # === Construct the base polygon
    # 1. Construct a polygon in UHP
    center = CC(I)
    polygon = HyperbolicRegularPolygon(num_sides, i_angle, center, {})

    # 2. Conformally transform: UHP -> PD
    points = list()
    for p in polygon._pts:
        _p = moebius_transform(Matrix(2, [CC(1.0), CC(-I), CC(1.0), CC(I)]), p)
        _p = PD.get_point(CC(_p))  # Change the float-point precision
        # _p = p  # for UHP
        points.append(_p)

    # define sides
    sides = list()
    for i in range(n):
        if i + 1 == n:
            _side = PD.get_geodesic(points[i], points[0])
        else:
            _side = PD.get_geodesic(points[i], points[i + 1])
        sides.append(_side)

    # === Check if the base point is in the base polygon
    if base_pt_x == 0: return []  # we don't deal w/h the y-axis (infinity)
    l_base_origin = PD.get_geodesic(p_base, PD.get_point(0.0))
    _cnt = 0
    for l in sides:
        try:
            l_base_origin.intersection(l)[0]
        except:
            _cnt += 1
    if _cnt != len(sides):
        return []
    # === Check if the base point is in the base polygon

    # Get 1st reflection transformations
    reflection_1st = [l.reflection_involution() for l in sides]

    # reflect the sides; fix R and iterate sides; [R1-l1, R1-l2, R1-l3, ...]
    reflect_1st_sides = list()  # n^2
    for i, R in enumerate(reflection_1st):  # n
        for j, s in enumerate(sides):  # n
            if i == j:
                reflect_1st_sides.append(s)
            else:
                reflect_1st_sides.append(R * s)

    # reflect base point: n
    reflect_1st_pBase = [R * p_base for R in reflection_1st]

    # get 2nd reflection transformations: n^2
    reflection_2nd = [l.reflection_involution() for l in reflect_1st_sides]

    # reflect all sides once more
    reflect_2nd_pBase = list()
    # reflect_1st_sides: [R1-l1, R1-l2, R1-l3, ...] (n^2)

    # For Second reflection, eg., R_{i,j}, we compute Diff in i and j to differentiate colours in plot later.
    for h in range(n):
        for ind_i, i in enumerate(range(h * n, (h + 1) * n)):  # Apply reflection
            for j in range(h * n, (h + 1) * n):  # side
                R = reflection_2nd[i]
                if i == j and ind_i != h:  # ind_i != h to avoid reflecting back to p-base
                    # Transform point
                    _p = reflect_1st_pBase[h]
                    _p = R * _p
                    reflect_2nd_pBase.append(_p)
                else:
                    # Transform sides
                    _s = reflect_1st_sides[j]

    # Sort by complex-argument
    reflect_2nd_pBase = sorted(reflect_2nd_pBase, key=lambda x: arg(x.coordinates()))

    # P-bisectors b/w 2nd reflections and base point
    list_perp_bisec = list()
    for p in reflect_2nd_pBase:
        _l = PD.get_geodesic(p_base, p).perpendicular_bisector().complete()
        list_perp_bisec.append(_l)

    # Compute all the intersections of perp-bisectors
    intersect_p = list()
    for i in range(len(list_perp_bisec)):
        for j in range(i + 1, len(list_perp_bisec)):  # avoid the replicates
            try:
                _p = list_perp_bisec[i].intersection(list_perp_bisec[j])[0]
                # _p = list_perp_bisec[j].intersection(list_perp_bisec[i])[0]

                intersect_p.append(_p)
            except Exception as e:
                msg = str(e)
                # print(msg)

    ## avoid the replicates; the above workaround sometime doesn't work so manually double-check
    new_intersect_p = list()
    # print([p.coordinates() for p in intersect_p])
    # print([arg(p.coordinates()) for p in intersect_p])
    for i in range(len(intersect_p)):
        if_exist = False
        p = intersect_p[i]
        for j in range(i + 1, len(intersect_p)):
            pp = intersect_p[j]
            # print(p.dist(pp))
            # if bool(_p.dist(__p) < 0.04):
            if bool(p.dist(pp) < 0.01):  # this value corresponds to float-pt precision
            # if bool(p.dist(pp) < 10**-6):  # this value corresponds to float-pt precision
                if_exist = True
                break
        if not if_exist:
            new_intersect_p.append(p)
    intersect_p = new_intersect_p

    table_if_exterior = np.zeros((len(intersect_p), len(reflect_2nd_pBase))).astype(bool)
    for i, p in enumerate(intersect_p):
        for j, p_2nd in enumerate(reflect_2nd_pBase):
            _d_p = p.dist(p_2nd)
            _d_p_base = p.dist(p_base)
            # print(i, j, _d_p, _d_p_base, bool(_d_p_base <= _d_p))
            _d_p = RR(_d_p)
            _d_p_base = RR(_d_p_base)
            table_if_exterior[i, j] = bool(_d_p_base <= _d_p)
    # print(table_if_exterior)
    ind = [intersect_p[i].coordinates() for i, row in enumerate(table_if_exterior) if np.all(row)]
    # ind = [intersect_p[i].coordinates() for i, row in enumerate(table_if_exterior)]
    return ind

def process_data_wrapper(args):
    num_sides, i_angle, base_pt_x, base_pt_y = args
    return process_data(num_sides, i_angle, base_pt_x, base_pt_y)

if __name__ == '__main__':
    # caching the computation outcomes!
    num_search_pt = 100
    num_sides = 5
    i_angle = pi / 4
    base_pt_x = 0.01
    base_pt_y = 0.01

    # ind = process_data(num_sides, i_angle, base_pt_x, base_pt_y)
    # print(len(ind))
    # asdf

    # === Grid search
    import time, datetime
    start = time.time()
    print(datetime.datetime.now())
    pool = multiprocessing.Pool()
    results = {}
    search_space = np.linspace(-0.5, 0.5, num_search_pt)
    for base_pt_x in search_space:
        for base_pt_y in search_space:
            args = (num_sides, i_angle, base_pt_x, base_pt_y)
            _key = f"({base_pt_x}, {base_pt_y})"
            results[_key] = pool.apply_async(process_data_wrapper, (args,))

    pool.close()
    pool.join()

    for k, v in results.items():
        ind = v.get()
        results[k] = len(ind)
        # print(k, len(ind))
    print(f"{datetime.datetime.now()}: Took {time.time() - start}")

    # Organise results and store them in CSV
    import pandas as pd

    df = pd.DataFrame([results])
    df = df.T.reset_index()
    df.columns = ['Coordinates', '#sides']
    df.to_csv(f"{num_sides}-gon.csv")
    print(df)

