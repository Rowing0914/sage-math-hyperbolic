""" this version includes the computation of Dirichlet domain, which is super-expensive...."""

import numpy as np
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
    p_base = PD.get_point(base_pt_x + base_pt_y * I)

    # === Construct the base polygon
    # 1. Construct a polygon in UHP
    center = CC(I)
    polygon = HyperbolicRegularPolygon(num_sides, i_angle, center, {})

    # 2. Conformally transform: UHP -> PD
    points = list()
    for p in polygon._pts:
        _p = moebius_transform(Matrix(2, [1.0, -I, 1.0, I]), p)
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

    """ we don't use this anymore..
    # === Check if the base point is in the base polygon
    table_if_interior = [False] * len(reflect_1st_pBase)
    _d_origin = abs(p_base.coordinates())
    _d_origin = RR(_d_origin)
    for i, p_reflect in enumerate(reflect_1st_pBase):
        _d_p_reflect = p_base.dist(p_reflect)
        print(i, _d_origin, _d_p_reflect, bool(_d_origin <= _d_p_reflect))
        _d_p_reflect = RR(_d_p_reflect)
        table_if_interior[i] = bool(_d_origin <= _d_p_reflect)  # p-base is closer to Origin -> Interior!
    if not all(table_if_interior):
        return []
    # === Check if the base point is in the base polygon
    """

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
    for i in range(len(intersect_p)):
        if_exist = False
        p = intersect_p[i]
        for j in range(i + 1, len(intersect_p)):
            pp = intersect_p[j]
            print(p.dist(pp))
            # if bool(_p.dist(__p) < 0.04):
            if bool(p.dist(pp) < 0.01):  # this value corresponds to float-pt precision
                if_exist = True
                break
        if not if_exist:
            new_intersect_p.append(p)
    intersect_p = new_intersect_p
    # print(len(new_intersect_p))
    # adfs

    table_if_exterior = np.zeros((len(intersect_p), len(reflect_2nd_pBase))).astype(bool)
    for i, p in enumerate(intersect_p):
        for j, p_2nd in enumerate(reflect_2nd_pBase):
            _d_p = p.dist(p_2nd)
            _d_p_base = p.dist(p_base)
            print(i, j, _d_p, _d_p_base, bool(_d_p_base <= _d_p))
            _d_p = RR(_d_p)
            _d_p_base = RR(_d_p_base)
            table_if_exterior[i, j] = bool(_d_p_base <= _d_p)
    print(table_if_exterior)
    ind = [intersect_p[i].coordinates() for i, row in enumerate(table_if_exterior) if np.all(row)]
    # ind = [intersect_p[i].coordinates() for i, row in enumerate(table_if_exterior)]
    return ind


# caching the computation outcomes!
num_search_pt = 20
num_sides = 3
i_angle = pi / 4
base_pt_x = 0.01
base_pt_y = 0.01

ind = process_data(num_sides, i_angle, base_pt_x, base_pt_y)
print(len(ind))
asdf

# Grid search
res = dict()
for base_pt_x in np.linspace(-0.35, 0.35, num_search_pt):
    for base_pt_y in np.linspace(-0.35, 0.35, num_search_pt):
        ind = process_data(num_sides, i_angle, base_pt_x, base_pt_y)
        _key = f"({base_pt_x}, {base_pt_y})"
        print(_key, len(ind))
        res[_key] = len(ind)

# Organise results and store them in CSV
import pandas as pd

df = pd.DataFrame([res])
df = df.T.reset_index()
df.columns = ['Coordinates', '#sides']
df.to_csv("3-gon.csv")
print(df)

# Visualisation
import matplotlib.pyplot as plt

df[['X', 'Y']] = df['Coordinates'].str.strip('()').str.split(', ', expand=True).astype(float)

fig, ax = plt.subplots(figsize=(10, 8))
scatter = ax.scatter(df['X'], df['Y'], c=df['#sides'], cmap='viridis', s=50)
plt.colorbar(scatter, label='Number of Sides')
plt.title('2D Scatter Plot of Polygon Sides')
c = circle((0, 0), 1)
c_matplotlib = c.matplotlib(figure=fig, sub=ax)
ax.axis('equal')  # Set aspect ratio to be equal

plt.show()
