import matplotlib.colors as mcolors
import matplotlib.pyplot as plt
import numpy as np
from sage.all import *
from sage.geometry.hyperbolic_space.hyperbolic_interface import HyperbolicPlane
from sage.plot.hyperbolic_regular_polygon import HyperbolicRegularPolygon

PD = HyperbolicPlane().PD()


@cached_function
def process_data(num_sides, i_angle, base_pt_x, base_pt_y):
    n = int(num_sides)

    # === Construct the base polygon
    # 1. Construct a polygon in UHP
    center = I
    polygon = HyperbolicRegularPolygon(num_sides, i_angle, center, {})

    # 2. Conformally transform: UHP -> PD
    from sage.geometry.hyperbolic_space.hyperbolic_model import moebius_transform
    points = list()
    for p in polygon._pts:
        _p = moebius_transform(Matrix(2, [1, -I, 1, I]), p)
        points.append(_p)

    # define sides
    sides = list()
    for i in range(n):
        if i + 1 == n:
            _side = PD.get_geodesic(points[i], points[0])
        else:
            _side = PD.get_geodesic(points[i], points[i + 1])
        sides.append(_side)

    # Get 1st reflection transformations
    reflection_1st = [l.reflection_involution() for l in sides]

    # reflect the sides; fix R and iterate sides; [R1-l1, R1-l2, R1-l3, ...]
    reflect_1st_sides = list()  # n^2
    diff_reflection_index = list()
    for i, R in enumerate(reflection_1st):  # n
        for j, s in enumerate(sides):  # n
            reflect_1st_sides.append(R * s)
            _diff = abs(i - j)
            _diff = 1 if _diff == (n - 1) else _diff  # modulo operation for cyclic ordering
            diff_reflection_index.append(_diff)

    # base point
    x, y = round(base_pt_x, 2), round(base_pt_y, 2)
    p_base = PD.get_point(x + y * I)

    # reflect base point: n
    reflect_1st_pBase = [R * p_base for R in reflection_1st]

    # get 2nd reflection transformations: n^2
    reflection_2nd = [l.reflection_involution() for l in reflect_1st_sides]

    # reflect all sides once more
    reflect_2nd_sides = list()
    reflect_2nd_pBase = list()
    # reflection_2nd: [R1R1, R2R1, R3R1, ...] (n^2)
    # reflect_1st_sides: [R1-l1, R1-l2, R1-l3, ...] (n^2)

    # For Second reflection, eg., R_{i,j}, we compute Diff in i and j to differentiate colours in plot later.
    diff_index = list()
    for h in range(n):
        for ind_i, i in enumerate(range(h * n, (h + 1) * n)):  # reflection transf
            for j in range(h * n, (h + 1) * n):  # side
                _diff = diff_reflection_index[i]
                R = reflection_2nd[i]
                if i == j and ind_i != h:  # ind_i != h to avoid reflecting back to p-base
                    # Transform point
                    _p = reflect_1st_pBase[h]
                    _p = R * _p
                    reflect_2nd_pBase.append(_p)
                    diff_index.append(_diff)
                else:
                    # Transform sides
                    _s = reflect_1st_sides[j]
                    reflect_2nd_sides.append(R * _s)

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
                # avoid the replicates; the above workaround sometime doesn't work so manually double-check
                if_exist = False
                for __p in intersect_p:
                    if bool(abs(_p.dist(__p)) < 10 ** -9):
                        if_exist = True
                        break
                if not if_exist:
                    intersect_p.append(_p)
            except:
                continue

    ind = list()
    table_if_exterior = np.zeros((len(reflect_2nd_pBase), len(intersect_p))).astype(bool)
    for i, p_2nd in enumerate(reflect_2nd_pBase):
        for j, p in enumerate(intersect_p):
            _d_p = p_2nd.dist(p)
            _d_p_base = p_base.dist(p_base)
            table_if_exterior[i, j] = _d_p > _d_p_base
    ind = [intersect_p[i].coordinates() for i, row in enumerate(table_if_exterior) if np.alltrue(row)]
    print(table_if_exterior)
    print(ind)

    # ==== Method 1: Intersection w/h the base point
    # # Compute all the intersections of perp-bisectors
    # intersect_p_l = list()
    # for i in range(len(list_perp_bisec)):
    #     for j in range(i + 1, len(list_perp_bisec)):  # avoid the replicates
    #         try:
    #             _p = list_perp_bisec[i].intersection(list_perp_bisec[j])[0]
    #             _l = PD.get_geodesic(p_base, _p)
    #
    #             # avoid the replicates; the above workaround sometime doesn't work so manually double-check
    #             if_exist = False
    #             for (__p, __l) in intersect_p_l:
    #                 if bool(abs(_p.dist(__p)) < 10 ** -9):
    #                     if_exist = True
    #                     break
    #             if not if_exist:
    #                 intersect_p_l.append((_p, _l))
    #         except:
    #             continue
    #
    # ind = list()
    # debug_pt = []
    # debug_line = []
    #
    # # p, l = intersect_p_l[1]
    # # debug_pt = [p]
    # # debug_line = [l]
    # # cnt_not_intersect = 0
    # # for i in range(len(list_perp_bisec)):
    # #     try:
    # #         # _p = list_perp_bisec[i].angle(l)
    # #         _p = l.intersection(list_perp_bisec[i])[0]
    # #         # _p = list_perp_bisec[i].intersection(l)[0]
    # #         debug_pt.append(_p)
    # #         debug_line.append(list_perp_bisec[i])
    # #         print(p, _p, _p.dist(p), bool(abs(_p.dist(p)) < 10 ** -9), arg(p.coordinates()))
    # #         if bool(abs(_p.dist(p)) < 10 ** -9):  # if intersect by itself at the same point
    # #             cnt_not_intersect += 1
    # #     except:
    # #         cnt_not_intersect += 1
    # #         pass
    # # print(p, cnt_not_intersect)
    # # # asdf
    # #
    # for _intersect in intersect_p_l:
    #     p, l = _intersect
    #     cnt_not_intersect = 0
    #     for i in range(len(list_perp_bisec)):
    #         try:
    #             # _p = l.intersection(list_perp_bisec[i])[0]  # bug: this has a different behaviour than the below...
    #             _p = list_perp_bisec[i].intersection(l)[0]
    #             if bool(abs(_p.dist(p)) < 10 ** -9):  # if intersect by itself at the same point
    #                 cnt_not_intersect += 1
    #         except:
    #             cnt_not_intersect += 1
    #             pass
    #     # print(cnt_not_intersect, p)
    #     if cnt_not_intersect == len(list_perp_bisec):
    #         ind.append(p.coordinates())
    return p_base, sides, reflect_1st_sides, reflect_1st_pBase, reflect_2nd_sides, reflect_2nd_pBase, list_perp_bisec, diff_index, ind
    # return p_base, sides, reflect_1st_sides, reflect_1st_pBase, reflect_2nd_sides, reflect_2nd_pBase, list_perp_bisec, diff_index, ind, debug_pt, debug_line


prev_num_sides = None
prev_i_angle = None
prev_base_pt_x = None
prev_base_pt_y = None

# caching the computation outcomes!
num_sides = 3
i_angle = pi / 4
base_pt_x = 0.0
base_pt_y = 0.0
p_base, sides, reflect_1st_sides, reflect_1st_pBase, reflect_2nd_sides, reflect_2nd_pBase, list_perp_bisec, diff_index, ind = process_data(
    num_sides, i_angle, base_pt_x, base_pt_y)


# p_base, sides, reflect_1st_sides, reflect_1st_pBase, reflect_2nd_sides, reflect_2nd_pBase, list_perp_bisec, diff_index, ind, debug_pt, debug_line = process_data(num_sides, i_angle, base_pt_x, base_pt_y)


@interact
def _(num_sides=3, i_angle=pi / 4, base_pt_x=0.0, base_pt_y=0.0, auto_update=False,
      if_plot_sides=True, if_plot_reflect_1st_sides=False, if_plot_reflect_1st_pBase=False,
      if_plot_reflect_2nd_sides=False, if_plot_reflect_2nd_pBase=False, if_plot_perp_bisec=True,
      if_show_dirichletDomain=False,
      ):
    global prev_num_sides, prev_i_angle, prev_base_pt_x, prev_base_pt_y
    global p_base, sides, reflect_1st_sides, reflect_1st_pBase, reflect_2nd_sides, reflect_2nd_pBase, list_perp_bisec, diff_index, ind
    # global p_base, sides, reflect_1st_sides, reflect_1st_pBase, reflect_2nd_sides, reflect_2nd_pBase, list_perp_bisec, diff_index, ind, debug_pt, debug_line

    if_rerun = num_sides != prev_num_sides or i_angle != prev_i_angle or base_pt_x != prev_base_pt_x or base_pt_y != prev_base_pt_y
    if if_rerun:  # Check if num_sides has changed
        # Update the previous value of num_sides
        prev_num_sides = num_sides
        prev_i_angle = i_angle
        prev_base_pt_x = base_pt_x
        prev_base_pt_y = base_pt_y

        p_base, sides, reflect_1st_sides, reflect_1st_pBase, reflect_2nd_sides, reflect_2nd_pBase, list_perp_bisec, diff_index, ind = process_data(
            num_sides, i_angle, base_pt_x, base_pt_y)
        # p_base, sides, reflect_1st_sides, reflect_1st_pBase, reflect_2nd_sides, reflect_2nd_pBase, list_perp_bisec, diff_index, ind, debug_pt, debug_line = process_data(num_sides, i_angle, base_pt_x, base_pt_y)

    # plot base point
    P = p_base.show(size=50, color="blue", legend_label='p-base')

    # if len(debug_pt) > 0:
    #     for i in debug_pt:
    #         P += i.show(color="green")
    # if len(debug_line) > 0:
    #     for i in debug_line:
    #         P += i.plot()

    # Plot sides
    if if_plot_sides:
        for i in sides:
            P += i.plot()

    # plot reflected sides
    if if_plot_reflect_1st_sides:
        for i in reflect_1st_sides:
            P += i.plot()

    # plot 1st reflection of points
    if if_plot_reflect_1st_pBase:
        for i in reflect_1st_pBase:
            P += i.show(color="red")

    # plot 2nd reflection of sides
    if if_plot_reflect_2nd_sides:
        for i in reflect_2nd_sides:
            P += i.plot()

    # plot 2nd reflection of points
    if if_plot_reflect_2nd_pBase:
        for i in reflect_2nd_pBase:
            P += i.show(color="red")

    # plot p-bisectors for 2nd reflections
    if if_plot_perp_bisec:
        used_colors = dict()
        cmap = plt.colormaps['tab10']
        for i, _diff in zip(list_perp_bisec, diff_index):
            _c = cmap(_diff)
            _c = mcolors.to_hex(_c)
            P += i.plot(color=_c)
            _name = f"L_i_i+{_diff}"
            if _name not in used_colors:
                used_colors[_name] = _c

    if if_show_dirichletDomain:
        ind = sorted(ind, key=arg)  # sort the vertices by complex-argument for correct instantiation of polygon!
        for p in ind:
            p = PD.get_point(coordinates=p)
            P += p.show(color="red", size=30)
        g = hyperbolic_polygon(pts=ind, model="PD", fill=True, alpha=0.2)
        P += g.plot()

    P.show(axes=True)

    if if_plot_perp_bisec:
        num_colors = len(used_colors)
        fig, ax = plt.subplots(figsize=(6, 1))

        # Display the custom colors in a table
        for i, (key, color) in enumerate(used_colors.items()):
            ax.add_patch(plt.Rectangle((i - 0.5, -0.5), 1, 1, color=color, ec='black'))
            ax.text(i, 0, key, ha='center', va='center', color='black', fontsize=8)

        # Remove axis ticks and set axis limits
        ax.set_xticks([])
        ax.set_yticks([])
        ax.set_xlim(-0.5, num_colors - 0.5)
        ax.set_ylim(-0.5, 0.5)

        # Show plot
        plt.show()
