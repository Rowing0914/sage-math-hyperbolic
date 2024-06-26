import matplotlib.colors as mcolors
import matplotlib.pyplot as plt
import numpy as np

PD = HyperbolicPlane().PD()
UHP = HyperbolicPlane().UHP()


@cached_function
def process_data(free_pt_x, free_pt_y, base_pt_x, base_pt_y):
    """ (free_pt_x, free_pt_y) are the coordinate in UHP! """
    # base point
    p_base = PD.get_point(CC(base_pt_x + base_pt_y * I))

    # Cayley transform: UHP -> PD
    z = CC(free_pt_x, free_pt_y)
    z = (z - I) / (z + I)

    # we need this as SageMath proceeds w/h a symbolic expression thus fraction keeps getting complex...
    # and that makes the computation time incredibly long.
    free_pt_x, free_pt_y = float(z.real()), float(z.imag())

    points = [
        PD.get_point(1.0),
        PD.get_point(free_pt_x + free_pt_y * I),
        PD.get_point(-1.0),
        PD.get_point(-I),
    ]
    print("Vertices of base polygon", points)

    # print(f"Vertices of Base Polygon: {[p.coordinates() for p in points]}")
    # n = int(num_sides)
    n = int(len(points))

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
            if i == j:
                reflect_1st_sides.append(s)
            else:
                reflect_1st_sides.append(R * s)
            _diff = abs(i - j)
            _diff = 1 if _diff == (n - 1) else _diff  # modulo operation for cyclic ordering
            diff_reflection_index.append(_diff)

    # reflect base point: n
    reflect_1st_pBase = [R * p_base for R in reflection_1st]

    """
    _sorted = sorted(reflect_1st_pBase, key=lambda x: arg(x.coordinates()))
    _sorted = [_p.coordinates() for _p in _sorted]
    cr = ((_sorted[0] - _sorted[2]) / (_sorted[0] - _sorted[3])) * ((_sorted[1] - _sorted[3]) / (_sorted[1] - _sorted[2]))
    cross_ratio = cr.imag()
    print(cross_ratio, cr)
    asdf
    """

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
        for ind_i, i in enumerate(range(h * n, (h + 1) * n)):  # Apply reflection
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
                    if i != j:
                        # Transform sides
                        _s = reflect_1st_sides[j]
                        reflect_2nd_sides.append(R * _s)

    # Sort by complex-argument
    # reflect_2nd_pBase = sorted(reflect_2nd_pBase, key=lambda x: arg(x.coordinates()))

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

    # avoid the replicates; the above workaround sometime doesn't work so manually double-check
    new_intersect_p = list()
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
            _check = bool(_d_p_base < _d_p) or bool((_d_p_base - _d_p) <= 10 ** -9)
            # print(i, j, _d_p, _d_p_base, _check)
            table_if_exterior[i, j] = _check
    # print(table_if_exterior)
    ind = [intersect_p[i].coordinates() for i, row in enumerate(table_if_exterior) if np.all(row)]
    # ind = [intersect_p[i].coordinates() for i, row in enumerate(table_if_exterior)]

    # Add the diagonals
    sides += [PD.get_geodesic(points[0], points[2]), PD.get_geodesic(points[1], points[3])]

    return p_base, sides, reflect_1st_sides, reflect_1st_pBase, reflect_2nd_sides, reflect_2nd_pBase, list_perp_bisec, diff_index, ind


prev_free_pt_x = None
prev_free_pt_y = None
prev_base_pt_x = None
prev_base_pt_y = None

# caching the computation outcomes!
# free_pt_x, free_pt_y = sqrt(2) / 2, sqrt(2) / 2
free_pt_x, free_pt_y = -1.0, 0.0
base_pt_x, base_pt_y = 0.2, 0.2

p_base, sides, reflect_1st_sides, reflect_1st_pBase, reflect_2nd_sides, reflect_2nd_pBase, list_perp_bisec, diff_index, ind = process_data(
    free_pt_x, free_pt_y, base_pt_x, base_pt_y)


@interact
def _(free_pt_x=free_pt_x, free_pt_y=free_pt_y, base_pt_x=base_pt_x, base_pt_y=base_pt_y, auto_update=False,
      if_plot_sides=True, if_plot_reflect_1st_sides=False, if_plot_reflect_1st_pBase=False,
      if_plot_reflect_2nd_sides=False, if_plot_reflect_2nd_pBase=False, if_plot_perp_bisec=True,
      if_show_dirichletDomain=False,
      ):
    global prev_free_pt_x, prev_free_pt_y, prev_base_pt_x, prev_base_pt_y
    global p_base, sides, reflect_1st_sides, reflect_1st_pBase, reflect_2nd_sides, reflect_2nd_pBase, list_perp_bisec, diff_index, ind

    if_rerun = free_pt_x != prev_free_pt_x or free_pt_y != prev_free_pt_y or base_pt_x != prev_base_pt_x or base_pt_y != prev_base_pt_y
    if if_rerun:  # Check if num_sides has changed
        # Update the previous value of num_sides
        prev_free_pt_x = free_pt_x
        prev_free_pt_y = free_pt_y
        prev_base_pt_x = base_pt_x
        prev_base_pt_y = base_pt_y

        p_base, sides, reflect_1st_sides, reflect_1st_pBase, reflect_2nd_sides, reflect_2nd_pBase, list_perp_bisec, diff_index, ind = process_data(
            free_pt_x, free_pt_y, base_pt_x, base_pt_y)

    # print(f"Found {len(ind)}-gon")

    # plot base point
    P = p_base.show(size=30, color="blue", legend_label='p-base')

    # Plot sides
    if if_plot_sides:
        for i in sides:
            P += i.plot()

    # plot reflected sides
    if if_plot_reflect_1st_sides:
        for i in reflect_1st_sides:
            if i not in sides:
                P += i.plot(color="magenta")

    # plot 1st reflection of points
    if if_plot_reflect_1st_pBase:
        for i in reflect_1st_pBase:
            P += i.show(color="red")

    # plot 2nd reflection of sides
    if if_plot_reflect_2nd_sides:
        for i in reflect_2nd_sides:
            if i not in sides:  # todo: this is not working...
                P += i.plot(color="purple")

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
            P += p.show(color="red", size=20)
        g = hyperbolic_polygon(pts=ind, model="PD", fill=True, alpha=0.0)
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
