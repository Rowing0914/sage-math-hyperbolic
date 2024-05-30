import numpy as np

KM = HyperbolicPlane().KM()
UHP = HyperbolicPlane().UHP()
PD = HyperbolicPlane().PD()


@cached_function
def process_data(r, num):
    _x, _y = (r - I) / (r + I)  # Cayley transform: UHP -> PD
    points = [
        PD.get_point(1.0),
        # we can't do this as SageMath proceeds w/h a symbolic expression thus fraction keeps getting complex...
        # and this makes the computation time incredibly long.
        PD.get_point(CC(_x, _y)),
        PD.get_point(-1.0),
        PD.get_point(-I),
    ]
    print("Vertices of base polygon", points)
    n = int(len(points))

    # define sides
    sides = list()
    for i in range(n):
        if i + 1 == n:
            _side = PD.get_geodesic(points[i], points[0])
        else:
            _side = PD.get_geodesic(points[i], points[i + 1])
        sides.append(_side)

    # ============ Visualisation in R^2
    x, y = var('x, y')
    eqn = -r * x + x ^ 2 + r * x ^ 2 - x ^ 3 + y ^ 2 + r * y ^ 2 - x * y ^ 2 == 0

    # q[r] == 0 for Cross-ratio
    plots = [implicit_plot(eqn, (x, -2, 2), (y, -2, 2), title=f"CR (r={r})")]

    # w/h other plots
    canvas = list()
    eqs = [(x - 1, "orange"),
           (x - r, "orange"),
           ((x - 0.5) ^ 2 + y ^ 2 - 0.5 ^ 2, "green"),
           ((x + 0.5) ^ 2 + y ^ 2 - 0.5 ^ 2, "pink"),
           (eqn, "red")]
    for eq, c in eqs:
        canvas.append(implicit_plot(eq, (x, -2, 2), (y, -2, 2), color=c))
    plots.append(sum(canvas))  # https://sagecell.sagemath.org/?q=ynscjh

    # ============ Visualisation in Poincare Disk
    # Ys = np.concatenate([np.linspace(start=-1 - r, stop=1, num=num // 2), np.linspace(start=1.01, stop=10, num=num // 2)])
    # Ys = np.concatenate([np.linspace(start=-1 - r, stop=1, num=num // 2), np.linspace(start=1.01, stop=10, num=num // 2)])
    Ys = np.linspace(start=0.0, stop=10, num=num)

    # for _x in np.linspace(start=-1 - r, stop=1, num=num):
    #     print(solve(eqn.subs({y: _x}), x))
    #     print(solve(eqn.subs({x: _x}), y))
    # asdf

    samples = []
    # for _x in Xs:
    #     samples += [CC(_x, s.rhs().n().imag()) for s in solve(eqn.subs({x: _x}), y)]

    for _y in Ys:
        # _samples = [CC(s.rhs().n().real(), _y) for s in solve(eqn.subs({y: _y}), x)[:2]]
        solutions = solve(eqn.subs({y: _y}), x)
        _samples = [CC(s.rhs().n().real(), _y) for s in [solutions[0], solutions[1]]]
        # print(_y, solve(eqn.subs({y: _y}), x))
        # print(_y, [s.rhs().n() for s in solve(eqn.subs({y: _y}), x)])
        print(_y, _samples)
        samples += _samples

    pts = list()

    for p in samples:
        if UHP.point_in_model(p):
            p = (p - I) / (p + I)  # Cayley transform: UHP -> PD
            pts.append(PD.get_point(p))

    return sides, pts, plots


prev_r = None

# caching the computation outcomes!
r = -1.25
num = 30
sides, pts, plots = process_data(r, num)


@interact
def _(r=r, num=num, auto_update=False):
    global prev_r, prev_num
    global sides, pts, plots

    if_rerun = r != prev_r or num != prev_num
    if if_rerun:  # Check if num_sides has changed
        prev_r = r
        prev_num = num
        sides, pts, plots = process_data(r, num)

    P = Graphics()
    for i in sides:
        P += i.plot()

    for i in pts:
        P += i.show(color="red")

    show(graphics_array([plots[0], plots[1]], 1, 2))
    P.show(axes=True)
