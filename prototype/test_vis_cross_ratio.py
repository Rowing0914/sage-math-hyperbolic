# x, y, r = var('x, y, r')
# # poly = ((x - y * I) ^ 2 - 2 * (x - y * I) + 1) * ((x - y * I) ^ 2 - 2 * r * (x - y * I) + r ^ 2) * ((x + y * I) ^ 2 - 2 * (x + y * I) + r) * ((x + y * I) ^ 2 - 2 * r * (x + y * I) + r)
#
# # irreducibility check
# p = -2 * r ^ 2 * x * y + 4 * r ^ 3 * x * y - 2 * r ^ 4 * x * y + 2 * r * x ^ 2 * y - 2 * r ^ 2 * x ^ 2 * y - 2 * r ^ 3 * x ^ 2 * y + 2 * r ^ 4 * x ^ 2 * y - 2 * x ^ 4 * y + 2 * r * x ^ 4 * y + 2 * r ^ 2 * x ^ 4 * y - 2 * r ^ 3 * x ^ 4 * y + 2 * x ^ 5 * y - 4 * r * x ^ 5 * y + 2 * r ^ 2 * x ^ 5 * y + 2 * r * y ^ 3 - 2 * r ^ 2 * y ^ 3 - 2 * r ^ 3 * y ^ 3 + 2 * r ^ 4 * y ^ 3 - 4 * x ^ 2 * y ^ 3 + 4 * r * x ^ 2 * y ^ 3 + 4 * r ^ 2 * x ^ 2 * y ^ 3 - 4 * r ^ 3 * x ^ 2 * y ^ 3 + 4 * x ^ 3 * y ^ 3 - 8 * r * x ^ 3 * y ^ 3 + 4 * r ^ 2 * x ^ 3 * y ^ 3 - 2 * y ^ 5 + 2 * r * y ^ 5 + 2 * r ^ 2 * y ^ 5 - 2 * r ^ 3 * y ^ 5 + 2 * x * y ^ 5 - 4 * r * x * y ^ 5 + 2 * r ^ 2 * x * y ^ 5
# factor(p)
#
# # irreducibility check
# q = -r * x + x ^ 2 + r * x ^ 2 - x ^ 3 + y ^ 2 + r * y ^ 2 - x * y ^ 2
# factor(q)

x, y, r = var('x, y, r')


def q(r):
    # define the polynomial w.r.t. r
    return -r * x + x ^ 2 + r * x ^ 2 - x ^ 3 + y ^ 2 + r * y ^ 2 - x * y ^ 2


print(q(r=-1) == 0)
# mathematica: ContourPlot[-x ^ 3 - x * y ^ 2 + x == 0, {x, -2, 2}, {y, -2, 2}]
# implicit_plot(q(r=-1) == 0, (x, -2, 2), (y, -2, 2))

canvas = list()
for eq in [x - 1, x + 1, (x - 0.5) ^ 2 + y ^ 2 - 0.5 ^ 2, (x + 0.5) ^ 2 + y ^ 2 - 0.5 ^ 2, q(r=-1) == 0]:
    canvas.append(implicit_plot(eq, (x, -2, 2), (y, -2, 2)))
sum(canvas)  # https://sagecell.sagemath.org/?q=ynscjh

# ContourPlot[{x - 1, x + 1, (x - 0.5)^2 + y^2 - 0.5^2, (x + 0.5)^2 + y^2 - 0.5^2, -x^3 - x*y^2 + x == 0}, {x, -2, 2}, {y, -2, 2}]

# ======= Visualisation
x, y, r = var('x, y, r')


def q(r):
    # define the polynomial w.r.t. r
    return -r * x + x ^ 2 + r * x ^ 2 - x ^ 3 + y ^ 2 + r * y ^ 2 - x * y ^ 2


# eqn = q(r=-1) == 0
# eqn = q(r=-1.02) == 0
# eqn = q(r=-1.25) == 0
eqn = q(r=-2.0) == 0

import numpy as np

_num = 50
Xs = np.linspace(start=-1, stop=1, num=_num)
Ys = np.concatenate([np.linspace(start=-1, stop=1, num=_num // 2), np.linspace(start=1.01, stop=10, num=_num // 2)])

samples = []
for _x in Xs:
    samples += [CC(_x, s.rhs().n().imag()) for s in solve(eqn.subs({x: _x}), y)]

for _y in Ys:
    samples += [CC(s.rhs().n().real(), _y) for s in solve(eqn.subs({y: _y}), x)]

print(samples)
KM = HyperbolicPlane().KM()
UHP = HyperbolicPlane().UHP()
PD = HyperbolicPlane().PD()

pts = list()

""" R^2 -> Klein Model -> PD
Ref: https://en.wikipedia.org/wiki/Coordinate_systems_for_the_hyperbolic_plane#Beltrami_coordinates

for p in samples:
    _x, _y = tanh(p.real()), tanh(p.imag())  # R^2 -> Klein Model
    if KM.point_in_model(CC(_x, _y)):
        # https://en.wikipedia.org/wiki/Poincar%C3%A9_disk_model#Relation_to_other_models_of_hyperbolic_geometry
        p = CC(_x / (1 + sqrt(1 - _x ^ 2 - _y ^ 2)), _y / (1 + sqrt(1 - _x ^ 2 - _y ^ 2)))  # KM -> PD
        pts.append(PD.get_point(p))
"""

for p in samples:
    if UHP.point_in_model(p):
        p = (p - I) / (p + I)  # Cayley transform: UHP -> PD
        pts.append(PD.get_point(p))

print(pts)
g = Graphics()
for p in pts:
    g += p.show()
g.show(axes=True)
