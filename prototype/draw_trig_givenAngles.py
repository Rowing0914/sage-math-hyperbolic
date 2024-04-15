import numpy as np

UHP = HyperbolicPlane().UHP()  # UHP is the upper half plane IH
HM = HyperbolicPlane().HM()  # HM  is the hyperboloid model
PD = HyperbolicPlane().PD()

n1, n2, n3 = 4, 4, 8
a1, a2, a3 = pi / n1, pi / n2, pi / n3  # given angles, we draw a hyperbolic triangle with these angles
N = 2 * lcm([n1, n2, n3])
F.<u> = CyclotomicField(N)


def c(a1, a2, a3):
    return (cos(a1) + cos(a2) * cos(a3)) / sin(a2) / sin(a3)  # !!! sin(a3) is in denominator !!!


c1, c2, c3 = c(a1, a2, a3), c(a2, a3, a1), c(a3, a1, a2)  # algebraic in the given example


# sides are CC(arccosh(c1)), CC(arccosh(c2)), CC(arccosh(c3))

def toUHP(M):
    """M is a point in the HM model, we construct an ad-hoc version in UHP
    """
    a, b, c = M.coordinates()
    a, b, c = QQbar(a), QQbar(b), QQbar(c)
    return UHP(-(a + i) * (c + 1) / ((b - 1) * c - a ^ 2 - b ^ 2 + b - 1))


def s(v, w):
    """v, w are vectors with three entries, we return the Minkowski product with signature ++-"""
    return v * diagonal_matrix([1, 1, -1]) * w  # no need of a transposition


# a, b, p, q, r are used in the following coordinates in Hyperboloid model
# as a parametrization of points
R.<a,b,p,q,r> = PolynomialRing(F)

V1 = vector([0, 0, 1])
V2 = vector([0, a, b])
V3 = vector([p, q, r])

J = R.ideal([s(V2, V2) + 1, s(V3, V3) + 1,
             s(V1, V2) + c3, s(V2, V3) + c1, s(V3, V1) + c2])
V = J.variety(ring=QQbar)
sols = [sol for sol in V if sol[a] > 0 and sol[q] > 0]  # so V2, V3 maps to IH
sol = sols[0]  # first solution

a0, b0, p0, q0, r0 = [sol[v].real() for v in R.gens()]

S1, S2, S3 = vector([0, 0, 1]), vector([0, a0, b0]), vector([p0, q0, r0])
M1, M2, M3 = HM.get_point(S1), HM.get_point(S2), HM.get_point(S3)
H1, H2, H3 = toUHP(M1), toUHP(M2), toUHP(M3)  # using the ad-hoc coercion from HM to UHP
P1, P2, P3 = PD(H1), PD(H2), PD(H3)
Q1, Q2, Q3 = H1.coordinates(), H2.coordinates(), H3.coordinates()
p1, p2, p3 = P1.coordinates(), P2.coordinates(), P3.coordinates()

l1 = UHP.get_geodesic(H1, H2)
l2 = UHP.get_geodesic(H2, H3)
l3 = UHP.get_geodesic(H3, H1)

_a1, _a2, _a3 = np.rad2deg(float(l1.angle(l2))), np.rad2deg(float(l2.angle(l3))), np.rad2deg(float(l3.angle(l1)))
print("Angles: ", round(_a1, ndigits=2), round(_a2, ndigits=2), round(_a3, ndigits=2))

# p = hyperbolic_polygon(pts=[Q1, Q2, Q3], model="UHP", fill=True, alpha=0.3)
p = hyperbolic_polygon(pts=[p1, p2, p3], model="PD", fill=True, alpha=0.3)

g = Graphics()
g += p.plot()
g.show(axes=True, aspect_ratio=1)
