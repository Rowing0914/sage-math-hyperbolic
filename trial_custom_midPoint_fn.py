from sage.all import *
from sage.geometry.hyperbolic_space.hyperbolic_constants import EPSILON
from sage.geometry.hyperbolic_space.hyperbolic_geodesic import HyperbolicGeodesic
from sage.geometry.hyperbolic_space.hyperbolic_interface import HyperbolicPlane
from sage.matrix.matrix_symbolic_dense import Matrix_symbolic_dense
from sage.misc.lazy_import import lazy_import

lazy_import('sage.geometry.hyperbolic_space.hyperbolic_isometry',
            'moebius_transform')


# sage.rings.complex_mpfr.set_global_complex_round_mode(n=20)


def find_midpoint(line: HyperbolicGeodesic, n_iter: int = 1):
    UHP = line._model.realization_of().a_realization()
    l_UHP = line.to_model(model=UHP)

    if l_UHP.length() == infinity:
        raise ValueError("the length must be finite")

    start = l_UHP._start.coordinates()
    end = l_UHP._end.coordinates()
    d = l_UHP._model._dist_points(start, end) / 2
    S = l_UHP.complete()._to_std_geod(start)

    # If the matrix is symbolic then needs to be simplified in order to
    #    make the calculations easier for the symbolic calculus module.
    if isinstance(S, Matrix_symbolic_dense):
        S = S.simplify_full().simplify_full()
    S_1 = S.inverse()
    T = matrix([[exp(d), 0], [0, 1]])
    M = S_1 * T * S
    if ((real(start - end) < EPSILON) or
            (abs(real(start - end)) < EPSILON and
             imag(start - end) < EPSILON)):
        end_p = start
    else:
        end_p = end

    for i in range(n_iter):
        P_3 = moebius_transform(M, end_p)
        P = l_UHP._model.get_point(P_3)
        end_p = P_3
    return line._model(P)


PD = HyperbolicPlane().PD()
UHP = HyperbolicPlane().UHP()

pts = [
    (-0.4409894081700523 + 0.8975123073706717 * I, 0.7311976448002898 + 0.6821656721343501 * I),
    (0.9173343623719807 + 0.3981176554884141 * I, 0.39811765407683164 + 0.9173343629845988 * I),
    (0.8929547612469968 - 0.45014641436572456 * I, 0.6942481764060302 + 0.71973569423567 * I),
    (-0.948267839912836 + 0.317471422000539 * I, 0.275282020775022 + 0.961363515553831 * I),
    (-0.637265030365423 + 0.770644717800204 * I, -0.911931414471390 - 0.410342655959882 * I),
    (-0.271916134746048 + 0.962320952522997 * I, -0.962320953102423 + 0.271916132695439 * I),
    (-0.997981068784061 - 0.0635120960811705 * I, -0.0635120979596448 - 0.997981068664514 * I),
    (-0.953185397195065 + 0.302386505277743 * I, -0.402744794622099 - 0.915312313041185 * I),
    (-0.903077261621167 - 0.429478124638280 * I, 0.311836556591047 - 0.950135759759330 * I),
    (0.2902399360181517 - 0.9569539066957087 * I, 0.9635728586831255 + 0.26744596839217694 * I),
    (0.24535690078936162 - 0.969432819351109 * I, 0.9694328190091052 - 0.24535690214065609 * I),
    (-0.400019793494238 - 0.916506500147614 * I, 0.784270329432078 - 0.620419253708732 * I),
]

ll = list()
for p in pts:
    _s, _e = p
    s = PD.get_point(_s)
    e = PD.get_point(_e)
    print(s)
    s = s.to_model(model=UHP)
    print(s)
    asdf
    e = e.to_model(model=UHP)
    l1 = PD.get_geodesic(s, e)
    p_mid = l1.midpoint()
    print(p_mid)
    asdf
    for i in range(3):
        _l = PD.get_geodesic(s, p_mid)
        p_mid = _l.midpoint()
        print(s, p_mid)
    ll.append(l1)
    asdf
#
# pp = l1.midpoint()
# ppp = find_midpoint(line=l1)
# print(pp, ppp)
#
# mid_pts = list()
# for i in range(3):
#     ll = PD.get_geodesic(p1, pp)
#     pp = ll.midpoint()
#     mid_pts.append(pp)
#
# pppp = find_midpoint(line=l1, n_iter=3)

p = PD.get_point(pts[0][0])
print(p)
P = p.show(color="green")
for p in pts:
    s, e = p
    s = PD.get_point(s)
    e = PD.get_point(e)
    P += s.show(color="green")
    P += e.show(color="red")

P = ll[0].plot()
for l in ll:
    P += l.plot()
# P += pp.show()
# for pp in mid_pts:
#     P += pp.show()
# P += pppp.show(color="red")
# P += p2.show(color="green")
# P += l1.plot()
P.show(axes=True)
