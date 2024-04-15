""" See the note of mine but this doesn't simply scale to n-gon as Hyperbolic law of cosine for n-gon is research topic.
"""

from sage.plot.hyperbolic_polygon import HyperbolicPolygon


def i_rotation(z, alpha):
    _a = alpha / 2
    _c = cos(_a)
    _s = sin(_a)
    G = matrix([[_c, _s], [-_s, _c]])
    return (G[0][0] * z + G[0][1]) / (G[1][0] * z + G[1][1])


def generate_polygon(center, sides, i_angle):
    beta = 2 * pi / sides  # compute the rotation angle to be used ahead
    alpha = i_angle / Integer(2)
    I = CC(0, 1)
    # compute using cosine theorem the radius of the circumscribed circle
    # using the triangle formed by the radius and the three known angles
    r = arccosh(cot(alpha) * (1 + cos(beta)) / sin(beta))

    # The first point will be always on the imaginary axis limited
    # to 8 digits for efficiency in the subsequent calculations.
    z_0 = I * (e ** r).n(digits=8)

    # Compute the dilation isometry used to move the center
    # from I to the imaginary part of the given center.
    scale = center.imag()

    # Compute the parabolic isometry to move the center to the
    # real part of the given center.
    h_disp = center.real()

    d_z_k = [z_0 * scale + h_disp]  # d_k has the points for the polygon in the given center
    z_k = z_0  # z_k has the Re(z)>0 vertices for the I centered polygon
    r_z_k = []  # r_z_k has the Re(z)<0 vertices
    if is_odd(sides):
        vert = (sides - 1) // 2
    else:
        vert = sides // 2 - 1
    for k in range(vert):
        # Compute with 8 digits to accelerate calculations
        z_k = i_rotation(z_k, beta).n(digits=8)
        d_z_k += [z_k * scale + h_disp]
        r_z_k = [-(z_k).conjugate() * scale + h_disp] + r_z_k
        print(z_k, -(z_k).conjugate())

    if is_odd(sides):
        return HyperbolicPolygon(d_z_k + r_z_k, "UHP", options)
    else:
        z_opo = [I * (e ** (-r)).n(digits=8) * scale + h_disp]
        return HyperbolicPolygon(d_z_k + z_opo + r_z_k, "UHP", options)


options = {"alpha": 1, "fill": False}
sides = 3
i_angle = pi / 4
center = CC(0.0 + I)
polygon = generate_polygon(center=center, sides=sides, i_angle=i_angle)
print(polygon._pts)

from sage.plot.hyperbolic_regular_polygon import HyperbolicRegularPolygon

p = HyperbolicRegularPolygon(sides=sides, i_angle=i_angle, center=center, options=options)
print(p._pts)
