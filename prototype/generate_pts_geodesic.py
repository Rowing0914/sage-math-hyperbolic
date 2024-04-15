def rotation(z, alpha, c=0 + I * 0):
    _a = alpha / 2
    _c = cos(_a)
    _s = sin(_a)
    G = matrix([[_c, _s], [-_s, _c]])
    translated_z = z - c
    rotated_z = (G[0][0] * translated_z + G[0][1]) / (G[1][0] * translated_z + G[1][1])
    return rotated_z + c


UHP = HyperbolicPlane().UHP()

g = Graphics()

z = zz = 1 + I
p = UHP.get_point(z)
g += p.show(color="red")

z = rotation(z=z, alpha=2 * pi / 3)
p1 = UHP.get_point(z)
g += p1.show(color="green")

z = rotation(z=z, alpha=2 * pi / 3)
p1 = UHP.get_point(z)
g += p1.show(color="green")

z = rotation(z=z, alpha=2 * pi / 3)
p1 = UHP.get_point(z)
g += p1.show(color="green")

z = zz
z = rotation(z=z, alpha=2 * pi / 3, c=z)
p1 = UHP.get_point(z)
g += p1.show(color="blue")

z = rotation(z=z, alpha=2 * pi / 3, c=z)
p1 = UHP.get_point(z)
g += p1.show(color="blue")

z = rotation(z=z, alpha=2 * pi / 3, c=z)
p1 = UHP.get_point(z)
g += p1.show(color="blue")

g.show(axes=True)
