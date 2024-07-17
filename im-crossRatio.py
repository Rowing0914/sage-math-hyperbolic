""" Want to compute Im(Cross-Ratio) == 0
Vertex positions are taken from https://sagecell.sagemath.org/?q=chjfof
"""

# ==== SageMath
# points = [-0.202808201418345 - 0.351274108164389 * I, 0.405616401426393,
#           0.330298764996358 + 0.572094243649821 * I, -0.202808201418345 + 0.351274108164389 * I]
# points = [I * (1 + p) / (1 - p) for p in points]
# z1, z2, z3, z4 = points
#
# l = [(z1, z2), (z2, z3), (z3, z4), (z1, z4)]
# # l = [(z1, z4)]
# a_r = list()
# for z in l:
#     p, r = var('p, r')  # (p: centre, r: radius)
#     _z1, _z2 = z
#     print(z)
#     sols = solve([(_z1.real() - p) ^ 2 + (_z1.imag()) ^ 2 == r ^ 2,
#                   (_z2.real() - p) ^ 2 + (_z2.imag()) ^ 2 == r ^ 2], [p, r])
#
#     for sol in sols:
#         if_pass = False
#         for a in sol:
#             if a.lhs() == r and a.rhs().n() > 0.0:
#                 if_pass = True
#         if if_pass:
#             print(sol[0].rhs().n(), float(sol[1].rhs().n()))
#             a_r.append((sol[0].rhs().n(), float(sol[1].rhs().n())))  # (p, r)
#
# z = var("z", domain="complex")
# r_1_z = a_r[0][0] + (a_r[0][1] ** 2) / (conjugate(z) - a_r[0][0])
# r_2_z = a_r[1][0] + (a_r[1][1] ** 2) / (conjugate(z) - a_r[1][0])
# r_3_z = a_r[2][0] + (a_r[2][1] ** 2) / (conjugate(z) - a_r[2][0])
# r_4_z = a_r[3][0] + (a_r[3][1] ** 2) / (conjugate(z) - a_r[3][0])
# eqn = imag((r_1_z - r_3_z) / (r_1_z - r_4_z) * (r_2_z - r_4_z) / (r_2_z - r_3_z)) == 0.0
# # print(r_1_z, r_2_z, r_3_z, r_4_z, eqn)
# # print(expand(eqn))
# print(factor(eqn))
#
# -5.70920562610768, 6.179599638171012
# 0.979544827702941, 2.559673978533838
# -1.08022232145605, 0.826766370105546
# 0.000000000000000, 0.6952248772854238
#
# z = var("z", domain="complex")
# A = -5.70920562610768 + (6.179599638171012) ^ 2 / (conjugate(z) - (-5.70920562610768));
# B = 0.979544827702941 + (2.559673978533838) ^ 2 / (conjugate(z) - (0.979544827702941));
# C = -1.08022232145605 + (0.826766370105546) ^ 2 / (conjugate(z) - (-1.08022232145605));
# D = 0.0 + (0.6952248772854238) ^ 2 / (conjugate(z) - (0.0));
# factor(im( ( (A - C) / (A - D) ) / ( (B - C) / (B - D) ) ) == 0)
#
# eqn = imag(((A - C) / (A - D)) / ((B - C) / (B - D))) == 0.0
#
# print(factor(eqn))

""" Sympy: https://sagecell.sagemath.org/?q=tpvqan """
from sympy import *

z = symbols("z", complex=True)
A = -5.70920562610768 + (6.179599638171012) ** 2 / (conjugate(z) - (-5.70920562610768))
B = 0.979544827702941 + (2.559673978533838) ** 2 / (conjugate(z) - (0.979544827702941))
C = -1.08022232145605 + (0.826766370105546) ** 2 / (conjugate(z) - (-1.08022232145605))
D = 0.0 + (0.6952248772854238) ** 2 / (conjugate(z) - (0.0))

A = -5 + (6) ** 2 / (conjugate(z) - (-5))
B = 0 + (2) ** 2 / (conjugate(z) - (0))
C = -1 + (0) ** 2 / (conjugate(z) - (-1))
D = 0.0 + (0) ** 2 / (conjugate(z) - (0.0))

eqn = im(((A - C) / (A - D)) / ((B - C) / (B - D)))

print(expand(eqn))
