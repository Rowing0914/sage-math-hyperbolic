n = 3
PD = HyperbolicPlane().PD()

points = list()
for i in range(n):
	points.append(0.5*cos((i+1)*2*pi/n) + 0.5*sin((i+1)*2*pi/n)*I)

# define sides
sides = list()
for i in range(n):
	if i+1 == n:
		_side = PD.get_geodesic(points[i], points[0])
	else:
		_side = PD.get_geodesic(points[i], points[i+1])
	sides.append(_side)

# Get 1st reflection transformations
reflection_1st = [l.reflection_involution() for l in sides]

# reflect the sides
reflect_1st_sides = list()
for i in range(n):
	for ii in range(n):
		print(i, ii)
		reflect_1st_sides.append(reflection_1st[i] * sides[ii])

# define a polygon
polygon = hyperbolic_polygon(points, model="PD", fill=True, color='purple')


@interact
def _(x=0.0, y=0.0):
	# base point
	p_base = PD.get_point(x + y*I)

	# reflect base point
	reflect_1st_pBase = [R*p_base for R in reflection_1st]

	# get 2nd reflection transformations
	reflection_2nd = [l.reflection_involution() for l in reflect_1st_sides]

	# reflect base point once more
	reflect_2nd_pBase = list()
	for R in reflection_2nd:
		for _point in reflect_1st_pBase:
			_p = R * _point
			reflect_2nd_pBase.append(_p)

	# reflect all sides once more
	reflect_1st_sides = list()
	for i in range(n):
		for ii in range(n):
			print(i, ii)
			reflect_1st_sides.append(reflection_1st[i] * sides[ii])
	R2R1_l1 = R2R1*R1_l1
	R2R1_l2 = R2R1*R1_l2
	R2R1_l3 = R2R1*R1_l3
	R2R1_l4 = R2R1*R1_l4
	R3R1_l1 = R3R1*R1_l1
	R3R1_l2 = R3R1*R1_l2
	R3R1_l3 = R3R1*R1_l3
	R3R1_l4 = R3R1*R1_l4
	R4R1_l1 = R4R1*R1_l1
	R4R1_l2 = R4R1*R1_l2
	R4R1_l3 = R4R1*R1_l3
	R4R1_l4 = R4R1*R1_l4

	R1R2_l1 = R1R2*R2_l1
	R1R2_l2 = R1R2*R2_l2
	R1R2_l3 = R1R2*R2_l3
	R1R2_l4 = R1R2*R2_l4
	R3R2_l1 = R3R2*R2_l1
	R3R2_l2 = R3R2*R2_l2
	R3R2_l3 = R3R2*R2_l3
	R3R2_l4 = R3R2*R2_l4
	R4R2_l1 = R4R2*R2_l1
	R4R2_l2 = R4R2*R2_l2
	R4R2_l3 = R4R2*R2_l3
	R4R2_l4 = R4R2*R2_l4

	R1R3_l1 = R1R3*R3_l1
	R1R3_l2 = R1R3*R3_l2
	R1R3_l3 = R1R3*R3_l3
	R1R3_l4 = R1R3*R3_l4
	R2R3_l1 = R2R3*R3_l1
	R2R3_l2 = R2R3*R3_l2
	R2R3_l3 = R2R3*R3_l3
	R2R3_l4 = R2R3*R3_l4
	R4R3_l1 = R4R3*R3_l1
	R4R3_l2 = R4R3*R3_l2
	R4R3_l3 = R4R3*R3_l3
	R4R3_l4 = R4R3*R3_l4

	R1R4_l1 = R1R4*R4_l1
	R1R4_l2 = R1R4*R4_l2
	R1R4_l3 = R1R4*R4_l3
	R1R4_l4 = R1R4*R4_l4
	R2R4_l1 = R2R4*R4_l1
	R2R4_l2 = R2R4*R4_l2
	R2R4_l3 = R2R4*R4_l3
	R2R4_l4 = R2R4*R4_l4
	R3R4_l1 = R3R4*R4_l1
	R3R4_l2 = R3R4*R4_l2
	R3R4_l3 = R3R4*R4_l3
	R3R4_l4 = R3R4*R4_l4

	# P-bisectors b/w 2nd reflections and base point
	l_pBase_R2R1base = PD.get_geodesic(p_base, R2R1_pBase).perpendicular_bisector().complete()
	l_pBase_R3R1base = PD.get_geodesic(p_base, R3R1_pBase).perpendicular_bisector().complete()
	l_pBase_R4R1base = PD.get_geodesic(p_base, R4R1_pBase).perpendicular_bisector().complete()
	l_pBase_R1R2base = PD.get_geodesic(p_base, R1R2_pBase).perpendicular_bisector().complete()
	l_pBase_R3R2base = PD.get_geodesic(p_base, R3R2_pBase).perpendicular_bisector().complete()
	l_pBase_R4R2base = PD.get_geodesic(p_base, R4R2_pBase).perpendicular_bisector().complete()
	l_pBase_R1R3base = PD.get_geodesic(p_base, R1R3_pBase).perpendicular_bisector().complete()
	l_pBase_R2R3base = PD.get_geodesic(p_base, R2R3_pBase).perpendicular_bisector().complete()
	l_pBase_R4R3base = PD.get_geodesic(p_base, R4R3_pBase).perpendicular_bisector().complete()
	l_pBase_R1R4base = PD.get_geodesic(p_base, R1R4_pBase).perpendicular_bisector().complete()
	l_pBase_R2R4base = PD.get_geodesic(p_base, R2R4_pBase).perpendicular_bisector().complete()
	l_pBase_R3R4base = PD.get_geodesic(p_base, R3R4_pBase).perpendicular_bisector().complete()

	# Plot sides
	P = l1.plot() + l2.plot() + l3.plot() + l4.plot()
	# plot base point
	P += p_base.show(size=50, color="blue", legend_label='p-base')
	# plot reflected sides
	P += R1_l1.plot() + R1_l2.plot() + R1_l3.plot() + R1_l4.plot()
	P += R2_l1.plot() + R2_l2.plot() + R2_l3.plot() + R2_l4.plot()
	P += R3_l1.plot() + R3_l2.plot() + R3_l3.plot() + R3_l4.plot()
	P += R4_l1.plot() + R4_l2.plot() + R4_l3.plot() + R4_l4.plot()

	# plot 1st reflection of points
	P += R1_pBase.show(color="red") + R2_pBase.show(color="red") + R3_pBase.show(color="red") + R4_pBase.show(color="red")
	# # plot 2nd reflection of sides
	# P += R2R1_l1.plot() + R2R1_l2.plot() + R2R1_l3.plot() + R2R1_l4.plot() + R3R1_l1.plot() + R3R1_l2.plot() + R3R1_l3.plot() + R3R1_l4.plot() + R4R1_l1.plot() + R4R1_l2.plot() + R4R1_l3.plot() + R4R1_l4.plot()
	# P += R1R2_l1.plot() + R1R2_l2.plot() + R1R2_l3.plot() + R1R2_l4.plot() + R3R2_l1.plot() + R3R2_l2.plot() + R3R2_l3.plot() + R3R2_l4.plot() + R4R2_l1.plot() + R4R2_l2.plot() + R4R2_l3.plot() + R4R2_l4.plot()
	# P += R1R3_l1.plot() + R1R3_l2.plot() + R1R3_l3.plot() + R1R3_l4.plot() + R2R3_l1.plot() + R2R3_l2.plot() + R2R3_l3.plot() + R2R3_l4.plot() + R4R3_l1.plot() + R4R3_l2.plot() + R4R3_l3.plot() + R4R3_l4.plot()
	# P += R1R4_l1.plot() + R1R4_l2.plot() + R1R4_l3.plot() + R1R4_l4.plot() + R2R4_l1.plot() + R2R4_l2.plot() + R2R4_l3.plot() + R2R4_l4.plot() + R3R4_l1.plot() + R3R4_l2.plot() + R3R4_l3.plot() + R3R4_l4.plot()
	# # # plot 2nd reflection of points
	# P += R2R1_pBase.show(color="red") + R3R1_pBase.show(color="red")
	# P += R1R2_pBase.show(color="red") + R3R2_pBase.show(color="red")
	# P += R1R3_pBase.show(color="red") + R2R3_pBase.show(color="red")
	# # plot p-bisectors for 2nd reflections
	P += l_pBase_R2R1base.plot(thickness=2, color="orange") + l_pBase_R3R1base.plot(thickness=2, color="orange") + l_pBase_R4R1base.plot(thickness=2, color="orange")
	P += l_pBase_R1R2base.plot(thickness=2, color="orange") + l_pBase_R3R2base.plot(thickness=2, color="orange") + l_pBase_R4R2base.plot(thickness=2, color="orange")
	P += l_pBase_R1R3base.plot(thickness=2, color="orange") + l_pBase_R2R3base.plot(thickness=2, color="orange") + l_pBase_R4R3base.plot(thickness=2, color="orange")
	P += l_pBase_R1R4base.plot(thickness=2, color="orange") + l_pBase_R2R4base.plot(thickness=2, color="orange") + l_pBase_R3R4base.plot(thickness=2, color="orange")

	P.show(axes=True)
