PD = HyperbolicPlane().PD()

@interact
def _(num_sides=3, x=0.0, y=0.0, auto_update=False,
	if_plot_sides=False, if_plot_reflect_1st_sides=False, if_plot_reflect_1st_pBase=False,
	if_plot_reflect_2nd_sides=False, if_plot_reflect_2nd_pBase=False, if_plot_perp_bisec=False,
	):

	n = int(num_sides)

	points = list()
	for i in range(n):
		_x, _y = 0.5*cos((i+1)*2*pi/n), 0.5*sin((i+1)*2*pi/n)
		_x, _y = round(_x, 2), round(_y, 2)
		points.append(_x + _y * I)

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

	# reflect the sides; fix R and iterate sides; [R1-l1, R1-l2, R1-l3, ...]
	reflect_1st_sides = list()  # n^2
	for R in reflection_1st:  # n
		for s in sides:  # n
			reflect_1st_sides.append(R * s)

	# base point
	x, y = round(x, 2), round(y, 2)
	p_base = PD.get_point(x + y*I)

	# reflect base point: n
	reflect_1st_pBase = [R*p_base for R in reflection_1st]

	# get 2nd reflection transformations: n^2
	reflection_2nd = [l.reflection_involution() for l in reflect_1st_sides]

	# reflect all sides once more
	reflect_2nd_sides = list()
	reflect_2nd_pBase = list()
	# reflect_1st_sides: [R1-l1, R1-l2, R1-l3, ...] (n^2)
	# reflection_2nd: [R1R1, R2R1, R3R1, ...] (n^2)
	
	# For Second reflection, eg., R_{i,j}, we compute Diff in i and j to differentiate colours in plot later.
	diff_index = list()
	for h in range(n):
		for ind_i, i in enumerate(range(h*n, (h+1)*n)):  # reflection transf
			for j in range(h*n, (h+1)*n):  # side
				R = reflection_2nd[i]
				if i == j and ind_i != h:  # ind_i != h to avoid reflecting back to p-base
					# Transform point
					_p = reflect_1st_pBase[h]
					_p = R * _p
					reflect_2nd_pBase.append(_p)
					diff_index.append(ind_i)  # 'i' already represents the diff
				else:
					# Transform sides
					_s = reflect_1st_sides[j]
					reflect_2nd_sides.append(R * _s)

	# P-bisectors b/w 2nd reflections and base point
	list_perp_bisec = list()
	for p in reflect_2nd_pBase:
		_l = PD.get_geodesic(p_base, p).perpendicular_bisector().complete()
		list_perp_bisec.append(_l)

	# plot base point
	P = p_base.show(size=50, color="blue", legend_label='p-base')

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
		color_keys = list(colors.keys())
		used_colors = dict()
		dict_if_plot_done = {k : False for k in set(diff_index)}
		for i, _diff in zip(list_perp_bisec, diff_index):
			# _c = cmap[_diff]
			_c = colors[color_keys[_diff]]
			P += i.plot(thickness=1.5, color=_c)
			if _diff not in used_colors:
				used_colors[_diff] = color_keys[_diff]
		res = list()
		for k, v in used_colors.items():
			_label = f"R_i_i+{k}: {v}"
			res.append(_label)
		print("=== Colour Legend ===")
		print(res)

	P.show(axes=True)
