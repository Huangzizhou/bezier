### Bezier subdivision pre-computations ###

from itertools import product
from sympy import \
	symbols, prod, Matrix, MatrixSymbol, det, \
	collect, Rational, expand, pprint, numbered_symbols, cse
from sympy.printing.c import C99CodePrinter

## Custom printer for rationals ##
class MyCodePrinter(C99CodePrinter):
	def _print_Rational(self, expr):
		return f'R({expr.p}, {expr.q})'

## C print ##
def C_print(expr, n, s, p):
	CSE_results = cse(expr, numbered_symbols('tmp_'), optimizations='basic')
	lines = [
		'template<>\n' +
		f'void lagrangeVectorT<{n}, {s}, {p}>' +
		'(const std::vector<fp_t> &cpFP, std::vector<Interval> &out) {'
	]
	L = len(CSE_results[1])
	lines.append(f'out.resize({L});')
	lines.append(f'std::vector<Interval> cp(cpFP.size());')
	lines.append('for (uint i = 0; i < cpFP.size(); ++i) cp[i] = cpFP[i];')
	my_ccode = MyCodePrinter().doprint
	for helper in CSE_results[0]:
		lines.append('const Interval ' + my_ccode(helper[1], helper[0]))
	for i,result in enumerate(CSE_results[1]):
		lines.append(my_ccode(result, f'out[{i}]'))
	return '\n\t'.join(lines) + '\n}'

## Subscripting ##
def subscripts(b, i):
	return symbols([f'{b}[{j}]' for j in i])

## Index set ##
def index_set(n, s, p):
	return [tup for tup in product(range(p+1), repeat=n)
		if sum(tup[:s]) <= p and all(x <= p for x in tup[s:])]

## Corner indices set ##
def corner_indices_set(n, s, p):
	simplex_ord = n * p - s
	tensor_ord = n * p - 1
	assert(simplex_ord >= 0 and tensor_ord >= 0)
	time_ord = n
	en = enumerate(index_set_J(n,s,p))
	return [ind for ind,tup in en
		if all(x == simplex_ord or x == 0 for x in tup[:s]) \
		and all(x == tensor_ord or x == 0 for x in tup[s:n]) \
		and tup[-1] == time_ord]

## Univariate Lagrange polynomial ##
def lagrange_uni(p, i, x):
	k_values = [k for k in range(0, p+1) if k != i]
	return prod((p * x - k) / (i - k) for k in k_values)

## Multivariate Lagrange polynomial ##
def lagrange(n, s, p, i, x_name):
	assert(len(i) == n)
	x = subscripts(x_name, range(n))
	i_slack = p - sum(i[:s])
	assert(i_slack >= 0)
	x_slack = 1 - sum(x[:s])
	return \
		lagrange_uni(i_slack, i_slack, x_slack) * \
		prod([lagrange_uni(i[k], i[k], x[k]) for k in range(s)]) * \
		prod([lagrange_uni(p, i[k], x[k]) for k in range(s,n)])

## Geometric map component ##
def geo_map_component(n, s, p, t, x_name, c_name, ni):
	indices = index_set(n, s, p)
	cp = subscripts(c_name, range(t + 2*ni, len(indices)*2*n, 2*n))
	lag = [lagrange(n, s, p, k, x_name) for k in indices]
	return sum(cp[k] * lag[k] for k in range(len(indices)))
	
## Geometric map ##
def geo_map(n, s, p, x_name, c_name, t_name):
	T = symbols(t_name)
	return [
		(1 - T) * geo_map_component(n, s, p, 0, x_name, c_name, i) + \
		T * geo_map_component(n, s, p, 1, x_name, c_name, i)
		for i in range(n)] + [T]

## Jacobian determinant ##
def jac_det(n, s, p, x_name, c_name, t_name):
	xt = subscripts(x_name, range(n)) + [symbols(t_name)]
	gmap = geo_map(n, s, p, x_name, c_name, t_name)
	res = det(Matrix([[poly.diff(v) for v in xt] for poly in gmap]))
	res = collect(res, symbols(t_name))
	return res

## Domain points for the Jacobian determinant ##
def index_set_J(n, s, p):
	simplex_ord = n * p - s
	tensor_ord = n * p - 1
	assert(simplex_ord >= 0 and tensor_ord >= 0)
	time_ord = n
	return [tup for tup in product(range(n*p+1), repeat=n+1)
		if sum(tup[:s]) <= simplex_ord \
		and all(x <= tensor_ord for x in tup[s:n]) \
		and tup[-1] <= time_ord]

def domain_pts_J(n, s, p):
	simplex_den = max(1, n * p - s)
	tensor_den = max(1, n * p - 1)
	time_den = n
	return [
		tuple(
			Rational(elem, simplex_den) if i < s else
			Rational(elem, time_den) if i == n else
			Rational(elem, tensor_den)
			for i, elem in enumerate(tup)
		)
		for tup in index_set_J(n, s, p)]

## Lagrange vector ##
def lagrange_vector(n, s, p, c_name):
	x_name = 'u'
	t_name = 't'
	poly = jac_det(n, s, p, x_name, c_name, t_name)
	pts = domain_pts_J(n, s, p)
	xt = subscripts(x_name, range(n)) + [symbols(t_name)]
	rule = lambda pt: {xt[k]: pt[k] for k in range(n+1)}
	lag_vec = [poly.subs(rule(pt)) for pt in pts]
	return lag_vec


## Time subdivision maps ##
def time_subdiv_map(pt, shift):
	h = 1 if shift else 0
	return(pt[:-1] + ((pt[-1] + h) / 2,))

## Space subdivision maps ##
def space_subdiv_map(pt, n, s, q):
	assert(1 <= s <= n)
	assert(len(pt) == n+1)
	assert(0 <= q < 2**(n+1))
	res = list(pt)
	qs = q % (2**s)
	# Simplex part
	if s == 1:
		res[0] = pt[0] + qs
	elif s == 2:
		if qs == 0:
			res[0] = pt[0]
			res[1] = pt[1]
		elif qs == 1:
			res[0] = pt[0] + 1
			res[1] = pt[1]
		elif qs == 2:
			res[0] = pt[0]
			res[1] = pt[1] + 1
		elif qs == 3:
			res[0] = pt[0] + 1
			res[1] = pt[1] + 1
	elif s == 3:
		if qs == 0:
			res[0] = pt[0]
			res[1] = pt[1]
			res[2] = pt[2]
		elif qs == 1:
			res[0] = pt[0] + 1
			res[1] = pt[1]
			res[2] = pt[2]
		elif qs == 2:
			res[0] = pt[0]
			res[1] = pt[1] + 1
			res[2] = pt[2]
		elif qs == 3:
			res[0] = pt[0]
			res[1] = pt[1]
			res[2] = pt[2] + 1
		elif qs == 4:
			res[0] = 1 - pt[1] - pt[2]
			res[1] = pt[1]
			res[2] = pt[0] + pt[1] + pt[2]
		elif qs == 5:
			res[0] = 1 - pt[1]
			res[1] = pt[0] + pt[1]
			res[2] = pt[1] + pt[2]
		elif qs == 6:
			res[0] = pt[0] + pt[1]
			res[1] = 1 - pt[0]
			res[2] = pt[2]
		elif qs == 7:
			res[0] = pt[0]
			res[1] = pt[1] + pt[2]
			res[2] = 1 - pt[0] - pt[1]
	else:
		raise Exception('Unsupported simplex dimension')

	# Tensor part
	for k in range(s,n+1):
		qt = (q >> k) & 1
		res[k] = pt[k] + qt

	# Divide and return
	return tuple(r/2 for r in res)

## Compress matrix ##
def mat_compress(m):
	l = m.rows
	assert(m.cols == l)
	s = str(l)
	for i in range(l):
		for j in range(l):
			if m[i, j] != 0:
				s += f', {i}, {j}, {m[i,j].p}, {m[i,j].q}'
	return s

## Binomial and multinomial coefficients ##
def multinomial(ii):
	res = 1
	k = sum(ii)
	i0 = ii.index(max(ii))
	for a in ii[:i0] + ii[i0+1:]:
		for j in range(1,a+1):
			res *= k
			res //= j
			k -= 1
	return res

## Subdivision matrices ##
def subdiv_matrices(n, s, p):
	simplex_ord = n * p - s
	tensor_ord = n * p - 1
	time_ord = n
	xt = subscripts('u', range(n)) + [symbols('t')]
	xt_full = xt + [1 - sum(xt[:s])] + [1 - xt[k] for k in range(s, n+1)]
	exponents = index_set_J(n, s, p)
	exponents_full = [list(e) + [simplex_ord - sum(e[:s])] + \
		[tensor_ord - e[k] for k in range(s, n)] + [time_ord - e[n]]
		for e in exponents]
	r_full = range(len(xt_full))
	basis = [
		multinomial(list(e[:s]) + [simplex_ord - sum(e[:s])]) *
		prod(multinomial((ee, tensor_ord - ee)) for ee in e[s:n]) *
		multinomial((e[n], time_ord - e[n])) *
		prod([xt_full[k]**exponents_full[i][k] for k in r_full])
		for i,e in enumerate(exponents)]
	assert(expand(sum(basis)) == 1)
	if (__debug__): print('Basis computed')
	pts = domain_pts_J(n, s, p)

	rule = lambda e: {xt[k]: e[k] for k in range(n+1)}
	b2l = Matrix([[b.subs(rule(pt)) for b in basis] for pt in pts])
	if (__debug__): print('B2L computed')
	l2b = b2l.inv()
	if (__debug__): print('L2B computed')
	tsd = [
		l2b * Matrix([[
			b.subs(rule(time_subdiv_map(pt, t)))
			for b in basis] for pt in pts])
		for t in [False, True]]
	if (__debug__): print('Time subdivision matrices computed')
	ssd = [
		l2b * Matrix([[
			b.subs(rule(space_subdiv_map(pt, n, s, q)))
			for b in basis] for pt in pts])
		for q in range(2**(n+1))]
	if (__debug__): print('Space subdivision matrices computed')
	return [l2b, tsd, ssd]

## Format matrices ##
def matrices_formatted(n, s, p):
	lines = [
		'template<>\n' +
		f'void initMatricesT<{n}, {s}, {p}> ' +
		'(Matrix<Interval> &l2b, ' +
		'std::pair<Matrix<Interval>, Matrix<Interval>> &tsd, ' +
		f'std::array<Matrix<Interval>, {2<<n}> &ssd) ' +
		'{'
	]
	
	mm = subdiv_matrices(n, s, p)
	lines.append('l2b.fill({' + mat_compress(mm[0]) + '});')
	lines.append('tsd.first.fill({' + mat_compress(mm[1][0]) + '});')
	lines.append('tsd.second.fill({' + mat_compress(mm[1][1]) + '});')
	for q,a in enumerate(mm[2]):
		lines.append(f'ssd[{q}].fill(' + '{' + mat_compress(a) + '});')
	return '\n\t'.join(lines) + '\n}'

## Format lag vector ##
def lag_vec_formatted(n, s, p):
	return C_print(lagrange_vector(n,s,p,'cp'), n, s, p)

## Format corner indices
def corners_formatted(n, s, p):
	ind = map(str,corner_indices_set(n, s, p))
	return f'template<>\nvoid cornerIndicesT<{n}, {s}, {p}>(std::vector<uint> &v) ' + \
		'{ v = {' + ','.join(ind) + '}; }'

## Write ##
combinations = [
	(1,1,1),
	(1,1,2),
	(1,1,3),
	(1,1,4),
	(1,1,5),
	(2,2,1),
	(2,2,2),
	(2,2,3),
	(3,3,1),
	(3,3,2),
]

WRITE_MATRICES = False
WRITE_LAGVEC = False
WRITE_CORNERS = False
INFO_ORDER = False
INFO_JAC_ORDER = True

if WRITE_MATRICES:
	print('Writing matrices...')
	with open('src/validity/transMatrices.cpp', 'w') as f:
		f.write('#include "transMatrices.hpp"\n\n')
		f.write('namespace element_validity {\n')

		for n,s,p in combinations:
			f.write(matrices_formatted(n,s,p))
			f.write('\n\n')

		f.write('\n}\n')
	print('Done writing matrices.')
else:
	print('Not writing matrices.')

if WRITE_LAGVEC:
	print('Writing Lagrange vectors...')
	with open('src/validity/lagrangeVector.cpp', 'w') as f:
		f.write('#include "lagrangeVector.hpp"\n\n')
		f.write('#define R(p, q) Interval(p) / q\n\n')
		f.write('namespace element_validity {\n')

		for n,s,p in combinations:
			f.write(lag_vec_formatted(n,s,p))
			f.write('\n\n')

		f.write('}\n')
	print('Done writing Lagrange vectors.')
else:
	print('Not writing Lagrange vectors.')

if WRITE_CORNERS:
	print('Writing corner indices...')
	with open('src/validity/cornerIndices.cpp', 'w') as f:
		f.write('#include "cornerIndices.hpp"\n\n')
		f.write('namespace element_validity {\n')

		for n,s,p in combinations:
			f.write(corners_formatted(n,s,p))
			f.write('\n\n')

		f.write('\n}\n')
	print('Done writing corner indices.')
else:
	print('Not writing corner indices.')

if INFO_ORDER:
	for n,s,p in combinations:
		print(f'This is the order of control points for elements of type ({n},{s},{p}):')
		ind = index_set(n, s, p)
		for i,t in enumerate(ind):
			print(f'\t{i}: {t}')
		print()

if INFO_JAC_ORDER:
	for n,s,p in combinations:
		print(f'This is the order of jacobian control points for elements of type ({n},{s},{p}):')
		ind = index_set_J(n, s, p)
		for i,t in enumerate(ind):
			print(f'\t{i}: {t}')
		print()