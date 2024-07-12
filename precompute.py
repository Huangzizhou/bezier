### Bezier subdivision pre-computations ###
from polyfem_indices import index_set_polyfem

## Setup ##
COMBINATIONS = [
	(1,1,1), (1,1,2), (1,1,3), (1,1,4), (1,1,5), # segments
	(2,1,1), (2,1,2), # quadrilaterals
	(2,2,1), (2,2,2), (2,2,3), (2,2,4), # triangles
	(3,1,1), (3,1,2), # hexahedra
	(3,3,1), (3,3,2), (3,3,3), (3,3,4), # tetrahedra
]

OPTIONAL_OPTIMIZE = [(3,3,4),]

DRY_RUN = False
POLYFEM_ORDER = True

WRITE_CMAKE = True
WRITE_MATRICES = False
WRITE_LAGVEC = False
WRITE_CORNERS = False

INFO_ORDER = False
INFO_JAC_ORDER = False
INFO_LAGBASIS = False
INFO_LAGVECTOR = False

TO_ORIGIN = False

## Imports

from itertools import product
from sympy import \
	symbols, prod, Matrix, MatrixSymbol, det, Poly, \
	collect, collect_const, expand, \
	pprint, numbered_symbols, cse, QQ, ZZ
from sympy.simplify.ratsimp import ratsimp
from sympy.printing.cxx import CXX11CodePrinter
from sympy.polys.matrices import DomainMatrix
from concurrent.futures import ProcessPoolExecutor

## Custom printer for rationals ##
class MyCodePrinter(CXX11CodePrinter):
	def _print_Rational(self, expr):
		return f'R({expr.p}, {expr.q})'

## C print ##
def C_print(expr, n, s, p):
	if (__debug__): print('Simplification/CSE start')
	CSE_results = cse(expr, numbered_symbols('tmp_'), order='canonical')
	if (__debug__): print('Simplification/CSE end')
	lines = [
		'template<>\n' +
		f'void lagrangeVectorT<{n}, {s}, {p}>' +
		'(const std::span<const fp_t> cpFP, const std::span<Interval> out) {'
	]
	# lines.append(f'out.resize({len(expr)});')
	lines.append(f'assert(out.size() == {len(expr)});')
	assert len(expr) == len(CSE_results[1])
	lines.append('const uint S = cpFP.size();')
	lines.append('std::vector<Interval> cp(S);')
	lines.append('for (uint i = 0; i < S; ++i) cp[i] = cpFP[i];')
	if TO_ORIGIN:
		lines.append(f'for (int i = S-1; i >= 0; --i) cp[i] -= cp[i%{n}];')
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
	assert(0 < s <= n)
	if POLYFEM_ORDER: return index_set_polyfem(n,s,p)
	else: return [tup for tup in product(range(p+1), repeat=n)
		if sum(tup[:s]) <= p and all(x <= p for x in tup[s:])]

## Simplification
def my_simplify(expr):
	return expand(expr)

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
def lagrange_uni(m, p, i, x):
	k_values = [k for k in range(0, p+1) if k != i]
	return prod((m * x - k) / (i - k) for k in k_values)

## Multivariate Lagrange polynomial ##
def lagrange(n, s, p, i, x_name):
	assert(len(i) == n)
	x = subscripts(x_name, range(n))
	i_slack = p - sum(i[:s])
	assert(i_slack >= 0)
	x_slack = 1 - sum(x[:s])
	return my_simplify(
		lagrange_uni(p, i_slack, i_slack, x_slack) * \
		prod([lagrange_uni(p, i[k], i[k], x[k]) for k in range(s)]) * \
		prod([lagrange_uni(p, p, i[k], x[k]) for k in range(s,n)])
	)

## Geometric map component ##
def geo_map_component(n, s, p, t, x_name, c_name, ni):
	if (__debug__): print(f"Geometric map ({n} {s} {p}), time {t}, component {ni}")
	indices = index_set(n, s, p)
	cp = subscripts(c_name, range(t + 2*ni, len(indices)*2*n, 2*n))
	if TO_ORIGIN: cp[0] = 0
	lag = [lagrange(n, s, p, k, x_name) for k in indices]
	return my_simplify(sum(cp[k] * lag[k] for k in range(len(indices))))
	
## Geometric map ##
def aux_map_component(args):
	n, s, p, t, x_name, c_name, i = args
	return geo_map_component(n, s, p, t, x_name, c_name, i)

def aux_combine_results(args):
	i, T, results_0, results_1 = args
	return (1 - T) * results_0[i] + T * results_1[i]

def geo_map(n, s, p, x_name, c_name, t_name):
	T = symbols(t_name)

	with ProcessPoolExecutor() as executor:
		args_0 = [(n, s, p, 0, x_name, c_name, i) for i in range(n)]
		args_1 = [(n, s, p, 1, x_name, c_name, i) for i in range(n)]

		# Submit tasks for t=0 and t=1 components
		futures_0 = [executor.submit(aux_map_component, arg) for arg in args_0]
		futures_1 = [executor.submit(aux_map_component, arg) for arg in args_1]

		results_0 = [future.result() for future in futures_0]
		results_1 = [future.result() for future in futures_1]

	with ProcessPoolExecutor() as executor:
		# Prepare arguments for combining results
		args_combined = [(i, T, results_0, results_1) for i in range(n)]

		futures_combined = [executor.submit(aux_combine_results, arg)
			for arg in args_combined]
		results_combined = [future.result() for future in futures_combined]

	return results_combined + [T]


## Jacobian determinant ##
def jacobian(n, s, p, x_name, c_name, t_name):
	xt = subscripts(x_name, range(n)) + [symbols(t_name)]
	gmap = geo_map(n, s, p, x_name, c_name, t_name)
	J =	Matrix([[my_simplify(poly.diff(v)) for v in xt] for poly in gmap])
	return J

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
			QQ(elem, simplex_den) if i < s else
			QQ(elem, time_den) if i == n else
			QQ(elem, tensor_den)
			for i, elem in enumerate(tup)
		)
		for tup in index_set_J(n, s, p)]

## Lagrange vector ##
def aux_evaluate_det(args):
	pt, J, xt, n = args
	if (__debug__): print(f"\t Det eval {pt} start")
	d = J.subs({xt[k]: pt[k] for k in range(n+1)}).det(method='berkowitz')
	if (__debug__): print(f"\t Det eval {pt} end")
	return d

def lagrange_vector(n, s, p, c_name):
	x_name = 'u'
	t_name = 't'
	if (__debug__): print("Jac start")
	J = jacobian(n, s, p, x_name, c_name, t_name)
	if (__debug__): print("Jac end")
	pts = domain_pts_J(n, s, p)
	xt = subscripts(x_name, range(n)) + [symbols(t_name)]

	with ProcessPoolExecutor() as executor:
		futures = [
			executor.submit(aux_evaluate_det, [pt, J, xt, n]) for pt in pts
		]
		lag_vec = [future.result() for future in futures]

	return lag_vec



## Time subdivision maps ##
def time_subdiv_map(pt, shift):
	h = 1 if shift else 0
	return(pt[:-1] + ((pt[-1] + h) * QQ(1,2),))

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
			res[0] = 1 - pt[0]
			res[1] = 1 - pt[1]
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
	return tuple(r * QQ(1,2) for r in res)

def time_subdiv_maps(pts, shift):
	return [time_subdiv_map(pt, shift) for pt in pts]
def space_subdiv_maps(pts, n, s, q):
	return [space_subdiv_map(pt, n, s, q) for pt in pts]

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
def eval_basis(args):
	xt, basis, pts = args
	return DomainMatrix.from_Matrix(Matrix(
	[[b.eval({xt[k]: pt[k] for k in range(len(xt))})
	for b in basis] for pt in pts])).to_field()
def aux_matmul(args): return args[0] * args[1]

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
	basis = [Poly(
			multinomial(list(e[:s]) + [simplex_ord - sum(e[:s])]) *
			prod(multinomial((ee, tensor_ord - ee)) for ee in e[s:n]) *
			multinomial((e[n], time_ord - e[n])) *
			prod([xt_full[k]**exponents_full[i][k] for k in r_full]),
		xt, domain='QQ') for i,e in enumerate(exponents)]
	assert(sum(basis) == 1)
	if (__debug__): print('Basis computed')
	pts = domain_pts_J(n, s, p)

	with ProcessPoolExecutor() as executor:
		args_combined = [(xt, basis, pts)] + \
			[(xt, basis, time_subdiv_maps(pts, t)) for t in (False, True)] + \
			[(xt, basis, space_subdiv_maps(pts, n, s, q)) for q in range(2**(n+1))]

		futures_combined = [executor.submit(eval_basis, args)
			for args in args_combined]
		results_combined = [future.result() for future in futures_combined]
	if (__debug__): print('Prepared subdivision matrices')

	l2b = results_combined[0].inv()
	tsd = [l2b * m for m in results_combined[1:3]]
	ssd = [l2b * m for m in results_combined[3:]]
	if (__debug__): print('Computed subdivision matrices')
	return [
		l2b.to_Matrix(),
		[m.to_Matrix() for m in tsd],
		[m.to_Matrix() for m in ssd]
	]

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
if WRITE_CMAKE and not DRY_RUN:
	with open('src/validity/CMakeLists.txt', 'w') as f:
		f.write('set(HEADERS\n')
		f.write('\telement_validity.hpp\n')
		f.write('\tcornerIndices.hpp\n')
		f.write('\tlagrangeVector.hpp\n')
		f.write('\ttransMatrices.hpp\n')
		f.write(')\n\n')

		f.write('set(SOURCES\n')
		for n,s,p in COMBINATIONS:
			f.write(f'\ttransMatrices_{n}_{s}_{p}.cpp\n')
			f.write(f'\tlagrangeVector_{n}_{s}_{p}.cpp\n')
		f.write('\tcornerIndices.cpp\n')
		f.write(')\n\n')

		f.write('source_group(TREE "${CMAKE_CURRENT_SOURCE_DIR}" PREFIX "Source Files" FILES ${SOURCES})\n')
		f.write('target_sources(bezier PRIVATE ${SOURCES})')

if WRITE_MATRICES:
	print('Writing matrices...')
	path = lambda n,s,p: f'src/validity/transMatrices_{n}_{s}_{p}.cpp'
	if DRY_RUN: path = lambda n,s,p: '/dev/null'
	for n,s,p in COMBINATIONS:
		with open(path(n,s,p), 'w') as f:
			f.write('#include "transMatrices.hpp"\n\n')
			f.write('namespace element_validity {\n')
			if (__debug__): print(f'({n} {s} {p})')
			f.write(matrices_formatted(n,s,p))
			f.write('\n}\n')
	print()
	print('Done writing matrices.')
else:
	print('Not writing matrices.')

if WRITE_LAGVEC:
	print('Writing Lagrange vectors...')
	path = lambda n,s,p: f'src/validity/lagrangeVector_{n}_{s}_{p}.cpp'
	if DRY_RUN: path = lambda n,s,p: '/dev/null'
	for n,s,p in COMBINATIONS:
		optopt = (n,s,p) in OPTIONAL_OPTIMIZE
		with open(path(n,s,p), 'w') as f:
			f.write('#include "lagrangeVector.hpp"\n\n')
			f.write('#define R(p, q) (Interval(p) / q)\n\n')
			f.write('namespace element_validity {\n')
			if optopt:
				f.write('#ifdef LAGVEC_GCC_O0\n')
				f.write('#pragma GCC push_options\n')
				f.write('#pragma GCC optimize ("-O0")\n')
				f.write('#endif\n')
			f.write(lag_vec_formatted(n,s,p))
			f.write('}\n')
			if optopt:
				f.write('#ifdef LAGVEC_GCC_O0\n')
				f.write('#pragma GCC pop_options\n')
				f.write('#endif\n')
			f.write('#undef R')
	print('Done writing Lagrange vectors.')
else:
	print('Not writing Lagrange vectors.')

if WRITE_CORNERS:
	print('Writing corner indices...')
	path = 'src/validity/cornerIndices.cpp'
	if DRY_RUN: path = '/dev/null'
	with open(path, 'w') as f:
		f.write('#include "cornerIndices.hpp"\n\n')
		f.write('namespace element_validity {\n')

		for n,s,p in COMBINATIONS:
			f.write(corners_formatted(n,s,p))
			f.write('\n\n')

		f.write('\n}\n')
	print('Done writing corner indices.')
else:
	print('Not writing corner indices.')

if INFO_ORDER:
	for n,s,p in COMBINATIONS:
		print(f'This is the order of geometric map control points for elements of type ({n},{s},{p}):')
		ind = index_set(n, s, p)
		for i,t in enumerate(ind):
			print(f'\t{i}: {t}')
		print()

if INFO_JAC_ORDER:
	for n,s,p in COMBINATIONS:
		print(f'This is the order of jacobian control points for elements of type ({n},{s},{p}):')
		ind = index_set_J(n, s, p)
		for i,t in enumerate(ind):
			print(f'\t{i}: {t}')
		print()

if INFO_LAGBASIS:
	for n,s,p in COMBINATIONS:
		print(f'This is the Lagrange basis for elements of type ({n},{s},{p}):')
		indices = index_set(n, s, p)
		lag = [lagrange(n, s, p, k, 'x') for k in indices]
		for i,poly in enumerate(lag):
			print(f'\t{i}: {poly}')
		print()

if INFO_LAGVECTOR:
	for n,s,p in COMBINATIONS:
		print(f'This is the Lagrange vector for elements of type ({n},{s},{p}):')
		lag = lagrange_vector(n, s, p, 'c')
		for i,poly in enumerate(lag):
			print(f'\t{i}: {poly}')
		print()