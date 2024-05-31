### Bezier subdivision pre-computations ###

from itertools import product
from sympy import symbols, prod, Matrix, det, collect, Rational

## Subscripting ##
def subscripts(b, i):
	return symbols([b+'_'+str(j) for j in range(i)])

def subscripts_arr(b, ii):
    return symbols([b + '_' + '_'.join(map(str, j)) for j in ii])

## Index set ##
def index_set(n, s, p):
	all_indices = product(range(p+1), repeat=n)
	valid_indices = [tup for tup in all_indices
		if sum(tup[:s]) <= p and all(x <= p for x in tup[s:])
	]
	return valid_indices

## Univariate Lagrange polynomial ##
def lagrange_uni(p, i, x):
	k_values = [k for k in range(0, p+1) if k != i]
	return prod((p * x - k) / (i - k) for k in k_values)

## Multivariate Lagrange polynomial ##
def lagrange(n, s, p, i, x_name):
	assert(len(i) == n)
	x = subscripts(x_name, n)
	i_slack = p - sum(i[:s])
	assert(i_slack >= 0)
	x_slack = 1 - sum(x[:s])
	return \
		lagrange_uni(i_slack, i_slack, x_slack) * \
		prod([lagrange_uni(i[k], i[k], x[k]) for k in range(s)]) * \
		prod([lagrange_uni(p, i[k], x[k]) for k in range(s,n)])

## Geometric map component ##
def geo_map_component(n, s, p, t, x_name, c_name):
	indices = index_set(n, s, p)
	cp_indices = [(tup + (t,)) for tup in indices]
	cp = subscripts_arr(c_name, cp_indices)
	lag = [lagrange(n, s, p, k, x_name) for k in indices]
	return sum(cp[k] * lag[k] for k in range(len(indices)))
	
# print(geo_map_component(2,1,1,1,'u','X'))
## Geometric map ##
def geo_map(n, s, p, x_name, c_names, t_name):
	assert(len(c_names) == n)
	T = symbols(t_name)
	return [
		(1 - T) * geo_map_component(n, s, p, 0, x_name, c) + \
		T * geo_map_component(n, s, p, 1, x_name, c) \
		for c in c_names
	] + [T]

## Jacobian determinant ##
def jac_det(n, s, p, x_name, c_names, t_name):
	xt = subscripts(x_name, n) + [symbols(t_name)]
	gmap = geo_map(n, s, p, x_name, c_names, t_name)
	res = det(Matrix([[poly.diff(v) for v in xt] for poly in gmap]))
	res = collect(res, symbols(t_name))
	return res

## Domain points for the Jacobian determinant ##
def domain_pts_J(n, s, p):
	simplex_den = max(1, n * p - s)
	tensor_den = max(1, n * p - 1)
	time_den = n
	all_indices = product(range(n*p+1), repeat=n+1)
	valid_indices = [tup for tup in all_indices
		if sum(tup[:s]) <= simplex_den \
		and all(x <= tensor_den for x in tup[s:n]) \
		and tup[-1] <= time_den
	]
	return [
        tuple(
            Rational(elem, simplex_den) if i < s else
            Rational(elem, time_den) if i == n else
            Rational(elem, tensor_den)
            for i, elem in enumerate(tup)
        )
        for tup in valid_indices
    ]

## Lagrange vector ##
def lagrange_vector(n, s, p, c_names):
	x_name = 'u'
	t_name = 't'
	poly = jac_det(n, s, p, x_name, c_names, t_name)
	pts = domain_pts_J(n, s, p)
	xt = subscripts(x_name, n) + [symbols(t_name)]
	rule = lambda pt: {xt[k]: pt[k] for k in range(n+1)}
	lag_vec = [poly.subs(rule(pt)) for pt in pts]
	return lag_vec


## Time subdivision maps ##
def time_sub_map(pt, shift):
	h = 1 if shift else 0
	return(pt[:-1] + ((pt[-1] + h) / 2,))

def space_sub_map(pt, n, s, q):
	assert(1 <= s <= n)
	assert(len(pt) == n+1)
	assert(0 <= q < (2<<(n+1)))
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


## Tests ## 
for p in domain_pts_J(2, 2, 1):
	print(p)
	for q in range(8):
		print(f'=={q}==> ', space_sub_map(p, 2, 2, q))
	print()