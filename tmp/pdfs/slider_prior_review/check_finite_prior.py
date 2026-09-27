"""Independent checks of the finite-grid derivation; not package implementation."""
from pathlib import Path
import re
import numpy as np

def logsumexp(values):
    peak = np.max(values)
    return peak + np.log(np.exp(values-peak).sum())

rng = np.random.default_rng(20260927)
d = np.linspace(-1, 1, 17)
assert np.array_equal(d[[0, 4, 8, 12, 16]], [-1, -.5, 0, .5, 1])
assert np.all(np.diff(d) == .125)
n, p, L = 24, 4, 3
raw = rng.binomial(2, [.25, .35, .45, .6], (n, p)).astype(float)
U = np.linalg.qr(np.column_stack([np.ones(n), rng.normal(size=(n, n-1))]))[0][:, 1:]
X = U.T @ raw
H = U.T @ (raw == 1).astype(float)
N = n-1
y = U.T @ (1.2*raw[:, 0] - .5*(raw[:, 0] == 1) + rng.normal(size=n))
Z = X[:, :, None] + H[:, :, None]*d
S_direct = np.sum(Z**2, axis=0)
A, B, D = (X*X).sum(0), (X*H).sum(0), (H*H).sum(0)
S = A[:, None] + 2*B[:, None]*d + D[:, None]*d**2
np.testing.assert_allclose(S, S_direct, rtol=1e-13, atol=1e-13)
pi = np.array([.1, .2, .3, .4])
w = np.ones(17)/17
V = np.array([.3, .8, .5])
s2 = .9

def ser(r, prior_var, noise, weights):
    T = (X.T@r)[:, None] + (H.T@r)[:, None]*d
    v = 1/(1/prior_var + S/noise)
    m = v*T/noise
    bf = -.5*np.log1p(prior_var*S/noise) + m*m/(2*v)
    logp = np.log(pi)[:, None] + np.log(weights)[None, :] + bf
    norm = logsumexp(logp)
    return np.exp(logp-norm), m, v, bf, norm

rho, m, v, bf, norm = ser(y, V[0], s2, w)
max_bf_error = 0.
for j in range(p):
    for k in range(17):
        z = Z[:, j, k]
        cov = s2*np.eye(N) + V[0]*np.outer(z, z)
        direct = -.5*(np.linalg.slogdet(cov)[1] - N*np.log(s2)
                        + y@np.linalg.solve(cov, y) - y@y/s2)
        max_bf_error = max(max_bf_error, abs(bf[j, k]-direct))
assert max_bf_error < 1e-10
compressed = X@(rho*m).sum(1) + H@(rho*m*d).sum(1)
expanded = np.einsum('jk,njk->n', rho*m, Z)
np.testing.assert_allclose(compressed, expanded, rtol=1e-13, atol=1e-13)

# The grid evidence agrees with an ordinary expanded 17p-column SER,
# and its zero-slider column agrees with the additive Gaussian formulas.
z_flat = Z.reshape(N, -1)
flat_v = 1/(1/V[0] + (z_flat*z_flat).sum(0)/s2)
flat_m = flat_v*(z_flat.T@y)/s2
flat_bf = -.5*np.log1p(V[0]*(z_flat*z_flat).sum(0)/s2) + flat_m**2/(2*flat_v)
np.testing.assert_allclose(flat_bf.reshape(p, 17), bf, rtol=1e-13, atol=1e-13)
add_v = 1/(1/V[0] + A/s2)
add_m = add_v*(X.T@y)/s2
np.testing.assert_allclose(m[:, 8], add_m, rtol=1e-13, atol=1e-13)
np.testing.assert_allclose(v[:, 8], add_v, rtol=1e-13, atol=1e-13)

# Start with prior variational factors and check each coordinate update,
# the common-w update, and both variance updates against the same ELBO.
rho = np.broadcast_to(pi[None, :, None]*w, (L, p, 17)).copy()
m = np.zeros((L, p, 17))
v = np.broadcast_to(V[:, None, None], (L, p, 17)).copy()

def objective():
    means = np.einsum('ljk,njk->ln', rho*m, Z)
    second = np.einsum('ljk,jk->l', rho*(v+m*m), S)
    erss = np.sum((y-means.sum(0))**2) + np.sum(second - (means*means).sum(1))
    kl_cat = np.sum(rho*(np.log(rho)-np.log(pi)[None, :, None]-np.log(w)[None, None, :]))
    kl_normal = .5*np.sum(rho*((v+m*m)/V[:, None, None]-1+np.log(V[:, None, None]/v)))
    return -.5*N*np.log(2*np.pi*s2)-erss/(2*s2)-kl_cat-kl_normal, means, erss

increments = []
for sweep in range(12):
    for ell in range(L):
        old, means, _ = objective()
        residual = y - means.sum(0) + means[ell]
        rho[ell], m[ell], v[ell], _, _ = ser(residual, V[ell], s2, w)
        increments.append(objective()[0]-old)
    old = objective()[0]
    w = rho.sum((0, 1))/L
    increments.append(objective()[0]-old)
    assert abs(w.sum()-1) < 1e-12
    old = objective()[0]
    V = (rho*(v+m*m)).sum((1, 2))
    increments.append(objective()[0]-old)
    old, _, erss = objective()
    s2 = erss/N
    increments.append(objective()[0]-old)
assert min(increments) > -1e-10, min(increments)

# The quadratic ERSS formula agrees with summing full covariance matrices.
_, means, erss = objective()
full_cov = np.zeros((N, N))
for ell in range(L):
    coeff = (rho[ell]*(v[ell]+m[ell]**2)).reshape(-1)
    full_cov += (z_flat*coeff)@z_flat.T - np.outer(means[ell], means[ell])
direct_erss = np.sum((y-means.sum(0))**2) + np.trace(full_cov)
np.testing.assert_allclose(erss, direct_erss, rtol=1e-12, atol=1e-12)

# Count-forced additive SNPs do not contribute to learned prior counts.
free = np.array([True, False, True, False])
counts = rho[:, free, :].sum((0, 1))
w_free = counts/counts.sum()
old_q = counts@np.log(w)
new_q = counts@np.log(w_free)
assert new_q >= old_q - 1e-12

source = Path(__file__).resolve().parents[3]/'output/latex/susie_slide_prior.tex'
tex = source.read_text(encoding='utf-8')
labels = re.findall(r'\\label\{([^}]+)\}', tex)
refs = re.findall(r'\\(?:eqref|ref)\{([^}]+)\}', tex)
assert len(labels) == len(set(labels))
assert set(refs) <= set(labels)
stack = []
for match in re.finditer(r'\\(begin|end)\{([^}]+)\}', tex):
    op, env = match.groups()
    if op == 'begin':
        stack.append(env)
    else:
        assert stack and stack.pop() == env
assert not stack
depth = 0
for token in re.findall(r'\\[\s\S]|[{}]', tex):
    if token == '{': depth += 1
    if token == '}': depth -= 1
    assert depth >= 0
assert depth == 0
print(f'Grid and sufficient statistics: PASS')
print(f'68 Bayes factors vs direct multivariate normals: max error {max_bf_error:.3g}')
print('Compressed fits, expanded SER, additive reduction, full-covariance ERSS: PASS')
print(f'{len(increments)} variational/EM coordinate updates: smallest ELBO change {min(increments):.6g}')
print('Count-filter prior update: PASS')
print(f'LaTeX source structure: {len(labels)} unique labels; references, environments and braces PASS')
