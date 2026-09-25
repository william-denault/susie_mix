"""Independent numerical checks of the manuscript's VEB identities."""
from itertools import product
import re
from pathlib import Path
import numpy as np

ROOT = Path('C:/Document/Serieux/Travail/Data_analysis_and_papers/susie_mix')
tex = (ROOT / 'output/latex/slider_model.tex').read_text(encoding='utf-8')
clean = re.sub(r'(?<!\\)%[^\n]*', '', tex)
clean = re.sub(r'\\begin\{lstlisting\}.*?\\end\{lstlisting\}', '', clean, flags=re.S)
labels = re.findall(r'\\label\{([^}]+)\}', clean)
refs = re.findall(r'\\(?:eqref|ref)\{([^}]+)\}', clean)
assert len(labels) == len(set(labels)), 'Duplicate labels'
assert not set(refs) - set(labels), 'Undefined references'
bibs = re.findall(r'\\bibitem\{([^}]+)\}', clean)
cites = [key for group in re.findall(r'\\cite\{([^}]+)\}', clean) for key in group.split(',')]
assert not set(cites) - set(bibs), 'Undefined citations'
stack = []
for match in re.finditer(r'\\(begin|end)\{([^}]+)\}', clean):
    op, env = match.groups()
    if op == 'begin':
        stack.append(env)
    else:
        assert stack and stack.pop() == env, ('Mismatched environment', env)
assert not stack
brace_depth = 0
for token in re.findall(r'\\[\s\S]|[{}]', clean):
    if token == '{': brace_depth += 1
    if token == '}': brace_depth -= 1
    assert brace_depth >= 0, 'Unbalanced brace'
assert brace_depth == 0
section = (ROOT / 'output/latex/slider_veb_section.tex').read_text(encoding='utf-8')
assert section.strip() in tex
assert '% SLIDER_VEB_SECTION_INSERT' not in tex
print(f'LaTeX source: {len(labels)} unique labels; all references, citations, braces and environments valid.')

def lse(x):
    m = np.max(x)
    return m + np.log(np.sum(np.exp(x - m)))

def lognormal(y, covariance):
    sign, ld = np.linalg.slogdet(covariance)
    assert sign > 0
    return -.5 * (len(y) * np.log(2*np.pi) + ld + y @ np.linalg.solve(covariance,y))

def candidate_fit(x, h, r, V, s2, delta):
    z = x + delta*h
    S, T = z@z, z@r
    v = 1/(1/V + S/s2)
    m = v*T/s2
    bf = -.5*np.log1p(V*S/s2) + V*T*T/(2*s2*(s2+V*S))
    return bf, m, v

def maximize_slider(x, h, r, V, s2):
    A, B, D = 1+V*(x@x)/s2, V*(x@h)/s2, V*(h@h)/s2
    C, E = np.sqrt(V)*(x@r)/s2, np.sqrt(V)*(h@r)/s2
    poly = [-D*D, B*E*E-C*D*E-3*B*D, A*E*E-D*(C*C+A)-2*B*B, A*C*E-B*(C*C+A)]
    roots = np.roots(np.trim_zeros(poly,'f')) if np.any(poly) else []
    candidates = [-1., 0., 1.]
    candidates += [float(v.real) for v in roots if abs(v.imag)<1e-8 and -1<v.real<1]
    return max(candidates,key=lambda d: candidate_fit(x,h,r,V,s2,d)[0])

def posterior(X,H,r,V,s2,delta,pi):
    vals = np.array([candidate_fit(X[:,j],H[:,j],r,V,s2,delta[j]) for j in range(X.shape[1])])
    weights = np.log(pi)+vals[:,0]
    evidence_ratio = lse(weights)
    return np.exp(weights-evidence_ratio), vals[:,1], vals[:,2], evidence_ratio

def kl_component(a,m,v,V,pi):
    return np.sum(a*np.log(a/pi)) + .5*np.sum(a*((v+m*m)/V-1+np.log(V/v)))

def full_objective(X,H,y,delta,a,m,v,V,s2,pi):
    Z = X[None,:,:]+H[None,:,:]*delta[:,None,:]
    means = np.einsum('lp,lnp->ln',a*m,Z)
    seconds = np.einsum('lp,lnp,lnp->l',a*(v+m*m),Z,Z)
    erss = np.sum((y-means.sum(axis=0))**2) + np.sum(seconds-np.sum(means**2,axis=1))
    KL = np.array([kl_component(a[l],m[l],v[l],V[l],pi) for l in range(len(V))])
    F = -.5*len(y)*np.log(2*np.pi*s2)-erss/(2*s2)-KL.sum()
    return F, means, erss, KL

rng=np.random.default_rng(73229)
max_bf_error=0.
max_block_error=0.
min_step=float('inf')
min_gap=float('inf')
updates=0
for case in range(30):
    n,p,L=18,3,2
    Xraw=rng.binomial(2,.4,(n,p)).astype(float)
    Hraw=(Xraw==1).astype(float)
    # An orthonormal representation of the centered subspace.
    U=np.linalg.qr(np.column_stack([np.ones(n),rng.normal(size=(n,n-1))]))[0][:,1:]
    X,H=U.T@Xraw,U.T@Hraw
    y=U.T@(Xraw[:,0]*.7+Hraw[:,1]*.4+rng.normal(size=n))
    V=np.array([.3,.8]); s2=.9; pi=np.array([.2,.3,.5])
    delta=rng.uniform(-1,1,(L,p))
    a=rng.dirichlet(np.ones(p),L)
    m=rng.normal(0,.2,(L,p)); v=rng.uniform(.05,.3,(L,p))
    for iteration in range(6):
        for l in range(L):
            oldF,means,_,KL=full_objective(X,H,y,delta,a,m,v,V,s2,pi)
            r=y-(means.sum(axis=0)-means[l])
            oldZ=X+H*delta[l]
            old_erss=np.dot(r,r)-2*np.dot(r,means[l])+np.sum(a[l]*(v[l]+m[l]**2)*np.sum(oldZ**2,axis=0))
            old_local=-.5*len(y)*np.log(2*np.pi*s2)-old_erss/(2*s2)-KL[l]
            delta[l]=[maximize_slider(X[:,j],H[:,j],r,V[l],s2) for j in range(p)]
            a[l],m[l],v[l],logratio=posterior(X,H,r,V[l],s2,delta[l],pi)
            newF,_,_,_=full_objective(X,H,y,delta,a,m,v,V,s2,pi)
            local_evidence=lognormal(r,s2*np.eye(len(y)))+logratio
            max_block_error=max(max_block_error,abs((newF-oldF)-(local_evidence-old_local)))
            min_step=min(min_step,newF-oldF)
            assert newF>=oldF-1e-10
            assert abs((newF-oldF)-(local_evidence-old_local))<1e-10
            # Dense multivariate Gaussian integration independently verifies each BF.
            for j in range(p):
                z=X[:,j]+delta[l,j]*H[:,j]
                dense=lognormal(r,s2*np.eye(len(y))+V[l]*np.outer(z,z))-lognormal(r,s2*np.eye(len(y)))
                formula=candidate_fit(X[:,j],H[:,j],r,V[l],s2,delta[l,j])[0]
                max_bf_error=max(max_bf_error,abs(dense-formula))
                assert abs(dense-formula)<1e-10
            # Enumerate all p^L identities to integrate the entire model exactly.
            terms=[]
            for indices in product(range(p),repeat=L):
                cov=s2*np.eye(len(y)); lp=0.
                for k,j in enumerate(indices):
                    z=X[:,j]+delta[k,j]*H[:,j]
                    cov=cov+V[k]*np.outer(z,z)
                    lp+=np.log(pi[j])
                terms.append(lp+lognormal(y,cov))
            gap=lse(terms)-newF
            min_gap=min(min_gap,gap)
            assert gap>=-1e-10
            updates+=1
print(f'Numerical checks: {updates} component updates across 30 simulated designs.')
print(f'Maximum BF error against dense Gaussian integration: {max_bf_error:.3g}')
print(f'Maximum block-ELBO identity error: {max_block_error:.3g}')
print(f'Minimum ELBO step: {min_step:.3g}')
print(f'Minimum exact log-evidence minus ELBO: {min_gap:.3g}')
print('All checks passed. These checks do not establish calibration or global optimality.')
