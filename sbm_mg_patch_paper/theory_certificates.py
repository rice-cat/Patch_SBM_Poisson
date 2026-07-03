"""Single-cell SBM block spectra: flat face vs staircase step corner vs chamfer.

Cell K = (0,1)^2, Q_p Lagrange (Chebyshev-Lobatto nodes). SBM form on K:
  a(u,v) = int_K grad u . grad v
           - int_F (dn u) v  - int_F (dn v) Eu  + sigma int_F Eu Ev
with Eu(x) = u(x + d) (exact polynomial extrapolation, k>=p Taylor).

Configs:
  flat(s):    face {y=1},          d=(0,s)
  corner(s):  faces {y=1},{x=1},   d=(0,s),(s,0)      [staircase step cell]
  diag(s):    faces {y=1},{x=1},   d=s(1,1)/sqrt2 both [pocket-directed proj.]
  chamfer(s): volume=triangle x+y<=1, face=hypotenuse, d=s(1,1)/sqrt2, Q_p cut
Outputs per (p,s,config):
  lmin(H) full block, lmin(H) on zero-extension subspace (funcs vanishing on
  interior edges x=0,y=0 -> global indefiniteness certificate), lmax(H),
  max |Im eig(A)|.
"""
import numpy as np
from numpy.polynomial import polynomial as npoly


def cheb_lobatto(p):
    if p == 0:
        return np.array([0.5])
    k = np.arange(p + 1)
    return 0.5 * (1 - np.cos(np.pi * k / p))


def basis_coeffs(p):
    x = cheb_lobatto(p)
    V = np.vander(x, p + 1, increasing=True)
    return np.linalg.inv(V), x  # C[:, i] = monomial coeffs of phi_i


def eval1d(C, pts, deriv=0):
    cols = []
    for i in range(C.shape[1]):
        c = C[:, i]
        for _ in range(deriv):
            c = npoly.polyder(c)
        cols.append(npoly.polyval(pts, c))
    return np.stack(cols, axis=1)  # [npts, p+1]


def eval2d(C, pts, dx=0, dy=0):
    Fx = eval1d(C, pts[:, 0], dx)
    Fy = eval1d(C, pts[:, 1], dy)
    n = C.shape[1]
    return (Fx[:, :, None] * Fy[:, None, :]).reshape(len(pts), n * n)


def gauss01(n):
    x, w = np.polynomial.legendre.leggauss(n)
    return 0.5 * (x + 1), 0.5 * w


def assemble(p, s, config, gamma):
    C, nodes = basis_coeffs(p)
    n1 = p + 1
    N = n1 * n1
    nq = p + 3
    xq, wq = gauss01(nq)

    A = np.zeros((N, N))

    # ---- volume ----
    if config == "chamfer":
        # Duffy on triangle x+y<=1: x=u, y=v(1-u), J=(1-u)
        U, V = np.meshgrid(xq, xq, indexing="ij")
        W = np.outer(wq, wq) * (1 - U)
        pts = np.stack([U.ravel(), (V * (1 - U)).ravel()], axis=1)
        Wv = W.ravel()
    else:
        U, V = np.meshgrid(xq, xq, indexing="ij")
        pts = np.stack([U.ravel(), V.ravel()], axis=1)
        Wv = np.outer(wq, wq).ravel()
    Gx = eval2d(C, pts, dx=1)
    Gy = eval2d(C, pts, dy=1)
    A += Gx.T @ (Wv[:, None] * Gx) + Gy.T @ (Wv[:, None] * Gy)

    # ---- faces ----
    sigma = gamma * max(p, 1) ** 2  # h=1
    faces = []
    if config == "flat":
        fp = np.stack([xq, np.ones(nq)], axis=1)
        faces.append((fp, wq, (0.0, 1.0), np.array([0.0, s])))
    elif config in ("corner", "diag"):
        d_top = np.array([0.0, s]) if config == "corner" else s / np.sqrt(2) * np.array([1.0, 1.0])
        d_rgt = np.array([s, 0.0]) if config == "corner" else s / np.sqrt(2) * np.array([1.0, 1.0])
        fp = np.stack([xq, np.ones(nq)], axis=1)
        faces.append((fp, wq, (0.0, 1.0), d_top))
        fp = np.stack([np.ones(nq), xq], axis=1)
        faces.append((fp, wq, (1.0, 0.0), d_rgt))
    elif config == "chamfer":
        t = xq
        fp = np.stack([t, 1 - t], axis=1)
        faces.append((fp, np.sqrt(2) * wq, (1 / np.sqrt(2), 1 / np.sqrt(2)),
                      s / np.sqrt(2) * np.array([1.0, 1.0])))

    for fp, fw, nvec, d in faces:
        Phi = eval2d(C, fp)
        dPhin = nvec[0] * eval2d(C, fp, dx=1) + nvec[1] * eval2d(C, fp, dy=1)
        PhiE = eval2d(C, fp + d[None, :])
        Wf = fw[:, None]
        # a(u,v): A[v_idx, u_idx]
        A += -Phi.T @ (Wf * dPhin)          # -int (dn u) v
        A += -dPhin.T @ (Wf * PhiE)         # -int (dn v) Eu
        A += sigma * PhiE.T @ (Wf * PhiE)   # +sigma int Eu Ev
    return A, nodes


def zero_ext_mask(p, nodes):
    # drop nodes on interior edges x=0 or y=0 -> v extendable by zero in CG
    n1 = p + 1
    keep = []
    for i in range(n1):
        for j in range(n1):
            if nodes[i] > 1e-12 and nodes[j] > 1e-12:
                keep.append(i * n1 + j)
    return np.array(keep)


def analyze(p, s, config, gamma=4.0):
    A, nodes = assemble(p, s, config, gamma)
    H = 0.5 * (A + A.T)
    ev = np.linalg.eigvalsh(H)
    keep = zero_ext_mask(p, nodes)
    evS = np.linalg.eigvalsh(H[np.ix_(keep, keep)])
    lam = np.linalg.eigvals(A)
    return ev[0], evS[0], ev[-1], np.max(np.abs(lam.imag))


if __name__ == "__main__":
    svals = [0.25, 0.5, 0.75, 1.0, np.sqrt(2)]
    for gamma in (4.0,):
        print(f"=== sigma = {gamma} p^2 ===")
        for config in ("flat", "corner", "diag", "chamfer"):
            print(f"-- {config} --")
            hdr = "p  " + "".join(f"| s={s:5.3f}: lminH lminH_S lmaxH maxIm " for s in svals)
            print("p    s      lminH      lminH_S    lmaxH     maxImEig")
            for p in (1, 2, 3, 4):
                for s in svals:
                    l0, l0S, lM, mi = analyze(p, s, config, gamma)
                    print(f"{p}  {s:5.3f}  {l0:10.3e} {l0S:10.3e} {lM:9.3e} {mi:9.3e}")
            print()


# ---------- P_p triangle (remeshed chamfer) and inward shifts ----------
def tri_basis(p):
    nodes = []
    for i in range(p + 1):
        for j in range(p + 1 - i):
            nodes.append((i / p if p else 1 / 3, j / p if p else 1 / 3))
    nodes = np.array(nodes)
    exps = [(a, b) for a in range(p + 1) for b in range(p + 1 - a)]
    V = np.array([[x ** a * y ** b for (a, b) in exps] for (x, y) in nodes])
    return np.linalg.inv(V), nodes, exps


def tri_eval(Ci, exps, pts, dx=0, dy=0):
    x, y = pts[:, 0], pts[:, 1]
    cols = []
    for (a, b) in exps:
        ca, cb = 1.0, 1.0
        aa, bb = a, b
        for _ in range(dx):
            ca *= aa; aa -= 1
        for _ in range(dy):
            cb *= bb; bb -= 1
        if aa < 0 or bb < 0 or ca == 0 or cb == 0:
            cols.append(np.zeros(len(pts)))
        else:
            cols.append(ca * cb * x ** aa * y ** bb)
    M = np.stack(cols, axis=1)  # monomial evals
    return M @ Ci


def assemble_tri(p, s, gamma):
    Ci, nodes, exps = tri_basis(p)
    nq = p + 3
    xq, wq = gauss01(nq)
    U, V = np.meshgrid(xq, xq, indexing="ij")
    W = (np.outer(wq, wq) * (1 - U)).ravel()
    pts = np.stack([U.ravel(), (V * (1 - U)).ravel()], axis=1)
    Gx = tri_eval(Ci, exps, pts, dx=1)
    Gy = tri_eval(Ci, exps, pts, dy=1)
    A = Gx.T @ (W[:, None] * Gx) + Gy.T @ (W[:, None] * Gy)
    sigma = gamma * max(p, 1) ** 2
    t = xq
    fp = np.stack([t, 1 - t], axis=1)
    fw = np.sqrt(2) * wq
    d = s / np.sqrt(2) * np.array([1.0, 1.0])
    Phi = tri_eval(Ci, exps, fp)
    dPhin = (tri_eval(Ci, exps, fp, dx=1) + tri_eval(Ci, exps, fp, dy=1)) / np.sqrt(2)
    PhiE = tri_eval(Ci, exps, fp + d[None, :])
    Wf = fw[:, None]
    A += -Phi.T @ (Wf * dPhin) - dPhin.T @ (Wf * PhiE) + sigma * PhiE.T @ (Wf * PhiE)
    return A, nodes


def analyze_tri(p, s, gamma=4.0):
    A, nodes = assemble_tri(p, s, gamma)
    H = 0.5 * (A + A.T)
    ev = np.linalg.eigvalsh(H)
    keep = [k for k, (x, y) in enumerate(nodes) if x > 1e-12 and y > 1e-12]
    evS = np.linalg.eigvalsh(H[np.ix_(keep, keep)]) if keep else np.array([np.nan])
    lam = np.linalg.eigvals(A)
    return ev[0], evS[0], ev[-1], np.max(np.abs(lam.imag))


def extra():
    print("=== P_p TRIANGLE (remeshed chamfer), sigma=4p^2 ===")
    print("p    s      lminH      lminH_S    lmaxH     maxImEig")
    for p in (1, 2, 3, 4):
        for s in (0.25, 0.5, 0.707, 1.0):
            l0, l0S, lM, mi = analyze_tri(p, s)
            print(f"{p}  {s:5.3f}  {l0:10.3e} {l0S:10.3e} {lM:9.3e} {mi:9.3e}")
    print()
    print("=== INWARD shift, flat face, d=(0,-s), sigma=4p^2 ===")
    print("p    s      lminH      lminH_S    lmaxH     maxImEig")
    for p in (1, 2, 3, 4):
        for s in (0.25, 0.5, 0.75, 0.9):
            A, nodes = assemble(p, s, "flat", 4.0)
            # rebuild with inward shift
            A2, _ = assemble_inward(p, s, 4.0)
            H = 0.5 * (A2 + A2.T)
            ev = np.linalg.eigvalsh(H)
            keep = zero_ext_mask(p, nodes)
            evS = np.linalg.eigvalsh(H[np.ix_(keep, keep)])
            lam = np.linalg.eigvals(A2)
            print(f"{p}  {s:5.3f}  {ev[0]:10.3e} {evS[0]:10.3e} {ev[-1]:9.3e} {np.max(np.abs(lam.imag)):9.3e}")
    print()
    print("=== penalty sensitivity, staircase corner, s=1.0 ===")
    print("p  gamma   lminH_S")
    for p in (3, 4):
        for gamma in (2.0, 4.0, 16.0, 64.0):
            _, l0S, _, _ = analyze(p, 1.0, "corner", gamma)
            print(f"{p}  {gamma:5.1f}  {l0S:10.3e}")


def assemble_inward(p, s, gamma):
    C, nodes = basis_coeffs(p)
    n1 = p + 1
    nq = p + 3
    xq, wq = gauss01(nq)
    U, V = np.meshgrid(xq, xq, indexing="ij")
    pts = np.stack([U.ravel(), V.ravel()], axis=1)
    Wv = np.outer(wq, wq).ravel()
    Gx = eval2d(C, pts, dx=1)
    Gy = eval2d(C, pts, dy=1)
    A = Gx.T @ (Wv[:, None] * Gx) + Gy.T @ (Wv[:, None] * Gy)
    sigma = gamma * max(p, 1) ** 2
    fp = np.stack([xq, np.ones(nq)], axis=1)
    d = np.array([0.0, -s])
    Phi = eval2d(C, fp)
    dPhin = eval2d(C, fp, dy=1)
    PhiE = eval2d(C, fp + d[None, :])
    Wf = wq[:, None]
    A += -Phi.T @ (Wf * dPhin) - dPhin.T @ (Wf * PhiE) + sigma * PhiE.T @ (Wf * PhiE)
    return A, nodes


extra()
"""Corrected certificates + truncated-Taylor extrapolation order k experiment."""
import numpy as np
from math import comb, factorial
# (definitions above)

def phiE_truncated(C, fp, d, k):
    out = 0.0
    for j in range(k + 1):
        term = 0.0
        for a in range(j + 1):
            b = j - a
            term = term + comb(j, a) * (d[0] ** a) * (d[1] ** b) * eval2d(C, fp, dx=a, dy=b)
        out = out + term / factorial(j)
    return out

def assemble_k(p, s, config, gamma, k, inward=False):
    C, nodes = basis_coeffs(p)
    nq = p + 3
    xq, wq = gauss01(nq)
    U, V = np.meshgrid(xq, xq, indexing="ij")
    pts = np.stack([U.ravel(), V.ravel()], axis=1)
    Wv = np.outer(wq, wq).ravel()
    Gx = eval2d(C, pts, dx=1); Gy = eval2d(C, pts, dy=1)
    A = Gx.T @ (Wv[:, None] * Gx) + Gy.T @ (Wv[:, None] * Gy)
    sigma = gamma * p ** 2
    sgn = -1.0 if inward else 1.0
    faces = []
    fp = np.stack([xq, np.ones(nq)], axis=1)
    faces.append((fp, wq, (0.0, 1.0), sgn * np.array([0.0, s])))
    if config == "corner":
        fp = np.stack([np.ones(nq), xq], axis=1)
        faces.append((fp, wq, (1.0, 0.0), sgn * np.array([s, 0.0])))
    for fp, fw, nvec, d in faces:
        Phi = eval2d(C, fp)
        dPhin = nvec[0] * eval2d(C, fp, dx=1) + nvec[1] * eval2d(C, fp, dy=1)
        PhiE = phiE_truncated(C, fp, d, k)
        Wf = fw[:, None]
        A += -Phi.T @ (Wf * dPhin) - dPhin.T @ (Wf * PhiE) + sigma * PhiE.T @ (Wf * PhiE)
    return A, nodes

def cert(p, s, config, gamma, k, inward=False):
    A, nodes = assemble_k(p, s, config, gamma, k, inward)
    H = 0.5 * (A + A.T)
    n1 = p + 1
    keep = []
    for i in range(n1):
        for j in range(n1):
            x, y = nodes[i], nodes[j]
            if config == "flat":
                ok = x > 1e-12 and x < 1 - 1e-12 and y > 1e-12  # 3 interior edges
            else:
                ok = x > 1e-12 and y > 1e-12                     # 2 interior edges
            if ok: keep.append(i * n1 + j)
    evS = np.linalg.eigvalsh(H[np.ix_(keep, keep)])
    return evS[0]

print("=== corrected certificate lminH_S (flat: 3 interior edges) ===")
print("cfg     p  s      k=p     k=3     k=2     k=1")
for config in ("flat", "corner"):
    for p in (3, 4):
        for s in (0.5, 0.75, 1.0, np.sqrt(2)):
            row = [cert(p, s, config, 4.0, k) for k in (p, 3, 2, 1)]
            print(f"{config:7s} {p}  {s:5.3f} " + " ".join(f"{v:8.4f}" for v in row))
print()
print("=== inward shift (lambda>0), flat, corrected cert ===")
print("p  s      k=p     k=2     k=1")
for p in (2, 3, 4):
    for s in (0.25, 0.5, 0.75):
        row = [cert(p, s, "flat", 4.0, k, inward=True) for k in (p, 2, 1)]
        print(f"{p}  {s:5.3f} " + " ".join(f"{v:8.4f}" for v in row))
"""PDE-substituted extension: replace normal derivatives >=2 via -Lap u = f (f=0 part):
 d_y^{2m} -> (-1)^m d_x^{2m},  d_y^{2m+1} -> (-1)^m d_x^{2m} d_y  (face normal = y).
All terms evaluated ON the face -> no extrapolation beyond first normal derivative."""
import numpy as np
from math import factorial
# (definitions above)

def phiE_pde(C, fp, delta, k, normal):  # normal: 'y' or 'x'
    out = 0.0
    for j in range(k + 1):
        m, r = divmod(j, 2)  # j = 2m + r
        sgn = (-1) ** m
        if normal == 'y':
            T = eval2d(C, fp, dx=2 * m, dy=r)
        else:
            T = eval2d(C, fp, dx=r, dy=2 * m)
        out = out + sgn * (delta ** j) / factorial(j) * T
    return out

def assemble_pde(p, s, config, gamma, k):
    C, nodes = basis_coeffs(p)
    nq = p + 3
    xq, wq = gauss01(nq)
    U, V = np.meshgrid(xq, xq, indexing="ij")
    pts = np.stack([U.ravel(), V.ravel()], axis=1)
    Wv = np.outer(wq, wq).ravel()
    Gx = eval2d(C, pts, dx=1); Gy = eval2d(C, pts, dy=1)
    A = Gx.T @ (Wv[:, None] * Gx) + Gy.T @ (Wv[:, None] * Gy)
    sigma = gamma * p ** 2
    faces = [("y", np.stack([xq, np.ones(nq)], axis=1), (0.0, 1.0))]
    if config == "step":
        faces.append(("x", np.stack([np.ones(nq), xq], axis=1), (1.0, 0.0)))
    for nrm, fp, nvec in faces:
        Phi = eval2d(C, fp)
        dPhin = nvec[0] * eval2d(C, fp, dx=1) + nvec[1] * eval2d(C, fp, dy=1)
        PhiE = phiE_pde(C, fp, s, k, nrm)
        Wf = wq[:, None]
        A += -Phi.T @ (Wf * dPhin) - dPhin.T @ (Wf * PhiE) + sigma * PhiE.T @ (Wf * PhiE)
    return A, nodes

def cert(p, s, config, gamma=4.0, k=None):
    k = p if k is None else k
    A, nodes = assemble_pde(p, s, config, gamma, k)
    H = 0.5 * (A + A.T)
    n1 = p + 1
    keep = []
    for i in range(n1):
        for j in range(n1):
            x, y = nodes[i], nodes[j]
            ok = (x > 1e-12 and x < 1 - 1e-12 and y > 1e-12) if config == "flat" \
                 else (x > 1e-12 and y > 1e-12)
            if ok: keep.append(i * n1 + j)
    evS = np.linalg.eigvalsh(H[np.ix_(keep, keep)])
    ev = np.linalg.eigvalsh(H)
    lam = np.linalg.eigvals(A)
    return evS[0], ev[0], ev[-1], np.max(np.abs(lam.imag))

print("=== PDE-substituted extension (f-part omitted), k=p, sigma=4p^2 ===")
print("config  p   s      lminH_S    lminH     lmaxH     maxIm")
for config in ("flat", "step"):
    for p in (3, 4):
        for s in (0.5, 0.75, 1.0, np.sqrt(2)):
            l0S, l0, lM, mi = cert(p, s, config)
            print(f"{config:5s}  {p}  {s:5.3f}  {l0S:9.4f}  {l0:9.4f} {lM:9.3e} {mi:8.4f}")
print()
print("=== amplification check: lmaxH growth vs s (was ~sigma*T_p^2 for Taylor) ===")
for p in (4,):
    for s in (0.5, 1.0, np.sqrt(2)):
        _, _, lM, _ = cert(p, s, "flat")
        print(f"p={p} s={s:5.3f}  lmaxH={lM:9.3e}")
