"""1D shifted interface method (SIM) test: -(kappa u')' = f on (0,1),
kappa jumps at x_G (not a mesh node). Surrogate interface = nearest node x_t.
Duplicated DoF at x_t. Consistent form:
  a = sum int k u'v' - {kEu'}_G [v]_xt - {kEv'}_G [Eu]_G + sig [Eu]_G [Ev]_G
'naive' variant uses [Ev]_G in the first flux term (symmetric-looking) -> fails.
E = exact polynomial extrapolation from the element adjacent to x_t."""
import numpy as np
from numpy.polynomial import polynomial as npoly

def lagr(p):
    x = 0.5 * (1 - np.cos(np.pi * np.arange(p + 1) / p)) if p > 0 else np.array([0.5])
    return np.linalg.inv(np.vander(x, p + 1, increasing=True)), x

def ev(C, t, der=0):
    out = []
    for i in range(C.shape[1]):
        c = C[:, i]
        for _ in range(der): c = npoly.polyder(c)
        out.append(npoly.polyval(t, c))
    return np.array(out)  # [ndofs] for scalar t

def solve(p, n, km, kp, xG, f, uex, form="consistent", gamma=4.0):
    C, ref = lagr(p)
    h = 1.0 / n
    nodes_elem = lambda e: e * h + ref * h
    # surrogate node: nearest mesh node to xG
    it = int(round(xG / h))          # element boundary index owning x_t
    xt = it * h
    d = xG - xt
    # global CG dofs with duplicate at xt: left dofs: elements 0..it-1 -> it*p+1
    Nl = it * p + 1
    Nr = (n - it) * p + 1
    N = Nl + Nr
    A = np.zeros((N, N)); b = np.zeros(N)
    def dofs(e):
        if e < it:  return np.arange(e * p, e * p + p + 1)
        else:       return Nl + np.arange((e - it) * p, (e - it) * p + p + 1)
    nq = p + 2
    gx, gw = np.polynomial.legendre.leggauss(nq)
    gx = 0.5 * (gx + 1); gw = 0.5 * gw
    for e in range(n):
        k = km if e < it else kp
        idx = dofs(e)
        P1 = np.stack([ev(C, t, 1) for t in gx])   # [nq, ndof] d/dref
        P0 = np.stack([ev(C, t, 0) for t in gx])
        Ae = (k / h) * P1.T @ (gw[:, None] * P1)
        A[np.ix_(idx, idx)] += Ae
        xg = e * h + gx * h
        b[idx] += h * P0.T @ (gw * f(xg))
    # interface functionals (global coefficient vectors)
    eL, eR = it - 1, it                      # elements left/right of xt
    iL, iR = dofs(eL), dofs(eR)
    tL = (xG - eL * h) / h                   # ref coord of Gamma in left elem (>1)
    tR = (xG - eR * h) / h                   # in right elem (may be <0 or in (0,1))
    def func(elem_idx, t, der):
        v = np.zeros(N); v[elem_idx] = ev(C, t, der) / h ** der
        return v
    JEu = func(iL, tL, 0) - func(iR, tR, 0)          # [Eu]_Gamma
    Jv_t = func(iL, 1.0, 0) - func(iR, 0.0, 0)       # [v]_xt (plain)
    bm, bp = kp / (km + kp), km / (km + kp)          # harmonic weights
    FE = bm * km * func(iL, tL, 1) + bp * kp * func(iR, tR, 1)   # {kEu'}_Gamma
    khm = 2 * km * kp / (km + kp)
    sig = gamma * khm * max(p, 1) ** 2 / h
    Jv1 = JEu if form == "naive" else Jv_t
    A += -np.outer(Jv1, FE) - np.outer(FE, JEu) + sig * np.outer(JEu, JEu)
    # shifted-source correction: flux at xt = flux at Gamma + int_strip f
    strip = d * np.sum(gw * f(xt + gx * d))
    b += strip * Jv_t
    # strong BCs u(0)=uex(0,'m'), u(1)=uex(1,'p')
    for dof, val in ((0, uex(0.0, "m")), (N - 1, uex(1.0, "p"))):
        A[dof, :] = 0; A[dof, dof] = 1; b[dof] = val
    u = np.linalg.solve(A, b)
    # physical error: sample; strip (min(xt,xG),max(xt,xG)) uses the CORRECT side's extension
    err = 0.0
    for e in range(n):
        idx = dofs(e)
        for t in np.linspace(0.01, 0.99, 9):
            x = e * h + t * h
            side = "m" if x < xG else "p"
            # physical value: if x in strip, represented by the other field's ext
            if xG > xt and xt <= x < xG:      # strip owned by right field, physical=left
                uh = u[iL] @ (ev(C, (x - eL * h) / h))
            elif xG < xt and xG <= x < xt:    # strip owned by left field, physical=right
                uh = u[iR] @ (ev(C, (x - eR * h) / h))
            else:
                uh = u[idx] @ ev(C, t)
            err = max(err, abs(uh - uex(x, side)))
    return err

# --- Test 1: broken line, f=0, u(0)=0, u(1)=1 ---
km, kp, xG = 1.0, 10.0, 0.5 + 0.023  # interface off-node
F = 1.0 / ((xG / km) + (1 - xG) / kp)
def uex_line(x, side):
    return F * x / km if side == "m" else F * xG / km + F * (x - xG) / kp
zero = lambda x: 0.0 * x
print("Broken line (f=0), kappa: 1 -> 10, xG=0.523, P1..P3, n=8:")
for form in ("consistent", "naive"):
    for p in (1, 2, 3):
        e = solve(p, 8, km, kp, xG, zero, uex_line, form)
        print(f"  {form:10s} p={p}:  max err = {e:.2e}")

# --- Test 2: convergence, f=1 ---
def make_uex(km, kp, xG):
    # u = -x^2/(2k) + a x + c per side; u(0)=0, u(1)=1, [u]=[ku']=0
    # unknowns: am, ap, cp (cm=0). ku' = -x + k a
    # flux: -xG + km am = -xG + kp ap -> km am = kp ap =: t
    # value: -xG^2/(2km) + am xG = -xG^2/(2kp) + ap xG + cp
    # u(1)=1: -1/(2kp) + ap + cp = 1
    M = np.array([[km, -kp, 0], [xG, -xG, -1], [0, 1, 1]])
    r = np.array([0, xG**2 / 2 * (1 / km - 1 / kp), 1 + 1 / (2 * kp)])
    am, ap, cp = np.linalg.solve(M, r)
    def uex(x, side):
        return -x**2 / (2 * km) + am * x if side == "m" else -x**2 / (2 * kp) + ap * x + cp
    return uex
uex2 = make_uex(km, kp, xG)
one = lambda x: np.ones_like(x)
print("\nConvergence, f=1 (max nodal-sample err / rate):")
for p in (1, 2, 3):
    prev = None
    row = f"  p={p}: "
    for n in (8, 16, 32, 64):
        e = solve(p, n, km, kp, xG, one, uex2)
        rate = "" if prev is None else f"({np.log2(prev/e):4.1f})"
        row += f" {e:.2e}{rate}"
        prev = e
    print(row)

# --- Test 3: extreme jump ---
print("\nExtreme jump kappa: 1 -> 1e6, broken line, p=2, n=8:")
km2, kp2 = 1.0, 1e6
F2 = 1.0 / ((xG / km2) + (1 - xG) / kp2)
def uex3(x, side):
    return F2 * x / km2 if side == "m" else F2 * xG / km2 + F2 * (x - xG) / kp2
e = solve(2, 8, km2, kp2, xG, zero, uex3)
print(f"  max err = {e:.2e}")
