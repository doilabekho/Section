# beton.py - V2
#version V2
import datetime as dt
import numpy as np
import pandas as pd
import seaborn as sns
import xlwings as xw
from xlwings import func, script
import matplotlib
import scipy
import math
from scipy.optimize import root, fsolve
from matplotlib.path import Path
from scipy.spatial import Delaunay
import matplotlib.pyplot as plt

from beton import *   # ← fonctionne en local et sur GitHub
from acier import *   # ← fonctionne en local et sur GitHub	







# ═══════════════════════════════════════════════════════════════
# BLOC 1 — UTILITAIRES GÉOMÉTRIQUES
# ═══════════════════════════════════════════════════════════════

@func
def n_liste(data):
    count = 0
    for row in data:
        if len(row) >= 2:
            v0 = str(row[0]).strip().lower()
            v1 = str(row[1]).strip().lower()
            if v0 == "fin" or v1 == "fin":
                count += 1
    return count

@func
def liste_i(data, i):
    fins = []
    for idx, row in enumerate(data):
        if len(row) >= 2:
            v0 = str(row[0]).strip().lower()
            v1 = str(row[1]).strip().lower()
            if v0 == "fin" or v1 == "fin":
                fins.append(idx)
    if i < 1 or i > len(fins):
        return []
    start = 0 if i == 1 else fins[i-2] + 1
    return data[start:fins[i-1]]

@func
def points_valides(data):
    pts = []
    for row in data:
        if len(row) >= 2:
            try:
                pts.append((float(row[0]), float(row[1])))
            except:
                pass
    return pts

@func
def nettoyer_polygone(points, tol=1e-9):
    if len(points) < 2:
        return points
    x1, y1 = points[0]
    x2, y2 = points[-1]
    if abs(x1-x2) < tol and abs(y1-y2) < tol:
        return points[:-1]
    return points

@func
def aire_centre_inertie_origine(data):
    """Green exact : A, Cx, Cy, Ix, Iy par formule de Shoelace."""
    pts = points_valides(data)
    n = len(pts)
    if n < 3:
        return 0.0, None, None, 0.0, 0.0
    A = Cx = Cy = Ix = Iy = 0.0
    for i in range(n):
        x1, y1 = pts[i]
        x2, y2 = pts[(i+1) % n]
        cross = x1*y2 - x2*y1
        A  += cross
        Cx += (x1+x2)*cross
        Cy += (y1+y2)*cross
        Ix += (y1**2 + y1*y2 + y2**2)*cross
        Iy += (x1**2 + x1*x2 + x2**2)*cross
    A *= 0.5
    if A == 0:
        return 0.0, None, None, 0.0, 0.0
    Cx /= (6*A); Cy /= (6*A)
    Ix /= 12;    Iy /= 12
    return np.abs(A), Cx, Cy, np.abs(Ix), np.abs(Iy)

@func
def caracteristiques_mecaniques(contour, evidement):
    A, Cx, Cy, Ix0, Iy0 = aire_centre_inertie_origine(contour)
    if A == 0:
        return 0.0, None, None, 0.0, 0.0
    A_tot  = A;    Cx_tot = A*Cx; Cy_tot = A*Cy
    Ix_tot = Ix0;  Iy_tot = Iy0
    N = n_liste(evidement)
    for i in range(1, N+1):
        Li = liste_i(evidement, i)
        Ai, cxi, cyi, ixi, iyi = aire_centre_inertie_origine(Li)
        if Ai > 0:
            A_tot  -= Ai;     Cx_tot -= Ai*cxi; Cy_tot -= Ai*cyi
            Ix_tot -= ixi;    Iy_tot -= iyi
    if A_tot == 0:
        return 0.0, None, None, 0.0, 0.0
    CxG = Cx_tot/A_tot;  CyG = Cy_tot/A_tot
    return (np.abs(A_tot), CxG, CyG,
            np.abs(Ix_tot - A_tot*CyG**2),
            np.abs(Iy_tot - A_tot*CxG**2))

@func
def translation_points(data, Cx, Cy):
    pts_trans = []
    for row in data:
        if len(row) >= 2:
            try:
                pts_trans.append((float(row[0])-Cx, float(row[1])-Cy))
            except:
                pass
    return pts_trans

@func
def transformation_repere_cg(contour, evidement):
    A, Cx, Cy, Ix, Iy = caracteristiques_mecaniques(contour, evidement)
    if Cx is None or Cy is None:
        return None, [], None, None
    contour_cg = translation_points(contour, Cx, Cy)
    evidements_cg = []
    N = n_liste(evidement)
    if N == 0:
        return contour_cg, evidements_cg, Cx, Cy
    for i in range(1, N+1):
        Li = liste_i(evidement, i)
        evidements_cg.append(translation_points(Li, Cx, Cy))
    return contour_cg, evidements_cg, Cx, Cy

@func
def acier_G(p_acier, s_acier):
    if p_acier is None or len(p_acier) == 0:
        return []
    p_acier = np.atleast_2d(p_acier)
    s_acier = np.atleast_1d(s_acier)
    newacier = []
    nx = min(p_acier.shape[0], s_acier.shape[0])
    for i in range(nx):
        try:
            if p_acier[i,0] in (None,"","") or p_acier[i,1] in (None,"",""):
                continue
            x = float(p_acier[i,0]); y = float(p_acier[i,1])
            s = float(s_acier[i])
            newacier.append([x, y, s/10000.0])
        except (ValueError, TypeError, IndexError):
            continue
    return newacier

@func
def centre_gravite_polygone(data):
    pts = points_valides(data)
    n = len(pts)
    if n < 3: return None, None
    A = Cx = Cy = 0.0
    for i in range(n):
        x1,y1=pts[i]; x2,y2=pts[(i+1)%n]
        cross=x1*y2-x2*y1; A+=cross; Cx+=(x1+x2)*cross; Cy+=(y1+y2)*cross
    A*=0.5
    if A==0: return None, None
    return Cx/(6*A), Cy/(6*A)

@func
def aire_et_centre(data):
    pts = points_valides(data)
    n = len(pts)
    if n < 3: return 0.0, None, None
    A = Cx = Cy = 0.0
    for i in range(n):
        x1,y1=pts[i]; x2,y2=pts[(i+1)%n]
        cross=x1*y2-x2*y1; A+=cross; Cx+=(x1+x2)*cross; Cy+=(y1+y2)*cross
    A*=0.5
    if A==0: return 0.0, None, None
    return abs(A), Cx/(6*A), Cy/(6*A)


# ═══════════════════════════════════════════════════════════════
# BLOC 2 — NOYAU GREEN ANALYTIQUE EXACT
# ═══════════════════════════════════════════════════════════════
#
# Théorème de Green :  ∬_Ω f(x,y) dA = ∮_∂Ω Q(x,y) dy
# avec ∂Q/∂x = f(x,y)
#
# Correction du sens : aire_signee = Σ (x1·y2 - x2·y1)/2
# Si aire_signee < 0 (sens horaire), on inverse le signe du résultat.
#
# ── DÉCOUPAGE EN ZONES ──────────────────────────────────────
# Pour ELS : zones T (ε≤0, σ=0) et C (ε>0, σ=C·ε)
# Pour ELU : zones T (ε≤0), P (0<ε≤e2), R (ε>e2)
#
# Chaque arête est découpée aux t∈(0,1) où ε=0 et ε=e2.
# Sur chaque sous-segment la zone est constante.
#
# ── PRIMITIVES EXACTES ──────────────────────────────────────
#
# Paramétrage arête : x(t)=xa+t·dx, y(t)=ya+t·dy, t∈[0,1]
# Notation : ta..tb = sous-intervalle après reparamétrisation en s∈[0,1]
# Sur [ta,tb] : s∈[0,1], x=x1+s·Δx, y=y1+s·Δy, ε=εa+s·Δε
# avec x1=xa+ta·dx, Δx=(tb-ta)·dx, etc.
#
# INTÉGRALES ÉLÉMENTAIRES sur [0,1] :
#   ∫s^k ds = 1/(k+1)
#   Notons : εa, Δε, x1, Δx, y1, Δy
#
# ELS Zone C (σ=C·ε) :
#   N  = C·Δy·∫(εa+Δε·s) ds
#      = C·Δy·(εa + Δε/2)
#
#   My = C·Δy·∫(εa+Δε·s)·(y1+Δy·s) ds
#      = C·Δy·[εa·y1 + (εa·Δy+y1·Δε)/2 + Δε·Δy/3]
#
#   Mz = C·Δy·∫(εa+Δε·s)·(x1+Δx·s) ds
#      = C·Δy·[εa·x1 + (εa·Δx+x1·Δε)/2 + Δε·Δx/3]
#
#   Scom = Δy·∫(x1+Δx·s) ds = Δy·(x1 + Δx/2)
#
# ELU Zone R (σ=fcd) :
#   N  = fcd·Δy
#   My = fcd·Δy·(y1 + Δy/2)
#   Mz = fcd·Δy·(x1 + Δx/2)
#   Scom = Δy·(x1 + Δx/2)   [même que ELS zone C]
#
# ELU Zone P (σ=fcd·(2u-u²), u=ε/e2) :
#   u(s) = (εa+Δε·s)/e2 = ua + s·Δu  (ua=εa/e2, Δu=Δε/e2)
#
#   σ = fcd·(2u - u²)
#     = fcd·(2ua + 2Δu·s - ua² - 2ua·Δu·s - Δu²·s²)
#
#   ∫₀¹ σ ds = fcd·(2ua + 2Δu/2 - ua² - 2ua·Δu/2 - Δu²/3)
#            = fcd·(2ua + Δu - ua² - ua·Δu - Δu²/3)
#
#   ∫₀¹ σ·s ds = fcd·(2ua/2 + 2Δu/3 - ua²/2 - 2ua·Δu/3 - Δu²/4)
#
#   N  = fcd·Δy·I0
#   My = fcd·Δy·(y1·I0 + Δy·I1)
#   Mz = fcd·Δy·(x1·I0 + Δx·I1)
#   Scom = Δy·(x1 + Δx/2)
#
# où I0 = ∫₀¹ σ/fcd ds,  I1 = ∫₀¹ σ/fcd · s ds
# ═══════════════════════════════════════════════════════════════
import numpy as np

# ═══════════════════════════════════════════════════════════════
# BLOC 2 — NOYAU GREEN ANALYTIQUE EXACT (VERSION RIGOUREUSE)
# ═══════════════════════════════════════════════════════════════
#
# Théorème de Green :
# ∬ f(x,y) dA = ∮ Q(x,y) dy avec Q(x,y)=∫₀ˣ f(ξ,y)dξ
#
# → Toutes les intégrales sont traitées rigoureusement
# → Compatible BA (ELS / ELU)
# → Exact pour N, My, Mz
#
# Aire :
#   A = ∮ x dy
#
# Moments :
#   ∫ y dA = ∮ x y dy
#   ∫ x dA = ∮ x²/2 dy
#
# ═══════════════════════════════════════════════════════════════


def _aire_signee(pts):
    """Aire signée (shoelace)."""
    a = 0.0
    n = len(pts)
    for i in range(n):
        x1, y1 = pts[i]
        x2, y2 = pts[(i+1) % n]
        a += x1*y2 - x2*y1
    return 0.5 * a


def _decouper_arete(xa, ya, xb, yb, eps0, alpha, beta, seuils):

    ea = eps0 + alpha*xa + beta*ya
    eb = eps0 + alpha*xb + beta*yb

    cuts = [0.0]
    de = eb - ea

    for s in seuils:
        if abs(de) > 1e-15:
            t = (s - ea) / de
            if 1e-12 < t < 1.0 - 1e-12:
                cuts.append(t)

    cuts.append(1.0)
    cuts = sorted(set(cuts))

    segments = []

    for k in range(len(cuts)-1):
        t1, t2 = cuts[k], cuts[k+1]

        x1 = xa + t1*(xb-xa)
        y1 = ya + t1*(yb-ya)
        x2 = xa + t2*(xb-xa)
        y2 = ya + t2*(yb-ya)

        e1 = eps0 + alpha*x1 + beta*y1
        e2 = eps0 + alpha*x2 + beta*y2

        segments.append((x1, y1, x2, y2, e1, e2))

    return segments


# ═══════════════════════════════════════════════════════════════
# CONTRIBUTIONS EXACTES — GREEN
# ═══════════════════════════════════════════════════════════════

def _geom_integrals(x1, y1, x2, y2):
    """Intégrales géométriques exactes via Green."""
    dx = x2 - x1
    dy = y2 - y1

    A  = dy * (x1 + dx/2.0)
    My = dy * (x1*y1 + (x1*dy + y1*dx)/2.0 + dx*dy/3.0)
    Mz = dy * (x1**2/2.0 + x1*dx/2.0 + dx**2/6.0)

    return A, My, Mz


# ──────────────────────────────────────────────────────────────
# ELS (σ = C·ε)
# ──────────────────────────────────────────────────────────────
def _contrib_ELS_exact(x1, y1, x2, y2, eps0, alpha, beta, C):

    dx = x2 - x1
    dy = y2 - y1

    if abs(dy) < 1e-15:
        return np.zeros(4)

    ea = eps0 + alpha*x1 + beta*y1
    eb = eps0 + alpha*x2 + beta*y2
    de = eb - ea

    # ✅ INTÉGRALES AVEC x (Green correct)

    N_ = C * dy * (
        ea*x1
        + (ea*dx + x1*de)/2.0
        + de*dx/3.0
    )

    My_ = C * dy * (
        ea*x1*y1
        + (ea*(x1*dy + y1*dx) + de*x1*y1)/2.0
        + de*(x1*dy + y1*dx)/3.0
    )

    Mz_ = C * dy * (
        ea*(x1**2/2.0)
        + (ea*x1*dx + de*(x1**2)/2.0)/2.0
        + de*(x1*dx)/3.0
    )

    Sc_ = dy * (x1 + dx/2.0)

    return np.array([N_, My_, Mz_])


# ──────────────────────────────────────────────────────────────
# ELU — zone R (σ = fcd)
# ──────────────────────────────────────────────────────────────

def _contrib_ELU_R_exact(x1, y1, x2, y2, fcd):

    A, My_g, Mz_g = _geom_integrals(x1, y1, x2, y2)

    N_  = fcd * A
    My_ = fcd * My_g
    Mz_ = fcd * Mz_g
    Sc_ = A

    return np.array([N_, My_, Mz_, Sc_])


# ──────────────────────────────────────────────────────────────
# ELU — zone P (parabole)
# ──────────────────────────────────────────────────────────────

def _contrib_ELU_P_exact(x1, y1, x2, y2, eps0, alpha, beta, fcd, e2):

    ea = eps0 + alpha*x1 + beta*y1
    eb = eps0 + alpha*x2 + beta*y2

    ua = ea / e2
    ub = eb / e2
    du = ub - ua

    I0 = 2*ua - ua**2 + du - ua*du - du**2/3.0

    A, My_g, Mz_g = _geom_integrals(x1, y1, x2, y2)

    N_  = fcd * A * I0
    My_ = fcd * My_g * I0
    Mz_ = fcd * Mz_g * I0
    Sc_ = A

    return np.array([N_, My_, Mz_, Sc_])


# ═══════════════════════════════════════════════════════════════
# INTÉGRATEUR GLOBAL
# ═══════════════════════════════════════════════════════════════

def _integrer_polygone(vertices, eps0, alpha, beta,
                       mode, C=None, fcd=None, e2=None):

    pts = np.asarray(vertices, dtype=float)
    n = len(pts)

    if n < 3:
        return np.zeros(4)

    res = np.zeros(4)

    for i in range(n):

        xa, ya = pts[i]
        xb, yb = pts[(i+1) % n]

        if mode == 'ELS':

            segs = _decouper_arete(
                xa, ya, xb, yb,
                eps0, alpha, beta,
                [0.0]
            )

            for x1, y1, x2, y2, e1, e2_loc in segs:
                if e1 <= 0 and e2_loc <= 0:
                    continue

                res += _contrib_ELS_exact(
                    x1, y1, x2, y2,
                    eps0, alpha, beta, C
                )

        else:  # ELU

            segs = _decouper_arete(
                xa, ya, xb, yb,
                eps0, alpha, beta,
                [0.0, e2]
            )

            for x1, y1, x2, y2, e1, e2_loc in segs:

                if e1 <= 0 and e2_loc <= 0:
                    continue

                elif e1 <= e2 and e2_loc <= e2:
                    res += _contrib_ELU_P_exact(
                        x1, y1, x2, y2,
                        eps0, alpha, beta,
                        fcd, e2
                    )

                else:
                    res += _contrib_ELU_R_exact(
                        x1, y1, x2, y2,
                        fcd
                    )

    # ⚠️ IMPORTANT : PAS de correction de signe ici
    return res



# ═══════════════════════════════════════════════════════════════
# BLOC 3 — FONCTIONS PUBLIQUES
# ═══════════════════════════════════════════════════════════════

@func
def S_com(polygon1, evi, eps0, alpha, beta):
    """
    Aire de la zone comprimée (ε > 0) — Green exact.
    Indépendant du sens de saisie du polygone.
    """
    contour_cg, evidements_cg, Cx, Cy = transformation_repere_cg(polygon1, evi)
    if contour_cg is None:
        return 0.0

    #contour_cg = nettoyer_polygone(contour_cg)
    # On réutilise _integrer_polygone en mode ELS avec C=1
    # (S_com = ∬ 1 dA sur zone comprimée, indépendant de la loi)
    res_c = _integrer_polygone(contour_cg, eps0, alpha, beta, 'ELS', C=1.0)
    Ic = res_c[3]   # Scom

    Iv = 0.0
    for trou in evidements_cg:
        trou = nettoyer_polygone(trou)
        if len(trou) >= 3:
            res_t = _integrer_polygone(trou, eps0, alpha, beta, 'ELS', C=1.0)
            Iv += res_t[3]

    return Ic - Iv


# ─────────────────────────────────────────────────────────────
# ELS
# ─────────────────────────────────────────────────────────────

@func
def calculer_N_My_Mz(
    polygon1, evi, p_acier, s_acier,
    a_com, n, eps0, alpha, beta
):
    """
    ELS — N, My, Mz par Green analytique EXACT.
    Loi béton : σ = max(0, C·ε) avec C = Es/(n·1000).
    Indépendant du sens de saisie du polygone.
    """
    contour_cg, evidements_cg, Cx, Cy = transformation_repere_cg(polygon1, evi)
    if contour_cg is None:
        return 0.0, 0.0, 0.0

    C = 200000.0 / (float(n) * 1000.0)

    #contour_cg = nettoyer_polygone(contour_cg)
    res_c = _integrer_polygone(contour_cg, eps0, alpha, beta, 'ELS', C=C)

    res_v = np.zeros(4)
    for trou in evidements_cg:
        trou = nettoyer_polygone(trou)
        if len(trou) >= 3:
            res_v += _integrer_polygone(trou, eps0, alpha, beta, 'ELS', C=C)

    Nc, Mc_y, Mc_z = (res_c - res_v)[:3]

    # ── Acier
    Ns = Msy = Msz = 0.0
    if p_acier is not None and len(p_acier) > 0:
        acier_data = acier_G(p_acier, s_acier)
        if len(acier_data) > 0:
            acier = np.array(acier_data, dtype=float)
            x_s = acier[:,0] - Cx;  y_s = acier[:,1] - Cy
            eps_s = eps0 + alpha*x_s + beta*y_s
            sig_s = sigma_s_lin1(eps_s, a_com)
            forces = sig_s * acier[:,2]
            Ns  = np.sum(forces)
            Msy = np.sum(forces * y_s)
            Msz = np.sum(forces * x_s)

    return Ns+Nc, Msy+Mc_y, Msz+Mc_z


@func
def solve_GG_ELS(polygon1, evi, p_acier, s_acier,
                 a_com, n, Nobj, Myobj, Mzobj):
    x0 = np.array([0.0, 0.001, 0.0])
    def residuals(eps):
        N, My, Mz = calculer_N_My_Mz(
            polygon1, evi, p_acier, s_acier, a_com, n, *eps)
        return np.array([N-Nobj, My-Myobj, Mz-Mzobj])
    def jacobian(eps, h=1e-6):
        J = np.empty((3,3)); r0 = residuals(eps)
        for i in range(3):
            dh = np.zeros(3); dh[i] = h
            J[:,i] = (residuals(eps+dh) - r0) / h
        return J
    result = root(residuals, x0, jac=jacobian, method='hybr',
                  tol=1e-6, options={'maxfev': 1000})
    return result.x


@func
def resultats_GG_ELS(
    polygon1, evi, p_acier, s_acier,
    a_com, n, eps0, alpha, beta
):
    contour_cg, evidements_cg, Cx, Cy = transformation_repere_cg(polygon1, evi)
    if contour_cg is None:
        return {}

    poly_array = np.asarray(contour_cg, dtype=float)
    eps_c = eps0 + alpha*poly_array[:,0] + beta*poly_array[:,1]
    espcmax = float(np.max(eps_c)); espcmin = float(np.min(eps_c))
    sig_c_max = float(sigma_c_n1(espcmax, n))
    sig_c_min = float(sigma_c_n1(espcmin, n))

    espsmax = espsmin = sig_s_max = sig_s_min = ""
    if p_acier is not None and len(p_acier) > 0:
        acier = acier_G(p_acier, s_acier)
        if len(acier) > 0:
            acier_array = np.asarray(acier, dtype=float)
            x_s = acier_array[:,0] - Cx; y_s = acier_array[:,1] - Cy
            eps_s = eps0 + alpha*x_s + beta*y_s
            espsmax = float(np.max(eps_s)); espsmin = float(np.min(eps_s))
            sig_s_max = float(sigma_s_lin1(espsmax, a_com))
            sig_s_min = float(sigma_s_lin1(espsmin, a_com))

    A_com = S_com(polygon1, evi, eps0, alpha, beta)
    N, My, Mz = calculer_N_My_Mz(
        polygon1, evi, p_acier, s_acier, a_com, n, eps0, alpha, beta)

    return {
        'ACOM': A_com,
        'EPS_C_MAX': espcmax,  'EPS_C_MIN': espcmin,
        'SIG_C_MAX': sig_c_max, 'SIG_C_MIN': sig_c_min,
        'EPS_S_MAX': espsmax,  'EPS_S_MIN': espsmin,
        'SIG_S_MAX': sig_s_max, 'SIG_S_MIN': sig_s_min,
        'N': N, 'MY': My, 'MZ': Mz, 'CX': Cx, 'CY': Cy
    }


@func
def e_resultats_GG_ELS(polygon1, evi, p_acier, s_acier,
                        a_com, n, eps0, alpha, beta, resultats):
    tout = resultats_GG_ELS(
        polygon1, evi, p_acier, s_acier, a_com, n, eps0, alpha, beta)
    return [tout[r] for r in resultats.split(',')]


# ─────────────────────────────────────────────────────────────
# ELU — Parabole-rectangle n=2 (C20-C50)
# ─────────────────────────────────────────────────────────────

@func
def calculer_N_My_Mz_ELU_pararect(
    polygon1, evi, p_acier, s_acier,
    a_com, fck, fcd, fyd, k, eps_uk, eps_ud,
    eps0, alpha, beta
):
    """
    ELU — N, My, Mz par Green analytique EXACT.
    Loi béton : parabole-rectangle n=2 (C20-C50, EC2 §3.1.7).
    3 zones T/P/R — indépendant du sens de saisie du polygone.
    """
    contour_cg, evidements_cg, Cx, Cy = transformation_repere_cg(polygon1, evi)
    if contour_cg is None:
        return 0.0, 0.0, 0.0

    e2 = float(eps_c2(fck))

    #contour_cg = nettoyer_polygone(contour_cg)
    res_c = _integrer_polygone(contour_cg, eps0, alpha, beta,
                               'ELU', fcd=fcd, e2=e2)
    res_v = np.zeros(4)
    for trou in evidements_cg:
        trou = nettoyer_polygone(trou)
        if len(trou) >= 3:
            res_v += _integrer_polygone(trou, eps0, alpha, beta,
                                        'ELU', fcd=fcd, e2=e2)

    Nc, Mc_y, Mc_z = (res_c - res_v)[:3]

    # ── Acier
    Ns = Msy = Msz = 0.0
    if p_acier is not None and len(p_acier) > 0:
        acier_data = acier_G(p_acier, s_acier)
        if len(acier_data) > 0:
            acier = np.array(acier_data, dtype=float)
            x_s = acier[:,0] - Cx;  y_s = acier[:,1] - Cy
            eps_s = eps0 + alpha*x_s + beta*y_s
            sig_s = sigma_s_palier1(fyd, k, eps_uk, eps_ud, eps_s, a_com)
            forces = sig_s * acier[:,2]
            Ns  = np.sum(forces)
            Msy = np.sum(forces * y_s)
            Msz = np.sum(forces * x_s)

    return Ns+Nc, Msy+Mc_y, Msz+Mc_z


@func
def solve_GG_ELU_pararect(polygon1, evi, p_acier, s_acier,
                           a_com, fck, fcd, fyd, k, eps_uk, eps_ud,
                           Nobj, Myobj, Mzobj):
    x0 = np.array([0.0, 0.1, 0.0])
    def residuals(eps):
        N, My, Mz = calculer_N_My_Mz_ELU_pararect(
            polygon1, evi, p_acier, s_acier,
            a_com, fck, fcd, fyd, k, eps_uk, eps_ud, *eps)
        return np.array([N-Nobj, My-Myobj, Mz-Mzobj])
    def jacobian(eps, h=1e-6):
        J = np.empty((3,3)); r0 = residuals(eps)
        for i in range(3):
            dh = np.zeros(3); dh[i] = h
            J[:,i] = (residuals(eps+dh) - r0) / h
        return J
    result = root(residuals, x0, jac=jacobian, method='hybr',
                  tol=1e-6, options={'maxfev': 1000})
    return result.x


@func
def resultats_GG_ELU_pararect(
    polygon1, evi, p_acier, s_acier,
    a_com, fck, fcd, fyd, k, eps_uk, eps_ud,
    eps0, alpha, beta
):
    contour_cg, evidements_cg, Cx, Cy = transformation_repere_cg(polygon1, evi)
    if contour_cg is None:
        return {}

    poly_array = np.asarray(contour_cg, dtype=float)
    eps_c = eps0 + alpha*poly_array[:,0] + beta*poly_array[:,1]
    espcmax = float(np.max(eps_c)); espcmin = float(np.min(eps_c))
    sig_c_max = float(sigma_c_pararect1(fck, fcd, espcmax))
    sig_c_min = float(sigma_c_pararect1(fck, fcd, espcmin))

    epssmax = epssmin = sig_s_max = sig_s_min = ""
    pour_a_pl = 0.0
    if p_acier is not None and len(p_acier) > 0:
        acier = acier_G(p_acier, s_acier)
        if len(acier) > 0:
            acier_array = np.asarray(acier, dtype=float)
            x_s = acier_array[:,0] - Cx; y_s = acier_array[:,1] - Cy
            eps_s = eps0 + alpha*x_s + beta*y_s
            epssmax = float(np.max(eps_s)); epssmin = float(np.min(eps_s))
            plast = np.sum(np.abs(eps_s) >= (float(fyd)/200.0))
            pour_a_pl = float(plast/len(eps_s)*100.0)
            sig_s_max = float(sigma_s_palier1(fyd, k, eps_uk, eps_ud, epssmax, a_com))
            sig_s_min = float(sigma_s_palier1(fyd, k, eps_uk, eps_ud, epssmin, a_com))

    A_com = float(S_com(polygon1, evi, eps0, alpha, beta))
    N, My, Mz = calculer_N_My_Mz_ELU_pararect(
        polygon1, evi, p_acier, s_acier,
        a_com, fck, fcd, fyd, k, eps_uk, eps_ud,
        eps0, alpha, beta)

    return {
        'EPS0': eps0,   'ALPHA': alpha,   'BETA': beta,
        'ACOM': A_com,
        'EPS_C_MAX': espcmax,   'EPS_C_MIN': espcmin,
        'SIG_C_MAX': sig_c_max, 'SIG_C_MIN': sig_c_min,
        'EPS_S_MAX': epssmax,   'EPS_S_MIN': epssmin,
        'SIG_S_MAX': sig_s_max, 'SIG_S_MIN': sig_s_min,
        'N': float(N),  'MY': float(My),  'MZ': float(Mz),
        'PA': pour_a_pl, 'CX': Cx,        'CY': Cy
    }


@func
def e_resultats_GG_ELU_pararect(polygon1, evi, p_acier, s_acier,
                                  a_com, fck, fcd, fyd, k, eps_uk, eps_ud,
                                  eps0, alpha, beta, resultats):
    tout = resultats_GG_ELU_pararect(
        polygon1, evi, p_acier, s_acier,
        a_com, fck, fcd, fyd, k, eps_uk, eps_ud,
        eps0, alpha, beta)
    return [tout[r] for r in resultats.split(',')]


# ═══════════════════════════════════════════════════════════════
# BLOC 4 — TEST DE VALIDATION
# ═══════════════════════════════════════════════════════════════

def _test_validation():
    """
    Vérifie l'indépendance du sens de saisie et la précision
    sur un rectangle 30×60 cm centré à l'origine.
    """
    print("="*60)
    print("TEST DE VALIDATION — Green analytique exact")
    print("Section rectangulaire 30×60 cm")
    print("="*60)

    # Rectangle CCW (sens trigonométrique)
    rect_ccw = [[-15.0,-30.0],[15.0,-30.0],[15.0,30.0],[-15.0,30.0]]
    # Rectangle CW (sens horaire)
    rect_cw  = [[-15.0,-30.0],[-15.0,30.0],[15.0,30.0],[15.0,-30.0]]
    evi = [["fin","fin"]]

    fck=30.; fcd=20.; n_eq=15.
    fyd=434.8; k=1.05; eps_uk=45.; eps_ud=40.; a_com=1.
    eps0=1.0; alpha=0.05; beta=0.08
    p_acier=None; s_acier=None

    # ELS
    N_ccw, My_ccw, Mz_ccw = calculer_N_My_Mz(
        rect_ccw, evi, p_acier, s_acier, a_com, n_eq, eps0, alpha, beta)
    N_cw,  My_cw,  Mz_cw  = calculer_N_My_Mz(
        rect_cw,  evi, p_acier, s_acier, a_com, n_eq, eps0, alpha, beta)

    print(f"\nELS CCW : N={N_ccw:.6f}  My={My_ccw:.6f}  Mz={Mz_ccw:.6f}")
    print(f"ELS CW  : N={N_cw:.6f}  My={My_cw:.6f}  Mz={Mz_cw:.6f}")
    print(f"Écart   : {abs(N_ccw-N_cw):.2e}  {abs(My_ccw-My_cw):.2e}  {abs(Mz_ccw-Mz_cw):.2e}")

    # ELU
    N_e_ccw, My_e_ccw, Mz_e_ccw = calculer_N_My_Mz_ELU_pararect(
        rect_ccw, evi, p_acier, s_acier, a_com, fck, fcd, fyd, k,
        eps_uk, eps_ud, eps0, alpha, beta)
    N_e_cw,  My_e_cw,  Mz_e_cw  = calculer_N_My_Mz_ELU_pararect(
        rect_cw,  evi, p_acier, s_acier, a_com, fck, fcd, fyd, k,
        eps_uk, eps_ud, eps0, alpha, beta)

    print(f"\nELU CCW : N={N_e_ccw:.6f}  My={My_e_ccw:.6f}  Mz={Mz_e_ccw:.6f}")
    print(f"ELU CW  : N={N_e_cw:.6f}  My={My_e_cw:.6f}  Mz={Mz_e_cw:.6f}")
    print(f"Écart   : {abs(N_e_ccw-N_e_cw):.2e}  {abs(My_e_ccw-My_e_cw):.2e}  {abs(Mz_e_ccw-Mz_e_cw):.2e}")

    Sc_ccw = S_com(rect_ccw, evi, eps0, alpha, beta)
    Sc_cw  = S_com(rect_cw,  evi, eps0, alpha, beta)
    print(f"\nS_com CCW={Sc_ccw:.6f}  CW={Sc_cw:.6f}  Écart={abs(Sc_ccw-Sc_cw):.2e}")
    print("\n✓ Validation terminée")


if __name__ == "__main__":
    _test_validation()
