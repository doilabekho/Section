"""
════════════════════════════════════════════════════════════════════════
 CALCUL BP POST-TENSION — MÉTHODE DU RETOUR À L'ÉTAT 0 (CDS IX.1.2)
 xlwings Lite
════════════════════════════════════════════════════════════════════════

 Méthode INSTANTANE — résumé du manuel CDS §IX.1.2 :

   Étape 1 — ÉTAT PERMANENT (QP)
     Équilibre sous G+P avec loi différée du béton.
     Calcule le plan de déformation (ε_i, α_i, β_i).

   Étape 2 — ÉTAT ZÉRO (décompression)
     On remet béton et aciers passifs à zéro.
     La déformation résiduelle des câbles devient :
                       σ_pm − n · σ_b_pm
       ε_p^(0) = ─────────────────────────
                          E_p · (df_p/dε)|0

     Formule simplifiée (acier linéaire au voisinage de 0) :
       ε_p^(0) = (σ_pm + σ_b_pm) / E_p × 1000
     avec σ_b_pm = n · f_b(ε(y_k, z_k))  contrainte béton au droit du câble

   Étape 3 — ÉTAT FINAL (G+P+Q)
     Équilibre avec loi instantanée sur l'incrément.
     Loi câble modifiée :
       σ_p = ρ · [f_p(ε_p^(0) + ε+) − f_p(ε_p^(0))] + f_p(ε_p^(0))

 Notations (cohérentes avec CDS) :
   ε_i, α_i, β_i : plan de déformation à l'état permanent QP
   ε_p^(0)        : déformation câble à l'état zéro
   ε+             : déformation totale à l'état final (eps0 + α·x + β·y)
   ρ              : facteur de participation (0 ≤ ρ ≤ 1)

════════════════════════════════════════════════════════════════════════
"""

import numpy as np
import scipy.integrate
from scipy.integrate import quad
from scipy.optimize  import root, fsolve, least_squares
from scipy.spatial   import Delaunay
from matplotlib.path import Path
from xlwings         import func

_ES = 200_000.0   # MPa


# ════════════════════════════════════════════════════════════════════════════
# HELPERS INTERNES — intégration et données acier
# ════════════════════════════════════════════════════════════════════════════

def _gauss7():
    pts = np.array([
        [1/3, 1/3],
        [0.1012865073234099, 0.1012865073234099],
        [0.7974269853530873, 0.1012865073234099],
        [0.1012865073234099, 0.7974269853530873],
        [0.4701420641051151, 0.0597158717897698],
        [0.4701420641051151, 0.4701420641051151],
        [0.0597158717897698, 0.4701420641051151],
    ], dtype=float)
    w = np.array([
        0.225,
        0.1259391805448272, 0.1259391805448272, 0.1259391805448272,
        0.1323941527885062, 0.1323941527885062, 0.1323941527885062,
    ], dtype=float)
    return pts.T, w


def _acier_array(p_acier, s_acier):
    # Conversion forcée en tableaux NumPy pour manipuler les dimensions facilement
    p_acier = np.atleast_2d(p_acier)
    s_acier = np.atleast_1d(s_acier)
    
    newacier = []
    # On boucle sur le nombre d'armatures (lignes)
    nx = p_acier.shape[0] 
    
    for i in range(nx):
        # p_acier[i, 0] est le X, p_acier[i, 1] est le Y
        # s_acier[i] est la section en mm² (divisée par 10000 pour cm²)
        newacier.append([p_acier[i, 0], p_acier[i, 1], s_acier[i] / 10000])
        
    return newacier



@func

def acier_G(p_acier, s_acier):
    """
    Construit la liste des armatures acier [x, y, section],
    en supprimant UNIQUEMENT les lignes non remplies.
    Les sections nulles ou négatives sont conservées.
    """

    if p_acier is None or len(p_acier) == 0:
        return None, None, None

    p_acier = np.atleast_2d(p_acier)
    s_acier = np.atleast_1d(s_acier)

    newacier = []

    nx = min(p_acier.shape[0], s_acier.shape[0])

    for i in range(nx):
        try:
            # nettoyage UNIQUEMENT des coordonnées vides
            if p_acier[i, 0] in (None, "", " ") or p_acier[i, 1] in (None, "", " "):
                continue

            x = float(p_acier[i, 0])
            y = float(p_acier[i, 1])
            s = float(s_acier[i])   # peut être 0 ou négatif → CONSERVÉ ✅

            newacier.append([x, y, s / 10000.0])  # mm² → cm²

        except (ValueError, TypeError, IndexError):
            # ignore uniquement les lignes vraiment invalides
            continue

    return newacier


###########################################################"##"
# travailler avec les liste pour évidement et transformer donnée
#######################################################

@func
def n_liste(data):
    """
    Compte le nombre de séparateurs 'fin / fin',
    insensible à la casse et aux espaces.
    """
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
    """
    Retourne la i-ème sous-liste séparée par 'fin / fin'
    (insensible à la casse et aux espaces).
    """
    fins = []

    for idx, row in enumerate(data):
        if len(row) >= 2:
            v0 = str(row[0]).strip().lower()
            v1 = str(row[1]).strip().lower()
            if v0 == "fin" or v1 == "fin":
                fins.append(idx)

    # sécurité
    if i < 1 or i > len(fins):
        return []

    start = 0 if i == 1 else fins[i - 2] + 1
    end = fins[i - 1]

    return data[start:end]

    return data[start:end]
@func
def points_valides(data):
    pts = []
    for row in data:
        if len(row) >= 2:
            try:
                x = float(row[0])
                y = float(row[1])
                pts.append((x, y))
            except:
                pass
    return pts
@func
def centre_gravite_polygone(data):
    pts = points_valides(data)
    n = len(pts)

    if n < 3:
        return None, None

    A = 0.0
    Cx = 0.0
    Cy = 0.0

    for i in range(n):
        x1, y1 = pts[i]
        x2, y2 = pts[(i + 1) % n]
        cross = x1 * y2 - x2 * y1
        A += cross
        Cx += (x1 + x2) * cross
        Cy += (y1 + y2) * cross

    A *= 0.5
    if A == 0:
        return None, None

    Cx /= (6 * A)
    Cy /= (6 * A)

    return Cx, Cy
@func
def aire_et_centre(data):
    pts = points_valides(data)
    n = len(pts)

    if n < 3:
        return 0.0, None, None

    A = 0.0
    Cx = 0.0
    Cy = 0.0

    for i in range(n):
        x1, y1 = pts[i]
        x2, y2 = pts[(i + 1) % n]
        cross = x1 * y2 - x2 * y1
        A += cross
        Cx += (x1 + x2) * cross
        Cy += (y1 + y2) * cross

    A *= 0.5
    if A == 0:
        return 0.0, None, None

    Cx /= (6 * A)
    Cy /= (6 * A)

    return abs(A), Cx, Cy
@func
def aire_centre_inertie_origine(data):
    pts = points_valides(data)
    n = len(pts)

    if n < 3:
        return 0.0, None, None, 0.0, 0.0

    A = 0.0
    Cx = 0.0
    Cy = 0.0
    Ix = 0.0
    Iy = 0.0

    for i in range(n):
        x1, y1 = pts[i]
        x2, y2 = pts[(i + 1) % n]
        cross = x1 * y2 - x2 * y1

        A += cross
        Cx += (x1 + x2) * cross
        Cy += (y1 + y2) * cross

        Ix += (y1**2 + y1*y2 + y2**2) * cross
        Iy += (x1**2 + x1*x2 + x2**2) * cross

    A *= 0.5
    if A == 0:
        return 0.0, None, None, 0.0, 0.0

    Cx /= (6 * A)
    Cy /= (6 * A)
    Ix /= 12
    Iy /= 12

    return np.abs(A), Cx, Cy, np.abs(Ix), np.abs(Iy)
@func
def caracteristiques_mecaniques(contour, evidement):
    A, Cx, Cy, Ix0, Iy0 = aire_centre_inertie_origine(contour)

    if A == 0:
        return 0.0, None, None, 0.0, 0.0

    A_tot = A
    Cx_tot = A * Cx
    Cy_tot = A * Cy
    Ix_tot = Ix0
    Iy_tot = Iy0

    N = n_liste(evidement)

    for i in range(1, N + 1):
        Li = liste_i(evidement, i)
        Ai, cxi, cyi, ixi, iyi = aire_centre_inertie_origine(Li)

        if Ai > 0:
            A_tot -= Ai
            Cx_tot -= Ai * cxi
            Cy_tot -= Ai * cyi
            Ix_tot -= ixi
            Iy_tot -= iyi

    if A_tot == 0:
        return 0.0, None, None, 0.0, 0.0

    CxG = Cx_tot / A_tot
    CyG = Cy_tot / A_tot

    IxG = Ix_tot - A_tot * CyG**2
    IyG = Iy_tot - A_tot * CxG**2

    return np.abs(A_tot), CxG, CyG, np.abs(IxG), np.abs(IyG)
@func
def translation_points(data, Cx, Cy):
    """
    Transforme une liste de points dans le repère centré en (Cx, Cy)
    """
    pts_trans = []

    for row in data:
        if len(row) >= 2:
            try:
                x = float(row[0])
                y = float(row[1])
                pts_trans.append((x - Cx, y - Cy))
            except:
                pass  # ignore lignes vides ou non numériques

    return pts_trans

@func
def transformation_repere_cg(contour, evidement):
    """
    Calcule le centre de gravité global puis transforme les coordonnées
    du contour et des évidements dans le repère centré (Cx, Cy).

    Retour :
    - contour_cg        : liste de points transformés
    - evidements_cg     : liste de listes (vide si N = 0)
    - Cx, Cy            : centre de gravité global
    """

    # --- Calcul du centre de gravité global
    A, Cx, Cy, Ix, Iy = caracteristiques_mecaniques(contour, evidement)

    if Cx is None or Cy is None:
        return None, [], None, None

    # --- Contour transformé
    contour_cg = translation_points(contour, Cx, Cy)

    # --- Évidements transformés (séparés)
    evidements_cg = []
    N = n_liste(evidement)

    if N == 0:
        # Aucun évidement
        return contour_cg, evidements_cg, Cx, Cy

    for i in range(1, N + 1):
        Li = liste_i(evidement, i)
        Li_cg = translation_points(Li, Cx, Cy)
        evidements_cg.append(Li_cg)

    return  contour_cg, evidements_cg, Cx, Cy

@func
def nettoyer_polygone(points, tol=1e-9):
    """
    Supprime le dernier point s'il est égal au premier
    à une tolérance numérique près.
    """
    if len(points) < 2:
        return points

    x1, y1 = points[0]
    x2, y2 = points[-1]

    if abs(x1 - x2) < tol and abs(y1 - y2) < tol:
        return points[:-1]

    return points


# ════════════════════════════════════════════════════════════════════════════
# INTÉGRATEUR POLYGONAL
# ════════════════════════════════════════════════════════════════════════════

@func
def polygone_integrate(f, vertices, tol=1e-9, rtol=1e-9, max_depth=4):
    """Intégration adaptative Gauss 7 pts sur polygone (Delaunay)."""
    if vertices is None or len(vertices) < 3:
        return 0.0

    gpts_T, gw = _gauss7()

    def _gauss(v0, v1, v2, area):
        J0x, J0y = v1[0]-v0[0], v1[1]-v0[1]
        J1x, J1y = v2[0]-v0[0], v2[1]-v0[1]
        X = v0[0] + J0x*gpts_T[0] + J1x*gpts_T[1]
        Y = v0[1] + J0y*gpts_T[0] + J1y*gpts_T[1]
        v = np.asarray(f(X, Y))
        return area * (np.dot(gw, v.T) if v.ndim > 1 else np.dot(gw, v))

    def _tri(v0, v1, v2, atol):
        stack = [(v0, v1, v2, atol, 0)]; tot = None
        while stack:
            v0, v1, v2, atol, d = stack.pop()
            area = 0.5*abs((v1[0]-v0[0])*(v2[1]-v0[1])-(v2[0]-v0[0])*(v1[1]-v0[1]))
            if area < 1e-15: continue
            c = _gauss(v0, v1, v2, area)
            if tot is None: tot = np.zeros_like(c)
            if d >= max_depth: tot += c; continue
            m01, m12, m02 = 0.5*(v0+v1), 0.5*(v1+v2), 0.5*(v0+v2); a4 = area*0.25
            fine = (_gauss(v0,m01,m02,a4)+_gauss(v1,m01,m12,a4)+
                    _gauss(v2,m12,m02,a4)+_gauss(m01,m12,m02,a4))
            err = float(np.max(np.abs(fine-c))); sc = max(float(np.max(np.abs(fine))),1e-30)
            if err <= atol and err <= rtol*sc: tot += fine
            else:
                atol4 = atol*0.25
                stack += [(v0,m01,m02,atol4,d+1),(v1,m01,m12,atol4,d+1),
                          (v2,m12,m02,atol4,d+1),(m01,m12,m02,atol4,d+1)]
        return tot if tot is not None else 0.0

    verts = np.asarray(vertices, float)
    tri = Delaunay(verts); path = Path(verts); tot = None
    for s in tri.simplices:
        v0,v1,v2 = verts[s]
        if path.contains_point((v0+v1+v2)/3.):
            r = _tri(v0,v1,v2,tol)
            if tot is None: tot = np.zeros_like(r)
            tot += r
    return tot if tot is not None else 0.0


# ════════════════════════════════════════════════════════════════════════════
# LOIS DE COMPORTEMENT
# ════════════════════════════════════════════════════════════════════════════

@func
def eps_c2(fck):
    fck = np.asarray(fck, float)
    return np.where(fck <= 50, 2.0, 2.0 + 0.085*(fck-50)**0.53)

@func
def eps_cu2(fck):
    fck = np.asarray(fck, float)
    return np.where(fck <= 50, 3.5, 2.6 + 35*((90-fck)/100)**4)

@func
def eps_n(fck):
    fck = np.asarray(fck, float)
    return np.where(fck <= 50, 2.0, 1.4 + 23.4*((90-fck)/100)**4)

@func
def sigma_c_n1(eps_c, n):
    """Béton fissuré ELS — loi linéaire, traction = 0."""
    return np.maximum(0.0, (_ES/n) * np.asarray(eps_c, float) / 1000.0)

@func
def sigma_c_n2(eps_c, n):
    """Béton fissuré ELS — loi linéaire, sans fissurée"""
    return (_ES/n) * np.asarray(eps_c, float) / 1000.0


@func
def sigma_c_pararect1(fck, fcd, eps_c):
    """Parabole-rectangle EC2 — formulation continue (sans np.where)."""
    eps = np.asarray(eps_c, float)
    e2, n = eps_c2(fck), eps_n(fck)
    t  = np.abs(1.0 - eps/e2)
    mask = (eps + np.abs(eps)) / (2.0*np.abs(eps) + 1e-30)
    return np.minimum(fcd*(1.0 - t**n)*mask, fcd*mask)


@func
def S_com(polygon1, evi, eps0, alpha, beta):

    contour_cg, evidements_cg, Cx, Cy = transformation_repere_cg(polygon1, evi)

    if contour_cg is None:
        return 0.0


    def filter_comprime(x, y):
        eps = eps0 + alpha * x + beta * y
        abs_eps = np.abs(eps)
        return (eps + abs_eps) / (2 * abs_eps + 1e-15)

    # Béton plein
    Ic = polygone_integrate(filter_comprime, contour_cg)

    # Évidements
    Iv = 0.0
    for trou in evidements_cg:
        trou = nettoyer_polygone(trou)
        if len(trou) >= 3:
            Iv += polygone_integrate(filter_comprime, trou)

    return Ic - Iv



# ════════════════════════════════════════════════════════════════════════════
# ÉTAPE 1 — ÉTAT PERMANENT QP
# ════════════════════════════════════════════════════════════════════════════
#
# Équilibre sous G+P (loi béton différée via n = n_diff).
# Les câbles participent avec leur contrainte initiale sig_p (fournie).
# Pas de loi de câble modifiée ici : σ_câble = sig_p (imposé).
#
# Formule CDS (§IX.1.2) :
#   G = ∫ f_b(ε) dAb  +  Σ f_s(ε_si) Asi  +  Σ sig_pk Apk
#
# Le solveur trouve (ε_i, α_i, β_i) tel que N=NQP, My=MYQP, Mz=MZQP.
@func
def calculer_P_iso(polygon1, evi, p_acier, s_acier, a_com, n,
                p_pre, s_pre, sig_p,
                fpd, Ep, kp, eps_ukp, eps_udp,
):
    # --------------------------------------------------
    # 0. Transformation géométrique vers le repère CG
    # --------------------------------------------------
    # transformation_repere_cg calcule Cx, Cy AUTOMATIQUEMENT
    contour_cg, evidements_cg, Cx, Cy = transformation_repere_cg(polygon1, evi)    

    # ── Câbles (contrainte initiale imposée) ───────────────────────────────
    acierp = np.array(acier_G(p_pre, s_pre), dtype=float)
    x_p = acierp[:, 0] - Cx
    y_p = acierp[:, 1] - Cy
    A_p = acierp[:, 2]

    sig_pk = np.asarray(
        [v for v in np.atleast_1d(sig_p) if v not in ("", None)],
        dtype=float
    )


    # IMPORTANT : à l'état QP, les câbles exercent leur force initiale.
    # On ne calcule PAS une déformation variable ; on utilise directement sig_p.
    F_p     = sig_pk * A_p
    Nsp     = F_p.sum()
    Msyp    = (F_p * y_p).sum()
    Mszp    = (F_p * x_p).sum()

    return (float(-Nsp ),
            float(-Msyp ),
            float(-Mszp))

@func
def calculer_QP(polygon1, evi, p_acier, s_acier, a_com, n,
                p_pre, s_pre, sig_p,
                fpd, Ep, kp, eps_ukp, eps_udp,
                eps0, alpha, beta):
    """
    Calcule (N, My, Mz) à l'état permanent QP.

    Béton    : loi linéaire fissuré (module différé via n = Es/Ec_diff)
    Acier passif : loi linéaire
    Câbles   : contrainte imposée sig_p (pas de variation)

    CORRECTION vs code original :
      - Les câbles sont traités avec leur contrainte initiale imposée sig_p
        et non via une loi incrémentale (qui n'a pas de sens à l'état QP).
      - sig_p est un vecteur de contraintes (MPa), positif = traction câble.
    """
    # --------------------------------------------------
    # 0. Transformation géométrique vers le repère CG
    # --------------------------------------------------
    # transformation_repere_cg calcule Cx, Cy AUTOMATIQUEMENT
    contour_cg, evidements_cg, Cx, Cy = transformation_repere_cg(polygon1, evi)
  

    # ── Béton ─────────────────────────────────────────────────────────────
    def f_complet(x, y):
        sig = sigma_c_n1(eps0 + alpha*x + beta*y, n)  #  fissuré
        return np.array([sig, sig*y, sig*x])

    # ---- Béton plein
    res_c = polygone_integrate(f_complet, contour_cg)

    # ---- Béton évidements
    if not evidements_cg:      # cas N = 0
        res_v = np.array([0.0, 0.0, 0.0])
    else:
        res_v = np.zeros(3)
        for trou in evidements_cg:
            if len(trou) >= 3:
                res_v += polygone_integrate(f_complet, trou)

    # Résultantes béton
    Nc, Mc_y, Mc_z = res_c - res_v

    # ── Aciers passifs ─────────────────────────────────────────────────────
    if p_acier is None or len(p_acier) == 0:
        Ns, Msy, Msz = 0.0, 0.0, 0.0
    else:
        acier_data = np.array(acier_G(p_acier, s_acier), dtype=float)

    if acier_data.size == 0:
            Ns, Msy, Msz = 0.0, 0.0, 0.0

    else:
        acier = np.array(acier_G(p_acier, s_acier), dtype=float)

        # ⚠️ translation des aciers vers le repère CG
        x_s = acier[:, 0] - Cx
        y_s = acier[:, 1] - Cy
        areas_s = acier[:, 2]

        eps_s = eps0 + alpha * x_s + beta * y_s
        sig_s = sigma_s_lin1(eps_s, a_com)

        forces_s = sig_s * areas_s
        Ns = np.sum(forces_s)
        Msy = np.sum(forces_s * y_s)
        Msz = np.sum(forces_s * x_s)



    return (float(Ns + Nc ),
            float(Msy + Mc_y ),
            float(Msz  + Mc_z  ))


@func
def solve_GG_QP(polygon1, evi, p_acier, s_acier, a_com, n,
                p_pre, s_pre, sig_p,
                fpd, Ep, kp, eps_ukp, eps_udp,
                NQP, MYQP, MZQP):
    """
    Résout l'équilibre à l'état permanent QP → (ε_i, α_i, β_i).

    CORRECTIONS vs code original :
      - x0 adapté : on initialise avec un état de compression légère
        cohérent avec un état permanent sous précontrainte.
      - Cascade fsolve → root/hybr pour robustesse.
    """
    
    targets = np.array([NQP, MYQP, MZQP], float)

    def resid(ep):
        N, My, Mz = calculer_QP(
            polygon1, evi, p_acier, s_acier, a_com, n,
            p_pre, s_pre, sig_p,
            fpd, Ep, kp, eps_ukp, eps_udp,
            ep[0], ep[1], ep[2])
        return np.array([N, My, Mz]) - targets

    # Point de départ : compression légère + gradient
    x0 = np.array([0.5, 0.0, 1e-4])

    # Essai 1 : fsolve
    try:
        x, info, ier, _ = fsolve(resid, x0, full_output=True)
        if ier == 1 and np.max(np.abs(resid(x))) < 1e-3:
            return x
    except Exception:
        pass

    # Essai 2 : root/hybr
    try:
        sol = root(resid, x0, method='hybr', tol=1e-6)
        if sol.success and np.max(np.abs(sol.fun)) < 1e-3:
            return sol.x
    except Exception:
        pass

    # Essai 3 : root/lm
    sol = root(resid, x0, method='lm')
    return sol.x


# ════════════════════════════════════════════════════════════════════════════
# CALCUL DE ε_p^(0) — DÉFORMATION CÂBLE À L'ÉTAT ZÉRO
# ════════════════════════════════════════════════════════════════════════════
#
# Formule CDS §IX.1.2 :
#
#       σ_pm − n · σ_b_pm
# ε_p^(0) = ──────────────────────
#               df_p/dε |_0
#
# Pour une loi linéaire au voisinage de 0 (Ep constant) :
#   df_p/dε |_0 = Ep/1000
#
# Donc :
#   ε_p^(0) = (σ_pm − n · σ_b_pm) / Ep × 1000
#            = eps_s_lin_p(Ep, σ_pm − n·σ_b_pm)
#
# σ_pm     : contrainte initiale dans le câble (fournie, = sig_p[k])
# σ_b_pm   : contrainte béton au droit du câble à l'état QP
#            = f_b(ε(y_k, z_k)) = sigma_c_n1(ε_QP, n)
#
# NOTE IMPORTANTE sur le signe de σ_b_pm :
#   La formule du manuel est : σ_pm - n·σ_b_pm
#   σ_b_pm est la CONTRAINTE DE COMPRESSION du béton (positive en compression EC2).
#   On la soustrait pour RÉDUIRE la déformation du câble lors de la décompression.
#   Plus le béton est comprimé au droit du câble, plus la déformation résiduelle
#   du câble à l'état 0 est petite.

def _eps_cable_etat_zero(
    eps_i, alpha_i, beta_i,
    x_p, y_p,
    sig_p, n, Ep,
    Cx, Cy
):
    """
    Calcule ε_p^(0) pour tous les câbles (retour à zéro instantané).

    ε_QP = déformation béton au droit du câble à l'état QP
    σ_b_pm = contrainte béton au droit du câble
    ε_p^(0) = (σ_pm − n·σ_b_pm) / Ep
    """

    # coordonnées câbles dans le repère CG
    xpc = np.asarray(x_p, float) - Cx
    ypc = np.asarray(y_p, float) - Cy

    # déformation béton au câble (état QP)
    eps_QP = eps_i + alpha_i * xpc + beta_i * ypc

    # contrainte béton correspondante
    sig_b_pm = sigma_c_n1(eps_QP, n)

    sig_pk = np.asarray(
    [v for v in np.atleast_1d(sig_p) if v not in ("", None)],
    dtype=float
    )

    sig_pm = np.asarray(sig_pk, float)  # tension câble < 0

    # déformation initiale câble ε_p^(0)
    return eps_s_lin_p(Ep, sig_pm - n * sig_b_pm) # sig_pm <0 car tension de cable => -n*sig_b_pm


# ════════════════════════════════════════════════════════════════════════════
# ÉTAPE 3 — ÉTAT FINAL ELS  (G + P + Q)
# ════════════════════════════════════════════════════════════════════════════
#
# Équilibre sous charges totales avec loi instantanée.
# Loi câble modifiée (CDS §IX.1.2) :
#   σ_p = ρ · [f_p(ε_p^(0) + ε+) − f_p(ε_p^(0))] + f_p(ε_p^(0))
#
# où ε+ = eps0 + α·x + β·y  est la déformation totale à l'état final.
#
# CORRECTION vs code original :
#   Le code original utilise eps_inc_p = eps_final - eps_QP (incrément QP→final)
#   Le manuel dit : ε+ = eps_final TOTAL (pas un incrément depuis QP).
#   La loi câble opère sur ε_p^(0) + ε+ (total depuis état 0),
#   pas sur ε_p^(0) + Δε (incrément depuis QP).

@func
def calculer_ELS(polygon1, evi, p_acier, s_acier, a_com, n,
                 p_pre, s_pre, sig_p,
                 fpd, Ep, kp, eps_ukp, eps_udp,
                 NQP, MYQP, MZQP, roh,
                 eps0, alpha, beta):
    """
    Calcule (N, My, Mz) à l'état final ELS — méthode retour à zéro INSTANTANE.

    Béton et aciers passifs : loi instantanée sur déformation TOTALE
    (pas incrémentale depuis QP, car on part de l'état zéro).
    Câbles : loi modifiée CDS §IX.1.2.
    """
    
    # --------------------------------------------------
    # 0. Transformation géométrique vers le repère CG
    # --------------------------------------------------
    contour_cg, evidements_cg, Cx, Cy = transformation_repere_cg(polygon1, evi)

    if contour_cg is None:
        return 0.0, 0.0, 0.0

    # ── Étape 1 : état permanent QP ───────────────────────────────────────
    res_qp   = solve_GG_QP(polygon1, evi, p_acier, s_acier, a_com, n,
                            p_pre, s_pre, sig_p, fpd, Ep, kp,
                            eps_ukp, eps_udp, NQP, MYQP, MZQP)
    eps_i, alpha_i, beta_i = float(res_qp[0]), float(res_qp[1]), float(res_qp[2])

    acierp = np.array(acier_G(p_pre, s_pre), dtype=float)
    x_p = acierp[:, 0]
    y_p = acierp[:, 1]
    A_p = acierp[:, 2]


    # ── Étape 2 : déformation câble à l'état zéro ε_p^(0) ────────────────
    eps_p0 = _eps_cable_etat_zero(eps_i, alpha_i, beta_i, x_p, y_p, sig_p, n, Ep, Cx, Cy)

    # ── Étape 3a : béton — loi instantanée sur déformation totale ─────────
    # À l'état final, le béton part de zéro (décompression faite) et
    # subit la déformation totale eps+ = eps0 + α·x + β·y.
    # ── Béton ─────────────────────────────────────────────────────────────
   
    def f_complet(x, y):
        sig = sigma_c_n1(eps0 + alpha*x + beta*y, n)  #  fissuré
        return np.array([sig, sig*y, sig*x])

    # ---- Béton plein
    res_c = polygone_integrate(f_complet, contour_cg)

    # ---- Béton évidements
    if not evidements_cg:      # cas N = 0
        res_v = np.array([0.0, 0.0, 0.0])
    else:
        res_v = np.zeros(3)
        for trou in evidements_cg:
            if len(trou) >= 3:
                res_v += polygone_integrate(f_complet, trou)

    # Résultantes béton
    Nc, Mc_y, Mc_z = res_c - res_v

    # ── Étape 3b : aciers passifs — loi instantanée sur déformation totale ─
    if p_acier is None or len(p_acier) == 0:
        Ns, Msy, Msz = 0.0, 0.0, 0.0
    else:
        acier_data = np.array(acier_G(p_acier, s_acier), dtype=float)

    if acier_data.size == 0:
            Ns, Msy, Msz = 0.0, 0.0, 0.0

    else:
        acier = np.array(acier_G(p_acier, s_acier), dtype=float)
        
        x_s = acier[:, 0] - Cx
        y_s = acier[:, 1] -  Cy
        A_s = acier[:, 2]
        eps_s_final = eps0 + alpha*x_s + beta*y_s
        sig_s       = sigma_s_lin1(eps_s_final, a_com)
        F_s  = sig_s * A_s
        Ns   = F_s.sum();  Msy = (F_s*y_s).sum();  Msz = (F_s*x_s).sum()

    # ── Étape 3c : câbles — loi modifiée CDS §IX.1.2 ──────────────────────
    roh = np.asarray(
    [v for v in np.atleast_1d(roh) if v not in ("", None)],
    dtype=float)

    x_pc = x_p - Cx
    y_pc = y_p - Cy

    eps_plus   = eps0 + alpha*x_pc + beta*y_pc   # déformation totale à l'état final

    # f_p(ε_p^(0))              : contrainte câble à l'état zéro
    sig_p0_val = sigma_s_lin1_p(Ep, eps_p0, a_com)

    # f_p(ε_p^(0) + ε+)         : contrainte câble état final
    sig_p1_val = sigma_s_lin1_p(Ep, eps_p0 + eps_plus, a_com)

    # Loi modifiée : σ_p = ρ·(f_p1 − f_p0) + f_p0
    sig_pk = roh * (sig_p1_val - sig_p0_val) + sig_p0_val

    F_p   = sig_pk * A_p
    #F_p   = roh * sig_p1_val * A_p # test
    Nsp   = F_p.sum();  Msyp = (F_p*y_pc).sum();  Mszp = (F_p*x_pc).sum()
    # Effet isostatique
    Niso, Myiso, Mziso = calculer_P_iso(polygon1, evi, p_acier, s_acier, a_com, n,
            p_pre, s_pre, sig_p,
            fpd, Ep, kp, eps_ukp, eps_udp)
    coef =1
    return (float(Ns +  Nsp + Nc + coef * Niso),
            float(Msy +  Msyp + Mc_y + coef * Myiso),
            float(Msz +   Mszp + Mc_z + coef * Mziso ))


@func
def solve_GG_ELS(polygon1, evi, p_acier, s_acier, a_com, n,
                 p_pre, s_pre, sig_p, fpd, Ep, kp, eps_ukp, eps_udp,
                 NQP, MYQP, MZQP, roh,
                 Nobj, Myobj, Mzobj):
    """
    Résout (N=Nobj, My=Myobj, Mz=Mzobj) → (ε₀, α, β) — ELS.
    Cascade : fsolve → hybr → lm.
    """
    targets = np.array([Nobj, Myobj, Mzobj], float)
    x0      = np.array([0.0, 1e-6, 1e-6], float)

    def resid(ep):
        N, My, Mz = calculer_ELS(
            polygon1, evi, p_acier, s_acier, a_com, n,
            p_pre, s_pre, sig_p, fpd, Ep, kp, eps_ukp, eps_udp,
            NQP, MYQP, MZQP, roh, ep[0], ep[1], ep[2])
        return np.array([N, My, Mz]) - targets

    for solver_fn in (
        lambda: fsolve(resid, x0, full_output=True)[0],
        lambda: root(resid, x0, method='hybr', tol=1e-6).x,
        lambda: root(resid, x0, method='lm').x,
    ):
        try:
            x = solver_fn()
            if np.max(np.abs(resid(x))) < 1e-3:
                return x
        except Exception:
            pass

    return root(resid, x0, method='lm').x


# ════════════════════════════════════════════════════════════════════════════
# RÉSULTATS ELS
# ════════════════════════════════════════════════════════════════════════════

@func
def resultats_GG_ELS(polygon1, evi, p_acier, s_acier, a_com, n,
                     p_pre, s_pre, sig_p, fpd, Ep, kp, eps_ukp, eps_udp,
                     NQP, MYQP, MZQP, roh,
                     eps0, alpha, beta):
    """Retourne un dict de résultats ELS pour le plan (ε₀, α, β)."""
    # --------------------------------------------------
    # 0. Transformation géométrique vers le repère CG
    # --------------------------------------------------
    contour_cg, evidements_cg, Cx, Cy = transformation_repere_cg(polygon1, evi)

    if contour_cg is None:
        return {}

    # --------------------------------------------------
    # 1. Vérification contraintes béton (repère CG)
    # --------------------------------------------------
    if len(contour_cg) > 0:
        poly_array = np.asarray(contour_cg, dtype=float)

        eps_c = eps0 + alpha * poly_array[:, 0] + beta * poly_array[:, 1]
        espcmax = float(np.max(eps_c))
        espcmin = float(np.min(eps_c))

        sig_c_max = float(sigma_c_n1(espcmax, n))
        sig_c_min = float(sigma_c_n1(espcmin, n))
    else:
        espcmax = espcmin = 0.0
        sig_c_max = sig_c_min = 0.0

    # --------------------------------------------------
    # 2. Vérification contraintes acier (repère CG)
    # --------------------------------------------------
    if p_acier is None or len(p_acier) == 0:
        espsmax = espsmin = ""
        sig_s_max = sig_s_min = ""
    else:
        acier = acier_G(p_acier, s_acier)

        if len(acier) == 0:
            espsmax = espsmin = ""
            sig_s_max = sig_s_min = ""
        else:
            acier_array = np.asarray(acier, dtype=float)

            # coordonnées déjà TRANSLATÉES vers CG
            x_s = acier_array[:, 0] - Cx
            y_s = acier_array[:, 1] - Cy

            eps_s = eps0 + alpha * x_s + beta * y_s
            espsmax = float(np.max(eps_s))
            espsmin = float(np.min(eps_s))

            sig_s_max = float(sigma_s_lin1(espsmax, a_com))
            sig_s_min = float(sigma_s_lin1(espsmin, a_com))

    # Câbles
    res_qp = solve_GG_QP(polygon1, evi, p_acier, s_acier, a_com, n,
                          p_pre, s_pre, sig_p, fpd, Ep, kp,
                          eps_ukp, eps_udp, NQP, MYQP, MZQP)
    eps_i, alpha_i, beta_i = float(res_qp[0]), float(res_qp[1]), float(res_qp[2])

    acierp = np.array(acier_G(p_pre, s_pre), dtype=float)
    x_p = acierp[:, 0]
    y_p = acierp[:, 1]
    A_p = acierp[:, 2]

    eps_p0    = _eps_cable_etat_zero(eps_i, alpha_i, beta_i, x_p, y_p, sig_p, n, Ep, Cx, Cy)

    x_pc = x_p - Cx
    y_pc = y_p - Cy
    eps_plus  = eps0 + alpha*x_pc + beta*y_pc
    sig_p0_v  = sigma_s_lin1_p(Ep, eps_p0, a_com)
    sig_p1_v  = sigma_s_lin1_p(Ep, eps_p0 + eps_plus, a_com)
    contrai = sig_p1_v
    #contrai   = roh*(sig_p1_v - sig_p0_v) + sig_p0_v
    espp      = eps_s_palier_p(fpd, Ep, kp, eps_ukp, eps_udp, contrai)

    A_com = S_com(polygon1, evi, eps0, alpha, beta)
    N, My, Mz = calculer_ELS(
        polygon1, evi, p_acier, s_acier, a_com, n,
        p_pre, s_pre, sig_p, fpd, Ep, kp, eps_ukp, eps_udp,
        NQP, MYQP, MZQP, roh, eps0, alpha, beta)

    return {
        "ACOM":      A_com,
        "EPS_C_MAX": espcmax,
        "EPS_C_MIN": espcmin,
        "SIG_C_MAX": sig_c_max,
        "SIG_C_MIN": sig_c_min,
        "EPS_P_MAX": float(espp.max()),  "EPS_P_MIN": float(espp.min()),
        "EPS_S_MAX": espsmax,
        "EPS_S_MIN": espsmin,
        "SIG_S_MAX": sig_s_max,
        "SIG_S_MIN": sig_s_min,
        "SIG_P_MAX": float(contrai.max()), "SIG_P_MIN": float(contrai.min()),
        "SIG_P":     contrai.tolist(),     "EPS_P":     espp.tolist(),
        "N": float(N), "MY": float(My), "MZ": float(Mz),
    }


@func
def e_resultats_GG_ELS(polygon1, evi, p_acier, s_acier, a_com, n,
                        p_pre, s_pre, sig_p, fpd, Ep, kp, eps_ukp, eps_udp,
                        NQP, MYQP, MZQP, roh,
                        eps0, alpha, beta, resultats):
    tout = resultats_GG_ELS(polygon1, evi, p_acier, s_acier, a_com, n,
                            p_pre, s_pre, sig_p, fpd, Ep, kp, eps_ukp, eps_udp,
                            NQP, MYQP, MZQP, roh, eps0, alpha, beta)
    return [tout[r] for r in resultats.split(',')]


# ════════════════════════════════════════════════════════════════════════════
# ÉTAPE 3 — ÉTAT FINAL ELU  (G + P + Q)
# ════════════════════════════════════════════════════════════════════════════
#
# Identique à ELS mais avec :
#   Béton    : loi parabole-rectangle (sigma_c_pararect1)
#   Acier passif : loi palier (sigma_s_palier1)
#   Câbles   : même loi modifiée CDS §IX.1.2, avec f_p = loi palier

@func
def calculer_N_My_Mz_ELU_pararect(
        polygon1, evi, p_acier, s_acier, a_com,
        fck, fcd, fyd, k, eps_uk, eps_ud,
        p_pre, s_pre, sig_p, fpd, Ep, kp, eps_ukp, eps_udp,
        n, NQP, MYQP, MZQP, roh,
        eps0, alpha, beta):
    """
    Calcule (N, My, Mz) ELU — méthode retour à zéro INSTANTANE.

    CORRECTIONS vs code original :
      1. ε+ est la déformation TOTALE (pas un incrément depuis QP).
      2. La loi câble utilise ε_p^(0) + ε+ (pas ε_p^(0) + Δε_QP→final).
      3. σ_b à l'état QP est corrigé : signe − dans la formule ε_p^(0).
    """
    contour_cg, evidements_cg, Cx, Cy = transformation_repere_cg(polygon1, evi)

    if contour_cg is None:
        return 0.0, 0.0, 0.0
    # ── État QP ──────────────────────────────────────────────────────────
    res_qp   = solve_GG_QP(polygon1, evi, p_acier, s_acier, a_com, n,
                            p_pre, s_pre, sig_p, fpd, Ep, kp,
                            eps_ukp, eps_udp, NQP, MYQP, MZQP)
    eps_i, alpha_i, beta_i = float(res_qp[0]), float(res_qp[1]), float(res_qp[2])

    acierp = np.array(acier_G(p_pre, s_pre), dtype=float)
    x_p = acierp[:, 0]
    y_p = acierp[:, 1]
    A_p = acierp[:, 2]

    # ── ε_p^(0) ──────────────────────────────────────────────────────────
    eps_p0 = _eps_cable_etat_zero(eps_i, alpha_i, beta_i, x_p, y_p, sig_p, n, Ep, Cx, Cy)

    # ── Béton ELU — loi parabole-rectangle ───────────────────────────────
    def f_complet(x, y):
        sig = sigma_c_pararect1(fck, fcd, eps0 + alpha*x + beta*y)
        return np.array([sig, sig*y, sig*x])

    # ---- Béton plein
    res_c = polygone_integrate(f_complet, contour_cg)

    # ---- Béton évidements
    if not evidements_cg:      # cas N = 0
        res_v = np.array([0.0, 0.0, 0.0])
    else:
        res_v = np.zeros(3)
        for trou in evidements_cg:
            if len(trou) >= 3:
                res_v += polygone_integrate(f_complet, trou)

    # Résultantes béton
    Nc, Mc_y, Mc_z = res_c - res_v

    # ── Aciers passifs ELU — loi palier ──────────────────────────────────
    if p_acier is None or len(p_acier) == 0:
        Ns, Msy, Msz = 0.0, 0.0, 0.0
    else:
        acier_data = np.array(acier_G(p_acier, s_acier), dtype=float)

    if acier_data.size == 0:
            Ns, Msy, Msz = 0.0, 0.0, 0.0

    else:
        acier = np.array(acier_G(p_acier, s_acier), dtype=float)
        x_s = acier[:, 0] - Cx
        y_s = acier[:, 1] - Cy
        A_s = acier[:, 2]
        eps_s_final = eps0 + alpha*x_s + beta*y_s
        sig_s       = sigma_s_palier1(fyd, k, eps_uk, eps_ud, eps_s_final, a_com)
        F_s  = sig_s * A_s
        Ns   = F_s.sum();  Msy = (F_s*y_s).sum();  Msz = (F_s*x_s).sum()

    # ── Câbles ELU — loi modifiée CDS §IX.1.2 ────────────────────────────
    roh = np.asarray(
    [v for v in np.atleast_1d(roh) if v not in ("", None)],
    dtype=float)
    x_pc = x_p - Cx
    y_pc = y_p - Cy
    eps_plus   = eps0 + alpha*x_pc + beta*y_pc

    sig_p0_val = sigma_s_palier_p1(fpd, Ep, kp, eps_ukp, eps_udp, eps_p0, a_com)
    sig_p1_val = sigma_s_palier_p1(fpd, Ep, kp, eps_ukp, eps_udp, eps_p0 + eps_plus, a_com)
    sig_pk     = roh*(sig_p1_val - sig_p0_val) + sig_p0_val

    F_p   = sig_pk * A_p
    Nsp   = F_p.sum();  Msyp = (F_p*y_pc).sum();  Mszp = (F_p*x_pc).sum()

    Niso, Myiso, Mziso = calculer_P_iso(polygon1, evi, p_acier, s_acier, a_com, n,
            p_pre, s_pre, sig_p,
            fpd, Ep, kp, eps_ukp, eps_udp)
    coef = 1
    return (float(Ns +  Nsp + Nc + coef * Niso),
            float(Msy +  Msyp + Mc_y + coef * Myiso),
            float(Msz +   Mszp + Mc_z + coef * Mziso ))
@func
def solve_GG_ELU_pararect(
        polygon1, evi, p_acier, s_acier, a_com,
        fck, fcd, fyd, k, eps_uk, eps_ud,
        p_pre, s_pre, sig_p, fpd, Ep, kp, eps_ukp, eps_udp,
        n, NQP, MYQP, MZQP, roh,
        Nobj, Myobj, Mzobj):
    """
    Résout (N=Nobj, My=Myobj, Mz=Mzobj) → (ε₀, α, β) — ELU.
    """
    targets = np.array([Nobj, Myobj, Mzobj], float)
    x0      = np.array([0.5, 0.0, 1e-4], float)

    def resid(ep):
        N, My, Mz = calculer_N_My_Mz_ELU_pararect(
            polygon1, evi, p_acier, s_acier, a_com,
            fck, fcd, fyd, k, eps_uk, eps_ud,
            p_pre, s_pre, sig_p, fpd, Ep, kp, eps_ukp, eps_udp,
            n, NQP, MYQP, MZQP, roh, ep[0], ep[1], ep[2])
        return np.array([N, My, Mz]) - targets

    for solver_fn in (
        lambda: fsolve(resid, x0, full_output=True)[0],
        lambda: root(resid, x0, method='hybr', tol=1e-6).x,
        lambda: root(resid, x0, method='lm').x,
    ):
        try:
            x = solver_fn()
            if np.max(np.abs(resid(x))) < 1e-3:
                return x
        except Exception:
            pass

    return root(resid, x0, method='lm').x


# ════════════════════════════════════════════════════════════════════════════
# RÉSULTATS ELU
# ════════════════════════════════════════════════════════════════════════════

@func
def resultats_GG_ELU_pararect(
        polygon1, evi, p_acier, s_acier, a_com,
        fck, fcd, fyd, k, eps_uk, eps_ud,
        p_pre, s_pre, sig_p, fpd, Ep, kp, eps_ukp, eps_udp,
        n, NQP, MYQP, MZQP, roh,
        eps0, alpha, beta):
    """Retourne un dict de résultats ELU pour le plan (ε₀, α, β)."""
    contour_cg, evidements_cg, Cx, Cy = transformation_repere_cg(polygon1, evi)

     # --------------------------------------------------
    # 1. Vérification contraintes béton (repère CG)
    # --------------------------------------------------
    if len(contour_cg) > 0:
        poly_array = np.asarray(contour_cg, dtype=float)

        eps_c = eps0 + alpha * poly_array[:, 0] + beta * poly_array[:, 1]
        espcmax = float(np.max(eps_c))
        espcmin = float(np.min(eps_c))

        sig_c_max = float(sigma_c_pararect1(fck, fcd, espcmax))
        sig_c_min = float(sigma_c_pararect1(fck, fcd, espcmin))
    else:
        espcmax = espcmin = 0.0
        sig_c_max = sig_c_min = 0.0

    # --------------------------------------------------
    # 2. Vérification contraintes acier (repère CG)
    # --------------------------------------------------
    if p_acier is None or len(p_acier) == 0:
        epssmax = epssmin = ""
        sig_s_max = sig_s_min = ""
        pour_a_pl = 0.0
    else:
        acier = acier_G(p_acier, s_acier)

        if len(acier) == 0:
            epssmax = epssmin = ""
            sig_s_max = sig_s_min = ""
            pour_a_pl = 0.0
        else:
            acier_array = np.asarray(acier, dtype=float)

            # coordonnées déjà TRANSLATÉES vers CG
            x_s = acier_array[:, 0] - Cx
            y_s = acier_array[:, 1] - Cy

            eps_s = eps0 + alpha * x_s + beta * y_s
            epssmax = float(np.max(eps_s))
            epssmin = float(np.min(eps_s))

            # Pourcentage plastifié (Es ≈ 200 000 MPa)
            plast = np.sum(np.abs(eps_s) >= (float(fyd) / 200.0))
            pour_a_pl = float(plast / len(eps_s) * 100.0)

            sig_s_max = float(sigma_s_palier1(
                fyd, k, eps_uk, eps_ud, epssmax, a_com
            ))
            sig_s_min = float(sigma_s_palier1(
                fyd, k, eps_uk, eps_ud, epssmin, a_com
            ))

    # Câbles
    res_qp = solve_GG_QP(polygon1, evi, p_acier, s_acier, a_com, n,
                          p_pre, s_pre, sig_p, fpd, Ep, kp,
                          eps_ukp, eps_udp, NQP, MYQP, MZQP)
    eps_i, alpha_i, beta_i = float(res_qp[0]), float(res_qp[1]), float(res_qp[2])

    acierp = np.array(acier_G(p_pre, s_pre), dtype=float)
    x_p = acierp[:, 0]
    y_p = acierp[:, 1]
    A_p = acierp[:, 2]

    eps_p0    = _eps_cable_etat_zero(eps_i, alpha_i, beta_i, x_p, y_p, sig_p, n, Ep, Cx, Cy)

    x_pc = x_p - Cx
    y_pc = y_p - Cy
    eps_plus = eps0 + alpha*x_pc + beta*y_pc
    sig_p0_v = sigma_s_palier_p1(fpd, Ep, kp, eps_ukp, eps_udp, eps_p0, a_com)
    sig_p1_v = sigma_s_palier_p1(fpd, Ep, kp, eps_ukp, eps_udp, eps_p0 + eps_plus, a_com)
    contrai  = sig_p1_v
    #contrai  = roh*(sig_p1_v - sig_p0_v) + sig_p0_v
    espp     = eps_s_palier_p(fpd, Ep, kp, eps_ukp, eps_udp, contrai)

    # Plastification acier passif


    A_com = S_com(polygon1, evi, eps0, alpha, beta)
    N, My, Mz = calculer_N_My_Mz_ELU_pararect(
        polygon1, evi, p_acier, s_acier, a_com,
        fck, fcd, fyd, k, eps_uk, eps_ud,
        p_pre, s_pre, sig_p, fpd, Ep, kp, eps_ukp, eps_udp,
        n, NQP, MYQP, MZQP, roh, eps0, alpha, beta)

    return {
        "EPS0": eps0, "ALPHA": alpha, "BETA": beta, "ACOM": A_com,
        "EPS_C_MAX": espcmax,
        "EPS_C_MIN": espcmin,
        "SIG_C_MAX": sig_c_max,
        "SIG_C_MIN": sig_c_min,
        "EPS_P_MAX": float(espp.max()),  "EPS_P_MIN": float(espp.min()),
        "EPS_S_MAX": epssmax,
        "EPS_S_MIN": epssmin,
        "SIG_S_MAX": sig_s_max,
        "SIG_S_MIN": sig_s_min,
        "SIG_P_MAX": float(contrai.max()), "SIG_P_MIN": float(contrai.min()),
        "SIG_P":     contrai.tolist(),     "EPS_P":     espp.tolist(),
        "N": float(N), "MY": float(My), "MZ": float(Mz),
        "PA": pour_a_pl,
    }


@func
def e_resultats_GG_ELU_pararect(
        polygon1, evi, p_acier, s_acier, a_com,
        fck, fcd, fyd, k, eps_uk, eps_ud,
        p_pre, s_pre, sig_p, fpd, Ep, kp, eps_ukp, eps_udp,
        n, NQP, MYQP, MZQP, roh,
        eps0, alpha, beta, resultats):
    tout = resultats_GG_ELU_pararect(
        polygon1, evi, p_acier, s_acier, a_com,
        fck, fcd, fyd, k, eps_uk, eps_ud,
        p_pre, s_pre, sig_p, fpd, Ep, kp, eps_ukp, eps_udp,
        n, NQP, MYQP, MZQP, roh, eps0, alpha, beta)
    return [tout[r] for r in resultats.split(',')]
