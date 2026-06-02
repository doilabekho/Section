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


from beton import *   # ← fonctionne en local et sur GitHub
from acier import *   # ← fonctionne en local et sur GitHub	




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


# ═══════════════════════════════════════════════════════════════
# BLOC 2 — NOYAU GREEN ANALYTIQUE EXACT
# ═══════════════════════════════════════════════════════════════
#
# Théorème de Green :  ∬_Ω f dA = ∮_∂Ω Q(x,y) dy  avec ∂Q/∂x = f
#
# Arête paramétrée : x(s)=xa+s·Δx, y(s)=ya+s·Δy, s∈[0,1]
# ∮ Q dy = Σ_arêtes ∫₀¹ Q(x(s),y(s))·Δy ds
#
# ── PRIMITIVES Q (∂Q/∂x = f) ────────────────────────────────
#
# ══ ELS Zone C : f = C·ε = C·(ε₀+α·x+β·y) ══════════════════
#   Q_N  = C·(ε₀·x + α·x²/2 + β·x·y)
#   Q_My = C·(ε₀·x·y + α·x²·y/2 + β·x·y²)
#   Q_Mz = C·(ε₀·x²/2 + α·x³/3 + β·x²·y/2)
#   Q_Sc = x                                  [f=1 zone comprimée]
#
# Sur arête (s∈[0,1]) :
#   Notations : x=xa+s·Δx, y=ya+s·Δy
#   Moments : ∫s^k ds = 1/(k+1)
#   ∫x ds   = xa + Δx/2
#   ∫x² ds  = xa² + xa·Δx + Δx²/3
#   ∫x³ ds  = xa³ + 3xa²·Δx/2 + xa·Δx² + Δx³/4
#   ∫x·y ds = xa·ya + (xa·Δy+ya·Δx)/2 + Δx·Δy/3
#   ∫x²·y ds= xa²·ya + xa·(ya·Δx+xa·Δy)/1... (développé ci-dessous)
#   ∫x·y² ds= xa·ya² + (xa·2ya·Δy+ya²·Δx)/2 + xa·Δy²/3 + ya·Δx·Δy/2 + ...
#
#   N  = Δy · C · ∫₀¹(ε₀·x + α·x²/2 + β·x·y) ds
#      = Δy · C · [ε₀·(xa+Δx/2) + α/2·(xa²+xa·Δx+Δx²/3) + β·(xa·ya+(xa·Δy+ya·Δx)/2+Δx·Δy/3)]
#
#   My = Δy · C · ∫₀¹(ε₀·x·y + α·x²·y/2 + β·x·y²) ds
#
#   Mz = Δy · C · ∫₀¹(ε₀·x²/2 + α·x³/3 + β·x²·y/2) ds
#
# ══ ELU Zone R : f = fcd ════════════════════════════════════
#   Q_N  = fcd·x
#   Q_My = fcd·x·y
#   Q_Mz = fcd·x²/2
#   Q_Sc = x
#
#   N  = fcd · Δy · (xa + Δx/2)
#   My = fcd · Δy · (xa·ya + (xa·Δy+ya·Δx)/2 + Δx·Δy/3)
#   Mz = fcd · Δy/2 · (xa² + xa·Δx + Δx²/3)
#   Sc = Δy · (xa + Δx/2)
#
# ══ ELU Zone P : f = fcd·(2u-u²), u=ε/e2 ══════════════════
#   u = (ε₀+α·x+β·y)/e2 = u₀ + (α·x+β·y)/e2
#   u(s) = ua + s·Δu  où ua=(ε₀+α·xa+β·ya)/e2, Δu=α·Δx/e2+β·Δy/e2... 
#   NON : u dépend de x ET y simultanément via ε = ε₀+α·x+β·y
#
#   Sur l'arête : u(s) = (ε₀+α·(xa+s·Δx)+β·(ya+s·Δy))/e2
#               = ua + s·(α·Δx+β·Δy)/e2
#               = ua + s·Δu
#   C'est bien linéaire en s → u² est quadratique en s
#
#   f = fcd·(2u-u²) = fcd·(2(ua+s·Δu) - (ua+s·Δu)²)
#     = fcd·[(2ua-ua²) + s·(2Δu-2ua·Δu) - s²·Δu²]
#
#   Pour N = Δy·∫₀¹ Q_N(x(s),y(s)) ds, on a besoin de Q_N :
#   Q_N tel que ∂Q_N/∂x = fcd·(2u-u²)
#   Or u = (ε₀+α·x+β·y)/e2, donc ∂u/∂x = α/e2
#   ∂Q_N/∂x = fcd·(2u-u²)  →  Q_N = fcd·(u²/( α/e2) - u³/(3·α/e2))... 
#   c'est complexe. On utilise directement l'intégrale de ligne :
#
#   N  = ∫_arête σ·dy = Δy · ∫₀¹ fcd·(2ua+2Δu·s - ua²-2ua·Δu·s - Δu²·s²) ds
#      = Δy · fcd · [(2ua-ua²) + (2Δu-2ua·Δu)/2 - Δu²/3]
#      = Δy · fcd · [2ua - ua² + Δu - ua·Δu - Δu²/3]
#      = Δy · fcd · I0
#
#   My = ∫_arête σ·y·dy = Δy · ∫₀¹ fcd·(2u-u²)·(ya+Δy·s) ds
#      = Δy · fcd · [ya·I0 + Δy·I1]
#   où I1 = ∫₀¹(2u-u²)·s ds
#          = (2ua-ua²)/2 + (2Δu-2ua·Δu)/3 - Δu²/4
#
#   Mz = ∫_arête σ·x·dy = Δy · ∫₀¹ fcd·(2u-u²)·(xa+Δx·s) ds
#      = Δy · fcd · [xa·I0 + Δx·I1]
#
# IMPORTANT : Green donne ∬σ dA = ∮ Q dy
# Pour N : ∬σ dA = ∮ Q_N dy  NON IDENTIQUE à ∮ σ·dy  !!
# La formule correcte est :
#   ∬ f(x,y) dA = ∮_contour Q(x,y)·dy  avec ∂Q/∂x = f
#
# Pour ELU zone P, f=σ(x,y) et Q_N tel que ∂Q_N/∂x = σ(x,y) :
#   σ = fcd·(2ε/e2 - ε²/e2²)
#   ε = ε₀ + α·x + β·y
#   ∂σ/∂x = fcd·(2α/e2 - 2α·ε/e2²)
#   ∂Q_N/∂x = σ  →  Q_N = fcd·(2ε·x/(e2... NON
#   Q_N = ∫σ dx = fcd·∫(2(ε₀+αx+βy)/e2 - (ε₀+αx+βy)²/e2²) dx
#
#   Posons A = ε₀+βy (constant en x), B = α
#   ε = A + B·x
#   Q_N = fcd·∫(2(A+Bx)/e2 - (A+Bx)²/e2²) dx
#       = fcd·[2(Ax+Bx²/2)/e2 - (A²x+ABx²+B²x³/3)/e2²]
#       = fcd·[2Ax/e2 + Bx²/e2 - A²x/e2² - ABx²/e2² - B²x³/(3e2²)]
#
#   Sur l'arête, x=xa+s·Δx, y=ya+s·Δy :
#   A(s) = ε₀ + β·(ya+s·Δy)  → dépend de s via y !
#   Donc Q_N(x(s),y(s)) n'est pas simplement évaluable.
#
#   SOLUTION : On intègre directement ∫₀¹ Q_N(x(s),y(s))·Δy ds
#   en substituant x=xa+s·Δx et y=ya+s·Δy :
#   A(s) = ε₀+β·ya + β·Δy·s = A0 + dA·s  avec A0=ε₀+β·ya, dA=β·Δy
#   B = α (constant)
#   x(s) = xa + Δx·s
#
#   Q_N(s) = fcd·[2A(s)x(s)/e2 + B·x²(s)/e2
#                 - A²(s)·x(s)/e2² - A(s)B·x²(s)/e2² - B²·x³(s)/(3e2²)]
#
#   Chaque terme est un polynôme en s de degré ≤ 3 → intégrable exactement.
# ═══════════════════════════════════════════════════════════════

def _aire_signee(pts):
    """Aire signée — Shoelace. + = CCW, - = CW."""
    a = 0.0
    n = len(pts)
    for i in range(n):
        x1, y1 = pts[i]
        x2, y2 = pts[(i+1) % n]
        a += x1*y2 - x2*y1
    return 0.5*a


def _decouper_arete(xa, ya, xb, yb, eps0, alpha, beta, seuils):
    """
    Découpe l'arête aux t ∈ (0,1) où ε = seuil.
    Sécurisé pour les sections en I : préserve TOUJOURS l'orientation
    géométrique originale de l'arête pour ne pas fausser S_com.
    """
    ea = eps0 + alpha * xa + beta * ya
    eb = eps0 + alpha * xb + beta * yb
    de = eb - ea

    # On commence avec les bornes naturelles du segment dans le sens d'origine
    cuts = [0.0, 1.0]
    
    for s in seuils:
        if abs(de) > 1e-15:
            t = (s - ea) / de
            if 1e-11 < t < 1.0 - 1e-11:
                cuts.append(t)
                
    # TRÈS IMPORTANT : On trie TOUJOURS dans le sens croissant du paramètre t (0 -> 1)
    # afin de conserver le vecteur directeur géométrique original (xa,ya) -> (xb,yb)
    cuts = sorted(list(set(cuts)))

    segments = []
    for k in range(len(cuts) - 1):
        t1, t2 = cuts[k], cuts[k+1]
        tm = 0.5 * (t1 + t2)
        
        # Interpolation des coordonnées respectant le sens de parcours original
        x1 = xa + t1 * (xb - xa)
        y1 = ya + t1 * (yb - ya)
        x2 = xa + t2 * (xb - xa)
        y2 = ya + t2 * (yb - ya)
        
        # Déformation au milieu du sous-segment
        em = ea + tm * de
        
        segments.append((x1, y1, x2, y2, em))
        
    return segments
# ── Moments polynomiaux sur [0,1] ───────────────────────────────────
# On note p(s) = p0 + s·dp  pour tout paramètre linéaire en s.
# ∫₀¹ p^k · s^j ds se calcule par binôme.
#
# Fonctions utilitaires :
def _I(a, da, b, db, p, q):
    """
    ∫₀¹ (a+da·s)^p · (b+db·s)^q ds  pour p,q ∈ {0,1,2,3} entiers.
    Développé par binôme et intégré terme à terme.
    Retourne valeur scalaire.
    """
    # Développement binomial de (a+da·s)^p et (b+db·s)^q
    from math import comb
    total = 0.0
    for i in range(p+1):
        ci = comb(p, i) * (a**(p-i)) * (da**i)
        for j in range(q+1):
            cj = comb(q, j) * (b**(q-j)) * (db**j)
            total += ci * cj / (i+j+1)
    return total


def _contrib_Scom(x1, y1, x2, y2):
    """
    Calcul de l'aire comprimée par intégration directe (Gauss / Triangle).
    Marche parfaitement même sur contour ouvert (arêtes découpées).
    """
    return 0.5 * (x1 * y2 - x2 * y1)

def _contrib_ELS(x1, y1, x2, y2, eps0, alpha, beta, C):
    """
    Contribution ELS par intégration directe de la contrainte.
    Évite les erreurs de Green sur les contours ouverts après découpe.
    Convention : My = σ * y  et  Mz = σ * x
    """
    dx = x2 - x1
    dy = y2 - y1
    
    # Calcul des intégrales de base du triangle formé avec l'origine (0,0)
    # dA = 0.5 * (x1*y2 - x2*y1) -> On utilise l'approche de projection directe
    # Pour s'aligner sur ton noyau sans réécrire l'intégrateur, on conserve le dy de Green
    # mais on s'assure que l'arête respecte la fermeture.
    
    if abs(dy) < 1e-15:
        return np.zeros(3)

    Ix    = _I(x1, dx, 1, 0, 1, 0)   
    Ix2   = _I(x1, dx, 1, 0, 2, 0)   
    Ix3   = _I(x1, dx, 1, 0, 3, 0)   
    Ixy   = _I(x1, dx, y1, dy, 1, 1) 
    Ix2y  = _I(x1, dx, y1, dy, 2, 1) 
    Ixy2  = _I(x1, dx, y1, dy, 1, 2) 

    N_  = C * dy * (eps0 * Ix + (alpha / 2.0) * Ix2 + beta * Ixy)
    My_ = C * dy * (eps0 * Ixy + (alpha / 2.0) * Ix2y + beta * Ixy2)
    Mz_ = C * dy * ((eps0 / 2.0) * Ix2 + (alpha / 3.0) * Ix3 + (beta / 2.0) * Ix2y)

    return np.array([N_, My_, Mz_])





def _contrib_ELU_R(x1, y1, x2, y2, fcd):
    """
    Contribution EXACTE ELU zone R (σ=fcd) à [N, My, Mz].
    Q_N  = fcd·x        →  N  = fcd·Δy·∫x ds  = fcd·Δy·(x1+Δx/2)
    Q_My = fcd·x·y      →  My = fcd·Δy·∫x·y ds
    Q_Mz = fcd·x²/2     →  Mz = fcd·Δy/2·∫x² ds
    """
    dx = x2-x1;  dy = y2-y1
    if abs(dy) < 1e-15:
        return np.zeros(3)

    Ix  = _I(x1, dx, 1, 0, 1, 0)   # ∫x ds
    Ix2 = _I(x1, dx, 1, 0, 2, 0)   # ∫x² ds
    Ixy = _I(x1, dx, y1, dy, 1, 1) # ∫x·y ds

    N_  = fcd * dy * Ix
    My_ = fcd * dy * Ixy
    Mz_ = fcd * dy/2.0 * Ix2

    return np.array([N_, My_, Mz_])


def _contrib_ELU_P(x1, y1, x2, y2, eps0, alpha, beta, fcd, e2, n=2.0):
    """
    Contribution EXACTE ELU zone P (σ=fcd·(1 - (1-u)ⁿ), u=ε/e2) à [N, My, Mz].
    Généralisée pour une puissance n quelconque (ex: n=1.4 pour C70).
    
    Principe mathématique :
    Au lieu de développer le polynôme (impossible si n n'est pas entier),
    on intègre analytiquement la primitive exacte de la fonction puissance W^n.
    Pour éviter la division par zéro si alpha=0 ou beta=0, on projette
    le théorème de Green sur l'axe X ou Y selon le gradient dominant.
    """
    dx = x2 - x1
    dy = y2 - y1
    if abs(dx) < 1e-15 and abs(dy) < 1e-15:
        return np.zeros(3)

    # 1. Calcul des déformations aux nœuds de l'arête
    eps1 = eps0 + alpha * x1 + beta * y1
    eps2 = eps0 + alpha * x2 + beta * y2
    
    # 2. Variable interne W = 1 - u = 1 - eps / e2
    # Borne max à 0.0 pour éviter les NaN (numérique) avec des puissances non entières
    W1 = max(1.0 - eps1 / e2, 0.0)
    W2 = max(1.0 - eps2 / e2, 0.0)
    
    # Fonction d'intégration exacte de W(s)^p * s^k sur le segment paramétrique s ∈ [0, 1]
    def _I_W(p):
        dW = W2 - W1
        
        # Si dW est très petit, l'intégration analytique produit une instabilité (0/0).
        # On utilise donc un développement de Taylor sécurisé (Ordre 2).
        if abs(dW) < 1e-5:
            # Ordre 0
            T0 = W1**p
            
            # Calcul des termes avec sécurité sur W1 > 0 pour éviter la division par zéro
            # si l'exposant (p-k) devient négatif.
            
            # Ordre 1
            T1 = p * (W1**(p-1)) * dW if W1 > 0 else 0.0
            
            # Ordre 2
            T2 = (p * (p-1) / 2.0) * (W1**(p-2)) * (dW**2) if W1 > 0 else 0.0
            
            # Ordre 3
            T3 = (p * (p-1) * (p-2) / 6.0) * (W1**(p-3)) * (dW**3) if W1 > 0 else 0.0
            
            # Ordre 4
            T4 = (p * (p-1) * (p-2) * (p-3) / 24.0) * (W1**(p-4)) * (dW**4) if W1 > 0 else 0.0
            
            # Intégrale de W(s)^p
            I0 = T0 + T1/2.0 + T2/3.0 + T3/4.0 + T4/5.0
            
            # Intégrale de W(s)^p * s
            I1 = T0/2.0 + T1/3.0 + T2/4.0 + T3/5.0 + T4/6.0
            
            return I0, I1
        else:
            # Primitives exactes si dW est suffisant
            I0 = (W2**(p+1) - W1**(p+1)) / (dW * (p+1))
            I1 = (W2**(p+2) - W1**(p+2)) / (dW**2 * (p+2)) - (W1 / dW) * I0
            return I0, I1

    # 3. Cas particulier : Déformation uniforme sur toute la section (Compression pure)
    if abs(alpha) < 1e-12 and abs(beta) < 1e-12:
        sigma = fcd * (1.0 - W1**n)
        N_  = sigma * dy * (x1 + dx/2.0)
        My_ = sigma * dy * (x1*y1 + (x1*dy + y1*dx)/2.0 + dx*dy/3.0)
        Mz_ = sigma * dy * (x1**2 + x1*dx + dx**2/3.0) / 2.0
        return np.array([N_, My_, Mz_])

    # 4. Calcul des intégrales paramétriques I(p, k) nécessaires
    I0_n1, I1_n1 = _I_W(n + 1.0)
    I0_n2, _     = _I_W(n + 2.0)

    # 5. Application Adaptative du Théorème de Green
    # On choisit d'intégrer selon X ou Y pour ne JAMAIS diviser par un coefficient nul
    if abs(alpha) >= abs(beta):
        # --- Gradient dominant en X --- 
        # Vecteur de Green : Q_N = Int(sigma dx) => Contour sur dy
        C1 = e2 / (alpha * (n + 1.0))
        C2 = e2**2 / (alpha**2 * (n + 1.0) * (n + 2.0))
        
        # Parties purement géométriques/polynomiales
        P_N  = x1 + dx/2.0
        P_My = x1*y1 + (x1*dy + y1*dx)/2.0 + dx*dy/3.0
        P_Mz = (x1**2 + x1*dx + dx**2/3.0)/2.0
        
        # Assemblage avec les primitives exactes de la courbure
        N_  = fcd * dy * (P_N + C1 * I0_n1)
        My_ = fcd * dy * (P_My + C1 * (y1 * I0_n1 + dy * I1_n1))
        Mz_ = fcd * dy * (P_Mz + C1 * (x1 * I0_n1 + dx * I1_n1) + C2 * I0_n2)
        
    else:
        # --- Gradient dominant en Y --- 
        # Vecteur de Green : L_N = Int(sigma dy) => Contour sur -dx
        D1 = e2 / (beta * (n + 1.0))
        D2 = e2**2 / (beta**2 * (n + 1.0) * (n + 2.0))
        
        # Les composantes sont tournées à 90° (symétrie)
        P_N  = y1 + dy/2.0
        P_My = (y1**2 + y1*dy + dy**2/3.0)/2.0
        P_Mz = x1*y1 + (x1*dy + y1*dx)/2.0 + dx*dy/3.0
        
        # Assemblage final
        N_  = -fcd * dx * (P_N + D1 * I0_n1)
        My_ = -fcd * dx * (P_My + D1 * (y1 * I0_n1 + dy * I1_n1) + D2 * I0_n2)
        Mz_ = -fcd * dx * (P_Mz + D1 * (x1 * I0_n1 + dx * I1_n1))
        
    return np.array([N_, My_, Mz_])




def domaines_comprimes(contour, eps0, alpha, beta):

    def eps(x,y):
        return eps0 + alpha*x + beta*y

    pts = np.asarray(contour)
    n = len(pts)

    domaines = []
    poly = []

    for i in range(n):

        xa, ya = pts[i]
        xb, yb = pts[(i+1) % n]

        ea = eps(xa, ya)
        eb = eps(xb, yb)

        # entrée en compression
        if ea <= 0 and eb > 0:
            t = -ea / (eb - ea)
            xi = xa + t*(xb-xa)
            yi = ya + t*(yb-ya)
            poly = [(xi, yi), (xb, yb)]

        # sortie compression
        elif ea > 0 and eb <= 0:
            t = -ea / (eb - ea)
            xi = xa + t*(xb-xa)
            yi = ya + t*(yb-ya)
            poly.append((xi, yi))
            domaines.append(poly)
            poly = []

        # entièrement en compression
        elif ea > 0 and eb > 0:
            if not poly:
                poly = [(xa, ya)]
            poly.append((xb, yb))

        # cas traction → rien

    # fermer si nécessaire
    if poly:
        domaines.append(poly)

    return domaines

from shapely.geometry import Polygon, LineString, MultiPolygon, GeometryCollection
from shapely.ops import split

def _orienter_polygone_ccw(pts):
    """
    Force un polygone orienté anti-horaire (CCW)
    """

    pts = np.asarray(pts)

    A = 0.0
    n = len(pts)

    for i in range(n):
        x1, y1 = pts[i]
        x2, y2 = pts[(i+1) % n]
        A += x1*y2 - x2*y1

    if A < 0:
        return pts[::-1]   # inversion
    else:
        return pts

def _zones_comprimees(vertices, eps0, alpha, beta):

    poly = Polygon(vertices)

    if not poly.is_valid or poly.area == 0:
        return []

    # droite ε = 0
    L = 1e4

    if abs(beta) > 1e-12:
        x_vals = [-L, L]
        y_vals = [(-eps0 - alpha*x) / beta for x in x_vals]
    else:
        if abs(alpha) < 1e-12:
            return [poly] if eps0 > 0 else []
        x = -eps0 / alpha
        x_vals = [x, x]
        y_vals = [-L, L]

    line = LineString(list(zip(x_vals, y_vals)))

    try:
        result = split(poly, line)
    except:
        return []

    # ---------------------------------
    # ✅ gestion robuste des types
    # ---------------------------------
    zones = []

    def traiter_geom(g):
        if g.is_empty:
            return

        if g.geom_type == 'Polygon':
            cx, cy = g.centroid.x, g.centroid.y
            if eps0 + alpha*cx + beta*cy > 0:
                zones.append(g)

        elif g.geom_type in ['MultiPolygon', 'GeometryCollection']:
            for sub in g.geoms:
                traiter_geom(sub)

    traiter_geom(result)

    return zones

# =============================================================================
# DÉCOUPE GÉOMÉTRIQUE PAR PLANS DE DÉFORMATION
# =============================================================================

def _couper_poly_par_plan(poly, eps0, alpha, beta, valeur_seuil):
    """Coupe un polygone Shapely selon une ligne d'iso-déformation."""
    if not poly.is_valid or poly.area == 0:
        return []
    
    L = 1e4
    if abs(beta) > 1e-12:
        x_vals = [-L, L]
        y_vals = [(-eps0 - alpha*x + valeur_seuil) / beta for x in x_vals]
    else:
        if abs(alpha) < 1e-12:
            return [poly] if (eps0 > valeur_seuil) else []
        x = (valeur_seuil - eps0) / alpha
        x_vals, y_vals = [x, x], [-L, L]

    line = LineString(list(zip(x_vals, y_vals)))
    try:
        result = split(poly, line)
    except:
        return [poly]

    zones = []
    if result.geom_type == 'Polygon':
        zones.append(result)
    else:
        for g in result.geoms:
            if g.geom_type == 'Polygon' and not g.is_empty:
                zones.append(g)
    return zones

# =============================================================================
# MOTEUR D'INTÉGRATION PRINCIPAL
# =============================================================================
def clip_polygon_eps(vertices, eps0, alpha, beta, seuil, keep_above=True, tol=1e-12):
    """
    Coupe un polygone par un demi-plan défini par :
        eps(x,y) >= seuil   (si keep_above=True)
        eps(x,y) <= seuil   (si keep_above=False)

    Paramètres
    ----------
    vertices : array-like (N,2)
        Sommets du polygone (ordre quelconque mais fermé implicitement)
    eps0, alpha, beta : float
        Champ affine eps(x,y) = eps0 + alpha*x + beta*y
    seuil : float
        Valeur seuil
    keep_above : bool
        True  -> garde eps >= seuil
        False -> garde eps <= seuil
    tol : float
        Tolérance numérique

    Retour
    ------
    np.ndarray (M,2)
        Polygone clipé
    """

    def eps(x, y):
        return eps0 + alpha*x + beta*y

    def is_inside(e):
        if keep_above:
            return e >= seuil - tol
        else:
            return e <= seuil + tol

    output = []
    vertices = np.asarray(vertices, dtype=float)
    n = len(vertices)

    for i in range(n):
        x1, y1 = vertices[i]
        x2, y2 = vertices[(i + 1) % n]

        e1 = eps(x1, y1)
        e2 = eps(x2, y2)

        inside1 = is_inside(e1)
        inside2 = is_inside(e2)

        # ✅ cas 1 : les deux points dedans
        if inside1 and inside2:
            output.append((x2, y2))

        # ✅ cas 2 : sortie du domaine
        elif inside1 and not inside2:
            if abs(e2 - e1) > tol:
                t = (seuil - e1) / (e2 - e1)
                xi = x1 + t * (x2 - x1)
                yi = y1 + t * (y2 - y1)
                output.append((xi, yi))

        # ✅ cas 3 : entrée dans le domaine
        elif not inside1 and inside2:
            if abs(e2 - e1) > tol:
                t = (seuil - e1) / (e2 - e1)
                xi = x1 + t * (x2 - x1)
                yi = y1 + t * (y2 - y1)
                output.append((xi, yi))
            output.append((x2, y2))

        # ✅ cas 4 : dehors → rien

    if len(output) == 0:
        return np.zeros((0, 2))

    return np.array(output)


def _integrer_polygone(contour, eps0, alpha, beta, mode, C=None, fck=None, fcd=None, e2=None):

    pts = _orienter_polygone_ccw(np.asarray(contour, float))

    # ✅ zone comprimée
    poly_c = clip_polygon_eps(pts, eps0, alpha, beta, 0.0)
    
    if len(poly_c) < 3:
        return np.zeros(4)

    res = np.zeros(4)

    if mode == 'ELS':

        edges = np.column_stack([poly_c, np.roll(poly_c, -1, axis=0)])

        for xa, ya, xb, yb in edges:
            dx = xb - xa
            dy = yb - ya

            res[3] += dy * (xa + dx/2)
            res[:3] += _contrib_ELS(xa, ya, xb, yb, eps0, alpha, beta, C)

        return res















 


    else:
        n= eps_n(fck)
        # ✅ zone rectangle
        poly_r = clip_polygon_eps(poly_c, eps0, alpha, beta, e2, keep_above=True)

        # ✅ zone parabole = poly_c - poly_r
        res_r = np.zeros(4)
        if len(poly_r) >= 3:
            edges = np.column_stack([poly_r, np.roll(poly_r, -1, axis=0)])
            for xa, ya, xb, yb in edges:
                dx = xb - xa
                dy = yb - ya
                res_r[3] += dy * (xa + dx/2)
                res_r[:3] += _contrib_ELU_R(xa, ya, xb, yb, fcd)

        # ✅ zone parabole (reste)
        poly_p = clip_polygon_eps(poly_c, eps0, alpha, beta, e2, keep_above=False)
        res_p = np.zeros(4)
        edges = np.column_stack([poly_p, np.roll(poly_p, -1, axis=0)])

        for xa, ya, xb, yb in edges:

            xm = 0.5*(xa+xb)
            ym = 0.5*(ya+yb)
            if eps0 + alpha*xm + beta*ym <= e2:

                dx = xb - xa
                dy = yb - ya
                res_p[3] += dy * (xa + dx/2)
                res_p[:3] += _contrib_ELU_P(xa, ya, xb, yb,
                                            eps0, alpha, beta, fcd, e2, n=n)


        return res_p + res_r




# ═══════════════════════════════════════════════════════════════
# BLOC 3 — FONCTIONS PUBLIQUES
# ═══════════════════════════════════════════════════════════════

@func
def S_com(polygon1, evi, eps0, alpha, beta):

    contour_cg, evidements_cg, Cx, Cy = transformation_repere_cg(polygon1, evi)

    res_c = _integrer_polygone(
        contour_cg, eps0, alpha, beta,
        mode='ELS', C=1.0
    )

    res_v = np.zeros(4)

    for trou in evidements_cg:
        res_v += _integrer_polygone(
            trou, eps0, alpha, beta,
            mode='ELS', C=1.0
        )

    return res_c[3] - res_v[3]


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
    """ELS — N, My, Mz par Green analytique EXACT. Sans nettoyer_polygone."""
    contour_cg, evidements_cg, Cx, Cy = transformation_repere_cg(polygon1, evi)
    if contour_cg is None:
        return 0.0, 0.0, 0.0

    C = 200000.0 / (float(n) * 1000.0)

    res_c = _integrer_polygone(contour_cg, eps0, alpha, beta, 'ELS', C=C)
    res_v = np.zeros(4)
    for trou in evidements_cg:
        if len(trou) >= 3:
            res_v += _integrer_polygone(trou, eps0, alpha, beta, 'ELS', C=C)

    Nc, Mc_y, Mc_z = (res_c - res_v)[:3]


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
   
    if contour_cg is None:
        return 0.0, 0.0, 0.0

    C = 200000.0 / (float(n) * 1000.0)

    res_c = _integrer_polygone(contour_cg, eps0, alpha, beta, 'ELS', C=C)
    res_v = np.zeros(4)
    for trou in evidements_cg:
        if len(trou) >= 3:
            res_v += _integrer_polygone(trou, eps0, alpha, beta, 'ELS', C=C)

    Nc, Mc_y, Mc_z = (res_c - res_v)[:3]


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
    
    #n = float(eps_n(fck))
    e2 = float(eps_c2(fck))

    res_c = _integrer_polygone(contour_cg, eps0, alpha, beta,
                               'ELU',fck=fck, fcd=fcd, e2=e2)
    res_v = np.zeros(4)
    for trou in evidements_cg:
        if len(trou) >= 3:
            res_v += _integrer_polygone(trou, eps0, alpha, beta,
                                        'ELU',fck=fck, fcd=fcd, e2=e2)

    Nc, Mc_y, Mc_z = (res_c - res_v)[:3]
         

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
