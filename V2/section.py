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


	
@func
def sigma_s_palier(fyd, k, eps_uk,eps_ud, eps_s_array):
    """renvoie la contrainte acier (fe500 classe B) selon 3.2.7 de NF EN 1992-1-1"""
    eps_s =np.array(eps_s_array)
    eps_sd = fyd / 200000 * 1000 

    
    sigma1 = 200000 * eps_s / 1000 # if abs(eps_s) <= eps_sd:
   
    sigma2 = sgn(eps_s)* (fyd + (k * fyd - fyd)/(eps_uk - eps_sd) * (np.abs(eps_s) - eps_sd)) #  else:   
    return  np.where(np.abs(eps_s) > eps_sd, sigma2, sigma1)	
	
@func
def sigma_s_palier1(fyd, k, eps_uk,eps_ud, eps_s_array,a_com):
    """renvoie la contrainte acier (fe500 classe B) selon 3.2.7 de NF EN 1992-1-1"""
    eps_s =np.array(eps_s_array)

    
    sigma1 =sigma_s_palier(fyd, k, eps_uk,eps_ud, eps_s)*a_com # if eps_s > 0.0:
    sigma2 = sigma_s_palier(fyd, k, eps_uk,eps_ud, eps_s)   # else:
    return  np.where(eps_s > 0, sigma1, sigma2)
@func
def eps_s_palier(fyd, k, eps_uk, eps_ud, sig_s_array):
    """
    Version vectorisée analytique : calcule la déformation eps_s à partir de sig_s.
    sig_s_array : tableau des contraintes (MPa)
    """
    sig = np.asarray(sig_s_array)
    Es = 200000
    eps_sd = (fyd / Es) * 1000
    
    # Préparation des paramètres de la pente plastique
    # Attention : si k=1 (palier horizontal), la déformation n'est pas unique.
    # On gère le cas k > 1
    pente_inv = (eps_uk - eps_sd) / (k * fyd - fyd) if k > 1 else 0
    
    abs_sig = np.abs(sig)
    
    # 1. Calcul pour la zone élastique
    eps_elaste = (sig / Es) * 1000
    
    # 2. Calcul pour la zone plastique (branche inclinée)
    eps_plast = np.sign(sig) * (eps_sd + (abs_sig - fyd) * pente_inv)
    
    # 3. Assemblage des conditions
    # Note : Si sig_s == 0, on renvoie 0.0 au lieu de None pour la compatibilité
    conditions = [
        (sig == 0),                          # Cas sigma nul
        (abs_sig <= fyd),                    # Zone élastique
        (abs_sig > fyd) & (abs_sig <= k*fyd) # Zone plastique
    ]
    
    choix = [
        0.0,
        eps_elaste,
        eps_plast
    ]
    
    # On renvoie 0.0 par défaut si sigma est hors limites
    return np.select(conditions, choix, default=0.0)
	
@func
def sigma_s_lin(eps_s_array):
    """
    Version vectorisée de la loi acier linéaire infinie.
    E = 200 000 MPa (pas de palier plastique).
    eps_s_array : déformation en mm/m.
    """
    # On s'assure que l'entrée est un tableau
    eps = np.asarray(eps_s_array)
    
    # sigma = E * epsilon
    # Le /1000 convertit les mm/m en déformation relative (sans unité)
    return 200000 * eps / 1000
	
@func
def eps_s_lin(sig_s_array):
    """
    Version vectorisée : renvoie la déformation eps_s (en mm/m).
    Loi élastique linéaire infinie (E = 200 000 MPa).
    sig_s_array : contrainte acier (MPa).
    """
    # Conversion en array NumPy
    sig = np.asarray(sig_s_array)
    
    # epsilon = (sigma / E) * 1000 pour avoir des mm/m
    return 1000 * sig / 200000	
	
@func
def sigma_s_lin1(eps_s, a_com):
    """
    Loi acier linéaire avec coefficient de pondération a_com.
    Appliqué si eps_s > 0 (traction).
    """
    eps = np.asarray(eps_s)
    # Calcul de la contrainte de base (E = 200 000 MPa)
    sig_base = 200000 * eps / 1000
    
    # On multiplie par a_com uniquement là où eps_s > 0
    return np.where(eps > 0.0, sig_base * a_com, sig_base)


@func
def sigma_c_n1(eps_c, n):
    """
    Loi béton fissuré. 
    Utilise une fonction de transfert pour annuler la traction sans np.where.
    """
    eps = np.asarray(eps_c)
    eps_abs = np.abs(eps)
    
    # Le terme (eps + eps_abs) / (2 * eps_abs + 1e-15) agit comme un filtre :
    # Si eps > 0  => (eps + eps) / (2 * eps) = 1
    # Si eps < 0  => (eps - eps) / (2 * eps) = 0
    # Le 1e-15 est le "garde-fou" contre la division par zéro.
    
    filtre = (eps + eps_abs) / (2 * eps_abs + 1e-30)
    sigma = (200000 / n * eps / 1000) * filtre
    return sigma



@func
def sigma_c_pararect1(fck, fcd, eps_c_array):
    eps_c=np.array(eps_c_array)
    eps_c2_val = eps_c2(fck)
    eps_n_val = eps_n(fck)
    term1 = np.abs(1 - eps_c / eps_c2_val)
    term2 = (1 - eps_c / eps_c2_val)
    sigma = fcd * (1 - (term1 + term2) / (2 * term1 + 1e-30) * term1 ** eps_n_val) * (eps_c + np.abs(eps_c)) / (2 * np.abs(eps_c) + 1e-30)
    return sigma


# -------------------------
# Intégration adaptative
# -------------------------
@func
def polygone_integrate(f, vertices, tol=1e-9, rtol=1e-9, max_depth=3):
    """
    Intégrateur universel (Scalaire ou Vectoriel).
    Renvoie 0.0 ou un vecteur de 0.0 au lieu de None en cas d'entrée vide.
    """
    # -- Constantes de Gauss (cache) -----------------------------------------
    if not hasattr(polygone_integrate, '_gauss'):
        pts = np.array([
            [0.3333333333333333, 0.3333333333333333],
            [0.1012865073234099, 0.1012865073234099],
            [0.7974269853530873, 0.1012865073234099],
            [0.1012865073234099, 0.7974269853530873],
            [0.4701420641051151, 0.0597158717897698],
            [0.4701420641051151, 0.4701420641051151],
            [0.0597158717897698, 0.4701420641051151],
        ], dtype=np.float64)

        w = np.array([0.225, 0.1259391805448272, 0.1259391805448272, 0.1259391805448272, 
                      0.1323941527885062, 0.1323941527885062, 0.1323941527885062], dtype=np.float64)
        polygone_integrate._gauss = (pts.T, w)

    gpts_T, gw = polygone_integrate._gauss

    # -- Quadrature de Gauss adaptée -----------------------------------------
    def gauss_integral(v0, v1, v2, area):
        J0x, J0y = v1[0] - v0[0], v1[1] - v0[1]
        J1x, J1y = v2[0] - v0[0], v2[1] - v0[1]

        X = v0[0] + J0x * gpts_T[0] + J1x * gpts_T[1]
        Y = v0[1] + J0y * gpts_T[0] + J1y * gpts_T[1]
        
        valeurs = np.asarray(f(X, Y)) 
        
        if valeurs.ndim > 1:
            return area * np.dot(gw, valeurs.T)
        return area * np.dot(gw, valeurs)

    # -- Version itérative ---------------------------------------------------
    def integrate_triangle_iter(v0, v1, v2, atol):
        stack = [(v0, v1, v2, atol, 0)]
        total_tri = None 

        while stack:
            v0, v1, v2, atol, depth = stack.pop()
            
            area = 0.5 * abs((v1[0]-v0[0])*(v2[1]-v0[1]) - (v2[0]-v0[0])*(v1[1]-v0[1]))
            if area < 1e-15: continue

            coarse = gauss_integral(v0, v1, v2, area)
            
            if total_tri is None:
                total_tri = np.zeros_like(coarse)

            if depth >= max_depth:
                total_tri += coarse
                continue

            m01, m12, m02 = 0.5*(v0+v1), 0.5*(v1+v2), 0.5*(v0+v2)
            a4 = area * 0.25
            
            g0 = gauss_integral(v0, m01, m02, a4); g1 = gauss_integral(v1, m01, m12, a4)
            g2 = gauss_integral(v2, m12, m02, a4); g3 = gauss_integral(m01, m12, m02, a4)
            
            fine = g0 + g1 + g2 + g3

            err = np.max(np.abs(fine - coarse))
            scale = np.max(np.abs(fine))
            if scale < 1e-30: scale = 1e-30

            if err <= atol and err <= rtol * scale:
                total_tri += fine
            else:
                atol4 = atol * 0.25
                stack.extend([(v0, m01, m02, atol4, depth+1), (v1, m01, m12, atol4, depth+1),
                              (v2, m12, m02, atol4, depth+1), (m01, m12, m02, atol4, depth+1)])
        
        # Sécurité : si le triangle était dégénéré, on renvoie une valeur neutre
        return total_tri if total_tri is not None else 0.0

    # -- Triangulation et boucle principale ----------------------------------
    if vertices is None or len(vertices) < 3:
        return 0.0

    vertices = np.asarray(vertices, dtype=np.float64)
    tri = Delaunay(vertices)
    polygon_path = Path(vertices)
    
    total_final = None

    for simplex in tri.simplices:
        v0, v1, v2 = vertices[simplex]
        if polygon_path.contains_point((v0 + v1 + v2) / 3.0):
            res = integrate_triangle_iter(v0, v1, v2, tol)
            
            if total_final is None:
                total_final = np.zeros_like(res)
            total_final += res

    # Renvoi final sécurisé
    return total_final if total_final is not None else 0.0
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
            s = float(s_acier[i])   # peut être 0 ou négatif ? CONSERVÉ ?

            newacier.append([x, y, s / 10000.0])  # mm² ? cm²

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

   

@func
def calculer_N_My_Mz(
    polygon1, evi,
    p_acier, s_acier,
    a_com, n,
    eps0, alpha, beta
):
    """
    Calcul de N, My, Mz dans le repère du centre de gravité,
    avec gestion automatique du contour et des évidements
    (0 ou plusieurs).
    """

    # --------------------------------------------------
    # 0. Transformation géométrique vers le repère CG
    # --------------------------------------------------
    # transformation_repere_cg calcule Cx, Cy AUTOMATIQUEMENT
    contour_cg, evidements_cg, Cx, Cy = transformation_repere_cg(polygon1, evi)

    # --------------------------------------------------
    # 1. FONCTION D’INTÉGRATION VECTORIELLE (BÉTON)
    # --------------------------------------------------
    def f_complet(x, y):
        eps = eps0 + alpha * x + beta * y
        sig = sigma_c_n1(eps, n)
        return np.array([sig, sig * y, sig * x])

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

    # --------------------------------------------------
    # 2. PARTIE ACIER (vecteur)
    # --------------------------------------------------

    if p_acier is None or len(p_acier) == 0:
        Ns, Msy, Msz = 0.0, 0.0, 0.0
    else:
        acier_data = np.array(acier_G(p_acier, s_acier), dtype=float)

    if acier_data.size == 0:
            Ns, Msy, Msz = 0.0, 0.0, 0.0

    else:
        acier = np.array(acier_G(p_acier, s_acier), dtype=float)

        # ?? translation des aciers vers le repère CG
        x_s = acier[:, 0] - Cx
        y_s = acier[:, 1] - Cy
        areas_s = acier[:, 2]

        eps_s = eps0 + alpha * x_s + beta * y_s
        sig_s = sigma_s_lin1(eps_s, a_com)

        forces_s = sig_s * areas_s
        Ns = np.sum(forces_s)
        Msy = np.sum(forces_s * y_s)
        Msz = np.sum(forces_s * x_s)

    # --------------------------------------------------
    # 3. TOTAL
    # --------------------------------------------------
    return Ns + Nc, Msy + Mc_y, Msz + Mc_z

@func
def solve_GG_ELS(polygon1, evi, p_acier, s_acier, a_com, n, Nobj, Myobj, Mzobj):
    x0 = np.array([0., 0.001, 0.0])  # Valeurs initiales [eps0, alpha, beta]
    def residuals(eps):
        N, My, Mz = calculer_N_My_Mz(polygon1, evi, p_acier, s_acier, a_com, n, *eps)
        return np.array([N - Nobj, My - Myobj, Mz - Mzobj])
    #result = fsolve(residuals, x0)
    #return result
    def jacobian(eps, h=1e-6):
        """Jacobien par différences finies centrées (3×3)."""
        J = np.empty((3, 3))
        r0 = residuals(eps)
        for i in range(3):
            dh = np.zeros(3); dh[i] = h
            J[:, i] = (residuals(eps + dh) - r0) / h
        return J

    result = root(
        residuals,
        x0,
        jac=jacobian,
        method= "hybr", #"hybr",  # ou "lm"
        tol=1e-5,
        options={"maxfev": 800}
        )

    return result.x

#-------------------------------------------------------------------------------
@func
def resultats_GG_ELS(
    polygon1, evi,
    p_acier, s_acier,
    a_com, n,
    eps0, alpha, beta
):
    """
    Résultats GG ELS calculés dans le repère du centre de gravité.
    Cohérent avec calculer_N_My_Mz et S_com.
    """

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

    # --------------------------------------------------
    # 3. Calculs globaux (cohérents CG)
    # --------------------------------------------------
    A_com = S_com(polygon1, evi, eps0, alpha, beta)

    N, My, Mz = calculer_N_My_Mz(
        polygon1, evi,
        p_acier, s_acier,
        a_com, n,
        eps0, alpha, beta
    )

    # --------------------------------------------------
    # 4. Résultats
    # --------------------------------------------------
    return {
        "ACOM": A_com,

        "EPS_C_MAX": espcmax,
        "EPS_C_MIN": espcmin,
        "SIG_C_MAX": sig_c_max,
        "SIG_C_MIN": sig_c_min,

        "EPS_S_MAX": espsmax,
        "EPS_S_MIN": espsmin,
        "SIG_S_MAX": sig_s_max,
        "SIG_S_MIN": sig_s_min,

        "N": N,
        "MY": My,
        "MZ": Mz,

        "CX": Cx,
        "CY": Cy
    }
@func
def e_resultats_GG_ELS(polygon1,evi,p_acier,s_acier,a_com, n, eps0,alpha, beta,resultats):
    resultats_Excel = []
    resultats_tout = resultats_GG_ELS(polygon1,evi,p_acier,s_acier,a_com, n, eps0,alpha, beta)                      
    resultats_list = resultats.split(',')
    for r in resultats_list:       
        resultats_Excel.append(resultats_tout[r])
    return resultats_Excel   # résultats en ligne

#--------------(eps0,béta)--->(N,M)


@func
def calculer_N_My_Mz_ELU_pararect(
    polygon1, evi,
    p_acier, s_acier,
    a_com, fck, fcd, fyd, k,
    eps_uk, eps_ud,
    eps0, alpha, beta
):
    """
    Calcul ELU N, My, Mz avec loi béton parabole-rectangle,
    calculé dans le repère du centre de gravité,
    avec gestion de 0 ou plusieurs évidements.
    """

    # ==================================================
    # 0. Transformation géométrique vers le repère CG
    # ==================================================
    contour_cg, evidements_cg, Cx, Cy = transformation_repere_cg(polygon1, evi)

    if contour_cg is None:
        return 0.0, 0.0, 0.0

    # ==================================================
    # 1. PARTIE BÉTON (plein + évidements)
    # ==================================================
    def f_ELU_complet(x, y):
        eps = eps0 + alpha * x + beta * y
        sig = sigma_c_pararect1(fck, fcd, eps)
        return np.array([sig, sig * y, sig * x])

    # --- Béton plein
    res_c = polygone_integrate(f_ELU_complet, contour_cg)

    # --- Béton évidements (ABS obligatoire)
    res_v = np.zeros(3)
    for trou in evidements_cg:
        if len(trou) >= 3:
            res_v += abs(polygone_integrate(f_ELU_complet, trou))

    # Résultantes béton
    Nc, Mc_y, Mc_z = res_c - res_v

    # ==================================================
    # 2. PARTIE ACIER (repère CG)
    # ==================================================
    if p_acier is None or len(p_acier) == 0:
        Ns, Msy, Msz = 0.0, 0.0, 0.0
    else:
        acier_data = np.array(acier_G(p_acier, s_acier), dtype=float)

        if acier_data.size == 0:
            Ns, Msy, Msz = 0.0, 0.0, 0.0
        else:
            # Coordonnées recentrées
            x_s = acier_data[:, 0] - Cx
            y_s = acier_data[:, 1] - Cy
            areas_s = acier_data[:, 2]

            # Déformations
            eps_s = eps0 + alpha * x_s + beta * y_s

            # Contraintes acier (loi palier ELU)
            sig_s = sigma_s_palier1(
                fyd, k, eps_uk, eps_ud,
                eps_s, a_com
            )

            forces = sig_s * areas_s
            Ns = np.sum(forces)
            Msy = np.sum(forces * y_s)
            Msz = np.sum(forces * x_s)

    # ==================================================
    # 3. RÉSULTATS TOTAUX
    # ==================================================
    return Ns + Nc, Msy + Mc_y, Msz + Mc_z

@func
def solve_GG_ELU_pararect(polygon1, evi, p_acier, s_acier, a_com, fck, fcd, fyd, k, eps_uk, eps_ud, Nobj, Myobj, Mzobj):
    x0 = np.array([0.0, 0.1, 0.0])  # [eps0, alpha, beta] valeurs initiales

    def residuals(eps):
        N, My, Mz = calculer_N_My_Mz_ELU_pararect(
            polygon1, evi, p_acier, s_acier, a_com, fck, fcd, fyd, k, eps_uk, eps_ud, *eps
        )
        return np.array([N - Nobj, My - Myobj, Mz - Mzobj])

    #result = root(residuals, x0, method='hybr')
    #return result.x

    
    #result = fsolve(residuals, x0)
    def jacobian(eps, h=1e-6):
        """Jacobien par différences finies centrées (3×3)."""
        J = np.empty((3, 3))
        r0 = residuals(eps)
        for i in range(3):
            dh = np.zeros(3); dh[i] = h
            J[:, i] = (residuals(eps + dh) - r0) / h
        return J

    result = root(
        residuals,
        x0,
        jac=jacobian,
        method= "hybr", #"hybr",  # ou "lm"
        tol=1e-5,
        options={"maxfev": 800}
    )

    return result.x

@func    
def resultats_GG_ELU_pararect(
    polygon1, evi,
    p_acier, s_acier,
    a_com, fck, fcd, fyd, k,
    eps_uk, eps_ud,
    eps0, alpha, beta
):
    """
    Résultats GG ELU avec loi béton parabole-rectangle,
    calculés dans le repère du centre de gravité,
    cohérents avec calculer_N_My_Mz_ELU_pararect et S_com.
    """

    # ==================================================
    # 0. Transformation géométrique vers le repère CG
    # ==================================================
    contour_cg, evidements_cg, Cx, Cy = transformation_repere_cg(polygon1, evi)

    if contour_cg is None:
        return {}

    # ==================================================
    # 1. BÉTON – déformations et contraintes (repère CG)
    # ==================================================
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

    # ==================================================
    # 2. ACIER – déformations, contraintes, plastification
    # ==================================================
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

            # coordonnées recentrées
            x_s = acier_array[:, 0] - Cx
            y_s = acier_array[:, 1] - Cy

            eps_s = eps0 + alpha * x_s + beta * y_s
            epssmax = float(np.max(eps_s))
            epssmin = float(np.min(eps_s))

            # Pourcentage plastifié (Es ˜ 200 000 MPa)
            plast = np.sum(np.abs(eps_s) >= (float(fyd) / 200.0))
            pour_a_pl = float(plast / len(eps_s) * 100.0)

            sig_s_max = float(sigma_s_palier1(
                fyd, k, eps_uk, eps_ud, epssmax, a_com
            ))
            sig_s_min = float(sigma_s_palier1(
                fyd, k, eps_uk, eps_ud, epssmin, a_com
            ))

    # ==================================================
    # 3. SURFACE COMPRIMÉE (cohérente CG)
    # ==================================================
    A_com = float(S_com(polygon1, evi, eps0, alpha, beta))

    # ==================================================
    # 4. EFFORTS GLOBAUX ELU (N, My, Mz)
    # ==================================================
    N, My, Mz = calculer_N_My_Mz_ELU_pararect(
        polygon1, evi,
        p_acier, s_acier,
        a_com, fck, fcd, fyd, k,
        eps_uk, eps_ud,
        eps0, alpha, beta
    )

    # ==================================================
    # 5. RÉSULTATS
    # ==================================================
    return {
        "EPS0": eps0,
        "ALPHA": alpha,
        "BETA": beta,

        "ACOM": A_com,

        "EPS_C_MAX": espcmax,
        "EPS_C_MIN": espcmin,
        "SIG_C_MAX": sig_c_max,
        "SIG_C_MIN": sig_c_min,

        "EPS_S_MAX": epssmax,
        "EPS_S_MIN": epssmin,
        "SIG_S_MAX": sig_s_max,
        "SIG_S_MIN": sig_s_min,

        "N": float(N),
        "MY": float(My),
        "MZ": float(Mz),

        "PA": pour_a_pl,

        "CX": Cx,
        "CY": Cy
    }


@func
def e_resultats_GG_ELU_pararect(polygon1,evi,p_acier,s_acier,a_com,fck,fcd,fyd,k,eps_uk,eps_ud, eps0,alpha, beta,resultats):
    resultats_Excel = []
    resultats_tout = resultats_GG_ELU_pararect(polygon1,evi,p_acier,s_acier,a_com,fck,fcd,fyd,k,eps_uk,eps_ud, eps0,alpha, beta)                      
    resultats_list = resultats.split(',')
    for r in resultats_list:       
        resultats_Excel.append(resultats_tout[r])
    return resultats_Excel   # résultats en ligne


 

