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




# section.py — VERSION FINALE COMPLETE
# Green analytique exact — ELS + ELU EC2 (n=2)

import numpy as np
from scipy.optimize import root
from xlwings import func

from beton import eps_c2, sigma_c_n1, sigma_c_pararect1
from acier import sigma_s_lin1, sigma_s_palier1

# ═══════════════════════════════════════════════════════════════
# BLOC 1 — UTILITAIRES
# ═══════════════════════════════════════════════════════════════

def _aire_signee(pts):
    a = 0.0
    for i in range(len(pts)):
        x1,y1 = pts[i]
        x2,y2 = pts[(i+1)%len(pts)]
        a += x1*y2 - x2*y1
    return 0.5*a


def _orienter_ccw(pts):
    if _aire_signee(pts) < 0:
        return pts[::-1]
    return pts


def points_valides(data):
    pts=[]
    for row in data:
        try:
            pts.append((float(row[0]), float(row[1])))
        except: pass
    return pts


def nettoyer_polygone(points):
    if len(points)<2:
        return points
    if np.linalg.norm(np.array(points[0])-np.array(points[-1]))<1e-9:
        return points[:-1]
    return points


def translation_points(data, Cx, Cy):
    return [(x-Cx,y-Cy) for x,y in data]


def transformation_repere_cg(contour,evi):

    pts=points_valides(contour)
    if len(pts)<3:
        return None,[],None,None

    A=Cx=Cy=0
    for i in range(len(pts)):
        x1,y1=pts[i]; x2,y2=pts[(i+1)%len(pts)]
        cross=x1*y2-x2*y1
        A+=cross; Cx+=(x1+x2)*cross; Cy+=(y1+y2)*cross

    A*=0.5
    if A==0:
        return None,[],None,None

    Cx/=(6*A); Cy/=(6*A)

    return translation_points(pts,Cx,Cy),[],Cx,Cy


def acier_G(p_acier,s_acier):
    if p_acier is None:
        return []
    res=[]
    for i in range(min(len(p_acier),len(s_acier))):
        try:
            res.append([float(p_acier[i][0]),
                        float(p_acier[i][1]),
                        float(s_acier[i])/10000])
        except: pass
    return res


# ═══════════════════════════════════════════════════════════════
# BLOC 2 — GREEN EXACT CORRIGÉ
# ═══════════════════════════════════════════════════════════════

def _decouper_arete(xa,ya,xb,yb,eps0,a,b,seuils):

    ea=eps0+a*xa+b*ya
    eb=eps0+a*xb+b*yb
    de=eb-ea

    cuts=[0.0]

    for s in seuils:
        if abs(de)>1e-15:
            t=(s-ea)/de
            if 1e-12<t<1-1e-12:
                cuts.append(t)

    cuts.append(1.0)
    cuts=sorted(set(cuts))

    segs=[]
    for i in range(len(cuts)-1):
        t1,t2=cuts[i],cuts[i+1]

        x1=xa+t1*(xb-xa)
        y1=ya+t1*(yb-ya)
        x2=xa+t2*(xb-xa)
        y2=ya+t2*(yb-ya)

        ea1=eps0+a*x1+b*y1
        eb1=eps0+a*x2+b*y2

        segs.append((x1,y1,x2,y2,ea1,eb1))

    return segs


def _contrib_ELS(x1,y1,x2,y2,eps0,a,b,C):

    dx=x2-x1; dy=y2-y1
    if abs(dy)<1e-15:
        return np.zeros(4)

    ea=eps0+a*x1+b*y1
    eb=eps0+a*x2+b*y2
    de=eb-ea

    N=C*dy*(ea+de/2)
    My=C*dy*(ea*y1+(ea*dy+y1*de)/2+de*dy/3)
    Mz=C*dy*(ea*x1+(ea*dx+x1*de)/2+de*dx/3)
    Sc=dy*(x1+dx/2)

    return np.array([N,My,Mz,Sc])


def _contrib_ELU_R(x1,y1,x2,y2,fcd):

    dx=x2-x1; dy=y2-y1
    if abs(dy)<1e-15:
        return np.zeros(4)

    N=fcd*dy
    My=fcd*dy*(y1+dy/2)
    Mz=fcd*dy*(x1+dx/2)
    Sc=dy*(x1+dx/2)

    return np.array([N,My,Mz,Sc])


def _contrib_ELU_P(x1,y1,x2,y2,eps0,a,b,fcd,e2):

    dx=x2-x1; dy=y2-y1
    if abs(dy)<1e-15:
        return np.zeros(4)

    ea=eps0+a*x1+b*y1
    eb=eps0+a*x2+b*y2

    ua=ea/e2; ub=eb/e2
    du=ub-ua

    I0=2*ua-ua**2+du-ua*du-du**2/3
    I1=(ua-ua**2/2)+(2*du/3-2*ua*du/3)-du**2/4

    N=fcd*dy*I0
    My=fcd*dy*(y1*I0+dy*I1)
    Mz=fcd*dy*(x1*I0+dx*I1)
    Sc=dy*(x1+dx/2)

    return np.array([N,My,Mz,Sc])


def _integrer_polygone(vertices,eps0,a,b,mode,C=None,fcd=None,e2=None):

    pts=_orienter_ccw(np.asarray(vertices,float))
    res=np.zeros(4)

    for i in range(len(pts)):
        xa,ya=pts[i]; xb,yb=pts[(i+1)%len(pts)]

        if mode=='ELS':
            segs=_decouper_arete(xa,ya,xb,yb,eps0,a,b,[0])

            for x1,y1,x2,y2,ea,eb in segs:
                if ea<=0 and eb<=0: continue
                res+=_contrib_ELS(x1,y1,x2,y2,eps0,a,b,C)

        else:
            segs=_decouper_arete(xa,ya,xb,yb,eps0,a,b,[0,e2])

            for x1,y1,x2,y2,ea,eb in segs:

                if ea<=0 and eb<=0:
                    continue
                elif ea>=e2 and eb>=e2:
                    res+=_contrib_ELU_R(x1,y1,x2,y2,fcd)
                else:
                    res+=_contrib_ELU_P(x1,y1,x2,y2,eps0,a,b,fcd,e2)

    return res


# ═══════════════════════════════════════════════════════════════
# BLOC 3 — ELS + ELU
# ═══════════════════════════════════════════════════════════════

def calculer_N_My_Mz(polygon1,evi,p_acier,s_acier,a_com,n,eps0,alpha,beta):

    contour,_,Cx,Cy=transformation_repere_cg(polygon1,evi)
    if contour is None:
        return 0,0,0

    C=200000/(n*1000)
    contour=nettoyer_polygone(contour)

    res=_integrer_polygone(contour,eps0,alpha,beta,'ELS',C=C)
    Nc,Myc,Mzc=res[:3]

    Ns=Mys=Mzs=0
    acier=acier_G(p_acier,s_acier)

    if len(acier)>0:
        ac=np.array(acier)
        x=ac[:,0]-Cx; y=ac[:,1]-Cy
        eps=eps0+alpha*x+beta*y
        sig=sigma_s_lin1(eps,a_com)
        F=sig*ac[:,2]
        Ns=np.sum(F); Mys=np.sum(F*y); Mzs=np.sum(F*x)

    return Ns+Nc, Mys+Myc, Mzs+Mzc


def calculer_N_My_Mz_ELU_pararect(polygon1,evi,p_acier,s_acier,
                                 a_com,fck,fcd,fyd,k,eps_uk,eps_ud,
                                 eps0,alpha,beta):

    contour,_,Cx,Cy=transformation_repere_cg(polygon1,evi)
    if contour is None:
        return 0,0,0

    contour=nettoyer_polygone(contour)
    e2=float(eps_c2(fck))

    res=_integrer_polygone(contour,eps0,alpha,beta,'ELU',fcd=fcd,e2=e2)
    Nc,Myc,Mzc=res[:3]

    Ns=Mys=Mzs=0
    acier=acier_G(p_acier,s_acier)

    if len(acier)>0:
        ac=np.array(acier)
        x=ac[:,0]-Cx; y=ac[:,1]-Cy
        eps=eps0+alpha*x+beta*y
        sig=sigma_s_palier1(fyd,k,eps_uk,eps_ud,eps,a_com)
        F=sig*ac[:,2]
        Ns=np.sum(F); Mys=np.sum(F*y); Mzs=np.sum(F*x)

    return Ns+Nc, Mys+Myc, Mzs+Mzc


def solve_GG_ELS(polygon1,evi,p_acier,s_acier,
                 a_com,n,Nobj,Myobj,Mzobj):

    x0=np.array([0.0,0.001,0.0])

    def R(eps):
        N,My,Mz=calculer_N_My_Mz(
            polygon1,evi,p_acier,s_acier,a_com,n,
            eps[0],eps[1],eps[2])
        return np.array([N-Nobj,My-Myobj,Mz-Mzobj])

    return root(R,x0).x


def solve_GG_ELU_pararect(polygon1,evi,p_acier,s_acier,
                         a_com,fck,fcd,fyd,k,eps_uk,eps_ud,
                         Nobj,Myobj,Mzobj):

    x0=np.array([0.0,0.05,0.0])

    def R(eps):
        N,My,Mz=calculer_N_My_Mz_ELU_pararect(
            polygon1,evi,p_acier,s_acier,
            a_com,fck,fcd,fyd,k,eps_uk,eps_ud,
            eps[0],eps[1],eps[2])
        return np.array([N-Nobj,My-Myobj,Mz-Mzobj])

    return root(R,x0).x

def resultats_GG_ELU_pararect(
    polygon1, evi,
    p_acier, s_acier,
    a_com, fck, fcd, fyd, k,
    eps_uk, eps_ud,
    eps0, alpha, beta
):
    """
    Résultats complets ELU — intégration analytique Green.
    """
    from section import (transformation_repere_cg, S_com,
                         acier_G, sigma_c_pararect1, sigma_s_palier1)

    contour_cg, evidements_cg, Cx, Cy = transformation_repere_cg(polygon1, evi)
    if contour_cg is None:
        return {}

    # ── Béton ────────────────────────────────────────────────────────
    poly_array = np.asarray(contour_cg, dtype=float)
    eps_c    = eps0 + alpha * poly_array[:, 0] + beta * poly_array[:, 1]
    espcmax  = float(np.max(eps_c))
    espcmin  = float(np.min(eps_c))
    sig_c_max = float(sigma_c_pararect1(fck, fcd, espcmax))
    sig_c_min = float(sigma_c_pararect1(fck, fcd, espcmin))

    # ── Acier ────────────────────────────────────────────────────────
    epssmax = epssmin = ""
    sig_s_max = sig_s_min = ""
    pour_a_pl = 0.0

    if p_acier is not None and len(p_acier) > 0:
        acier = acier_G(p_acier, s_acier)
        if len(acier) > 0:
            acier_array = np.asarray(acier, dtype=float)
            x_s  = acier_array[:, 0] - Cx
            y_s  = acier_array[:, 1] - Cy
            eps_s = eps0 + alpha * x_s + beta * y_s
            epssmax = float(np.max(eps_s))
            epssmin = float(np.min(eps_s))
            plast   = np.sum(np.abs(eps_s) >= (float(fyd) / 200.0))
            pour_a_pl = float(plast / len(eps_s) * 100.0)
            sig_s_max = float(sigma_s_palier1(fyd, k, eps_uk, eps_ud, epssmax, a_com))
            sig_s_min = float(sigma_s_palier1(fyd, k, eps_uk, eps_ud, epssmin, a_com))

    # ── Efforts ──────────────────────────────────────────────────────
    A_com = float(S_com(polygon1, evi, eps0, alpha, beta))
    N, My, Mz = calculer_N_My_Mz_ELU_pararect(
        polygon1, evi, p_acier, s_acier,
        a_com, fck, fcd, fyd, k, eps_uk, eps_ud,
        eps0, alpha, beta
    )

    return {
        "EPS0": eps0, "ALPHA": alpha, "BETA": beta,
        "ACOM": A_com,
        "EPS_C_MAX": espcmax,  "EPS_C_MIN": espcmin,
        "SIG_C_MAX": sig_c_max, "SIG_C_MIN": sig_c_min,
        "EPS_S_MAX": epssmax,  "EPS_S_MIN": epssmin,
        "SIG_S_MAX": sig_s_max, "SIG_S_MIN": sig_s_min,
        "N": float(N), "MY": float(My), "MZ": float(Mz),
        "PA": pour_a_pl, "CX": Cx, "CY": Cy
    }


def e_resultats_GG_ELU_pararect(
    polygon1, evi, p_acier, s_acier,
    a_com, fck, fcd, fyd, k, eps_uk, eps_ud,
    eps0, alpha, beta, resultats
):
    tout = resultats_GG_ELU_pararect(
        polygon1, evi, p_acier, s_acier,
        a_com, fck, fcd, fyd, k, eps_uk, eps_ud,
        eps0, alpha, beta
    )
    return [tout[r] for r in resultats.split(',')]

