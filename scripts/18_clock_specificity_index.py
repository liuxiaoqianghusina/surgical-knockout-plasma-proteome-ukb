"""器官时钟特异性指数：被切器官自身时钟的年龄差变化 − 其余 11 个器官时钟（不含 Organismal、Conventional）变化的平均值。
同一设计矩阵下多结局 OLS 的系数协方差 = 残差协方差 / (e_r'e_r)，据此得到对比的标准误。另给出调整 eGFR 后的结果。"""
import pandas as pd, numpy as np, sys
from scipy import stats
sys.path.insert(0,"."); from b_common import *
X=pd.read_pickle(f"{C}/atlas_cohort.pkl"); G=pd.read_pickle(f"{C}/organ_agegap_oh.pkl")
organs=[c for c in G.columns if c not in ("Organismal","Conventional")]
scr=X.creat/88.4; fem=(X.sex==0); k=np.where(fem,0.7,0.9); a=np.where(fem,-0.241,-0.302)
X["egfr"]=142*np.minimum(scr/k,1)**a*np.maximum(scr/k,1)**-1.2*0.9938**X.age*np.where(fem,1.012,1)
X["egfr_z"]=(X.egfr-X.egfr.mean())/X.egfr.std(); X["egfr_z2"]=X.egfr_z**2; X["egfr_z3"]=X.egfr_z**3
TGT={"tonsil":"Immune","appendix":"Intestine","spleen":"Immune","stomach":"Intestine","colon":"Intestine","colon_right":"Intestine","small_bowel":"Intestine","kidney":"Kidney","kidney_tx":"Kidney","lung":"Lung"}
def si(n,adj):
    sx=OPS[n][2]; g=X["g_"+n].values; m=np.isin(g,[0,2])&(X.egfr.notna().values if adj else True)
    S=X[m]; Cm=design(S,sx is not None,["egfr_z","egfr_z2","egfr_z3"] if adj else ())
    e=(g==2)[m].astype(float); Y=G[organs].values[m]
    Q,_=np.linalg.qr(Cm); er=e-Q@(Q.T@e); Yr=Y-Q@(Q.T@Y); den=er@er; B=(er@Yr)/den
    res=Yr-np.outer(er,B); Sig=(res.T@res)/(len(e)-Cm.shape[1]-1); V=Sig/den
    c=np.full(len(organs),-1/(len(organs)-1)); c[organs.index(TGT[n])]=1.0
    est=c@B; se=np.sqrt(c@V@c); return int(e.sum()),B[organs.index(TGT[n])],np.mean(np.delete(B,organs.index(TGT[n]))),est,se,2*stats.norm.sf(abs(est/se))
rows=[]
for n in TGT:
    for adj in (False,True):
        ne,bt,bo,est,se,p=si(n,adj); rows.append([n,TGT[n],"调整 eGFR" if adj else "主模型",ne,bt,bo,est,est-1.96*se,est+1.96*se,p])
D=pd.DataFrame(rows,columns=["op","target_clock","model","n_removed","beta_target","mean_beta_other_organs","specificity_index","lo","hi","p"])
D.to_csv("../results/m09_clock_specificity_index.csv",index=False); print(D.round(3).to_string(index=False))
