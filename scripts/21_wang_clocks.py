"""第二套器官时钟：按 Wang 等（Nat Aging 2025）补充表 ST4 的各器官时钟蛋白，在本样本中重新估计权重（年龄对蛋白的线性回归，10 折交叉拟合），
与其 ST5 报告的 UKB 相关系数比较；再按与 Oh 时钟相同的流程计算年龄差、切除效应和特异性指数（被切器官时钟 − 其余 9 个器官时钟的平均）。
用法：python3 -P 21_wang_clocks.py <Wang SuppTables.xlsx>"""
import pandas as pd, numpy as np, sys
from scipy import stats
sys.path.insert(0,"."); from b_common import *
fw=sys.argv[1]
X=pd.read_pickle(f"{C}/atlas_cohort.pkl")
ol=load_olink().set_index("eid").loc[X.eid]
s4=pd.read_excel(fw,sheet_name="ST4",header=0)
st5=pd.read_excel(fw,sheet_name="ST5",header=None).iloc[1:,:2].set_index(0)[1].astype(float)
Z=(ol-ol.mean())/ol.std()
rng=np.random.default_rng(2025); fold=rng.integers(0,10,len(X)); age=X.age.values.astype(float)
ages={}; rows=[]
for org in s4.columns:
    ps=s4[org].dropna().astype(str).str.lower().tolist()
    M=np.column_stack([np.ones(len(X)),Z[ps].fillna(0).values])          # 缺失按均值（0）处理，与 Oh 时钟一致
    pred=np.empty(len(X))
    for k in range(10):
        tr=fold!=k; w=np.linalg.lstsq(M[tr],age[tr],rcond=None)[0]; pred[~tr]=M[~tr]@w
    ages[org]=pred; rows.append([org,len(ps),np.corrcoef(pred,age)[0,1],st5.get(org,np.nan)])
A=pd.DataFrame(ages,index=X.index)
RC=pd.DataFrame(rows,columns=["clock","n_proteins","r_age_reconstructed","r_age_reported_UKB"]); print(RC.round(3).to_string(index=False))
RC.to_csv("../results/m13_wang_clock_r.csv",index=False)
# 年龄差（与 Oh 时钟相同：对年龄、年龄²、性别、性别×年龄回归的残差，标准化）
Dm=np.column_stack([np.ones(len(X)),X.age,(X.age-57)**2,X.sex,X.sex*(X.age-57)])
Q,_=np.linalg.qr(Dm); G=A.values-Q@(Q.T@A.values); G=pd.DataFrame(G/G.std(0),columns=A.columns,index=X.index); G.to_pickle(f"{C}/organ_agegap_wang.pkl")
# 切除效应（采血前切除 vs 从未；采血后切除 vs 从未）
rows=[]
for n,(sr,og,sx) in OPS.items():
    g=X["g_"+n].values; base=np.ones(len(X),bool) if sx is None else (X.sex==sx).values
    if (base&(g==2)).sum()<30: continue
    for cmp,(m,e) in {"before_vs_never":(base&np.isin(g,[0,2]),g==2),"after_vs_never":(base&np.isin(g,[0,1]),g==1)}.items():
        if (m&e).sum()<20: continue
        B,SE,P,N,NE=ols_many(e[m].astype(float),design(X[m],sx is not None),G.values[m])
        for org,b,se,p in zip(G.columns,B,SE,P): rows.append([n,cmp,org,int((m&e).sum()),b,se,p])
T=pd.DataFrame(rows,columns=["op","contrast","clock","n_exp","beta_sd","se","p"])
T["fdr"]=np.nan; k=T.contrast=="before_vs_never"; T.loc[k,"fdr"]=stats.false_discovery_control(T.loc[k,"p"])
T.to_csv("../results/m13_wang_clock_effects.csv",index=False)
# 特异性指数（与 m09 相同的方法）
organs=[c for c in G.columns if c!="Organismal"]
scr=X.creat/88.4; fem=(X.sex==0); kk=np.where(fem,0.7,0.9); a=np.where(fem,-0.241,-0.302)
X["egfr"]=142*np.minimum(scr/kk,1)**a*np.maximum(scr/kk,1)**-1.2*0.9938**X.age*np.where(fem,1.012,1)
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
SI=pd.DataFrame(rows,columns=["op","target_clock","model","n_removed","beta_target","mean_beta_other_organs","specificity_index","lo","hi","p"])
SI.to_csv("../results/m13_wang_specificity_index.csv",index=False); print(SI.round(3).to_string(index=False))
piv=T[k].pivot(index="op",columns="clock",values="beta_sd")
pd.set_option("display.width",250); print(piv.loc[[o for o in TGT if o in piv.index]].round(2).to_string())
