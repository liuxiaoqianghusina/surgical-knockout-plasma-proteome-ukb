"""补齐补充数据中尚未保存的结果：器官时钟（调整 eGFR）、贡献分解、器官富集蛋白比例；扁桃体手术年龄梯度与重复测定复现（正式定义）。
用法：python3 -P 17_clock_and_tonsil_supplementary.py <Oh xlsx> <Wang xlsx>"""
import pandas as pd, numpy as np, sys
from scipy import stats
sys.path.insert(0,"."); from b_common import *
fo,fw=sys.argv[1:3]
X=pd.read_pickle(f"{C}/atlas_cohort.pkl"); R=pd.read_pickle(f"{C}/atlas_main.pkl"); G=pd.read_pickle(f"{C}/organ_agegap_oh.pkl")
scr=X.creat/88.4; fem=(X.sex==0); k=np.where(fem,0.7,0.9); a=np.where(fem,-0.241,-0.302)
X["egfr"]=142*np.minimum(scr/k,1)**a*np.maximum(scr/k,1)**-1.2*0.9938**X.age*np.where(fem,1.012,1)
X["egfr_z"]=(X.egfr-X.egfr.mean())/X.egfr.std(); X["egfr_z2"]=X.egfr_z**2; X["egfr_z3"]=X.egfr_z**3; ok=X.egfr.notna().values
rows=[]
for n,(sr,og,sx) in OPS.items():
    g=X["g_"+n].values; base=ok&(np.ones(len(X),bool) if sx is None else (X.sex==sx).values)
    if (base&(g==2)).sum()<30: continue
    s=base&np.isin(g,[0,2]); B,SE,P,N,NE=ols_many((g==2)[s].astype(float),design(X[s],sx is not None,["egfr_z","egfr_z2","egfr_z3"]),G.values[s])
    for c,b,se,p in zip(G.columns,B,SE,P): rows.append([n,c,int(NE.max()),b,se,p])
pd.DataFrame(rows,columns=["op","clock","n_removed","beta_sd","se","p"]).to_csv("../results/m07_clock_egfr_adjusted.csv",index=False)
W=pd.read_excel(fo,sheet_name="ST3").set_index("organ")
TGT={"colon":"Intestine","colon_right":"Intestine","small_bowel":"Intestine","appendix":"Intestine","stomach":"Intestine","spleen":"Immune","tonsil":"Immune","kidney":"Kidney","kidney_tx":"Kidney","lung":"Lung"}
dec=[]
for o,c in TGT.items():
    w=W.loc[c].drop("intercept").dropna(); w=w[w!=0]; r=R[R.op==o].set_index("protein").reindex(w.index)
    con=(w*r.beta).dropna()
    for p,v in con.items(): dec.append([o,c,p,w[p],r.beta[p],v])
pd.DataFrame(dec,columns=["op","clock","protein","clock_weight","protein_beta_sd","contribution_years"]).to_csv("../results/m07_clock_decomposition.csv",index=False)
s2=pd.read_excel(fo,sheet_name="ST2"); w3=pd.read_excel(fw,sheet_name="ST3"); w3["prot"]=w3.GeneName.str.lower()
fr=[]
for o,c in {**TGT,"prostate_radical":"Male","testis":"Male","uterus":"Female","ovary_bilat":"Female","thyroid":"Thyroid"}.items():
    r=R[R.op==o].set_index("protein"); down=(r.fdr_global<0.05)&(r.beta<=-0.1); up=(r.fdr_global<0.05)&(r.beta>=0.1)
    so=[p for p in s2[(s2.organ==c)&(s2.enriched==True)].prot if p in r.index]; sw=[p for p in w3[w3.Tissue==c].prot if p in r.index]
    fr.append([o,c,len(so),down.reindex(so).mean() if so else np.nan,up.reindex(so).mean() if so else np.nan,len(sw),down.reindex(sw).mean() if sw else np.nan,up.reindex(sw).mean() if sw else np.nan])
pd.DataFrame(fr,columns=["op","organ","n_Oh_enriched","frac_down_Oh","frac_up_Oh","n_Wang_clock_proteins","frac_down_Wang","frac_up_Wang"]).to_csv("../results/m07_clock_enriched_fraction.csv",index=False)
# ---- 扁桃体：手术年龄梯度、早年协变量、重复测定（正式定义）
d=pd.read_pickle(f"{C}/ops.pkl"); d["eid"]=d.eid.astype(int)
codes=[c for c in d.columns if c.startswith("p20004_i0_")]
L=pd.concat([pd.DataFrame({"eid":d.eid,"code":d[c],"age":pd.to_numeric(d[c.replace("p20004","p20011")],errors="coerce")}) for c in codes]); L=L[L.code.notna()]
fa=lambda cs: L[L.code.isin(cs)].groupby("eid").age.min()
ton=fa(["1478","1561"]); ade=fa(["1562"]); app=fa(["1458"])
X["ton_age"]=X.eid.map(ton); X["ade"]=X.eid.isin(ade.index)&~X.eid.isin(ton.index); X["app_child"]=(X.eid.map(app)<20).astype(float)
for b in [1940,1945,1950,1955,1960,1965]: X[f"bc{b}"]=(X.byear>=b).astype(float)
extra=[f"bc{b}" for b in [1940,1945,1950,1955,1960,1965]]+["app_child"]
KEY=["tmprss11d","ccl20","fcrl5","ly96","fcar","slamf7","mzb1","cd22","cr2","fcrl1","tnfrsf13c","fcer2","tcl1a"]
ol=load_olink().set_index("eid").loc[X.eid]; Y=ol[KEY].values.astype(float)
never=X.ton_age.isna().values&~X.ade.values
grad=[]
for lab,lo,hi in [("<5",0,5),("5–9",5,10),("10–15",10,16),("16–19",16,20),("20–39",20,40),("≥40",40,200)]:
    e=(X.ton_age>=lo).values&(X.ton_age<hi).values; m=e|never
    B,SE,P,N,NE=ols_many(e[m].astype(float),design(X[m],False,extra),Y[m])
    for p,b,se,pv in zip(KEY,B,SE,P): grad.append([lab,int(NE.max()),p,b,se,pv])
GR=pd.DataFrame(grad,columns=["age_at_tonsillectomy","n_removed","protein","beta_sd","se","p"]); GR.to_csv("../results/m07_tonsil_age_gradient.csv",index=False)
print(GR[GR.protein.isin(["tmprss11d","ccl20","fcrl5"])].pivot(index="protein",columns="age_at_tonsillectomy",values="beta_sd").round(2))
# 早年协变量敏感性
e10=pd.read_pickle(f"{C}/ops.pkl")[["eid","p1687_i0","p1697_i0","p1647_i0"]]; e10["eid"]=e10.eid.astype(int); Z=X.merge(e10,on="eid",how="left")
Z["s10"]=pd.to_numeric(Z.p1687_i0,errors="coerce"); Z["h10"]=pd.to_numeric(Z.p1697_i0,errors="coerce"); Z["cob"]=Z.p1647_i0
child=(Z.ton_age<16).values; mm=(child|never)&(Z.s10>0).values&(Z.h10>0).values&Z.cob.notna().values
dum=pd.get_dummies(Z[["s10","h10","cob"]].astype(str),drop_first=True).astype(float); ex2=extra+list(dum.columns); Z=pd.concat([Z,dum],axis=1)
B,SE,P,N,NE=ols_many(child[mm].astype(float),design(Z[mm],False,ex2),Y[mm])
EL=pd.DataFrame({"protein":KEY,"beta_sd_early_life_adjusted":B,"se":SE,"p":P,"n_removed":int(NE.max())}); EL.to_csv("../results/m07_tonsil_early_life.csv",index=False); print(EL.round(3).head(4))
# 重复测定（i2）
o2=pd.read_csv(f"{M}/olink_instance2.csv").rename(columns={"Participant.ID":"eid"}).set_index("eid")
k2=[p for p in KEY if p in o2.columns]; S=X[X.eid.isin(o2.index)].copy(); child2=(S.ton_age<16).values; nev2=S.ton_age.isna().values&~S.ade.values; m=child2|nev2
Y2=o2.loc[S.eid,k2].values.astype(float)
B,SE,P,N,NE=ols_many(child2[m].astype(float),design(S[m],False,extra),Y2[m])
RP=pd.DataFrame({"protein":k2,"beta_sd_repeat":B,"se":SE,"p":P,"n_removed":int(NE.max()),"n_total":int(N.max())}); RP.to_csv("../results/m07_tonsil_repeat.csv",index=False); print(RP.round(3))
