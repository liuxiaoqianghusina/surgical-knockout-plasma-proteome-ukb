"""儿童期扁桃体切除(<16 岁) → 基线后新发结局，Cox 预试验。"""
import pandas as pd, numpy as np, warnings, time; warnings.filterwarnings("ignore")
from lifelines import CoxPHFitter
t0=time.time()
d=pd.read_pickle("../cache/ops.pkl"); d["eid"]=d.eid.astype(int)
codes=[c for c in d.columns if c.startswith("p20004_i0_")]
L=pd.concat([pd.DataFrame({"eid":d.eid,"code":d[c],"age":pd.to_numeric(d[c.replace("p20004","p20011")],errors="coerce")}) for c in codes]); L=L[L.code.notna()]
def fa(cs): return L[L.code.isin(cs)].groupby("eid").age.min()
ton=fa(["1478","1561"]); ade=fa(["1562"]); app=fa(["1458"])
X=pd.DataFrame({"eid":d.eid,"sex":pd.to_numeric(d.p31,errors="coerce"),"byear":pd.to_numeric(d.p34,errors="coerce"),"age":pd.to_numeric(d.p21003_i0,errors="coerce"),
 "tdi":pd.to_numeric(d.p189,errors="coerce"),"bmi":pd.to_numeric(d.p21001_i0,errors="coerce"),"smk":pd.to_numeric(d.p20116_i0,errors="coerce"),
 "white":(d.p21000_i0.fillna("").str[:1]=="1").astype(int),"centre":d.p54_i0,"bdate":pd.to_datetime(d.p53_i0,errors="coerce")})
X["ton_age"]=X.eid.map(ton); X["ade_only"]=(X.eid.isin(ade.index)&~X.eid.isin(ton.index)).astype(int); X["app_child"]=(X.eid.map(app)<20).astype(int)
X=X.dropna(subset=["sex","byear","age","tdi","bmi","bdate"]); X=X[(X.smk>=0)&(X.ade_only==0)]
X["tchild"]=(X.ton_age<16).astype(int); X["tadult"]=(X.ton_age>=20).astype(int)
X=X[(X.ton_age.isna())|(X.ton_age<16)|(X.ton_age>=20)]
o=pd.read_pickle("../cache/outcomes.pkl"); X=X.merge(o,on="eid",how="left")
END=pd.Timestamp("2022-12-31")
cov=pd.concat([X[["sex","age","tdi","bmi","white","app_child"]],pd.get_dummies(X.smk.astype(int),prefix="s",drop_first=True),
   pd.get_dummies(pd.cut(X.byear,[1930,1940,1945,1950,1955,1960,1965,1975]).astype(str),prefix="b",drop_first=True),
   pd.get_dummies(X.centre,prefix="c",drop_first=True)],axis=1).astype(float)
outs=[c for c in o.columns if c.startswith(("h_","c_"))]
print("N",len(X),"child",X.tchild.sum(),"adult",X.tadult.sum(),flush=True)
rows=[]
for oc in outs:
    ev=X[oc]; prev=ev<=X.bdate
    stop=pd.concat([ev.where(ev>X.bdate),X.death,pd.Series(END,index=X.index)],axis=1).min(axis=1)
    e=(ev.notna()&(ev>X.bdate)&(ev<=stop)).astype(int)
    T=(stop-X.bdate).dt.days/365.25
    m=(~prev)&(T>0)
    for expo in ["tchild","tadult"]:
        mm=m&((X.tchild==0)&(X.tadult==0)|(X[expo]==1))
        df=pd.concat([T[mm].rename("T"),e[mm].rename("E"),X.loc[mm,[expo]],cov[mm]],axis=1)
        if df.E[df[expo]==1].sum()<5: rows.append([oc,expo,int(df.E.sum()),int(df.E[df[expo]==1].sum()),np.nan,np.nan,np.nan,np.nan]); continue
        try:
            cph=CoxPHFitter(penalizer=1e-4).fit(df,"T","E"); s=cph.summary.loc[expo]
            rows.append([oc,expo,int(df.E.sum()),int(df.E[df[expo]==1].sum()),s["exp(coef)"],s["exp(coef) lower 95%"],s["exp(coef) upper 95%"],s["p"]])
        except Exception as ex: rows.append([oc,expo,int(df.E.sum()),-1,np.nan,np.nan,np.nan,np.nan])
    print(oc,rows[-2][2:],round(time.time()-t0),flush=True)
r=pd.DataFrame(rows,columns=["outcome","exposure","events","events_exp","HR","lo","hi","p"])
r.to_csv("../cache/a05_cox.csv",index=False)
w=r.pivot(index="outcome",columns="exposure",values=["HR","p","events_exp"]); print(w.round(3).to_string())
