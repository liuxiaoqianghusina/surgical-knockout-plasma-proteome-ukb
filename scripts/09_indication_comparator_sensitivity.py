"""对主图谱显著(全局 FDR<0.05)的 手术×蛋白，用'采血后才切除'组分离器官效应与指征效应，并做敏感性分析与时间动力学。"""
import pandas as pd, numpy as np, time
from scipy import stats
from b_common import *
t=time.time()
X=pd.read_pickle(f"{C}/atlas_cohort.pkl"); R=pd.read_pickle(f"{C}/atlas_main.pkl")
ol=load_olink().set_index("eid").loc[X.eid]; prot=np.array(ol.columns)
sig=R[(R.fdr_global<0.05)&(R.op.isin(OPS.keys()))]
out=[]
for n,(sr,og,sx) in OPS.items():
    s=sig[sig.op==n]
    if s.empty: continue
    ps=s.protein.values; Y=ol[ps].values.astype(float)
    base=np.ones(len(X),bool) if sx is None else (X.sex==sx).values
    g=X["g_"+n].values
    def fit(mask,e,extra=()):
        S=X[mask]; Cm=design(S,sx is not None,extra); return ols_many(e[mask].astype(float),Cm,Y[mask])
    r=pd.DataFrame({"op":n,"protein":ps,"beta_before":s.beta.values,"p_before":s.p.values,"n_before":int((base&(g==2)).sum()),"n_after":int((base&(g==1)).sum())})
    if (base&(g==1)).sum()>=20:
        b,se,p,_,_=fit(base&np.isin(g,[0,1]),g==1); r["beta_after"]=b; r["p_after"]=p
        b,se,p,_,_=fit(base&np.isin(g,[1,2]),g==2); r["beta_b_vs_a"]=b; r["p_b_vs_a"]=p
    # 敏感性：+用药+胱抑素C
    Xm=X.copy(); Xm["cysc_f"]=Xm.cysc.fillna(Xm.cysc.median())
    meds=[c for c in X.columns if c.startswith("med_")]
    S=Xm[base&np.isin(g,[0,2])]; Cm=design(S,sx is not None,meds+["cysc_f"])
    b,_,p,_,_=ols_many((g==2)[base&np.isin(g,[0,2])].astype(float),Cm,Y[base&np.isin(g,[0,2])]); r["beta_adj_meds_cysc"]=b; r["p_adj"]=p
    # 敏感性：仅 HES 确认（暴露=HES 采血前；对照=从未）
    h=(X["hes_"+n]==1).values
    if (base&h).sum()>=20:
        m=base&((g==0)|h); b,_,p,_,_=fit(m,h); r["beta_hes_only"]=b; r["p_hes_only"]=p
    # 敏感性：去掉 UKB-PPP 联盟入选者
    m=base&np.isin(g,[0,2])&(X.ppp_sel.values==0); b,_,p,_,_=fit(m,g==2); r["beta_random_subset"]=b
    # 动力学：距手术年数
    yrs=X["yrs_"+n].values
    for lo,hi in [(0,5),(5,15),(15,30),(30,90)]:
        e=(g==2)&(yrs>=lo)&(yrs<hi)
        if (base&e).sum()>=20:
            m=base&((g==0)|e); b,_,p,_,_=fit(m,e); r[f"beta_y{lo}_{hi}"]=b
    out.append(r); print(n,len(ps),round(time.time()-t),flush=True)
T=pd.concat(out,ignore_index=True)
def classify(x):
    if pd.isna(x.get("beta_b_vs_a",np.nan)): return "no comparator"
    nb=x.beta_before; ba=x.beta_after; bba=x.beta_b_vs_a; pba=x.p_b_vs_a
    if pba<0.05 and np.sign(bba)==np.sign(nb) and (abs(ba)<abs(nb)/2 or np.sign(ba)!=np.sign(nb)): return "organ removal"
    if pba>=0.05 and np.sign(ba)==np.sign(nb) and abs(ba)>=abs(nb)/2: return "indication"
    if pba<0.05 and np.sign(bba)==np.sign(nb): return "organ removal + indication"
    return "undetermined"
T["class"]=T.apply(classify,axis=1)
T.to_pickle(f"{C}/atlas_comparator.pkl"); T.to_csv("../results/b03_comparator.csv",index=False)
print(T.groupby(["op","class"]).size().unstack(fill_value=0).to_string())
print("time",round(time.time()-t))
