"""在影像随访 i2 的重复 Olink（约 1,100 人，2014 年后）中复现主图谱信号（|β|≥0.2 且 FDR<0.05）。"""
import pandas as pd, numpy as np
from b_common import *
X=pd.read_pickle(f"{C}/atlas_cohort.pkl"); R=pd.read_pickle(f"{C}/atlas_main.pkl")
o2=pd.read_csv(f"{M}/olink_instance2.csv").rename(columns={"Participant.ID":"eid"})
Z=X.merge(o2,on="eid",suffixes=("","_i2")); print("i2 n",len(Z))
sig=R[(R.fdr_global<0.05)&(R.beta.abs()>=0.2)&R.op.isin(OPS.keys())]
rows=[]
for n,(sr,og,sx) in OPS.items():
    s=sig[(sig.op==n)&sig.protein.isin(o2.columns)]
    if s.empty: continue
    g=Z["g_"+n].values; base=np.ones(len(Z),bool) if sx is None else (Z.sex==sx).values
    m=base&np.isin(g,[0,2]); ne=int((m&(g==2)).sum())
    if ne<10: continue
    Y=Z[s.protein.values].values.astype(float)
    B,SE,P,N,NE=ols_many((g==2)[m].astype(float),design(Z[m],sx is not None),Y[m])
    for pr,b0,b,p in zip(s.protein,s.beta,B,P): rows.append([n,pr,ne,b0,b,p])
T=pd.DataFrame(rows,columns=["op","protein","n_exp_i2","beta_i0","beta_i2","p_i2"])
T["same_dir"]=np.sign(T.beta_i0)==np.sign(T.beta_i2)
T.to_csv("../results/b04_i2_replication.csv",index=False)
print(T.groupby("op").agg(n_pairs=("protein","size"),n_exp=("n_exp_i2","max"),same_dir=("same_dir","mean"),p05=("p_i2",lambda s:(s<0.05).mean()),
   r=("beta_i2",lambda s: np.corrcoef(T.loc[s.index,"beta_i0"],s)[0,1] if len(s)>2 else np.nan)).round(2).to_string())
