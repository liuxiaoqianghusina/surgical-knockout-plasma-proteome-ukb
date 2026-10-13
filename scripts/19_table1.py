"""表 1：各手术的人数、住院记录确认比例、距采血年数和采血后切除人数。"""
import pandas as pd, numpy as np, sys
sys.path.insert(0,"."); from b_common import *
X=pd.read_pickle(f"{C}/atlas_cohort.pkl")
print("N",len(X),"age %.1f (%.1f)"%(X.age.mean(),X.age.std()),"female %.1f%%"%(100*(X.sex==0).mean()),"white %.1f%%"%(100*X.white.mean()))
rows=[]
for k,(sr,og,sx) in OPS.items():
    base=np.ones(len(X),bool) if sx is None else (X.sex==sx).values
    b=base&(X["g_"+k]==2).values; a=base&(X["g_"+k]==1).values
    if b.sum()<30: continue
    y=X.loc[b,"yrs_"+k]; ya=X.loc[a,"yafter_"+k]
    rows.append([k,LABEL[k],{None:"both",1:"men",0:"women"}[sx],int(b.sum()),int((X.loc[b,"hes_"+k]==1).sum()),f"{y.median():.0f} ({y.quantile(.25):.0f}–{y.quantile(.75):.0f})",
                 int(a.sum()),f"{ya.median():.1f}" if a.sum() else "-",int((base&(X['g_'+k]==0).values).sum())])
T=pd.DataFrame(rows,columns=["op","label","sex","n_before","n_hes_confirmed","yrs_since_median_IQR","n_after","yrs_to_surgery_after","n_never"])
T["pct_hes"]=(100*T.n_hes_confirmed/T.n_before).round(0); T.to_csv("../results/table1_operations.csv",index=False); print(T.to_string(index=False))
