"""主图谱：每种手术（采血前切除 vs 从未切除）× 2,923 个 Olink 蛋白。"""
import pandas as pd, numpy as np, time
from scipy import stats
from b_common import *
t=time.time()
X=build(); ol=load_olink(); X=X.merge(ol[["eid"]],on="eid").reset_index(drop=True)
P=ol.set_index("eid").loc[X.eid]; prot=np.array(P.columns); Yall=P.values.astype(float)
X.to_pickle(f"{C}/atlas_cohort.pkl"); print("cohort",len(X),round(time.time()-t),flush=True)
res=[]; summ=[]
def run(name,S_mask,e,sexfix,label):
    S=X[S_mask]; Cm=design(S,sexfix); B,SE,Pv,N,NE=ols_many(e[S_mask].astype(float),Cm,Yall[S_mask])
    r=pd.DataFrame({"op":name,"label":label,"protein":prot,"beta":B,"se":SE,"p":Pv,"n":N,"n_exp":NE}); res.append(r)
    print(f"{name:22s} n_exp={int(e[S_mask].sum()):6d}  min p={np.nanmin(Pv):.1e}  top: "+", ".join(f"{a}({b:+.2f})" for a,b in r.sort_values('p').head(6)[['protein','beta']].values),flush=True)
for n,(sr,og,sx) in OPS.items():
    g=X["g_"+n]; m=g.isin([0,2]).values
    if sx is not None: m&=(X.sex==sx).values
    if ((g==2)&m).sum()<30: print("skip",n,((g==2)&m).sum()); continue
    run(n,m,(g==2).values,sx is not None,LABEL[n])
# 衍生对比
u=X.g_uterus==2; o=X.g_ovary_bilat==2; f=(X.sex==0)
run("uterus_ovaries_kept",(f&((u&~o)|((X.g_uterus==0)&(X.g_ovary_bilat==0)))).values,(u&~o).values,True,"Uterus (ovaries kept)")
run("ovary_given_uterus",(f&u).values,o.values,True,"Ovaries | hysterectomy")
cr=X.g_colon_right==2; cl=X.g_colon_left==2
run("right_vs_left_colon",((cr&~cl)|(cl&~cr)).values,(cr&~cl).values,False,"Right vs left colectomy")
pr=X.g_prostate_radical==2; tu=X.g_prostate_turp==2
run("radical_vs_turp",((pr&~tu)|(tu&~pr)).values&(X.sex==1).values,(pr&~tu).values,True,"Radical prostatectomy vs TURP")
R=pd.concat(res,ignore_index=True)
main=~R.op.str.startswith("nc_")
R["fdr_global"]=np.nan; ok=R.p.notna()
R.loc[ok,"fdr_global"]=stats.false_discovery_control(R.loc[ok,"p"])
ntest=ok.sum(); R["bonf"]=R.p<0.05/ntest
R.to_parquet(f"{C}/atlas_main.parquet") if False else R.to_pickle(f"{C}/atlas_main.pkl")
S=R.groupby(["op","label"]).agg(n_exp=("n_exp","max"),fdr05=("fdr_global",lambda s:(s<0.05).sum()),bonf=("bonf","sum"),max_abs_beta=("beta",lambda s:s.abs().max())).reset_index()
S=S.sort_values("fdr05",ascending=False); print(S.to_string(index=False)); S.to_csv("../results/b02_summary.csv",index=False)
print("tests",ntest,"time",round(time.time()-t))
