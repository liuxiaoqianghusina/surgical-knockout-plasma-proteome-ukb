"""器官时钟检验：用 Oh 2025 Nat Med 的 UKB Olink 器官时钟权重计算器官年龄，检验切除器官后对应器官年龄是否改变。
用法：python3 -P 13_organ_clocks.py <Oh SuppTables.xlsx> <Wang SuppTables.xlsx>"""
import pandas as pd, numpy as np, sys
from scipy import stats
sys.path.insert(0,".")
from b_common import *
fo,fw=sys.argv[1:3]
X=pd.read_pickle(f"{C}/atlas_cohort.pkl"); R=pd.read_pickle(f"{C}/atlas_main.pkl")
ol=load_olink().set_index("eid").loc[X.eid]
W=pd.read_excel(fo,sheet_name="ST3").set_index("organ"); icpt=W.pop("intercept")
Z=(ol-ol.mean())/ol.std()
ages={}
for org in W.index:
    w=W.loc[org].dropna(); w=w[w!=0]; w=w[w.index.isin(Z.columns)]
    ages[org]=icpt[org]+Z[w.index].fillna(0).values@w.values
A=pd.DataFrame(ages,index=X.index)
st4=pd.read_excel(fo,sheet_name="ST4").set_index("organ").iloc[:,0]
print("r(pred, age) ours vs paper:"); print(pd.DataFrame({"ours":A.corrwith(X.age),"paper":st4}).round(3).to_string())
# 年龄差：对年龄(样条近似：年龄+年龄²)与性别回归的残差
D=np.column_stack([np.ones(len(X)),X.age,(X.age-57)**2,X.sex,X.sex*(X.age-57)])
Q,_=np.linalg.qr(D); G=A.values-Q@(Q.T@A.values); G=G/G.std(0)   # 年龄差（SD）
G=pd.DataFrame(G,columns=A.columns,index=X.index); G.to_pickle(f"{C}/organ_agegap_oh.pkl")
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
T.to_csv("../results/b07_clock_effects.csv",index=False)
piv=T[k].pivot(index="op",columns="clock",values="beta_sd")
order=["Organismal","Conventional","Intestine","Immune","Kidney","Lung","Liver","Pancreas","Heart","Artery","Brain","Muscle","Adipose"]
pd.set_option("display.width",250); print("\n器官年龄差变化（SD；采血前切除 vs 从未）"); print(piv[order].round(2).to_string())
# 目标对：手术 → 对应器官时钟
TGT={"colon":"Intestine","colon_right":"Intestine","small_bowel":"Intestine","appendix":"Intestine","spleen":"Immune","tonsil":"Immune","kidney":"Kidney","kidney_tx":"Kidney","lung":"Lung","stomach":"Intestine"}
print("\n目标器官时钟（before vs never / after vs never）")
for o,c in TGT.items():
    r=T[(T.op==o)&(T.clock==c)].set_index("contrast")
    if r.empty: continue
    b=r.loc["before_vs_never"]; a=r.loc["after_vs_never"] if "after_vs_never" in r.index else None
    print(f"{o:12s} → {c:9s}  before {b.beta_sd:+.2f} (p={b.p:.1e})  after {a.beta_sd if a is not None else np.nan:+.2f}")
# 分解：每个蛋白对'器官预测年龄'变化的贡献 = 权重 × 蛋白效应（SD）
print("\n贡献分解（年，前 6 个蛋白）")
for o,c in TGT.items():
    w=W.loc[c].dropna(); w=w[w!=0]; r=R[R.op==o].set_index("protein").reindex(w.index)
    contrib=(w*r.beta).dropna().sort_values(key=np.abs,ascending=False)
    if contrib.empty: continue
    print(f"{o:12s} {c:9s} total {contrib.sum():+.2f} y | "+", ".join(f"{p.upper()} {v:+.2f}" for p,v in contrib.head(6).items()))
# Oh 与 Wang 的器官富集蛋白在切除后下降的比例
s2=pd.read_excel(fo,sheet_name="ST2"); w3=pd.read_excel(fw,sheet_name="ST3"); w3["prot"]=w3.GeneName.str.lower()
print("\n器官富集蛋白中切除后显著下降(FDR<0.05, β≤-0.1)的比例")
for o,c in {**TGT,"prostate_radical":"Male","testis":"Male","uterus":"Female","ovary_bilat":"Female","thyroid":"Thyroid"}.items():
    r=R[R.op==o].set_index("protein"); down=(r.fdr_global<0.05)&(r.beta<=-0.1); up=(r.fdr_global<0.05)&(r.beta>=0.1)
    so=set(s2[(s2.organ==c)&(s2.enriched==True)].prot); sw=set(w3[w3.Tissue==c].prot)
    fo_=down.reindex(list(so)).fillna(False).mean() if so else np.nan; fw_=down.reindex(list(sw)).fillna(False).mean() if sw else np.nan
    uo=up.reindex(list(so)).fillna(False).mean() if so else np.nan
    print(f"{o:16s} {c:9s} Oh enriched n={len(so):3d} down {fo_:.0%} up {uo:.0%} | Wang n={len(sw):3d} down {fw_:.0%}")
