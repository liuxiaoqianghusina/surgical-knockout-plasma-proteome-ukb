"""儿童期扁桃体切除 × 蛋白组（正式版，与方法一致）：<16 岁切除 vs 从未；排除 16–19 岁手术、成年手术、仅腺样体；
协变量 = 主分析协变量 + 出生年代 + 儿童期阑尾切除；阴性对照暴露：儿童期疝修补、仅腺样体。"""
import pandas as pd, numpy as np, sys
from scipy import stats
sys.path.insert(0,"."); from b_common import *
X=pd.read_pickle(f"{C}/atlas_cohort.pkl")
d=pd.read_pickle(f"{C}/ops.pkl"); d["eid"]=d.eid.astype(int)
codes=[c for c in d.columns if c.startswith("p20004_i0_")]
L=pd.concat([pd.DataFrame({"eid":d.eid,"code":d[c],"age":pd.to_numeric(d[c.replace("p20004","p20011")],errors="coerce")}) for c in codes]); L=L[L.code.notna()]
L.loc[L.age<0,"age"]=np.nan   # 20011：-1 不知道，-3 不愿回答 → 手术年龄未知
fa=lambda cs: L[L.code.isin(cs)].groupby("eid").age.min()
ton=fa(["1478","1561"]); ade=fa(["1562"]); app=fa(["1458"]); her=fa(["1403","1563","1564","1402","1404"])
X["ton_age"]=X.eid.map(ton); X["ade"]=X.eid.isin(ade.index)&~X.eid.isin(ton.index); X["app_child"]=(X.eid.map(app)<20).astype(float)
X["her_child"]=(X.eid.map(her)<16).astype(float); her_unk=(X.eid.isin(her.index)&X.eid.map(her).isna()).values
for b in [1940,1945,1950,1955,1960,1965]: X[f"bc{b}"]=(X.byear>=b).astype(float)
extra=[f"bc{b}" for b in [1940,1945,1950,1955,1960,1965]]+["app_child"]
ol=load_olink().set_index("eid").loc[X.eid]; prot=np.array(ol.columns); Y=ol.values.astype(float)
def run(mask,e,ex):
    B,SE,P,N,NE=ols_many(e[mask].astype(float),design(X[mask],False,ex),Y[mask])
    r=pd.DataFrame({"protein":prot,"beta":B,"se":SE,"p":P,"n":N,"n_exp":NE}).dropna(); r["fdr"]=stats.false_discovery_control(r.p); return r.sort_values("p")
child=(X.ton_age<16).values; never=~X.eid.isin(ton.index).values&~X.ade.values   # 手术年龄未知者两组都不纳入
main=run(child|never,child,extra)
rng=np.random.default_rng(7); half=rng.integers(0,2,len(X)).astype(bool)
h1=run((child|never)&half,child,extra); h2=run((child|never)&~half,child,extra)
neg_her=run(never&~her_unk,(X.her_child==1).values,extra)
neg_ade=run(X.ade.values|never,X.ade.values,extra)
print("child n_exp",int(child.sum()),"control",int(never.sum()))
print("main FDR<.05",(main.fdr<.05).sum(),"Bonf",(main.p<0.05/len(main)).sum())
print("neg hernia FDR<.05",(neg_her.fdr<.05).sum(),"n_exp",int(neg_her.n_exp.max()),"| neg adenoid FDR<.05",(neg_ade.fdr<.05).sum(),"n_exp",int(neg_ade.n_exp.max()))
k=["tmprss11d","ccl20","fcrl5","ly96","fcar","slamf7","mzb1","tnfrsf17","cd22","cr2","fcrl1","tnfrsf13c","fcer2","tcl1a"]
m=main.set_index("protein"); a=h1.set_index("protein"); b=h2.set_index("protein"); na=neg_ade.set_index("protein")
print(pd.DataFrame({"beta":m.beta.reindex(k),"p":m.p.reindex(k),"half1":a.beta.reindex(k),"half2":b.beta.reindex(k),"adenoid_only":na.beta.reindex(k)}).round(3).to_string())
sig=main[main.fdr<.05]; print("split-half r among FDR proteins: %.3f"%np.corrcoef(a.beta.reindex(sig.protein),b.beta.reindex(sig.protein))[0,1])
with pd.ExcelWriter("../results/m03_tonsil_child_final.xlsx") as w:
    main.to_excel(w,sheet_name="tonsil_child",index=False); neg_her.to_excel(w,sheet_name="neg_hernia_child",index=False); neg_ade.to_excel(w,sheet_name="neg_adenoid_only",index=False)
