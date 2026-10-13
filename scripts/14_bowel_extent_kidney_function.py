"""A：肠切除范围细分；B：肾功能（肌酐 eGFR，避免与胱抑素 C 循环）对肾/甲旁腺/脾/小肠等广泛效应的解释。"""
import pandas as pd, numpy as np, sys
from scipy import stats
sys.path.insert(0,"."); from b_common import *
X=pd.read_pickle(f"{C}/atlas_cohort.pkl"); R=pd.read_pickle(f"{C}/atlas_main.pkl"); G=pd.read_pickle(f"{C}/organ_agegap_oh.pkl")
gx=pd.read_pickle(f"{C}/opcs_gut.pkl"); X=X.merge(gx,on="eid",how="left"); assert len(X)==len(G)
ol=load_olink().set_index("eid").loc[X.eid]; prot=np.array(ol.columns)
# ---------- A 肠切除范围 ----------
bef=lambda c: (X[c]<=X.bdate).values
tot=bef("x_col_total"); rig=bef("x_col_right")&~tot; lef=bef("x_col_left")&~tot&~rig; rec=bef("x_rectum")&~tot&~rig&~lef
anycol=bef("x_col_total")|bef("x_col_right")|bef("x_col_left")|bef("x_col_transverse")|bef("x_col_other")|bef("x_rectum")
ile=bef("x_ileum")&~anycol; jej=bef("x_jejunum")&~anycol&~ile
never=(X.g_colon.values==0)&(X.g_small_bowel.values==0)&~np.column_stack([X[c].notna().values for c in gx.columns if c.startswith("x_")]).any(1)
grp={"Total/subtotal colectomy":tot,"Right hemicolectomy (incl. terminal ileum)":rig,"Left/sigmoid colectomy":lef,"Rectum only (anterior resection)":rec,"Ileal resection (no colon)":ile,"Jejunal resection (no colon)":jej}
P=["fgf19","fabp6","slc51b","fgfr4","guca2a","gpa33","cdhr2","lgals4","reg4","insl5","pyy","mln","gcg","ca1","cdh17","fabp2","apoa4"]
P=[p for p in P if p in ol.columns]; Y=np.column_stack([ol[P].values.astype(float),G[["Intestine","Organismal"]].values])
rows=[]
for k,m in grp.items():
    if m.sum()<10: rows.append([k,int(m.sum())]+[np.nan]*(len(P)+2)); continue
    s=never|m; B,SE,Pv,N,NE=ols_many(m[s].astype(float),design(X[s],False),Y[s])
    rows.append([k,int(m.sum())]+[f"{b:+.2f}{'*' if p<0.05/len(P) else ''}" for b,p in zip(B,Pv)])
A=pd.DataFrame(rows,columns=["group","n"]+[p.upper() for p in P]+["IntestineAge","OrganismalAge"])
pd.set_option("display.width",320); print("A. 肠切除范围（β, SD；* Bonferroni）"); print(A.to_string(index=False)); A.to_csv("../results/b09_gut_extent.csv",index=False)
# ---------- B 肾功能 ----------
scr=X.creat/88.4; fem=(X.sex==0); k=np.where(fem,0.7,0.9); a=np.where(fem,-0.241,-0.302)
X["egfr"]=142*np.minimum(scr/k,1)**a*np.maximum(scr/k,1)**-1.2*0.9938**X.age*np.where(fem,1.012,1)
X["egfr_z"]=(X.egfr-X.egfr.mean())/X.egfr.std(); X["egfr_z2"]=X.egfr_z**2; X["egfr_z3"]=X.egfr_z**3
ok=X.egfr.notna().values
Yall=ol.values.astype(float)
nev=ok&(X.g_kidney.values==0)&(X.g_kidney_tx.values==0)&(X.g_parathyroid.values==0)
Bg,_,_,_,_=ols_many(-X.egfr_z.values[nev],design(X[nev],False),Yall[nev])   # 每降低 1 SD eGFR 的蛋白效应
fil=pd.Series(Bg,index=prot); fil.to_csv("../results/b09_filtration_signature.csv")
print("\nB1. 各手术蛋白效应与'肾滤过签名'(每降 1 SD eGFRcr)的跨蛋白相关")
res=[]
for op in ["kidney","kidney_tx","parathyroid","spleen","small_bowel","colon","stomach","lung","gallbladder","uterus","tonsil","nc_cataract","nc_hip"]:
    r=R[R.op==op].set_index("protein").beta.reindex(prot)
    m=r.notna()&fil.notna(); rr=np.corrcoef(r[m],fil[m])[0,1]
    # eGFR 调整后重估
    sx=OPS[op][2]; g=X["g_"+op].values; base=ok&(np.ones(len(X),bool) if sx is None else (X.sex==sx).values); s=base&np.isin(g,[0,2])
    B2,_,P2,_,_=ols_many((g==2)[s].astype(float),design(X[s],sx is not None,["egfr_z","egfr_z2","egfr_z3"]),Yall[s])
    q2=stats.false_discovery_control(np.nan_to_num(P2,nan=1)); r0=R[R.op==op]
    n0=((r0.fdr_global<0.05)&(r0.beta.abs()>=0.2)).sum(); n1=((q2<0.05)&(np.abs(B2)>=0.2)).sum()
    res.append([op,round(rr,2),n0,n1]); 
    if op in ["kidney","kidney_tx","parathyroid","spleen","small_bowel"]:
        t=pd.Series(B2,index=prot); tq=pd.Series(q2,index=prot)
        top=t[(tq<0.05)].sort_values(key=np.abs,ascending=False).head(8)
        print(f"   {op:12s} r={rr:+.2f}  显著且|β|≥0.2: {n0} → 调整 eGFR 后 {n1} | 调整后最强: "+", ".join(f"{p.upper()}({v:+.2f})" for p,v in top.items()))
print(pd.DataFrame(res,columns=["op","r_with_filtration","n_sig_big_raw","n_sig_big_adj_egfr"]).to_string(index=False))
pd.DataFrame(res,columns=["op","r_with_filtration","n_sig_big_raw","n_sig_big_adj_egfr"]).to_csv("../results/b09_kidney_function.csv",index=False)
# 肾切除：eGFR≥60 分层 + 肾特异蛋白
s=ok&np.isin(X.g_kidney.values,[0,2])&(X.egfr.values>=60)
KP=["umod","nphs1","egf","dpep1","ren","cst3","b2m","tnfrsf1a","igfbp6"]; KP=[p for p in KP if p in ol.columns]
B3,_,P3,_,NE=ols_many((X.g_kidney.values==2)[s].astype(float),design(X[s],False,["egfr_z","egfr_z2"]),ol[KP].values[s].astype(float))
print("\nB2. 肾切除（仅 eGFRcr≥60，且调整 eGFR），n_exp=",NE.max()); print("   "+", ".join(f"{p.upper()} {b:+.2f}(p={pv:.0e})" for p,b,pv in zip(KP,B3,P3)))
# 器官时钟：调整 eGFR
print("\nB3. 器官年龄差（调整 eGFRcr）")
for op in ["kidney","kidney_tx","spleen","small_bowel","colon","lung","stomach"]:
    sx=OPS[op][2]; g=X["g_"+op].values; s=ok&np.isin(g,[0,2])
    B4,_,P4,_,_=ols_many((g==2)[s].astype(float),design(X[s],False,["egfr_z","egfr_z2","egfr_z3"]),G.values[s])
    print(f"   {op:12s} "+", ".join(f"{c} {b:+.2f}" for c,b in zip(G.columns,B4) if c in ["Intestine","Immune","Kidney","Lung","Pancreas","Artery","Organismal"]))
