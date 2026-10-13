"""UKB 个体内前后对照：在两次 Olink 采血之间做手术的人，其蛋白变化 vs 从未做该手术者的变化。
配对：i0→i2（约 2014–2019）或 i2→i3（约 2019–2022），取包住手术日期的那一对。"""
import pandas as pd, numpy as np, sys, statsmodels.api as sm, warnings; warnings.filterwarnings("ignore")
sys.path.insert(0,"."); from b_common import *
X=pd.read_pickle(f"{C}/atlas_cohort.pkl"); v=pd.read_pickle(f"{C}/visit_dates.pkl"); v["eid"]=v.eid.astype(int)
v["d2"]=pd.to_datetime(v.p53_i2,errors="coerce"); v["d3"]=pd.to_datetime(v.p53_i3,errors="coerce")
X=X.merge(v[["eid","d2","d3"]],on="eid",how="left")
o0=load_olink().set_index("eid"); o2=pd.read_csv(f"{M}/olink_instance2.csv").rename(columns={"Participant.ID":"eid"}).set_index("eid")
o3=pd.read_csv(f"{M}/olink_instance3.csv").rename(columns={"Participant.ID":"eid"}).set_index("eid")
TEST={"prostate_radical":(["msmb","tspan1"],1,"down"),"prostate_turp":(["msmb","tspan1"],1,"?"),"uterus":(["chrdl2","lefty2","klk4","paep","rspo1"],0,"down"),
 "ovary_bilat":(["wfikkn2","prok1","chrdl2"],0,"down"),"gallbladder":(["agr2","scgn","anxa10","krt19"],None,"up"),"colon":(["guca2a","fgf19","fabp6"],None,"down"),
 "thyroid":(["calca"],None,"down"),"nc_knee":(["cdcp1","klk4","col9a1"],None,"up"),"nc_hip":(["cdcp1","klk4","col9a1"],None,"up"),
 "nc_hernia":(["cdcp1","agr2","chrdl2"],None,"none"),"nc_varicose":(["cdcp1","agr2","chrdl2"],None,"none"),"nc_cataract":(["cdcp1","agr2"],None,"none")}
op=pd.read_pickle(f"{C}/opcs_full.pkl")
rows=[]
for n,(ps,sx,exp) in TEST.items():
    og=OPS[n][1]; od=pd.concat([op.set_index("eid")["opcs_"+g] for g in og],axis=1).min(axis=1)
    X["od"]=X.eid.map(od); before=X["g_"+n]==2   # 基线前已切除者排除
    recs=[]
    for a,b,oa,ob,da,db in [("i0","i2",o0,o2,"bdate","d2"),("i2","i3",o2,o3,"d2","d3")]:
        S=X[X.eid.isin(oa.index)&X.eid.isin(ob.index)&~before].copy()
        if sx is not None: S=S[S.sex==sx]
        S["ta"]=S[da]; S["tb"]=S[db]; S=S.dropna(subset=["ta","tb"])
        S["between"]=(S.od>S.ta)&(S.od<=S.tb-pd.Timedelta(days=30))
        S["never"]=S.od.isna()|(S.od>S.tb)
        S["yrs"]=(S.tb-S.ta).dt.days/365.25; S["pair"]=f"{a}->{b}"
        for p in ps:
            if p not in oa.columns or p not in ob.columns: continue
            S["d_"+p]=(ob.loc[S.eid,p].values-oa.loc[S.eid,p].values)/np.nanstd(oa.loc[S.eid,p].values)
        recs.append(S[(S.between)|(S.never)])
    R=pd.concat(recs).drop_duplicates("eid",keep="first")   # 每人用一对
    for p in ps:
        c="d_"+p
        if c not in R: continue
        Z=R.dropna(subset=[c]); ne=int(Z.between.sum())
        if ne<3: rows.append([n,p,exp,ne,np.nan,np.nan,np.nan]); continue
        Z=Z.assign(pair2=(Z.pair=="i2->i3").astype(int)); f=f"{c} ~ between + age + yrs + pair2"+("" if sx is not None else " + sex")
        Z=Z.assign(between=Z.between.astype(int))
        r=sm.OLS.from_formula(f,Z).fit()
        rows.append([n,p,exp,ne,r.params["between"],r.bse["between"],r.pvalues["between"]])
T=pd.DataFrame(rows,columns=["operation","protein","expected","n_operated_between","delta_change_SD","se","p"])
cs=pd.read_pickle(f"{C}/atlas_main.pkl").set_index(["op","protein"]).beta
T["cross_sectional_beta"]=[cs.get((o,p),np.nan) for o,p in zip(T.operation,T.protein)]
T.to_csv("../results/l05_ukb_within_person.csv",index=False)
pd.set_option("display.width",200); print(T.round(3).to_string(index=False))
