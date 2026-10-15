"""手术敲除图谱：公共函数（队列构建 + 按缺失模式分组的精确 OLS）。"""
import pandas as pd, numpy as np
from scipy import stats
M="/path/to/ukb_data"
C="/path/to/project/cache"
# name: (自报码, OPCS组, 限定性别 1男/0女/None)
OPS={"tonsil":(["1478","1561"],["tonsil"],None),"adenoid_only":(["1562"],["adenoid"],None),"appendix":(["1458"],["appendix"],None),
 "gallbladder":(["1455"],["gallbladder"],None),"spleen":(["1453"],["spleen"],None),"stomach":(["1448"],["stomach"],None),"oesophagus":([],["oesophagus"],None),
 "colon":(["1464","1461","1465"],["colon_any"],None),"colon_right":([],["colon_right"],None),"colon_left":([],["colon_left"],None),
 "small_bowel":(["1462"],["jejunum","ileum"],None),"pancreas":([],["pancreas"],None),"liver_resect":([],["liver_resect"],None),
 "kidney":(["1487"],["kidney"],None),"kidney_tx":(["1195"],["kidney_tx"],None),"adrenal":(["1232"],["adrenal"],None),
 "thyroid":(["1432"],["thyroid"],None),"thyroid_total":([],["thyroid_total"],None),"thyroid_partial":([],["thyroid_partial"],None),
 "parathyroid":(["1229"],["parathyroid"],None),"prostate_radical":(["1208"],["prostate_radical"],1),"prostate_turp":(["1209"],["prostate_turp"],1),
 "testis":(["1217"],["testis_any"],1),"vasectomy":(["1218"],["vasectomy"],1),"uterus":(["1357","1358","1359"],["uterus"],0),
 "ovary_bilat":(["1355","1358"],["adnexa_bilat"],0),"breast":(["1368"],["breast"],0),"lung":(["1444"],["lung"],None),
 "nc_hernia":(["1403","1563","1564","1402"],["nc_hernia"],None),"nc_cataract":(["1435"],["nc_cataract"],None),"nc_hip":(["1318"],["nc_hip"],None),
 "nc_knee":(["1319"],["nc_knee"],None),"nc_carpal":(["1501"],["nc_carpal"],None),"nc_varicose":(["1479"],["nc_varicose"],None)}
LABEL={"tonsil":"Tonsils","adenoid_only":"Adenoids only","appendix":"Appendix","gallbladder":"Gallbladder","spleen":"Spleen","stomach":"Stomach","oesophagus":"Oesophagus",
 "colon":"Colon (any)","colon_right":"Right colon","colon_left":"Left colon","small_bowel":"Small bowel","pancreas":"Pancreas","liver_resect":"Liver (partial)",
 "kidney":"Kidney","kidney_tx":"Kidney transplant","adrenal":"Adrenal","thyroid":"Thyroid","thyroid_total":"Thyroid (total)","thyroid_partial":"Thyroid (partial)",
 "parathyroid":"Parathyroid","prostate_radical":"Prostate (radical)","prostate_turp":"Prostate (TURP)","testis":"Testis","vasectomy":"Vasectomy","uterus":"Uterus",
 "ovary_bilat":"Ovaries (bilateral)","breast":"Breast (mastectomy)","lung":"Lung","nc_hernia":"Inguinal hernia repair","nc_cataract":"Cataract",
 "nc_hip":"Hip replacement","nc_knee":"Knee replacement","nc_carpal":"Carpal tunnel","nc_varicose":"Varicose veins"}
MEDPAT={"ppi":"omeprazole|lansoprazole|pantoprazole|esomeprazole|rabeprazole","thyroxine":"^thyroxine|sodium thyroxine","hrt":"oestrogen|^oestradiol|premarin|tibolone|conjugated|polyestradiol",
 "five_ari":"finasteride|dutasteride","adt":"bicalutamide|goserelin|leuprorelin|cyproterone acetate$","b12":"cobalamin|b12 ","pancreatin":"pancreatin|creon",
 "bile_seq":"cholestyramine|colestyramine","pred":"^prednisolone|methylprednisolone$","endo_ca":"tamoxifen|anastrozole|letrozole"}

def build():
    d=pd.read_pickle(f"{C}/ops.pkl"); d["eid"]=d.eid.astype(int)
    codes=[c for c in d.columns if c.startswith("p20004_i0_")]
    L=pd.concat([pd.DataFrame({"eid":d.eid,"code":d[c],"age":pd.to_numeric(d[c.replace("p20004","p20011")],errors="coerce")}) for c in codes]); L=L[L.code.notna()]
    L.loc[L.age<0,"age"]=np.nan   # 20011：-1 不知道，-3 不愿回答 → 手术年龄未知
    X=pd.DataFrame({"eid":d.eid,"sex":pd.to_numeric(d.p31,errors="coerce"),"byear":pd.to_numeric(d.p34,errors="coerce"),"age":pd.to_numeric(d.p21003_i0,errors="coerce"),
      "tdi":pd.to_numeric(d.p189,errors="coerce"),"bmi":pd.to_numeric(d.p21001_i0,errors="coerce"),"smk":pd.to_numeric(d.p20116_i0,errors="coerce"),
      "white":(d.p21000_i0.fillna("").str[:1]=="1").astype(float),"centre":d.p54_i0,"bdate":pd.to_datetime(d.p53_i0,errors="coerce")})
    e=pd.read_pickle(f"{C}/extra_cov.pkl"); e["eid"]=e.eid.astype(int)
    for c in ["p74_i0","p30720_i0","p30700_i0","p30903_i0","p845_i0"]: e[c]=pd.to_numeric(e[c],errors="coerce")
    X=X.merge(e.rename(columns={"p74_i0":"fast","p30720_i0":"cysc","p30700_i0":"creat","p30903_i0":"ppp_sel","p30901_i0":"plate","p2814_i0":"hrt_ever","p845_i0":"edu_age"})[
        ["eid","fast","cysc","creat","ppp_sel","plate","hrt_ever","edu_age"]],on="eid",how="left")
    med=pd.read_pickle(f"{C}/meds.pkl"); med["eid"]=med.eid.astype(int)
    cod=pd.read_csv(f"{M}/ukb_codings.csv",dtype=str); c4=cod[cod.coding_name=="data_coding_4"]
    mc=[c for c in med.columns if c.startswith("p20003_i0_a")]; long=med.melt("eid",mc).dropna(); long=long[long.value!=""]
    for k,p in MEDPAT.items():
        cs=set(c4[c4.meaning.str.contains(p,case=False,regex=True,na=False)].code); X["med_"+k]=X.eid.isin(long[long.value.isin(cs)].eid).astype(float)
    op=pd.read_pickle(f"{C}/opcs_full.pkl"); X=X.merge(op,on="eid",how="left")
    X["age_base"]=X.age
    for n,(sr,og,sx) in OPS.items():
        a=L[L.code.isin(sr)].groupby("eid").age.min() if sr else pd.Series(dtype=float)
        od=pd.concat([X["opcs_"+g] for g in og],axis=1).min(axis=1)
        srm=X.eid.isin(a.index); hes_before=od<=X.bdate
        before=srm|hes_before; after=(~before)&(od>X.bdate+pd.Timedelta(days=30))
        X["g_"+n]=np.where(before,2,np.where(after,1,0))   # 0 从未 1 采血后切除 2 采血前切除
        yrs_sr=X.age-X.eid.map(a); yrs_hes=(X.bdate-od).dt.days/365.25
        X["yrs_"+n]=np.where(before,np.fmax(yrs_sr.where(srm),yrs_hes.where(hes_before)),np.nan)
        X["hes_"+n]=hes_before.astype(int); X["yafter_"+n]=np.where(after,(od-X.bdate).dt.days/365.25,np.nan)
    X=X.dropna(subset=["sex","age","tdi","bmi","bdate"]); X=X[X.smk>=0]
    X["fast"]=X.fast.fillna(X.fast.median()).clip(0,24); X["ppp_sel"]=X.ppp_sel.fillna(0)
    return X.reset_index(drop=True)

def design(S, sexfix, extra=()):
    S=S.assign(age2=(S.age-57)**2)
    parts=[S[["age","age2","tdi","bmi","white","fast","ppp_sel"]+list(extra)],pd.get_dummies(S.smk.astype(int),prefix="s",drop_first=True),pd.get_dummies(S.centre.astype(str),prefix="c",drop_first=True)]
    if not sexfix: parts.append(pd.DataFrame({"sex":S.sex,"sex_age":S.sex*(S.age-57),"sex_age2":S.sex*S.age2,"sex_bmi":S.sex*(S.bmi-27)}))
    Z=pd.concat(parts,axis=1).astype(float).values
    return np.column_stack([np.ones(len(Z)),Z])

def ols_many(e, Cm, Y):
    """e: n 暴露；Cm: n×k 协变量；Y: n×p 蛋白（含 NaN）。按缺失模式分组，逐组 QR，精确完整病例 OLS。"""
    n,p=Y.shape; miss=np.isnan(Y)
    keys=pd.Series([hash(miss[:,j].tobytes()) for j in range(p)])
    B=np.full(p,np.nan); SE=np.full(p,np.nan); N=np.zeros(p,int); NE=np.zeros(p,int)
    for key,idx in keys.groupby(keys).groups.items():
        idx=np.array(list(idx)); ok=~miss[:,idx[0]]
        if ok.sum()<100: continue
        Q,_=np.linalg.qr(Cm[ok]); ee=e[ok]; er=ee-Q@(Q.T@ee); den=er@er
        if den<1e-8: continue
        Yo=Y[ok][:,idx]; Yo=(Yo-Yo.mean(0))/Yo.std(0); Yr=Yo-Q@(Q.T@Yo)
        b=(er@Yr)/den; res=Yr-np.outer(er,b); df=ok.sum()-Cm.shape[1]-1
        B[idx]=b; SE[idx]=np.sqrt((res**2).sum(0)/df/den); N[idx]=ok.sum(); NE[idx]=int(ee.sum())
    P=2*stats.norm.sf(np.abs(B/SE))
    return B,SE,P,N,NE

def load_olink():
    ol=pd.read_csv(f"{M}/olink_instance0.csv").rename(columns={"Participant.ID":"eid"})
    return ol
