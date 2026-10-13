"""HES ICD10 首次日期 + 癌症登记 + 死亡，用于扁桃体切除 → 结局预试验。分块读取。"""
import pandas as pd, numpy as np, time
t0=time.time()
M="/path/to/ukb_data"
F=f"{M}/ukb_health_outcomes.csv"
HES={"pta":("J36","J35"),"igan":("N02",),"pneu":("J12","J13","J14","J15","J16","J17","J18"),"sepsis":("A40","A41"),"uri":("J00","J01","J02","J03","J04","J05","J06"),
 "lri":("J20","J21","J22"),"asthma":("J45","J46"),"copd":("J44",),"crs":("J32","J33"),"cd":("K50",),"uc":("K51",),"coeliac":("K900",),"ibs":("K58",),
 "ms":("G35",),"ra":("M05","M06"),"pso":("L40",),"sle":("M32",),"sjog":("M350",),"t1d":("E10",),"thyai":("E050","E063"),"sarc":("D86",),"as":("M45",),
 "ami":("I21",),"ihd":("I20","I21","I22","I23","I24","I25"),"stroke":("I60","I61","I62","I63","I64"),"dep":("F32","F33"),"anx":("F41",),
 "cataract":("H25","H26"),"koa":("M17",),"vte":("I26","I801","I802","I803"),"gord":("K21",),"divert":("K57",),"gallst":("K80",),"appx":("K35","K36","K37"),
 "hzos":("B02",),"cellul":("L03",),"uti":("N390",),"ckd":("N18",),"ap":("K85",),"hl_h":("C81",),"nhl_h":("C82","C83","C84","C85"),"mm_h":("C90",)}
CAN={"c_tonsil":("C09","C10","C01","C024"),"c_hl":("C81",),"c_nhl":("C82","C83","C84","C85","C86","C88"),"c_cll":("C911",),"c_mm":("C90",),"c_crc":("C18","C19","C20"),
     "c_lung":("C34",),"c_hn_other":("C00","C02","C03","C04","C05","C06","C11","C12","C13","C14","C32"),"c_any":("C",)}
cols=pd.read_csv(F,nrows=0).columns
dcol=[c for c in cols if c.startswith("p41280_a")]
can=[f"p40006_i{i}" for i in range(22)]+[f"p40005_i{i}" for i in range(22)]
use=["eid","p40000_i0","p40001_i0","p41270"]+can+dcol
parts=[]
for k,ch in enumerate(pd.read_csv(F,usecols=use,dtype=str,chunksize=25000)):
    ch=ch.reset_index(drop=True); o=pd.DataFrame({"eid":ch.eid.astype(int)})
    o["death"]=pd.to_datetime(ch.p40000_i0,errors="coerce"); o["death_cause"]=ch.p40001_i0
    for n,pre in CAN.items():
        best=pd.Series(pd.NaT,index=ch.index,dtype="datetime64[ns]")
        for i in range(22):
            m=ch[f"p40006_i{i}"].fillna("").str.startswith(pre); d=pd.to_datetime(ch[f"p40005_i{i}"],errors="coerce")
            best=best.where(~(m&(best.isna()|(d<best))),d)
        o[n]=best
    ex=ch["p41270"].fillna("").str.split("|").explode()
    ex=pd.DataFrame({"row":ex.index,"code":ex.values}); ex["pos"]=ex.groupby("row").cumcount(); ex=ex[ex.code.notna()&(ex.code!="")]
    D=ch[dcol]; D.columns=range(D.shape[1]); st=D.stack()
    ex["date"]=pd.to_datetime(st.reindex(pd.MultiIndex.from_arrays([ex.row.values,ex.pos.values])).values,errors="coerce")
    for n,pre in HES.items():
        s=ex[ex.code.str.startswith(pre)]; o["h_"+n]=s.groupby("row").date.min().reindex(ch.index)
    parts.append(o); print(k,len(o),round(time.time()-t0),flush=True)
pd.concat(parts,ignore_index=True).to_pickle("../cache/outcomes.pkl"); print("done",round(time.time()-t0))
