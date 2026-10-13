"""HES OPCS4：各手术组首次日期（含亚码，用于剂量/部位对比）。"""
import pandas as pd, time
t0=time.time()
F="/path/to/ukb_data/ukb_health_outcomes.csv"
G={"tonsil":("F34",),"adenoid":("E20",),"appendix":("H01","H02"),"gallbladder":("J18",),"spleen":("J69",),"stomach":("G27","G28"),"oesophagus":("G01","G02","G03"),
 "colon_right":("H06","H07"),"colon_left":("H09","H10"),"colon_any":("H04","H05","H06","H07","H08","H09","H10","H11","H29","H33"),"rectum":("H33",),
 "jejunum":("G58",),"ileum":("G69",),"pancreas":("J55","J56","J57"),"liver_resect":("J02",),"kidney":("M02",),"kidney_partial":("M03",),"kidney_tx":("M01",),
 "adrenal":("B22",),"thyroid":("B08",),"thyroid_total":("B081",),"thyroid_partial":("B082","B083","B084","B085","B086","B088","B089"),"parathyroid":("B14",),
 "prostate_radical":("M61",),"prostate_turp":("M65",),"testis_bilat":("N05",),"testis_any":("N05","N06"),"vasectomy":("N17",),"uterus":("Q07","Q08"),
 "adnexa_bilat":("Q22",),"adnexa_uni":("Q23",),"breast":("B27",),"lung":("E54",),
 "nc_hernia":("T20",),"nc_cataract":("C71","C72","C75"),"nc_hip":("W37","W38","W39"),"nc_knee":("W40","W41","W42"),"nc_carpal":("A65",),"nc_varicose":("L85","L86","L87","L88")}
cols=pd.read_csv(F,nrows=0).columns; dcol=[c for c in cols if c.startswith("p41282_a")]
parts=[]
for k,ch in enumerate(pd.read_csv(F,usecols=["eid","p41272"]+dcol,dtype=str,chunksize=25000)):
    ch=ch.reset_index(drop=True)
    ex=ch["p41272"].fillna("").str.split("|").explode(); ex=pd.DataFrame({"row":ex.index,"code":ex.values}); ex["pos"]=ex.groupby("row").cumcount(); ex=ex[ex.code.notna()&(ex.code!="")]
    D=ch[dcol]; D.columns=range(D.shape[1]); st=D.stack()
    ex["date"]=pd.to_datetime(st.reindex(pd.MultiIndex.from_arrays([ex.row.values,ex.pos.values])).values,errors="coerce")
    o=pd.DataFrame({"eid":ch.eid.astype(int)})
    for n,pre in G.items(): o["opcs_"+n]=ex[ex.code.str.startswith(pre)].groupby("row").date.min().reindex(ch.index)
    parts.append(o)
pd.concat(parts,ignore_index=True).to_pickle("../cache/opcs_full.pkl"); print("done",round(time.time()-t0))
