"""肠切除范围细分的 OPCS 首次日期。"""
import pandas as pd, time
t0=time.time()
F="/path/to/ukb_data/ukb_health_outcomes.csv"
G={"col_total":("H04","H05","H29"),"col_right":("H06","H07"),"col_transverse":("H08",),"col_left":("H09","H10"),"col_other":("H11",),"rectum":("H33",),
   "jejunum":("G58",),"ileum":("G69",),"ileocaecal":("H07",),"stoma":("G74","H14","H15"),"gastrectomy_total":("G27",),"gastrectomy_partial":("G28",),"nephrectomy_partial":("M03",)}
cols=pd.read_csv(F,nrows=0).columns; dcol=[c for c in cols if c.startswith("p41282_a")]
parts=[]
for ch in pd.read_csv(F,usecols=["eid","p41272"]+dcol,dtype=str,chunksize=25000):
    ch=ch.reset_index(drop=True)
    ex=ch["p41272"].fillna("").str.split("|").explode(); ex=pd.DataFrame({"row":ex.index,"code":ex.values}); ex["pos"]=ex.groupby("row").cumcount(); ex=ex[ex.code.notna()&(ex.code!="")]
    D=ch[dcol]; D.columns=range(D.shape[1]); st=D.stack()
    ex["date"]=pd.to_datetime(st.reindex(pd.MultiIndex.from_arrays([ex.row.values,ex.pos.values])).values,errors="coerce")
    o=pd.DataFrame({"eid":ch.eid.astype(int)})
    for n,pre in G.items(): o["x_"+n]=ex[ex.code.str.startswith(pre)].groupby("row").date.min().reindex(ch.index)
    parts.append(o)
pd.concat(parts,ignore_index=True).to_pickle("../cache/opcs_gut.pkl"); print("done",round(time.time()-t0))
