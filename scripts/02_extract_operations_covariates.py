"""抽取基线自报手术(20004/20010/20011) + 出生地/出生年/性别/Townsend 等，用于扁桃体/阑尾切除可行性。"""
import pandas as pd, time
F="/path/to/ukb_data/ukb_assessment_center.csv"
cols=pd.read_csv(F,nrows=0).columns
want=[c for c in cols if c.startswith(("p20004_i0_","p20010_i0_","p20011_i0_"))]
extra=["p31","p34","p52","p21003_i0","p53_i0","p189","p129_i0","p130_i0","p1647_i0","p20022_i0","p1687_i0","p1697_i0","p21001_i0","p20116_i0","p6138_i0_a0","p21000_i0","p54_i0"]
use=["eid"]+want+[c for c in extra if c in cols]
print(len(use),[c for c in extra if c not in cols])
t=time.time(); parts=[]
for k,ch in enumerate(pd.read_csv(F,usecols=use,dtype=str,chunksize=25000)):
    parts.append(ch); print(k,round(time.time()-t),flush=True)
d=pd.concat(parts)
pc=pd.read_csv("/path/to/ukb_data/ukb_population_characteristics.csv",usecols=["eid","p31","p34","p52","p22189"],dtype=str).rename(columns={"p22189":"p189"})
d=d.drop(columns=[c for c in ["p31","p34","p52","p189"] if c in d.columns]).merge(pc,on="eid",how="left")
d.to_pickle("../cache/ops.pkl"); print("done",round(time.time()-t))
