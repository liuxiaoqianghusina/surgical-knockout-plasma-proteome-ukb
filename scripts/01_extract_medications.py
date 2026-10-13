"""基线用药（字段 20003，48 个数组位）。"""
import pandas as pd, time
F="/path/to/ukb_data/ukb_assessment_center.csv"
use=["eid"]+[f"p20003_i0_a{i}" for i in range(48)]
t=time.time(); parts=[]
for ch in pd.read_csv(F,usecols=use,dtype=str,chunksize=20000): parts.append(ch)
pd.concat(parts).to_pickle("../cache/meds.pkl"); print("done",round(time.time()-t))
