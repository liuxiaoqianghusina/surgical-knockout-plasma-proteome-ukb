import pandas as pd, time
F="/path/to/ukb_data/ukb_assessment_center.csv"
t=time.time(); parts=[]
for ch in pd.read_csv(F,usecols=["eid","p53_i2","p53_i3","p21003_i2","p21003_i3","p21001_i2","p21001_i3"],dtype=str,chunksize=25000): parts.append(ch)
pd.concat(parts).to_pickle("../cache/visit_dates.pkl"); print("done",round(time.time()-t))
