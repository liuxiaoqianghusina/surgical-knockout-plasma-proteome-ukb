"""补充协变量：空腹时间、HRT、绝经、用药大类、教育年龄、胱抑素C/肌酐、Olink 板号/联盟入选标记。"""
import pandas as pd, time
M="/path/to/ukb_data"
t=time.time()
bs=pd.read_csv(f"{M}/ukb_biological_samples.csv",usecols=["eid","p30720_i0","p30700_i0","p30710_i0","p30901_i0","p30903_i0","p30900_i0"],dtype=str)
print("bio",round(time.time()-t),flush=True)
use=["eid","p74_i0","p2814_i0","p2724_i0","p6153_i0","p6177_i0","p845_i0"]
parts=[]
for k,ch in enumerate(pd.read_csv(f"{M}/ukb_assessment_center.csv",usecols=use,dtype=str,chunksize=25000)): parts.append(ch); print(k,round(time.time()-t),flush=True)
a=pd.concat(parts).merge(bs,on="eid",how="left"); a.to_pickle("../cache/extra_cov.pkl"); print("done",round(time.time()-t))
