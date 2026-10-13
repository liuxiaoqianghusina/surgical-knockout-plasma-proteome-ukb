"""HPA 组织特异性验证：切除器官 X 后下降的蛋白，是否富集于 X 的 HPA 组织富集/增强基因？反之 X 富集蛋白中有多少真的下降？
用法：python3 -P 12_hpa_tissue_specificity.py /path/to/hpa_tissue.tsv"""
import pandas as pd, numpy as np, sys
from scipy import stats
sys.path.insert(0,".")
from b_common import C, LABEL
h=pd.read_csv(sys.argv[1],sep="\t")
def parse(s):
    if not isinstance(s,str): return {}
    out={}
    for kv in s.split(";"):
        if ":" in kv: k,v=kv.rsplit(":",1); out[k.strip().lower()]=float(v)
    return out
h["tis"]=h["RNA tissue specific nTPM"].map(parse); h["spec"]=h["RNA tissue specificity"]; h["cell"]=h["RNA single cell type specific nCPM"].map(parse); h["cspec"]=h["RNA single cell type specificity"]
R=pd.read_pickle(f"{C}/atlas_main.pkl"); T=pd.read_pickle(f"{C}/atlas_comparator.pkl")
olp=R.protein.unique()
gmap={p:[g.upper() for g in p.split("_")] for p in olp}
H=h.set_index(h.Gene.str.upper())
def tissue_set(tc,levels):
    tissues,cells=tc; s=set()
    for p,gs in gmap.items():
        for g in gs:
            if g in H.index:
                r=H.loc[g]; r=r.iloc[0] if isinstance(r,pd.DataFrame) else r
                if r.spec in levels and any(t in r.tis for t in tissues): s.add(p)
                if cells and r.cspec in {"Cell type enriched","Group enriched"} and any(c in r.cell for c in cells) and "Tissue enriched" in levels: s.add(p)
                if cells and r.cspec=="Cell type enhanced" and "Tissue enhanced" in levels and any(c in r.cell for c in cells): s.add(p)
    return s
SQ=["esophagus","tongue","vagina","skin 1"]
MAP={"tonsil":(["lymphoid tissue"]+SQ,["b-cells","plasma cells","basal keratinocytes","suprabasal keratinocytes","esophageal basal cells","esophageal suprabasal cells","esophageal apical cells"]),
 "appendix":(["intestine","lymphoid tissue"],["enterocytes","colonocytes","goblet cells","b-cells"]),"gallbladder":(["gallbladder"],["cholangiocytes"]),
 "spleen":(["lymphoid tissue"],["b-cells","macrophages"]),"stomach":(["stomach 1"],["gastric chief cells","parietal cells","foveolar cells","mucous neck cells","gastric progenitor cells"]),
 "colon":(["intestine"],["enterocytes","colonocytes","goblet cells","paneth cells","tuft cells","enteric stem cells"]),
 "colon_right":(["intestine"],["enterocytes","colonocytes","goblet cells","paneth cells","tuft cells","enteric stem cells"]),
 "colon_left":(["intestine"],["enterocytes","colonocytes","goblet cells","tuft cells"]),"small_bowel":(["intestine"],["enterocytes","goblet cells","paneth cells","tuft cells","enteric stem cells"]),
 "kidney":(["kidney"],["proximal tubule cells","podocytes","distal convoluted tubule cells","loop of henle epithelial cells","renal collecting duct principal cells"]),
 "thyroid":(["thyroid gland"],[]),"thyroid_total":(["thyroid gland"],[]),"parathyroid":(["parathyroid gland"],[]),
 "prostate_radical":(["prostate"],["prostatic glandular cells","basal prostatic cells","prostatic club cells","prostatic hillock cells"]),
 "testis":(["testis"],["leydig cells","sertoli cells","early spermatids","late spermatids","early primary spermatocytes","late primary spermatocytes","peritubular myoid cells"]),
 "vasectomy":(["epididymis","seminal vesicle","testis"],["epididymal principal cells","epididymal clear cells","late spermatids","early spermatids"]),
 "uterus":(["endometrium 1","cervix","smooth muscle"],["endometrial stromal cells","endometrial glandular cells","endometrial luminal cells","endometrial secretory cells","endometrial ciliated cells","decidual stromal cells"]),
 "uterus_ovaries_kept":(["endometrium 1","cervix","smooth muscle"],["endometrial stromal cells","endometrial glandular cells","endometrial luminal cells","endometrial secretory cells","endometrial ciliated cells","decidual stromal cells"]),
 "ovary_bilat":(["ovary","fallopian tube"],["granulosa cells","ovarian stromal cells","oocytes","fallopian secretory cells","fallopian tube ciliated cells"]),
 "ovary_given_uterus":(["ovary","fallopian tube"],["granulosa cells","ovarian stromal cells","oocytes","fallopian secretory cells","fallopian tube ciliated cells"]),
 "breast":(["breast"],["breast secretory cells","breast hormone-responsive cells","breast myoepithelial cells","breast lactating cells"]),
 "lung":(["lung"],["alveolar cells type 1","alveolar cells type 2","respiratory secretory cells","respiratory ciliated cells","transitional alveolar cells"])}
print("HPA tissues sample:",sorted({t for d in h.tis for t in d})[:60])
rows=[]
for op,tis in MAP.items():
    r=R[R.op==op].dropna(subset=["beta"])
    if r.empty: continue
    strict=tissue_set(tis,{"Tissue enriched","Group enriched"}); loose=strict|tissue_set(tis,{"Tissue enhanced"})
    r=r.assign(inS=r.protein.isin(strict),inL=r.protein.isin(loose))
    down=(r.fdr_global<0.05)&(r.beta<=-0.1); up=(r.fdr_global<0.05)&(r.beta>=0.1)
    # 富集：下降蛋白中器官富集比例 vs 背景
    a=(down&r.inL).sum(); b=(down&~r.inL).sum(); c=(~down&r.inL).sum(); d=(~down&~r.inL).sum()
    OR,p=stats.fisher_exact([[a,b],[c,d]]) if a+b>0 else (np.nan,np.nan)
    # 敏感度：严格富集蛋白中下降比例；秩检验
    sens=(down&r.inS).sum()/max(r.inS.sum(),1)
    w=stats.mannwhitneyu(r.loc[r.inS,"beta"],r.loc[~r.inS,"beta"],alternative="less").pvalue if r.inS.sum()>=3 else np.nan
    top=", ".join(r[down&r.inL].sort_values("beta").protein.head(8).str.upper())
    rows.append([op,LABEL.get(op,op),"/".join(tis[0]),r.inS.sum(),r.inL.sum(),int(down.sum()),int(up.sum()),a,round(OR,2) if a+b>0 else np.nan,p,round(sens,2),r.loc[r.inS,"beta"].median(),w,top])
D=pd.DataFrame(rows,columns=["op","label","HPA tissue","n_enriched(Olink)","n_enriched+enhanced","n_down","n_up","down∩tissue","OR","p_fisher","frac_enriched_down","median_beta_enriched","p_rank","top_down_tissue_proteins"])
pd.set_option("display.width",300); print(D.drop(columns=["top_down_tissue_proteins"]).round(3).to_string(index=False))
for _,x in D.iterrows(): print(f"{x.op:20s} {x.top_down_tissue_proteins}")
D.to_csv("../results/b06_hpa_enrichment.csv",index=False)
# 列出每个器官"HPA 认为富集但切除后不降"的蛋白
miss=[]
for op,tis in MAP.items():
    r=R[R.op==op].dropna(subset=["beta"]); strict=tissue_set(tis,{"Tissue enriched"})
    x=r[r.protein.isin(strict)][["protein","beta","p","fdr_global"]].assign(op=op); miss.append(x)
pd.concat(miss).to_csv("../results/b06_tissue_enriched_proteins_by_op.csv",index=False)
