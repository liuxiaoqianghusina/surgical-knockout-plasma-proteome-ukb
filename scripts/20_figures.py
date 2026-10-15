"""投稿版图（Nature Communications）：Arial、6–7 pt、180/136 mm 宽、PDF + 600 dpi TIFF。
用法：python3 -P 20_figures.py <Oh 2025 SuppTables.xlsx>"""
import pandas as pd, numpy as np, sys, io, matplotlib; matplotlib.use("Agg")
import matplotlib.pyplot as plt
from matplotlib.patches import FancyBboxPatch, Rectangle
from matplotlib import gridspec
from PIL import Image
from scipy import stats
sys.path.insert(0,"."); from b_common import *
OH=sys.argv[1]
OUT="../figures/"; A="../results/"
MM=1/25.4
plt.rcParams.update({"font.family":"Arial","font.size":6.5,"axes.titlesize":7,"axes.labelsize":6.5,"xtick.labelsize":6,"ytick.labelsize":6,"legend.fontsize":6,
  "mathtext.fontset":"custom","mathtext.rm":"Arial","mathtext.it":"Arial:italic","mathtext.bf":"Arial:bold","axes.linewidth":0.5,"xtick.major.width":0.5,"ytick.major.width":0.5,
  "xtick.major.size":2,"ytick.major.size":2,"axes.spines.top":False,"axes.spines.right":False,"pdf.fonttype":42,"lines.linewidth":0.8,"patch.linewidth":0.5})
RED,BLUE,GREY,LGREY="#B2182B","#2166AC","#7F7F7F","#D0D0D0"
def save(fig,name):
    fig.savefig(OUT+name+".pdf",bbox_inches="tight")
    buf=io.BytesIO(); fig.savefig(buf,format="png",dpi=600,bbox_inches="tight",facecolor="white"); buf.seek(0)
    im=Image.open(buf).convert("RGB"); im.save(OUT+name+".tiff",format="TIFF",compression="tiff_lzw",dpi=(600,600))
    im.resize((im.size[0]//2,im.size[1]//2)).save(OUT+name+"_preview.png"); plt.close(fig)
def plab(ax,t,x=-0.12,y=1.04): ax.text(x,y,t,transform=ax.transAxes,fontsize=8,fontweight="bold",va="bottom",ha="left")
R=pd.read_pickle(f"{C}/atlas_main.pkl"); T=pd.read_pickle(f"{C}/atlas_comparator.pkl")
def se_from(b,p): z=stats.norm.isf(np.clip(np.asarray(p,float),1e-300,1)/2); return np.abs(b)/np.maximum(z,1e-6)
LAB={**LABEL,"uterus_ovaries_kept":"Hysterectomy (ovaries kept)","ovary_given_uterus":"Oophorectomy | hysterectomy","nc_hernia":"Inguinal hernia repair","nc_varicose":"Varicose vein surgery",
     "nc_carpal":"Carpal tunnel release","nc_cataract":"Cataract extraction","nc_hip":"Hip replacement","nc_knee":"Knee replacement","stomach":"Gastric surgery","colon":"Colectomy (any)",
     "breast":"Mastectomy","lung":"Lung resection","prostate_radical":"Radical prostatectomy","prostate_turp":"TURP","thyroid":"Thyroidectomy","uterus":"Hysterectomy",
     "ovary_bilat":"Bilateral oophorectomy","small_bowel":"Small-bowel resection","kidney":"Nephrectomy","spleen":"Splenectomy","gallbladder":"Cholecystectomy","appendix":"Appendicectomy",
     "tonsil":"Tonsillectomy","testis":"Orchidectomy","vasectomy":"Vasectomy","thyroid_total":"Total thyroidectomy","thyroid_partial":"Partial thyroidectomy",
     "colon_right":"Right hemicolectomy","colon_left":"Left/sigmoid colectomy"}
# ================= 图 1 =================
fig,ax=plt.subplots(figsize=(180*MM,88*MM)); ax.axis("off"); ax.set_xlim(0,180); ax.set_ylim(0,88)
def box(x,y,w,h,t,fc="#F2F2F2",ec="#595959",bold=False,fs=6.5):
    ax.add_patch(FancyBboxPatch((x,y),w,h,boxstyle="round,pad=0.6,rounding_size=1.5",fc=fc,ec=ec,lw=0.6))
    ax.text(x+w/2,y+h/2,t,ha="center",va="center",fontsize=fs,fontweight="bold" if bold else "normal",linespacing=1.3)
def arr(x1,y1,x2,y2): ax.annotate("",xy=(x2,y2),xytext=(x1,y1),arrowprops=dict(arrowstyle="-|>",color="#595959",lw=0.6,mutation_scale=6))
box(2,64,56,18,"UK Biobank baseline, 2006–2010\nOlink Explore 3072, 2,923 proteins\nn = 52,054",fc="#DCE9F5",bold=True)
box(2,36,56,20,"Surgical history\nSelf-reported operations with age\nHospital procedures (OPCS-4) with dates\n24 organ removals; 6 other operations")
box(2,6,56,20,"Repeat Olink at imaging visits\n2014–2022, 1,463 proteins\nOperation between two blood draws",fc="#DCE9F5")
arr(30,64,30,57)
box(62.5,64,50,14,"Removed before blood draw\n(organ absent at sampling)",fc="#F6D5D5",bold=True)
box(62.5,40,50,14,"Removed after blood draw\n(same indication, organ present)",fc="#E6E6E6",bold=True)
box(62.5,16,50,14,"Never operated",fc="white",bold=True)
for y in (71,47,23): arr(58,46,62,y)
items=["Atlas: 30 operations × 2,923 proteins","Organ loss versus indication;\nanatomy within organs","Within-person change\nbefore and after surgery","Tissue origin: HPA RNA specificity","Proteomic organ clocks\n(13 published clocks)"]
for i,t in enumerate(items):
    y=70-i*15.5; box(122,y,56,11,f"{i+1}   "+t.replace(chr(10),chr(10)+"     "),fs=6.3); arr(113,47,121.5,y+5.5)
save(fig,"Fig1")
# ================= 图 2 =================
organ=["tonsil","appendix","gallbladder","stomach","colon","small_bowel","spleen","kidney","thyroid","prostate_radical","testis","vasectomy","uterus","ovary_bilat","breast","lung","nc_cataract","nc_knee","nc_hernia","nc_varicose","nc_carpal"]
sel=[]
for o in organ:
    s=R[(R.op==o)&(R.fdr_global<0.05)].assign(a=lambda d:d.beta.abs()).sort_values("a",ascending=False)
    sel+=[p for p in s.protein.head(3) if p not in sel]
H=R[R.op.isin(organ)&R.protein.isin(sel)].pivot(index="op",columns="protein",values="beta").loc[organ,sel]
Q=R[R.op.isin(organ)&R.protein.isin(sel)].pivot(index="op",columns="protein",values="fdr_global").loc[organ,sel]
fig,ax=plt.subplots(figsize=(180*MM,95*MM))
im=ax.imshow(H.clip(-1,1),cmap="RdBu_r",vmin=-1,vmax=1,aspect="auto",interpolation="nearest")
yy,xx=np.where((Q.values<0.05)&(np.abs(H.values)>=0.2)); ax.scatter(xx,yy,s=1.5,color="black",lw=0)
ax.set_yticks(range(len(organ))); ax.set_yticklabels([LAB[o] for o in organ])
ax.set_xticks(range(len(sel))); ax.set_xticklabels([p.upper() for p in sel],rotation=90,fontsize=5)
ax.axhline(15.5,color="black",lw=0.6); ax.text(len(sel)+0.6,7.5,"Organ removal",rotation=-90,va="center",fontsize=6); ax.text(len(sel)+0.6,18,"Other operations",rotation=-90,va="center",fontsize=6)
for s_ in ax.spines.values(): s_.set_visible(False)
ax.tick_params(length=0)
cb=fig.colorbar(im,ax=ax,shrink=0.5,pad=0.04,aspect=15); cb.set_label("β (SD), removed before blood draw vs never",fontsize=6); cb.ax.tick_params(labelsize=5.5,width=0.5,length=2); cb.outline.set_linewidth(0.5)
save(fig,"Fig2")
# ================= 图 3 =================
pairs=[("prostate_radical","klk3"),("thyroid","tg"),("thyroid","tshb"),("uterus","chrdl2"),("uterus","paep"),("colon","guca2a"),("colon","fgf19"),("spleen","wfikkn1"),("gallbladder","agr2"),("tonsil","tmprss11d"),("stomach","pga4")]
fig,ax=plt.subplots(figsize=(110*MM,85*MM)); ys=[]
for i,(o,p) in enumerate(pairs[::-1]):
    r=T[(T.op==o)&(T.protein==p)].iloc[0]
    for b,pp,c,dy,lab in [(r.beta_after,r.p_after,GREY,0.17,"Removed after blood draw (organ present)"),(r.beta_before,r.p_before,RED,-0.17,"Removed before blood draw (organ absent)")]:
        ax.errorbar(b,i+dy,xerr=1.96*se_from(b,pp),fmt="o",color=c,ms=3,capsize=1.5,elinewidth=0.6,capthick=0.6,label=lab if i==0 else None)
    ys.append(f"{LAB[o]}: {p.upper()}")
ax.axvline(0,color="black",lw=0.5); ax.set_yticks(range(len(pairs))); ax.set_yticklabels(ys); ax.set_xlabel("β (SD) versus never operated")
ax.legend(frameon=False,loc="lower left",bbox_to_anchor=(0,1.0),ncol=1,handletextpad=0.3)
save(fig,"Fig3")
# ================= 图 4 =================
fig=plt.figure(figsize=(180*MM,125*MM)); gs=gridspec.GridSpec(2,4,height_ratios=[1,1.15],hspace=0.7,wspace=0.5)
panels=[("a",["thyroid_total","thyroid_partial"],["tg","tshb","calca","pth"]),("b",["prostate_radical","prostate_turp"],["klk3","msmb","tspan1"]),
        ("c",["uterus_ovaries_kept","ovary_given_uterus"],["chrdl2","lefty2","klk4","rln2","prok1","wfikkn2"]),("d",["colon_right","colon_left"],["fgf19","fabp6","guca2a","pyy","insl5"])]
for k,(L_,ops,ps) in enumerate(panels):
    ax=fig.add_subplot(gs[0,k]); x=np.arange(len(ps)); w=0.38
    for j,(o,c) in enumerate(zip(ops,[RED,BLUE])):
        r=R[R.op==o].set_index("protein").reindex(ps)
        ax.bar(x+(j-0.5)*w,r.beta,w,yerr=1.96*r.se,color=c,error_kw=dict(lw=0.5,capsize=1,capthick=0.5),label=f"{LAB[o]} (n={int(r.n_exp.max())})")
    ax.axhline(0,color="black",lw=0.5); ax.set_xticks(x); ax.set_xticklabels([p.upper() for p in ps],rotation=45,ha="right",fontsize=5.5)
    ax.legend(frameon=False,fontsize=5,loc="lower left",bbox_to_anchor=(-0.05,1.0),handlelength=1); plab(ax,L_,-0.28,1.22)
    if k==0: ax.set_ylabel("β (SD) vs never operated")
G=pd.read_csv(A+"b09_gut_extent.csv").dropna(subset=["FGF19"])
cols=["FGF19","FABP6","GUCA2A","GPA33","INSL5","PYY","REG4","APOA4"]
V=G[cols].apply(lambda s: s.astype(str).str.replace("*","",regex=False).astype(float)); S=G[cols].apply(lambda s: s.astype(str).str.contains(r"\*"))
names={"Total/subtotal colectomy":"Total or subtotal colectomy","Right hemicolectomy (incl. terminal ileum)":"Right hemicolectomy (with terminal ileum)",
       "Left/sigmoid colectomy":"Left or sigmoid colectomy","Rectum only (anterior resection)":"Rectal resection only","Ileal resection (no colon)":"Ileal resection without colectomy"}
ax=fig.add_subplot(gs[1,1:4]); im=ax.imshow(V.clip(-1.2,1.2),cmap="RdBu_r",vmin=-1.2,vmax=1.2,aspect="auto")
for i in range(V.shape[0]):
    for j in range(V.shape[1]): ax.text(j,i,f"{V.iloc[i,j]:+.2f}"+("*" if S.iloc[i,j] else ""),ha="center",va="center",fontsize=5.5,color="white" if abs(V.iloc[i,j])>0.8 else "black")
ax.set_xticks(range(len(cols))); ax.set_xticklabels(cols); ax.set_yticks(range(len(G))); ax.set_yticklabels([f"{names[g]} (n={n})" for g,n in zip(G.group,G.n)])
for s_ in ax.spines.values(): s_.set_visible(False)
ax.tick_params(length=0); plab(ax,"e",-0.62,1.02)
cb=fig.colorbar(im,ax=ax,shrink=0.8,pad=0.02,aspect=12); cb.set_label("β (SD)",fontsize=6); cb.ax.tick_params(labelsize=5.5,width=0.5,length=2); cb.outline.set_linewidth(0.5)
save(fig,"Fig4")
# ================= 图 5 =================
W5=pd.read_csv(A+"l05_ukb_within_person.csv").dropna(subset=["delta_change_SD"]); W5=W5[W5.operation!="prostate_turp"]
order=["prostate_radical","uterus","ovary_bilat","gallbladder","thyroid","nc_knee","nc_hernia","nc_varicose","nc_cataract"]
W5["o"]=W5.operation.map({o:i for i,o in enumerate(order)}); W5=W5.sort_values(["o","protein"]).reset_index(drop=True)
fig,ax=plt.subplots(figsize=(120*MM,110*MM)); y=np.arange(len(W5))[::-1]
ax.errorbar(W5.delta_change_SD,y,xerr=1.96*W5.se,fmt="o",color=RED,ms=3,capsize=1.5,elinewidth=0.6,capthick=0.6,label="Within-person change, operated between two blood draws")
ax.scatter(W5.cross_sectional_beta,y,marker="D",s=9,color=BLUE,zorder=3,label="Cross-sectional estimate, removed before baseline")
ax.axvline(0,color="black",lw=0.5)
ax.set_yticks(y); ax.set_yticklabels([f"{LAB[o]}: {p.upper()} (n={int(n)})" for o,p,n in zip(W5.operation,W5.protein,W5.n_operated_between)],fontsize=5.8)
ax.set_xlabel("Protein change (SD)"); ax.legend(frameon=False,loc="lower left",bbox_to_anchor=(0,1.0),fontsize=5.8)
nc_start=(W5.operation.str.startswith("nc_hernia")).idxmax(); ax.axhline(len(W5)-nc_start-0.5,color=LGREY,lw=0.6,ls="--")
save(fig,"Fig5")
# ================= 图 6 =================
fil=pd.read_csv(A+"b09_filtration_signature.csv",index_col=0).iloc[:,0]
D=pd.read_csv(A+"b06_hpa_enrichment.csv")
keep=["tonsil","spleen","stomach","small_bowel","colon_right","kidney","thyroid","prostate_radical","testis","uterus_ovaries_kept","ovary_given_uterus","lung"]
D=D.set_index("op").loc[keep].reset_index()
fig=plt.figure(figsize=(180*MM,68*MM)); go=gridspec.GridSpec(1,2,width_ratios=[1,1.75],wspace=0.75); gs=gridspec.GridSpecFromSubplotSpec(1,2,subplot_spec=go[1],width_ratios=[1,0.85],wspace=0.08)
ax=fig.add_subplot(go[0]); r=R[R.op=="kidney"].set_index("protein").beta.reindex(fil.index); m=r.notna()&fil.notna()
ax.scatter(fil[m],r[m],s=1.5,alpha=0.4,color=GREY,lw=0,rasterized=True)
for p,c in [("umod",RED),("nphs1",RED),("cst3",BLUE),("b2m",BLUE),("tnfrsf1a",BLUE)]:
    ax.scatter(fil[p],r[p],s=8,color=c,zorder=3); ax.annotate(p.upper(),(fil[p],r[p]),xytext=(3,2),textcoords="offset points",fontsize=5.5,color=c)
ax.axhline(0,color="black",lw=0.4); ax.axvline(0,color="black",lw=0.4)
ax.set_xlabel("Protein change per 1-SD lower eGFR\n(participants without kidney surgery)"); ax.set_ylabel("Protein change after nephrectomy (SD)")
ax.text(0.03,0.97,f"r = {np.corrcoef(fil[m],r[m])[0,1]:.2f}",transform=ax.transAxes,va="top",fontsize=6); plab(ax,"a",-0.22)
ax2=fig.add_subplot(gs[0]); ax3=fig.add_subplot(gs[1],sharey=ax2); yy=np.arange(len(D))[::-1]
nm={"uterus_ovaries_kept":"Hysterectomy (ovaries kept)","ovary_given_uterus":"Oophorectomy | hysterectomy"}
ax2.barh(yy,np.log2(D.OR),color=RED,height=0.65)
for y_,o in zip(yy,D.OR): ax2.text(np.log2(o)+0.08,y_,f"{o:.1f}",va="center",fontsize=5.5)
ax2.set_yticks(yy); ax2.set_yticklabels([nm.get(o,LAB.get(o,o)) for o in D.op]); ax2.set_xlabel("log$_2$ odds ratio\n(organ genes among proteins that fall)"); ax2.set_xlim(0,4.2); plab(ax2,"b",-0.62)
ax2.yaxis.set_label_position("right")
ax3.barh(yy,100*D.frac_enriched_down,color=BLUE,height=0.65)
for y_,f,n in zip(yy,D.frac_enriched_down,D["n_enriched(Olink)"]): ax3.text(100*f+0.8,y_,f"{100*f:.0f}% of {n}",va="center",fontsize=5.5)
ax3.set_xlim(0,42); ax3.set_xlabel("% of organ-enriched proteins\nthat fall after removal"); plt.setp(ax3.get_yticklabels(),visible=False); ax3.tick_params(axis="y",length=0); plab(ax3,"c",-0.08)
save(fig,"Fig6")
# ================= 图 7 =================
CE=pd.read_csv(A+"b07_clock_effects.csv"); CE=CE[CE.contrast=="before_vs_never"]
ops=["tonsil","appendix","spleen","stomach","colon","colon_right","small_bowel","kidney","lung","gallbladder","uterus","ovary_bilat","prostate_radical","thyroid","nc_hernia","nc_varicose"]
clocks=["Intestine","Immune","Kidney","Lung","Liver","Pancreas","Heart","Artery","Brain","Muscle","Adipose","Organismal"]
TG={"tonsil":"Immune","appendix":"Intestine","spleen":"Immune","stomach":"Intestine","colon":"Intestine","colon_right":"Intestine","small_bowel":"Intestine","kidney":"Kidney","lung":"Lung"}
P_=CE.pivot(index="op",columns="clock",values="beta_sd").loc[ops,clocks]; F_=CE.pivot(index="op",columns="clock",values="fdr").loc[ops,clocks]
fig=plt.figure(figsize=(180*MM,95*MM)); go=gridspec.GridSpec(1,2,width_ratios=[1.45,1],wspace=0.72); gr=gridspec.GridSpecFromSubplotSpec(2,1,subplot_spec=go[1],hspace=0.75,height_ratios=[1,1.1])
ax=fig.add_subplot(go[0]); im=ax.imshow(P_.clip(-0.6,0.6),cmap="RdBu_r",vmin=-0.6,vmax=0.6,aspect="auto")
yy,xx=np.where(F_.values<0.05); ax.scatter(xx,yy,s=1.5,color="black",lw=0)
for i,o in enumerate(ops):
    if o in TG: j=clocks.index(TG[o]); ax.add_patch(Rectangle((j-0.5,i-0.5),1,1,fill=False,lw=1.0,ec="black"))
ax.set_xticks(range(len(clocks))); ax.set_xticklabels(clocks,rotation=45,ha="right"); ax.set_yticks(range(len(ops))); ax.set_yticklabels([LAB[o] for o in ops])
for s_ in ax.spines.values(): s_.set_visible(False)
ax.tick_params(length=0); plab(ax,"a",-0.42,1.01)
cb=fig.colorbar(im,ax=ax,shrink=0.5,pad=0.02,aspect=15); cb.set_label("Organ age gap (SD)",fontsize=6); cb.ax.tick_params(labelsize=5.5,width=0.5,length=2); cb.outline.set_linewidth(0.5)
W=pd.read_excel(OH,sheet_name="ST3").set_index("organ"); w=W.loc["Intestine"].drop("intercept").dropna(); w=w[w!=0]
ax=fig.add_subplot(gr[0]); ps=["guca2a","fgf19","mln","cdhr2","reg4","lgals4","ccl25"]
for k,(o,c) in enumerate([("colon_right",RED),("small_bowel",BLUE)]):
    r=R[R.op==o].set_index("protein").reindex(w.index); con=(w*r.beta).dropna()
    ax.barh(np.arange(len(ps))+(k-0.5)*0.38,con.reindex(ps).fillna(0).values,0.38,color=c,label=f"{LAB[o]} (total {con.sum():+.2f} y)")
ax.set_yticks(range(len(ps))); ax.set_yticklabels([p.upper() for p in ps],fontsize=5.5); ax.invert_yaxis(); ax.axvline(0,color="black",lw=0.5)
ax.set_xlabel("Contribution to predicted intestine age (years)"); ax.legend(frameon=False,fontsize=5.3,loc="lower left",bbox_to_anchor=(0,1.0)); plab(ax,"b",-0.3,1.3)
SI=pd.read_csv(A+"m09_clock_specificity_index.csv"); SI=SI[SI.op!="kidney_tx"]
order_si=["tonsil","appendix","spleen","stomach","colon","colon_right","small_bowel","kidney","lung"]
ax=fig.add_subplot(gr[1]); yy=np.arange(len(order_si))[::-1]
for k,(mdl,c,dy) in enumerate([("主模型","black",0.15),("调整 eGFR",GREY,-0.15)]):
    d=SI[SI.model==mdl].set_index("op").loc[order_si]
    ax.errorbar(d.specificity_index,yy+dy,xerr=[d.specificity_index-d.lo,d.hi-d.specificity_index],fmt="o",color=c,ms=2.5,capsize=1.2,elinewidth=0.6,capthick=0.6,
                label="Main model" if k==0 else "Adjusted for eGFR")
ax.axvline(0,color="black",lw=0.5); ax.set_yticks(yy); ax.set_yticklabels([f"{LAB[o]} → {SI[SI.op==o].target_clock.iloc[0]}" for o in order_si],fontsize=5.5)
ax.set_xlabel("Specificity index (SD): own-organ clock\nminus mean of other organ clocks"); ax.legend(frameon=False,fontsize=5.3,loc="lower left",bbox_to_anchor=(0,1.0),ncol=2); plab(ax,"c",-0.62,1.12)
save(fig,"Fig7")
# ================= 补充图 1 扁桃体 =================
P1=pd.read_excel(A+"m03_tonsil_child_final.xlsx",sheet_name="tonsil_child"); H1=pd.read_excel(A+"m03_tonsil_child_final.xlsx",sheet_name="neg_hernia_child")
O=pd.read_csv("../cache/a05_cox.csv"); O=O[O.exposure=="tchild"].set_index("outcome")
fig=plt.figure(figsize=(180*MM,78*MM)); gs=gridspec.GridSpec(1,2,width_ratios=[1,1.15],wspace=0.95)
ax=fig.add_subplot(gs[0]); sig=P1.fdr<0.05
ax.scatter(H1.beta,-np.log10(H1.p),s=1.5,color=LGREY,lw=0,rasterized=True,label="Hernia repair in childhood (control)")
ax.scatter(P1.beta[~sig],-np.log10(P1.p[~sig]),s=1.5,color=GREY,lw=0,rasterized=True)
ax.scatter(P1.beta[sig],-np.log10(P1.p[sig]),s=2.5,color=RED,lw=0,rasterized=True,label="Tonsillectomy in childhood, FDR < 0.05")
off={"cd22":(3,4),"cr2":(3,-6),"fcer2":(3,-4),"tnfrsf13c":(3,3),"ccl20":(-24,3),"fcrl5":(3,-5),"tmprss11d":(3,-2),"slamf7":(-27,-3),"fcar":(-20,2)}
for p in ["tmprss11d","ccl20","fcrl5","ly96","fcar","slamf7","mzb1","cd22","cr2","fcrl1","tnfrsf13c","fcer2","tcl1a"]:
    r=P1[P1.protein==p].iloc[0]; ax.annotate(p.upper(),(r.beta,-np.log10(r.p)),xytext=off.get(p,(3,1)),textcoords="offset points",fontsize=5)
ax.set_xlim(-0.3,0.3); ax.set_xlabel("β (SD), tonsillectomy before age 16 vs never"); ax.set_ylabel("−log$_{10}$ $P$"); ax.legend(frameon=False,fontsize=5.5,loc="lower left",bbox_to_anchor=(0,1.0),markerscale=3); plab(ax,"a",-0.18,1.18)
outs=[("h_pta","Peritonsillar abscess or chronic tonsil disease"),("c_tonsil","Tonsil or oropharyngeal cancer"),("c_nhl","Non-Hodgkin lymphoma"),("c_hl","Hodgkin lymphoma"),
      ("h_pso","Psoriasis"),("h_asthma","Asthma"),("h_ibs","Irritable bowel syndrome"),("h_cd","Crohn's disease"),("h_uc","Ulcerative colitis"),("h_ms","Multiple sclerosis"),
      ("h_pneu","Pneumonia"),("h_koa","Knee osteoarthritis (negative control)"),("h_cataract","Cataract (negative control)")]
ax=fig.add_subplot(gs[1]); yy=np.arange(len(outs))[::-1]
for k,(o,lab) in zip(yy,outs):
    r=O.loc[o]; ax.errorbar(r.HR,k,xerr=[[r.HR-r.lo],[r.hi-r.HR]],fmt="o",color="black",ms=2.5,capsize=1.2,elinewidth=0.6,capthick=0.6)
    ax.text(2.35,k,f"{r.HR:.2f} ({r.lo:.2f}–{r.hi:.2f})   {int(r.events_exp):,}",va="center",fontsize=5.5)
ax.text(2.35,len(outs)-0.3,"HR (95% CI)   Events",fontsize=5.5,fontweight="bold")
ax.axvline(1,color="black",lw=0.5); ax.set_xscale("log"); ax.set_xlim(0.4,2.2); ax.set_xticks([0.5,0.75,1,1.5,2]); ax.set_xticklabels(["0.5","0.75","1","1.5","2"])
ax.minorticks_off(); ax.set_yticks(yy); ax.set_yticklabels([l for _,l in outs]); ax.set_xlabel("Hazard ratio, tonsillectomy before age 16 vs never"); plab(ax,"b",-0.75)
save(fig,"SuppFig1")
# ================= 补充图 2 动力学 =================
kp=[("tonsil","tmprss11d"),("uterus","chrdl2"),("uterus","paep"),("gallbladder","agr2"),("colon","fgf19"),("spleen","wfikkn1")]
fig,axs=plt.subplots(2,3,figsize=(150*MM,85*MM)); fig.subplots_adjust(hspace=0.7,wspace=0.5); lab=["<5","5–14","15–29","≥30"]
for i,(ax,(o,p)) in enumerate(zip(axs.flat,kp)):
    r=T[(T.op==o)&(T.protein==p)].iloc[0]; v=[r.get(c,np.nan) for c in ["beta_y0_5","beta_y5_15","beta_y15_30","beta_y30_90"]]
    ax.plot(range(4),v,"o-",color=RED,ms=3); ax.axhline(0,color="black",lw=0.5); ax.set_xticks(range(4)); ax.set_xticklabels(lab,fontsize=5.5)
    ax.set_title(f"{LAB[o]}: {p.upper()}",fontsize=6); plab(ax,"abcdef"[i],-0.3,1.08)
    if i%3==0: ax.set_ylabel("β (SD) vs never operated")
fig.supxlabel("Years from surgery to blood draw",fontsize=6.5,y=0.0)
save(fig,"SuppFig2")
# ================= 补充图 3 第二套时钟 =================
RCw=pd.read_csv(A+"m13_wang_clock_r.csv"); SIw=pd.read_csv(A+"m13_wang_specificity_index.csv"); SIo=pd.read_csv(A+"m09_clock_specificity_index.csv")
fig=plt.figure(figsize=(170*MM,72*MM)); gs=gridspec.GridSpec(1,2,width_ratios=[0.75,1.45],wspace=0.95)
ax=fig.add_subplot(gs[0]); ax.plot([0.1,1],[0.1,1],color=GREY,lw=0.5,ls="--"); ax.scatter(RCw.r_age_reported_UKB,RCw.r_age_reconstructed,s=7,color="black",zorder=3)
off3={"Lung":(3,1),"Pancreas":(-5,-7),"Heart":(-5,3),"Muscle":(3,-5),"Intestine":(3,-4),"Kidney":(3,-2),"Immune":(-5,3),"Liver":(3,-5),"Artery":(-5,3),"Brain":(3,-5),"Organismal":(-5,3)}
for r in RCw.itertuples():
    dx,dy=off3.get(r.clock,(3,0)); ax.annotate(r.clock,(r.r_age_reported_UKB,r.r_age_reconstructed),xytext=(dx,dy),textcoords="offset points",fontsize=5,ha="right" if dx<0 else "left")
ax.set_xlim(0.1,1); ax.set_ylim(0.1,1); ax.set_xlabel("Correlation with age, reported"); ax.set_ylabel("Correlation with age, rebuilt"); plab(ax,"a",-0.32,1.2)
order3=["tonsil","appendix","spleen","stomach","colon","colon_right","small_bowel","kidney","lung","kidney_tx"]
ax=fig.add_subplot(gs[1]); yy=np.arange(len(order3))[::-1]
for k,(D_,c,dy,lb) in enumerate([(SIo,"black",0.15,"Oh and colleagues (published weights)"),(SIw,RED,-0.15,"Wang and colleagues (rebuilt)")]):
    d=D_[D_.model=="主模型"].set_index("op").loc[order3]
    ax.errorbar(d.specificity_index,yy+dy,xerr=[d.specificity_index-d.lo,d.hi-d.specificity_index],fmt="o",color=c,ms=2.5,capsize=1.2,elinewidth=0.6,capthick=0.6,label=lb)
ax.axvline(0,color="black",lw=0.5); ax.set_yticks(yy); ax.set_yticklabels([f"{LAB[o]} → {SIo[SIo.op==o].target_clock.iloc[0]}" for o in order3],fontsize=5.5)
ax.set_xlabel("Specificity index (SD): own-organ clock\nminus mean of other organ clocks"); ax.legend(frameon=False,fontsize=5.3,loc="lower left",bbox_to_anchor=(0,1.0),ncol=1); plab(ax,"b",-0.62,1.2)
save(fig,"SuppFig3")
print("done")
