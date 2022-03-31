#import sys
import matplotlib.pyplot as plt
import seaborn as sns
import pandas
import numpy
sns.set()

from Bio import pairwise2
from Bio.pairwise2 import format_alignment

def getinfo(INFO,KEY,caster=float):
    splitter=KEY+"="
    if splitter in INFO:
        value=INFO.split(splitter)[1].split(";")[0]
        return caster(value)
    else:
        return None

def get_scores(ids_intersect,base,call):
    base=parse_vcf(base, ids=ids_intersect)
    call=parse_vcf(call)
    pairs=[]
    scores=[]
    for matchid in base:
        pairs.append((base[matchid],call[matchid]))
        curr_score=pairwise2.align.globalxx(base[matchid], call[matchid], one_alignment_only=True, score_only=True)
        scores.append(curr_score/float(len(base[matchid])))
    return scores

def parse_vcf(filename,ids=None):
    matchtab={}
    for line in open(filename):
        if line[0]=="#":
            continue
        CHROM,POS,ID,REF,ALT,QUAL,FILTER,INFO,FORMAT,SAMPLE=line.split("\t")
        svtype=getinfo(INFO,"SVTYPE",caster=str)
        if svtype!="INS":
            continue
        if ids!=None and not ID in ids:
            continue
        matchtab[getinfo(INFO,"MatchId")]=ALT
    return matchtab

def do_plot(tech,coverage,title,csv_handle):
    id_sets=list()

    #Get the intersection of all caller's TP INS calls (ensures a comparable benchmark of INS sequence accuracy on the same set of SVs for all callers)
    for caller in ["sniffles2","sniffles1","cutesv","svim","pbsv"]:
        #param="def" if caller=="sniffles2" else "sns"
        param="def"
        ids=set()
        for line in open(f"data/pctsim/gb{tech}_{coverage}_{caller}_{param}_base.vcf","r"):
            if line[0]=="#":
                continue
            CHROM,POS,ID,REF,ALT,QUAL,FILTER,INFO,FORMAT,SAMPLE=line.split("\t")
            svtype=getinfo(INFO,"SVTYPE",caster=str)
            if svtype!="INS":
                continue
            ids.add(ID)
        id_sets.append(ids)

    ids_intersect=id_sets[0]
    for curr_set in id_sets:
        ids_intersect=ids_intersect&curr_set

    rows=[]
    callers=["sniffles2","sniffles1","cutesv","svim","pbsv"]
    for caller in callers:
        param="def"
        scores=get_scores(ids_intersect,f"data/pctsim/gb{tech}_{coverage}_{caller}_{param}_base.vcf",f"data/pctsim/gb{tech}_{coverage}_{caller}_{param}_call.vcf")
        for s in scores:
            rows.append([caller,s])

    plt.clf()
    plt.close()
    df=pandas.DataFrame(rows,columns=["caller","score"])
    sns.boxplot(data=df,x="caller",y="score")
    plt.title(title)
    plt.savefig(f"out_pdf/ins_seq_accuracy_{tech}_{coverage}.pdf")
    plt.savefig(f"out_png/ins_seq_accuracy_{tech}_{coverage}.png")

    tech_title={"ont":"ONT","hifi":"HiFi"}

    for caller in callers:
        caller_scores=[row[1] for row in rows if row[0]==caller]
        csv_handle.write(f"HG002,GIAB_Tier1,{tech_title[tech]},{coverage},{caller},{numpy.mean(caller_scores):.4f},{numpy.median(caller_scores):.4f}\n")

out_csv=open(f"out_csv/ins_seq_accuracy.csv","w")
out_csv.write("#Sample,GoldStandard,SequencingTechnology,Coverage,Caller,MeanNormalizedScore,MedianNormalizedScore\n")
do_plot(tech="hifi",coverage="30x", title="Insertion sequence accuracy (HiFi 30x, HG002 GIAB)", csv_handle=out_csv)
do_plot(tech="ont",coverage="30x", title="Insertion sequence accuracy (ONT 30x, HG002 GIAB)", csv_handle=out_csv)
