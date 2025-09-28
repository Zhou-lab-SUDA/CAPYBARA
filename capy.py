"""Author: Naclist
Contact: Zxzzzhu@gmail.com
Update time: 2025/9/20
Version: 1.1.3 Beta

This file is part of Capybara, a Core-snp Assignment PYthon tool for Acinetobacter baumannii, which is totally free to use and redevelop based on the GNU General Public License. Capybara is designed for pathogen detection and is not for profit. For more detailed information on the usage of the program, please refer to the original GNU license text: http://www.gnu.org/licenses/

This version is modified to downsample for metagenomic data.

"""

#!/usr/bin/env python3
from __future__ import annotations
import argparse,json,os,subprocess,sys,tempfile
from collections import defaultdict
from typing import Dict,Iterable,List,Tuple
from configure import exe,esl_list,lineage_snp,variant_snp
try:
    from configure import hc1030_snp,hc_ref
    HAS_HC=True
except Exception:
    hc1030_snp={}
    hc_ref=None
    HAS_HC=False
MINIMAP2=exe.get("minimap2","minimap2")
SAMTOOLS=exe.get("samtools","samtools")
BCFTOOLS=exe.get("bcftools","bcftools")
ESL_REF=esl_list.get("esl_ref")
def run(cmd:str,check:bool=True)->subprocess.CompletedProcess:
    proc=subprocess.run(cmd,shell=True,stdout=subprocess.PIPE,stderr=subprocess.PIPE,text=True)
    if check and proc.returncode!=0: raise RuntimeError(f"Command failed: {cmd}\nSTDERR:\n{proc.stderr}")
    return proc
def first_fasta_header(fasta_path:str)->str:
    with open(fasta_path,"rt") as fh:
        for line in fh:
            if line.startswith(">"): return line[1:].strip().split()[0]
    raise ValueError(f"No FASTA header in {fasta_path}")
def write_bed_for_positions(bed_path:str,chrom:str,positions:Iterable[int])->None:
    with open(bed_path,"wt") as out:
        for p in sorted(set(int(x) for x in positions)): out.write(f"{chrom}\t{p-1}\t{p}\n")
def depth_for_positions(bam:str,ref:str,bed_path:str)->Dict[int,int]:
    with tempfile.TemporaryDirectory() as td:
        out_tsv=os.path.join(td,"dp.tsv")
        cmd=f"{BCFTOOLS} mpileup -a DP -f {ref} -R {bed_path} {bam} -Ou | {BCFTOOLS} query -f '%CHROM\t%POS\t%DP\n' > {out_tsv}"
        run(cmd);dp={}
        with open(out_tsv,"rt") as fh:
            for line in fh:
                _,pos,dp_str=line.rstrip().split("\t")
                try: dp[int(pos)]=int(dp_str)
                except: dp[int(pos)]=0
        return dp
def targeted_call_with_af(bam:str,ref:str,bed:str,out_vcf:str,threads:int)->None:
    cmd=f"{BCFTOOLS} mpileup -Ou -f {ref} -R {bed} -a FORMAT/AD,FORMAT/DP {bam} | {BCFTOOLS} call -mv -Ov -o {out_vcf}"
    run(cmd)
def parse_vcf_with_af(vcf_path:str)->Dict[int,Tuple[str,str,int,int,float]]:
    calls={}
    with open(vcf_path,"rt") as fh:
        for line in fh:
            if line.startswith("#"): continue
            parts=line.rstrip().split("\t");pos=int(parts[1]);ref,alt=parts[3],parts[4];alt0=alt.split(",")[0]
            fmt=parts[8].split(":") if len(parts)>8 else [];smp=parts[9].split(":") if len(parts)>9 else []
            ad_ref=ad_alt=dp=0
            if fmt and smp:
                idx={k:i for i,k in enumerate(fmt)}
                if "AD" in idx and idx["AD"]<len(smp):
                    try: ad=smp[idx["AD"]].split(",");ad_ref=int(ad[0]);ad_alt=int(ad[1]) if len(ad)>1 else 0
                    except: ad_ref,ad_alt=0,0
                if "DP" in idx and idx["DP"]<len(smp):
                    try: dp=int(smp[idx["DP"]])
                    except: dp=ad_ref+ad_alt
            denom=max(1,ad_ref+ad_alt);af=ad_alt/denom
            calls[pos]=(ref,alt0,ad_ref,ad_alt,af)
    return calls
def map_to_bam(query:str,ref:str,threads:int,subsample:float,out_bam:str)->None:
    subs=f"-s {subsample}" if subsample<1.0 else ""
    cmd=f"{MINIMAP2} -a {ref} {query} -t {threads} | {SAMTOOLS} view -bS {subs} - | {SAMTOOLS} sort -@ {threads} -o {out_bam} -"
    run(cmd);run(f"{SAMTOOLS} index {out_bam}")
def samtools_flagstat(bam:str)->Dict[str,int]:
    out=run(f"{SAMTOOLS} flagstat {bam}").stdout;total=mapped=0
    for line in out.splitlines():
        if " in total " in line: total=int(line.split(" ")[0])
        if " mapped (" in line and " primary " not in line: mapped=int(line.split(" ")[0])
    return {"total":total,"mapped":mapped}
def samtools_breadth(bam:str)->float:
    out=run(f"{SAMTOOLS} coverage {bam}").stdout.strip().splitlines()
    for line in out:
        if line.startswith("#"): continue
        parts=line.split()
        if len(parts)>=7:
            try: return float(parts[6])/100.0
            except: pass
    return 0.0
def decide_lineage(lineage_calls:Dict[int,Tuple[str,str,str,int,float]],min_dp:int,min_af:float)->str:
    votes=defaultdict(int)
    for pos,(label,exp_ref,exp_alt,dp,af) in lineage_calls.items():
        if dp>=min_dp and af>=min_af and exp_alt not in (".","N",None): votes[label]+=1
    if not votes: return "NA"
    return max(votes.items(),key=lambda x:x[1])[0]
def decide_sublineage(best_lineage:str,variant_calls:Dict[int,Tuple[str,str,str,int,float]],min_dp:int,min_af:float)->str:
    hits=[];prefix=best_lineage+"."
    for pos,(label,exp_ref,exp_alt,dp,af) in variant_calls.items():
        if not label.startswith(prefix): continue
        if dp>=min_dp and af>=min_af and exp_alt not in (".","N",None): hits.append((label,af))
    if not hits: return "NA"
    hits.sort(key=lambda x:x[1],reverse=True);top_label,top_af=hits[0]
    if len(hits)>1 and hits[1][1]>=top_af*0.95: return "NA"
    return top_label
def evaluate_panel(bam:str,ref:str,panel:Dict[int,Tuple[str,str,str]],min_dp:int,min_af:float)->Tuple[Dict[int,Tuple[str,str,str,int,float]],Dict[int,Tuple[str,str,int,int,float]]]:
    chrom=first_fasta_header(ref)
    with tempfile.TemporaryDirectory() as td:
        bed=os.path.join(td,"panel.bed");write_bed_for_positions(bed,chrom,panel.keys())
        dp_map=depth_for_positions(bam,ref,bed);vcf=os.path.join(td,"panel.vcf");targeted_call_with_af(bam,ref,bed,vcf,threads=1);vcf_calls=parse_vcf_with_af(vcf)
    enriched={}
    for pos,(label,exp_ref,exp_alt) in panel.items():
        ref,alt,ad_ref,ad_alt,af=vcf_calls.get(pos,(".","",0,0,0.0));dp=dp_map.get(pos,ad_ref+ad_alt)
        enriched[pos]=(label,exp_ref,exp_alt,dp,af)
    return enriched,vcf_calls
def main():
    p=argparse.ArgumentParser(description="ESL barcode caller with optional HC1030")
    p.add_argument("-i","--input",required=True);p.add_argument("-o","--output",default="Capy");p.add_argument("-t","--threads",type=int,default=4)
    p.add_argument("-m","--metagenomic",action="store_true");p.add_argument("--subsample",type=float,default=1.0)
    p.add_argument("--min-snp-depth",type=int,default=5);p.add_argument("--min-allele-frac",type=float,default=0.8)
    args=p.parse_args()
    if not os.path.exists(ESL_REF): sys.exit(f"ESL reference not found: {ESL_REF}")
    sample=os.path.basename(args.input)
    with tempfile.TemporaryDirectory() as td:
        bam_esl=os.path.join(td,"esl.sorted.bam");map_to_bam(args.input,ESL_REF,args.threads,args.subsample,bam_esl)
        flag_esl=samtools_flagstat(bam_esl);breadth=samtools_breadth(bam_esl)
        lineage_calls,_=evaluate_panel(bam_esl,ESL_REF,lineage_snp,args.min_snp_depth,args.min_allele_frac)
        variant_calls,_=evaluate_panel(bam_esl,ESL_REF,variant_snp,args.min_snp_depth,args.min_allele_frac)
        lineage_call=decide_lineage(lineage_calls,args.min_snp_depth,args.min_allele_frac)
        sublineage_call="NA"
        if lineage_call!="NA": sublineage_call=decide_sublineage(lineage_call,variant_calls,args.min_snp_depth,args.min_allele_frac)
        gc1,gc2="NA","NA";hc_markers=[]
        if not args.metagenomic and HAS_HC and hc_ref:
            bam_hc=os.path.join(td,"hc1030.sorted.bam");map_to_bam(args.input,hc_ref,args.threads,args.subsample,bam_hc)
            hc_calls,hc_vcf=evaluate_panel(bam_hc,hc_ref,hc1030_snp,args.min_snp_depth,args.min_allele_frac)
            gc1_cnt=gc2_cnt=0
            for pos,(label,exp_ref,exp_alt,dp,af) in hc_calls.items():
                status="lowcov"
                if dp>=args.min_snp_depth and af>=args.min_allele_frac:
                    ref,alt,_ar,_aa,_af=hc_vcf.get(pos,(".","",0,0,0.0));alt0=alt.split(",")[0]
                    if ref==exp_ref and alt0==exp_alt:
                        status="match"
                        if label=="GC1": gc1_cnt+=1
                        elif label=="GC2": gc2_cnt+=1
                hc_markers.append((label,pos,exp_ref,exp_alt,dp,af,status))
            gc1,gc2=str(gc1_cnt),str(gc2_cnt)
        cov={"sample":sample,"mode":"metagenomic" if args.metagenomic else "isolate","lineage_call":lineage_call,"sublineage_call":sublineage_call,"gc1_votes":gc1,"gc2_votes":gc2,"reads_total":flag_esl.get("total",0),"reads_mapped":flag_esl.get("mapped",0),"breadth_ge1x":round(breadth,6),"variant_total":len(variant_snp),"variant_cov":sum(1 for (_lab,_r,_a,dp,_af) in variant_calls.values() if dp>=args.min_snp_depth)}
        with open(f"{args.output}.coverage.json","wt") as outj: json.dump(cov,outj,indent=2)
        with open(f"{args.output}.summary.tsv","wt") as outs:
            outs.write("sample\tmode\tlineage\tsublineage\tGC1_votes\tGC2_votes\treads_total\treads_mapped\tbreadth\tvariant_cov\tvariant_total\n")
            outs.write(f"{sample}\t{cov['mode']}\t{lineage_call}\t{sublineage_call}\t{gc1}\t{gc2}\t{cov['reads_total']}\t{cov['reads_mapped']}\t{cov['breadth_ge1x']}\t{cov['variant_cov']}\t{cov['variant_total']}\n")
        with open(f"{args.output}.markers.tsv","wt") as outm:
            outm.write("panel\tlabel\tpos\tref\talt\tdp\taf\tstatus\n")
            for pos,(label,exp_ref,exp_alt,dp,af) in lineage_calls.items(): outm.write(f"lineage\t{label}\t{pos}\t{exp_ref}\t{exp_alt}\t{dp}\t{af:.3f}\t.\n")
            for pos,(label,exp_ref,exp_alt,dp,af) in variant_calls.items(): outm.write(f"variant\t{label}\t{pos}\t{exp_ref}\t{exp_alt}\t{dp}\t{af:.3f}\t.\n")
            if hc_markers:
                for (label,pos,ref,alt,dp,af,status) in hc_markers: outm.write(f"hc1030\t{label}\t{pos}\t{ref}\t{alt}\t{dp}\t{af:.3f}\t{status}\n")
    print(f"[OK] Wrote: {args.output}.summary.tsv, {args.output}.markers.tsv, {args.output}.coverage.json")
if __name__=="__main__":
    try: main()
    except Exception as e: sys.stderr.write(f"[ERROR] {e}\n");sys.exit(1)
