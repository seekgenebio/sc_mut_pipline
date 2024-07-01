import re
import click
import math

def filter_depth(value, cutoff):
    if float(value) >= float(cutoff):
       return True
    else:
       return False
def fra(value,cutoff):
    if value=="":
        return True
    elif float(value)<=cutoff:
        return True
    else:
        return False  
def getindex(lines, sep, strs):
    line = lines.strip().split(sep)
    return line.index(strs)
def sor_f(rf, rr, af, ar, sorcut):
    rf=float(rf)+0.1
    rr=float(rr)+0.1
    af=float(af)+0.1
    ar=float(ar)+0.1
    R=(rf*ar)/(rr*af)
    refratio=min(rf,rr)/max(rf,rr)
    altratio=min(ar,ar)/max(af,ar)
    syratio=R+1.0/R
    sor=math.log(syratio)+math.log(refratio)-math.log(altratio)
    return True if sor <= float(sorcut) else False
def read_vcf(vcf):
    posdic={}
    with open(vcf,'r') as infile:
        for line in infile:
            if line.startswith("##"):
                if "CSQ" not in line:continue
                formate=re.search(r'Format: (.*?)">', line).group(1)
                infolist=formate.strip().split('|')
                continue
            elif line.startswith("#CHROM"):
                headlist = line.strip().split('\t')
                continue
            ids = line.strip().split('\t')
            chr,pos,ref,alt = ids[0:2]+ids[3:5]
            infoindex = headlist.index('INFO')
            infos = ids[infoindex]
            info = infos.strip().split(';')
            form=info[-1].split("CSQ=")[-1]
            formlist=form.split('|')
            infoname = [i.split('=')[0] for i in info[:-1]]
            infovalue = [i.split('=')[-1] for i in info[:-1]]
            valuelist=[]
            for m in (infoname.index(i) for i in ['DP', 'AO',"SRF","SRR","SAF","SAR"]): 
                valuelist.append(infovalue[m])
            ro, ao, rf, rr, af, ar = valuelist
            eas_1k = formlist[infolist.index('EAS_AF')]
            region=formlist[infolist.index('Feature_type')]
            types= formlist[infolist.index('Consequence')]
            needtype=["frameshift_variant","inframe_deletion","inframe_insertion","missense_variant","start_lost","stop_gained"]
            if ',' in ro or ',' in ao: continue
            if filter_depth(ro, 10) and filter_depth(ao, 3) and fra(eas_1k, 0.01) and sor_f(rf, rr, af, ar, 3) and region=="Transcript" and types in needtype:
                posdic.update({(chr,pos,ref,alt): ""})
    return posdic, infolist, infoindex
@click.command("filter vcf.")
@click.option("--starvcf", help="vcf file obtained using star")
@click.option("--hisatvcf", help="vcf file obtained using hisat2")
@click.option("--finalvcf",  help="candidate variation sites")
def final_vcf(starvcf, hisatvcf, finalvcf):
    minipos,infolist, infoindex = read_vcf(starvcf)
    hisatpos,infolist, infoindex = read_vcf(hisatvcf)
    with open(starvcf, 'r') as infile:
        with open(finalvcf, 'w') as out:
         with open(finalvcf+'.txt', 'w') as outtxt:
           outtxt.write('\t'.join(["chr",'pos','ref']+infolist)+'\n')
           for lines in infile:
               if lines.startswith("#"):
                   out.write(lines)
                   continue
               ids = lines.strip().split('\t')
               chr,pos,ref,alt = ids[0:2]+ids[3:5]
               if (chr,pos,ref,alt) in minipos and (chr,pos,ref,alt) in hisatpos:
                   out.write(lines)
                   outlist=[chr,pos,ref]
                   infos = ids[infoindex]
                   info = infos.strip().split(';')
                   form=info[-1].split("CSQ=")[-1]
                   formlist=form.split('|')
                   outlist+=formlist
                   outtxt.write('\t'.join(outlist)+'\n')
if __name__ == "__main__":
    final_vcf()
