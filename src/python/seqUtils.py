import math
import subprocess
import logging

def make_chunks(w,b,vcf):
    pos = []
    chrom = set()
    cmd="bcftools query -f '%CHROM %POS\\n' "+vcf
#    print(cmd)
    regions=subprocess.check_output(cmd,shell=True)
    regions=regions.decode("utf-8").strip().split('\n')
    logging.info("%d positions found"%len(regions))
    for line in regions:
        c,p = line.strip().split()
        pos.append(int(p))
        chrom.add(c)
    assert len(chrom)==1
    chrom=chrom.pop()
    if (pos[-1]-pos[0]) < w:
        return [(chrom,pos[0],pos[-1])]
    
    r=float(len(pos))/(pos[-1] - pos[0])
    w_m=int(round( r * w))
    b_m=int(round( r * b))
    nwin=round(len(pos)/float(w_m))
    w_m=int(math.ceil(len(pos)/nwin))
    chunks=[(chrom,pos[max(0,val-b_m)],pos[min(len(pos)-1,val+w_m+b_m)]) for val in range(0,len(pos),w_m)]
    if len(chunks)>1 and chunks[-1][2]==chunks[-2][2]:
        del chunks[-1]
    logging.info("split into %d chunks"%len(chunks))
    return chunks
