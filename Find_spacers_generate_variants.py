#####*********Find PAMs and extract flanking sequences**********

def GC_check(seq):
    gct=0
    cct=0
    spacer=seq[15:35]
    for i in range(len(spacer)):
        if spacer[i]=="G":
            gct=gct+1
        if spacer[i]=="C":
            cct=cct+1
    GC=(gct+cct)/len(spacer)
    if GC<0.80 and GC>0.40:
        return(seq)

from Bio import SeqIO
out1=open("spacer_list_mismatch.txt", "w")
out2=open("spacer_list.txt", "w")
with open("Genes.fasta", "r") as handle:
    for record in SeqIO.parse(handle, "fasta"):
        seq=record.seq.rstrip()
        le=len(seq)/10
        count=0
        mark=0
        for i in range(15,len(seq)-38):
            if seq[i]=="T" and mark==0 and count==0:
                if seq[i+1]=="T" or seq[i+1]=="G" or seq[i+1]=="C":
                    new2=seq[i-12:i+38]+"\n"
                    if GC_check(new2) is not None:
                        new2=GC_check(new2)
                        out1.write(str(new2))
                        out2.write(str(new2))
                        mark=i
                        count=count+1
            elif seq[i]=="T" and i>=mark+le and count==1:
                if seq[i+1]=="T" or seq[i+1]=="G" or seq[i+1]=="C":
                    new2=seq[i-12:i+38]+"\n"
                    if GC_check(new2) is not None:
                        new2=GC_check(new2)
                        out1.write(str(new2))
                        out2.write(str(new2))
                        mark=i
                        count=count+1
            elif seq[i]=="T" and i>=mark+le and count>1 and count<5:
                if seq[i+1]=="T" or seq[i+1]=="G" or seq[i+1]=="C":
                    new2=seq[i-12:i+38]+"\n"
                    if GC_check(new2) is not None:
                        new2=GC_check(new2)
                        out2.write(str(new2))
                        mark=i
                        count=count+1

        seq=record.seq.reverse_complement()
        le=len(seq)/10
        count=0
        mark=0
        for i in range(15,len(seq)-38):
            if seq[i]=="T" and mark==0 and count==0:
                if seq[i+1]=="T" or seq[i+1]=="G" or seq[i+1]=="C":
                    new2=seq[i-12:i+38]+"\n"
                    if GC_check(new2) is not None:
                        new2=GC_check(new2)
                        out1.write(str(new2))
                        out2.write(str(new2))
                        mark=i
                        count=count+1
            elif seq[i]=="T" and i>=mark+le and count==1:
                if seq[i+1]=="T" or seq[i+1]=="G" or seq[i+1]=="C":
                    new2=seq[i-12:i+38]+"\n"
                    if GC_check(new2) is not None:
                        new2=GC_check(new2)
                        out1.write(str(new2))
                        out2.write(str(new2))
                        mark=i
                        count=count+1
            elif seq[i]=="T" and i>=mark+le and count>1 and count<5:
                if seq[i+1]=="T" or seq[i+1]=="G" or seq[i+1]=="C":
                    new2=seq[i-12:i+38]+"\n"
                    if GC_check(new2) is not None:
                        new2=GC_check(new2)
                        out2.write(str(new2))
                        mark=i
                        count=count+1
out2.close()
out1.close()

#####*********Generate mismatiched targets**********

def mut(seq):
    import random
    A=list('TGC')
    T=list('AGC')
    G=list('ATC')
    C=list('ATG')
    if seq=="A":
        seq2=random.choice(A)
        return(seq2)
    if seq=="T":
        seq2=random.choice(T)
        return(seq2)
    if seq=="G":
        seq2=random.choice(G)
        return(seq2)
    if seq=="C":
        seq2=random.choice(C)
        return(seq2)
file = 'spacer_list_mismatch.txt'
wr = open("Mut_spacer_file.txt","w")
with open(file) as fp:
    for sp in fp:
        spacer=sp[15:35]
        full=spacer+"TTTTNNNNNNNNNNNNNNNN"+sp
        wr.write(full)
        for i in range(15,35):
            if str(mut(sp[i])) !="None":
                kseq2=str(mut(sp[i]))
                kseq3=sp[0:i]+kseq2+sp[i+1:int(len(sp))]
                full=spacer+"TTTTNNNNNNNNNNNNNNNN"+kseq3
                wr.write(full)
        for i in range(15,34):
            if str(mut(sp[i])) and str(mut(sp[i+1])) !="None":
                kseq2=str(mut(sp[i]))
                kseq21=str(mut(sp[i+1]))
                kseq3=sp[0:i]+kseq2+kseq21+sp[i+2:int(len(sp))]
                full=spacer+"TTTTNNNNNNNNNNNNNNNN"+kseq3
                wr.write(full)
        for i in range(15,32):
            if str(mut(sp[i])) and str(mut(sp[i+1])) and str(mut(sp[i+2])) !="None":
                kseq2=str(mut(sp[i]))
                kseq21=str(mut(sp[i+1]))
                kseq22=str(mut(sp[i+2]))
                kseq23=str(mut(sp[i+3]))
                kseq3=sp[0:i]+kseq2+kseq21+kseq22+kseq23+sp[i+4:int(len(sp))]
                full=spacer+"TTTTNNNNNNNNNNNNNNNN"+kseq3
                wr.write(full)
        for i in range(15,24):
            if str(mut(sp[i])) and str(mut(sp[i+1])) and str(mut(sp[i+10])) and str(mut(sp[i+11])) !="None":
                kseq2=str(mut(sp[i]))
                kseq21=str(mut(sp[i+1]))
                kseq22=str(mut(sp[i+10]))
                kseq23=str(mut(sp[i+11]))
                kseq3=sp[0:i]+kseq2+kseq21+sp[i+2:i+10]+kseq22+kseq23+sp[i+12:len(sp)]
                full=spacer+"TTTTNNNNNNNNNNNNNNNN"+kseq3
                wr.write(full)
wr.close()
fp.close()

#####*********Generate variant Pams**********

import itertools
file = 'spacer_list.txt'
wr = open("Pam_variant_file.txt","w")
with open(file) as fp:
    for sp in fp:
        spacer=sp[15:35]
        for st in itertools.product('ATGC', repeat=2):
            npam="T"+''.join(st)
            spam=sp[0:12]+npam+sp[15:len(sp)]
            fullC=spacer+"TTTTNNNNNNNNNNNNNNNN"+spam
            wr.write(fullC)
wr.close()
fp.close()

###******Concatenating and replacing with barcodes******
cat=open("Final_noncode_list.txt","w")
fnames=["Pam_variant_file.txt", "Mut_spacer_file.txt"]
for fn in fnames:
    with open(fn) as cur:
        cat.write(cur.read())
    cur.close()
cat.close()

import random
with open("barcodes16-2.txt") as bc:
    bcodes=bc.read().splitlines()
    bcodes2=random.sample(bcodes, len(bcodes))
fin = open("Final_oligo_list.txt","w")
cnt=0
with open("Final_noncode_list.txt") as f:
    for oli in f:
        oli = oli.strip("\n")
        Primer1="AACACCGTAATTTCTACTCTTGT"
        Primer2="GCTTGGCGTAACTAGATCTTG"
        repeat="GTCGGAACGCTCAACGATTGCCCCTCACGAGGGGAC"
        new = Primer1 + repeat + oli[0:24] + bcodes2[cnt] + oli[40:len(oli)] + Primer2 + "\n"
        #new2=new+Prime2
        cnt+=1
        fin.write(new)

bc.close()
fin.close()
f.close()

###Error checks
def GC_check(seq):
    gct=0
    cct=0
    spacer=seq[59:81]
    for i in range(len(spacer)):
        if spacer[i]=="G":
            gct=gct+1
        if spacer[i]=="C":
            cct=cct+1
    GC=(gct+cct)/len(spacer)
    if GC<0.80 and GC>0.40:
        return(seq)
error_out = open("Error_out.txt", "w")
count=0
lenct=0
with open("Final_oligo_list.txt") as fp:
    for seq in fp:
        seq=seq.rstrip()
        if GC_check(seq) is None:
            error_out.write("GC content off for:\n")
            error_out.write(seq+"\n")
            count=count+1
    if count==0:
        error_out.write("GC content all OK\n")
print("GC off for "+str(count)+" spacers")
with open("Final_oligo_list.txt") as fp:
    for seq in fp:
        seq=seq.rstrip()
        if len(seq)!= 170:
            lenct=lenct+1
            error_out.write("Length off for:\n")
            error_out.write(seq+"\n")
            lenct=lenct+1
    if lenct==0:
        error_out.write("Lengths all OK\n")
print("Length off for "+str(lenct)+" spacers")
print("Script successfully completed. Check for errors in Error_out.txt")
error_out.close()
