import subprocess as subp
import time

t1=time.time()
fq='INDEX3.fasta'
fqfile=open(fq,'w')

def CreateFq(read,fq,gap,length,start,end,num):
    global fqfile
    DNA=open(read,'r')
    DNAlist=[]
    DNAstr=''
    for bp in DNA:
        DNAlist.append(bp.strip())
    for bp in range(len(DNAlist)):
        if bp!=0:
            DNAstr+=str(DNAlist[bp])
    if end=='END':
        end=len(DNAstr)+1
    for bp in range(start//gap-1,end//gap-1):
        if bp==start//gap-1:
            fqfile.write('>@'+str(num)+'_bp_'+str(start)+'_to_'+str(end)+'\n')
        fqfile.write(DNAstr[int(bp*gap+(start-1)%gap):int(bp*gap+(start-1)%gap+length)])
        fqfile.write('\n')
        #fqfile.write('~'*length+'\n')

CreateFq('sequence.fasta',fq,gap=60,length=60,start=1,end=120567768,num=1)
CreateFq('sequence.fasta',fq,gap=60,length=60,start=127687768,end='END',num=1)
CreateFq('chr3.fasta',fq,gap=60,length=60,start=1,end='END',num=3)
CreateFq('sequence2.fasta',fq,gap=60,length=60,start=1,end='END',num=2)
CreateFq('sequence4.fasta',fq,gap=60,length=60,start=1,end='END',num=4)
CreateFq('sequence5.fasta',fq,gap=60,length=60,start=1,end='END',num=5)
CreateFq('sequence6.fasta',fq,gap=60,length=60,start=1,end='END',num=6)
CreateFq('sequence7.fasta',fq,gap=60,length=60,start=1,end='END',num=7)
CreateFq('sequence8.fasta',fq,gap=60,length=60,start=1,end='END',num=8)
CreateFq('sequence9.fasta',fq,gap=60,length=60,start=1,end='END',num=9)
CreateFq('sequence10.fasta',fq,gap=60,length=60,start=1,end='END',num=10)
CreateFq('sequence11.fasta',fq,gap=60,length=60,start=1,end='END',num=11)
CreateFq('sequence12.fasta',fq,gap=60,length=60,start=1,end='END',num=12)
CreateFq('sequence13.fasta',fq,gap=60,length=60,start=1,end='END',num=13)
CreateFq('sequence14.fasta',fq,gap=60,length=60,start=1,end='END',num=14)
CreateFq('sequence15.fasta',fq,gap=60,length=60,start=1,end='END',num=15)
CreateFq('sequence16.fasta',fq,gap=60,length=60,start=1,end='END',num=16)
CreateFq('sequence17.fasta',fq,gap=60,length=60,start=1,end='END',num=17)
CreateFq('sequence18.fasta',fq,gap=60,length=60,start=1,end='END',num=18)
CreateFq('sequence19.fasta',fq,gap=60,length=60,start=1,end='END',num=19)
CreateFq('sequence20.fasta',fq,gap=60,length=60,start=1,end='END',num=20)
CreateFq('sequence21.fasta',fq,gap=60,length=60,start=1,end='END',num=21)
CreateFq('sequence22.fasta',fq,gap=60,length=60,start=1,end='END',num=22)
CreateFq('sequenceX.fasta',fq,gap=60,length=60,start=1,end='END',num=23)

print('\n======= Finished! Time consumed:',time.time()-t1,'s\n>>> New file saved in',fq,':')
a=subp.call(['head',fq])
fqfile.close()