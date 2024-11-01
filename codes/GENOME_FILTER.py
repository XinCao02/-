import subprocess as subp
import os
import time
DNAlist=[]
t1=time.time()
bplist=[]
counter=0

def FetchRepeat(read,index):
    global counter
    global bplist
    checkerlist=[]
    s1=subp.Popen('/home/xcao/ncbi-blast-2.12.0+/bin/blastn -query '+read+' -out seq.blast -db '+index+' -outfmt "6 qseqid sseqid qcovs length mismatch gapopen qstart qend sstart send evalue bitscore pident" -evalue 1e-12 -num_threads 42',shell=True)
    sleep1=0.07
    d=0
    while True:
        a=[]
        try:
            time.sleep(sleep1)
            a=[]
            test1=open('seq.blast','r')
            a=test1.readlines(1)
            b=len(a)
            if d>2:
                e=subp.call(['cp','blank','seq.blast'])
                break
            else:
                c=1/b
                break
        except:
            a=[]
            sleep1+=0.05
            print('++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++')
            d+=1
    s2=subp.Popen("awk '$3 * $4 >=5700 "+'{'+'print '+'}'+"' seq.blast > test.blast",shell=True)
    sleep2=0.002
    while True:
        try:
            time.sleep(sleep2)
            test2=open('test.blast','r')
            test2.close() 
            break
        except:
            sleep2+=0.001
            print('----------------------------------------------------------------------')
    s3=subp.Popen('wc -l test.blast > checker2',shell=True)
    while True:
        s4=subp.call(['wc','-l','checker2'])
        s5=subp.call(['head','checker2'])
        try:
            sleep3=0.0019
            time.sleep(sleep3)
            #s5=subp.call(['head','checker2'])
            checker=open('checker2','r')
            for i in checker:
                checkerlist.append(i)
            checker.close()
            # checknum=eval(checkerlist[0][0:4])
            # print('checknum is:',checknum)
            if len(checkerlist[0])>14:
                print('len:',len(checkerlist[0]))
                checknum=eval(checkerlist[0][0:4])
                print('checknum is:',checknum)
                if 300<=checknum<=750:
                    goodDNA=open(read,'r')
                    for i in goodDNA:
                        if i[0]=='>':
                            i='>'+str(checknum)+i
                        bplist.append(i)
                    counter=0
            else:
                print('--------------------------------------------------------------------------------------------------------------------------------------------')
                counter+=1
            break
        except:
            sleep3+=0.00005
            print('==================================================================')
    s4=subp.call(['rm','-f','checker2'])
    s5=subp.call(['rm','-f','seq.blast'])
    s6=subp.call(['rm','-f','test.blast'])

def Narrow(output,INDEX,OUTPUT):
    print('>>> Narrow start! Original double line count:')
    s0=os.system('wc -l '+output+'')
    ratefilter="awk '$3 * $4 > 5100 "+"{"+"print}' "+output+"bp1 > "+output+"bp2"
    print('>>> BLASTing large index ... ...')
    s1=os.system('/home/xcao/ncbi-blast-2.12.0+/bin/blastn -query '+output+' -out '+output+'bp1 -db '+INDEX+' -outfmt "6 qseqid sseqid qcovs length mismatch gapopen qstart qend sstart send evalue bitscore pident" -evalue 1e-1 -num_threads 42 && '+ratefilter)
    print('>>> cover rate > 80 found! ['+output+'bp2]:')
    s2=os.system('wc -l '+output+'bp2 && head '+output+'bp2')
    print('>>> sorting out CR > 80 bps ... ...')
    s3=os.system("cat "+output+"bp2 |awk '{print $1"+'" "'+"$2}'|sort |uniq -c|awk '$1 >100 "+"{"+"print}' > "+output+"bp3")
    print('>>> frequent CR > 80 bps & reference INDEXs ['+output+'bp3]:')
    s4=os.system("wc -l "+output+"bp3 && head "+output+"bp3 && cat "+output+"bp3| awk '{print $2}'|sort|uniq -c > "+output+"bp4")
    print('>>> CR > 80 bps sorted! ['+output+'bp4]:')
    s5=os.system('wc -l '+output+'bp4 && head '+output+'bp4')
    s6=os.system("cp "+output+" "+OUTPUT+" && cat "+output+"bp4 |awk '{print $2}' > badbp")
    badbp=open('badbp','r')
    for bp in badbp:
        s7=os.system("sed -i '/"+bp.strip()+"/,+1d' "+OUTPUT)
    badbp.close()
    print('>>> FILTERing finished! Qualified bps:')
    os.system('wc -l '+OUTPUT+' && head '+OUTPUT)
    
def Aligner(read,index,INDEX,fa,gap,length,start,end,output,OUTPUT):
    DNA=open(read,'r')
    DNAlist=[]
    DNAstr=''
    for bp in DNA:
        DNAlist.append(bp.strip())
    for bp in range(len(DNAlist)):
        if bp!=0:
            DNAstr+=str(DNAlist[bp])
    for i in range((end-start-length)//gap):
        bpstart=start+i*gap
        bpend=bpstart+length
        s1=subp.call(['cp','blank','READS.fa'])
        bp=open(fa,'w')
        bp.write('>bp_'+str(bpstart)+'_'+str(bpend)+'\n')
        bp.write(DNAstr[bpstart:bpend])
        bp.close()
        print('Good bp:',len(bplist)/2)
        if counter>700000:
            print("======= PROGRAM CUT OFF BY COUNTER LIMIT =======")
            break
        FetchRepeat('READS.fa','dbname')
        s2=subp.call(['rm','-f','READS.fa'])
        estimated=(time.time()-t1)/(i+1)*((end-start-length)//gap)
        print('Succeeded:',i+1,'out of',(end-start-length)//gap+1,', >>> ESTIMATED TIME:',time.time()-t1,'out of',int(estimated),'s')
    out=open(output,'w')
    for result in range(len(bplist)):
        out.write(bplist[result])
        if result%2==1:
            out.write('\n')
    out.close()
    t2=time.time()
    print('Good bp numbers:',len(bplist)//2)
    print('============ STEP 1 of 2: CHOOSING COMPLETE! ============\n>>> Time consumed:',t2-t1,'s\n>>>>>>>>>>  Initial BPs saved in:',output,'\n')
    Narrow(output,INDEX,OUTPUT)

read='sequence.fasta'
index='index1'
INDEX='INDEX1'
fa='READS.fa'
gap=25000
length=60
start=122567768
end=127687768
output='TESTchr_1'
OUTPUT='OUTPUTtest_1.txt'

Aligner(read,index,INDEX,fa,gap,length,start,end,output,OUTPUT)
t3=time.time()
print('============ STEP 2 of 2: READS COMPLETE! ============\n>>> Time consumed: ',t3-t1,'s\n>>>>>>>>>>  Qualified BPs saved in:',OUTPUT)
