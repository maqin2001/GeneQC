import re
import sys,os


ref_genome = sys.argv[1]
gff_file = sys.argv[2]
output = sys.argv[3]
gff_out = sys.argv[4]

f1 = open(ref_genome,'r')
length={}
genome={}
for line in f1.readlines():
    if(re.match("^>",line)):
        array1=re.split(r'\s+',line)
        array1[0]=array1[0].replace(">","",1)
        idnumber=array1[0]
        number=0
        check=0
    else:
        line=line.rstrip('\n')
        if(check==0):
            check=line
            length[idnumber]=str(len(check))
            print(idnumber+"\t"+length[idnumber])
        genome[idnumber]=genome.get(idnumber,{})
        genome[idnumber][number]=line
        number+=1
f1.close()

use = 0
f2 = open(gff_file,'r')
f3 = open(output,"w+")
f4 = open(gff_out,"w+")
exon={}
exon_num=0
for line in f2.readlines():
    line=line.rstrip('\n')
    if not line.startswith("#"):
        array1=re.split(r'\t',line)
        exon[exon_num]=""
        if(re.match("transcript|primary_transcript|mRNA",array1[2])):
            arrayid = re.split(r'transcript_id=',array1[8])
            arraygff = re.split(r'\t',line,5)
            gfflen=int(array1[4]) - int(array1[3])       
            if(use==1):
                if(re.match("-",array1[6])):
                    i=exon_num
                    while(i>=0):
                        l=l+exon[i]
                        i-=1
                else:
                    for i in range(0,exon_num+1):
                        l=l+exon[i]
                exon = None
                exon_num=0
                exon={}
                exon[exon_num]=""
                print(l,file=f3)
                print(">"+arrayid[1],file=f3)
                l=""
                print(str(gffend)+"\t"+temp+"\n"+arrayid[1]+"\t"+arraygff[1]+"\t"+arraygff[2]+"\t"+"1\t",file=f4,end="")
                gffend=0
                temp=arraygff[5]
            else:
                print(">"+arrayid[1],file=f3)
                l=""
                gffend=0
                print(arrayid[1]+"\t"+arraygff[1]+"\t"+arraygff[2]+"\t"+"1\t",file=f4,end="")
                temp=arraygff[5]
                use=1
                exon_num=0
                exon={}
                exon[exon_num]=""
        if(re.match("exon",array1[2])):
            start=int(array1[3])
            end=int(array1[4])
            gffend=gffend+end-start+1
            ps=int((start-1)/int(length[array1[0]]))
            ps2=(start-1) % int(length[array1[0]])
            pe=int((end/int(length[array1[0]])))
            pe2=end % int(length[array1[0]])
            if(pe2==0):
                pe+=1
            arraystart=genome[array1[0]][ps]
            arrayend=genome[array1[0]][pe]
            if(ps==pe):
                for i in range(ps2,pe2):
                    exon[exon_num]=exon[exon_num]+arraystart[i]
            else:
                for i in range(ps2,int(length[array1[0]])):
                    exon[exon_num]=str(exon[exon_num])+arraystart[i]
                for i in range(ps+1,pe):
                    exon[exon_num]=exon[exon_num]+genome[array1[0]][i]
                for i in range(0,pe2):
                    exon[exon_num]=exon[exon_num]+arrayend[i]
        exon_num+=1
f2.close()
f3.close()
f4.close()

os.mkdir("hisatindex")
