import sys,os
import shutil
import re
import math
import numpy as np
import pandas as pd
from sklearn.linear_model import ElasticNet

# find the path
Script_dir = os.path.dirname(os.path.abspath(sys.argv[0]))
# read files
option = sys.argv[1]
# GeneQC for Plant
if(option == "1"):
    print("GeneQC for Plant")
    ref = sys.argv[2]
    gff = sys.argv[3]
    mapping = sys.argv[4]

    # load package
    samtools = Script_dir + "/samtools-1.2/samtools"

    # read name of bam or 
    indexname=mapping[0:-4]
    bam_un=indexname+".bam"
    bam_sort=indexname+"_sorted"
    bam=bam_sort+".bam"

    # covert to bam for sam file
    if mapping.endswith('.sam'):
        print("Convert Sam to Bam")
        os.system(samtools+" view -bS " +mapping+" > "+bam_un)

    # use samtools building index
    print("Building index")
    os.system(samtools+" "+"sort"+" "+bam_un+" "+bam_sort)
    os.system(samtools+" index "+bam)

    # create output folder
    out=Script_dir+"/"+indexname+"_out"
    if os.path.exists(out):
        shutil.rmtree(out)
    os.mkdir(out)
    geneoutput=out+"/temp/"
    os.mkdir(geneoutput)

    # Step1: extract the genes' sequences for self blast
    f1 = open(ref,'r')
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

    gene=geneoutput+"gene.fastq"
    f2 = open(gff,"r")
    f3 = open(gene,"w+")
    for line in f2.readlines():
        array1=re.split(r'\t',line)
        if(len(array1)>1):
            if(re.match("gene",array1[2])):
                start=array1[3]
                end=array1[4]
                ps= int(int(start)/ int(length[array1[0]]))
                ps2=int(start) % int(length[array1[0]])
                pe=int(int(end)/ int(length[array1[0]]))
                pe2=int(end) % int(length[array1[0]])
                arrayid=re.split(";",array1[8])
                arrayid2=re.split("=",arrayid[0])
                print(">"+arrayid2[1],file=f3)
                arraystart=genome[array1[0]][ps]
                arrayend=genome[array1[0]][pe]
                if(ps==pe):
                    for i in range(ps2,pe2):
                        print(arraystart[i],file=f3,end="")
                else:
                    for i in range(ps2,int(length[array1[0]])):
                        print(arraystart[i],file=f3,end="")
                    for i in range(ps+1,pe):
                        print(genome[array1[0]][i],file=f3,end="")
#                   for i in range(0,pe2):
#                       print(arrayend[i],file=f3,end="")
#                       i+=1
                print("\n",file=f3,end="")
    f2.close()
    f3.close()

    blast = Script_dir + "/blast/bin/"
    makeblastdb = blast+"makeblastdb"

    os.system(makeblastdb+" -in "+gene+" -out "+geneoutput+"self"+" -dbtype nucl")
    blastn=blast+"blastn"
    os.mkdir(geneoutput+"blast")
    os.system(blastn+" -query "+gene+" -db "+geneoutput+"self -out "+ geneoutput+"blast/self -outfmt 6 -num_threads 12 -evalue 0.000000001")
    # Step1 Self blast done

    # step2: parse blast results and find groupA gene output
    # create tier1, tierA and group
    tier1={}
    group={}
    tierA_parter={}
    file = open(geneoutput+"blast/self",'r')
    line = file.readline()
    while line:
        array1=re.split(r'\t+',line)
        if((array1[0]!=array1[1]) and (float(array1[2])*float(array1[3])>10000) and (int(array1[4])<=5) and (int(array1[5])<=5)):
            if(array1[0] not in tier1):
                tier1[array1[0]]=array1[2]+"\t"+array1[3]
            if not (array1[1] in group.get(array1[0],{})):
	            group[array1[0]]=group.get(array1[0],{})
	            group[array1[0]][array1[1]]=array1[2]+"\t"+array1[3]
        else:
            if not (re.match(array1[1],array1[0])) and (array1[0] not in tierA_parter):
                tierA_parter[array1[0]]=array1[1]+"\t"+array1[2]+"\t"+array1[3]
        line=file.readline()
    file.close()

    # get genelist
    genelist={}
    file = open(gff,"r")
    line = file.readline()
    while line:
        line=line.rstrip('\n')
        array1=re.split(r'\t+',line)
        if(len(array1)>1):
            if(re.match("gene",array1[2])):
                arrayid=re.split(";",array1[8])
                arrayid2=re.split("=",arrayid[0])
                genelist[arrayid2[1]]=1
        line = file.readline()
    file.close()

    # write tierA
    file1= open(out+ "/tierA.txt","w+")
    for Key in genelist:
        if (Key not in tier1):
            if(Key in tierA_parter):
                print(Key+"\t"+"A\t"+tierA_parter[Key],file=file1)
            else:
                print(Key+"\t"+"A\tUNIQUE\t0",file=file1)
    file1.close()

    # write tier1
    file2= open(out+ "/tier1.txt","w+")
    for Key in tier1:
        print(Key+"\t"+tier1[Key],file=file2)
    file2.close()

    # write group
    file3 = open(out+ "/group.txt","w+")
    for Key in genelist:
        for key2 in group.get(Key,{}):
            print(Key+"\t"+key2+"\t"+group[Key][key2],file=file3)
    file3.close()
    # step2 done

    # Step3: extract short reads
    # samtools extract
    samout = geneoutput+"samtools/"
    os.mkdir(samout)
    gene={}
    tier1 = out+"/tier1.txt"
    file=open(tier1,"r")
    line=file.readline()
    while line:
        line=line.rstrip('\n')
        array1=re.split(r'\s+',line)
        gene[array1[0]]=1
        line=file.readline()
    file.close()

    file=open(gff,"r")
    line=file.readline()
    while line:
        line=line.rstrip('\n')
        array1=re.split(r'\t+',line)
        if(len(array1)>1):
            if(re.match("gene",array1[2])):
                start=array1[3]
                end=array1[4]
                arrayid=re.split(";",array1[8])
                arrayid2=re.split("=",arrayid[0])
                if(arrayid2[1] in gene):
                    loc=array1[0]+":"+array1[3]+"-"+array1[4]
                    os.system(samtools+" view "+bam+" "+loc+" >"+samout+arrayid2[1])
        line=file.readline()
    file.close()

    # b-c-d
    file1=open(out+"/bcd.txt","w+")
    num1=num2=num3=0
    array1=[]
    gene={}
    list=os.listdir(samout)
    for filename in list:
        files=open(samout+filename,"r")
        line=files.readline()
        while line:
            array1=re.split(r'\s+',line)
            gene[filename]=gene.get(filename,{})
            gene[filename][array1[0]]=1
            line=files.readline()
        files.close()

    group=open(out+ "/group.txt","r")
    line=group.readline()
    while line:
        array1=re.split(r'\s+',line)
        mul={}
        for key in gene.get(array1[0],{}):
            if key in gene.get(array1[1],{}):
                mul[key]=1
        num1 = len(gene.get(array1[0],{}))           
        num2 = len(gene.get(array1[1],{}))
        num3 = len(mul)      
        if(num1>0):
            num4=num3/num1
        else:
            num4=0
        print(array1[0]+"\t"+array1[1]+"\t",num1,"\t",num2,"\t",num3,"\t",num4,file=file1)
        num3=0
        line=group.readline()
    group.close()
    file1.close()

    # unique count
    num=0
    gene={}
    group=open(out+ "/group.txt","r")
    record="start gene"
    file1=open(out+"/unique_reads.txt","w+")
    line=group.readline()
    while line:
        array1=re.split(r'\s+',line)
        if(re.match(record,array1[0])):
            inp = samout+array1[1]
            file=open(inp,"r")
            line2=file.readline()
            while line2:
                array2=re.split(r'\s+',line2)
                if array2[0] in gene:
                    del gene[array2[0]]
                line2=file.readline()
            file.close()
        else:
            num=len(gene)
            print(record+"\t",num,file=file1)
            gene={}
            record=array1[0]
            inp = samout+array1[0]
            file=open(inp,"r")
            line2=file.readline()
            while line2:
                array2=re.split(r'\s+',line2)
                gene[array2[0]]=1
                line2=file.readline()
            file.close()
            inp = samout+array1[1]
            file=open(inp,"r")
            line2=file.readline()
            while line2:
                array2=re.split(r'\s+',line2)
                if array2[0] in gene:
                    del gene[array2[0]]
                line2=file.readline()
            file.close
        line=group.readline()
    print(record+"\t",num,file=file1)
    group.close()
    file1.close()

    # score
    group=open(out+"/group.txt","r")
    bcd=open(out+ "/bcd.txt","r")
    tiera=open(out+ "/tierA.txt","r")
    uniquereads=open(out+ "/unique_reads.txt","r")
    gene={}
    sim={}
    gene2={}
    line1=group.readline()
    while line1:
        line1=line1.rstrip('\n')
        array1=re.split(r'\s+',line1)
        if array1[0] not in gene:
            sim[array1[0]]=1
            gene[array1[0]]=array1[2]
            gene2[array1[0]]=array1[3]
        else:
            sim[array1[0]]=float(sim[array1[0]])+1
            if(gene[array1[0]]<array1[2]):
                gene[array1[0]]=array1[2]
                gene2[array1[0]]=array1[3]
        line1=group.readline()
    group.close()

    mmr={}
    per={}
    total={}
    line1=bcd.readline()
    while line1:
        line1=line1.rstrip('\n')
        array1=re.split(r'\s+',line1)    
        if array1[0] not in mmr:
            per[array1[0]]=0
            mmr[array1[0]]=array1[5]
            total[array1[0]]=array1[2]
            if(float(array1[5])>0):
                per[array1[0]]+=1
        else:
            if(float(array1[5])>0):
                per[array1[0]]+=1
            if(mmr[array1[0]]<array1[5]):
                mmr[array1[0]]=array1[5]            
        line1=bcd.readline()
    bcd.close()

    mmrp={}
    line1=uniquereads.readline()
    while line1:
        if(re.match("start",line1)):
            line1=uniquereads.readline()
        line1=line1.rstrip('\n')
        array1=re.split(r'\s+',line1)
        if(float(total[array1[0]])>0):
            mmrp[array1[0]]=(float(total[array1[0]])-float(array1[1]))/float(total[array1[0]])
        else:
            mmrp[array1[0]]=0
        line1=uniquereads.readline()
    uniquereads.close()    
    
    mf={}
    line1=tiera.readline()
    while line1:    
        line1=line1.rstrip('\n')
        array1=re.split(r'\s+',line1)    
        mf[array1[0]]=0
        if(re.match("UNIQUE",array1[2])):
            gene[array1[0]]=0
        else:
            gene[array1[0]]=array1[3]
        sim[array1[0]]=per[array1[0]]=0
        mmr[array1[0]]=mmrp[array1[0]]=0
        line1=tiera.readline()
    tiera.close()

    temp=open(out+"/temp.txt","w+")
    min=1
    max=0
    for key in gene:
        gene[key]=float(gene[key])/100
        num=math.log(sim[key]+per[key]+1)/math.log(10)
        d=num*(float(gene[key])+3*float(mmr[key]))
        if(d<min):
            min=d
        if(d>max):
            max=d
        if key in mf:
            print(key+"\t",mf[key],"\t",mmr[key],"\t",mmrp[key],"\t",num,"\t",d,"\t0",file=temp)
        else:
            p=gene[key]*float(gene2[key])
            print(key+"\t",gene[key],"\t",mmr[key],"\t",mmrp[key],"\t",num,"\t",d,"\t",p,file=temp)
        num=0
    temp.close()
    temp=open(out+"/temp.txt","r")
    results=open(Script_dir+"/"+indexname+"_out.txt","w+")
    print("Gene_ID\tMAX_Similarity_persentage:D1\tMAX_MMR_persentage:D2\tDegree_weight:D3\tUnion_MMR_persentage",file=results)
    line1=temp.readline()
    while line1:
        line1=line1.rstrip('\n')
        array1=re.split(r'\s+',line1)
        print(array1[0]+"\t",array1[6],"\t",array1[2],"\t",array1[4],"\t",array1[3],file=results)
        line1=temp.readline()
    results.close()
    temp.close()
    shutil.rmtree(out)

    # Modeling parts
    data = pd.read_csv(Script_dir+"/"+indexname+"_out.txt",delim_whitespace=True, skiprows=1,  index_col=0, names=['D1','D2','D3','D4'],float_precision='round_trip')
    data.D4 = (data.D4) /2 #modify D3
    data=data[data.D2 !=0 & pd.notnull(data.D4)]

    X=np.array(data.iloc[:,0:3])
    y=np.array(data.D4)

    # choose the best alpha
    scores = []
    alpha_range = np.arange(0.1,1,0.1)
    for a in alpha_range:
        regr = ElasticNet(alpha=a)
        regr.fit(X,y)
        score = regr.score(X,y)
        scores.append(score)
    bestalpha = alpha_range[np.argmax(scores)]
    regr = ElasticNet(alpha=bestalpha)
    regr.fit(X,y)
    Dscore = regr.predict(X)
    Dscore = Dscore.reshape(-1,1)
    data['Dscore']=Dscore
    pd.DataFrame.to_csv(data,Script_dir+"/"+indexname+"_out.csv")

    os.remove(Script_dir+"/"+bam)
    os.remove(Script_dir+"/"+bam+".bai")

    print("Job Done!")

# GeneQC for animal
elif(option=="2"):
    print("GeneQC for Animal")
    ref = sys.argv[2]
    gff = sys.argv[3]
    mapping = sys.argv[4]
    
    # load package
    samtools = Script_dir + "/samtools-1.2/samtools"
    # read name of bam or 
    indexname=mapping[0:-4]
    bam_un=indexname+".bam"
    bam_sort=indexname+"_sorted"
    bam=bam_sort+".bam"

    # covert to bam for sam file
    if mapping.endswith('.sam'):
        print("Convert Sam to Bam")
        os.system(samtools+" view -bS " + mapping+" > "+bam_un)
    # use samtools building index
    print("Building index")
    os.system(samtools+" "+"sort"+" "+bam_un+" "+bam_sort)
    os.system(samtools+" index "+bam)
    # create output folder
    out=Script_dir+"/"+indexname+"_out"
    if os.path.exists(out):
        shutil.rmtree(out)
    os.mkdir(out)
    geneoutput=out+"/temp/"
    os.mkdir(geneoutput)
    blast = Script_dir + "/blast/bin/"
    makeblastdb = blast+"makeblastdb"
    query=geneoutput+"gene_fastq"
    shutil.copyfile(ref,query)
    os.system(makeblastdb+" -in "+query+" -out "+geneoutput+"self"+" -dbtype nucl")
    blastn=blast+"blastn"
    os.mkdir(geneoutput+"blast")
    os.system(blastn+" -query "+query+" -db "+geneoutput+"self -out "+geneoutput+"blast/self -outfmt 6 -num_threads 12 -evalue 0.000000001")
    # Step1 Self blast done
    
    # step2: parse blast results and find groupA gene output
    # create tier1, tierA and group
    tier1={}
    group={}
    tierA_parter={}
    file = open(geneoutput+"blast/self",'r')
    line = file.readline()
    while line:
        array1=re.split(r'\t+',line)
        if((array1[0]!=array1[1]) and (float(array1[2])*float(array1[3])>10000) and (int(array1[4])<=5) and (int(array1[5])<=5)):
            if(array1[0] not in tier1):
                tier1[array1[0]]=array1[2]+"\t"+array1[3]
            if not (array1[1] in group.get(array1[0],{})):
	            group[array1[0]]=group.get(array1[0],{})
	            group[array1[0]][array1[1]]=array1[2]+"\t"+array1[3]
        else:
            if not (re.match(array1[1],array1[0])) and (array1[0] not in tierA_parter):
                tierA_parter[array1[0]]=array1[1]+"\t"+array1[2]+"\t"+array1[3]
        line=file.readline()
    file.close()
    
    # get genelist
    genelist={}
    file = open(gff,"r")
    line = file.readline()
    while line:
        line=line.rstrip('\n')
        array1=re.split(r'\s+',line)
        genelist[array1[0]]=1
        line = file.readline()
    file.close()
    
    # write tierA
    file1= open(out+ "/tierA.txt","w+")
    for Key in genelist:
        if (Key not in tier1):
            if(Key in tierA_parter):
                print(Key+"\t"+"A\t"+tierA_parter[Key],file=file1)
            else:
                print(Key+"\t"+"A\tUNIQUE\t0",file=file1)
    file1.close()

    # write tier1
    file2= open(out+"/tier1.txt","w+")
    for Key in tier1:
        print(Key+"\t"+tier1[Key],file=file2)
    file2.close()

    # write group
    file3 = open(out+ "/group.txt","w+")
    for Key in group:
        for key2 in group.get(Key,{}):
            print(Key+"\t"+key2+"\t"+group[Key][key2],file=file3)
    file3.close()
    # step2 done

    # Step3: extract short reads
    # samtools extract
    samout = geneoutput+"samtools/"
    os.mkdir(samout)
    gene={}
    tier1 = out+"/tier1.txt"
    file=open(tier1,"r")
    line=file.readline()
    while line:
        line=line.rstrip('\n')
        array1=re.split(r'\s+',line)
        gene[array1[0]]=1
        line=file.readline()
    file.close()

    file=open(gff,"r")
    line=file.readline()
    while line:
        array1=re.split(r'\s+',line)
        if "transcript" in array1[2]:
            start=array1[3]
            end=array1[4]
            if(array1[0] in gene):
                loc=array1[0]+":"+array1[3]+"-"+array1[4]
                os.system(samtools+" view "+bam+" "+loc+" >"+samout+array1[0])
        line=file.readline()
    file.close()
    
    # b-c-d
    file1=open(out+"/bcd.txt","w+")
    num1=num2=num3=0
    gene={}
    list=os.listdir(samout)
    for filename in list:
        files=open(samout+filename,"r")
        line=files.readline()
        while line:
            array1=re.split(r'\s+',line)
            gene[filename]=gene.get(filename,{})
            gene[filename][array1[0]]=1
            line=files.readline()
        files.close()

    group=open(out+ "/group.txt","r")
    line=group.readline()
    while line:
        array1=re.split(r'\s+',line)
        mul={}
        for key in gene.get(array1[0],{}):
            if key in gene.get(array1[1],{}):
                mul[key]=1
        num1 = len(gene.get(array1[0],{}))           
        num2 = len(gene.get(array1[1],{}))
        num3 = len(mul)      
        if(num1>0):
            num4=num3/num1
        else:
            num4=0
        print(array1[0]+"\t"+array1[1]+"\t",num1,"\t",num2,"\t",num3,"\t",num4,file=file1)
        num3=0
        line=group.readline()
    group.close()
    file1.close()
    
    # unique count
    num=0
    gene={}
    group=open(out+ "/group.txt","r")
    record="start gene"
    file1=open(out+"/unique_reads.txt","w+")
    line=group.readline()
    while line:
        array1=re.split(r'\s+',line)
        if(re.match(record,array1[0])):
            inp = samout+array1[1]
            if os.path.exists(inp):
                file=open(inp,"r")
                line2=file.readline()
                while line2:
                    array2=re.split(r'\s+',line2)
                    if array2[0] in gene:
                        del gene[array2[0]]
                    line2=file.readline()
                file.close()
        else:
            num=len(gene)
            print(record+"\t",num,file=file1)
            gene={}
            record=array1[0]
            inp = samout+array1[0]
            if os.path.exists(inp):
                file=open(inp,"r")
                line2=file.readline()
                while line2:
                    array2=re.split(r'\s+',line2)
                    gene[array2[0]]=1
                    line2=file.readline()
                file.close()
            inp = samout+array1[1]
            if os.path.exists(inp):
                file=open(inp,"r")
                line2=file.readline()
                while line2:
                    array2=re.split(r'\s+',line2)
                    if array2[0] in gene:
                        del gene[array2[0]]
                    line2=file.readline()
                file.close
        line=group.readline()
    print(record+"\t",num,file=file1)
    group.close()
    file1.close()
    
    # score
    group=open(out+"/group.txt","r")
    bcd=open(out+ "/bcd.txt","r")
    tiera=open(out+ "/tierA.txt","r")
    uniquereads=open(out+ "/unique_reads.txt","r")
    gene={}
    sim={}
    gene2={}
    line1=group.readline()
    while line1:
        line1=line1.rstrip('\n')
        array1=re.split(r'\s+',line1)
        if array1[0] not in gene:
            sim[array1[0]]=1
            gene[array1[0]]=array1[2]
            gene2[array1[0]]=array1[3]
        else:
            sim[array1[0]]=float(sim[array1[0]])+1
            if(gene[array1[0]]<array1[2]):
                gene[array1[0]]=array1[2]
                gene2[array1[0]]=array1[3]
        line1=group.readline()
    group.close()

    mmr={}
    per={}
    total={}
    line1=bcd.readline()
    while line1:
        line1=line1.rstrip('\n')
        array1=re.split(r'\s+',line1)    
        if array1[0] not in mmr:
            per[array1[0]]=0
            mmr[array1[0]]=array1[5]
            total[array1[0]]=array1[2]
            if(float(array1[5])>0):
                per[array1[0]]+=1
        else:
            if(float(array1[5])>0):
                per[array1[0]]+=1
            if(mmr[array1[0]]<array1[5]):
                mmr[array1[0]]=array1[5]            
        line1=bcd.readline()
    bcd.close()

    mmrp={}
    line1=uniquereads.readline()
    while line1:
        if(re.match("start",line1)):
            line1=uniquereads.readline()
        line1=line1.rstrip('\n')
        array1=re.split(r'\s+',line1)
        if(float(total[array1[0]])>0):
            mmrp[array1[0]]=(float(total[array1[0]])-float(array1[1]))/float(total[array1[0]])
        else:
            mmrp[array1[0]]=0
        line1=uniquereads.readline()
    uniquereads.close()    
    
    mf={}
    line1=tiera.readline()
    while line1:    
        line1=line1.rstrip('\n')
        array1=re.split(r'\s+',line1)    
        mf[array1[0]]=0
        if(re.match("UNIQUE",array1[2])):
            gene[array1[0]]=0
        else:
            gene[array1[0]]=array1[3]
        sim[array1[0]]=per[array1[0]]=0
        mmr[array1[0]]=mmrp[array1[0]]=0
        line1=tiera.readline()
    tiera.close()

    temp=open(out+"/temp.txt","w+")
    min=1
    max=0
    for key in gene:
        gene[key]=float(gene[key])/100
        num=math.log(sim[key]+per[key]+1)/math.log(10)
        d=num*(float(gene[key])+3*float(mmr[key]))
        if(d<min):
            min=d
        if(d>max):
            max=d
        if key in mf:
            print(key+"\t",mf[key],"\t",mmr[key],"\t",mmrp[key],"\t",num,"\t",d,"\t0",file=temp)
        else:
            p=gene[key]*float(gene2[key])
            print(key+"\t",gene[key],"\t",mmr[key],"\t",mmrp[key],"\t",num,"\t",d,"\t",p,file=temp)
        num=0
    temp.close()
    temp=open(out+"/temp.txt","r")
    results=open(Script_dir+"/"+indexname+"_out.txt","w+")
    print("Gene_ID\tMAX_Similarity_persentage:D1\tMAX_MMR_persentage:D2\tDegree_weight:D3\tUnion_MMR_persentage",file=results)
    line1=temp.readline()
    while line1:
        line1=line1.rstrip('\n')
        array1=re.split(r'\t',line1)
        print(array1[0]+"\t",array1[6],"\t",array1[2],"\t",array1[4],"\t",array1[3],file=results)
        line1=temp.readline()
    results.close()
    temp.close()

    shutil.rmtree(out)
    
    # Modeling parts
    data = pd.read_csv(Script_dir+"/"+indexname+"_out.txt",delim_whitespace=True, skiprows=1,  index_col=0, names=['D1','D2','D3','D4'],float_precision='round_trip')
    data.D4 = (data.D4) /2 #modify D3
    data=data[data.D2 !=0 & pd.notnull(data.D4)]

    X=np.array(data.iloc[:,0:3])
    y=np.array(data.D4)

    # choose the best alpha
    scores = []
    alpha_range = np.arange(0.1,1,0.1)
    for a in alpha_range:
        regr = ElasticNet(alpha=a)
        regr.fit(X,y)
        score = regr.score(X,y)
        scores.append(score)
    bestalpha = alpha_range[np.argmax(scores)]

    regr = ElasticNet(alpha=bestalpha)
    regr.fit(X,y)
    Dscore = regr.predict(X)
    Dscore = Dscore.reshape(-1,1)
    data['Dscore']=Dscore
    pd.DataFrame.to_csv(data,Script_dir+"/"+indexname+"_out.csv")

    os.remove(Script_dir+"/"+bam)
    os.remove(Script_dir+"/"+bam+".bai")

    print("Job Done!")
    
else:
    print("Invalid Inputs. 1 for Plant, 2 for Animal")