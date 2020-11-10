#script to replace amino acid change at perticular position in fasta file
#python MutantDatabaseCreation.py <Fasta file> <mutation file>
import sys
import re
from BasicSequenceUtilities import *

class MutantDatabaseCreation():
    def __init__(self,arg):
        refdb_ = open(arg[1])
        aachange_=open(arg[2])
        U=BasicUtilities()
        gnmutmap={}
        pdbmap={}
        bname=re.split(".txt",arg[2])[0]
        fo=open(bname+"_mutantdb.txt","w")
        fo1 = open(bname+"_unmached_cord.txt", "w")
        #self.mutdb=open()

        for i in aachange_:
            gene = i.strip().split("\t")[0]
            achg=i.strip().split("\t")[1]
            ref=re.split("[0-9]+", achg)[0]
            alt=re.split("[0-9]+", achg)[1]
            pos=re.split("[A-z]+", achg)[1]
            #print(ref,alt,pos)

            if ref != "Ter":
                if alt != "?":
                    if len(ref)>3:
                        if len(ref)==6:
                            r1=ref[0:3]
                            r2 = ref[3:]
                            a1=alt[0:3]
                            a2=alt[3:]
                            if r1==a1:
                                gnmutmap.setdefault(gene, []).append(U.aaTriLToSlConvert(r2.upper()) + "\t" + str(int(pos)+1) + "\t" +U.aaTriLToSlConvert(a2.upper()))
                            elif r1!=a1:
                                gnmutmap.setdefault(gene, []).append(U.aaTriLToSlConvert(r1.upper()) + "\t" + pos+ "\t" + U.aaTriLToSlConvert(a1.upper()))
                                gnmutmap.setdefault(gene, []).append(U.aaTriLToSlConvert(r2.upper()) + "\t" + str(int(pos)+1) + "\t" +U.aaTriLToSlConvert(a2.upper()))

                    else:

                        rfsn = U.aaTriLToSlConvert(ref.upper())
                        atsn = U.aaTriLToSlConvert(alt.upper())
                        gnmutmap.setdefault(gene,[]).append(rfsn+"\t"+pos+"\t"+atsn)


        for rc in refdb_:
            header=rc.strip().split("\t")[0]
            seq=rc.strip().split("\t")[1]
            if "GN" in header:
                dbgene=header.split("GN=")[1].split("PE=")[0]
                pdbmap[dbgene.strip()]=header+"\t"+seq


        for m,p in gnmutmap.items():
            if m in pdbmap.keys():
                
                seq=pdbmap[m].split("\t")[1]
                head=pdbmap[m].split("\t")[0].replace(">sp|",">sp|m_").replace("_HUMAN",'')
                Gn=head.split("GN=")[1].split("PE=")[0].strip()
                mseq=seq
                mcord=""
                for ind in range(1,len(list(seq))):
                    for mv in p:

                        if ind==int(mv.split("\t")[1]):
                            if seq[ind-1]==mv.split("\t")[0]:
                                mcord = mcord + mv.split("\t")[0] + mv.split("\t")[1] + mv.split("\t")[2] + ","
                                mseq=mseq[:ind-1] +mv.split("\t")[2]+mseq[1 + ind-1:]
                            elif seq[ind-1]!=mv.split("\t")[0]:
                                fo1.write(m+"\t"+seq[ind-1]+'\t'+mv+"\n")
                if mseq!=seq:
                    fo.write(head.split("|")[0]+"|"+head.split("|")[1]+"_MU="+mcord[0:len(mcord)-1].replace(",","_")+"_Gn="+Gn+"|"+head.split("|")[2] +"\n"+mseq+"\n")
if __name__=="__main__":
    MutantDatabaseCreation(sys.argv)

