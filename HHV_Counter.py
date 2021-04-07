import argparse
import os
import sys
import re
from Bio import SeqIO
from dataclasses import dataclass
import pickle
import pandas as pd
import matplotlib, random, matplotlib.cm

@dataclass
class SeqComp:
    alignpos: int
    refpos: int
    ref: str
    sam: str
    dtype: str

@dataclass
class ComRef:
    atype: str
    dtype: str
    ctype: []
    ref: str
    sam: str
    insertion: []
    refpos: []
    mut: []

@dataclass
class gene:
    id : str
    fstrand : int
    exonstart : []
    exonend : []  #+1?
    exonseq : []
    aa : str

@dataclass
class sampleAA:
    non : bool = False
    fs : bool = False
    cfs: str = ''
    aa : str = ''
    error : bool = False
    fail : bool = False
    # refpos : int = 0
    # outcome : str = ''

@dataclass
class genemut:
    gene : str = ''
    refpos : int = 0
    ref : str = ''
    sam : str = ''
    change : str = '' 
    insertion : bool = False
    inserted : str = ''
    frameshift : str = ''
    splicemut : bool = False
    rstrand: bool = False


class SampleComp:
    def __init__(self,ref,samp,mdel=100000,genelocs = None):
        self.comps = {}
        self.barcols = {}
        for s in samp.getsamples():
            print(s)
            self.comps[s] = self.CompareSeqs(samp.getalignseq(ref, masked=True),samp.getalignseq(s, masked=True), mdel)
            self.comps[s].ctype = list('n'*len(self.comps[s].dtype))
        self.reflen = len(self.comps[ref].dtype)

    def CompareSeqs(self,r,s,mdel):
        coded = ComRef('','',[],'','',[],[],[])
        for base_index in range(len(r)):
            rb = r[base_index]
            sb = s[base_index]
            dtype = ''
            if rb == '-' and sb == '-':
                coded.atype += '0'
                continue
            elif rb == sb:
                coded.atype += 'n'
                coded.dtype += 'n'
                coded.ref += rb
                coded.sam += sb
                coded.insertion.append('')
            elif sb == '-':
                coded.atype += 'd'
                coded.dtype += 'd'
                coded.ref += rb
                coded.sam += sb
                coded.insertion.append('')                
            elif rb == '-':
                coded.atype += 'i'
                if len(coded.insertion) !=0:
                    coded.insertion[-1] += sb
            elif rb == 'E' or sb == 'E':
                coded.atype += 'E'
                coded.dtype += 'E'
                coded.ref += rb
                coded.sam += sb
                coded.insertion.append('')
            else:
                coded.atype += 's'
                coded.dtype += 's'
                coded.ref += rb
                coded.sam += sb
                coded.insertion.append('')
 
        coded.atype = re.sub(r"[d]{"+str(mdel)+",}",self.repld,coded.atype)
        coded.dtype = re.sub(r"[d]{"+str(mdel)+",}",self.repld,coded.dtype)
        coded.sam = re.sub(r"[-]{"+str(mdel)+",}",self.repld,coded.sam)
        ps = 1
        for x in coded.ref:
            coded.refpos.append(ps)
            coded.mut.append('')
            ps += 1
        return coded

    def getcomp(self,sample, start = None, end = None):
        if start == None:
            return self.comps[sample]
        else:
            return ComRef('',self.comps[sample].dtype[start-1:end-1],self.comps[sample].ctype[start-1:end-1],self.comps[sample].ref[start-1:end-1],self.comps[sample].sam[start-1:end-1],self.comps[sample].insertion[start-1:end-1],self.comps[sample].refpos[start-1:end-1],self.comps[sample].mut[start-1:end-1])

    def repld(self,m):
        return 'F' * len(m.group())
    
    def getsamples(self):
        return list(self.comps.keys())

    def initialise_barcols(self):
        hex_colors_only = []
        cmap=matplotlib.cm.get_cmap(matplotlib.cm.tab20c)
        for i in range(cmap.N):
            rgb = cmap(i)[:3]
            hex_colors_only.append(matplotlib.colors.rgb2hex(rgb))
        cmap=matplotlib.cm.get_cmap(matplotlib.cm.tab20b)
        for i in range(cmap.N):
            rgb = cmap(i)[:3]
            hex_colors_only.append(matplotlib.colors.rgb2hex(rgb))
        cmap=matplotlib.cm.get_cmap(matplotlib.cm.tab20)
        for i in range(cmap.N):
            rgb = cmap(i)[:3]
            hex_colors_only.append(matplotlib.colors.rgb2hex(rgb))
        cmap=matplotlib.cm.get_cmap(matplotlib.cm.Set3)
        for i in range(cmap.N):
            rgb = cmap(i)[:3]
            hex_colors_only.append(matplotlib.colors.rgb2hex(rgb))
        random.shuffle(hex_colors_only)
        for s in self.getsamples():
            self.barcols[s] = hex_colors_only.pop(0)

    def getgenes(self,w,o):
        genes = {}
        genedir = {}
        for s in self.getsamples():
            genes[s] = {}
            for m in self.comps[s].mut:
                if m != '':
                    for key in [*m]:
                        genedir[key] = m[key].rstrand
                        if key not in [*genes[s]]:
                            genes[s][key] = {}
                            genes[s][key]['rpos'] = []
                            genes[s][key]['ref'] = []
                            genes[s][key]['sam'] = []
                            # print(key,m[key].rstrand)
                        if m[key].rstrand and m[key].ref != '':
                            genes[s][key]['rpos'] = [m[key].refpos] + genes[s][key]['rpos']
                            genes[s][key]['ref'] = [m[key].ref[::-1]] + genes[s][key]['ref']
                            genes[s][key]['sam'] = [m[key].sam[::-1]] + genes[s][key]['sam']
                        elif m[key].ref != '':
                            genes[s][key]['rpos'].append(m[key].refpos)
                            genes[s][key]['ref'].append(m[key].ref)
                            genes[s][key]['sam'].append(m[key].sam)

        samps = self.getsamples()

        galignref = {}
        galign = {}
        for g in list(genes[samps[0]].keys()):
            galignref[g] = {}
            galign[g] = {}
            gstart = 200000
            gend = 0
            for pos in range(len(genes[samps[0]][g]['rpos'])):
                if genes[samps[0]][g]['rpos'][pos] < gstart:
                    gstart = genes[samps[0]][g]['rpos'][pos]
                if genes[samps[0]][g]['rpos'][pos] > gend:
                    gend = genes[samps[0]][g]['rpos'][pos]                
                maxlen = 0
                for s in samps:
                    if len(genes[s][g]['ref'][pos]) > maxlen:
                        maxlen = len(genes[s][g]['ref'][pos])
                for s in samps:
                    if s not in galignref[g].keys():
                        galignref[g][s] = ''
                    m = maxlen - len(genes[s][g]['ref'][pos])
                    galignref[g][s] += genes[s][g]['ref'][pos] + m *'-'
                    if s not in galign[g].keys():
                        galign[g][s] = ''
                    m = maxlen - len(genes[s][g]['sam'][pos])
                    galign[g][s] += genes[s][g]['sam'][pos] + m *'-'
            if not os.path.exists('./'+ref+'/genes'):
                os.makedirs('./'+ref+'/genes')

            gcount = self.tiler(w,o,gstart-10,gend+10,g)
            gpick = './'+ref+'/genes/'+g+'|'+str(w)+'|'+str(o)+'.pickle'
            outfile = open(gpick,'wb')
            with pd.ExcelWriter('./'+ref+'/genes/'+g+'|'+str(w)+'|'+str(o)+'.xlsx') as writer:  
                for samp,df in gcount.items():
                    df.to_excel(writer, sheet_name=samp)
            pickle.dump(gcount,outfile)
            outfile.close()
            
            f = open('./'+ref+'/genes/'+g+'.fasta','w')
            print('>ref',file=f)
            print(galignref[g][samps[0]],file=f)
            for s in samps:
                print('>'+s,file=f)
                print(galign[g][s],file=f)
            f.close()
        gdpick = './'+ref+'/genes/genes.direction'
        outfile = open(gdpick,'wb')
        pickle.dump(genedir,outfile)
        outfile.close()
        return genes
                        

    def tiler(self,win, ol, start = 0, stop = False, selgene=None):
        if len(self.barcols)==0:
            self.initialise_barcols()
        if not stop or stop > self.reflen:
            stop = self.reflen
        count = {}
        print(win,ol,start,stop,selgene)
        for s in self.getsamples():        
            count[s] = pd.DataFrame(columns=['refpos',
                                'var_substitutionsAll',
                                'var_substitutions',
                                'var_deletions',
                                'var_coding_substitutionsAll',
                                'var_coding_substitutions',
                                'var_coding_deletions',                                
                                'flag_err',
                                'flag_gap',
                                'var_insertions',
                                'var_coding_insertions',
                                'aa_frameshift',
                                'aa_splice-site',
                                'aa_insertions',
                                'aa_missenseAll',
                                'aa_missense',
                                'aa_nonsense',
                                'aa_gene-insertion',
                                'aa_frameshift-missense',
                                'flag_aa_gap',
                                'flag_aa_err',
                                'aa_deletions',
                                'aa_truncation',
                                'aa_nocds',
                                'aa_run-on',
                                'color'
                                ])
        end = start + win
        while start <= stop:
            # print(start,end)
            for s in self.getsamples():
                count[s] = count[s].append(self.wcounts(s,start,end,selgene),ignore_index=False)
            start = end - ol
            end = end + win - ol
        return count



    def wcounts(self,s,start, stop, selgene=None):
        if stop > self.comps[s].refpos[-1]-1:
            stop = self.comps[s].refpos[-1]-1
        rpos = int((self.comps[s].refpos[stop] + self.comps[s].refpos[start])/2)

        dseq = self.comps[s].dtype[start:stop]
        subsall = dseq.count('s')
        subs = re.sub(r'(.)\1+',r'\1',dseq).count('s')
        dels = re.sub(r'(.)\1+',r'\1',dseq).count('d')
        err = dseq.count('E')
        gap = dseq.count('F')

        coding = self.comps[s].ctype[start:stop]  #still need to check right position
        cdseq = ''.join([val for is_coding, val in zip(coding,dseq) if is_coding=='c'])
        csubsall = cdseq.count('s')
        csubs = re.sub(r'(.)\1+',r'\1',cdseq).count('s')
        cdels = re.sub(r'(.)\1+',r'\1',cdseq).count('d')      

        inslist = self.comps[s].insertion[start:stop]
        ins = len(inslist) - inslist.count('')
        cinslist = [val for is_coding, val in zip(coding,inslist) if is_coding=='c']
        cins = len(cinslist) - cinslist.count('')

        muts = {} #per gene in mut, added together, prob need to return separate too
        aafs = 0
        aasp = 0
        aain = 0
        aamissall = 0
        aamiss = 0
        aanons = 0
        gins = 0
        aafmiss = 0
        aagap = 0
        aaerr = 0
        aadel = 0
        aatrun = 0
        aanocds = 0
        aarun = 0
        for mutpos in self.comps[s].mut[start:stop]:
            if mutpos != '':
                for gene,mut in mutpos.items():
                    if selgene == None or gene == selgene:
                        if mut.frameshift != '':
                            aafs += 1
                        if mut.insertion == True:
                            aain += 1
                        if mut.splicemut == True:
                            aasp += 1
                        if gene in muts.keys():
                            muts[gene] = muts[gene] + mut.change
                        else:
                            muts[gene] = mut.change
        for geneaa in muts:
            gene = muts[geneaa]
            aamissall += gene.count('M')
            gene = re.sub(r'(.)\1+',r'\1',gene)
            aamiss += gene.count('M')
            aanons += gene.count('N')
            gins += gene.count('I') #not including nonsense insertions counted above
            aafmiss += gene.count('F')
            aagap += gene.count('G')
            aaerr += gene.count('E')
            aadel += gene.count('D')
            aatrun += gene.count('T')
            aanocds += gene.count('X')
            aarun += gene.count('R')

            # if e or g then blank window? add flag
        bcol = self.barcols[s]
        if err > 0 or gap > 0:
            bcol = 'crimson'
        df = pd.DataFrame([[rpos,
                            subsall, 
                            subs, 
                            dels,
                            csubsall, 
                            csubs, 
                            cdels,                            
                            err,
                            gap,
                            ins,
                            cins,
                            aafs,
                            aasp,
                            aain,
                            aamissall,
                            aamiss,
                            aanons,
                            gins,
                            aafmiss,
                            aagap,
                            aaerr,
                            aadel,
                            aatrun,
                            aanocds,
                            aarun,
                            bcol
                            ]], index=[rpos],
                        columns=['refpos',
                                'var_substitutionsAll',
                                'var_substitutions',
                                'var_deletions',
                                'var_coding_substitutionsAll',
                                'var_coding_substitutions',
                                'var_coding_deletions',                                
                                'flag_err',
                                'flag_gap',
                                'var_insertions',
                                'var_coding_insertions',
                                'aa_frameshift',
                                'aa_splice-site',
                                'aa_insertions',
                                'aa_missenseAll',
                                'aa_missense',
                                'aa_nonsense',
                                'aa_gene-insertion',
                                'aa_frameshift-missense',
                                'flag_aa_gap',
                                'flag_aa_err',
                                'aa_deletions',
                                'aa_truncation',
                                'aa_nocds',
                                'aa_run-on',
                                'color'
                                ])
        return df 

class Samples:
    def __init__(self, file):
        self.alignseq = {}

        with open(file) as f:
            index = ''
            for line in f:
                if line[0] == '>':
                    index = line[1:].rstrip()
                    self.alignseq[index] = ''
                else:
                    self.alignseq[index] += line.lower().rstrip() 
                   
        self.samples = self.alignseq.keys()

    def getsamples(self):
        return self.samples
    
    def getrefpos(self, sample, alignpos):
        return len(self.alignseq[sample][0:alignpos].replace('-',''))

    def getalignpos(self, sample, refpos):
        tmp = len(self.alignseq[sample][0:refpos].replace('-',''))
        x = 0
        while tmp != refpos:
            x += 1
            tmp = len(self.alignseq[sample][0:refpos + x].replace('-',''))
        return refpos + x

    def getrefseq(self,sample, start=None, end=None, coord="ref", masked=False):
        if start and coord == "ref":
            start = self.getalignpos(sample,start)-1
            end = self.getalignpos(sample,end)
        elif start and coord == 'align':
            start = start - 1
        if masked:
            return re.sub(r'[^acgt-]','E',self.alignseq[sample][start:end].replace('-',''))
        return self.alignseq[sample][start:end].replace('-','')

    def getalignseq(self,sample,start=None, end=None, coord="align", masked=False):
        if start and coord == "align":
            start = start - 1
        elif start and coord == 'ref':
            start = self.getalignpos(sample,start)-1
            end = self.getalignpos(sample,end)   
        if masked:
            return re.sub(r'[^acgt-]','E',self.alignseq[sample][start:end])                
        return self.alignseq[sample][start:end]

class GenesFromRef:
    def __init__(self, gbfile, reference=None):
        self.tab = str.maketrans("ACTGactg", "TGACtgac")
        self.genes = self.parseGenes(gbfile)

    def parseGenes (self,gbfile):
        genes = {}
        for rec in SeqIO.parse(gbfile, "genbank"):
            if rec.features:
                for feature in rec.features:
                    if feature.type == "CDS":
                        aa = feature.qualifiers['translation'][0]
                        gname = feature.qualifiers['gene'][0].replace('/','_')
                        exs = []
                        exe = []
                        exseq = []
                        for p in feature.location.parts:
                            exs.append(int(p.start))
                            exe.append(int(p.end))
                            if (feature.strand == 1):
                                exseq.append(str(rec.seq[p.start:p.end]).lower())
                            else:
                                exseq.append(str(rec.seq[p.start:p.end].reverse_complement()).lower())
                        while gname in [*genes]:
                            gname+='_r'
                        genes[gname] = gene(gname,feature.strand,exs,exe,exseq,aa)
        return genes

    def getgenenames(self):
        return self.genes.keys()

    def getgene(self,gname):
        return self.genes[gname]

def introncheck (testgene,recomp,cursam,reverse):
    if len(testgene.exonstart)>1: #checks intron exon bundaries
        for i,v in enumerate(testgene.exonstart):
            if i == 0:
                if reverse:
                    introns = recomp.getcomp(cursam,testgene.exonstart[i]-1, testgene.exonstart[i]+3)
                else:
                    introns = recomp.getcomp(cursam,testgene.exonend[i]-1, testgene.exonend[i]+3)
                # print(introns.ref, introns.sam)
                for p,r,s in zip (introns.refpos, introns.ref, introns.sam):
                    if r != s:
                        if type(recomp.comps[cursam].mut[p-1]) != dict:
                            recomp.comps[cursam].mut[p-1] = {}
                        if testgene.id not in recomp.comps[cursam].mut[p-1].keys():
                            recomp.comps[cursam].mut[p-1][testgene.id] = genemut(gene = testgene.id, refpos = p, splicemut=True, rstrand=reverse)
                        # print(p,r,s,'splice mutation')
            elif i < len(testgene.exonstart)-1:
                if reverse:
                    introns = recomp.getcomp(cursam,testgene.exonstart[i]-1, testgene.exonstart[i]+3)
                    introne = recomp.getcomp(cursam,testgene.exonend[i]-1, testgene.exonend[i]+3)
                else:
                    introne = recomp.getcomp(cursam,testgene.exonstart[i]-1, testgene.exonstart[i]+3)
                    introns = recomp.getcomp(cursam,testgene.exonend[i]-1, testgene.exonend[i]+3)
                for p,r,s in zip (introne.refpos, introne.ref, introne.sam):
                    if r != s:
                        if type(recomp.comps[cursam].mut[p-1]) != dict:
                            recomp.comps[cursam].mut[p-1] = {}
                        if testgene.id not in recomp.comps[cursam].mut[p-1].keys():
                            recomp.comps[cursam].mut[p-1][testgene.id] = genemut(gene = testgene.id, refpos = p, splicemut=True, rstrand=reverse)
                        # print(p,r,s,'splice mutation')
                for p,r,s in zip (introns.refpos, introns.ref, introns.sam):
                    if r != s:
                        if type(recomp.comps[cursam].mut[p-1]) != dict:
                            recomp.comps[cursam].mut[p-1] = {}
                        if testgene.id not in recomp.comps[cursam].mut[p-1].keys():
                            recomp.comps[cursam].mut[p-1][testgene.id] = genemut(gene = testgene.id, refpos = p, splicemut=True, rstrand=reverse)
                        # print(p,r,s,'splice mutation')
                # print(introne.ref,introne.sam)
                # print(introns.ref,introns.sam)  
            elif i == len(testgene.exonstart)-1:
                if reverse:
                    introne = recomp.getcomp(cursam,testgene.exonend[i]-1, testgene.exonend[i]+3)
                else:
                    introne = recomp.getcomp(cursam,testgene.exonstart[i]-1, testgene.exonstart[i]+3)
                for p,r,s in zip (introne.refpos, introne.ref, introne.sam):
                    if r != s:
                        if type(recomp.comps[cursam].mut[p-1]) != dict:
                            recomp.comps[cursam].mut[p-1] = {}
                        if testgene.id not in recomp.comps[cursam].mut[p-1].keys():
                            recomp.comps[cursam].mut[p-1][testgene.id] = genemut(gene = testgene.id, refpos = p, splicemut=True,rstrand=reverse)
                        # print(p,r,s,'splice mutation')
                # print(introne.ref,introne.sam)

def processSamAA(sc,ins,cfs,fs,flag,saseq,non,sout):
    while len(sc) >= 3:
        psco = ''
        sout += sc + ' '
        while sc[:3] == '---':
            psco += '-'
            sc = sc[3:]
        sco = '-'*int(sc.count('-')/3) #number of deleted codons
        tfs = sc.count('-')%3 * -1 #frameshift from deletions
        sc = sc.replace('-','') #remove deleted bases

        if sc[:3].isupper(): #check inserted codon
            ins = True
        elif any(l.isupper() for l in sc[:3]): # check partial codon insertion
            tfs += sum(1 for u in sc[:3] if u.isupper()) # count inserted bases taking deletions into account
            ins = False
        else:
            ins = False

        if tfs < 0:
            cfs = "FSD"
        elif tfs > 0:
            cfs = "FSI"
        else:
            cfs = ''

        fs += tfs
        if fs < 0:
            fs = (fs%3)-3
        if fs > 0:
            if int(fs/3)>0:
                ins = True
            fs = fs%3
        if fs != 0:
            flag = True
        else:
            flag = False


        for d in psco:
            saseq.append([sampleAA(non,flag,cfs,'-')])

        if len(sc) >= 3:
            sout += psco + ' ' + codontable[sc[:3].upper()] + ' ' + sco + ' ' + str(fs) + ' ' + str(flag) + ' ' + str(ins)

            if ins:
                if non:
                    saseq[-1].append(sampleAA(non,flag,cfs,'-'))
                else:
                    saseq[-1].append(sampleAA(non,flag,cfs,codontable[sc[:3].upper()]))
            else:
                if non:
                    saseq.append([sampleAA(non,flag,cfs,'-')])
                else:   
                    saseq.append([sampleAA(non,flag,cfs,codontable[sc[:3].upper()])])
            for d in sco:
                saseq.append([sampleAA(non,flag,cfs,'-')]) 

            if codontable[sc[:3].upper()] == '_':
                non = True

            sc = sc[3:] #remove codon and repeat until done
    return sc,ins,cfs,fs,flag,saseq,non,sout

def getaaseqs(testgene,recomp,cursam,reverse):
    # if testgene.id == 'HN2':
    #     print('pause')
    dt = ''
    gc = ''
    rc = ''
    sc = ''
    saseq = []
    raseq = []
    aapos = []
    non = False
    flag = False
    fs = 0
    cfs = ''
    error = False
    fail = False
    pre_error = False
    old_chars = "ACGTacgt"
    replace_chars = "TGCAtgca"
    tab = str.maketrans(old_chars,replace_chars)
    for s, e, seq in zip(testgene.exonstart, testgene.exonend, testgene.exonseq):
        t = recomp.getcomp(cursam, start=s+1, end=e+1)
        recomp.comps[cursam].ctype[s:e] = 'c'*(e-s) #need to check coordinates for this
        for x,y in zip(range(len(seq)),reversed(range(len(seq)))): #y for - strand
            if reverse:
                dt += t.dtype[y].translate(tab)[::-1]
                gc += seq[x]
                rc += t.ref[y].translate(tab)[::-1]
                sc += (t.sam[y]  + t.insertion[y].upper()).translate(tab)[::-1]                
            else:
                dt += t.dtype[x]
                gc += seq[x]
                rc += t.ref[x]
                sc += t.sam[x]  + t.insertion[x].upper() 

            if 'E' in rc or 'E' in sc:
                error = True
            if 'F' in rc or 'F' in sc:
                fail = True
            if (error or fail) and sc[0] != 'E' and sc[0] != 'F': #added to keep same length
                tempsc = ''
                bcount = 0
                while 'E' not in sc[:3] and 'F' not in sc[:3]:
                    pre_error = True
                    if sc[:3].isupper():
                        saseq[-1].append(sampleAA(non,flag,cfs,'Y'))
                    else:
                        saseq.append([sampleAA(non,flag,cfs,'X')]) 
                    sc = sc[3:]                    
                 
            ins = False
            sout = ''
            if sc == '---':
                saseq.append([sampleAA(non,flag,cfs,'-')])
                sc = ''
            if len(sc) - sc.count('-') >= 3 and not error and not fail:
                sc,ins,cfs,fs,flag,saseq,non,sout = processSamAA(sc,ins,cfs,fs,flag,saseq,non,sout)
            # while len(sc) - sc.count('-') >= 3 and not error and not fail:
            #     psco = ''
            #     sout += sc + ' '
            #     while sc[:3] == '---':
            #         psco += '-'
            #         sc = sc[3:]
            #     sco = '-'*int(sc.count('-')/3) #number of deleted codons
            #     tfs = sc.count('-')%3 * -1 #frameshift from deletions
            #     sc = sc.replace('-','') #remove deleted bases

            #     if sc[:3].isupper(): #check inserted codon
            #         ins = True
            #     elif any(l.isupper() for l in sc[:3]): # check partial codon insertion
            #         tfs += sum(1 for u in sc[:3] if u.isupper()) # count inserted bases taking deletions into account
            #         ins = False
            #     else:
            #         ins = False

            #     if tfs < 0:
            #         cfs = "FSD"
            #     elif tfs > 0:
            #         cfs = "FSI"
            #     else:
            #         cfs = ''

            #     fs += tfs
            #     if fs < 0:
            #         fs = (fs%3)-3
            #     if fs > 0:
            #         if int(fs/3)>0:
            #             ins = True
            #         fs = fs%3
            #     if fs != 0:
            #         flag = True
            #     else:
            #         flag = False

            #     sout += psco + ' ' + codontable[sc[:3].upper()] + ' ' + sco + ' ' + str(fs) + ' ' + str(flag) + ' ' + str(ins)

            #     for d in psco:
            #         saseq.append([sampleAA(non,flag,cfs,'-')])
            #     if ins:
            #         if non:
            #             saseq[-1].append(sampleAA(non,flag,cfs,'-'))
            #         else:
            #             saseq[-1].append(sampleAA(non,flag,cfs,codontable[sc[:3].upper()]))
            #     else:
            #         if non:
            #             saseq.append([sampleAA(non,flag,cfs,'-')])
            #         else:   
            #             saseq.append([sampleAA(non,flag,cfs,codontable[sc[:3].upper()])])
            #     for d in sco:
            #         saseq.append([sampleAA(non,flag,cfs,'-')]) 

            #     if codontable[sc[:3].upper()] == '_':
            #         non = True

            #     sc = sc[3:] #remove codon and repeat until done

            if len(rc) == 3:
                if gc != rc:
                    print("FAIL ref and gene do not match" + gc + rc)
                if (error or fail) and not pre_error:
                    # print(testgene.id,t.refpos[y],"ref",rc, 'EF')
                    saseq.append([sampleAA(error=error,fail=fail)])
             
                if reverse:
                    # if not error and not fail:
                        # print(testgene.id,t.refpos[y],"ref",rc, codontable[rc[:3].upper()], sout) #position is recorded at 3rd base of codon in case last pos in new exon
                    aapos.append(t.refpos[y])
                else:
                    # if not error and not fail:
                        # print(testgene.id,t.refpos[x],"ref",rc, codontable[rc[:3].upper()], sout)
                    aapos.append(t.refpos[x])
                try:
                    raseq.append(codontable[rc[:3].upper()])
                    if codontable[rc[:3].upper()] == '_' and len(saseq) < len(raseq) and len(sc) > 0 and not error and not fail:
                        saseq.append([sampleAA(non,flag,cfs,'-')])
                    # if codontable[rc[:3].upper()] == '_' and len(saseq) < len(raseq):
                        # saseq.append([sampleAA(non,flag,cfs,'-')])
                except:
                    raseq.append('X')
                if sout == '':
                    x=''
                gc = ''
                rc = ''
            pre_error=False
            print(str(len(raseq)),str((len(saseq))))
    if len(rc) > 0:
        raseq.append('X')     
    if len(raseq) != len(saseq):
        if len(sc) >= 3 and not error and not fail:
            sc,ins,cfs,fs,flag,saseq,non,sout = processSamAA(sc,ins,cfs,fs,flag,saseq,non,sout)
        # if (len(sc)%3 != 0 or error or fail) and len(raseq)-len(saseq) == 1:
        if (len(sc)%3 != 0 or error or fail or non):
            while len(saseq) < len(raseq):                          #or keep adding until stop...
                saseq.append([sampleAA(error=error,fail=fail, non=non)])
        while (len(saseq) > len(raseq)) and (error or fail or non):
            del saseq[-1] 
        if len(raseq) != len(saseq):
            print ("PROBLEM")
            sys.exit("Problem with amino acid lengths")
    return aapos,raseq,saseq,error,fail

def aaresult(aapos,raseq,saseq,error,fail,recomp,cursam,reverse,testgene):
    start = True
    for p,r,s in zip(aapos,raseq,saseq):
        if type(recomp.comps[cursam].mut[p-1]) != dict: 
            recomp.comps[cursam].mut[p-1] = {}
        if testgene.id not in recomp.comps[cursam].mut[p-1].keys(): 
            recomp.comps[cursam].mut[p-1][testgene.id] = genemut(gene=testgene.id,refpos = p,rstrand=reverse)
        recomp.comps[cursam].mut[p-1][testgene.id].ref = r
        insertion = False
        nocds = False

        for saa in s:
            recomp.comps[cursam].mut[p-1][testgene.id].sam += saa.aa
            if insertion == True:
                r = '-'
                recomp.comps[cursam].mut[p-1][testgene.id].ref += '-'
                recomp.comps[cursam].mut[p-1][testgene.id].insertion = True        
            if saa.cfs != '':
                # print(p,r,saa,cfs)
                recomp.comps[cursam].mut[p-1][testgene.id].frameshift = saa.cfs 
            if error:
                # print(p,r,saa,"E")
                recomp.comps[cursam].mut[p-1][testgene.id].change += 'E'
            elif fail:
                # print(p,r,saa,"G")
                recomp.comps[cursam].mut[p-1][testgene.id].change += 'G'
            elif (start and r != saa.aa) or nocds == True:
                nocds = True #no similar for E and F?
                # print(p,r,saa,"X") #no cds first aa disrupted
                recomp.comps[cursam].mut[p-1][testgene.id].change += 'X'
            elif saa.non:
                # print(p,r,saa,"T") #prior non so aa truncated
                recomp.comps[cursam].mut[p-1][testgene.id].change += 'T'
            elif r == '-' and saa.aa == '_':
                # print(p,r,saa,"N") #NONSENSE insertion
                recomp.comps[cursam].mut[p-1][testgene.id].change += 'N'
            elif r == '-':
                # print(p,r,saa,"I") #insertion aa
                recomp.comps[cursam].mut[p-1][testgene.id].change += 'I'
            elif r != '_' and saa.aa == '_':
                # print(p,r,saa,"N") #NONSENSE
                recomp.comps[cursam].mut[p-1][testgene.id].change += 'N'
            elif r == '_' and saa.aa != '_':
                # print(p,r,saa,"R") #Run on
                recomp.comps[cursam].mut[p-1][testgene.id].change += 'R'
            elif saa.aa == '-':
                # print(p,r,saa,"D") #Deletion (fs deletion separate?)
                recomp.comps[cursam].mut[p-1][testgene.id].change += 'D'
            elif r != saa.aa and not saa.fs:
                # print(p,r,saa,"M") #missense
                recomp.comps[cursam].mut[p-1][testgene.id].change += 'M'
            elif r !=  saa.aa and saa.fs:
                # print(p,r,saa,"F") #Frameshift missense
                recomp.comps[cursam].mut[p-1][testgene.id].change += 'F'
            elif r == saa.aa:
                # print(p,r,saa,".")
                recomp.comps[cursam].mut[p-1][testgene.id].change += '.'
            else:
                print('other')

            insertion = True
            start = False

def aachecks(genes,recomp):
    for cursam in recomp.comps.keys():
        print('GENES',cursam)
        for gn in genes.getgenenames():
            print(gn)
            testgene = genes.getgene(gn)
            reverse = False
            if testgene.fstrand != 1:
                reverse = True

            introncheck(testgene,recomp,cursam,reverse)
            aapos,raseq,saseq,error,fail = getaaseqs(testgene,recomp,cursam,reverse)
            aaresult(aapos,raseq,saseq,error,fail,recomp,cursam,reverse,testgene)

parser = argparse.ArgumentParser()
parser.add_argument("-i", "--infile", help="sample multiple alignment",default="HHV-6B_final_aligned.fa.fasta")
# parser.add_argument("-i", "--infile", help="sample multiple alignment",default="6Aalign_auto_no11DXy.txt")
parser.add_argument("-r", "--reference", help="reference sequence",default="Z29")
parser.add_argument("-gb", "--genbank", help="genbank file for reference",default="Z29.gb")
parser.add_argument("-w", "--winsize", help="window size",type=int,default=500)
parser.add_argument("-o", "--overlap", help="overlap between windows, must be less than winsize",type=int,default=0)
parser.add_argument("-mi", "--maxdel", help="determine cut off between deletion and gap",type=int,default=200)
args = parser.parse_args()


infile = args.infile
ref = args.reference
genefromgb = args.genbank
print(infile,ref,args.genbank)
# wsize = args.winsize
# osize = args.overlap
wsize = [100,100,300,300,500,500,1000,1000]
osize = [0,50,0,150,0,250,0,500]
maxdel = args.maxdel
dirName = ref +'_'+str(wsize)+'_'+str(osize)

tab = str.maketrans("ACTGactg", "TGACtgac")

codontable = {
'ATA':'I', 'ATC':'I', 'ATT':'I', 'ATG':'M',
'ACA':'T', 'ACC':'T', 'ACG':'T', 'ACT':'T',
'AAC':'N', 'AAT':'N', 'AAA':'K', 'AAG':'K',
'AGC':'S', 'AGT':'S', 'AGA':'R', 'AGG':'R',
'CTA':'L', 'CTC':'L', 'CTG':'L', 'CTT':'L',
'CCA':'P', 'CCC':'P', 'CCG':'P', 'CCT':'P',
'CAC':'H', 'CAT':'H', 'CAA':'Q', 'CAG':'Q',
'CGA':'R', 'CGC':'R', 'CGG':'R', 'CGT':'R',
'GTA':'V', 'GTC':'V', 'GTG':'V', 'GTT':'V',
'GCA':'A', 'GCC':'A', 'GCG':'A', 'GCT':'A',
'GAC':'D', 'GAT':'D', 'GAA':'E', 'GAG':'E',
'GGA':'G', 'GGC':'G', 'GGG':'G', 'GGT':'G',
'TCA':'S', 'TCC':'S', 'TCG':'S', 'TCT':'S',
'TTC':'F', 'TTT':'F', 'TTA':'L', 'TTG':'L',
'TAC':'Y', 'TAT':'Y', 'TAA':'_', 'TAG':'_',
'TGC':'C', 'TGT':'C', 'TGA':'_', 'TGG':'W',
}

print('reading alignment')
samp = Samples(infile)

print('loading genes')
genelist = GenesFromRef(genefromgb)

# for gn in genes.getgenenames():
#     print(genes.getgene(gn).id)
#     print(genes.getgene(gn).fstrand)

print('Performing comparison')
recomp = SampleComp(ref,samp,mdel=maxdel)

print('processing amino acids')
aachecks(genelist,recomp)

print('saving comparison')
if not os.path.exists('./'+ref):
    os.makedirs('./'+ref+'/data')
filepick = './'+ref+'/data/'+ref+'.pickle'
outfile = open(filepick,'wb')
pickle.dump(recomp,outfile)
outfile.close()

print("saved")

# ref = "GLA_25506"

print("starting counting")
# infile = open('./'+ref+'/data/'+ref+".pickle",'rb')
# test = pickle.load(infile)
recomp.barcols={}
# print("loaded file")
# infile.close()

for win,ov in zip(wsize,osize):
    genes = recomp.getgenes(win,ov)
for win,ov in zip(wsize,osize):
    output = recomp.tiler(win,ov)
    with pd.ExcelWriter('./'+ref+'/'+ref+'|'+str(win)+'|'+str(ov)+'.xlsx') as writer:  
        for samp,df in output.items():
            df.to_excel(writer, sheet_name=samp)
    filepick = './'+ref+'/'+ref+'|'+str(win)+'|'+str(ov)+'.pickle'
    outfile = open(filepick,'wb')
    pickle.dump(output,outfile)
    outfile.close()

print("output complete")
