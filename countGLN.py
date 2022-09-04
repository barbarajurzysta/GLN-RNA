#!/usr/bin/python
# -*- coding: utf-8 -*-


import sys
import math
import argparse
## Kolejne dwa tylko do testowych rysunkow
#from mpl_toolkits.mplot3d.axes3d import Axes3D
#import matplotlib.pyplot as plt

from HBonds import BondsRead, FindLoopsTails


EPS = 0.0001
g_DEFAULT=0.69 # TRESHOLDS for lasso classification
h_DEFAULT=0.55
sl_DEFAULT=1.5
R = 2 # precision for rounding an output


##########################

def CompareEq(a,b):
    if not type(a)==list and not type(a)==tuple:
        return abs(a-b)<EPS
    else:
        assert len(a)==len(b)
    for i in range(len(a)):
        if not CompareEq(a[i],b[i]): return False
    return True

#--------------------------------------------------------------------------------

def CompareGeq(a,b): return a>b-EPS
def CompareGt(a,b): return a>b+EPS

##########################

#--------------------------------------------------------------------------------

def ScalarProduct(a,b):
# return scalar product of two vectors (assumption: a,b vectors (lists) of three double numbers
    assert (type(a)==list and type(b)==list and len(a)==3 and len(b)==3)
    return (a[0]*b[0]+a[1]*b[1]+a[2]*b[2])

#--------------------------------------------------------------------------------

def VectorProduct(a,b):
# return vector product of two vectors (assumption: a,b vectors (lists) of three double numbers
    assert (type(a)==list and type(b)==list and len(a)==3 and len(b)==3)
    return [a[1]*b[2]-a[2]*b[1],a[2]*b[0]-a[0]*b[2],a[0]*b[1]-a[1]*b[0]]

#--------------------------------------------------------------------------------

def NormalizedVectorProduct(a,b):
# return vp := a x b / |a x b| (assumption: a,b vectors (lists) of three double numbers
    assert (type(a)==list and type(b)==list and len(a)==3 and len(b)==3)
    v = VectorProduct(a,b)
    d = (v[0]*v[0]+v[1]*v[1]+v[2]*v[2])**0.5
    if d==0: return 0
    return [el/d for el in v]

##########################
def FileCheck(fn):
    try:
      open(fn, "r")
      return 1
    except IOError:
      print("Error: File",fn,"does not appear to exist.\n")
      return 0

#--------------------------------------------------------------------------------

def ChainRead(filename, begin=False, end=False, limit=15):
#Read coordinates from .xyz file, cuts them from id "begin" to "end", and put them into the returned list; empty list -> error/problem -> end of the program; works with 4 and 5 columns file.  
    if (not FileCheck(filename)): return []
    comp = []
    x,y,z = -1,-1,-1
    indexes = []
    gaps = []
    started = False
    i = 0
    f = open(filename)
    for line in f.readlines():
        line = line.split()[:4]
        x2,y2,z2 = [float(x) for x in line[1:]]
        if started and CompareEq([x,y,z],[x2,y2,z2]):
            print("Error: There were two identical atoms with the index",k,"in the file. Program is done becouse of that.\n")
            return []
        dist = math.sqrt((x2-x)**2 + (y2-y)**2 + (z2-z)**2)
        if started and dist > limit:
            gaps.append(i-1)   # gap between i-1 and i
        i+=1
        k=line[0]
        indexes.append(k)
        x,y,z = x2,y2,z2
        started = True
        if (len(indexes)>begin and len(indexes)-1<=end) or (not begin and not end):
            comp.append((k,(x,y,z)))
    f.close()
    return comp, indexes, gaps

#--------------------------------------------------------------------------------

def PrepareFragments(chain,ids=[]):
    ''' Cuts parts of chain indicated in the list ids.
    For ids==[] it returns whole chain, for ids=[(a,b),(c,d)] it returns chain consisting atoms [a,a+1,...,b,c,c+1,...,d].
    '''
    if (ids==[]): return chain
    assert(type(ids)==list or type(ids)==tuple)
    ids = list(ids)
    if (type(ids[0])!=list and type(ids[0])!=tuple): ids = [ids]
    ids.sort()

    loop = []

    for piece in ids:
        loop+=chain[piece[0]:piece[1]+1]

    return loop

def PrepareLoop(chain,ids=[]):
    ''' Additionally closes structure
    '''
    loop = PrepareFragments(chain,ids)
    if len(loop)>0: loop.append(loop[0])
    return loop


#--------------------------------------------------------------------------------

def CompsOverlie(comp1,comp2):
#Return True if there are points on two comps that overlie
    for i in range(len(comp1)):
        for j in range(len(comp2)):
            if CompareEq(comp1[i][1],comp2[j][1]): return True
    return False

##########################

def linking_oneSegment(A1,A2,B1,B2):
# returns self linking number between A1A2 and B1B2
# Ken's notations

    a = [B1[i]-A1[i] for i in [0,1,2]]
    b = [B2[i]-A1[i] for i in [0,1,2]]
    c = [B2[i]-A2[i] for i in [0,1,2]]
    d = [B1[i]-A2[i] for i in [0,1,2]]

    n1 = NormalizedVectorProduct(a,b)
    n2 = NormalizedVectorProduct(b,c)
    n3 = NormalizedVectorProduct(c,d)
    n4 = NormalizedVectorProduct(d,a)
    if n1==0 or n2==0 or n3==0 or n4==0: 
        print("Problem with Vector Product...", A1, A2, B1, B2, a, b, c, d)
        return 0


    f = [B2[i]-B1[i] for i in [0,1,2]]
    h = [A2[i]-A1[i] for i in [0,1,2]]
    s = VectorProduct(f,h)

    x = ScalarProduct(s,a)
    sign = 0
    if CompareEq(x,0): x,sign = 0,0
    if x>0: sign = 1
    if x<0: sign = -1

    as1, as2, as3, as4 = ScalarProduct(n1,n2), ScalarProduct(n2,n3), ScalarProduct(n3,n4), ScalarProduct(n4,n1)
    if (CompareGeq(as1,1)): as1=1
    if (CompareGeq(as2,1)): as2=1
    if (CompareGeq(as3,1)): as3=1
    if (CompareGeq(as4,1)): as4=1
    if (CompareGeq(-1,as1)): as1=-1
    if (CompareGeq(-1,as2)): as2=-1
    if (CompareGeq(-1,as3)): as3=-1
    if (CompareGeq(-1,as4)): as4=-1

    res = sign*(math.asin(as1) + math.asin(as2) + math.asin(as3) + math.asin(as4)) / (4*math.pi)
    return res

#--------------------------------------------------------------------------------

def linking_loopTail(loop,tail):
    N1,N2 = len(loop), len(tail)
    links = [ [0 for i in range(N2-1)] for i in range(N1-1)]
    #links[i][j] = GLN ( comp1<i,i+1>, comp2<j,j+1> ) - GLN between each pair of segments: ONE SEGMENT with ONE SEGMENT

    for i in range(N1-1):
      for j in range(N2-1):
        links[i][j] = linking_oneSegment(loop[i][1], loop[i+1][1], tail[j][1], tail[j+1][1])


    vlinksLOOP1 = [0 for i in range(N2-1)]  #vlinksLOOP1[j] = GLN ( comp1, comp2<j,j+1> ) - GLN between WHOLE comp1 (as a LOOP1) and ONE SEGMENT of comp2
    for i in range(N1-1):
        for j in range(N2-1):
            vlinksLOOP1[j] += links[i][j]

    linksLOOP1 = [ [0 for i in range(N2)] for i in range(N2)]
    #linksLOOP1[i][j] = GLN ( comp1, comp2<i-j> ) - GLN between WHOLE comp1 (as a LOOP1) and PART of comp2 <i-j>

    for j1 in range(N2): 
      for j2 in range(j1+1,N2):
        linksLOOP1[j1][j2] = linksLOOP1[j1][j2-1] + vlinksLOOP1[j2-1]
    
    return linksLOOP1
    #return linksLOOP1


##########################

def give_gLassoType (min1,max1,whole, g=g_DEFAULT, sl=sl_DEFAULT, h=h_DEFAULT ):
    cl = ""
    amin1, awhole = abs(min1),abs(whole)
    if max(max1,amin1)>=sl: cl = "gLS"
    elif max(max1,amin1)<g: cl = "gL0"
    elif amin1<g<=max1 or amin1>=g>max1: cl = "gL1"
    elif awhole<h: cl = "gL2+"
    else: cl = "gL3+"
    return cl


def give_GLNstats(linksLoop):
    #linksLoop[i][j] = GLN ( loop, tail<i-j> ) - GLN between whole loop and PART of tail <i-j>
    results = {}
    N = len(linksLoop)-1
    whole = round(linksLoop[0][N],R)
    results['wh'] = whole
    
    max1, xmax1, ymax1, min1, xmin1, ymin1 = 0,0,0,0,0,0
    for j1 in range(N):
        for j2 in range(j1+1,N):
            if linksLoop[j1][j2] > max1: max1,xmax1,ymax1 = linksLoop[j1][j2],j1,j2
            if linksLoop[j1][j2] < min1: min1,xmin1,ymin1 = linksLoop[j1][j2],j1,j2

    max1, min1 = round(max1,R), round(min1,R)
    results['cl'] = give_gLassoType(min1,max1,whole)
    
    mmax = max(max1,-min1)
    if (mmax==max1):  l,xmax,ymax = 1,xmax1,ymax1
    if (mmax==-min1): l,mmax,xmax,ymax = 1,-mmax,xmin1,ymin1
    mmax = round(mmax,R)
    
    results['max']= [max1,(xmax1,ymax1)]
    results['min']= [min1,(xmin1,ymin1)]
    results['mmax']= [mmax,(xmax,ymax)]
    
    return results



##########################

def GLN(loop, tail):

    N = len(tail)
    linksLoop = linking_loopTail(loop, tail)  #linksLOOP: N*N
    
    results = {}
    whole = round(linksLoop[0][N-1],R)
    results['wh'] = whole
    
    max1, xmax1, ymax1, min1, xmin1, ymin1 = 0,0,0,0,0,0
    for j1 in range(N):
        for j2 in range(j1+1,N):
            if linksLoop[j1][j2] > max1: max1,xmax1,ymax1 = linksLoop[j1][j2],j1,j2
            if linksLoop[j1][j2] < min1: min1,xmin1,ymin1 = linksLoop[j1][j2],j1,j2

    max1, min1 = round(max1,R), round(min1,R)
    results['cl'] = give_gLassoType(min1,max1,whole)
    
    mmax = max(max1,-min1)
    if (mmax==max1):  xmax,ymax = xmax1,ymax1
    if (mmax==-min1): mmax,xmax,ymax = -mmax,xmin1,ymin1
    mmax = round(mmax,R)
    
    results['max']= [max1,(tail[xmax1][0],tail[ymax1][0])]
    results['min']= [min1,(tail[xmin1][0],tail[ymin1][0])]
    results['mmax']= [mmax,(tail[xmax][0],tail[ymax][0])]

    return results


##########################

def GLN_complexTail(loop, tails, joined=True):
    ''' loop is chain of atoms, tails is collection od chains of atoms - next pieces of tail, but not connected between each other 
    '''

    if (len(loop)==0 or len(tails)==0):
        #print("[WARNING] Empty structure given to GLN_complexTail.\n")
        return {},{}

    linkTails = []
    for tail in tails:
        linkTails.append(linking_loopTail(loop, tail))
    

    # najpierw miny i maxy etc. w każdym kawałku ogona oddzielnie
    local_minmax = {}
    if not joined: jmmax=0
    for i in range(len(tails)):
        if len(tails[i])==0: continue
        
        t, g = tails[i], linkTails[i]
        tends = (t[0][0],t[-1][0])
        local_minmax[tends] = {'wh': round(g[0][len(g)-1],R), 'min': [0,(0,0)], 'max': [0,(0,0)], 'mmax': [0,(0,0)], 'cl': ''}

        max1, xmax1, ymax1, min1, xmin1, ymin1 = 0,0,0,0,0,0
        N = len(t)
        for j1 in range(N):
            for j2 in range(j1+1,N):
                if g[j1][j2] > max1: max1,xmax1,ymax1 = g[j1][j2],j1,j2
                if g[j1][j2] < min1: min1,xmin1,ymin1 = g[j1][j2],j1,j2

        max1, min1 = round(max1,R), round(min1,R)
        local_minmax[tends]['cl'] = give_gLassoType(min1,max1,local_minmax[tends]['wh'])
    
        mmax = max(max1,-min1)
        if (mmax==max1):  xmax,ymax = xmax1,ymax1
        if (mmax==-min1): mmax,xmax,ymax = -mmax,xmin1,ymin1
        mmax = round(mmax,R)

        local_minmax[tends]['max'] = [max1,(t[xmax1][0],t[ymax1][0])]
        local_minmax[tends]['min'] = [min1,(t[xmin1][0],t[ymin1][0])]
        local_minmax[tends]['mmax']= [mmax,(t[xmax][0], t[ymax][0])]
        if not joined:
            if abs(mmax)>=abs(jmmax):
                jmmax=mmax
                total_minmax=local_minmax[tends]
                total_minmax['min'].append([total_minmax['min'][1]])
                total_minmax['max'].append([total_minmax['max'][1]])


    def tail_format(gln, x,y):
        tstart,start = x
        tend,end = y
        if (tstart==tend):
            skladanka = [(tails[tstart][start][0],tails[tend][end][0])]
        else:
            skladanka = [(tails[tstart][start][0],tails[tstart][-1][0])]
            skladanka += [(tails[i][0][0],tails[i][-1][0]) for i in range(tstart+1,tend)]
            skladanka += [(tails[tend][0][0], tails[tend][end][0])]
        return [gln, (tails[tstart][start][0], tails[tend][end][0]), skladanka]

    # teraz dla calosci
    if joined:
        whole = 0
        for g in linkTails: whole += g[0][len(g)-1]
        max1, xmax1, ymax1, min1, xmin1, ymin1 = 0,(0,0),(0,0),0,(0,0),(0,0) # xmax1 =(tail_start, start), indeks ogona, i który na nim atom
        for tail_start in range(len(tails)):
            t = tails[tail_start]
            for start in range(len(tails[tail_start])):
            #sprawdzamy wszystkie kawałki zaczynające się od t[start], czyli od kawałka ogona t, od atomu start
                gln = 0
                tail_end, end = tail_start, start+1
                while tail_end < len(tails):
                    while end < len(tails[tail_end]):
                        gln += linkTails[tail_end][end-1][end]
                        if gln > max1: max1, xmax1, ymax1 = gln, (tail_start,start), (tail_end,end)
                        if gln < min1: min1, xmin1, ymin1 = gln, (tail_start,start), (tail_end,end)
                        end += 1
                    end = 1
                    tail_end += 1

        max1, min1, whole = round(max1,R), round(min1,R), round(whole,R)
        total_minmax = {'cl': give_gLassoType(min1,max1,whole), 'wh': whole}
        total_minmax['max'] = tail_format(max1,xmax1,ymax1)
        total_minmax['min'] = tail_format(min1,xmin1,ymin1)

        mmax = max(max1,-min1)
        if (mmax==max1):  xmax,ymax = xmax1,ymax1
        if (mmax==-min1): mmax,xmax,ymax = -mmax,xmin1,ymin1
        mmax = round(mmax,R)
        total_minmax['mmax'] = tail_format(mmax,xmax,ymax)

    return total_minmax, local_minmax





##########################

def colorFromGLN(gln):
    if (gln<-1): return (int(255*1/(gln*gln)), 0, 0) 
    elif (gln<=0): return (255, int(255*(1+gln)), int(255*(1+gln)))
    elif (gln<=1): return (int(255*(1-gln)), int(255*(1-gln)), 255)
    else: return (0, 0, int(255*1/(gln*gln)))

##########################

##########################
##########################
##########################
def main():
# Reading arguments
    if (len(sys.argv) == 1):
        print("*\nUsage of the program: python "+sys.argv[0]+" <filename1> (<begin of a comp1> <end of a comp1>) <filename2> (<begin of a comp2> <end of a comp2>) (<-additional_option and argument>^n)")
        print("*\nAdditional options:\n   -> -close (00,10,01,11): if we close 1st/2nd comp; implicitly close=00;\n   -> -loop (0,1,2): if we calculate bots \"glns\", 0:both, 1:only 1st comp as a loop, 2:only 2nd comp as a loop; implicitly loop=0;\n   -> -out (0,1,2): type of output: 0: for each comp as a loop, 1: in total (for KnotGenom), 2: in total in csv; implicitly out=2.\n*")
        return 0

    b1,e1,b2,e2 = False,False,False,False
    ar=1
    filename1 = sys.argv[ar]
    ar=ar+1
    if (sys.argv[ar]).isdigit():
        b1,e1 = int(sys.argv[ar]), int(sys.argv[ar+1])
        ar=ar+2
    filename2 = sys.argv[ar]
    ar=ar+1
    if len(sys.argv)>ar and (sys.argv[ar]).isdigit():
        b2,e2 = int(sys.argv[ar]), int(sys.argv[ar+1])
        ar=ar+2

# Read optional arguments
    cl1, cl2 = False, False # if we threat 1/2 comp as closed
    loop = 0        # if we check both combinations - 0-both, 1-comp1 as a loop, 2-comp2 as a loop 
    out = 2

    while ar<len(sys.argv)-1: # -closed, -loop, -out
      if sys.argv[ar]=="-close":
        if len(sys.argv[ar+1])!=2:
          print("Error: argument after -close must consists of two digits (0-1).")
          return 0
        if sys.argv[ar+1][0]=="1": cl1 = True
        if sys.argv[ar+1][1]=="1": cl2 = True
        ar+=2
      elif sys.argv[ar]=="-loop":
        loop = int(sys.argv[ar+1])
        if loop<0 or loop>2: loop = 0
        ar+=2
      elif sys.argv[ar]=="-out":
        out = int(sys.argv[ar+1])
        if out<0 or out>2: out = 2
        ar+=2

# Creating components   
    comp1 = ChainRead(filename1,b1,e1)
    if cl1: comp1.append(comp1[0])
    comp2 = ChainRead(filename2,b2,e2)
    if cl2: comp2.append(comp2[0])
    N1, N2 = len(comp1), len(comp2)

    if N1==0 or N2==0:
        print("empty-component")
        return 0
    if N1==1 or N2==1:
        print("1-length-comp")
        return 0
    if CompsOverlie(comp1,comp2):
        print("components-overlie")
        return 0

#    print(GLN(comp1, comp2))
#    print("KONIEC")
#    return
                                                           
# Calculate linking 
##########################
    linksLOOP1,linksLOOP2 = linking_loopTail(comp1,comp2)  #linksLOOP1: N2*N2, linksLOOP2: N1*N1
##########################

# Preparing data for out
#    out =1: wh1 max1 min1 cl1
#    out =2: spacja N2 wh1: wh1+: max1: min1: max10: min10: cl 

    whole = round(linksLOOP2[0][N1-1],R)
    max1, xmax1, ymax1, min1, xmin1, ymin1, max2, xmax2, ymax2, min2, xmin2, ymin2 = 0,0,0,0,0,0,0,0,0,0,0,0

    if loop==0 or loop==2:
      for i1 in range(N1):
       for i2 in range(i1+1,N1):
        if linksLOOP2[i1][i2] > max2: max2,xmax2,ymax2 = linksLOOP2[i1][i2],i1,i2
        if linksLOOP2[i1][i2] < min2: min2,xmin2,ymin2 = linksLOOP2[i1][i2],i1,i2
    if loop==0 or loop==1:
      for j1 in range(N2):
       for j2 in range(j1+1,N2):
        if linksLOOP1[j1][j2] > max1: max1,xmax1,ymax1 = linksLOOP1[j1][j2],j1,j2
        if linksLOOP1[j1][j2] < min1: min1,xmin1,ymin1 = linksLOOP1[j1][j2],j1,j2
 
    cl1 = give_gLassoType(min1,max1,whole)
    cl2 = give_gLassoType(min2,max2,whole)

    absmax = max(max1,-min1,max2,-min2)
    if (absmax==max1):  l,xmax,ymax = 1,xmax1,ymax1
    if (absmax==max2):  l,xmax,ymax = 2,xmax2,ymax2
    if (absmax==-min1): l,absmax,xmax,ymax = 1,-absmax,xmin1,ymin1
    if (absmax==-min2): l,absmax,xmax,ymax = 2,-absmax,xmin2,ymin2
    absmax = round(absmax,R)

    if l==1: xmax_id, ymax_id = comp2[xmax][0], comp2[ymax][0] 
    if l==2: xmax_id, ymax_id = comp1[xmax][0], comp1[ymax][0] 

    if out==0:
       if loop==0 or loop==1: print("*LOOP1 wh: "+str(whole)+" max: "+str(max1)+" min: "+str(min1)+" cl: "+cl1) 
       if loop==0 or loop==2: print("*LOOP2 wh: "+str(whole)+" max: "+str(max2)+" min: "+str(min2)+" cl: "+cl2)
    if out==1:
       print("wh: "+str(whole)+" max|GLN|: "+str(absmax)+" compAsLoop: "+str(l)+" fragmentOfOther: ("+str(xmax_id)+","+str(ymax_id)+")")
    if out==2:
       print(str(whole)+";"+str(absmax)+";"+str(l)+";"+str(xmax_id)+";"+str(ymax_id))

##########################

# def rysuj(loop, tail):
#     fig, ax = plt.subplots(subplot_kw={'projection': '3d'})

#     datasets = [{"x":[P[1][0] for P in loop], "y":[P[1][1] for P in loop], "z":[P[1][2] for P in loop], "colour": "red"}, {"x":[P[1][0] for P in tail], "y":[P[1][1] for P in tail], "z":[P[1][2] for P in tail], "colour": "blue"}]

#     for dataset in datasets:
#         ax.plot(dataset["x"], dataset["y"], dataset["z"], color=dataset["colour"])

#     plt.show()


def main2():
    # output = "total" # min, middle or total

    # assert(len(sys.argv)>=3)
    parser = argparse.ArgumentParser(description='Welcome to CountGLN!')
    parser.add_argument('-x', '--xyz', required=True, type=str, help='file with RNA structure in .xyz format')
    parser.add_argument('-b', '--bonds', required=True, type=str, help='file with hydrogen bonds')
    parser.add_argument('-m', '--maxbonds', type=int, help='max number of bonds between loop and tail; set to -1 to not use this parameter (any number of bonds is valid)')
    parser.add_argument('--out', type=str, help='type of output format: min/middle/total/df, default: total')
    parser.add_argument('-n', '--name', type=str, help='name of the structure')
    parser.add_argument('--notjoined', action='store_true', help='count GLN only on continuous fragments of tails')
    args = parser.parse_args()

    filenameXYZ = args.xyz
    filenameBonds = args.bonds

    if args.maxbonds==None: maxbonds=0
    else: maxbonds=args.maxbonds

    if args.out==None: output = "total"
    else: output=args.out

    if args.name==None: name=""
    else: name=args.name

    if args.notjoined: joined=False
    else: joined=True
    
    chain,indexes,gaps = ChainRead(filenameXYZ)
    if len(chain)==0:
        print("[WARNING] An empty structure.\n")
        return

    B,types = BondsRead(filenameBonds, indexes)

    LT = FindLoopsTails(B, beg=0, end=len(indexes)-1, gaps=gaps, maxbonds=maxbonds)


    total_lid, total_mmax = 0,[0,0,0]
    for lids in LT:
        loop = PrepareLoop(chain,lids)
        if len(loop) < 3: continue
        tailsN, tailsC = [], []
        beg_loop = lids[0] if type(lids[0]!=list) else lids[0][0]
        end_loop = lids[1] if type(lids[0]!=list) else lids[-1][1]
        for tids in LT[lids]:
            if (tids[1] < beg_loop):
                t = PrepareFragments(chain,[tids])
                if len(t)>0: tailsN.append(t)
            elif (tids[0] > end_loop):
                t = PrepareFragments(chain,[tids])
                if len(t)>0: tailsC.append(t)
            else:
                print('whats up? is there a middle tail?')
        total_minmaxN = GLN_complexTail(loop,tailsN,joined)[0]
        total_minmaxC = GLN_complexTail(loop,tailsC,joined)[0]
        mmax = [0,0,0]

        if (len(tailsN)>0 and abs(total_minmaxN["mmax"][0])>abs(mmax[0])): mmax = total_minmaxN["mmax"]
        if (len(tailsC)>0 and abs(total_minmaxC["mmax"][0])>abs(mmax[0])): mmax = total_minmaxC["mmax"]
        if (abs(mmax[0])>abs(total_mmax[0])): 
            total_mmax = mmax
            total_lid = (indexes[lids[0]],indexes[lids[1]])

        full_lids = (indexes[lids[0]],indexes[lids[1]])

        if (output=="total"):
            print("\nLoop and tails: ",full_lids, LT[lids], ":")
            print("Tail N: ",total_minmaxN)
            print("Tail C: ",total_minmaxC)
        elif (output=="middle"):
            print("Loop:",full_lids," mmax: ",mmax)
        elif (output=="df"):    #as a dataframe
            #PDB;Loop;BTypes;Tails;Nmmax;Nclass;Nwhole;Nmax;NmaxTail;NmaxTails;Nmin;NminTail;NminTails;Cmmax;Cclass;Cwhole;Cmax;CmaxTail;CmaxTails;Cmin;CminTail;CminTails;WithGaps
            withGaps = gaps!=[]
            if LT[lids]:
                full_LT = []
                for piece in LT[lids]:
                    full_LT.append((indexes[piece[0]], indexes[piece[1]]))
                if not total_minmaxN: 
                    total_minmaxN={'mmax':['x'],'cl':'x','wh':'x','max':['x','x','x'],'min':['x','x','x'] }
                if not total_minmaxC: 
                    total_minmaxC={'mmax':['x'],'cl':'x','wh':'x','max':['x','x','x'],'min':['x','x','x'] }
                L=[full_lids,types[lids],full_LT,total_minmaxN['mmax'][0],total_minmaxN['cl'],total_minmaxN['wh'],total_minmaxN['max'][0],total_minmaxN['max'][1],total_minmaxN['max'][2],total_minmaxN['min'][0],total_minmaxN['min'][1],total_minmaxN['min'][2],total_minmaxC['mmax'][0],total_minmaxC['cl'],total_minmaxC['wh'],total_minmaxC['max'][0],total_minmaxC['max'][1],total_minmaxC['max'][2],total_minmaxC['min'][0],total_minmaxC['min'][1],total_minmaxC['min'][2],withGaps]
                L=[str(el) for el in L]
                if name: print(name+";") 
                print(";".join(L)+"\r")

    if output!="df": print("\nTOTAL mmax:", total_mmax[0], ", loop:", total_lid, ", tail:", total_mmax[1:],"\n")


    #rysuj(loop, tail3)


##########################
if __name__ == "__main__":

    main2()

