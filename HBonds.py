#!/usr/bin/python
# -*- coding: utf-8 -*-


def FileCheck(fn):
    try:
        open(fn, "r")
        return 1
    except IOError:
        print("Error: File",fn,"does not appear to exist.\n")
        return 0


def BondsRead(filename, indexes):
    if (not FileCheck(filename)): return {}
    bonds = {}
    f = open(filename)
    types = {}
    not_in_indexes = []
    for line in f.readlines():
        line = line.split('"')
        part1 = line[1].split("|")
        part2 = line[5].split("|")
        b1 = part1[4]
        b2 = part2[4]
        if len(part1) == 9 or len(part2) == 9:   # zly lancuch typu 1_555
            continue
        if len(part1) == 8:
            if part1[7] != '': 
                b1 += part1[7]
        if len(part2) == 8:
            if part2[7] != '': 
                b2 += part2[7]
        if b1 not in indexes or b2 not in indexes:
            if b1 not in indexes and b1 not in not_in_indexes:  # symetryczne = b2 nie trzeba
                not_in_indexes.append(b1)
            continue
        b1 = indexes.index(b1)
        b2 = indexes.index(b2)
        if b1 in bonds:
            if b2 not in bonds[b1]: bonds[b1].append(b2)
        else: bonds[b1] = [b2]
        types[tuple(sorted([b1,b2]))] = line[3]
    f.close()
    if not_in_indexes:
        print(filename, 'Warning! Indexes are present in hbonds file but not in .xyz file. Omitting them. \n', not_in_indexes)
    return bonds, types



def FindLoops(bonds, gaps):
    loops = []
    for (b,bc) in bonds.items():
        for c in bc:
            b1,b2 = min(b,c), max(b,c)
            if any([True for gap in gaps if gap>=b1 and gap<b2]): continue      # dziura w pętli
            if (b1+1 in bonds and (b2-1 in bonds[b1+1] or b2-2 in bonds[b1+1] or b2+1 in bonds[b1+1] or b2+2 in bonds[b1+1] or b2 in bonds[b1+1])): continue
            if (b1+2 in bonds and (b2-1 in bonds[b1+2] or b2-2 in bonds[b1+2] or b2+1 in bonds[b1+2] or b2+2 in bonds[b1+2] or b2 in bonds[b1+2])): continue
            if (b2-1 in bonds and (b1 in bonds[b2-1])): continue
            if (b1,b2) not in loops: loops.append((b1,b2))
    loops.sort()
    return loops


def FindLoopsTails(bonds, beg, end, gaps, minloop = 6, mintail = 4, maxbonds = 0):    
# maxbonds- ile dopuszczalnych wiazan, zeby uznac za nieistotne i polaczyc 2 ogony
# maxbonds = -1 - w ogóle nie dzielimy ogonów
    L = FindLoops(bonds, gaps)
    L = [x for x in L if x[1]-x[0] > minloop-2]
    LT = {}
    for loop in L:
        LT[loop] = []
        if maxbonds == -1:
            tail1 = (beg, loop[0]-1)
            tail2 = (loop[1]+1, end)
            if tail1[1] - tail1[0] > mintail - 2:
                LT[loop].append(tail1)
            if tail2[1] - tail2[0] > mintail - 2:
                LT[loop].append(tail2)
        elif maxbonds >= 0:
            t1 = beg
            nr = beg
            last = float('-inf')
            for nr in range(beg,loop[0]):
                if nr not in bonds: continue
                for b in bonds[nr]:
                    if b in range(loop[0]+1, loop[1]):
                    #wrzucam ten kawalek tail (jesli dlugosc wieksza niz mintail) i przesuwam t1 dalej,
                        if nr > t1 and nr-t1 > mintail-1:
                            if t1-last > maxbonds:
                                LT[loop].append((t1, nr-1))
                            else:        # lacze obecny ogon z poprzednim
                                LT[loop].append((LT[loop].pop()[0], nr-1))
                            last = nr
                        t1 = nr+1
                        break
            
            if nr-t1 > mintail-2:
                if t1-last > maxbonds:
                    LT[loop].append((t1,nr))
                else:        # lacze obecny ogon z poprzednim
                    LT[loop].append((LT[loop].pop()[0], nr))
            t1 = loop[1]+1
            nr = t1
            for nr in range(loop[1]+1, end+1):
                if nr not in bonds: continue
                for b in bonds[nr]:
                    if b in range(loop[0]+1, loop[1]):
                        #wrzucam ten kawalek tail (jesli dlugosc wieksza niz mintail) i przesuwam t1 dalej,
                        if nr > t1 and nr-t1 > mintail-1:
                            if t1-last > maxbonds or LT[loop][-1][0] < loop[0]:
                                LT[loop].append((t1,nr-1))
                            else:        # lacze obecny ogon z poprzednim
                                print('!', loop, nr)
                                LT[loop].append((LT[loop].pop()[0], nr-1))
                            last = nr
                        t1 = nr+1
                        break
            if nr-t1 > mintail-2:
                if t1-last > maxbonds:
                    LT[loop].append((t1,nr))
                else:        # lacze obecny ogon z poprzednim
                    LT[loop].append((LT[loop].pop()[0], nr))
        else:
            raise ValueError(f"Invalid value of maxbonds: maxbonds={maxbonds}.")

    if gaps == []: return LT 
    for (loop,tails) in LT.items():     # dzielenie ogonów w miejscach gaps
        T = []
        for i in range(len(tails)):
            tail = tails[i]
            for g in gaps:
                if tail[0] <= g and tail[1] > g:
                    if g-tail[0] > mintail-2:
                        T.append((tail[0], g))
                    tail = (g+1, tail[1])
            if tail[1]-tail[0] > mintail-2:
                T.append((tail[0], tail[1]))
        LT[loop] = T

    return LT     

def InLoop(nr,loop):    # czy indeks nr należy do ktoregos z przedzialow petli loop?
    for piece in loop:
        if nr >= piece[0] and nr <= piece[1]:
            return True
    return False

