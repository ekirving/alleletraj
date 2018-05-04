#!/usr/bin/env python
# -*- coding: ASCII -*-

"""

:Author: Martin Kircher
:Contact: Martin.Kircher@eva.mpg.de
:Date: *10.03.2010

"""

import os,sys
from optparse import OptionParser
import string
import math
table = string.maketrans('.N','AA')

# sys.path.append('/home/sbsuser/lib/python/')
#from bx.seq import twobit

parser = OptionParser()
parser.add_option("-p","--pp",dest="propPaired",help="Filter for proper pair flag",default=False,action='store_true')
#parser.add_option("-t", "--twobit", dest="twobit", help="Read genome from 2bit FILE",default="")
#parser.add_option("--maxlength",dest="maxlength",help="Maximum insert size to be extracted from the genome (default 1000)",default=1000,type="int")
parser.add_option("--length_cutoff",dest="length",help="Consider only sequences of at least X nt (default X = None)",default=None,type="int")
parser.add_option("--frequency_cutoff",dest="frequency",help="Consider only sequences of at least X coverage (default X = None)",default=None,type="int")
parser.add_option("--no_length",dest="nolength",help="Do not consider length of merged reads",default=False,action='store_true')
parser.add_option("-f","--5p",dest="fivePrime",help="Cluster single reads on five prime coordinate",default=False,action='store_true')
parser.add_option("--5p_max_length",dest="fivePrimeML",help="Maximum length of single reads (needed for five prime clustering, def 500)",default=500,type="int")
parser.add_option("-c","--cutoff",dest="cutoff",help="Entropy complexity cutoff (default = None)",type='float')
parser.add_option("-k","--keep",dest="keep",help="Keep unmapped sequences",default=False,action='store_true')
parser.add_option("--buffer",dest="rbuffer",help="Lowest number of PE reads buffered before write (def 50000)",default=50000,type="int")
parser.add_option("--quality",dest="quality",help="Quality score offset (def 33)",default=33,type="int")
parser.add_option("-m","--merged",dest="merged",help="Require SR reads to be merged",default=False,action='store_true')
parser.add_option("-v","--verbose",dest="verbose",help="Print progress messages",default=False,action='store_true')
(options, args) = parser.parse_args()

have_genome = False
genome = None
#if os.path.isfile(options.twobit):
#  genome = twobit.TwoBitFile( open( options.twobit ) )
#  have_genome = True
#  sys.stderr.write("Opened 2bit file.\n")

def calc_GC(seq):
  GC = 0
  for elem in seq: 
    if elem in "GC": GC += 1
  return GC/float(len(seq))

def is_complex_entropy(seq,cutoff=options.cutoff):
  global table
  if cutoff == None: return True
  seq=seq.translate(table)
  counts = []
  for x in 'ACGT':
    counts.append(seq.count(x))
  total = float(sum(counts))
  entropy = 0
  for elem in counts:
    if (total > 0) and (elem/total <> 0):
      entropy -= elem/total*math.log(elem/total,2)
  return (entropy >= cutoff)

def extractLengthCigar(seq):
  digit = ""
  length = 0
  for i in seq:
    if i.isdigit(): digit+=i
    else:
      if len(digit) > 0:
        if i not in "SHIP": length += int(digit)
        digit = ""
      else:
        sys.stderr.write("Unexpected cigar line: %s\n"%seq)
  return length

def split_cigar(seq):
  digit = ""
  cigar = []
  for i in seq:
    if i.isdigit(): digit+=i
    else:
      if len(digit) > 0:
        cigar.append((int(digit),i))
        digit = ""
      else:
        sys.stderr.write("Unexpected cigar line: %s\n"%seq)
  return cigar

def aln_length(cigar):
  res = 0
  for value,ctype in cigar:
    if ctype == "D" or ctype == "M": res += value
  return res

strand_flag = 16
first_read_flag = 64
second_read_flag = 128
proper_pair_flag = 2
paired_flag = 1
query_unmapped_flag = 4
mate_unmapped_flag = 8

fid = 0
fflag = 1
fchr = 2
fstart = 3
fmqual = 4
fcigar = 5
fschr = 6
fsstart = 7
fisize = 8
fseq = 9
fqual = 10

def calc_consensus(lines):
  count = len(lines)
  outline = lines[0][:fqual+1]
  seq2mqualofields = {} # KEEP MAPQ AND FIELDS FROM ENTRY WITH IDENTICAL SEQUENCE
  for elem in lines:
    ofields = []
    for field in elem[fqual+1:]: # LOOK FOR PREVIOUS PCR DUPLICATE COUNTS
      if field.startswith('XP:i:'):
        count += int(field.split(':')[-1])-1
      else:
        ofields.append(field.rstrip("\n"))
    seq2mqualofields[elem[fseq]] = (elem[fmqual],ofields)

  if len(lines) > 1:
    # DETERMINE CONSENSUS SEQUENCE FOR THIS VARIANT
    seqstring,qualstring = "",""
    for pos in range(len(outline[fseq])):
      bases = [0,0,0,0]
      bcount = 0
      base,qualchr = None,None
      for elem in lines:
        base = elem[fseq][pos]
        qualchr = elem[fqual][pos]
        if base == 'N': continue
        bcount += 1
        qual = (ord(elem[fqual][pos])-options.quality)/-10.0
        if qual == 0: qual = -0.1
        rev_qual = math.log10(1.0-10**qual)-math.log10(3.0)
        for i,b in enumerate('ACGT'):
          if b == base: bases[i]+=qual
          else: bases[i]+=rev_qual
      if bcount > 1:
        total_prob = math.log10(max(0.000000001,sum(map(lambda x:10**x,bases))))
        max_base,max_qual,min_val = 'N',chr(options.quality),0
        for i,b in enumerate('ACGT'):
          cval = bases[i]-total_prob
          if cval < min_val:
            min_val = cval
            max_base,max_qual = b, chr(int(round(min(60,-10.0*(bases[i]-total_prob))))+options.quality)
        seqstring+=max_base
        qualstring+=max_qual
      else:
        seqstring+=base
        qualstring+=qualchr

    outline[fseq]=seqstring
    outline[fqual]=qualstring
  else:
    seqstring = outline[fseq]

  if seqstring in seq2mqualofields:
    outline[fmqual] = seq2mqualofields[seqstring][0]
    outline.extend(seq2mqualofields[seqstring][1])
  outline.append("XP:i:%d\n"%count)
  return "\t".join(outline)

def get_consensus(lines):
  by_cigar = {}
  cigar_count = {}
  # DETERMINE MOST FREQUENT CIGAR LINE PAIR
  for (line1,line2) in lines:
    fields1 = line1.split("\t")
    fields2 = line2.split("\t")
    cigars = fields1[fcigar]+","+fields2[fcigar]
    if cigars in by_cigar:
      cigar_count[cigars]+=1
      by_cigar[cigars][0].append(fields1)
      by_cigar[cigars][1].append(fields2)
    else:
      cigar_count[cigars]=1
      by_cigar[cigars]=([fields1],[fields2])
  to_sort = map(lambda (y,x): (x,-len(y),y),cigar_count.iteritems())
  to_sort.sort()
  selcigar = to_sort[-1][-1]
  lines = by_cigar[selcigar]
  del by_cigar
  del cigar_count
  del to_sort
  return calc_consensus(lines[0]),calc_consensus(lines[1])

def get_consensus_SR(lines):
  # DETERMINE MOST FREQUENT CIGAR LINE
  by_cigar = {}
  cigar_count = {}
  for line in lines:
    fields = line.split("\t")
    if fields[fcigar] in by_cigar:
      cigar_count[fields[fcigar]]+=1
      by_cigar[fields[fcigar]].append(fields)
    else:
      cigar_count[fields[fcigar]]=1
      by_cigar[fields[fcigar]]=[fields]
  to_sort = map(lambda (y,x): (x,-len(y),y),cigar_count.iteritems())
  to_sort.sort()
  selcigar = to_sort[-1][-1]
  lines = by_cigar[selcigar]
  del by_cigar
  del cigar_count
  del to_sort
  return calc_consensus(lines)

variants = {}
incomplete_variants = {}
curpos = None
curvariants = {}
total_lines = 0
out_lines = 0
out_lines_SR = 0
firstAln = False
if True:
#try:
  for line in sys.stdin:
    if line.startswith("@") and not firstAln:
      sys.stdout.write(line)
    total_lines += 1
    if options.verbose and total_lines % 10000 == 0: sys.stderr.write("Lines in %d / Lines out %d / Lines out SR %d\n"%(total_lines,out_lines,out_lines_SR))
    fields = line.split("\t")
    if len(fields) > 10:
      firstAln = True
      if options.length != None and len(fields[9]) < options.length: continue
      if not is_complex_entropy(fields[fseq]): continue
      if (fields[fchr] == "*" or fields[fcigar] == "*") and not options.keep: continue
      if not (int(fields[fflag]) & paired_flag == paired_flag) and options.merged and not fields[fid].startswith("M_"): continue

      if len(variants) > options.rbuffer:
        if options.verbose: sys.stderr.write("Full buffer (%d)"%len(variants)+str(curpos)+" \n")
        hvariants = {}
        for (hchr,outpos,outpos_r2),lines in variants.iteritems():
          if type(hchr) != type(()) and ((hchr != curpos[0]) or ((hchr == curpos[0]) and (abs(outpos[1]-outpos_r2[1])+max(outpos[1],outpos_r2[1])) < curpos[1])):
            if (options.frequency != None) and (len(lines) < options.frequency): continue
            line1,line2 = get_consensus(lines)
            sys.stdout.write(line1)
            sys.stdout.write(line2)
            out_lines += 2
          elif type(hchr) == type(()) and (((hchr[0] != curpos[0]) and (hchr[1] != curpos[0])) or ((hchr[0] == curpos[0]) and (outpos[1]+options.fivePrimeML < curpos[1])) or ((hchr[1] == curpos[0]) and (outpos_r2[1]+options.fivePrimeML < curpos[1]))):
            if (options.frequency != None) and (len(lines) < options.frequency): continue
            line1,line2 = get_consensus(lines)
            sys.stdout.write(line1)
            sys.stdout.write(line2)
            out_lines += 2
          else:
            hvariants[(hchr,outpos,outpos_r2)]=lines
        variants = hvariants
        if options.verbose: sys.stderr.write("- Full buffer (%d)"%len(variants)+str(curpos)+" \n")
        count = 0

      if (int(fields[fflag]) & paired_flag == paired_flag): # PE DATA
        if not(int(fields[fflag]) & proper_pair_flag == proper_pair_flag) and options.propPaired: continue
        if (int(fields[fflag]) & query_unmapped_flag == query_unmapped_flag) or (int(fields[fflag]) & mate_unmapped_flag == mate_unmapped_flag): continue

        if  ((int(fields[fflag]) & first_read_flag == first_read_flag)): #FORWARD READ
          curpos = (fields[fchr],int(fields[fstart]))

          hchr = fields[fchr]
          cigar = split_cigar(fields[fcigar])
          start = 0
          end = len(fields[fqual])
          if cigar[0][1] == "S": start = cigar[0][0]
          if cigar[-1][1] == "S": end -= cigar[-1][0]
          outpos = (fields[fchr],int(fields[fstart]))
          #qualscore = sum(map(lambda x:ord(x),fields[fqual][start:end]))

          if (int(fields[fflag]) & strand_flag == strand_flag):
            outpos = (fields[fchr],curpos[1]+end-start)

          if fields[fid] not in incomplete_variants:
            incomplete_variants[fields[fid]] = [line,outpos]
          else:
            line_r2,outpos_r2 = incomplete_variants[fields[fid]]
            del incomplete_variants[fields[fid]]
            if outpos_r2[0] != hchr: hchr = hchr,outpos_r2[0]
            if (hchr,outpos,outpos_r2) not in variants:
              variants[(hchr,outpos,outpos_r2)] = [(line,line_r2)]
            else:
              variants[(hchr,outpos,outpos_r2)].append((line,line_r2))
        elif (int(fields[fflag]) & second_read_flag == second_read_flag):  #REVERSE READ
          curpos = (fields[fchr],int(fields[fstart]))

          hchr = fields[fchr]
          cigar = split_cigar(fields[fcigar])
          start = 0
          end = len(fields[fqual])
          if cigar[0][1] == "S": start = cigar[0][0]
          if cigar[-1][1] == "S": end -= cigar[-1][0]
          outpos = (fields[fchr],int(fields[fstart]))
          #qualscore = sum(map(lambda x:ord(x),fields[fqual][start:end]))

          if (int(fields[fflag]) & strand_flag == strand_flag):
            outpos = (fields[fchr],curpos[1]+end-start)

          if fields[fid] not in incomplete_variants:
            incomplete_variants[fields[fid]] = [line,outpos]
          else:
            line_r1,outpos_r1 = incomplete_variants[fields[fid]]
            del incomplete_variants[fields[fid]]
            if outpos_r1[0] != hchr: hchr = outpos_r1[0],hchr
            if (hchr,outpos_r1,outpos) not in variants:
              variants[(hchr,outpos_r1,outpos)] = [(line_r1,line)]
            else:
              variants[(hchr,outpos_r1,outpos)].append((line_r1,line))
        else:
          sys.stderr.write("Should not happen!")
      else: # SR DATA
        if (curpos != None) and ((fields[fchr],int(fields[fstart])) != curpos):
          if options.fivePrime and (fields[fchr] == curpos[0]):
            hpos = int(fields[fstart])-options.fivePrimeML
            hvariants = {}
            for key,value in curvariants.iteritems():
              if (key[1] < hpos):
                if (options.frequency != None) and (value[1] < options.frequency): continue
                #print key,value
                cline = get_consensus_SR(value[0])
                sys.stdout.write(cline)
                out_lines_SR += 1
              else:
                hvariants[key]=value
            curvariants = hvariants
          else:
            for key,value in curvariants.iteritems():
              if (options.frequency != None) and (value[1] < options.frequency): continue
              #print key,value
              cline = get_consensus_SR(value[0])
              sys.stdout.write(cline)
              out_lines_SR += 1
            curvariants = {}
        curpos = (fields[fchr],int(fields[fstart]))

        strand = int(fields[fflag]) & strand_flag == strand_flag
        cigar = split_cigar(fields[fcigar])
        outpos = curpos[1]
        if strand and options.fivePrime: outpos+=aln_length(cigar)

        ## DETERMINE SCORE AFTER CONSIDERING SOFT-CLIP
        #start = 0
        #end = len(fields[fqual])
        #if cigar[0][1] == "S": start = cigar[0][0]
        #if cigar[-1][1] == "S": end -= cigar[-1][0]
        #qualscore = sum(map(lambda x:ord(x),fields[fqual][start:end]))

        nkey = (strand,outpos)
        if not options.nolength and (fields[fid].startswith("M_") or fields[fid].startswith("C_M_")): nkey = (strand,outpos,aln_length(cigar))

        if nkey in curvariants:
          curvariants[nkey][0].append(line)
          curvariants[nkey][1]+=1
        else:
          curvariants[nkey] = [[line],1]
#except:
  #print "Error. Stopping iteration over PIPE..."

for key,value in curvariants.iteritems():
  if (options.frequency != None) and (value[1] < options.frequency): continue
  #print key,value
  line = get_consensus_SR(value[0])
  sys.stdout.write(line)
  out_lines_SR += 1
curvariants = {}

for key,value in variants.iteritems():
  if (options.frequency != None) and (len(value) < options.frequency): continue
  line1,line2 = get_consensus(value)
  sys.stdout.write(line1)
  sys.stdout.write(line2)
  out_lines += 2
