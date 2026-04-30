#!/usr/bin/env python

import sys, os, re, tempfile, csv, pysam, json, binascii, argparse, subprocess, gzip
import sqlite3
import pandas as pd
import pyranges as pr
import numpy as np
from time import gmtime, strftime
from cyvcf2 import VCF
from pathlib import Path
import code

__version__ = '2.0.0'

ACROCENTRICS = ['13','14','15','21','22']


def sort_chrompos(row,chrom='Chromosome',pos='Start'):
    # Extract the numeric part from the 'Chromosome' value
    if row[chrom].startswith('chr'):
        num_part = row[chrom][3:]
        if num_part.isdigit():
            chromosome_value = int(num_part)
        else:
            chromosome_value = float('inf')
    else:
        chromosome_value = float('inf')
    
    return (chromosome_value, row[pos])


#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
#
# Script
#
#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

parser = argparse.ArgumentParser(description='Filter CNVs')
parser.add_argument('--ploidy',required=True,type=str,help='ploidy tsv from wakhan')
parser.add_argument('-i','--wakhan_dir',required=True,type=str,help='wakhan output dir')
parser.add_argument('--genebed',type=str,help='Gene bed file')
parser.add_argument('--svgenebed',type=str,help='sv gene report')
parser.add_argument('--outfile',required=True,type=str,default=None,help='output filename')
parser.add_argument('--cb', type=str, help='cytobands bedfile')
parser.add_argument('-s', '--min_cnv_size', required=True, default=1000000, type=int, help='min cnv size')
parser.add_argument('-c', '--centromeres', required=True)

args = parser.parse_args()



# dataframe with CNVs
svs = pd.DataFrame(columns=['category','type','chrom1','pos1','chrom2','pos2','length','csyntax','psyntax','genes','filters','id','abundance','info'])


#########################################
#
# Collect coverage metrics
#
#########################################

#with open(covReport, 'r') as f:
#    lines = f.readlines()

#csgenes='/storage1/fs1/timley/Active/aml_ppg/tmp/abelhj/spencerlab/cle-chromoseq/assets/data/blood_cancer_files/chromoseq_genes.bed.gz'
#cssvgenes='/storage1/fs1/timley/Active/aml_ppg/tmp/abelhj/spencerlab/cle-chromoseq/assets/data/blood_cancer_files/chromoseq_sv_genes.bed.gz'
geneCovDf=pd.read_csv(args.genebed, sep="\t", header=None, names=['Chromosome','Start','End','Gene','Info'])
covDf=pd.read_csv(args.svgenebed, sep="\t", names=['Chromosome','Start','End','Gene','Info'])
covDf=pd.concat([geneCovDf, covDf], axis=0).drop_duplicates().reset_index()

#
# get transcript info for genes and sv genes
#
ids = geneCovDf['Info'].str.split('|',expand=True).replace('\s\+\s\S+','',regex=True).loc[:,2:]
ids.columns = ['Transcript','Region','cdsStart','cdsEnd','strand']
geneTrx = pd.concat([geneCovDf,ids],axis=1)[['Gene','Transcript','cdsStart','cdsEnd','strand']].drop_duplicates().reset_index()
geneTrx['cdsStart'] = geneTrx['cdsStart'].astype(int)
geneTrx['cdsEnd'] = geneTrx['cdsEnd'].astype(int)

#covDf['Region'] = covDf['Info'].str.split('|').apply(lambda x: x[3] if len(x) > 2 else '')

geneList = covDf['Gene'].unique().tolist()
reportableCnvGeneList = geneTrx['Gene'].drop_duplicates().tolist()

#cytobands=pd.read_table(args.cb, sep="\t", header=None)
#cbfile='/storage1/fs1/timley/Active/aml_ppg/tmp/abelhj/spencerlab/cle-chromoseq/assets/data/hg38.cytoBandIdeo.bed.gz'
cytobands=pd.read_table(args.cb, sep="\t", header=None, 
names=['Chromosome', 'Start', 'End', 'Band', 'Col'])

centromeres=pd.read_table(args.centromeres, sep="\t", header=None,
names=['Chromosome', 'Start', 'End'])


########################
#
# Get CNVs
#
########################

print("Gathering CNVs...",file=sys.stderr)
vcfs=[args.clonal_vcf, args.subclonal_vcf]
vcf_labs=['clonal', 'subclonal']


for jj in range(len(vcfs)):
  with open(vcfs[jj], 'rt') as f:
    ct=0
    for line in f:

      ll=line.strip().split('\t')
      if ll[0][0]=='#':
        continue
      
      ct=ct+1
      dd = dict(s.split('=') for s in ll[7].split(';'))
      [xx, vartype, chr1, ss]=ll[2].split(':')
      ss=ss.split('-')
      [pos1, pos2] = [int(ss[0]), int(ss[1])]
      chr2 = chr1
      svlen = pos2 - pos1 + 1
      filter=ll[6]
      
      # get cytobands
      bands = 'None'
      bandstring = 'None'
      bands = pr.PyRanges(cytobands).intersect(pr.PyRanges(chromosomes = chr1,starts = [pos1],ends = [pos2])).df
      bands=bands['Band'].tolist()
      if len(bands)>0:
          bandstring=bands[0]
          if len(bands)>1:
              bandstring=bands[0]+'-'+bands[-1]


      #code.interact(local=locals())
      cen = pr.PyRanges(centromeres).intersect(pr.PyRanges(chromosomes = chr1,starts = [pos1],ends = [pos2])).df
      if cen.shape[0]>0:
        cenlen=cen['End'][0]-cen['Start'][0]
#        code.interact(local=locals())
        if svlen<cenlen*1.5:
          continue
        #code.interact(local=locals())
    
      # gene by overlap between variant and genes covqcdf. This is all we need for this resolution.
      genestring = 'None'
      genes = pr.PyRanges(covDf).intersect(pr.PyRanges(chromosomes = chr1,starts = [pos1],ends = [pos2])).df
      if genes.shape[0] > 0:
          genes = genes['Gene'].unique().tolist()
          if len(genes) == 0:
              genestring = 'None'
          elif len(genes) > 10:
              genestring = str(len(genes)) + " genes"
          else:
              genestring = ",".join(genes)
              
      csyntax = '.'
      psyntax = '.'
      if vartype == 'LOSS':
          csyntax = chr1 + ":g." + str(pos1) + "_" + str(pos2) + "del"
          if bands[0].find('p') > -1 and bands[-1].find('q') > -1: # if the CNA spans the centromere then the whole chromosome is lost/gained
              psyntax = "seq[GRCh38] -" + chr1.replace('chr','')
              
          elif 'q11' in bands and 'qter' in bands and chr1 in ACROCENTRICS:
              psyntax = "seq[GRCh38] -" + chr1.replace('chr','')

          elif bands[0].find('p') > -1:
              bands.reverse()
              psyntax = "seq[GRCh38] del(" + chr1.replace('chr','') + ")(" + bands[0] + bands[-1] + ")"
              
          else:
              psyntax = "seq[GRCh38] del(" + chr1.replace('chr','') + ")(" + bands[0] + bands[-1] + ")"
          
      elif vartype == 'GAIN':
          csyntax = chr1 + ":g." + str(pos1) + "_" + str(pos2) + "dup"
          if bands[0].find('p') > -1 and bands[-1].find('q') > -1:
              psyntax = "seq[GRCh38] +" + chr1.replace('chr','')

          elif 'q11' in bands and 'qter' in bands and chr1 in ACROCENTRICS:
              psyntax = "seq[GRCh38] +" + chr1.replace('chr','')

          elif bands[0].find('p') > -1:
              bands.reverse()
              psyntax = "seq[GRCh38] dup(" + chr1.replace('chr','') + ")(" + bands[0] + bands[-1] + ")"
              
          else:
              psyntax = "seq[GRCh38] dup(" + chr1.replace('chr','') + ")(" + bands[0] + bands[-1] + ")"
          
      elif vartype == 'CNLOH':
          csyntax = chr1 + ":g." + str(pos1) + "_" + str(pos2) + "cnLOH"
          if bands[0].find('p') > -1 and bands[-1].find('q') > -1:
              psyntax = "seq[GRCh38] +" + chr1.replace('chr','')

          elif 'q11' in bands and 'qter' in bands and chr1 in ACROCENTRICS:
              psyntax = "seq[GRCh38] +" + chr1.replace('chr','')

          elif bands[0].find('p') > -1:
              bands.reverse()
              psyntax = "seq[GRCh38] cnLOH(" + chr1.replace('chr','') + ")(" + bands[0] + bands[-1] + ")"
              
          else:
              psyntax = "seq[GRCh38] cnLOH(" + chr1.replace('chr','') + ")(" + bands[0] + bands[-1] + ")"

      #code.interact(local=locals())
      format=ll[8].split(':')
      sample=ll[9].split(':')
      sample=dict([format[ii], sample[ii]] for ii in range(len(format)))
      copynumber=sample['TCN']
      abundance='NA'
      category = 'CNV'
      if filter!='PASS':
          category = 'OtherSv'
      if svlen < args.min_cnv_size:
          category = 'SubthresholdCNV'
      if vcf_labs[jj] == 'subclonal':
          category ='subclonal'
      
      
      svs = pd.concat([svs,pd.DataFrame([dict(zip(svs.columns,[category,vartype,chr1,str(pos1),chr1,str(pos2),
str(svlen),csyntax,psyntax,genestring,filter,ll[0]+'_'+ll[1]+'_'+ll[2],str(abundance)+"%",'CN='+str(copynumber)]))])])




svs.to_csv(args.outfile, sep="\t", header=True, index=False)
