#!/usr/bin/env python3
import argparse, sys, os, csv, re
import pysam
import code
import pandas as pd
import pyranges as pr

__version__ = '1.0.0'

def checkfile(file_path):
    """Check if a file exists at the given path."""
    if not os.path.exists(file_path):
        raise argparse.ArgumentTypeError(f"The file {file_path} does not exist.")
    return file_path

def fileexists(file_path):
    """Check if a file exists at the given path."""
    if os.path.exists(file_path):
        raise argparse.ArgumentTypeError(f"The outfile {file_path} exists!")
    return file_path

# write a function to parse the VEP CSQ VCF header
def vepHeader(header):
    fields = header.info.get('CSQ').record.get('Description').split(':')[1].replace('"',"").strip().split('|')
    return dict(zip(fields,range(0,len(fields))))

def svpackHeader(header):
    fields = header.info.get('BCSQ').record.get('Description').split(':')[1].replace('"',"").strip().split('|')
    return dict(zip(fields,range(0,len(fields))))


def parseVepCsq(csq, vep, field=None):
    if not isinstance(csq, tuple):
        raise TypeError("csq should be a string")
    if not isinstance(vep, list):
        raise TypeError("vep should be a list")
    out = []
    for val in csq:
        fields = val.split("|")    
        if len(fields) < len(vep):
            raise ValueError("VEP header has more fields than the CSQ record")
        csq_dict = dict(zip(vep, fields))
        out.append(csq_dict)
    if field:
        out = list(set([d[field] for d in out if field in d and d[field]!='']))
    return out

def main():
    
    parser = argparse.ArgumentParser(description="Filter VCF records based on criteria")
    parser.add_argument("vcffile", type=checkfile, help="Input VCF file")
    parser.add_argument("-g", "--geneList", type=checkfile, required=True, help="Gene list file")
    parser.add_argument("-c", "--cnvgeneList", type=checkfile, required=True, help="CNV Gene list file")
    parser.add_argument("-b", "--svgeneList", type=checkfile, required=True, help="SV Gene list file")
    parser.add_argument("--caller", type=str, default='severus', help='sv caller, sniffles or severus')
    parser.add_argument("-l", "--minSvLength", type=int, default=5000, help="Minimum SV length")
    parser.add_argument("-L", "--maxSvLength", type=int, default=1e9, help="Maximum SV length")
    parser.add_argument("-m", "--minSvReads", type=int, default=2, help="Minimum SR and PR alt-supporting reads")
    parser.add_argument("-s", "--slop", type=int, default=20000, help="slop around non-vep annotated variants")
    parser.add_argument("-a", "--minSvAbundance", type=float, default=5.0, help="Minimum SV abundance in percent")
    parser.add_argument("-o", "--outfile", type=str, help="Outfile")
    parser.add_argument('-v', '--version', action='version', version='%(prog)s: ' + __version__)
    parser.add_argument('-n', '--ins', type=checkfile, required=True, help='bed file of normal insertions')

    nonSynon = {"splice_acceptor_variant","splice_donor_variant","stop_gained","frameshift_variant","stop_lost","start_lost","transcript_ablation","transcript_amplification","inframe_insertion","inframe_deletion","missense_variant","protein_altering_variant"}
    nonSynonb = {'sv:cds', 'sv:bnd', 'sv:utr'}

    args = parser.parse_args()

    outfile = args.outfile
    if outfile is None:
        outfile = sys.stdout

    recurrentSvs = pd.read_csv(args.geneList,sep=',', header=None, names = ['gene1','strand1','gene2','strand2','type'])
    gl=recurrentSvs['gene1'].tolist()
    gl.extend(recurrentSvs['gene2'].tolist())

    chromoseqGenes=pd.read_csv(args.cnvgeneList, header=None, skiprows=0, names=['Chromosome','Start','End','Gene','Info'],sep="\t")
    chromoseqGenes=chromoseqGenes['Gene'].unique().tolist()
    gl.extend(chromoseqGenes)

    svBed = pd.read_csv(args.svgeneList, header=None, skiprows=0, names=['Chromosome','Start','End','Gene','Info'],sep="\t")
    gl.extend(svBed['Gene'].tolist())

    insBed = pd.read_csv(args.ins, header=None, skiprows=0, names=['Chromosome','Start','End', 'Variants'], sep="\t")

    knownGenes=set(gl)
    cseqGenes=set(chromoseqGenes)
    
    passingRecords = set()
    knownSvGenes = {}
    recSv={}
    recCnv={}
    recurrent={}
    pon={}

    # Open VCF file
    vcf_in = pysam.VariantFile(args.vcffile, "r")
    vep = vepHeader(vcf_in.header)
    svpack = svpackHeader(vcf_in.header)


    ct=0
    for record in vcf_in:
        # this gets genes that overlap SVs from both the standard VEP annotation to get upstream/downstream events and custom BED overlap annotations
        # so we dont miss anything. It also annotates variants with gene clusters, like IGH, IGL, TRA, etc.
        ct=ct+1
        sys.stderr.write(str(ct)+'\n')
        svtype = record.info.get("SVTYPE")
        genes=[]
        genesb=[]
        genes_itx=[]
        consequences=[]
        consequencesb=[]

        [chr2, pos2] = [record.chrom, record.stop]
        if svtype=='BND':
            gg=re.search('[\[\]](chr[^:]*):([0-9]+)', record.alts[0]).groups()
            [chr2, pos2]=[gg[0], int(gg[1])]
            
        if 'CSQ' in record.info.keys():
            genes = list(set(parseVepCsq(record.info.get("CSQ"),list(vep.keys()),'SYMBOL') + parseVepCsq(record.info.get("CSQ"),list(vep.keys()),'KnownSvGenes')))
            consequences = '&'.join(parseVepCsq(record.info.get("CSQ"),list(vep.keys()),'Consequence')).split('&')
        else:     #if variant was ignored by vep
            if svtype in ['BND', 'INV']:
                genes1a = pr.PyRanges(svBed).intersect(pr.PyRanges(chromosomes = str(record.chrom),starts = [record.start-args.slop],ends = [record.start+args.slop])).df
                genes1b = pr.PyRanges(svBed).intersect(pr.PyRanges(chromosomes = chr2,starts = [pos2-args.slop],ends = [pos2+args.slop])).df
                if genes1a.shape[0] > 0:
                    genes_itx=genes1a['Gene'].unique().tolist()
                elif genes1b.shape[0] > 0:
                    genes_itx.extend(genes1b['Gene'].unique().tolist())
            if svtype in ['DEL', 'DUP', 'INS']:
                genes1a=pr.PyRanges(svBed).intersect(pr.PyRanges(chromosomes = str(record.chrom),starts = [record.start-args.slop],ends = [pos2+args.slop])).df
                if genes1a.shape[0] > 0:
                    genes_itx=genes1a['Gene'].unique().tolist()

        if 'BCSQ' in record.info.keys():
            genesb = list(set(parseVepCsq(record.info.get("BCSQ"),list(svpack.keys()),'gene')))
            consequencesb = '&'.join(parseVepCsq(record.info.get("BCSQ"),list(svpack.keys()),'Consequence')).split('&')

        pon_filter = False
        if svtype in ['BND', 'INS']:
            ins1 = pr.PyRanges(insBed).intersect(pr.PyRanges(chromosomes = str(record.chrom),starts = [record.start-args.slop],ends = [record.start+args.slop])).df
            ins2=ins1.copy()
            if svtype=='BND':
                ins2=pr.PyRanges(insBed).intersect(pr.PyRanges(chromosomes = chr2,starts = [pos2-args.slop],ends = [pos2+args.slop])).df
            if ins1.shape[0]>0 or  ins2.shape[0]>0:
                pon_filter = True
        if 'CCDGID' in record.info.keys() or 'HPRC_ID' in record.info.keys():
            pon_filter = True

                
        # skip dels/dups that are too small
        if svtype in ['DEL','DUP', 'INS'] and record.info.get("SVLEN") is not None:
            if abs(record.info.get("SVLEN")) < args.minSvLength:
                continue

        DV = 0
        if 'DV' in record.format.keys():
            DV = record.samples[0]['DV']
            # some records have no read support. Skip these.
            if DV == 0:
                continue

        # csq has a known gene
        recordVariant = False
        geneHits = ( set(genes) | set(genes_itx) | set(genesb) )  & knownGenes
        cseqgeneHits = (set(genes) | set(genes_itx ) | set(genesb)  ) & cseqGenes
        
        isRecurrent=False
        if len(geneHits)>0:
            genes1a = pr.PyRanges(svBed).intersect(pr.PyRanges(chromosomes = str(record.chrom),starts = [record.start-args.slop],ends = [record.start+args.slop])).df
            genes1b = pr.PyRanges(svBed).intersect(pr.PyRanges(chromosomes = chr2,starts = [pos2-args.slop],ends = [pos2+args.slop])).df
            if genes1a.shape[0]>0 and genes1b.shape[0]>0:
                genes1a=genes1a['Gene'].unique().tolist()
                genes1b=genes1b['Gene'].unique().tolist()
                if len(recurrentSvs[(recurrentSvs['gene1'].isin(genes1a)) & (recurrentSvs['gene2'].isin(genes1b))])>0 or \
                   len(recurrentSvs[(recurrentSvs['gene1'].isin(genes1b)) & (recurrentSvs['gene2'].isin(genes1a))])>0:
                    isRecurrent = True
                    recurrent[record.id]=1

        # get records that overlap a known gene
        if len(geneHits) > 0 and svtype in ['BND', 'INV']:
            recordVariant = True

        elif len(geneHits) > 0 and svtype in ['INS','DEL','DUP'] and ( len(set(consequences) & nonSynon)>0 or len(set(consequencesb) & nonSynonb)>0):

            recordVariant = True
            
        # also get BND records where the mate has been recorded (this is the main purpose of this script) 
        if svtype=='BND' and (('MATEID' in record.info.keys() and record.info.get("MATEID")[0] in passingRecords) or record.id in passingRecords ):
            recordVariant = True

        # finally get passing records with a functional consequence if there is a contig and there are at least minReads alternate-supporting SR and PR reads
        if len(record.filter)==0 and (len(set(consequences) & nonSynon)>0 or len(set(consequencesb) & nonSynonb)>0) and DV > args.minSvReads and not pon_filter:
            recordVariant = True

        if recordVariant:
            passingRecords.add(record.id)
            if len(geneHits) > 0:
                recSv.setdefault(record.id, []).extend(list(geneHits))
                recCnv.setdefault(record.id, []).extend(list(cseqgeneHits))
            if pon_filter:
                pon[record.id]=True

    vcf_in.close()

    # Re-open the VCF and print header and records that are in passingRecords
    vcf_in = pysam.VariantFile(args.vcffile, "r")

    if 'KnownSvGenes' not in vcf_in.header.info.keys():
        vcf_in.header.info.add("KnownSvGenes", '.', 'String', 'List of recurrent SV genes that overlap this variant')

    if 'RecurrentSV' not in vcf_in.header.info.keys():
        vcf_in.header.info.add('RecurrentSV', '.', 'String', 'Recurrent SV')

    if 'ChromoseqGenes' not in vcf_in.header.info.keys():
        vcf_in.header.info.add("ChromoseqGenes", '.', 'String', 'List of Chromoseq genes that overlap this variant')

    if 'Imprecise' not in vcf_in.header.filters.keys():
        vcf_in.header.filters.add("Imprecise",None,None,'No contig found so breakend are imprecise')

    if 'MinSvReads' not in vcf_in.header.filters.keys():
        vcf_in.header.filters.add("MinSvReads",None,None,f'Fails minimum SR or PR alt-supporting reads (${args.minSvReads})')

    if 'MinSvAbundance' not in vcf_in.header.filters.keys():
        vcf_in.header.filters.add("MinSvAbundance",None,None,f'Fails minimum SV abundance (${args.minSvAbundance})')

    if 'PONFilter' not in vcf_in.header.info.keys():
        vcf_in.header.info.add("PONFilter", '.', 'String' , 'Overlaps insertion from PON.')


    header = vcf_in.header.copy()

    vcf_out = pysam.VariantFile(outfile, "w", header=header)

    for record in vcf_in:
        printRecord = False

        if record.id in passingRecords:
            printRecord = True
        elif record.info.get("SVTYPE") == 'BND' and ('MATEID' in record.info.keys() and  record.info.get("MATEID")[0] in passingRecords):
            printRecord = True

        if printRecord:

            [DV, DR] = [0,0]
            VAF=0.0
            if 'DV' in record.format.keys():
                DV = record.samples[0]['DV']
            if 'VAF' in record.format.keys():
                VAF = record.samples[0]['VAF']
            elif 'DR' in record.format.keys():
                DR = record.samples[0]['DR']
                VAF = DV/(DV+DR)

            if DV < int(args.minSvReads):
                record.filter.add('MinSvReads')

            if VAF * 100.0 < float(args.minSvAbundance):
                record.filter.add('MinSvAbundance')

            new_record = vcf_out.new_record()

            new_record.chrom = record.chrom
            new_record.pos = record.pos
            new_record.stop = record.stop
            new_record.id = record.id
            new_record.alleles = record.alleles

            for k in record.filter.keys():
                new_record.filter.add(k)

            for k in record.info.keys():
                new_record.info[k] = record.info.get(k)

            if new_record.id in recSv.keys():
                if len(recSv[new_record.id])>0:
                    new_record.info['KnownSvGenes'] = ','.join(list(set(recSv[new_record.id])))

            if new_record.id in recCnv.keys():
                if len(recCnv[new_record.id])>0:
                    new_record.info['ChromoseqGenes'] = ','.join(list(set(recCnv[new_record.id])))

            if new_record.id in recurrent.keys():
                new_record.info['RecurrentSV']="recurrent"

            if new_record.id in pon.keys():
                new_record.info['PONFilter']="True"
                
            for k in record.format.keys():
                new_record.samples[0][k] = record.samples[0][k]
            
            vcf_out.write(new_record)

    vcf_in.close()
    vcf_out.close()

if __name__ == "__main__":
    main()
