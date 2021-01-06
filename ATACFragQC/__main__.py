import os, sys, re, getopt, pysam
import pandas as pd
import numpy as np
from plotnine import *
from PIL import Image

from ATACFragQC import __version__

class ArgumentList:
    file_bam = ''
    file_ref = ''
    file_out = False
    quality = 50
    isize = 147
    def __init__(self):
        self.file_bam = ''
        self.file_ref = ''
        file_out = False
        self.quality = 50
        self.isize = 147

def bedScan(args):
    print("Processing the reference...")
    ref = pd.read_table(args.file_ref, comment='#', header=None)
    ref.columns = ['seq_id', 'source', 'type', 'start', 'end', 'score', 'strand', 'phase', 'attributes']
    ref = ref[(ref['type'] == 'transcript') & (ref['strand'] != '-') & ref['seq_id'].str.match('^chr[0-9XY]+$') & (ref['start'] > 1000)]
    ref = ref[['seq_id', 'strand', 'start', 'end']].sort_values(by=['seq_id', 'start', 'strand'])
    ref = ref.drop_duplicates(subset=['seq_id', 'start', 'strand'], keep='first')
    ref.loc[ref['strand'] == '+', 'start'] -= 1000
    ref.loc[ref['strand'] == '+', 'end'] = ref.loc[ref['strand'] == '+', 'start'] + 2000
    ref.index = list(range(ref.index.size))
    
    fs = pysam.AlignmentFile(args.file_bam, "rb")
    print("Scaning the distribution of fragments...")
    chr_count = {}
    len_count = [0] * 501
    chr_list = list(set(ref['seq_id']))
    chr_list.append('chrM')
    for chr in chr_list:
        count = 0
        for read in fs.fetch(chr):
            if read.flag == 99 and read.mapq > args.quality and read.isize < 501:
                len_count[read.isize] += 1
                count += 1
        chr_count[chr] = count
    len_count = pd.DataFrame({'V1': list(range(1, 501)), 'V2': len_count[1:]})
    chr_count = pd.DataFrame({'V1': list(chr_count.keys()), 'V2': list(chr_count.values())})
    chr_count = chr_count[chr_count['V1'].str.match('chr.*')]
    chr_list = list(map(lambda x: int(x.replace('chr', '')), list(chr_count.loc[chr_count['V1'].str.match('^chr[0-9]+$'), 'V1'])))
    chr_list.sort()
    chr_list += ['X', 'Y', 'M']
    chr_list = list(map(lambda x: 'chr'+str(x), chr_list))
    chr_count['V1'] = pd.Categorical(chr_count['V1'], categories=chr_list, ordered=True)
    
    print("Scaning the fragments around TSSs...")
    dist_count = []
    for index, row in ref.iterrows():
        count = np.zeros(2001, dtype=np.int64)
        for frag in fs.fetch(row['seq_id'], row['start'], row['end']):
            if frag.flag == 99 and frag.mapq > args.quality and frag.isize < args.isize:
                count[range(max(0, frag.pos - row['start']), min(2001, frag.pos + frag.isize - row['start']))] += 1
        if sum(count) > 0:
            dist_count.append(count)
    dist_count = np.array(dist_count)
    factors = np.mean(dist_count[:, list(range(0,100))+list(range(1901,2001))], axis=1)
    factors[factors == 0] = np.mean(factors)
    dist_count = pd.DataFrame({'V1': list(range(-900, 901)), 'V2': list(np.mean(dist_count[:, 100:1901] / factors.reshape(len(factors), 1), axis=0))})
    
    print("Saving the results...")
    (pathname, extension) = os.path.splitext(args.file_bam)
    (filepath, filename) = os.path.split(pathname)
    ggsave(plot=ggplot(chr_count, aes(x='V1', y='V2'))+geom_bar(stat="identity", width=0.8, fill="#80B1D3")+
        labs(x="Chromosome", y="Fragments")+
        theme(plot_title=element_blank(), panel_background=element_blank(), axis_line=element_line(colour="black"), 
        axis_text_y=element_text(colour="black"), axis_text_x=element_text(angle=270, hjust=0.3, vjust=1, colour="black")), 
        width=4, height=6, dpi=200, filename=pathname+'.tmp0.png', limitsize=False, verbose=False)
    ggsave(plot=ggplot(len_count, aes(x='V1', y='V2'))+geom_bar(stat="identity", colour="#80B1D3", fill="#80B1D3")+
        labs(x="\nInsert Size", y="Fragments")+
        theme(plot_title=element_blank(), panel_background=element_blank(), axis_line=element_line(colour="black"), axis_text=element_text(colour="black")), 
        width=6, height=6, dpi=200, filename=pathname+'.tmp1.png', limitsize=False, verbose=False)
    ggsave(plot=ggplot(dist_count, aes(x='V1', y='V2'))+geom_line(size=1, colour="#80B1D3")+
        labs(x="\nDistance from TSS (bp)", y="Mean TSS enrichment score")+
        scale_x_continuous(breaks=range(-800, 801, 200))+
        theme(plot_title=element_blank(), panel_background=element_blank(), axis_line=element_line(colour="black"), axis_text=element_text(colour="black")), 
        width=6, height=6, dpi=200, filename=pathname+'.tmp2.png', limitsize=False, verbose=False)
    width = 0
    height = 0
    imgs = [Image.open(pathname+'.tmp'+str(i)+'.png') for i in range(3)]
    for img in imgs:
        width += img.width
        height = max(height, img.height)
    result = Image.new('RGBA', (width, height), 'white')
    pos = 0
    for img in imgs:
        result.paste(img, (pos, 0))
        pos += img.width
        img.close()
    result.save(pathname+'_qc.png')
    [os.remove(pathname+'.tmp'+str(i)+'.png') for i in range(3)]
    if args.file_out:
        chr_count.to_csv(pathname+'_chr.tsv', sep='\t', index=False, header=False)
        len_count.to_csv(pathname+'_fl.tsv', sep='\t', index=False, header=False)
        dist_count.to_csv(pathname+'_tss.tsv', sep='\t', index=False, header=False)

def main():
    opts, args = getopt.getopt(sys.argv[1:], 
        'hoi:r:q:l:', 
        ['help', 'output', 'input=', 'reference=', 'quality=', 'length='])
    arguments = ArgumentList()
    help_flag = False
    help_info = 'ATACFragQC - Version: '+__version__+'\n'\
        +'Usage:\nATACFragQC [options] -i <input.bam> -r <reference.gtf>\nArguments:\n'\
        +'-i <file>\tA aligned & deduped BAM file\n'\
        +'-r <file>\tGTF genome annotation\n'\
        +'-o \t\tThe table of results would be saved if -o was set (default: False)\n'\
        +'-q [1-255]\tThe quality limit of alignment (default: 50)\n'\
        +'-l [50-500]\tThe length limit of nucleosome-free fragment (default: 147)\n'
    for opt, arg in opts:
        if opt in ('-h', '--help'):
            help_flag = True
        elif opt in ("-i", "--input"):
            arguments.file_bam = arg
        elif opt in ("-r", "--reference"):
            arguments.file_ref = arg
        elif opt in ("-o", "--output"):
            arguments.file_out = True
        elif opt in ("-q", "--quality"):
            if int(arg) >= 1 and int(arg) <= 255:
                arguments.quality = int(arg)
        elif opt in ("-l", "--length"):
            if int(arg) >= 50 and int(arg) <= 500:
                arguments.isize = int(arg)
    if help_flag or arguments.file_bam == '' or arguments.file_ref == '':
        print(help_info)
        sys.exit()
    bedScan(arguments)

if __name__ == '__main__':
    main()
    file_bam = ''
    file_ref = ''
    file_out = ''
    quality = 0