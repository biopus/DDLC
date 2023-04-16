#标准库
import argparse
import os
import sys
import shutil
import subprocess



#第三方库
from Bio import Phylo
import matplotlib


#自定义模块
from lib.build_database import fetch_gb_file, parse_gene_in_txt, parse_its_in_txt, parse_HA_in_txt
from lib.check import check_rf_marker, check_raw_gb, check_query_marker
from lib.trim import my_trim


def main(args):
    genus = args.genus
    marker = args.marker
    rfdb = args.rfdb
    query_ready = args.query_ready
    forward = args.fw
    reverse = args.rs
    cwd = os.getcwd()
    reference_path = os.path.join(cwd,'reference',marker+'.fasta')
    query_path = os.path.join(cwd,'query',marker+'.fasta')
    my_prog_path = os.path.abspath(__file__)
    my_prog_dir = os.path.dirname(my_prog_path)
    
    ###build reference marker###
    #如果有准备好的参考marker的fasta文件
    if rfdb == '1':
        pass
    #如果只有包含参考marker的gb文件
    elif rfdb == '2':
        if not os.path.exists('reference'):
            os.mkdir('reference')
        if check_raw_gb(marker):
            if marker in ['matk','rbcl']:
                parse_gene_in_txt(marker)
            elif marker == 'its':
                parse_its_in_txt(marker)
            elif marker == 'trnH-psbA':
                parse_HA_in_txt(marker)
            else:
                sys.stdout.write("We don't support this marker!")
    #如果什么也没有
    elif rfdb == '3':
        if not os.path.exists('reference'):
            os.mkdir('reference')
        if not os.path.exists('raw_gb'):
            os.mkdir('raw_gb')
        fetch_gb_file(genus,marker)    
        if marker in ['matk','rbcl']:
            parse_gene_in_txt(marker)
        elif marker == 'its':
            parse_its_in_txt(marker)
        elif marker == 'trnH-psbA':
            parse_HA_in_txt(marker)
        else:
            sys.stdout.write("We don't support this marker!")
    check_rf_marker(marker)

    ###build query marker###
    #find query marker in sequencing reads with geneminer
    if not os.path.exists('query'):
        os.mkdir('query')
    if not query_ready:
        geneminer_cmd = f"geneminer.py -1 {forward} -2 {reverse} -rtfa reference/{marker}.fasta -o {marker} -t 4"
        p = subprocess.run(geneminer_cmd,shell=True)
        with open(os.path.join(os.getcwd(),marker,'GM_results',f'{marker}.fasta'),'r') as output_handle:
            des = output_handle.readline()
            seq = output_handle.readline()
        with open(query_path,'w') as input_handle:
            input_handle.write(des)
            input_handle.write(seq)
    check_query_marker(marker)
    shutil.rmtree(marker)

    ###identify with phylogenetic relationship###
    #multipul alignments
    shutil.copy(reference_path,'data.fasta')
    with open('data.fasta','a') as in_handle, open(query_path,'r') as out_handle:
        in_handle.write(out_handle.read())
    muscle_cmd = f"{my_prog_dir}/muscle -align data.fasta -output aligned_data.fasta"
    subprocess.run(muscle_cmd,shell=True)

    #trim aligned result
    #设置裁剪分数
    Sg = 0.8
    Ss = 0.6
    my_trim(Sg,Ss)

    #build tree
    #way1：直接用raxml的GRTGAMMA模型
    if not os.path.exists('raxml_out'):
        os.mkdir('raxml_out')
    os.chdir('raxml_out')
    raxml_cmd = f"{my_prog_dir}/raxml -f a -m GRTGAMMA -p 12345 -x 12345 -# 1000 -s {cwd}/trimmed_data.fasta -n out" 
    subprocess.run(raxml_cmd,shell=True)
    os.chdir('..')
    os.remove('trimmed_data.fasta.reduced')

    #way2:使用modeltestng+raxml
    if not os.path.exists('raxmlng_out'):
        os.mkdir('raxmlng_out')
    

    #way3:使用raxml-ng

    #way4:使用iqtree

    #可视化
    Phylo.convert('RAxML_bipartitions.out','newick','test_tree','phyloxml')
    
    

    ###identify with piarwise alignment###


    ###identify with blast###
                




if __name__ == "__main__":
    parser = argparse.ArgumentParser(usage="%(prog)s <genus> -m [{its,matk,rbcl}]",formatter_class=argparse.RawDescriptionHelpFormatter,description="name:cmidentifier\n"
                                     "version:1.0")
    parser.add_argument("-g",metavar="<str>",dest="genus",type=str,required=True,help="specify the probable genus")
    parser.add_argument("-m",metavar="<str>",dest="marker",type=str,choices=["its","matk","rbcl","trnH-psbA"],required=True,help="Specify marker.")
    parser.add_argument("-rf",metavar="<str>",dest="rfdb",choices=['1','2','3'],required=True,help="""choose 1 if you have marker in fasta format in the reference directory.
    choose 2 if you have raw gb file in raw_gb directory. choose 3 if you want to use this prog to downlownd""")
    parser.add_argument("-q",dest="query_ready",action="store_true",help="if you already have query marker,choose this")
    parser.add_argument("-fw",metavar='<str>',help="if you don't have query marker, choose the file with forward paired-end reads (*.fq/*.fq.gz)")
    parser.add_argument("-rs",metavar='<str>',help="if you don't have query marker, choose the file with reverse paired-end reads (*.fq/*.fq.gz)")
    args = parser.parse_args()
    main(args)