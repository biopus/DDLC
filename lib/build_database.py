from Bio import Entrez,SeqIO
import os

Entrez.email = "biopus@163.com"

#查找搜索项在nc库中一共有多少条记录
def get_num(search_info):
    handle = Entrez.egquery(term = search_info)
    record = Entrez.read(handle)
    for row in record["eGQueryResult"]:
        if row["DbName"]=="nuccore":
             n = row["Count"]
    return n

#获取查询结果的gi号列表
def get_gi_lis(search_info,num):
    handle = Entrez.esearch(db="nucleotide",term=search_info,retmax = num)
    record = Entrez.read(handle)
    gi_lis = record["IdList"]
    return gi_lis

#下载gb文件到当前工作目录下的raw_gb目录
def fetch_gb_file(genus,marker):
    file_name = marker + '.gb'
    search_info = f"{genus}[Organism] AND {marker}[Gene]"
    term_num = get_num(search_info)
    gi_lis = get_gi_lis(search_info,term_num)
    
    with open(os.path.join(os.getcwd(),'raw_gb',file_name),'w') as f:
        for gi in gi_lis:
            handle = Entrez.efetch(db='nucleotide', id=gi, rettype='gb', retmode='text')
            record = SeqIO.read(handle, "genbank")
            SeqIO.write(record, f, "genbank")

#解析基因类型的marker
def parse_gene_in_txt(gene):
    gb_file_name = gene + '.gb'
    fa_file_name = gene + '.fasta'
    f = open(os.path.join(os.getcwd(),'reference',fa_file_name), 'w')
    with open(os.path.join(os.getcwd(),'raw_gb',gb_file_name), 'r') as gb_handle:
        for record in SeqIO.parse(gb_handle, 'genbank'):
            for feature in record.features:
                if feature.type == "gene":
                    gene_name = feature.qualifiers['gene'][0]
                    if gene_name.lower() == gene:
                        id = record.id
                        org = '_'.join(record.annotations['organism'].split(' '))
                        seq = feature.extract(record.seq)#自动判断在哪条链上再提取
                        f.write(f">{id}_orgnism_{org}\n{seq}\n")
    f.close()

def parse_its_in_txt():
    pass

def parse_HA_in_txt():
    pass


