from Bio import SeqIO
import sys
import os
import time

def is_gb(file_path):
    try:
        with open(file_path,'r') as gb_handle:
            records = SeqIO.parse(gb_handle,'genbank')
            return any(records)
    except:
        return False

def is_fasta(file_path):
    try:
        with open(file_path,'r') as gb_handle:
            records = SeqIO.parse(gb_handle,'fasta')
            return any(records)
    except:
        return False


#检查当前工作目录下的reference文件夹内容是否符合我们的要求
def check_rf_marker(marker_name):
    std_extension = ('.fa','.fasta','.fna','ffn','faa')
    if not os.path.exists('reference'):
        sys.stdout.write("You don't have the reference database!\n")#检查工作目录是否有reference目录
        exit(1)
    else:
        reference_files = os.listdir('reference')
        flag = 0
        if reference_files:
            for file in reference_files:
                filename, extension = os.path.splitext(os.path.basename(file))
                if filename.lower() == marker_name and extension in std_extension:
                    file_path = os.path.join(os.getcwd(),'reference',file)
                    flag = 1
                    break
            if flag:
                if is_fasta(file_path):
                    for i in range(3):
                        sys.stdout.write('\t  .\n')
                        time.sleep(0.5)
                    sys.stdout.write("reference check complete!\n")
                else:
                    sys.stdout.write("we can't parse your reference file,please make sure your file in fasta format!\n")#检查是否真的是fasta文件
            else:
                sys.stdout.write(f"we don't find {marker_name}.fasta in reference dir!\n")#检查reference文件夹是否包含目标marker的fasta文件
        else:
            sys.stdout.write("The reference dir is empty!\n")#检查reference文件夹是否为空
            exit(1)

#检查当前工作目录下的raw_gb文件夹内容是否符合我们的要求
def check_raw_gb(marker_name):
    std_extension =('gb','gbk')
    if not os.path.exists('raw_gb'):
        sys.stdout.write("You don't have the raw_gb !\n")#检查工作目录是否有raw_gb目录
        exit(1)
    else:
        gb_files = os.listdir('raw_gb')
        flag =0
        if gb_files:
            for file in gb_files:
                filename, extension = os.path.splitext(os.path.basename(file))
                if filename.lower() == marker_name and extension in std_extension:
                    file_path = os.path.join(os.getcwd(),'raw_gb',file)
                    flag = 1
                    break
            if flag:
                if is_gb(file_path):
                    pass
                else:
                    sys.stdout.write("we can't parse your reference file,please make sure your file in genbank format!\n")#检查是否真的是gb文件
                    exit(1)
            else:
                sys.stdout.write(f"we don't find {marker_name}.gb in raw_gb dir!\n")#检查raw_gb文件夹是否包含目标marker的genbank文件
                exit(1)
        else:
            sys.stdout.write("The raw_gb dir is empty!\n")#检查raw_gb文件夹是否为空
            exit(1)

#检查当前工作目录下的query文件夹内的目标marker
def check_query_marker(marker_name):
    std_extension = ('.fa','.fasta','.fna','ffn','faa')
    if not os.path.exists('query'):
        sys.stdout.write("You don't have the query directory !\n")#检查工作目录是否有query目录
        exit(1)
    else:
        query_files = os.listdir('query')
        flag =0
        if query_files:
            for file in query_files:
                filename, extension = os.path.splitext(os.path.basename(file))
                if filename.lower() == marker_name and extension in std_extension:
                    file_path = os.path.join(os.getcwd(),'query',file)
                    flag = 1
                    break
            if flag:
                if is_fasta(file_path):
                    for i in range(3):
                        sys.stdout.write('\t  .\n')
                        time.sleep(0.5)
                    sys.stdout.write("query check complete!\n")
                    for i in range(3):
                        sys.stdout.write('\t  .\n')
                        time.sleep(0.5)
                else:
                    sys.stdout.write("we can't parse your query file,please make sure your file in fasta format!\n")#检查是否真的是gb文件
                    exit(1)
            else:
                sys.stdout.write(f"we don't find {marker_name}.fasta in query dir!\n")#检查raw_gb文件夹是否包含目标marker的genbank文件
                exit(1)
        else:
            sys.stdout.write("The query dir is empty!\n")#检查raw_gb文件夹是否为空
            exit(1)

            

