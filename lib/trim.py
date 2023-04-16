import numpy as np
import os
from Bio import SeqIO

def get_array(file_path):#读取为矩阵
    seq_lis = []
    ID_lis = []

    for record in SeqIO.parse(file_path,'fasta'):
        seq_lis.append(list(str(record.seq)))
        ID_lis.append(record.description)
    seq_arr = np.array(seq_lis,dtype='U1')
    return seq_arr, ID_lis

def gap_score(seq_arr):#计算gap score Sg
    row_n = len(seq_arr[:,0])
    col_n = len(seq_arr[0,:])
    Sg_tup = ()
    for col in range(col_n):
        gap = 0
        for row in range(row_n):
            if seq_arr[row,col] == '-':
                gap += 1
        Sg = 1-(gap/row_n)
        Sg_tup = (*Sg_tup,Sg)
    return Sg_tup

def real_simi_score():
    pass

def simi_score(seq_arr):#计算similarity score Ss
    row_n = len(seq_arr[:,0])
    col_n = len(seq_arr[0,:])
    Ss_tuple = ()
    for col in range(col_n):
        base_dic = {'a':0,'t':0,'c':0,'g':0}
        for row in range(row_n):
            residue = seq_arr[row,col].lower()
            if residue in base_dic:
                base_dic[residue] += 1
        Ss = max(base_dic.values())/row_n
        Ss_tuple = (*Ss_tuple,Ss)
    return Ss_tuple

def trim(seq_arr,Sg_tup,Ss_tup,Sg_t,Ss_t):#根据Sg和Ss裁剪
    row_n = len(seq_arr[:,0])
    col_n = len(seq_arr[0,:])
    remain_tup = ()
    for col in range(col_n):
        if Sg_tup[col] > Sg_t and Ss_tup[col] > Ss_t:
            remain_tup = (*remain_tup,col)
    trimmed_arr = seq_arr[:,remain_tup]
    return trimmed_arr

def out_put(trimmed_arr,ID_lis,out_path):
    row_n = len(trimmed_arr[:,0])
    with open(out_path,'w') as out_handle:
        for row in range(row_n):
            out_handle.write(f">{ID_lis[row]}\n{''.join(trimmed_arr[row])}\n\n")
        
def my_trim(Sg_t,Ss_t):
    # in_path = os.path.join(os.path.dirname(os.path.abspath(__file__)),'aligned_data.fasta')
    in_path = os.path.join(os.getcwd(),'aligned_data.fasta')
    # out_path = os.path.join(os.path.dirname(os.path.abspath(__file__)),'out.fasta')
    out_path = os.path.join(os.getcwd(),'trimmed_data.fasta')
    seq_arr, ID_lis = get_array(in_path)
    Sg_tup = gap_score(seq_arr)
    Ss_tup = simi_score(seq_arr)
    trimmed_arr = trim(seq_arr,Sg_tup,Ss_tup,Sg_t,Ss_t)
    out_put(trimmed_arr,ID_lis,out_path)



        
