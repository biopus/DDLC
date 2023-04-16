from filter import * 
from collections import deque
from operator import itemgetter
import gc
from collections import Counter

pars.add_argument('-ta', metavar='<int>', type=int, help='''thread of assemble''', required=False, default= 1)
pars.add_argument('-ka', metavar='<int>', type=int, help='''kmer of assemble''',  default=41)
pars.add_argument('-limit_count', metavar='<int>', type=int, help='''limit of kmer count''', required=False, default= 0)

# 简并碱基字典
D_BASE_DICT = {'AG':'R','CT':'Y', 'GT':'K', 'GC':'S','AC':'M', 'AT':'W','GA':'R','TC':'Y','TG':'K', 'CG':'S','CA':'M', 'TA':'W',}

def consensus_sequence(multi_alignment):
    """
    根据排序后的序列生成一致序列
    """ 
    num_sequences = len(multi_alignment)
    seq_length = len(multi_alignment[0])
    consensus = []
    for column in range(seq_length):
        base_count = Counter(multi_alignment[row][column] for row in range(1, num_sequences) if multi_alignment[row][column] != '-')
        if base_count:
            consensus.append(base_count.most_common(1)[0][0])
        else:
            consensus.append('?')
    return ''.join(consensus)

def Reverse_Complement_ACGT(seq):
    """
    简化版反向互补
    """ 
    return seq.translate(str.maketrans('ACGT', 'TGCA'))[::-1]

def Make_Assemble_Dict(file_list, kmer_size, _kmer_dict, _ref_dict):
    """
    构建拼接用的字典
    :param file_list: 文件列表
    :param kmer_size: kmer的长度
    :param _kmer_dict: 待生成的字典value的格式为[深度，位置（1000以内的整数）]
    :param _ref_dict: 参考序列的字典
    :return: 返回kmer的总数量
    """
    MASK_BIN = (1<< (kmer_size<<1)) - 1 # kmer的掩码
    kmer_count = 0
    for file in file_list:
        infile = open(file, 'r', encoding='utf-8', errors='ignore')
        infile.readline()
        for line in infile:
            temp_str = [] # 为了支持多行的fasta文件作为源数据
            if Filted_File_Ext == '.fasta':
                while line and line[0] != '>':
                    temp_str.append(line)
                    line = infile.readline()
            else:
                temp_str.append(line)
                infile.readline()
                infile.readline()
                infile.readline()
            read_seq = ''.join(filter(str.isalpha, ''.join(temp_str).upper()))
            intseqs, read_len = Seq_To_Int(read_seq) # 序列转整数，获取长度
            intseqs.append(Seq_To_Int(read_seq, True)) # 加入反向互补序列
            for x in intseqs:
                kmer_count += read_len - kmer_size + 1
                for j in range(0, read_len - kmer_size):
                    temp_list, temp_pos,  kmer = [], 0, x >> (j<<1) & MASK_BIN
                    if kmer in _kmer_dict:
                        _kmer_dict[kmer][0] += 1
                    else:
                        if kmer in _ref_dict: # kmer的位置
                            if  _ref_dict[kmer] & 1073741824: # 判断是否为反向互补的序列
                                temp_pos = 1000 - (_ref_dict[kmer] & 1023)
                                temp_list = [1, temp_pos, 1] # 标记为反向的的kmer
                            else: 
                                temp_pos = _ref_dict[kmer] & 1023
                                temp_list = [1, temp_pos, 0] 
                        else:
                            temp_list = [1, 1023, 1] 
                        _kmer_dict[kmer] = temp_list
        infile.close()
    return kmer_count

def Get_Closest_Power_Of_Two(num):
    """
    获取小于正整数的最接近的2的整数次方
    """
    power = 1
    while power <= num:
        power <<= 1
    return power >> 1

def Get_Middle_Fragment(text, slice_len):
    """
    截取reads中间的高质量片段
    """
    start = (len(text) - slice_len) >> 1
    end = start + 128
    return text[start:end]

def Make_Reads_Dict(file_list, _reads_dict):
    """
    截取reads中间的片段，构建高质量的reads字典
    :param file_list: 文件列表
    :param _reads_dict: 待生成的字典value的格式为[seq, count]
    :return: 返回切片的长度
    """
    read_len = 0
    slice_len = 0
    for file in file_list:
        infile = open(file, 'r', encoding='utf-8', errors='ignore')
        infile.readline()
        for line in infile:
            temp_str = []
            if Filted_File_Ext == '.fasta':
                while line and line[0] != '>':
                    temp_str.append(line)
                    line = infile.readline()
            else:
                temp_str.append(line)
                infile.readline()
                infile.readline()
                infile.readline()
            read_seq = ''.join(filter(str.isalpha, ''.join(temp_str).upper()))
            if not read_len:
                read_len = len(read_seq)
                slice_len = Get_Closest_Power_Of_Two(read_len)
            intseqs = [read_seq, Reverse_Complement_ACGT(read_seq)] # 加入反向互补序列
            for x in intseqs:
                slice_reads = Get_Middle_Fragment(x, slice_len)
                if slice_reads in _reads_dict:
                    _reads_dict[slice_reads] += 1
                else:
                    _reads_dict[slice_reads] = 1
        infile.close()
    return slice_len

def Get_Similar_Ref(file_list, kmer_size, _reads_set_combine):
    """
    构建拼接用的字典
    :param file_list: 参考序列文件列表
    :param kmer_size: kmer的长度
    :param _reads_set_combine: reads的集合
    :return: 返回最近源的参考序列
    """
    MASK_BIN = (1<< (kmer_size<<1)) - 1 # kmer的掩码
    Similar_Ref, ref_kmer_count = ('',0,''),  0
    for file in file_list:
        infile = open(file, 'r', encoding='utf-8', errors='ignore')
        # 第一条参考序列的名字
        ref_name = infile.readline()
        for line in infile:
            # 读取每一条参考序列
            temp_str = []
            while line and line[0] != '>':
                temp_str.append(line)
                line = infile.readline()
            
            ref_seq = ''.join(filter(str.isalpha, ''.join(temp_str).upper()))
            intseqs, read_len = Seq_To_Int(ref_seq) # 序列转整数，获取长度
            intseqs.append(Seq_To_Int(ref_seq, True)) # 加入反向互补序列
            # 参考序列的kmer在集合中的计数
            kmer_count = 0
            for x in intseqs:
                for j in range(0, read_len-kmer_size):
                    kmer = x >> (j<<1) & MASK_BIN
                    if kmer in _reads_set_combine:
                        kmer_count += 1
            if ref_kmer_count < kmer_count:
                 ref_kmer_count = kmer_count
                 Similar_Ref = (intseqs[0], read_len, ref_name)
            # 参考序列的名字
            if line[0] == '>': ref_name = line
        infile.close()
    return Similar_Ref[2], Int_To_Seq(Similar_Ref[0], Similar_Ref[1])


def Median(x):
    """
    使用中位数分割列表
    :return: 左半边列表，右半边列表，中位数
    """ 
    x = sorted(x)
    length = len(x)
    mid, rem = divmod(length, 2)    # divmod函数返回商和余数
    if rem:
        return x[:mid], x[mid+1:], x[mid]
    else:
        return x[:mid], x[mid:], (x[mid-1]+x[mid])/2

def Quartile(x):
    """
    获取列表的四分位数
    :return: 左四分位数，中位数，右四分位数，最大值加1
    """ 
    lHalf, rHalf, q2 = Median(x)
    return Median(lHalf)[2], q2, Median(rHalf)[2], max(x) + 1


def Forward_Bin(seq_int, mask):
    """
    正向的迭代器
    """ 
    for x in (0,1,2,3): yield ((seq_int & mask) << 2) + x

def Get_Weight(_pos, new_pos, weight = 4):
    """
    距离和权重的关系模型，默认权重值为16，最高权重值为256，最低为0
    :param _pos: 当前kmer在参考序列的位置
    :param new_pos: 新的kmer位置
    :param weight: 默认权重
    :return: 返回计算后的权重
    """ 
    return int.bit_length((1024 - abs(_pos - new_pos)) >> 2) if (_pos and new_pos) else weight


def Get_Forward_Contig_v6(_dict, seed, kmer_size, iteration = 1024, weight = 4):
    """
    带权重的DBG贪婪拼接
    :param _pos: 当前kmer在参考序列的位置
    :param seed: 新的kmer位置
    :param kmer_size: kmer的大小
    :param iteration: 最大循环数量
    :param weight: 默认权重
    :return: best_seq, kmer_list, best_kmc, best_pos, best_snp
    """ 
    temp_list, kmer_list, stack_list, pos_list = deque([seed]), deque([seed]), deque(), deque()
    cur_kmc, cur_seq, contigs = deque(), deque(), deque()
    _pos, node_distance, best_kmc_sum  = 0, 0, 0
    MASK = (1 << ((kmer_size << 1) - 2)) - 1
    while True and iteration:
        node = [[i, _dict[i][1], _dict[i][0] << Get_Weight(_pos, _dict[i][1], weight), ACGT_DICT[i & 3]] for i in Forward_Bin(temp_list[-1], MASK) if i in _dict]
        node.sort(key = itemgetter(2), reverse=True)
        while node and node[0][0] in temp_list: node.pop(0) 
        if not node: 
            iteration -= 1
            cur_kmc_sum = sum(cur_kmc)
            contigs.append((cur_kmc_sum, cur_seq.copy()))
            if cur_kmc_sum > best_kmc_sum:
                best_kmc_sum = cur_kmc_sum
            for _ in range(node_distance):
                temp_list.pop()
                cur_kmc.pop()
                cur_seq.pop()
            if stack_list: 
                node, node_distance, _pos = stack_list.pop()
            else: 
                break
        if len(node) >= 2:
            stack_list.append((node[1:], node_distance, _pos))
            node_distance = 0
        if node[0][1] > 0: _pos = node[0][1]
        temp_list.append(node[0][0])
        kmer_list.append(node[0][0])
        pos_list.append(node[0][1])
        cur_kmc.append(node[0][2])
        cur_seq.append(node[0][0]&3)
        node_distance += 1
    return contigs, kmer_list, pos_list, best_kmc_sum

def Process_Contigs(contigs, max_weight, slice_len, reads_dict):
    """
    通过将contigs与reads进行map，来检测contig的可靠性
    :param contigs: 拼接过程获取的contigs
    :param max_weight: 最大的权重，只考虑大于最大权重一半的contigs
    :param slice_len: reads的高质量切片的长度
    :param reads_dict: reads的高质量切片的词典
    :return: 按照map上的reads的数量倒序排序过后的contigs
    """ 
    processed_contigs = sorted([[''.join(ACGT_DICT[k] for k in x[1]), x[0], 0] for x in contigs if x[0] > max_weight >> 1], key=itemgetter(1), reverse=True)
    for x in processed_contigs:
        contig_len = len(x[0])
        for j in range(contig_len - slice_len):
            if contig_len - slice_len - j >= 0:
                slice_str = x[0][contig_len - slice_len - j:contig_len - j]
                if slice_str in reads_dict:
                    x[2] += reads_dict[slice_str]
    processed_contigs.sort(key=itemgetter(2), reverse=True)
    return processed_contigs

def Get_Contig_v5(_reads_dict, slice_len, _dict, seed, kmer_size, iteration = 1024, weight = 4):
    """
    获取最优的contig
    :param _reads_dict: reads的高质量切片的词典
    :param slice_len: reads的高质量切片的长度
    :param _dict: 用于拼接的kmer字典
    :param seed: 拼接种子
    :param kmer_size: kmer的长度
    :param iteration: 构建contig时允许的最大路径分支数
    :param weight: 没有ref时的默认权重
    :return: 最优contigs，用到所有的kmer的集合，contig的大概位置，contig的权重，map到contig上的reads数量
    """ 

    contigs_1, kmer_list_1, pos_list_1, weight_1 = Get_Forward_Contig_v6(_dict, seed, kmer_size, iteration, weight)
    contigs_2, kmer_list_2, pos_list_2, weight_2 = Get_Forward_Contig_v6(_dict, Reverse_Int(seed, kmer_size), kmer_size, iteration, weight)
    # 清理位置列表
    pos_list = [x for x in pos_list_1+ pos_list_2 if x > 0 and x < 1000]
    # 获取位置中位数
    contig_pos = int(Quartile(pos_list)[1] if len(pos_list)>1 else -1)
    # 获取最可能的两侧的contig
    contigs_1_16 = Process_Contigs(contigs_1, weight_1, slice_len, _reads_dict)
    contigs_2_16 = Process_Contigs(contigs_2, weight_2, slice_len, _reads_dict)
    # 组合contig
    contig = (Reverse_Complement_ACGT(contigs_2_16[0][0]) if contigs_2_16 else '') + Int_To_Seq(seed, kmer_size) + (contigs_1_16[0][0] if contigs_1_16 else '')
    contig_weight = contigs_1_16[0][1] if contigs_1_16 else 0 + contigs_2_16[0][1] if contigs_2_16 else 0
    contig_count = contigs_1_16[0][2] if contigs_1_16 else 0 + contigs_2_16[0][2] if contigs_2_16 else 0
    return contig, set(kmer_list_1 + kmer_list_2), contig_pos, contig_weight, contig_count

if __name__ == '__main__':
    args = pars.parse_args()
    # 初始化文件夹
    if args.ta > 0:
        if not os.path.isdir(os.path.join(args.o, 'contigs')):
            os.mkdir(os.path.join(args.o, 'contigs'))
    time_stamp = args.ts if args.ta < 0 else str(time.time())

    # 通过多次运行脚本进行多进程，可改进为程序内的多进程
    for i in range(1, args.ta):
        command = [args.py,"assembler.py","-r",args.r,"-ka",args.ka,"-o",args.o,"-ta",-i,
        "-change_seed",args.change_seed,
        "-limit_count",args.limit_count,
        "-ts",time_stamp]
        _thread_subproc = threading.Thread(target = Execute_Subproc, args = (list(map(str,command)),))
        _thread_subproc.start()

    # 载入参考序列信息
    if not Get_Ref_Info(args.r, ref_path_list, ref_path_dict, ref_length_dict):
        print('Invaild reference!')
        sys.exit(0)
    t0 = time.time()

    if args.ta > 0: print('======================== Assemble =========================')
    # 定义默认权重值
    cur_weight, assemble_count, failed_count  = 4, 0, 0
    for key, value in ref_length_dict.items(): 
        limit = args.limit_count
        depth = 0
        assemble_count += 1
        contig_path = os.path.join(args.o, "contigs",key + ".contig.fasta")
        if os.path.isfile(contig_path) == False:
            with open(contig_path, 'w') as out: pass
            write_contig = False
            print('Assembling', key, assemble_count,'/',len(ref_length_dict)," " * 32)
            # 检查是哪种扩展名
            file_extensions = ['.fasta', '.fq']
            Filted_File_Ext = '.fq'
            filtered_file_path = None
            for ext in file_extensions:
                file_path = os.path.join(args.o,'filtered', key + ext)
                if os.path.exists(file_path):
                    filtered_file_path = file_path
                    Filted_File_Ext = ext
                    break
            # 清理文件
            if os.path.isfile(os.path.join(args.o, 'filtered',key + Filted_File_Ext )) == False:
                if os.path.isfile(contig_path): os.remove(contig_path)
                continue
            # 获取种子列表
            ref_dict, filtered_dict = {}, {}
            # 如果k大于参考序列长度的0.75，则减小k
            current_ka = args.ka 
            if args.ka > (value * 0.75):
                current_ka = int(value * 0.75)
                print(key, "change vaule of k to", current_ka)
            # 制作参考序列的kmer字典
            Make_Kmer_Dict_v6(ref_dict, [ref_path_dict[key]], current_ka, True, True)
            # 制作用于拼接的kmer字典
            Make_Assemble_Dict([filtered_file_path], current_ka, filtered_dict, ref_dict)
            # 缩减filtered_dict，保留大于limit和有位置信息的
            if limit > 0:
                filtered_dict = {k: v for k, v in filtered_dict.items() if v[0] > limit or v[1] > 0}
            # 必须有filtered_dict才能继续
            if not filtered_dict:
                failed_count += 1
                if os.path.isfile(contig_path): os.remove(contig_path)
                print('Could not get enough reads from filter. Estimate depth:', depth ,' '*16)
                continue
            # 处理ref_dict，标记不在filtered_dict中的kmer
            for i in ref_dict:
                if i not in filtered_dict:
                    ref_dict[i] |= 1023  #前10位全部置1，用来代表没有位置信息
            # 如果用ref_dict做seed_list，主要考虑参考序列的保守区
            # 用filtered_dict做seed_list，主要考虑测序的高丰度区
            # 长度位置在10~990之间，与参考序列方向一致v[2] == 0
            seed_list = [(k, v[0], v[1]) for k, v in filtered_dict.items() if v[1]>10 and v[1]<990 and not v[2]]
            list.sort(seed_list, key=itemgetter(1), reverse=True)
            # 必须有seed_list, 否则意味着跟参考序列差别过大
            if not seed_list:
                failed_count += 1
                if os.path.isfile(contig_path): os.remove(contig_path)
                print('Could not get enough seeds. Estimate depth:', depth ," "*16)
                continue
            # 获取seed集合，用来加速集合操作
            seed_list_len = len(seed_list)
            seed_set = set([i[0] for i in seed_list])
            # 获取contigs
            contigs = []
            reads_dict = {}
            # 获取最大切片长度，建立reads切片字典
            slice_len = Make_Reads_Dict([filtered_file_path], reads_dict)
            # 获取contigs
            while len(seed_list) > seed_list_len >> 1: # 该子集耗费了大于一半的seed就没必要再做了 
                contig, kmer_set, contig_pos, contig_weight, reads_count = Get_Contig_v5(reads_dict, slice_len, filtered_dict, seed_list[0][0], current_ka, weight = cur_weight)
                seed_list = [item for item in seed_list if (item[0] not in kmer_set) and (Reverse_Int(item[0], current_ka) not in kmer_set)]
                if reads_count * slice_len > len(contig): # 起码要有reads高质量切片能够覆盖contig，否则就是错误的拼接
                    contigs.append([contig, len(seed_set & kmer_set), contig_pos, reads_count, contig_weight])
            if not contigs:
                failed_count += 1
                if os.path.isfile(contig_path): os.remove(contig_path)
                print("No valid contig.")
                continue
            # 获取contigs中seed个数的最大值
            max_seed_count = max([i[1] for i in contigs])
            # 移除小于max_seed_count一半的contig
            if len(contigs) > 1:
                contigs = [i for i in contigs if i[1] > max_seed_count >> 1]
            if len(contigs) > 1:
                print(key, 'have mutiply contigs', ' '*10)
            # 保存contigs，格式为 长度_使用的种子数_粗略的千分比位置_跟参考序列重叠的切片数_权重
            with open(contig_path, 'w') as out:
                for x in contigs:
                    out.write('>contig_' + str(len(x[0])) + '_' + str(x[1]) + '_' + str(x[2]) + '_' + str(x[3]) + '_' + str(x[4]) + '\n')
                    out.write(x[0] + '\n')
            ref_dict, filtered_dict = {}, {}
            gc.collect()
    t1 = time.time()
    print('\nTime cost:', t1 - t0, " "*32,'\n') # 拼接所用的时间
    # 合并多进程和结果
    if args.ta > 0:
        sum_failed_count = failed_count
        wait_results = True 
        while wait_results:
            wait_results = False
            for i in range(1, args.ta):
                if os.path.isfile(time_stamp +str(-i)) == False:
                    wait_results = True
                    time.sleep(1)
                    break
        for i in range(1, args.ta):
            with open(time_stamp + str(-i) ,"r" ) as infile:
                sum_failed_count += int(infile.readline())
            os.remove(time_stamp +str(-i))
        print("Assemble Summary:", sum_failed_count, '/', assemble_count, 'failed', " "*32)
    else:
        with open(time_stamp + str(args.ta) ,"w" ) as outfile:
            outfile.write(str(failed_count))
    