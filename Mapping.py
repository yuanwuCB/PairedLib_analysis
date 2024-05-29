from Bio.Seq import Seq
import pandas as pd
import re

header = ['ID','gRNA','PAM','target_seq','edit_site_in_proto','context','clinvar_id','disease_name',
              'gene_info','chr','coordinate','depth',
              'pos_1_A_to_G','pos_1_C_to_T','pos_1_C_to_G','pos_1_C_to_A',
              'pos_2_A_to_G','pos_2_C_to_T','pos_2_C_to_G','pos_2_C_to_A',
              'pos_3_A_to_G','pos_3_C_to_T','pos_3_C_to_G','pos_3_C_to_A',
              'pos_4_A_to_G','pos_4_C_to_T','pos_4_C_to_G','pos_4_C_to_A',
              'pos_5_A_to_G','pos_5_C_to_T','pos_5_C_to_G','pos_5_C_to_A',
              'pos_6_A_to_G','pos_6_C_to_T','pos_6_C_to_G','pos_6_C_to_A',
              'pos_7_A_to_G','pos_7_C_to_T','pos_7_C_to_G','pos_7_C_to_A',
              'pos_8_A_to_G','pos_8_C_to_T','pos_8_C_to_G','pos_8_C_to_A',
              'pos_9_A_to_G','pos_9_C_to_T','pos_9_C_to_G','pos_9_C_to_A',
              'pos_10_A_to_G','pos_10_C_to_T','pos_10_C_to_G','pos_10_C_to_A',
              'pos_11_A_to_G','pos_11_C_to_T','pos_11_C_to_G','pos_11_C_to_A',
              'pos_12_A_to_G','pos_12_C_to_T','pos_12_C_to_G','pos_12_C_to_A',
              'pos_13_A_to_G','pos_13_C_to_T','pos_13_C_to_G','pos_13_C_to_A',
              'pos_14_A_to_G','pos_14_C_to_T','pos_14_C_to_G','pos_14_C_to_A',
              'pos_15_A_to_G','pos_15_C_to_T','pos_15_C_to_G','pos_15_C_to_A',
              'pos_16_A_to_G','pos_16_C_to_T','pos_16_C_to_G','pos_16_C_to_A',
              'pos_17_A_to_G','pos_17_C_to_T','pos_17_C_to_G','pos_17_C_to_A',
              'pos_18_A_to_G','pos_18_C_to_T','pos_18_C_to_G','pos_18_C_to_A',
              'pos_19_A_to_G','pos_19_C_to_T','pos_19_C_to_G','pos_19_C_to_A',
              'pos_20_A_to_G','pos_20_C_to_T','pos_20_C_to_G','pos_20_C_to_A']
PATH = '/home/yuanw/scratch-midway2/010823-nextseq/'
gRNA_backbone = 'GTTTTAGAGC'

def mapping(file_name):
    counter = [0,0,0,0,0] # 1.target matches with gRNA. 2.PAM exist
    edit_dict = {}
    gRNA_orders = pd.read_excel('/project2/weixintang/yuanwu/103023-nextseq/09082023_disease_C_lib.xlsx',sheet_name='summary_final')
    gRNA_order = gRNA_orders['gRNA'].to_list()
    gRNAs = []
    num = 0
    with open(PATH+'%s_R1_001.fastq'%file_name) as f1:
        with open(PATH+'%s_R2_001.fastq'%file_name) as f2:
            for gRNA_read,target_read in zip(f1.readlines(),f2.readlines()):
                num += 1
                if num % 4 != 2:
                    continue
                gRNA_read = gRNA_read.rstrip('\n')
                target_read = target_read.rstrip('\n')
                counter[0] += 1  # total read
                gRNA_posis = [_.start() for _ in re.finditer(str(gRNA_backbone).upper(),str(gRNA_read))]
                if len(gRNA_posis)==1:
                    counter [1] += 1 # gRNA backbone
                    gRNA = gRNA_read[gRNA_posis[0]-20:gRNA_posis[0]]
                    if gRNA in gRNA_order:
                        counter [2] += 1 # gRNA in order list
                        inform = list(gRNA_orders.loc[gRNA_orders.loc[gRNA_orders['gRNA'] == gRNA].index].values[0,:])
                    else:
                        continue
                    sequence = Seq(target_read).reverse_complement()
                    posis = [_.start() for _ in re.finditer(str(gRNA).upper().replace('A','(A|G)').replace('C','(C|T|A|G)'),str(sequence))]
                    if len(posis) == 1:
                        target = sequence[posis[0]:posis[0]+20]
                        counter[3] += 1 # gRNA match with target
                        gRNAs.append(gRNA)
                        if sequence[posis[0]+21:posis[0]+23] == 'GG':
                            counter[4] += 1 # PAM exist
                            if gRNA not in edit_dict.keys():
                                edit_dict[gRNA] = inform+[1]+[0]*80 # ID    gRNA    PAM target_seq  edit_site_in_proto  context clinvar_id  disease_name    gene_info   chr coordinate  read
                            elif gRNA in edit_dict.keys():
                                edit_dict[gRNA][11] += 1
                            for _ in range(len(gRNA)):
                                if gRNA[_] != target[_]:
                                    if gRNA[_] == 'A' and target[_] == 'G':
                                        edit_dict[gRNA][12+4*_+0] += 1
                                    elif gRNA[_] == 'C' and target[_] == 'T':
                                        edit_dict[gRNA][12+4*_+1] += 1
                                    elif gRNA[_] == 'C' and target[_] == 'G':
                                        edit_dict[gRNA][12+4*_+2] += 1
                                    elif gRNA[_] == 'C' and target[_] == 'A':
                                        edit_dict[gRNA][12+4*_+3] += 1
    pd.DataFrame(edit_dict).T.to_excel(PATH+'%s.xlsx'%file_name,index=False,header=header)
    print('\t'.join(map(str,counter))+'\n')
    return

if __name__ == '__main__':
    mapping(f'{PATH}050e_228lib_L001')
