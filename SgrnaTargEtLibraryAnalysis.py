"""
<STELA> -- Sgrna TargEt Library Analysis
Tang Lab Paired sgRNA Target Library Mapping Pacakge.
Updated on 05/06/2025 by Yuan Wu
"""
import os.path
import re
import time
import gzip
import math
import datetime
import warnings
import numpy as np
import pandas as pd
from Bio import SeqIO
from Bio.Seq import Seq
import seaborn as sns
import matplotlib.pyplot as plt
from matplotlib.colors import LinearSegmentedColormap

def fn_assemble(file_list):
    try:
        map_object = map(int,file_list.split('-'))
        file_tag = list(map_object)
    except AttributeError:
        print('error')
    fn_num_list = []
    for _ in range(0,len(file_tag),2):
        for fn_num in range(file_tag[_], file_tag[_+1]+1):
            fn_num_list.append(fn_num)
    return fn_num_list

# Hook is used to find the spacer.
Hook = 'GTTTTAGAGC' # the first 10 base of sgRNA scaffold.
# directory of raw reads
input_path = '/Users/yuanwu/Downloads/20250601_NextSeq/'
# directory of where you want the output file to
output_path = '/Users/yuanwu/Downloads/20250601_NextSeq/'
# All the barcodes of FASTQ files that you want to analyze
Barcode_list = fn_assemble('481-488-529-536-721-730-493-500-541-548-733-742-505-512-553-560-745-754')
# Mammalian cell library used here. [Library_A, Library_C]
Library_Type = 'Library_A'
Editing_Type = 'A'
# The region used for context selectivity heatmap plotting.
CONTEXT_WINDOW = [range(4,9),[6]]
REPLICATE = 3 # optional, only if you want to take the average number of your replicate samples.
DEPTH = 1
NAN_CUTOFF = 0.1

# https://matplotlib.org/stable/gallery/text_labels_and_annotations/font_family_rc.html
plt.rcParams['font.family'] = 'Arial'
# https://matplotlib.org/stable/users/explain/text/fonts.html
plt.rcParams["svg.fonttype"] = 'none'
pd.set_option('display.max_columns', None)
NOW = datetime.datetime.now()
DATE = NOW.strftime("%Y-%m-%d")
warnings.filterwarnings("ignore")

def CONTEXT_HEATMAP_AVG(File_list,WINDOW,REPLICATE):
    colors = ["white","green"]
    cmap = LinearSegmentedColormap.from_list("WhiteGreen", colors, N=256)
    fig, axs = plt.subplots(math.ceil(len(File_list)//3 /4) + 1, 4, constrained_layout=True, figsize=(8, 2*(math.ceil(len(File_list)//3 /4)+1)))
    for num,REP_GROUP in enumerate(zip(File_list[:len(File_list)//REPLICATE],
                                       File_list[len(File_list)//REPLICATE:-1*len(File_list)//REPLICATE],
                                       File_list[-1*len(File_list)//REPLICATE:])):
        EDIT_RATE_GROUP = []
        for File in REP_GROUP:
            if not os.path.exists(f'{output_path}TCDI{File}.xlsx'):
                print(f'{File} has nothing to map. Skipping...')
                continue
            try:
                (TOTAL_READS_DF, EDITED_READS_DF) = (pd.read_excel(f'{output_path}TCDI{File}.xlsx',sheet_name='TOTAL_READS',index_col=0,header=0),
                                                     pd.read_excel(f'{output_path}TCDI{File}.xlsx',sheet_name='EDITED_READS',index_col=0,header=0))
            except:
                print(f'{File} has nothing to map. Skipping...')
                continue
            EDIT_RATE_DF = EDITED_READS_DF/TOTAL_READS_DF*100
            if EDIT_RATE_DF.isna().sum().sum() > EDIT_RATE_DF.size * NAN_CUTOFF:
                print(f'{File} has more than {NAN_CUTOFF*100}% no-read context. Skipping...')
                continue
            else:
                EDIT_RATE_GROUP.append(EDIT_RATE_DF.loc[WINDOW].to_numpy())
        try:
            data = np.nanmean(EDIT_RATE_GROUP,axis=0).mean(axis=0).reshape((4,4))
            data = pd.DataFrame(data,index=['A','T','C','G'],columns=['A','T','C','G'])
            sns.heatmap(data, annot=True, fmt=".1f", cmap=cmap, ax=axs[num // 4, num % 4], linewidths=0.2,
                        linecolor='black', cbar_kws={"shrink": 0.65,'format':'%.1f'}, annot_kws={"color": "black"}, square=True)
        except:
            print(f'{REP_GROUP} has no correct reading file. Skipping...')
    for num in range(axs.size):
        if num < len(File_list)//3:
            if len(WINDOW) > 1:
                axs[num // 4, num % 4].set_title(f'{File_list[num]}_Pos_[{WINDOW[0]},{WINDOW[-1]}]')
            else:
                axs[num // 4, num % 4].set_title(f'{File_list[num]}_Pos_[{WINDOW[0]}]')
        else:
            axs[num // 4, num % 4].axis('off')
    if len(WINDOW) > 1:
        plt.savefig(f'{output_path}HEATMAP_CONTEXT_AVERAGE_{DATE}_Pos_[{WINDOW[0]},{WINDOW[-1]}].svg', dpi=300)
    else:
        plt.savefig(f'{output_path}HEATMAP_CONTEXT_AVERAGE_{DATE}_Pos_[{WINDOW[0]}].svg', dpi=300)
    plt.show()
    return

def CONTEXT_HEATMAP(File_list,WINDOW):
    colors = ["white","green"]
    cmap = LinearSegmentedColormap.from_list("WhiteGreen", colors, N=256)
    fig, axs = plt.subplots(math.ceil(len(File_list)/4) + 1, 4, constrained_layout=True, figsize=(8, 2*(math.ceil(len(File_list)/4)+1)))
    for num,File in enumerate(File_list):
        if not os.path.exists(f'{output_path}TCDI{File}.xlsx'):
            print(f'{File} has nothing to map. Skipping...')
            continue
        try:
            (TOTAL_READS_DF, EDITED_READS_DF) = (pd.read_excel(f'{output_path}TCDI{File}.xlsx',sheet_name='TOTAL_READS',index_col=0,header=0),
                                                 pd.read_excel(f'{output_path}TCDI{File}.xlsx',sheet_name='EDITED_READS',index_col=0,header=0))
        except:
            print(f'{File} has nothing to map. Skipping...')
            continue
        EDIT_RATE_DF = EDITED_READS_DF/TOTAL_READS_DF*100
        data = EDIT_RATE_DF.loc[WINDOW].mean(axis=0,skipna=True).to_numpy().reshape((4,4))
        data = pd.DataFrame(data,index=['A','T','C','G'],columns=['A','T','C','G'])
        sns.heatmap(data, annot=True, fmt=".1f", cmap=cmap, ax=axs[num // 4, num % 4], linewidths=0.2,
                    linecolor='black', cbar_kws={"shrink": 0.65,'format':'%.1f'}, annot_kws={"color": "black"}, square=True)
    for num in range(axs.size):
        if num < len(File_list):
            if len(WINDOW) > 1:
                axs[num // 4, num % 4].set_title(f'{File_list[num]}_Pos_[{WINDOW[0]},{WINDOW[-1]}]')
            else:
                axs[num // 4, num % 4].set_title(f'{File_list[num]}_Pos_[{WINDOW[0]}]')
        else:
            axs[num // 4, num % 4].axis('off')
    if len(WINDOW) > 1:
        plt.savefig(f'{output_path}HEATMAP_CONTEXT_{DATE}_Pos_[{WINDOW[0]},{WINDOW[-1]}].svg', dpi=300)
    else:
        plt.savefig(f'{output_path}HEATMAP_CONTEXT_{DATE}_Pos_[{WINDOW[0]}].svg', dpi=300)
    plt.show()
    return

def ANALYZER(MAPPING_DATA,FILE,Editing_Type,OutputExcel,LOG):
    global CONTEXT_WINDOW, DEPTH
    # Configuration
    INDEX = {'A': 0, 'T': 1, 'C': 2, 'G': 3}
    CONVERSION = {'C': 'T', 'A': 'G'}
    Context_Header = []
    BASE = Editing_Type
    for Five_Prime_Base in ['A','T','C','G']:
        for Three_Prime_Base in ['A','T','C','G']:
            Context_Header.append(f'{Five_Prime_Base}{BASE}{Three_Prime_Base}')

    # Start calculating the editing rate.
    OUTPUT_TOTAL_df,OUTPUT_EDITED_df = (pd.DataFrame(0,columns=Context_Header,index=range(1,21,1)),
                                        pd.DataFrame(0,columns=Context_Header,index=range(1,21,1)))
    for NUM,ROW in MAPPING_DATA.iterrows():
        DATA = ROW.to_dict()
        COVERAGE = DATA['depth']
        SPACER = DATA['gRNA']
        FULL_TARGET = DATA['target_seq']
        for POS in range(1,21,1):
            if SPACER[POS-1] == BASE:
                CONTEXT = FULL_TARGET[8 + POS:11 + POS].upper()
                OUTPUT_TOTAL_df.loc[POS,CONTEXT] += COVERAGE
                OUTPUT_EDITED_df.loc[POS,CONTEXT] += DATA[f'pos_{POS}_{BASE}_to_{CONVERSION[BASE]}']
        # statistical number
        AVERAGE_COVERAGE = MAPPING_DATA['depth'][MAPPING_DATA['depth'] >= DEPTH].mean()
        MEDIAN_COVERAGE = MAPPING_DATA['depth'][MAPPING_DATA['depth'] >= DEPTH].median()
    LOG.writelines(f"TCDI{FILE}\t{len(MAPPING_DATA['depth'][MAPPING_DATA['depth'] >= DEPTH])}\t{AVERAGE_COVERAGE}\t{MEDIAN_COVERAGE}\n")
    OUTPUT_TOTAL_df.to_excel(OutputExcel,sheet_name='TOTAL_READS')
    OUTPUT_EDITED_df.to_excel(OutputExcel, sheet_name='EDITED_READS')
    return

def ALIGNMENT(R1_FQ, R2_FQ, Library_Type, OutputExcel, File, LOG):
    global input_path,output_path,Hook
    COUNTER = {"Total_Reads":0,'gRNA_Scaffold_Found':0,'Protospacer_Match':0,'Spacer_Protospacer_Match':0,'Correct_PAM':0}

    # Reference Excel which contains all the information of library sites.
    if Library_Type == 'Library_A':
        REFERENCE = pd.read_excel(
            '/Users/yuanwu/YuanWu_Tang_lab/publication/TadA8r/Revision/oligo pool design/0_YW_oligo_design.xlsx',
            sheet_name='summary_final')
    elif Library_Type == 'Library_C':
        REFERENCE = pd.read_excel(
            '/Users/yuanwu/YuanWu_Tang_lab/publication/Topipotent adenine base editor/C_disease_lib/09082023_disease_C_lib.xlsx',
            sheet_name='summary_final')

    # CLUSTER all the READ base on the SPACER first.
    READS_DICT = {}
    for (READ1, READ2) in zip(SeqIO.parse(gzip.open(R1_FQ, 'rt'), "fastq"),
                              SeqIO.parse(gzip.open(R2_FQ, 'rt'), "fastq")):
        COUNTER['Total_Reads'] += 1
        SPACER_READ,PROTOSPACER_READ = str(READ1.seq), str(READ2.seq)
        SPACER_POS = [_.start() for _ in re.finditer(str(Hook).upper(), SPACER_READ)]
        if len(SPACER_POS) == 1:
            COUNTER['gRNA_Scaffold_Found'] += 1
            SPACER = SPACER_READ[SPACER_POS[0] - 20:SPACER_POS[0]]
            if SPACER in READS_DICT.keys():
                READS_DICT[SPACER].append((READ1, READ2))
            elif SPACER not in READS_DICT.keys():
                READS_DICT[SPACER] = [(READ1, READ2)]

    # MATCH Reads to gRNA reference list one by one.
    EDIT_DICT = {}
    SPACER_WHITE_LIST = REFERENCE['gRNA'].to_list()
    INDEX = {'AG':0,'CT':1,'CG':2,'CA':3}
    for SPACER in READS_DICT.keys():
        if SPACER in SPACER_WHITE_LIST:
            COUNTER['Protospacer_Match'] += len(READS_DICT[SPACER])  # gRNA in order list
            CLINVAR_RECORD = list(REFERENCE.loc[REFERENCE.loc[REFERENCE['gRNA'] == SPACER].index].values[0, :])
        else:
            continue
        for (SPACER_READ,PROTOSPACER_READ) in READS_DICT[SPACER]:
            SPACER_READ,PROTOSPACER_READ = SPACER_READ.seq,Seq(PROTOSPACER_READ.seq).reverse_complement().upper()
            PROTOSPACER_POS = [_.start() for _ in
                     re.finditer(str(SPACER).upper().replace('A', '(A|G)').replace('C', '(C|T|A|G)'),
                                 str(PROTOSPACER_READ))]
            if len(PROTOSPACER_POS) != 1:
                continue
            PROTOSPACER = PROTOSPACER_READ[PROTOSPACER_POS[0]:PROTOSPACER_POS[0] + 20]
            COUNTER['Spacer_Protospacer_Match'] += 1  # gRNA match with target
            if PROTOSPACER_READ[PROTOSPACER_POS[0] + 21:PROTOSPACER_POS[0] + 23] == 'GG':
                COUNTER['Correct_PAM'] += 1  # PAM exist
                if SPACER not in EDIT_DICT.keys():
                    # ID    gRNA    PAM target_seq  edit_site_in_proto  context clinvar_id  disease_name    gene_info   chr coordinate  read
                    EDIT_DICT[SPACER] = (CLINVAR_RECORD + [1] +
                                         [0] * 80)
                elif SPACER in EDIT_DICT.keys():
                    EDIT_DICT[SPACER][11] += 1
                for POS in range(len(SPACER)):
                    if PROTOSPACER[POS] != SPACER[POS]:
                        EDIT_DICT[SPACER][len(CLINVAR_RECORD)+1+20*INDEX[f'{SPACER[POS]}{PROTOSPACER[POS]}']+POS] += 1
    header = (['ID', 'gRNA', 'PAM', 'target_seq', 'edit_site_in_proto', 'context', 'clinvar_id', 'disease_name',
              'gene_info', 'chr', 'coordinate', 'depth'] +
              [f'pos_{pos}_A_to_G' for pos in range(1,21,1)] +
              [f'pos_{pos}_C_to_T' for pos in range(1,21,1)] +
              [f'pos_{pos}_C_to_G' for pos in range(1,21,1)] +
              [f'pos_{pos}_C_to_A' for pos in range(1,21,1)])
    # generate and save mapping output of each sample.
    MAPPING_df = pd.DataFrame(EDIT_DICT).T
    try:
        MAPPING_df.to_excel(OutputExcel, sheet_name='Mapping', index=False, header=header)
        MAPPING_df.columns = header
    except ValueError:
        print(f'{File} has nothing to map. Skipping...')
        return None
    # Output the filter process log.
    COUNTER_df = pd.DataFrame(COUNTER, index=['Counts'])
    print(COUNTER_df)
    LOG.writelines(COUNTER_df.to_string()+'\n')
    return MAPPING_df

def main():
    def TIME_PRINT(string):
        print(string)
        return time.time()

    global REFERENCE, input_path, output_path, Excelwriter, Editing_Type
    T1 = time.time()
    MAPLOG, ANALYZERLOG = (open(f'{output_path}MAPPING_{DATE}.log', 'w'),
                           open(f'{output_path}ANALYZING_{DATE}.log', 'w'))
    ANALYZERLOG.writelines(f'File\tCover_sites\tAVERAGE_Coverage\tMedian_Coverage\n') # Header of Analytical LOG file.
    # Start alignment.
    for file in Barcode_list:
        if os.path.isfile(f'{output_path}TCDI{file}.xlsx'):
            print(f'Mapping result of {file} already exist. Skipping...')
            continue
        Excelwriter = pd.ExcelWriter(f'{output_path}TCDI{file}.xlsx')
        print(f'Starting mapping of {file}...')
        MAPPING_df = ALIGNMENT(f'{input_path}TCDI{file}_S{file}_L001_R1_001.fastq.gz',
                               f'{input_path}TCDI{file}_S{file}_L001_R2_001.fastq.gz',
                               Library_Type,
                               Excelwriter,
                               file,
                               MAPLOG)
        if MAPPING_df is None: continue
        T1 = TIME_PRINT(f'Finished mapping of {file}...Using time: {(time.time()-T1):.2f} s')
        # Start editing analyzing.
        print(f'Starting analyzing of {file}...')
        ANALYZER(MAPPING_df,
                 file,
                 Editing_Type,
                 Excelwriter,
                 ANALYZERLOG)
        T1 = TIME_PRINT(f'Finished analyzing of {file}...Using time: {(time.time()-T1):.2f} s')
        Excelwriter.close()
    MAPLOG.close()
    ANALYZERLOG.close()
    # Plotting the 4x4 editing rate heatmap of every single file.
    for WINDOW in CONTEXT_WINDOW:
        print(f'Starting plotting of individual heatmap...')
        CONTEXT_HEATMAP(Barcode_list,
                        WINDOW)
        T1 = TIME_PRINT(f'Finish plotting of individual heatmap ...Using time: {(time.time() - T1):.2f} s')
    # Plotting the 4x4 editing rate heatmap of average editing rate from the replicate results.
    if REPLICATE != None and len(Barcode_list) % REPLICATE != 0:
        raise IndexError(f'The file number ({len(Barcode_list)}) from the input did not contain all your replicates.')
    else:
        for WINDOW in CONTEXT_WINDOW:
            print(f'Starting plotting of average heatmap...')
            CONTEXT_HEATMAP_AVG(Barcode_list,
                                WINDOW,
                                REPLICATE)
            T1 = TIME_PRINT(f'Finish plotting of average heatmap ...Using time: {(time.time() - T1):.2f} s')
    return

if __name__ == '__main__':
    main()