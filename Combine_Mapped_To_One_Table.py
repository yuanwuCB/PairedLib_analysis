import pandas as pd
from _Final_Name_List import *
import numpy as np

# CONSTANT PARAMETER
path = '/Users/yuanwu/BaseSpace/010824-nextseq/excel_output/' # OUTPUT FOLDER OF MAPPING.PY
ot_path = '/Users/yuanwu/BaseSpace/010824-nextseq/041424-aio_table/' # OUTPUT FOLDER OF COMBINED TABLE FILE
depth = 20 # COVERAGE DEPTH FILTER

# FILE LIST. Can also use "[file for file in os.listdir(path) if file.endswith('xlsx')]". I use list to avoid other xlsx file under same folder.
file_lists = ['050e_228lib_S1_L001.xlsx','050h_228lib_S2_L001.xlsx','156e_228lib_S3_L001.xlsx',
                 'BE4_228lib_S4_L001.xlsx',  'pYW093uc_228lib_S24_L001.xlsx', 'pYW094ub_228lib_S23_L001.xlsx',
                 'pYW095ua_228lib_S16_L001.xlsx', 'pYW095ub_228lib_S17_L001.xlsx', 'pYW096uc_228lib_S20_L001.xlsx',
                 'pYW097ua_228lib_S13_L001.xlsx', 'pYW097ub_228lib_S14_L001.xlsx',  'pYW105b_228lib_S51_L001.xlsx',
                 'pYW105ue_228lib_S26_L001.xlsx', 'pYW106ub_228lib_S21_L001.xlsx',
                 'pYW106uc_228lib_S22_L001.xlsx', 'pyx165f_228lib_S7_L001.xlsx', 'pYX203ua_228lib_S5_L001.xlsx',
                 'pYX204ua_228lib_S11_L001.xlsx', 'pyx205ub_228lib_S12_L001.xlsx',  'pYX206ua_228lib_S15_L001.xlsx',
                 'pYX207ua_228lib_S18_L001.xlsx', 'pYX207ud_228lib_S19_L001.xlsx', 'pyx208uc_228lib_S10_L001.xlsx',
                 'pyx220b_228lib_S8_L001.xlsx', 'pyx220c_228lib_S9_L001.xlsx',  'pYX217d_228lib_S6_L001.xlsx']

# Convert mapping output from PROTOSPACER per line to TARGET_C per line.
def format_convert(wtbase,edbase):
    for num,file in enumerate(file_lists):
        excelwriter = pd.ExcelWriter(
            f'/Users/yuanwu/BaseSpace/010824-nextseq/041424-aio_table/Adata/{final_rename_dict[file.split("_")[0]]}_{file.split("_")[0]}.xlsx')
        rf = pd.read_excel(path + file, index_col=False)
        column = ['ID', 'gRNA', 'PAM', 'target_seq', 'edit_site_in_proto', 'context', 'clinvar_id', 'disease_name',
                  'gene_info', 'chr', 'coordinate', 'position', 'target_context']
        aio_potency_ot_df_list, aio_purity_ot_df_list = [], []
        for _ in range(len(rf['gRNA'])):
            total_read = rf['depth'][_]
            gRNA = rf['gRNA'][_]
            sequence = rf['target_seq'][_]
            info = rf.iloc[_].to_list()[:11]
            total_As = 0
            if total_read <= depth or total_read == '':
                continue
            for pos in range(20):
                if gRNA[pos] == wtbase:
                    total_As += rf[f'pos_{str(pos + 1)}_{wtbase}_to_{edbase}'][_] / total_read
            for pos in range(20):
                if gRNA[pos] == wtbase:
                    context = sequence[9 + pos:12 + pos].upper()
                    if total_As == 0:
                        potency = 0
                        purity = np.nan
                    else:
                        potency = rf[f'pos_{str(pos + 1)}_{wtbase}_to_{edbase}'][_]/total_read
                        purity = potency/total_As
                    aio_potency_ot_df_list.append(info + [pos + 1, context, potency])
                    aio_purity_ot_df_list.append(info + [pos + 1, context, purity])
        pd.DataFrame(aio_potency_ot_df_list, columns=column + [final_rename_dict[file.split("_")[0]]]).to_excel(excelwriter, 'potency')
        pd.DataFrame(aio_purity_ot_df_list, columns=column + [final_rename_dict[file.split("_")[0]]]).to_excel(excelwriter, 'purity')
        excelwriter.close()
        print(f'{file} is finished.')
    return

# Combined converted table to one excel sharing same index for convenience of removal NAN.
def table_combine(wtbase):
    excelwriter = pd.ExcelWriter('/Users/yuanwu/BaseSpace/010824-nextseq/aio-table/010824_aio_table_041424_A.xlsx')
    rf = []
    list1 = []
    index = pd.read_excel(
        '/Users/yuanwu/YuanWu_Tang_lab/publication/Topipotent adenine base editor/C_disease_lib/09082023_disease_C_lib.xlsx',
        sheet_name='summary_final')
    for num,gRNA in enumerate(index['gRNA']):
        for pos in range(20):
            if gRNA[pos] == wtbase:
                list1.append(['_'.join([gRNA, str(pos + 1)]),index['target_seq'][num][9 + pos:12 + pos].upper()]+index.iloc[num].to_list()[1:])
    list1 = pd.DataFrame(list1, columns=['gRNA','target_context','proto','PAM','target_seq','edit_site_in_proto','context','clinvar_id','disease_name','gene_info','chr','coordinate']).drop_duplicates(subset=['gRNA'],keep='first').set_index('gRNA')
    concat = list1.copy(deep=True)
    for option in ['purity', 'potency']:
        for num, file in enumerate(file_lists):
            rf.append(pd.read_excel(f'/Users/yuanwu/BaseSpace/010824-nextseq/041424-aio_table/Adata/{final_rename_dict[file.split("_")[0]]}_{file.split("_")[0]}.xlsx',
                                    sheet_name=option, index_col=0))
            rf[-1]['index'] = rf[-1]['gRNA'] + '_' + rf[-1]['position'].map(str)
            rf[-1] = rf[-1].set_index('index')
            concat = pd.concat([concat, rf[-1].iloc[:, -1]], axis=1)
            print(len(concat))
        concat.to_excel(excelwriter,option)
    excelwriter.close()
    return

if __name__ == '__main__':
    format_convert('A','G')
    table_combine('A')



