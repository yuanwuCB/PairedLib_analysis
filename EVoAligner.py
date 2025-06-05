"""
<EVA> -- EVolution library NGS Alignment
Tang Lab Evolution Library Mapping Pacakge.
Updated on 05/06/2025 by Yuan Wu
"""

import re
import gzip
import time
import pandas as pd
import warnings
import logomaker
import numpy as np
from Bio.Seq import Seq
from Bio import SeqIO
import matplotlib.cm as cm
import matplotlib.pyplot as plt

# Configuration (Please label your randomized region with NNK).
#TadA REFERENCE_DNA = 'GgcaaaAcgTgcatgggaCgaaNNKNNKNNKcccgtGggCgcAgtTctggtgcacaacaaCCgGgtTatTggTgaAggCtggAATNNKNNKNNKNNKNNKcaTgaTccGacCgcAcaTgcAgaAatCatggcCctgagAcagggCggActggtTatgcagaattaTcgTctgatcgaCgcCacActgtaCgtgacCctggaAccTtgTgtgatgtgcgcCggCgcCatgatccacTCTCggatCggTCgTgtTgtTttcggCgcaNNKNNKNNKNNKNNKggTgcAgcAggTAGcctgatg'
#AaTadA REFERENCE_DNA = 'cgaaacgcgcgtttgaaaaaNNKNNKNNKccggtgggcgcgattattgtgaaagaaggcgaaattattagcaaagcgcatNNKNNKNNKNNKNNKNNKaaagatccgaccgcgcatgcggaaatgctggcgattaaagaagcgtgccgccgcctgaacaccaaatatctggaaggctgcgaactgtatgtgaccctggaaccgtgcattatgtgcagctatgcgctggtgctgagccgcattgaaaaagtgatttttagcgcgNNKNNKNNKNNKNNKggcggcgtggtgagcgtg'
#AaTadA_shuffle_N REFERENCE_DNA = 'gcaaagaatattttctgaaagtggcgctgcgcgaagcgaaacgcgcgNNKgaaaaaggcgaagtgccggtgggcgcgattattgtgNNKgaaggcgaaattattagcaaagcgcataacNNKNNKgaagaactgaaagatccgaccgcgcatgcggaaatgctggcgattaaagaagcgtgccgccgcctgaacaccaaatatctgNNKggctgcgaactgtatNNKaccNNKgaaccgtgcattatgtgcagctatgcgctggtgctgagccgcattgaaaaagtgatttttagc'
REFERENCE_DNA = 'aagcgtgccgccgcctgaacaccaaatatctgNNKggctgcgaactgtatNNKaccNNKgaaccgtgcattatgtgcagctatgcgctggtgctgagccgcattgaaaaagtgatttttagcNNKctgNNKNNKNNKNNKggcNNKgtggtgagcgtgtttNNKattctgNNKNNKccgaccNNKNNKcatcgcgtgaaatgggaatattatccgctggaagaagcgagcgaactgctgNNKNNKtttNNKaaaaaaNNKcgcNNKNNKNNKNNKtccggtagcgaaacaccggggact'
# for some circumstances (Cannot find one by the code), input a preset hook to the code
Hook = ['gaaccgtgcattatg'.upper(),59]
# directory of raw reads
input_path = '/Users/yuanwu/Downloads/20250601_NextSeq/'
# directory of where you want the output file to
output_path = '/Users/yuanwu/Downloads/20250601_NextSeq/'
# All the barcodes of FASTQ files that you want to analyze
Barcode_list = [278,279,280,281,282,283]
# If you want to have correct amino acid label in your LogoPlot.
Random_Pos = [71,77,79,101,103,104,105,106,108,114,117,118,121,122,140,141,143,146,148,149,150,151]
# Nominate a list that you want to exclude from your reads when specific amino acids show up. (e.g. Stop Codon/wild type)
Exclude_List = ['X','*']

# https://matplotlib.org/stable/gallery/text_labels_and_annotations/font_family_rc.html
plt.rcParams['font.family'] = 'Arial'
# https://matplotlib.org/stable/users/explain/text/fonts.html
plt.rcParams["svg.fonttype"] = 'none'
COLOR_SCHEME = {'A':'g','T':'r','C':'b','G':'y','N':'k'}
warnings.filterwarnings("ignore")

def ATCG_percentage_plot(SeqDF,QualityDF,REF,file):
    def myplot(x, y, bins):
        heatmap, xedges, yedges = np.histogram2d(x[np.isfinite(x)], y[np.isfinite(y)], bins=bins)
        extent = [xedges[0], xedges[-1], yedges[0], yedges[-1]]
        return heatmap.T, extent

    def dataframe_to_frequency(df):
        frequency_df = pd.DataFrame()
        for column in df.columns:
            frequency_df[column] = df[column].value_counts(normalize=True).reindex(['A','T','G','C','N'],fill_value=0)
        return frequency_df

    def normalize_plot(BASE,ax):
        perc = [base_percentages_by_pos[pos].get(BASE, 0) for pos in positions]
        ax.plot(positions, perc, label=BASE, marker='o', color=COLOR_SCHEME[BASE], markersize=4)
        return

    base_percentages_by_pos = dataframe_to_frequency(SeqDF)
    positions = sorted(base_percentages_by_pos.keys())
    fig,axs = plt.subplots(nrows=3, ncols=1, figsize = (20,12), constrained_layout=True)
    for ax in (axs[0],axs[1]):
        for BASE in ['A','T','C','G','N']:
            normalize_plot(BASE,ax)
    axs[1].set_ylim(0,0.2)
    TOPY = axs[0].get_ylim()[1]
    for num, _ in enumerate(REF):
        try:
            ax.plot(num, TOPY, COLOR_SCHEME[_] + 'o')
        except KeyError:
            ax.plot(num, TOPY, 'ko')
    QualityDF = QualityDF.sample(frac=0.05,random_state=1).melt().dropna()
    img, extent = myplot(QualityDF.variable.to_numpy(), QualityDF.value.to_numpy(), [len(REF),10])
    axs[2].imshow(img, extent=extent, origin='lower', cmap=cm.jet)
    for Title,ax in zip(['ATCG Base Percentage at Each Genomic Position','ATCG Base Percentage at Each Genomic Position (Zoom In to 20%)','PhreadQuality_Score_Heatmap']
            ,axs.flatten()):
        ax.set_xlabel('Genomic Position')
        ax.set_ylabel('Count')
        ax.grid(True)
        ax.set_title(Title)
        ax.set_xlim(0, len(REF))
    for ax in (axs[0],axs[1]):
        ax.legend()
    plt.gca()
    plt.savefig(f'{output_path}TCDI{file}_ATCG_percentage_plot.pdf')
    plt.show()
    return

def Logoplot(df,title):
    fig,ax = plt.subplots(figsize=(0.2*len(Random_Pos),1.9),constrained_layout=True)
    df = df.div(df.sum(axis=0), axis=1).astype(float).T
    df = df.reset_index(drop=True)
    logomaker.Logo(df,
                   vpad=.05,
                   width=.8,
                   stack_order='fixed',
                   color_scheme='skylign_protein',
                   ax=ax)
    # style using Logo methods
    try:
        if Random_Pos:
            ax.set_xticks(range(len(Random_Pos)))
            ax.set_xticklabels(Random_Pos, fontsize=7)
    except:
        ax.tick_params(axis='x', labelsize=7)
    ax.set_xlabel(title, size=7)
    ax.set_yticks([0, .5, 1])
    ax.set_yticklabels([0, .5, 1], fontsize=7)
    ax.set_ylabel('Probability', size=7)
    plt.savefig(f'{output_path}{title}_Logoplot.svg', dpi=300)
    plt.savefig(f'{output_path}{title}_Logoplot.pdf', dpi=300)
    plt.show()
    return

def RANDOMIZE_REGION_FINDER(REFERENCE):
    RANDOMIZE_REGIONS = list(re.finditer('NNK', REFERENCE))
    tag = RANDOMIZE_REGIONS[0].end()
    RANDOMIZE_REGION = [RANDOMIZE_REGIONS[0].start(), tag]
    for region in RANDOMIZE_REGIONS[1:]:
        if region.start() != tag:
            RANDOMIZE_REGION += [region.start(), region.end()]
            tag = region.end()
        else:
            tag = region.end()
            RANDOMIZE_REGION[-1] = tag
            continue
    return RANDOMIZE_REGION

def ALIGNMENT(R1_FQ,R2_FQ,REFERENCE,file):
    def Hook_finder(REF,List):
        if List[1] <= 146:
            Hook_st_index = 1
            Hook_st = List[Hook_st_index]
        elif List[1] > 146:
            Hook_st = 70
        if len(List) > 2:
            if Hook_st + 20 > 145 or Hook_st + 20 > List[Hook_st_index + 2]:
                raise ValueError("Cannot find a suitable hook.")
        elif Hook_st + 20 > 145:
            raise ValueError("Cannot find a suitable hook.")
        Hook_Len = 12
        while REFERENCE.count(REF[Hook_st:Hook_st + Hook_Len]) > 1:
            Hook_Len += 1
            if Hook_Len >= 20:
                raise ValueError("Cannot find a suitable hook.")
        return [REF[Hook_st:Hook_st + Hook_Len],Hook_st]

    def PE_Merger(Read1,Read2,Read1_Q,Read2_Q,Ref,Hook):
        if Read1.find(Hook[0])-Hook[1] < 0:
            return None
        if len(Ref) >= 300 - 4:
            return (Read1[Read1.find(Hook[0]) - Hook[1]:] +
                   Ref [len(Read1[Read1.find(Hook[0])-Hook[1]:]) : len(Ref) - 150] +
                   str(Seq(Read2).reverse_complement()),
                    Read1_Q[Read1.find(Hook[0]) - Hook[1]:] +
                    [0 for _ in range(len(Ref) - 150 - (len(Read1[Read1.find(Hook[0]) - Hook[1]:])))] +
                    Read2_Q[::-1])
        else:
            return (Read1[Read1.find(Hook[0]) - Hook[1]:] +
                    str(Seq(Read2[:len(Ref)-len(Read1[Read1.find(Hook[0])-Hook[1]:])]).reverse_complement()),
                    Read1_Q[Read1.find(Hook[0]) - Hook[1]:] +
                    Read2_Q[len(Ref) - len(Read1[Read1.find(Hook[0]) - Hook[1]:]) - 1::-1]
                    )

    def AA_convert(Sequence,Ref):
        ST = Ref - Ref // 3 * 3
        return Seq(Sequence[ST:]).translate()

    def Randomize_AA_Finder(Read,AA_Seq,List):
        return ('_'.join(Read[List[i]:List[i+1]] for i in range(0,len(List),2)),
                ''.join(map(str,[AA_Seq[List[i]//3:List[i+1]//3] for i in range(0,len(List),2)])))

    # Find the randomize region in amplicon and hook for mapping.
    Reads_Counter = [0,0,0]
    Randomize_Region = RANDOMIZE_REGION_FINDER(REFERENCE)
    print(f'The randomized region in your sequence are:{Randomize_Region}')
    try:
        global Hook
        print(f'Using input Hook: {Hook}')
    except:
        Hook = Hook_finder(REFERENCE,Randomize_Region)

    # Parse reads from fastq file and match randomized region one by one.
    Len_REF_AA = len(AA_convert(REFERENCE,Randomize_Region[0]))
    AA_data,Randomize_dna_data,Reads_data,Quality_data,Unmapped_Reads = [],dict(),[],[],[]
    for (READ1,READ2) in zip(SeqIO.parse(gzip.open(R1_FQ, 'rt'), "fastq"),
                             SeqIO.parse(gzip.open(R2_FQ, 'rt'), "fastq")):
        Reads_Counter[0] += 1
        READ1_QS,READ2_QS = READ1.letter_annotations["phred_quality"],READ2.letter_annotations["phred_quality"]
        READ1_SEQ,READ2_SEQ = str(READ1.seq),str(READ2.seq)
        try:
            Merge_Read,Merge_Quality = PE_Merger(READ1_SEQ,READ2_SEQ,READ1_QS,READ2_QS,REFERENCE,Hook)
        except:
            for READ in [READ1,READ2]:
                Unmapped_Reads.append(READ)
            continue
        Reads_Counter[1] += 1
        Amino_Acid_Sequence = AA_convert(Merge_Read,Randomize_Region[0])
        if len(Amino_Acid_Sequence) != Len_REF_AA:
            continue
        Reads_Counter[2] += 1
        Randomize_Region_DNA,Randomize_Region_AA = Randomize_AA_Finder(Merge_Read,Amino_Acid_Sequence,Randomize_Region)
        if any(Exclude_AA in Randomize_Region_AA for Exclude_AA in Exclude_List):
            Including_Tag = 0
        else:
            Including_Tag = 1
        Quality_data.append(Merge_Quality)
        Reads_data.append(list(Merge_Read) + [Including_Tag])
        AA_data.append(list(str(Amino_Acid_Sequence)) + [Including_Tag])
        if Randomize_Region_AA in Randomize_dna_data.keys():
            if Randomize_Region_DNA in Randomize_dna_data[Randomize_Region_AA][0]:
                Randomize_dna_data[Randomize_Region_AA][0][Randomize_Region_DNA] += 1
            else:
                Randomize_dna_data[Randomize_Region_AA][0][Randomize_Region_DNA] = 1
        else:
            Randomize_dna_data[Randomize_Region_AA] = [{Randomize_Region_DNA: 1},Including_Tag]

    # Output data to excel..
    List_df = pd.DataFrame(AA_data)
    Randomize_aa_data_df = pd.DataFrame.from_dict(Randomize_dna_data).T.reset_index(drop=False).rename(
        columns={'index': 'AA', 0: 'Codon_Dict', 1: 'Including_Tag'})
    Randomize_aa_data_df['Frequency'] = Randomize_aa_data_df['Codon_Dict'].apply(lambda x: sum(x.values()))
    for Including_Statu,Tag in zip(['no_filter','with_filter'],[0,1]): # 0 meangs include all reads. 1 means include only reads without AA in Exclude_List.
        AA_Output_pd = pd.DataFrame(0,
                                    columns=list(AA_convert(REFERENCE, Randomize_Region[0])),
                                    index=['A', 'C', 'D', 'E', 'F', 'G', 'H', 'I', 'K', 'L', 'M', 'N', 'P', 'Q', 'R',
                                           'S', 'T', 'V', 'W', 'Y', '*', 'X'])
        List_sub_df = List_df[List_df.iloc[:,-1] >= Tag].iloc[:,:-1]
        for num, Ref_aa in enumerate(AA_Output_pd.columns):
            AA_Output_pd.iloc[:, num] = List_sub_df.groupby(num).size().reindex(AA_Output_pd.index, fill_value=0)
        AA_Output_pd.to_excel(Excelwriter, sheet_name=f"Full_Length_Region_{Including_Statu}")
        Randomize_aa_data_sub_df = Randomize_aa_data_df[Randomize_aa_data_df['Including_Tag'] >= Tag].sort_values(by=['Frequency'], ascending=False)
        Randomize_aa_data_sub_df = Randomize_aa_data_sub_df.reset_index(drop=True)
        Randomize_aa_data_sub_df.index = Randomize_aa_data_sub_df.index + 1
        Randomize_aa_data_sub_df[['AA','Codon_Dict','Frequency']].to_excel(Excelwriter,sheet_name=f'Randomize_Region_{Including_Statu}')
    print(f'Total reads number of fastq is {Reads_Counter[0]};\n'
          f'The number of Successfully aligned reads is {Reads_Counter[1]};\n'
          f'The number of reads with correct length is {Reads_Counter[2]}.')
    # Plot the whole mapping result as curve plot.
    ATCG_percentage_plot(pd.DataFrame(Reads_data),pd.DataFrame(Quality_data),REFERENCE_DNA,file)
    # Plot the logo plot.
    Logoplot(AA_Output_pd.iloc[:,[pos for _ in [range(Randomize_Region[i]//3,Randomize_Region[i+1]//3) for i in range(0,len(Randomize_Region),2)] for pos in _]],f'TCDI{file}')
    # Output unmapped reads.
    with open(f'{output_path}TCDI{file}.unmapped.fastq','w') as file:
        SeqIO.write(Unmapped_Reads,file,'fastq')
    return

def main():
    global REFERENCE_DNA, input_path, output_path
    # Reformat input and find the hook position.
    REFERENCE_DNA = REFERENCE_DNA.upper()
    # Start alignment.
    for file in Barcode_list:
        global Excelwriter
        T1 = time.time()
        print(f'Starting mapping of {file}...')
        Excelwriter = pd.ExcelWriter(f'{output_path}TCDI{file}.xlsx')
        ALIGNMENT(f'{input_path}TCDI{file}_S{file}_L001_R1_001.fastq.gz',
                  f'{input_path}TCDI{file}_S{file}_L001_R2_001.fastq.gz',
                  REFERENCE_DNA,
                  file)
        Time = time.time() - T1
        print(f'Finished mapping of {file}...Using time: {Time:.2f} s')
        Excelwriter.close()
    return

if __name__ == '__main__':
    main()