import pandas as pd
import os
from _Final_Name_List.py import *
import seaborn as sns
import matplotlib.pyplot as plt
import numpy as np
from matplotlib import colors
from matplotlib.colors import LinearSegmentedColormap
import warnings

warnings.filterwarnings(action='ignore')

# color bar for protospacer plot
cmap_c = LinearSegmentedColormap.from_list('custom_c',['white',np.array((255,201,212,256))/256.0,np.array((256,45,90,256))/256.0],N=256)
# aio-table color bar for main figure
cmap_c_cbar = LinearSegmentedColormap.from_list('custom_a',['white',np.array((256,45,90,256))/256.0],N=256)
c_sm = plt.cm.ScalarMappable(cmap=cmap_c_cbar,norm=colors.Normalize(vmin=0,vmax=100))

# CONSTANT PARAMETER
path = '/Users/yuanwu/BaseSpace/010824-nextseq/aio-table/' #
aio_table = '010824_aio_table_041424_C_unique.xlsx'

# PLOTING FORMAT
plt.rcParams['font.family'] = ['Arial']
plt.rcParams.update({'font.size': 7})

def NCN_heatmap(): # For Fig. S16
    guide = pd.read_excel(path + aio_table, sheet_name='purity', index_col=0).dropna()
    index = guide.index
    pos = 4
    proto = guide['proto'].to_list()
    for file in [_ for _ in os.listdir(path) if _.endswith('C_unique.xlsx')]:
        potency_rf = pd.read_excel(path + aio_table, sheet_name='potency', index_col=0).dropna().loc[index]
        for editor in potency_rf.columns[12:]:
            potency_ot_df = pd.DataFrame(0,index=['A','T','C','G'],columns=['A','T','C','G']).fillna(0)
            site_ot_df = pd.DataFrame(0,index=['A','T','C','G'],columns=['A','T','C','G']).fillna(0)
            for context in contexts:
                potency_ot_df.loc[context[0],context[2]] = potency_rf[(potency_rf['position'] == pos) &
                                                                      (potency_rf['target_context'] == context)].loc[:,editor].mean()*100
                site_ot_df.loc[context[0],context[2]] = potency_rf[(potency_rf['position'] == pos) &
                                                                   (potency_rf['target_context'] == context)].loc[:,editor].size
            potency_ot_df.to_excel(excelwriter,f'{editor}_potency')
            site_ot_df.to_excel(excelwriter, f'{editor}_site')
    return


def C_Disease_Lib():
    aio_table = '010824_aio_table_041424_C_unique.xlsx'
    guide = pd.read_excel(path + aio_table, sheet_name='purity', index_col=0).dropna()
    index = guide.index
    proto = guide['gRNA'].to_list()
    pos = 6
    for file in [_ for _ in os.listdir(path) if _.endswith('C_unique.xlsx')]:
        purity_rf = pd.read_excel(path+file, sheet_name='purity',index_col=0).dropna().loc[index]
        potency_rf = pd.read_excel(path+file,sheet_name='potency',index_col=0).dropna().loc[index]

        # extract target C data
        # target_C_potency = potency_rf.copy(deep=True).iloc[:,:16]
        # target_C_potency['Tad-CBE(NCN)'] = 0
        # target_C_purity = purity_rf.copy(deep=True).iloc[:,:16]
        # target_C_purity['Tad-CBE(NCN)'] = 0
        # for _ in range(len(target_C_potency['proto'])):
        #     target_C_potency.loc[target_C_potency.index[_],'Tad-CBE(NCN)'] = potency_rf.loc[target_C_potency.index[_],context2editor[target_C_potency['target_context'][_]]]
        # for _ in range(len(target_C_purity['proto'])):
        #     target_C_purity.loc[target_C_purity.index[_],'Tad-CBE(NCN)'] = purity_rf.loc[target_C_purity.index[_],context2editor[target_C_purity['target_context'][_]]]
        # target_C_potency.to_excel(excelwriter,'potency')
        # target_C_purity.to_excel(excelwriter,'purity')

        print(len(purity_rf.drop_duplicates(subset=['clinvar_id'])))
        filename = '_'.join(file.rstrip('.xlsx').split('_')[-2:])
        # calculate mean potency or purity
        pd.concat([purity_rf[purity_rf['position']==pos].iloc[:,12:].mean().T for pos in range(1,21,1)],axis=1)\
            .set_axis(range(1,21,1),axis=1)\
            .rename({'Tad-CBE(ACG)':'Tad-CBE(ACG1)'})\
            .reindex(editor_order)\
            .to_excel(excelwriter,f'{filename}_purity')
        pd.concat([potency_rf[potency_rf['position']==pos].iloc[:,12:].mean().T for pos in range(1,21,1)],axis=1)\
            .set_axis(range(1,21,1),axis=1)\
            .rename({'Tad-CBE(ACG)':'Tad-CBE(ACG1)'})\
            .reindex(editor_order)\
            .to_excel(excelwriter,f'{filename}_potency')
        pd.DataFrame([len(potency_rf[potency_rf['position']==pos]) for pos in range(1,21,1)],index=range(1,21,1),columns=['n'])\
            .T.reindex(editor_order) \
            .rename({'Tad-CBE(ACG)': 'Tad-CBE(ACG1)'}) \
            .reindex(editor_order)\
            .to_excel(excelwriter,f'{filename}_site_number')

    for file in ['010824_aio_table_041424_A_unique.xlsx']:
        potency_rf = pd.read_excel(path + file, sheet_name='potency', index_col=0).dropna()
        filename = '_'.join(file.rstrip('.xlsx').split('_')[-2:])
        potency_rf = pd.concat([potency_rf.iloc[num,:] for num,pro in enumerate(potency_rf['proto']) if pro in proto],axis=1).T
        print(len(potency_rf.drop_duplicates(subset=['clinvar_id'])))
        pd.concat([potency_rf[potency_rf['position'] == pos].iloc[:, 12:].mean().T for pos in range(1, 21, 1)], axis=1) \
            .set_axis(range(1, 21, 1), axis=1)\
            .rename({'Tad-CBE(ACG)': 'Tad-CBE(ACG1)'}) \
            .reindex(editor_order)\
            .to_excel(excelwriter, f'{filename}_potency')
        pd.DataFrame([len(potency_rf[potency_rf['position'] == pos]) for pos in range(1, 21, 1)], index=range(1, 21, 1),
                     columns=['n']) \
            .rename({'Tad-CBE(ACG)': 'Tad-CBE(ACG1)'}) \
            .reindex(editor_order)\
            .T.to_excel(excelwriter, f'{filename}_site_number')

# For data in S15B
def A_Disease_Lib():
    for file in ['Alib-010824_aio_table_041424_C.xlsx']:
        potency_rf = pd.read_excel(path+file,sheet_name='potency',index_col=0).dropna()
        filename = '_'.join(file.rstrip('.xlsx').split('_')[-2:])
        print(len(potency_rf.drop_duplicates(subset=['clinvar_id'])))
        pd.concat([potency_rf[potency_rf['position'] == pos].iloc[:, 12:].mean().T for pos in range(1, 21, 1)], axis=1) \
            .set_axis(range(1, 21, 1), axis=1) \
            .to_excel(excelwriter, f'{filename}_potency')

# For Fig. S17
def purity_vs_potency():
    def label_point(x, y, val, ax, con):
        a = pd.concat({'x': x, 'y': y, 'val': val}, axis=1)
        for i, point in a.iterrows():
            # if str(point['val']).split('(')[-1][:3] == con:
            #     ax.text(point['x'] + .005, point['y'], str(point['val']),fontsize=10,weight='bold')
            # else:
            ax.text(point['x'] - .02, point['y'], str(point['val']),fontsize=6)

    purity_rf = pd.read_excel(f'{path}{aio_table}',
                              sheet_name='purity',
                              index_col=0).dropna()
    index = purity_rf.index
    potency_rf = pd.read_excel(f'{path}{aio_table}',
                               sheet_name='potency',
                               index_col=0).loc[index]
    excelwriter = pd.ExcelWriter(f'{path}potency_vs_purity_designate_6-042924.xlsx',engine="openpyxl")

    fig, axs = plt.subplots(4, 4, figsize=(8,6), constrained_layout=True,tight_layout=False)
    for pos in range(6,7,1):
        for num,context in enumerate(contexts):
            ax=axs[num//4,num%4]
            potency_ot, purity_ot, site_ot = pd.DataFrame(index=purity_rf.columns.to_list()[12:],columns=['potency']), \
                                             pd.DataFrame(index=purity_rf.columns.to_list()[12:],columns=['purity']), \
                                             pd.DataFrame(index=purity_rf.columns.to_list()[12:],columns=['site_num'])
            editor_list = control+[context2editor[context]]
            purity_series = purity_rf[(purity_rf['position'] == pos) & (purity_rf['target_context'] == context)]
            potency_series = potency_rf[(potency_rf['position'] == pos) & (potency_rf['target_context'] == context)]

            for editor in editor_list:
                if editor in control:
                    purity, potency, site_num = purity_series.loc[:,editor].mean(), potency_series.loc[:,editor].mean(), len(purity_series.loc[:,editor])
                else:
                    purity,potency,site_num = purity_series.loc[:,editor].mean(), potency_series.loc[:,editor].mean(), len(purity_series.loc[:,editor])
                potency_ot.loc[editor,'potency'] = potency
                purity_ot.loc[editor,'purity'] = purity
                site_ot.loc[editor,'site_num'] = site_num
            data = pd.concat([potency_ot,purity_ot,site_ot],axis=1).dropna()
            data.to_excel(excelwriter,f'{pos}_{context}')
            sns.scatterplot(x='potency',y='purity',data=data,s=6,ax=ax)
            ax.set_xlabel('C:G-to-T:A conversion (%)',fontsize=7)
            # ax.set_xlim(0,0.8)
            # ax.set_ylim(0.2,0.8)
            ax.set_ylabel('Target C editing/\nTotal C editing',fontsize=7)
            ax.set_title(context,fontsize=8)
            data['editor']=data.index.to_list()
            label_point(data['potency'],data['purity'],data['editor'],ax,context)
        plt.savefig(f'{path}purity_vs_potency_{pos}_designate.svg',dpi=600)
        plt.show()
    excelwriter.close()
    return

# For Fig. S18
def aio_heatmap(position,files,context,bool):
    def norm(raw):
        # return [float(i) / normalization * 100 for i in raw]
        return [float(i) / max(raw) * 100 for i in raw]

    share_column = []
    counter = [0,0]
    for file in files:
        share_column.append(rf.columns.to_list().index(file))
    fig, axs = plt.subplots(1,len(files),figsize=(8, 1.4), constrained_layout=True)
    for num,file in enumerate(files):
        column = rf.columns.to_list().index(file)
        ax = axs[num]
        if (num == 0 or bool == True) and context != '':
            gRNA_in_order = rf[(rf['position'] == position) & (rf['target_context'] == context)].iloc[:,
                            list(range(12)) + share_column].sort_values(by=[rf.columns.to_list()[column]],
                                                                  ascending=False)['gRNA']
        elif (num == 0 or bool == True) and context == '':
            gRNA_in_order = rf[rf['position'] == position].iloc[:,
                            list(range(12)) + share_column].sort_values(by=[rf.columns.to_list()[column]],
                                                                        ascending=False)['gRNA']
        df = pd.read_excel(f'/Users/yuanwu/BaseSpace/010824-nextseq/excel_output/{editor2file[file]}',index_col=False)
        datas,ytick = [],[]
        for i,gRNA in enumerate(gRNA_in_order):
            if num == 0: counter[0] += 1
            data = []
            tag = df[df['gRNA'] == gRNA].index.to_list()[0]
            total_read = df['depth'][tag]
            for pos in range(20):
                data.append(df['pos_' + str(pos + 1) + '_C_to_T'][tag] / total_read * 100)
            datas.append(norm(data))
            # if (max(datas[-1]) > 200 or df['pos_' + str(position) + '_C_to_T'][tag] / total_read * 100 == 0):
            #     if num==0:
            #         counter[1] += 1
            #     ytick.append(i)
        df = pd.DataFrame(datas, columns=list(range(1,21,1)), index=list(range(len(datas))))
        #df.to_excel(f'{ot_path}{final_rename_dict[file.split("_")[0]]}.xlsx')
        N = df.index.size
        value = np.array(df)
        c_norm = colors.Normalize(vmin=0, vmax=100)
        heatmap = []
        for m in range(20):
            heatmap.append(ax)
        # first heatmap
        x_label = []  # put value definition here if you want switch between black and white
        for m in range(len(df.columns)):
            column_name = df.columns[m]
            x_label.append(gRNA[m])
            heatmap[m] = ax.imshow(np.vstack([df[column_name], df[column_name]]).T, aspect='auto',
                                   extent=[0.5+m, m+1.5, -0.5, N-0.5], origin='lower', cmap=cmap_c,
                                   norm=c_norm)



        # for i in range(len(df.index)):
        #     for j in range(len(df.columns)):
        #         if value[i, j] != 'nan' and gRNA[j][0] == 'C' and value[i, j] > 10:
        #             text = ax.text(j - 9.5, i, format(value[i, j], '.1f'),
        #                            ha="center", va="center",
        #                            color='black', fontsize=2.3,
        #                            rotation='vertical')  # color=textcolors[heatmap[j].norm(value[i, j]) > threshold]
        ax.set_xlim(0.5, 20.5)
        ax.set_xticks(list(range(1,21,1)))
        ax.xaxis.set_ticks_position('top')
        ax.set_xticklabels(list(range(1,21,1)),fontsize=3)
        ax.set_yticks(ytick,labels=None)
        ax.set_yticklabels(['*']*len(ytick))
        ax.yaxis.set_ticks_position('none')
        ax.tick_params(axis='y', which='major', labelsize=5, pad=-2)
        ax.tick_params(axis='x', which='major', labelsize=5, pad=1, length=2)
        ax.invert_yaxis()
        ax.set_title(file, fontstyle='italic', fontsize=8)
    # fig.subplots_adjust(right=0.8)
    # cbar_ax = fig.add_axes([0.85,0.15,0.05,0.7])
    cbar1 = fig.colorbar(c_sm, ax=axs.ravel(),fraction=0.04/len(files), pad=0.01) # ticks=range(0,101,20), # (int(df.stack().max())//10+1)*10, int(((int(df.stack().max())//10+1)*10)/5)), fraction=0.04/len(files), pad=0.02)
    cbar1.set_label('Normalized C:G-to-T:A\nconversion (%)', size=6, labelpad=0.8, linespacing=1)
    cbar1.ax.tick_params(labelsize=6, length=2, pad=-0.01)
    # cbar2 = fig.colorbar(a_sm, ax=axs.ravel(), fraction=0.04 / len(files),
    #                      pad=0.02)  # ticks=range(0,101,20), # (int(df.stack().max())//10+1)*10, int(((int(df.stack().max())//10+1)*10)/5)), fraction=0.04/len(files), pad=0.02)
    # cbar2.set_label('Normalized A:T-to-C:G efficiency (%)', size=8, labelpad=0.8, linespacing=1)
    # cbar2.ax.tick_params(labelsize=8, length=2, pad=-0.01)
    plt.savefig(
        f'{path}heatmap_5D_shared/{files[0]}_{position}_{context}_shared.svg',
        dpi=300)
    plt.show()
    print([files[0]]+list(map(str,counter)))
    return

# For Fig S22 and Fig 5G
def aio_heatmap_Tad_CBEs(position,files,context,bool):
    def norm(raw):
        # return [float(i) / normalization * 100 for i in raw]
        return [float(i) / max(raw) * 100 for i in raw]
    rf = pd.read_excel(f'/Users/yuanwu/BaseSpace/010824-nextseq/aio-table/best-perform-Tad-CBEs(NCN)-winning_no_match_test.xlsx',index_col=False)
    counter = [0,0]
    fig, axs = plt.subplots(1,len(files),figsize=(8, 1.4), constrained_layout=True)
    for num,file in enumerate(files):
        ax = axs[num]
        if (num == 0 or bool == True) and context == '':
            gRNA_in_order = rf[rf['position'] == position].iloc[:,
                            list(range(12)) + [32,33,34,35,36]].sort_values(by=['max_purity'],ascending=False)['gRNA']
        df = pd.read_excel(f'/Users/yuanwu/BaseSpace/010824-nextseq/excel_output/{editor2file[file]}',index_col=False)
        datas,ytick = [],[]
        for i,gRNA in enumerate(gRNA_in_order):
            if num == 0: counter[0] += 1
            data = []
            tag = df[df['gRNA'] == gRNA].index.to_list()[0]
            total_read = df['depth'][tag]
            for pos in range(20):
                data.append(df['pos_' + str(pos + 1) + '_C_to_T'][tag] / total_read * 100)
            datas.append(norm(data))
            # if (max(datas[-1]) > 200 or df['pos_' + str(position) + '_C_to_T'][tag] / total_read * 100 == 0):
            #     if num==0:
            #         counter[1] += 1
            #     ytick.append(i)
        df = pd.DataFrame(datas, columns=list(range(1,21,1)), index=list(range(len(datas))))
        #df.to_excel(f'{ot_path}{final_rename_dict[file.split("_")[0]]}.xlsx')
        N = df.index.size
        value = np.array(df)
        c_norm = colors.Normalize(vmin=0, vmax=100)
        heatmap = []
        for m in range(20):
            heatmap.append(ax)
        # first heatmap
        x_label = []  # put value definition here if you want switch between black and white
        for m in range(len(df.columns)):
            column_name = df.columns[m]
            x_label.append(gRNA[m])
            heatmap[m] = ax.imshow(np.vstack([df[column_name], df[column_name]]).T, aspect='auto',
                                   extent=[0.5+m, m+1.5, -0.5, N-0.5], origin='lower', cmap=cmap_c,
                                   norm=c_norm)



        # for i in range(len(df.index)):
        #     for j in range(len(df.columns)):
        #         if value[i, j] != 'nan' and gRNA[j][0] == 'C' and value[i, j] > 10:
        #             text = ax.text(j - 9.5, i, format(value[i, j], '.1f'),
        #                            ha="center", va="center",
        #                            color='black', fontsize=2.3,
        #                            rotation='vertical')  # color=textcolors[heatmap[j].norm(value[i, j]) > threshold]
        ax.set_xlim(0.5, 20.5)
        ax.set_xticks(list(range(1,21,1)))
        ax.xaxis.set_ticks_position('top')
        ax.set_xticklabels(list(range(1,21,1)),fontsize=3)
        ax.set_yticks(ytick,labels=None)
        ax.set_yticklabels(['*']*len(ytick))
        ax.yaxis.set_ticks_position('none')
        ax.tick_params(axis='y', which='major', labelsize=5, pad=-2)
        ax.tick_params(axis='x', which='major', labelsize=5, pad=1, length=2)
        ax.invert_yaxis()
        ax.set_title(file if file != 'Tad-CBEs(NCN)' else 'Tad-CBE(NCN)', fontsize=8)
    # fig.subplots_adjust(right=0.8)
    # cbar_ax = fig.add_axes([0.85,0.15,0.05,0.7])
    cbar1 = fig.colorbar(c_sm, ax=axs.ravel(),fraction=0.04/len(files), pad=0.01) # ticks=range(0,101,20), # (int(df.stack().max())//10+1)*10, int(((int(df.stack().max())//10+1)*10)/5)), fraction=0.04/len(files), pad=0.02)
    cbar1.set_label('Normalized C:G-to-T:A\nconversion (%)', size=6, labelpad=0.8, linespacing=1)
    cbar1.ax.tick_params(labelsize=6, length=2, pad=-0.01)
    # cbar2 = fig.colorbar(a_sm, ax=axs.ravel(), fraction=0.04 / len(files),
    #                      pad=0.02)  # ticks=range(0,101,20), # (int(df.stack().max())//10+1)*10, int(((int(df.stack().max())//10+1)*10)/5)), fraction=0.04/len(files), pad=0.02)
    # cbar2.set_label('Normalized A:T-to-C:G efficiency (%)', size=8, labelpad=0.8, linespacing=1)
    # cbar2.ax.tick_params(labelsize=8, length=2, pad=-0.01)
    plt.savefig(
        f'{path}heatmap_5E_shared/{files[0]}_{position}_{context}_shared.svg',
        dpi=300)
    plt.show()
    print([files[0]]+list(map(str,counter)))
    return

# For Fig S14
def potency_heatmap():
    gRNA = 'C'*20
    df = pd.read_excel(f'{path}mean_by_position_order.xlsx', keep_default_na=False, sheet_name='C_unique_potency_percentage',index_col=0)
    fig = plt.figure(figsize=(4.2, 0.143 * len(df.index)))
    ax = fig.add_axes([0.2, 0.1, 0.8, 0.8])
    N = df.index.size
    value = np.array(df)
    c_norm = colors.Normalize(vmin=0, vmax=100)
    a_norm = colors.Normalize(vmin=0, vmax=100)
    heatmap = []
    threshold = [0 for _ in range(len(gRNA))]
    for m in range(len(gRNA)):
        heatmap.append(ax)
    # first heatmap
    x_label = []  # put value definition here if you want switch between black and white
    for m in range(len(df.columns)):
        base = 'C'
        x_label.append(base)
        column_name = df.columns[m]
        if base == 'A':
            heatmap[m] = ax.imshow(np.vstack([df[column_name], df[column_name]]).T, aspect='auto',
                                   extent=[-10 + m, -9 + m, -0.5, N - 0.5], origin='lower', cmap=cmap_a, norm=a_norm)
            threshold[m] = heatmap[m].norm(value.max()) / 2.
        elif base == 'C':
            heatmap[m] = ax.imshow(np.vstack([df[column_name], df[column_name]]).T, aspect='auto',
                                   extent=[-10 + m, -9 + m, -0.5, N - 0.5], origin='lower', cmap=cmap_c, norm=c_norm)
            threshold[m] = heatmap[m].norm(value.max()) / 2.
        else:
            heatmap[m] = ax.imshow(np.vstack([df[column_name], df[column_name]]).T, aspect='auto',
                                   extent=[-10 + m, -9 + m, -0.5, N - 0.5], origin='lower', cmap='Greys', vmin=0,
                                   vmax=100)
            threshold[m] = heatmap[m].norm(value.max()) / 2.

    print(x_label)
    try:
        A_posi, C_posi = x_label.index("A"), x_label.index('C')
        cbar1 = fig.colorbar(heatmap[C_posi], ax=ax, ticks=range(0, 101, 20),
                             pad=0.03, aspect=40)
        cbar2 = fig.colorbar(heatmap[A_posi], ax=ax, ticks=range(0, 101, 20),
                             pad=0.01, aspect=40)
        cbar1.set_label('C:G-to-T:A conversion (%)', size=7, labelpad=0)
        cbar1.ax.tick_params(labelsize=7, length=2)
        cbar2.set_label('A:T-to-G:C conversion (%)', size=7, labelpad=0)
        cbar2.ax.tick_params(labelsize=7, length=2)
    except:
        C_posi = x_label.index('C')
        cbar1 = fig.colorbar(heatmap[C_posi], ax=ax, ticks=range(0, 101, 20),
                             pad=0.03, aspect=40)
        cbar1.set_label('C:G-to-T:A conversion (%)', size=7, labelpad=0)
        cbar1.ax.tick_params(labelsize=7, length=2)
    for i in range(len(df.index)):
        for j in range(len(df.columns)):
            if value[i, j] != 'nan' and (gRNA[j][0] == 'C' or gRNA[j][0] == 'A') and value[i, j] > 5:
                text = ax.text(j - 9.5, i, format(value[i, j], '.1f'),
                               ha="center", va="center", color='black',
                               fontsize=3.2,
                               rotation='vertical')  # color=textcolors[heatmap[j].norm(value[i, j]) > threshold]
    x_label = ['C$_{'+str(_+1)+'}$' for _ in range(len(x_label))]
    ax.set_xlim(-10, 10)
    ax.set_xticks(np.arange(-9.5, 10.5, 1))
    ax.set_xticklabels(x_label)
    ax.set_yticks(range(N))
    ax.set_yticklabels(df.index, horizontalalignment='right', verticalalignment='center', x=-0.04)
    ax.tick_params(axis='y', which='major', labelsize=7, pad=-4)
    ax.tick_params(axis='x', which='major', labelsize=6, pad=1, length=2)
    ax.invert_yaxis()
    fig.tight_layout()
    plt.title('Pair sgRNA-target library', fontweight='bold', fontsize=10, pad=0, y=1.03)
    plt.savefig(path + 'competency.svg', dpi=300)
    plt.show()

# For extract data for Fig S19. Did not include ploting part, which was done in GraphPad.
def potency_purity_bell_curve():
    options = ['potency','purity']
    purity_rf = pd.read_excel(f'{path}{aio_table}',
                              sheet_name='purity',
                              index_col=0).dropna()
    print(len(purity_rf.drop_duplicates(subset='clinvar_id')))
    index = purity_rf.index
    potency_rf = pd.read_excel(f'{path}{aio_table}',
                               sheet_name='potency',
                               index_col=0).loc[index]
    excelwriter = pd.ExcelWriter(f'{path}/purity_potency_bell_curve/purity_potency_bell_curve-042924.xlsx', engine="openpyxl")

    for num, context in enumerate(contexts):
        fig, axs = plt.subplots(2, 1, figsize=(5, 4), constrained_layout=True)
        potency_ot, purity_ot, site_ot = pd.DataFrame(index=purity_rf.columns.to_list()[12:], columns=list(range(1,21,1))), \
                                         pd.DataFrame(index=purity_rf.columns.to_list()[12:], columns=list(range(1,21,1))), \
                                         pd.DataFrame(index=purity_rf.columns.to_list()[12:], columns=list(range(1,21,1)))
        editor_list = control + [context2editor[context]]

        for pos in range(1, 21, 1):
            purity_series = purity_rf[(purity_rf['position'] == pos) & (purity_rf['target_context'] == context)]
            potency_series = potency_rf[(potency_rf['position'] == pos) & (potency_rf['target_context'] == context)]
            for editor in editor_list:
                if editor in control:
                    purity, potency, site_num = purity_series.loc[:, editor].mean(), potency_series.loc[:,
                                                                                     editor].mean(), len(
                        purity_series.loc[:, editor])
                else:
                    purity, potency, site_num = purity_series.loc[:, editor].mean(), potency_series.loc[:,
                                                                                     editor].mean(), len(
                        purity_series.loc[:, editor])
                potency_ot.loc[editor, pos] = potency*100
                purity_ot.loc[editor, pos] = purity
                site_ot.loc[editor, pos] = site_num
        site_ot.dropna().to_excel(excelwriter, f'{context}_site')
        for ax_num,data in enumerate([potency_ot.dropna(),purity_ot.dropna()]):
            data=data.T
            ax = axs[ax_num]
            data.to_excel(excelwriter, f'{context}_{options[ax_num]}')
            data.plot(ax=ax,linewidth=1)
            ax.legend(data.columns, bbox_to_anchor=(1.0, 1.05), prop={'size': 6}, loc=2)
            ax.set_xlim(1, 20)
            ax.set_xticks(list(range(1, 21, 1)))
            if ax_num == 0:
                ax.set_yticks([_ * 10.0 for _ in range(11)])
                ax.set_ylabel('C:G-to-T:A conversion (%)',fontsize=7)
            elif ax_num == 1:
                ax.set_yticks([_ / 10.0 for _ in range(11)])
                ax.set_ylabel('Target C editing/Total C editing', fontsize=7)
            ax.tick_params(axis='y', which='major', labelsize=7)
            # ax.set_xlabel('C:G-to-T:A conversion (%)', fontsize=7)
            # # ax.set_xlim(0,0.8)
            # # ax.set_ylim(0.2,0.8)
            # ax.set_ylabel('Target C editing/\nTotal C editing', fontsize=7)
            # ax.set_title(context, fontsize=8)
        plt.title(context, fontweight='bold', fontsize=10, pad=0, y=1.03)
        plt.savefig(f'{path}/purity_potency_bell_curve/purity_potency_bell_curve_{context}.svg', dpi=600)
        plt.show()
    excelwriter.close()
    return

# Not used in the paper.
def purity_vs_potency_targetC_context_dissect():
    def label_point(x, y, val, ax, con):
        a = pd.concat({'x': x, 'y': y, 'val': val}, axis=1)
        for i, point in a.iterrows():
            # if str(point['val']).split('(')[-1][:3] == con:
            #     ax.text(point['x'] + .005, point['y'], str(point['val']),fontsize=10,weight='bold')
            # else:
            ax.text(point['x'] - .02, point['y'], str(point['val']),fontsize=6)

    purity_rf = pd.read_excel(f'{path}{aio_table}',
                              sheet_name='purity',
                              index_col=0).dropna()
    purity_rf = purity_rf[purity_rf['edit_site_in_proto'] == purity_rf['position']]
    print(len(purity_rf.drop_duplicates(subset=['clinvar_id'])))
    index = purity_rf.index
    potency_rf = pd.read_excel(f'{path}{aio_table}',
                               sheet_name='potency',
                               index_col=0).loc[index]
    excelwriter = pd.ExcelWriter(f'{path}potency_vs_purity_designate_5-7_TargetC-042924.xlsx')

    fig, axs = plt.subplots(4, 4, figsize=(8,6), constrained_layout=True,tight_layout=False)
    total_potency, total_purity, total_site = pd.DataFrame(index=control+['Tad-CBEs(NCN)'], columns=['potency']).fillna(0), \
                                              pd.DataFrame(index=control+['Tad-CBEs(NCN)'], columns=['purity']).fillna(0), \
                                              pd.DataFrame(index=control+['Tad-CBEs(NCN)'], columns=['site_num']).fillna(0)
    for num, context in enumerate(contexts):
        editor_list = control + [context2editor[context]]
        potency_ot, purity_ot, site_ot = pd.DataFrame(index=editor_list, columns=['potency']).fillna(0), \
                                         pd.DataFrame(index=editor_list, columns=['purity']).fillna(0), \
                                         pd.DataFrame(index=editor_list, columns=['site_num']).fillna(0)
        ax = axs[num // 4, num % 4]
        for pos in range(5,8,1):
            purity_series = purity_rf[(purity_rf['position'] == pos) & (purity_rf['target_context'] == context)]
            potency_series = potency_rf[(potency_rf['position'] == pos) & (potency_rf['target_context'] == context)]

            for editor_index,editor in enumerate(editor_list):
                purity, potency, site_num = purity_series.loc[:,editor].sum(), potency_series.loc[:,editor].sum(), len(purity_series.loc[:,editor])
                potency_ot.loc[editor,'potency'] += potency
                purity_ot.loc[editor,'purity'] += purity
                site_ot.loc[editor,'site_num'] += site_num
                total_potency.iloc[editor_index,0] += potency
                total_purity.iloc[editor_index, 0] += purity
                total_site.iloc[editor_index, 0] += site_num
        potency_ot=potency_ot/site_ot['site_num'].mean()
        purity_ot=purity_ot/site_ot['site_num'].mean()
        data = pd.concat([potency_ot,purity_ot,site_ot],axis=1).dropna()
        data.to_excel(excelwriter,f'5-7_{context}')
        sns.scatterplot(x='potency',y='purity',data=data,s=6,ax=ax)
        ax.set_xlabel('C:G-to-T:A conversion (%)',fontsize=7)
        # ax.set_xlim(0,0.8)
        # ax.set_ylim(0.2,0.8)
        ax.set_ylabel('Target C editing/\nTotal C editing',fontsize=7)
        ax.set_title(context,fontsize=8)
        data['editor']=data.index.to_list()
        label_point(data['potency'],data['purity'],data['editor'],ax,context)
    data = pd.concat([total_potency, total_purity, total_site], axis=1).dropna().to_excel(excelwriter,f'5-7_Tad-CBEs(NCN)')

    plt.savefig(f'{path}potency_vs_purity_designate_5-7_TargetC.svg',dpi=600)
    plt.show()
    excelwriter.close()
    return

# For Fig S21 B
def purity_vs_potency_targetC_position_combined():
    def label_point(x, y, val, ax):
        a = pd.concat({'x': x, 'y': y, 'val': val}, axis=1)
        for i, point in a.iterrows():
            # if str(point['val']).split('(')[-1][:3] == con:
            #     ax.text(point['x'] + .005, point['y'], str(point['val']),fontsize=10,weight='bold')
            # else:
            ax.text(point['x'] - .02, point['y'], str(point['val']),fontsize=6)

    purity_rf = pd.read_excel(f'{path}{aio_table}',
                              sheet_name='purity',
                              index_col=0).dropna()
    purity_rf = purity_rf[purity_rf['edit_site_in_proto'] == purity_rf['position']]
    print(len(purity_rf.drop_duplicates(subset=['clinvar_id'])))
    index = purity_rf.index
    potency_rf = pd.read_excel(f'{path}{aio_table}',
                               sheet_name='potency',
                               index_col=0).loc[index]

    fig, axs = plt.subplots(4, 4, figsize=(7.7,6.6), constrained_layout=True,tight_layout=False)
    total_potency, total_purity, total_site = pd.DataFrame(index=control+['Tad-CBEs(NCN)'], columns=['potency']).fillna(0), \
                                              pd.DataFrame(index=control+['Tad-CBEs(NCN)'], columns=['purity']).fillna(0), \
                                              pd.DataFrame(index=control+['Tad-CBEs(NCN)'], columns=['site_num']).fillna(0)
    for num, context in enumerate(contexts):
        editor_list = control + [context2editor[context]]
        potency_ot, purity_ot, site_ot = pd.DataFrame(index=editor_list, columns=['potency']).fillna(0), \
                                         pd.DataFrame(index=editor_list, columns=['purity']).fillna(0), \
                                         pd.DataFrame(index=editor_list, columns=['site_num']).fillna(0)
        for pos in range(5,8,1):
            purity_series = purity_rf[(purity_rf['position'] == pos) & (purity_rf['target_context'] == context)]
            potency_series = potency_rf[(potency_rf['position'] == pos) & (potency_rf['target_context'] == context)]

            for editor_index,editor in enumerate(editor_list):
                purity, potency, site_num = purity_series.loc[:,editor].sum(), potency_series.loc[:,editor].sum(), len(purity_series.loc[:,editor])
                potency_ot.loc[editor,'potency'] += potency
                purity_ot.loc[editor,'purity'] += purity
                site_ot.loc[editor,'site_num'] += site_num
                total_potency.iloc[editor_index,0] += potency
                total_purity.iloc[editor_index, 0] += purity
                total_site.iloc[editor_index, 0] += site_num
        potency_ot=potency_ot/site_ot['site_num'].mean()
        purity_ot=purity_ot/site_ot['site_num'].mean()
        data = pd.concat([potency_ot,purity_ot,site_ot],axis=1).dropna()
        data['editor']=data.index.to_list()
    data = pd.concat([total_potency/total_site['site_num'].mean(), total_purity/total_site['site_num'].mean(), total_site], axis=1).dropna()
    ax = axs[0,0]
    for ax in axs.reshape(-1):
        sns.scatterplot(x='potency', y='purity', data=data, s=6, ax=ax)
        ax.set_xlabel('C:G-to-T:A conversion (%)', fontsize=7)
        # ax.set_xlim(0,0.8)
        # ax.set_ylim(0.2,0.8)
        data['editor'] = data.index.to_list()
        label_point(data['potency'], data['purity'], data['editor'], ax)
        ax.set_ylabel('Target C editing/\nTotal C editing', fontsize=7)
        ax.set_title('Position 5-7', fontsize=8)
    plt.savefig(f'{path}potency_vs_purity_designate_5-7_combined_TargetC.svg',dpi=600)
    plt.show()
    return

# For Fig S21 A
def purity_vs_potency_targetC_position_dissect():
    def label_point(x, y, val, ax):
        a = pd.concat({'x': x, 'y': y, 'val': val}, axis=1)
        for i, point in a.iterrows():
            # if str(point['val']).split('(')[-1][:3] == con:
            #     ax.text(point['x'] + .005, point['y'], str(point['val']),fontsize=10,weight='bold')
            # else:
            ax.text(point['x'] - .02, point['y'], str(point['val']),fontsize=6)

    purity_rf = pd.read_excel(f'{path}010824_aio_table_041424_TargetC_unique.xlsx',
                              sheet_name='purity',
                              index_col=0).dropna().loc[rf.index]
    print(len(purity_rf.drop_duplicates(subset='clinvar_id')))
    index = purity_rf.index
    potency_rf = pd.read_excel(f'{path}010824_aio_table_041424_TargetC_unique.xlsx',
                               sheet_name='potency',
                               index_col=0).loc[index]
    excelwriter = pd.ExcelWriter(f'{path}potency_vs_purity_1-12_TargetC-042924.xlsx')

    fig, axs = plt.subplots(3, 4, figsize=(8, 5), constrained_layout=True, tight_layout=False)
    total_potency, total_purity, total_site = pd.DataFrame(index=control + ['Tad-CBEs(NCN)'],
                                                           columns=list(range(1,13,1))).fillna(0), \
                                              pd.DataFrame(index=control + ['Tad-CBEs(NCN)'],
                                                           columns=list(range(1,13,1))).fillna(0), \
                                              pd.DataFrame(index=control + ['Tad-CBEs(NCN)'],
                                                           columns=list(range(1,13,1))).fillna(0)
    purity_rf = purity_rf[purity_rf['edit_site_in_proto'] == purity_rf['position']]
    potency_rf = potency_rf[potency_rf['edit_site_in_proto'] == potency_rf['position']]
    for pos in range(1, 13, 1):
        potency_ot, purity_ot, site_ot = pd.DataFrame(index=control + ['Tad-CBEs(NCN)'], columns=['potency']).fillna(0), \
                                         pd.DataFrame(index=control + ['Tad-CBEs(NCN)'], columns=['purity']).fillna(0), \
                                         pd.DataFrame(index=control + ['Tad-CBEs(NCN)'], columns=['site_num']).fillna(0)
        ax = axs[(pos-1)//4,(pos-1) % 4]
        purity_series = purity_rf[purity_rf['position'] == pos]
        potency_series = potency_rf[potency_rf['position'] == pos]
        for editor_index, editor in enumerate(control + ['Tad-CBEs(NCN)']):
            purity, potency, site_num = purity_series.loc[:, editor].mean(), potency_series.loc[:,
                                                                            editor].mean(), len(
                purity_series.loc[:, editor])
            potency_ot.loc[editor, 'potency'] += potency
            purity_ot.loc[editor, 'purity'] += purity
            site_ot.loc[editor, 'site_num'] += site_num
            total_potency.loc[editor, pos] += potency
            total_purity.loc[editor, pos] += purity
            total_site.loc[editor, pos] += site_num
        data = pd.concat([potency_ot, purity_ot, site_ot], axis=1).dropna()
        sns.scatterplot(x='potency', y='purity', data=data, s=6, ax=ax)
        ax.set_xlabel('C:G-to-T:A conversion (%)', fontsize=7)
        # ax.set_xlim(0,0.8)
        # ax.set_ylim(0.2,0.8)
        ax.set_ylabel('Target C editing/\nTotal C editing', fontsize=7)
        ax.set_title(f'Position {pos}', fontsize=8)
        data['editor'] = data.index.to_list()
        label_point(data['potency'], data['purity'], data['editor'], ax)
    total_potency.to_excel(excelwriter,'potency')
    total_purity.to_excel(excelwriter,'purity')
    total_site.to_excel(excelwriter,'site_number')
    plt.savefig(f'{path}potency_vs_purity_designate_{pos}_TargetC.svg', dpi=600)
    plt.show()
    excelwriter.close()
    return

# Not used in the paper.
def targetC_lib_winning_match():
    purity_rf = pd.read_excel(f'{path}010824_aio_table_041424_TargetC_unique.xlsx',
                              sheet_name='purity',
                              index_col=0).dropna().loc[rf.index]
    print(len(purity_rf.drop_duplicates(subset='clinvar_id')))
    index = purity_rf.index
    potency_rf = pd.read_excel(f'{path}010824_aio_table_041424_TargetC_unique.xlsx',
                               sheet_name='potency',
                               index_col=0).loc[index]
    purity_rf = purity_rf[purity_rf['edit_site_in_proto'] == purity_rf['position']]
    potency_rf = potency_rf[potency_rf['edit_site_in_proto'] == potency_rf['position']]
    purity_series = purity_rf[(purity_rf['Tad-CBEs(NCN)'] > purity_rf['Tad-CBE2.4']) & (
            purity_rf['Tad-CBEs(NCN)'] > purity_rf['Tad-CBE3.1']) & (
                                    purity_rf['Tad-CBEs(NCN)'] > purity_rf['TadCBEd']) & (
                                    purity_rf['Tad-CBEs(NCN)'] > purity_rf['rAPOBEC1-BE4']) & (
                                    purity_rf['edit_site_in_proto'] == purity_rf['position'])]
    potency_series = potency_rf[(potency_rf['Tad-CBEs(NCN)'] > potency_rf['Tad-CBE2.4']) & (
            potency_rf['Tad-CBEs(NCN)'] > potency_rf['Tad-CBE3.1']) & (
                                      potency_rf['Tad-CBEs(NCN)'] > potency_rf['TadCBEd']) & (
                                      potency_rf['Tad-CBEs(NCN)'] > potency_rf['rAPOBEC1-BE4']) & (
                                      potency_rf['edit_site_in_proto'] == potency_rf['position'])]
    ot_df = pd.DataFrame(index=['potency','purity','sweetspot_NCN'],columns=list(range(1,13,1))).fillna(0)
    for pos in range(1,13,1):
        ot_df.loc['potency',pos]=potency_series[potency_series['position'] == pos]['Tad-CBEs(NCN)'].mean()
        ot_df.loc['purity', pos] = purity_series[purity_series['position'] == pos]['Tad-CBEs(NCN)'].mean()
        ot_df.loc['sweetspot_NCN', pos] = len(purity_series[purity_series['position'] == pos])
    ot_df.to_excel(f'{path}targetC_unique_sweetspot.xlsx')
    return

# For extract data for Fig S20. Did not include ploting part, which was done in GraphPad.
def targetC_lib_winning_no_match():
    excelwriter = pd.ExcelWriter(f'{path}winning_no_match.xlsx')
    purity_rf = pd.read_excel(f'{path}010824_aio_table_041424_C_unique.xlsx',
                              sheet_name='purity',
                              index_col=0).dropna().loc[rf.index]
    print(len(purity_rf.drop_duplicates(subset='clinvar_id')))
    index = purity_rf.index
    potency_rf = pd.read_excel(f'{path}010824_aio_table_041424_C_unique.xlsx',
                               sheet_name='potency',
                               index_col=0).loc[index]
    purity_rf = purity_rf[purity_rf['edit_site_in_proto'] == purity_rf['position']]
    potency_rf = potency_rf[potency_rf['edit_site_in_proto'] == potency_rf['position']]
    for (df,option) in zip([purity_rf,potency_rf],['purity','potency']):
        fulltable = {}
        for context in contexts:
            fulltable[context] = pd.DataFrame(index=df.columns.to_list()[12:], columns=range(1, 21, 1)).fillna(0)
        ot_df = pd.DataFrame(index=df.columns.to_list()[12:],columns=range(1,21,1)).fillna(0)

        df['max_purity'] = df.iloc[:, 12:].apply(max, axis=1)
        df['max_purity_editor'] = [[] for _ in range(len(df['max_purity']))]
        df['group'] = [0 for _ in range(len(df['max_purity']))]
        if option == 'purity':
            df['best_editor'] = ['' for _ in range(len(df['max_purity']))]
            df['potency_of_the_best_editor'] = [0 for _ in range(len(df['max_purity']))]
        site_df = pd.DataFrame(index=contexts,columns=range(1,21,1)).fillna(0)
        win_groups = [0,0,0] # 1.only control win 2.both NCN and control win. 3.only Tad-CBEs(NCN) win
        for _ in range(len(df)):
            site_df.loc[df['context'][_], df['edit_site_in_proto'][_]] += 1  # record site number
            for editor in df.columns[12:32]:
                if df.iloc[_].loc[editor] == df['max_purity'][_]:
                    ot_df.loc[editor, int(df['edit_site_in_proto'][_])] += 1
                    df['max_purity_editor'][_].append(editor)
                    fulltable[df['context'][_]].loc[editor, int(df['edit_site_in_proto'][_])] += 1
            assign = False
            #group sites by winning situations
            for editor in control:
                if editor in df['max_purity_editor'][_]:
                    for winner in df['max_purity_editor'][_]:
                        if not(winner in control):
                            win_groups[1] += 1
                            df['group'][_] = 2
                            assign = True
                            break
                    if assign == False and editor == df['max_purity_editor'][_][-1]:
                        win_groups[0] += 1
                        assign = True
                        df['group'][_] = 1
                        break
                    elif assign == True:
                        break
            if assign == False:
                win_groups[2] += 1
                df['group'][_] = 3
            #find potency of winning editors
            if option == 'purity':
                potency_of_max_purity = [potency_rf.loc[purity_rf.index[_],editor] for editor in purity_rf['max_purity_editor'][_]]
                df['best_editor'][_] = df['max_purity_editor'][_][np.argmax(potency_of_max_purity)]
                df['potency_of_the_best_editor'][_] = potency_of_max_purity[np.argmax(potency_of_max_purity)]
        print(win_groups)
        df.to_excel(excelwriter, f'unmatch_{option}')
        site_df.to_excel(excelwriter, f'unmatch_site_number_{option}')
        ot_df.to_excel(excelwriter, f'unmatch_winning_{option}')
        if option == 'purity':
            print([len(df[(df['position'] == pos) & ((df['group']==3)|(df['group']==2))]) for pos in range(1,13,1)])
        # for key in fulltable.keys():
        #     fulltable[key].to_excel(excelwriter, f'{key}_no_match_Tad-CBEs(NCN)')
    excelwriter.close()
    return

# For extract data for Fig S20. Did not include ploting part, which was done in GraphPad.
def targetC_lib_winning_no_match_position():
    def label_point(x, y, val, ax):
        a = pd.concat({'x': x, 'y': y, 'val': val}, axis=1)
        for i, point in a.iterrows():
            # if str(point['val']).split('(')[-1][:3] == con:
            #     ax.text(point['x'] + .005, point['y'], str(point['val']),fontsize=10,weight='bold')
            # else:
            ax.text(point['x'] - .02, point['y'], str(point['val']),fontsize=6)

    excelwriter = pd.ExcelWriter(f'{path}best-perform-Tad-CBEs(NCN)-winning_no_match_test.xlsx')
    purity_rf = pd.read_excel(f'{path}010824_aio_table_041424_C_unique.xlsx',
                              sheet_name='purity',
                              index_col=0).dropna().loc[rf.index]
    print(len(purity_rf.drop_duplicates(subset='clinvar_id')))
    index = purity_rf.index
    potency_rf = pd.read_excel(f'{path}010824_aio_table_041424_C_unique.xlsx',
                               sheet_name='potency',
                               index_col=0).loc[index]
    purity_rf = purity_rf[purity_rf['edit_site_in_proto'] == purity_rf['position']]
    potency_rf = potency_rf[potency_rf['edit_site_in_proto'] == potency_rf['position']]
    for (df,option) in zip([purity_rf,potency_rf],['purity','potency']):
        fulltable = {}
        for context in contexts:
            fulltable[context] = pd.DataFrame(index=df.columns.to_list()[12:], columns=range(1, 21, 1)).fillna(0)
        ot_df = pd.DataFrame(index=df.columns.to_list()[12:],columns=range(1,21,1)).fillna(0)
        if option == 'purity':
            df['max_purity'] = df.iloc[:, 16:].apply(max, axis=1)
        elif option == 'potency':
            df['max_purity'] = df.iloc[:, list(range(16,29))+[30]].apply(max, axis=1)
        df['max_purity_editor'] = [[] for _ in range(len(df['max_purity']))]
        df['group'] = [0 for _ in range(len(df['max_purity']))]
        if option == 'purity':
            df['best_editor'] = ['' for _ in range(len(df['max_purity']))]
            df['potency_of_the_best_editor'] = [0 for _ in range(len(df['max_purity']))]
        site_df = pd.DataFrame(index=contexts,columns=range(1,21,1)).fillna(0)
        win_groups = [0,0,0] # 1.only control win 2.both NCN and control win. 3.only Tad-CBEs(NCN) win
        for _ in range(len(df)):
            site_df.loc[df['context'][_], df['edit_site_in_proto'][_]] += 1  # record site number
            if option == 'purity':
                for editor in df.columns[16:32]:
                    if df.iloc[_].loc[editor] == df['max_purity'][_]:
                        ot_df.loc[editor, int(df['edit_site_in_proto'][_])] += 1
                        df['max_purity_editor'][_].append(editor)
                        fulltable[df['context'][_]].loc[editor, int(df['edit_site_in_proto'][_])] += 1
            elif option == 'potency':
                for editor in df.columns[list(range(16,29))+[30]]:
                    if df.iloc[_].loc[editor] == df['max_purity'][_]:
                        ot_df.loc[editor, int(df['edit_site_in_proto'][_])] += 1
                        df['max_purity_editor'][_].append(editor)
                        fulltable[df['context'][_]].loc[editor, int(df['edit_site_in_proto'][_])] += 1
            assign = False
            #group sites by winning situations
            for editor in control:
                if editor in df['max_purity_editor'][_]:
                    for winner in df['max_purity_editor'][_]:
                        if not(winner in control):
                            win_groups[1] += 1
                            df['group'][_] = 2
                            assign = True
                            break
                    if assign == False and editor == df['max_purity_editor'][_][-1]:
                        win_groups[0] += 1
                        assign = True
                        df['group'][_] = 1
                        break
                    elif assign == True:
                        break
            if assign == False:
                win_groups[2] += 1
                df['group'][_] = 3
            #find potency of winning editors
            if option == 'purity':
                potency_of_max_purity = [potency_rf.loc[purity_rf.index[_],editor] for editor in purity_rf['max_purity_editor'][_]]
                df['best_editor'][_] = df['max_purity_editor'][_][np.argmax(potency_of_max_purity)]
                df['potency_of_the_best_editor'][_] = potency_of_max_purity[np.argmax(potency_of_max_purity)]
        print(win_groups)
        df.to_excel(excelwriter, f'unmatch_{option}')
        site_df.to_excel(excelwriter, f'unmatch_site_number_{option}')
        ot_df.to_excel(excelwriter, f'unmatch_winning_{option}')
        # df[(df['position']==5)|(df['position']==6)|(df['position']==7)].iloc[:,[12,13,14,15,32,36]].mean() # best-perform Tad-CBEs(NCN) potency/purity
        # potency_rf.loc[df[(df['position']==5)|(df['position']==6)|(df['position']==7)].index].iloc[:,[12,13,14,15]].mean() #control potency

        # potency vs. purity plot
        if option == 'purity':
            df_series = df[df['position']==6]
            data = pd.DataFrame(index=control+['Tad-CBEs(NCN)'],columns=['potency','purity']).fillna(0)
            for num in [12,14,13,15]:
                data.iloc[num-12,1] = df_series.iloc[:,num].mean()
                data.iloc[num-12,0] = potency_rf.loc[df_series.index].iloc[:,
                                         num].mean()
            data.loc['Tad-CBEs(NCN)','purity'] = df_series.iloc[:,32].mean()
            data.loc['Tad-CBEs(NCN)','potency'] = df_series.iloc[:,36].mean()
            data['editor'] = data.index.to_list()
            fig, axs = plt.subplots(4, 4, figsize=(7.7, 6.6), constrained_layout=True, tight_layout=False)
            for ax in axs.reshape(-1):
                sns.scatterplot(x='potency', y='purity', data=data, s=6, ax=ax)
                ax.set_xlabel('C:G-to-T:A conversion (%)', fontsize=7)
                # ax.set_xlim(0,0.8)
                # ax.set_ylim(0.2,0.8)
                data['editor'] = data.index.to_list()
                label_point(data['potency'], data['purity'], data['editor'], ax)
                ax.set_ylabel('Target C editing/\nTotal C editing', fontsize=7)
                ax.set_title('Position 6', fontsize=8)
            plt.savefig(f'{path}potency_vs_purity_designate_6_combined_TargetC.svg', dpi=600)
            plt.show()
        # for key in fulltable.keys():
        #     fulltable[key].to_excel(excelwriter, f'{key}_no_match_Tad-CBEs(NCN)')
    excelwriter.close()
    return

def making_excel_output_from_list():
    dfs_list = {}
    ot_excel = []
    df = pd.read_excel(f'{path}best-perform-Tad-CBEs(NCN)-winning_no_match_test.xlsx',sheet_name='unmatch_purity',index_col=0)
    for _ in range(len(df['gRNA'])):
        editor = df['best_editor'][_]
        if editor in dfs_list.keys():
            excel_df = dfs_list[editor]
        else:
            excel_df = pd.read_excel(f'/Users/yuanwu/BaseSpace/010824-nextseq/excel_output/{editor2file[editor]}',index_col=False)
            dfs_list[editor] = excel_df.copy(deep='True')
        append_series = excel_df[excel_df['clinvar_id'] == str(df['clinvar_id'][_])]
        append_series['purity'] = df['max_purity'][_]
        ot_excel.append(append_series)
    pd.concat(ot_excel).to_excel(f'/Users/yuanwu/BaseSpace/010824-nextseq/excel_output/Tad-CBEs(NCN).xlsx')
    return

if __name__ == '__main__':
    # FOR aio-heatmap
    # for num,file in enumerate(['Tad-CBEs(NCN)']):
    #     for pos in range(5,8,1):
    #         for bool in [False]:
    #             # if num <= 3:
    #             #     continue
    #             context = file.split('(')[1][:3]
    #             aio_heatmap_Tad_CBEs(pos,[file]+control,'',bool) # ploting all-in-one heatmap
