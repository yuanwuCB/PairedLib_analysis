context2editor = {'CCG':'Tad-CBE(CCG3)','CCT':'Tad-CBE(CCT2)','GCT':'Tad-CBE(GCT1)',
                  'GCC':'Tad-CBE(GCC3)','TCC':'Tad-CBE(TCC1)','CCC':'Tad-CBE(CCC5)',
                  'CCA':'Tad-CBE(CCA2)','ACG':'Tad-CBE(ACG1)' ,'ACA':'Tad-CBE(ACA1)',
                  'TCT':'Tad-CBE(TCT1)','TCG':'Tad-CBE(TCG2)','GCA':'Tad-CBE(GCA1)',
                  'GCG':'Tad-CBE(GCG1)','TCA':'Tad-CBE(TCA3)','ACC':'Tad-CBE(ACC3)',
                  'ACT':'Tad-CBE(ACT4)'}

# editor_order = ['Tad-CBE2.4','Tad-CBE3.1','TadCBEd','rAPOBEC1-BE4','Tad-CBEs(NCN)']

editor_order = ['Tad-CBE2.4','Tad-CBE3.1','TadCBEd','rAPOBEC1-BE4',
                'Tad-CBE(ACA1)','Tad-CBE(ACT4)','Tad-CBE(ACC3)','Tad-CBE(ACG1)',
                'Tad-CBE(TCA3)','Tad-CBE(TCT1)','Tad-CBE(TCC1)','Tad-CBE(TCG2)',
                'Tad-CBE(CCA2)','Tad-CBE(CCT2)','Tad-CBE(CCC5)','Tad-CBE(CCG3)',
                'Tad-CBE(GCA1)','Tad-CBE(GCT1)','Tad-CBE(GCC3)','Tad-CBE(GCG1)']

contexts = ['ACA','ACT','ACC','ACG',
            'TCA','TCT','TCC','TCG',
            'CCA','CCT','CCC','CCG',
            'GCA','GCT','GCC','GCG']

editor2file = {'Tad-CBE2.4':'050e_228lib_S1_L001.xlsx','Tad-CBE3.1':'156e_228lib_S3_L001.xlsx',
				'TadCBEd':'050h_228lib_S2_L001.xlsx','rAPOBEC1-BE4':'BE4_228lib_S4_L001.xlsx',
                'Tad-CBE(ACA1)':'pYX203ua_228lib_S5_L001.xlsx','Tad-CBE(ACT4)':'pYX217d_228lib_S6_L001.xlsx',
                'Tad-CBE(ACC3)':'pyx220c_228lib_S9_L001.xlsx','Tad-CBE(ACG1)':'pyx165f_228lib_S7_L001.xlsx',
                'Tad-CBE(TCA3)':'pyx208uc_228lib_S10_L001.xlsx','Tad-CBE(TCT1)':'pYX204ua_228lib_S11_L001.xlsx',
                'Tad-CBE(TCC1)':'pYW097ua_228lib_S13_L001.xlsx','Tad-CBE(TCG2)':'pyx205ub_228lib_S12_L001.xlsx',
                'Tad-CBE(CCA2)':'pYW106uc_228lib_S22_L001.xlsx','Tad-CBE(CCT2)':'pYW094ub_228lib_S23_L001.xlsx',
                'Tad-CBE(CCC5)':'pYW105ue_228lib_S26_L001.xlsx','Tad-CBE(CCG3)':'pYW093uc_228lib_S24_L001.xlsx',
                'Tad-CBE(GCA1)':'pYX206ua_228lib_S15_L001.xlsx','Tad-CBE(GCT1)':'pYW095ua_228lib_S16_L001.xlsx',
                'Tad-CBE(GCC3)':'pYW096uc_228lib_S20_L001.xlsx','Tad-CBE(GCG1)':'pYX207ua_228lib_S18_L001.xlsx',
                'Tad-CBEs(NCN)':'Tad-CBEs(NCN).xlsx'}

# rf = pd.read_excel(f'{path}{aio_table}', sheet_name='purity',
#                    index_col=0).dropna()

control = ['Tad-CBE2.4','Tad-CBE3.1','TadCBEd','rAPOBEC1-BE4']

