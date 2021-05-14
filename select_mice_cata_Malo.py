from configuration import *
# os.system('pwd')
# os.system('which python')
# os.system('python --version')



# def get_mice(group):
#     print(data_dir)
#     print(os.listdir(data_dir))
#     path = '{}{}/'.format(data_dir, group)
#     print(path)
#     files = os.listdir(path)
#     if '.DS_Store' in files:
#         files.remove('.DS_Store')
#
#     return files
#
#
# def get_mice_for_spectrum(group):
#     path = '{}{}/'.format(data_spectrum_dir, group)
#     files = os.listdir(path)
#     if '.DS_Store' in files:
#         files.remove('.DS_Store')
#     files = [f[4:] for f in files]
#     return files

def get_mice(miRNA, genotype, exclude_mice=False):
    df_excel = pd.read_excel(work_dir + 'datetime_reference_miRNA.xlsx', index_col = 0)
    df_excel['miRNA'] = df_excel['miRNA'].astype('str')
    mylist = df_excel[(df_excel.genotype == genotype) & (df_excel.miRNA == miRNA)].index.to_list()
    files = [i[4:] for i in mylist]
    if exclude_mice:
        mice = []
        for mouse in files:
            if mouse not in RT_PCR_execption:
                mice.append(mouse)
    else :
        mice = files
    return mice

def get_mouse_info(mouse):
    df_excel = pd.read_excel(work_dir + 'datetime_reference_miRNA.xlsx', index_col = 0)
    miRNA = str(df_excel.at['MTA-{}'.format(mouse), 'miRNA'])
    genotype = df_excel.at['MTA-{}'.format(mouse), 'genotype']
    return miRNA, genotype



if __name__ == '__main__':

    # get_mouse_info('B2884')
    # print(get_mice('137', 'test'))
    print(get_mice('137', 'test'))
    print(get_mice('137', 'test', True))
    print(RT_PCR_execption)
    # print(get_mice_for_spectrum('Control'))
    # print(get_mice_for_spectrum('DCR-HCRT'))
