from configuration import *
from select_mice_cata_Malo import get_mice

def average_baseline_spectrum(spectrum_method = 'welch', exclude_mice = False, auc = False, agg_power = 'mean'):
    for miRNA in ['128', '137', '665']:
        groups = {'control' : get_mice(miRNA, 'control'), 'test' : get_mice(miRNA, 'test')}
        for genotype in groups:
            if exclude_mice:
                dirname = excel_dir + '/RT_PCR_excluded_mice/power_spectrum_baseline_auc_{}_agg_power_{}/{}/{}/'.format(auc, agg_power,miRNA,genotype)
            else:
                dirname = excel_dir + '/power_spectrum_baseline_auc_{}_agg_power_{}/{}/{}/'.format(auc, agg_power,miRNA,genotype)
            states = ['w', 'n', 'r','t']
            for state in states :
                period = 'light'
                df1 = pd.read_excel(dirname+'/{}_spectrum_{}_{}_bl1_{}_exclude_mice_{}.xlsx'.format(miRNA,spectrum_method, state, period,exclude_mice), index_col =0)
                df2 = pd.read_excel(dirname+'/{}_spectrum_{}_{}_bl2_{}_exclude_mice_{}.xlsx'.format(miRNA,spectrum_method, state, period,exclude_mice), index_col =0)
                mice1 = df1.index
                mice2 = df2.index
                freqs = df1.columns[4:].to_list()
                df = pd.DataFrame(index = get_mice(miRNA, genotype), columns = freqs )
                for ind in range(mice1.size) :
                    mouse1 = mice1[ind]
                    mouse2 = mice2[ind]
                    mouse = mouse1[:5]
                    values1 = df1.loc[mouse1, freqs].values
                    values2 = df2.loc[mouse2, freqs].values
                    if sum(np.isnan(np.array(values1, dtype = 'float64'))) == values1.size  or sum(np.isnan(np.array(values2, dtype = 'float64'))) == values2.size :
                        weighted_mean = np.nansum((values1,values2), axis = 0)
                    else :
                        weighted_mean = (12/24)*values1 + (12/24)*values2
                    df.loc[mouse,freqs] = weighted_mean
                df.to_excel(dirname + '/{}_{}_spectrum_{}_{}_BASELINE_{}_exclude_mice_{}.xlsx'.format(miRNA,genotype, spectrum_method, state, period, exclude_mice))
                #
                period = 'dark'
                df1 = pd.read_excel(dirname+'/{}_spectrum_{}_{}_bl1_{}_exclude_mice_{}.xlsx'.format(miRNA,spectrum_method, state, period, exclude_mice), index_col =0)
                df2 = pd.read_excel(dirname+'/{}_spectrum_{}_{}_bl2_{}_exclude_mice_{}.xlsx'.format(miRNA,spectrum_method, state, period, exclude_mice), index_col =0)
                mice1 = df1.index
                mice2 = df2.index
                freqs = df1.columns[4:].to_list()
                df = pd.DataFrame(index = get_mice(miRNA, genotype), columns = ['ref_values'] + freqs )
                for ind in range(mice1.size) :
                    mouse1 = mice1[ind]
                    mouse2 = mice2[ind]
                    mouse = mouse1[:5]
                    values1 = df1.loc[mouse1, freqs].values
                    values2 = df2.loc[mouse2, freqs].values
                    if sum(np.isnan(np.array(values1, dtype = 'float64'))) == values1.size  or sum(np.isnan(np.array(values2, dtype = 'float64'))) == values2.size :
                        weighted_mean = np.nansum((values1,values2), axis = 0)
                    else :
                        weighted_mean = (12/24)*values1 + (12/24)*values2
                    df.loc[mouse,freqs] = weighted_mean
                df.to_excel(dirname + '/{}_{}_spectrum_{}_{}_BASELINE_{}_exclude_mice_{}.xlsx'.format(miRNA,genotype, spectrum_method, state, period,exclude_mice))

def average_recovery_spectrum(spectrum_method = 'welch', exclude_mice = False, auc = False, agg_power = 'mean'):
    date_ref = pd.read_excel(work_dir + 'datetime_reference_miRNA.xlsx', index_col = 0)
    for miRNA in ['128', '137', '665']:
        groups = {'control' : get_mice(miRNA, 'control'), 'test' : get_mice(miRNA, 'test')}
        for genotype in groups:
            if exclude_mice:
                dirname = excel_dir + '/RT_PCR_excluded_mice/power_spectrum_baseline_auc_{}_agg_power_{}/{}/{}/'.format(auc, agg_power,miRNA,genotype)
            else :
                dirname = excel_dir + '/power_spectrum_baseline_auc_{}_agg_power_{}/{}/{}/'.format(auc, agg_power,miRNA,genotype)
            states = ['w', 'n', 'r', 't']
            for state in states :
                period = 'light'
                df1 = pd.read_excel(dirname+'/{}_spectrum_{}_{}_sd_{}_exclude_mice_{}.xlsx'.format(miRNA, spectrum_method, state, period,exclude_mice), index_col =0)
                df2 = pd.read_excel(dirname+'/{}_spectrum_{}_{}_r1_{}_exclude_mice_{}.xlsx'.format(miRNA, spectrum_method, state, period,exclude_mice), index_col =0)
                mice1 = df1.index
                mice2 = df2.index
                freqs = df1.columns[4:].to_list()
                df = pd.DataFrame(index = get_mice(miRNA, genotype), columns = freqs )
                for ind in range(mice1.size) :
                    mouse1 = mice1[ind]
                    mouse2 = mice2[ind]
                    mouse = mouse1[:5]
                    values1 = df1.loc[mouse1, freqs].values
                    values2 = df2.loc[mouse2, freqs].values
                    if sum(np.isnan(np.array(values1, dtype = 'float64'))) == values1.size  or sum(np.isnan(np.array(values2, dtype = 'float64'))) == values2.size :
                        weighted_mean = np.nansum((values1,values2), axis = 0)
                    else :
                        weighted_mean = (6/18)*values1 + (12/18)*values2
                    df.loc[mouse,freqs] = weighted_mean
                df.to_excel(dirname + '/{}_{}_spectrum_{}_{}_RECOVERY_{}_exclude_mice_{}.xlsx'.format(miRNA, genotype, spectrum_method, state, period,exclude_mice))
                #
                period = 'dark'
                df1 = pd.read_excel(dirname+'/{}_spectrum_{}_{}_sd_{}_exclude_mice_{}.xlsx'.format(miRNA,spectrum_method, state, period,exclude_mice), index_col =0)
                df2 = pd.read_excel(dirname+'/{}_spectrum_{}_{}_r1_{}_exclude_mice_{}.xlsx'.format(miRNA,spectrum_method, state, period,exclude_mice), index_col =0)
                mice1 = df1.index
                mice2 = df2.index
                freqs = df1.columns[4:].to_list()
                df = pd.DataFrame(index = get_mice(miRNA, genotype), columns =  freqs )
                for ind in range(mice1.size) :
                    mouse1 = mice1[ind]
                    mouse2 = mice2[ind]
                    mouse = mouse1[:5]
                    if  date_ref.at['MTA-{}'.format(mouse), 'short_end'] == 'yes':
                        weighted_mean = df1.loc[mouse1, freqs].values
                    else :
                        values1 = df1.loc[mouse1, freqs].values
                        values2 = df2.loc[mouse2, freqs].values
                        # if sum(np.isnan(np.array(values1, dtype = 'float64'))) == values1.size  or sum(np.isnan(np.array(values2, dtype = 'float64'))) == values2.size :
                        #     weighted_mean = np.nansum((values1,values2), axis = 0)
                        weighted_mean = (12/24)*values1 + (12/24)*values2
                    df.loc[mouse,freqs] = weighted_mean

                df.to_excel(dirname + '/{}_{}_spectrum_{}_{}_RECOVERY_{}_exclude_mice_{}.xlsx'.format(miRNA,genotype, spectrum_method, state, period,exclude_mice))
                #
                #

def average_baseline_fragmentation(exclude_mice = False):
    date_ref = pd.read_excel(work_dir + 'datetime_reference_miRNA.xlsx', index_col = 0)
    encoding_to_state = {'w':'wake', 'n':'nrem', 'r':'rem', 't':'tdw'}
    for folder in ['sleep_fragmentation','sleep_percent_fragmentation']:
        for miRNA in ['128', '137', '665']:
            groups = {'control' : get_mice(miRNA, 'control'), 'test' : get_mice(miRNA, 'test')}
            for genotype in groups:
                if exclude_mice:
                    dirname = excel_dir + '/RT_PCR_excluded_mice/{}/{}/{}/'.format(folder,miRNA, genotype)
                else :
                    dirname = excel_dir + '/{}/{}/{}/'.format(folder,miRNA, genotype)
                encoding = ['w', 'n', 'r', 't']
                for e in encoding :
                    state = encoding_to_state[e]
                    period = 'light'
                    df1 = pd.read_excel(dirname+'/{}/{}_{}_{}_{}1.xlsx'.format(state,miRNA, genotype, state, period), index_col =0)
                    df2 = pd.read_excel(dirname+'/{}/{}_{}_{}_{}2.xlsx'.format(state,miRNA, genotype, state, period), index_col =0)
                    df = pd.DataFrame(index = df1.index, columns = df1.columns)
                    mouse = df1.index.to_list()
                    bouts = df1.columns.to_list()
                    for mouse in mouse :
                        values1 = df1.loc[mouse, bouts].values
                        values2 = df2.loc[mouse, bouts].values
                        weighted_mean = (12/24)*values1 + (12/24)*values2
                        df.loc[mouse, bouts] = weighted_mean
                    df.to_excel(dirname+'/{}/{}_{}_{}_BASELINE_{}.xlsx'.format(state, miRNA, genotype,state, period))
                    period = 'dark'
                    df1 = pd.read_excel(dirname+'/{}/{}_{}_{}_{}1.xlsx'.format(state,miRNA, genotype, state, period), index_col =0)
                    df2 = pd.read_excel(dirname+'/{}/{}_{}_{}_{}2.xlsx'.format(state,miRNA, genotype, state, period), index_col =0)
                    df = pd.DataFrame(index = df1.index, columns = df1.columns)
                    mouse = df1.index.to_list()
                    bouts = df1.columns.to_list()
                    for mouse in mouse :
                        values1 = df1.loc[mouse, bouts].values
                        values2 = df2.loc[mouse, bouts].values
                        weighted_mean = (12/24)*values1 + (12/24)*values2
                        df.loc[mouse, bouts] = weighted_mean
                    df.to_excel(dirname+'/{}/{}_{}_{}_BASELINE_{}.xlsx'.format(state, miRNA, genotype,state, period))

def average_recovery_fragmentation(exclude_mice = False):
    date_ref = pd.read_excel(work_dir + 'datetime_reference_miRNA.xlsx', index_col = 0)
    encoding_to_state = {'a':'cata', 'w':'wake', 'n':'nrem', 'r':'rem', 't':'tdw'}
    for folder in ['sleep_fragmentation','sleep_percent_fragmentation']:
        for miRNA in ['128', '137', '665']:
            groups = {'control' : get_mice(miRNA, 'control'), 'test' : get_mice(miRNA, 'test')}
            for genotype in groups:
                if exclude_mice:
                    dirname = excel_dir + '/RT_PCR_excluded_mice/{}/{}/{}/'.format(folder,miRNA, genotype)
                else :
                    dirname = excel_dir + '/{}/{}/{}/'.format(folder,miRNA, genotype)
                encoding = ['w', 'n', 'r','t']
                for e in encoding :
                    state = encoding_to_state[e]
                    period = 'light'
                    df1 = pd.read_excel(dirname+'/{}/{}_{}_{}_{}3.xlsx'.format(state,miRNA, genotype, state, period), index_col =0)
                    df2 = pd.read_excel(dirname+'/{}/{}_{}_{}_{}4.xlsx'.format(state,miRNA, genotype, state, period), index_col =0)
                    df = pd.DataFrame(index = df1.index, columns = df1.columns)
                    mouse = df1.index.to_list()
                    bouts = df1.columns.to_list()
                    for mouse in mouse :
                        values1 = df1.loc[mouse, bouts].values
                        values2 = df2.loc[mouse, bouts].values
                        weighted_mean = (6/18)*values1 + (12/18)*values2
                        df.loc[mouse, bouts] = weighted_mean
                    df.to_excel(dirname+'/{}/{}_{}_{}_RECOVERY_{}.xlsx'.format(state,miRNA, genotype, state, period))
                    period = 'dark'
                    df1 = pd.read_excel(dirname+'/{}/{}_{}_{}_{}3.xlsx'.format(state,miRNA, genotype, state, period), index_col =0)
                    df2 = pd.read_excel(dirname+'/{}/{}_{}_{}_{}4.xlsx'.format(state,miRNA, genotype, state, period), index_col =0)
                    df = pd.DataFrame(index = df1.index, columns = df1.columns)
                    mouse = df1.index.to_list()
                    bouts = df1.columns.to_list()
                    for mouse in mouse :
                        if date_ref.at['MTA-{}'.format(mouse), 'short_end'] == 'yes':
                            weighted_mean= df1.loc[mouse, bouts].values
                        else :
                            values1 = df1.loc[mouse, bouts].values
                            values2 = df2.loc[mouse, bouts].values
                            weighted_mean = (12/24)*values1 + (12/24)*values2
                        df.loc[mouse, bouts] = weighted_mean
                    df.to_excel(dirname+'/{}/{}_{}_{}_RECOVERY_{}.xlsx'.format(state,miRNA, genotype, state, period))



def average_baseline_directory(directory = 'state_duration', exclude_mice = False):
    encoding_to_state = {'w':'wake', 'n':'nrem', 'r':'rem', 't':'tdw'}
    for miRNA in ['128', '137', '665']:
        groups = {'control' : get_mice(miRNA, 'control'), 'test' : get_mice(miRNA, 'test')}
        for genotype in groups:
            if exclude_mice:
                dirname = excel_dir + '/RT_PCR_excluded_mice/{}/{}/{}/'.format(directory,miRNA, genotype)
            else :
                dirname = excel_dir + '/{}/{}/{}/'.format(directory,miRNA, genotype)
            encoding = ['w', 'n', 'r', 't']
            for e in encoding :
                state = encoding_to_state[e]
                period = 'light'
                df1 = pd.read_excel(dirname+'/{}/{}_{}_{}_{}1.xlsx'.format(state, miRNA, genotype, state, period), index_col =0)
                df2 = pd.read_excel(dirname+'/{}/{}_{}_{}_{}2.xlsx'.format(state, miRNA, genotype, state, period), index_col =0)
                df = pd.DataFrame(index = df1.index, columns = df1.columns)
                mouse = df1.index.to_list()
                for mouse in mouse :
                    # values1 = df1.loc[mouse, 'mean_duration']
                    # values2 = df2.loc[mouse, 'mean_duration']
                    # weighted_mean = (12/24)*values1 + (12/24)*values2
                    # df.loc[mouse, 'mean_duration'] = weighted_mean
                    for col in df.columns:
                          values1 = df1.loc[mouse,col]
                          values2 = df2.loc[mouse,col]
                          weighted_mean = (12/24)*values1 + (12/24)*values2
                          df.loc[mouse, col] = weighted_mean
                df.to_excel(dirname+'/{}/{}_{}_{}_BASELINE_{}.xlsx'.format(state, miRNA, genotype, state, period))
                period = 'dark'
                df1 = pd.read_excel(dirname+'/{}/{}_{}_{}_{}1.xlsx'.format(state, miRNA, genotype, state, period), index_col =0)
                df2 = pd.read_excel(dirname+'/{}/{}_{}_{}_{}2.xlsx'.format(state, miRNA, genotype, state, period), index_col =0)
                df = pd.DataFrame(index = df1.index, columns = df1.columns)
                mouse = df1.index.to_list()
                for mouse in mouse :
                    # values1 = df1.loc[mouse, 'mean_duration']
                    # values2 = df2.loc[mouse, 'mean_duration']
                    # weighted_mean = (12/24)*values1 + (12/24)*values2
                    # df.loc[mouse, 'mean_duration'] = weighted_mean
                    for col in df.columns:
                          values1 = df1.loc[mouse,col]
                          values2 = df2.loc[mouse,col]
                          weighted_mean = (12/24)*values1 + (12/24)*values2
                          df.loc[mouse, col] = weighted_mean
                df.to_excel(dirname+'/{}/{}_{}_{}_BASELINE_{}.xlsx'.format(state, miRNA, genotype, state, period))

def average_recovery_directory(directory = 'state_duration', exclude_mice = False):
    date_ref = pd.read_excel(work_dir + 'datetime_reference_miRNA.xlsx', index_col = 0)
    encoding_to_state = {'a':'cata', 'w':'wake', 'n':'nrem', 'r':'rem', 't':'tdw'}
    for miRNA in ['128', '137', '665']:
        groups = {'control' : get_mice(miRNA, 'control'), 'test' : get_mice(miRNA, 'test')}
        for genotype in groups:
            if exclude_mice:
                dirname = excel_dir + '/RT_PCR_excluded_mice/{}/{}/{}/'.format(directory,miRNA, genotype)
            else:
                dirname = excel_dir + '/{}/{}/{}/'.format(directory,miRNA, genotype)
            encoding = ['w', 'n', 'r', 't']
            for e in encoding :
                state = encoding_to_state[e]
                period = 'light'
                df1 = pd.read_excel(dirname+'/{}/{}_{}_{}_{}3.xlsx'.format(state,miRNA, genotype, state, period), index_col =0)
                df2 = pd.read_excel(dirname+'/{}/{}_{}_{}_{}4.xlsx'.format(state,miRNA, genotype, state, period), index_col =0)
                df = pd.DataFrame(index = df1.index, columns = df1.columns)
                mouse = df1.index.to_list()
                for mouse in mouse :
                    # values1 = df1.loc[mouse, 'mean_duration']
                    # values2 = df2.loc[mouse, 'mean_duration']
                    # weighted_mean = (6/18)*values1 + (12/18)*values2
                    # df.loc[mouse, 'mean_duration'] = weighted_mean
                    for col in df.columns:
                          values1 = df1.loc[mouse,col]
                          values2 = df2.loc[mouse,col]
                          weighted_mean = (6/18)*values1 + (12/18)*values2
                          df.loc[mouse, col] = weighted_mean
                df.to_excel(dirname+'/{}/{}_{}_{}_RECOVERY_{}.xlsx'.format(state,miRNA, genotype, state, period))
                period = 'dark'
                df1 = pd.read_excel(dirname+'/{}/{}_{}_{}_{}3.xlsx'.format(state,miRNA, genotype, state, period), index_col =0)
                df2 = pd.read_excel(dirname+'/{}/{}_{}_{}_{}4.xlsx'.format(state,miRNA, genotype, state, period), index_col =0)
                df = pd.DataFrame(index = df1.index, columns = df1.columns)
                mouse = df1.index.to_list()
                for mouse in mouse :
                    if date_ref.at['MTA-{}'.format(mouse), 'short_end'] == 'yes':
                        # weighted_mean= df1.loc[mouse, 'mean_duration']
                        for col in df.columns:
                              weighted_mean = df1.loc[mouse,col]
                              df.loc[mouse, col] = weighted_mean
                    else :
                        # values1 = df1.loc[mouse, 'mean_duration']
                        # values2 = df2.loc[mouse, 'mean_duration']
                        # weighted_mean = (12/24)*values1 + (12/24)*values2
                        for col in df.columns:
                              values1 = df1.loc[mouse,col]
                              values2 = df2.loc[mouse,col]
                              weighted_mean = (12/24)*values1 + (12/24)*values2
                              df.loc[mouse, col] = weighted_mean
                    # df.loc[mouse, 'mean_duration'] = weighted_mean
                df.to_excel(dirname+'/{}/{}_{}_{}_RECOVERY_{}.xlsx'.format(state,miRNA, genotype, state, period))


def average_baseline_IRI(exclude_mice = False):
    directory = 'IRI'
    encoding_to_state = {'w':'wake', 'n':'nrem', 'r':'rem'}
    for miRNA in ['128', '137', '665']:
        groups = {'control' : get_mice(miRNA, 'control'), 'test' : get_mice(miRNA, 'test')}
        for genotype in groups:
            if exclude_mice :
                dirname = excel_dir + '/RT_PCR_excluded_mice/{}/{}/{}/'.format(directory,miRNA, genotype)
            else:
                dirname = excel_dir + '/{}/{}/{}/'.format(directory,miRNA, genotype)
            encoding = ['w', 'n', 'r']
            period = 'light'
            df1 = pd.read_excel(dirname+'/{}_{}_{}1.xlsx'.format(miRNA, genotype, period), index_col =0)
            df2 = pd.read_excel(dirname+'/{}_{}_{}2.xlsx'.format(miRNA, genotype, period), index_col =0)
            df = pd.DataFrame(index = df1.index, columns = df1.columns)
            mouse = df1.index.to_list()
            for mouse in mouse :
                # values1 = df1.loc[mouse, 'mean_duration']
                # values2 = df2.loc[mouse, 'mean_duration']
                # weighted_mean = (12/24)*values1 + (12/24)*values2
                # df.loc[mouse, 'mean_duration'] = weighted_mean
                for col in df.columns:
                      values1 = df1.loc[mouse,col]
                      values2 = df2.loc[mouse,col]
                      weighted_mean = (12/24)*values1 + (12/24)*values2
                      df.loc[mouse, col] = weighted_mean
            df.to_excel(dirname+'/{}_{}_BASELINE_{}.xlsx'.format(miRNA, genotype, period))
            period = 'dark'
            df1 = pd.read_excel(dirname+'/{}_{}_{}1.xlsx'.format(miRNA, genotype, period), index_col =0)
            df2 = pd.read_excel(dirname+'/{}_{}_{}2.xlsx'.format(miRNA, genotype, period), index_col =0)
            df = pd.DataFrame(index = df1.index, columns = df1.columns)
            mouse = df1.index.to_list()
            for mouse in mouse :
                # values1 = df1.loc[mouse, 'mean_duration']
                # values2 = df2.loc[mouse, 'mean_duration']
                # weighted_mean = (12/24)*values1 + (12/24)*values2
                # df.loc[mouse, 'mean_duration'] = weighted_mean
                for col in df.columns:
                      values1 = df1.loc[mouse,col]
                      values2 = df2.loc[mouse,col]
                      weighted_mean = (12/24)*values1 + (12/24)*values2
                      df.loc[mouse, col] = weighted_mean
            df.to_excel(dirname+'/{}_{}_BASELINE_{}.xlsx'.format(miRNA, genotype, period))


def average_recovery_IRI(exclude_mice = False):
    directory = 'IRI'
    date_ref = pd.read_excel(work_dir + 'datetime_reference_miRNA.xlsx', index_col = 0)
    encoding_to_state = {'a':'cata', 'w':'wake', 'n':'nrem', 'r':'rem', 't':'tdw'}
    for miRNA in ['128', '137', '665']:
        groups = {'control' : get_mice(miRNA, 'control'), 'test' : get_mice(miRNA, 'test')}
        for genotype in groups:
            if exclude_mice :
                dirname = excel_dir + '/RT_PCR_excluded_mice/{}/{}/{}/'.format(directory,miRNA, genotype)
            else :
                dirname = excel_dir + '/{}/{}/{}/'.format(directory,miRNA, genotype)
            period = 'light'
            df1 = pd.read_excel(dirname+'/{}_{}_{}3.xlsx'.format(miRNA, genotype, period), index_col =0)
            df2 = pd.read_excel(dirname+'/{}_{}_{}4.xlsx'.format(miRNA, genotype, period), index_col =0)
            df = pd.DataFrame(index = df1.index, columns = df1.columns)
            mouse = df1.index.to_list()
            for mouse in mouse :
                # values1 = df1.loc[mouse, 'mean_duration']
                # values2 = df2.loc[mouse, 'mean_duration']
                # weighted_mean = (6/18)*values1 + (12/18)*values2
                # df.loc[mouse, 'mean_duration'] = weighted_mean
                for col in df.columns:
                      values1 = df1.loc[mouse,col]
                      values2 = df2.loc[mouse,col]
                      weighted_mean = (6/18)*values1 + (12/18)*values2
                      df.loc[mouse, col] = weighted_mean
            df.to_excel(dirname+'/{}_{}_RECOVERY_{}.xlsx'.format(miRNA, genotype, period))
            period = 'dark'
            df1 = pd.read_excel(dirname+'/{}_{}_{}3.xlsx'.format(miRNA, genotype, period), index_col =0)
            df2 = pd.read_excel(dirname+'/{}_{}_{}4.xlsx'.format(miRNA, genotype, period), index_col =0)
            df = pd.DataFrame(index = df1.index, columns = df1.columns)
            mouse = df1.index.to_list()
            for mouse in mouse :
                if date_ref.at['MTA-{}'.format(mouse), 'short_end'] == 'yes':
                    # weighted_mean= df1.loc[mouse, 'mean_duration']
                    for col in df.columns:
                          weighted_mean = df1.loc[mouse,col]
                          df.loc[mouse, col] = weighted_mean
                else :
                    # values1 = df1.loc[mouse, 'mean_duration']
                    # values2 = df2.loc[mouse, 'mean_duration']
                    # weighted_mean = (12/24)*values1 + (12/24)*values2
                    for col in df.columns:
                          values1 = df1.loc[mouse,col]
                          values2 = df2.loc[mouse,col]
                          weighted_mean = (12/24)*values1 + (12/24)*values2
                          df.loc[mouse, col] = weighted_mean
                # df.loc[mouse, 'mean_duration'] = weighted_mean
            df.to_excel(dirname+'/{}_{}_RECOVERY_{}.xlsx'.format(miRNA, genotype, period))


def average_TDW_ratio(exclude_mice = False):
    directory = 'tdw_w_ratio'
    encoding_to_state = {'w':'wake', 'n':'nrem', 'r':'rem', 't':'tdw'}
    for miRNA in ['128', '137', '665']:
        groups = {'control' : get_mice(miRNA, 'control'), 'test' : get_mice(miRNA, 'test')}
        for genotype in groups:
            if exclude_mice :
                dirname = excel_dir + '/RT_PCR_excluded_mice/{}/{}/{}/'.format(directory,miRNA, genotype)
            else:
                dirname = excel_dir + '/{}/{}/{}/'.format(directory,miRNA, genotype)

            df = pd.read_excel(dirname+'/tdw_w_ratio_by_mouse_by_12hours.xlsx', index_col =0)
            columns = ['BASELINE_light', 'BASELINE_dark', 'RECOVERY_light', 'RECOVERY_dark']
            df_av = pd.DataFrame(index = df.index, columns = columns)
            mouse = df.index.to_list()
            pairs = [[0,2], [1,3],[4,6],[5,7]]
            for mouse in mouse :
                # values1 = df1.loc[mouse, 'mean_duration']
                # values2 = df2.loc[mouse, 'mean_duration']
                # weighted_mean = (12/24)*values1 + (12/24)*values2
                # df.loc[mouse, 'mean_duration'] = weighted_mean
                for c,col in enumerate(columns):
                    pair = pairs[c]
                    values1 = df.at[mouse,pair[0]]
                    values2 = df.at[mouse,pair[1]]
                    weighted_mean = (12/24)*values1 + (12/24)*values2
                    df_av.loc[mouse, col] = weighted_mean
            df_av.to_excel(dirname+'/tdw_12hours_AVERAGE.xlsx')


def fusion_baseline_recovery_REM_latency():
    for miRNA in ['128', '137', '665']:
        for m, merging in enumerate(['BASELINE', 'RECOVERY']):
            day1, day2 = 2*m+1, 2*(m+1)
            for genotype in ['control', 'test']:
                dirname = excel_dir + '/REM_latency/{}/{}/'.format(miRNA, genotype)
                for period in ['dark', 'light']:
                    df1 = pd.read_excel(dirname+'/REM_latency_{}_{}_{}{}.xlsx'.format(miRNA, genotype, period,day1), index_col =0)
                    df2 = pd.read_excel(dirname+'/REM_latency_{}_{}_{}{}.xlsx'.format(miRNA, genotype, period,day2), index_col =0)
                    df = pd.DataFrame(index = df1.index, columns = np.arange(4))
                    mice = np.unique(df1.index.to_list()+df2.index.to_list())
                    for mouse in mice:
                        list1 = df1.loc[mouse,:]
                        list2 = df2.loc[mouse,:]
                        list1 = list1[~np.isnan(list1)].tolist()
                        list2 = list2[~np.isnan(list2)].tolist()
                        latencies = list1 + list2
                        if len(df.columns.to_list())<len(latencies):
                            df = df.reindex(columns = np.arange(len(latencies)))
                        df.loc[mouse, np.arange(len(latencies))] = np.array(latencies)
                    filename = dirname + 'REM_latency_{}_{}_{}_{}.xlsx'.format(miRNA,genotype, merging, period)
                    df.to_excel(filename)

def average_all(exclude_mice=False):
    average_baseline_spectrum(spectrum_method = 'welch', exclude_mice = exclude_mice, auc = False, agg_power='mean')
    average_recovery_spectrum(spectrum_method = 'welch', exclude_mice = exclude_mice, auc = False, agg_power='mean')
    average_baseline_fragmentation(exclude_mice = exclude_mice)
    average_recovery_fragmentation(exclude_mice = exclude_mice)
    average_baseline_directory('state_time',exclude_mice = exclude_mice)
    average_recovery_directory('state_time',exclude_mice = exclude_mice)
    average_baseline_directory('state_duration',exclude_mice = exclude_mice)
    average_recovery_directory('state_duration',exclude_mice = exclude_mice)
    average_baseline_directory('period_amount',exclude_mice = exclude_mice)
    average_recovery_directory('period_amount',exclude_mice = exclude_mice)
    average_TDW_ratio(exclude_mice = exclude_mice)
    average_TDW_ratio(exclude_mice = exclude_mice)
    average_baseline_IRI(exclude_mice)
    average_recovery_IRI(exclude_mice)

if __name__ == '__main__':
    # average_baseline_spectrum('welch')
    # average_recovery_spectrum('welch')
    # average_baseline_fragmentation()
    # average_recovery_fragmentation()
    # average_baseline_fragmentation_mean()
    # average_recovery_fragmentation_mean()
    # fusion_baseline_recovery_REM_latency()
    # average_TDW_ratio(True)

    average_all(exclude_mice=False)
    average_all(exclude_mice=True)
    # average_baseline_directory('state_time')
    # average_recovery_directory('state_time')
    # average_baseline_directory('state_duration')
    # average_recovery_directory('state_duration')
    # average_baseline_directory('period_amount')
    # average_recovery_directory('period_amount')
    # average_baseline_IRI()
    # average_recovery_IRI()
