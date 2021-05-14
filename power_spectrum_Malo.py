from configuration import *
from select_mice_cata_Malo import get_mice, get_mouse_info
from sklearn.metrics import auc

local_path = os.path.dirname(os.path.realpath(__file__))
print(local_path)


# def get_clear_spectral_score_one_mouse_one_condition(mouse, condition, state):
#     debug = False
#     ds = xr.open_dataset(precompute_dir + '/spectrums/spectrum_scoring_{}.nc'.format(mouse))
#     times = ds.coords['times_somno'].values/3600
#     # score = ds['score'].values
#     dirname = precompute_dir + '/tdw_score/'
#     ds_tdw = xr.open_dataset(dirname + 'tdw_score_{}.nc'.format(mouse))
#     score = ds_tdw['new_score'].values
#
#     if condition == 'bl0':
#         t2 = 4
#         mask = (times<t2)
#     if condition == 'bl1':
#         t1, t2 = 4, 28
#         mask = (times>t1)&(times<t2)
#     elif condition == 'bl2':
#         t1, t2 = 24, 52
#         mask = (times>=t1) & (times<t2)
#     elif condition == 'sd':
#         t1, t2 = 52+6, 76
#         mask = (times>=t1) & (times<t2)
#     elif condition == 'r1':
#         t1 = 76
#         mask = (times>=t1)
#     else :
#         print('Condition does not exist !')
#
#     # mask = (times>=t1) & (times<t2)
#     score = score[mask]
#
#     if state == 'w':
#             score_no_artifact = (score == 'w' )|(score == 't')
#     else :
#         score_no_artifact = score == state
#
#     score_no_transition = score_no_artifact.copy()
#     count = 0
#     event_length = []
#     for num, value in enumerate(score_no_artifact) :
#         if value:
#             count +=1
#         elif value == False:
#             event_length.append(count)
#             if count == 1 :
#                 score_no_transition[num-1] = False
#             elif count == 2 :
#                 score_no_transition[num-1] = False
#                 score_no_transition[num-2] = False
#                 count = 0
#             elif count > 2 :
#                 score_no_transition[num-1] = False
#                 score_no_transition[num-count] = False
#             count = 0
#     if count == 1 :
#         score_no_transition[num-1] = False
#     elif count == 2 :
#         score_no_transition[num-1] = False
#         score_no_transition[num-2] = False
#         count = 0
#     elif count > 2 :
#         score_no_transition[num-1] = False
#         score_no_transition[num-count] = False
#
#
#     if debug :
#         fig, ax = plt.subplots()
#         s = ds.coords['times_somno'].values.size
#         for i in range(5):
#             ax.axvline(i*s/(3600))
#         h = ds.coords['times_somno'].values/3600
#         ax.plot(h, np.ones(s))
#         ax.plot(h[mask], np.ones(s)[mask])
#         ax.plot(h[mask], score == state)
#
#         fig, ax = plt.subplots()
#         ax.scatter(h[mask], score_no_artifact)
#         ax.plot(h[mask], score_no_transition, color = 'green')
#         print(count)
#         plt.show()
#
#     return score_no_transition, mask
#
#  def get_relative_spectrums_one_mouse_one_condition(mouse, condition, spectrum_method = 'somno'):
#     debug = 0
#     cleared_score_w, _= get_clear_spectral_score_one_mouse_one_condition(mouse, condition, 'w')
#     cleared_score_t, _= get_clear_spectral_score_one_mouse_one_condition(mouse, condition, 't')
#     cleared_score_r, _ = get_clear_spectral_score_one_mouse_one_condition(mouse, condition, 'r')
#     cleared_score_n, mask_condition = get_clear_spectral_score_one_mouse_one_condition(mouse, condition, 'n')
#
#     ds = xr.open_dataset(precompute_dir + '/spectrums/spectrum_scoring_{}.nc'.format(mouse))
#
#     spectrum = ds['{}_spectrum'.format(spectrum_method)].values
#     freqs = ds.coords['freqs_{}'.format(spectrum_method)].values
#     spectrum = spectrum[mask_condition, :]
#
#     total_epochs = sum(cleared_score_w) + sum(cleared_score_r) + sum(cleared_score_n)
#
#     proportion_w = sum(cleared_score_w)/total_epochs
#     proportion_r = sum(cleared_score_r)/total_epochs
#     proportion_n = sum(cleared_score_n)/total_epochs
#     proportion_t = sum(cleared_score_t)/total_epochs
#
#     absolute_spectrum_w = np.median(spectrum[cleared_score_w, :], axis = 0)
#     absolute_spectrum_r = np.median(spectrum[cleared_score_r, :], axis = 0)
#     absolute_spectrum_n = np.median(spectrum[cleared_score_n, :], axis = 0)
#     absolute_spectrum_t = np.median(spectrum[cleared_score_t, :], axis = 0)
#
#     #########CAUTION power spectrum is area under curv#########
#     auc_w = auc(freqs, absolute_spectrum_w)
#     auc_r = auc(freqs, absolute_spectrum_r)
#     auc_n = auc(freqs, absolute_spectrum_n)
#     auc_t = auc(freqs, absolute_spectrum_t)
#     ref_values = proportion_w*auc_w + proportion_r*auc_r + proportion_n*auc_n + proportion_t*auc_t
#
#
#     ####FALSE
#     # mean_w = np.mean(absolute_spectrum_w)
#     # mean_r = np.mean(absolute_spectrum_r)
#     # mean_n = np.mean(absolute_spectrum_n)
#     # ref_values = proportion_w*mean_w + proportion_r*mean_r + proportion_n*mean_n
#
#
#
#     #####ALi's
#     # sum_w = np.sum(absolute_spectrum_w)
#     # sum_r = np.sum(absolute_spectrum_r)
#     # sum_n = np.sum(absolute_spectrum_n)
#     # ref_values = (proportion_w*sum_w + proportion_r*sum_r + proportion_n*sum_n)/100
#     #
#     # relative_spectrum_w = 100*absolute_spectrum_w / ref_values
#     # relative_spectrum_r = 100*absolute_spectrum_r / ref_values
#     # relative_spectrum_n = 100*absolute_spectrum_n / ref_values
#
#     # print(ref_values)
#     # fig, ax = plt.subplots()
#     # global_spec = np.mean(spectrum[cleared_score_w + cleared_score_r + cleared_score_n], axis = 0)
#     # print(auc(freqs,global_spec))
#     # ax.plot(freqs,global_spec )
#     # ax.plot(freqs,global_spec/ref_values )
#     #
#     #
#     # fig, ax = plt.subplots()
#     # ax.plot(freqs, relative_spectrum_w)
#     # ax.plot(freqs, relative_spectrum_r)
#     # ax.plot(freqs, relative_spectrum_n)
#     # plt.show()
#
#     if debug == True:
#         fig, ax = plt.subplots()
#         ax.plot(absolute_spectrum_w, color = 'blue')
#         ax.plot(absolute_spectrum_r, color = 'blue')
#         ax.plot(absolute_spectrum_n, color = 'blue')
#
#         ax.plot(relative_spectrum_w, color = 'green')
#         ax.plot(relative_spectrum_r, color = 'green', alpha = .66)
#         ax.plot(relative_spectrum_n, color = 'green', alpha = .33)
#         plt.show()
#
#     return ref_values, {'w' : relative_spectrum_w, 'r':relative_spectrum_r, 'n' :relative_spectrum_n,'t' :relative_spectrum_t}
def compute_ref_values_baseline_Sha(miRNA, genotype, spectrum_method = 'welch',period = 'light', auc = False, agg_power = 'mean'):
    mice = get_mice(miRNA, genotype)
    df_ref = pd.DataFrame(index=mice)
    for mouse in mice:
        cleared_score_w, _= get_clear_spectral_score_one_mouse_one_condition_day_night(mouse, 'bl1', 'w', period = period)
        cleared_score_r, _ = get_clear_spectral_score_one_mouse_one_condition_day_night(mouse, 'bl1', 'r', period = period)
        cleared_score_t, _ = get_clear_spectral_score_one_mouse_one_condition_day_night(mouse, 'bl1', 't', period = period)
        cleared_score_n, mask_condition = get_clear_spectral_score_one_mouse_one_condition_day_night(mouse, 'bl1', 'n', period = period)

        ds = xr.open_dataset(precompute_dir + '/spectrums/spectrum_scoring_{}.nc'.format(mouse))

        spectrum = ds['{}_spectrum'.format(spectrum_method)].values
        freqs = ds.coords['freqs_{}'.format(spectrum_method)].values
        spectrum = spectrum[mask_condition, :]

        total_epochs = np.sum(cleared_score_w) + np.sum(cleared_score_r) + np.sum(cleared_score_n)+ np.sum(cleared_score_t)

        proportion_w = np.sum(cleared_score_w)
        proportion_r = np.sum(cleared_score_r)
        proportion_n = np.sum(cleared_score_n)
        proportion_t = np.sum(cleared_score_t)

        absolute_spectrum_w = spectrum[cleared_score_w, :]
        absolute_spectrum_r = spectrum[cleared_score_r, :]
        absolute_spectrum_n = spectrum[cleared_score_n, :]
        absolute_spectrum_t = spectrum[cleared_score_t, :]

        cleared_score_w, _= get_clear_spectral_score_one_mouse_one_condition_day_night(mouse, 'bl2', 'w', period = period)
        cleared_score_r, _ = get_clear_spectral_score_one_mouse_one_condition_day_night(mouse, 'bl2', 'r', period = period)
        cleared_score_t, _ = get_clear_spectral_score_one_mouse_one_condition_day_night(mouse, 'bl2', 't', period = period)
        cleared_score_n, mask_condition = get_clear_spectral_score_one_mouse_one_condition_day_night(mouse, 'bl2', 'n', period = period)

        ds = xr.open_dataset(precompute_dir + '/spectrums/spectrum_scoring_{}.nc'.format(mouse))

        spectrum = ds['{}_spectrum'.format(spectrum_method)].values
        freqs = ds.coords['freqs_{}'.format(spectrum_method)].values
        spectrum = spectrum[mask_condition, :]

        total_epochs = total_epochs + np.sum(cleared_score_w) + np.sum(cleared_score_r) + np.sum(cleared_score_n)+ np.sum(cleared_score_t)

        proportion_w = (proportion_w+ np.sum(cleared_score_w))/total_epochs
        proportion_r = (proportion_r+ np.sum(cleared_score_r))/total_epochs
        proportion_n = (proportion_n+ np.sum(cleared_score_n))/total_epochs
        proportion_t = (proportion_t+ np.sum(cleared_score_t))/total_epochs

        if agg_power == 'median':
            absolute_spectrum_w = np.nanmedian(np.concatenate([absolute_spectrum_w,spectrum[cleared_score_w, :]]), axis = 0)
            absolute_spectrum_r = np.nanmedian(np.concatenate([absolute_spectrum_r,spectrum[cleared_score_r, :]]), axis = 0)
            absolute_spectrum_n = np.nanmedian(np.concatenate([absolute_spectrum_n,spectrum[cleared_score_n, :]]), axis = 0)
            absolute_spectrum_t = np.nanmedian(np.concatenate([absolute_spectrum_t,spectrum[cleared_score_t, :]]), axis = 0)

        elif agg_power == 'mean':
            absolute_spectrum_w = np.nanmean(np.concatenate([absolute_spectrum_w,spectrum[cleared_score_w, :]]), axis = 0)
            absolute_spectrum_r = np.nanmean(np.concatenate([absolute_spectrum_r,spectrum[cleared_score_r, :]]), axis = 0)
            absolute_spectrum_n = np.nanmean(np.concatenate([absolute_spectrum_n,spectrum[cleared_score_n, :]]), axis = 0)
            absolute_spectrum_t = np.nanmean(np.concatenate([absolute_spectrum_t,spectrum[cleared_score_t, :]]), axis = 0)

        df = (freqs[11]-freqs[1])/10


        if np.sum(cleared_score_w) !=0:
            if auc == True:
                sum_w = np.sum(absolute_spectrum_w)*df
            elif auc == False :
                sum_w = np.sum(absolute_spectrum_w)
        else:
            sum_w =0

        if np.sum(cleared_score_t) !=0:
            if auc == True:
                sum_t = np.sum(absolute_spectrum_t)*df
            elif auc == False :
                sum_t = np.sum(absolute_spectrum_t)
        else:
            sum_t =0

        if np.sum(cleared_score_r) !=0:
            if auc == True:
                sum_r = np.sum(absolute_spectrum_r)*df
            elif auc == False :
                sum_r = np.sum(absolute_spectrum_r)
        else:
            sum_r =0

        if np.sum(cleared_score_n) !=0:
            if auc == True:
                sum_n = np.sum(absolute_spectrum_n)*df
            elif auc == False :
                sum_n = np.sum(absolute_spectrum_n)
        else:
            sum_n = 0


        #########CAUTION power spectrum is area under curv#########
        # auc_w = auc(freqs, absolute_spectrum_w)
        # auc_r = auc(freqs, absolute_spectrum_r)
        # auc_n = auc(freqs, absolute_spectrum_n)
        # ref_values = proportion_w*auc_w + proportion_r*auc_r + proportion_n*auc_n


        #####ALi's
        # sum_w = np.sum(absolute_spectrum_w)*df
        # sum_r = np.sum(absolute_spectrum_r)*df
        # sum_n = np.sum(absolute_spectrum_n)*df
        ref_values = (proportion_w*sum_w + proportion_r*sum_r + proportion_n*sum_n +  proportion_t*sum_t)
        # print(ref_values)
        #
        # ref_values = (proportion_w*sum_w + proportion_r*sum_r + proportion_n*sum_n)
        # ref, _ = get_relative_spectrums_one_mouse_one_condition_day_night(mouse, 'bl1', spectrum_method, period)
        # print(ref)
        # ref, _ = get_relative_spectrums_one_mouse_one_condition_day_night(mouse, 'bl2', spectrum_method, period)
        # print(ref)
        # ref_values = (proportion_w*sum_w + proportion_r*sum_r + proportion_n*sum_n)
        # ref, _ = get_relative_spectrums_one_mouse_one_condition_day_night(mouse, 'sd', spectrum_method, period)
        # print(ref)
        # ref, _ = get_relative_spectrums_one_mouse_one_condition_day_night(mouse, 'r1', spectrum_method, period)
        # print(ref)
        df_ref.at[mouse,'ref'] = ref_values
    dirname = excel_dir + '/power_spectrum_baseline_auc_{}_agg_posser_{}/'.format(auc, agg_power)
    if not os.path.exists(dirname):
        os.makedirs(dirname)
    df_ref.to_excel(dirname + 'ref_{}_{}_{}_{}.xlsx'.format(miRNA, genotype,spectrum_method,period))



def get_ref_values_Sha(mouse,spectrum_method, period = 'light', auc = False, agg_power = 'mean'):
    miRNA, genotype = get_mouse_info(mouse)
    dirname = excel_dir + '/power_spectrum_baseline_auc_{}_agg_posser_{}/'.format(auc, agg_power)
    filname = dirname + 'ref_{}_{}_{}_{}.xlsx'.format(miRNA, genotype,spectrum_method,period)
    if not os.path.exists(filname):
        compute_ref_values_baseline_Sha(miRNA, genotype, spectrum_method, period, auc, agg_power)
    df = pd.read_excel(filname, index_col =0)
    ref_values = df.at[mouse, 'ref']
    return ref_values
    # exit()

def get_relative_spectrums_one_mouse_one_condition_day_night_Sha(mouse, condition, spectrum_method = 'somno', period = 'light', auc = False, agg_power = 'mean'):

    cleared_score_w, _= get_clear_spectral_score_one_mouse_one_condition_day_night(mouse, condition, 'w', period = period)
    cleared_score_r, _ = get_clear_spectral_score_one_mouse_one_condition_day_night(mouse, condition, 'r', period = period)
    cleared_score_t, _ = get_clear_spectral_score_one_mouse_one_condition_day_night(mouse, condition, 't', period = period)
    cleared_score_n, mask_condition = get_clear_spectral_score_one_mouse_one_condition_day_night(mouse, condition, 'n', period = period)

    ds = xr.open_dataset(precompute_dir + '/spectrums/spectrum_scoring_{}.nc'.format(mouse))

    spectrum = ds['{}_spectrum'.format(spectrum_method)].values
    freqs = ds.coords['freqs_{}'.format(spectrum_method)].values
    spectrum = spectrum[mask_condition, :]

    # total_epochs = np.sum(cleared_score_w) + np.sum(cleared_score_r) + np.sum(cleared_score_n)+ np.sum(cleared_score_t)
    #
    # proportion_w = np.sum(cleared_score_w)/total_epochs
    # proportion_r = np.sum(cleared_score_r)/total_epochs
    # proportion_n = np.sum(cleared_score_n)/total_epochs
    # proportion_t = np.sum(cleared_score_t)/total_epochs


    if agg_power == 'median':
        absolute_spectrum_w = np.nanmedian(spectrum[cleared_score_w, :], axis = 0)
        absolute_spectrum_r = np.nanmedian(spectrum[cleared_score_r, :], axis = 0)
        absolute_spectrum_n = np.nanmedian(spectrum[cleared_score_n, :], axis = 0)
        absolute_spectrum_t = np.nanmedian(spectrum[cleared_score_t, :], axis = 0)
    elif agg_power=='mean':
        absolute_spectrum_w = np.nanmean(spectrum[cleared_score_w, :], axis = 0)
        absolute_spectrum_r = np.nanmean(spectrum[cleared_score_r, :], axis = 0)
        absolute_spectrum_n = np.nanmean(spectrum[cleared_score_n, :], axis = 0)
        absolute_spectrum_t = np.nanmean(spectrum[cleared_score_t, :], axis = 0)

    ref_values = get_ref_values_Sha(mouse,spectrum_method, period = 'light', auc = auc, agg_power = agg_power)

    relative_spectrum_w = 100*absolute_spectrum_w / ref_values
    relative_spectrum_r = 100*absolute_spectrum_r / ref_values
    relative_spectrum_n = 100*absolute_spectrum_n / ref_values
    relative_spectrum_t = 100*absolute_spectrum_t / ref_values

    return ref_values, {'w' : relative_spectrum_w, 'r':relative_spectrum_r, 'n' :relative_spectrum_n, 't' :relative_spectrum_t}



def plot_spectrum_compare_day_night_Sha(spectrum_method = 'somno', period = 'light', exclude_mice = False, auc = False, agg_power = 'mean'):
    date_ref = pd.read_excel(work_dir + 'datetime_reference_miRNA.xlsx', index_col = 0)
    if period == 'light':
        sessions = ['bl1', 'bl2', 'sd', 'r1' ]
    elif period == 'dark':
        sessions = ['bl0', 'bl1', 'bl2', 'sd', 'r1' ]

    for miRNA in ['128', '137', '665']:
        mice_control = get_mice(miRNA, 'control', exclude_mice)
        mice_test = get_mice(miRNA, 'test', exclude_mice)

        groups = {'control' :mice_control , 'test' : mice_test}
        ds = xr.open_dataset(precompute_dir + '/spectrums/spectrum_scoring_B2884.nc')   #B2576
        freqs = ds.coords['freqs_{}'.format(spectrum_method)].values
        freqs = np.round(freqs, 2).tolist()

        indices = []
        for state in ['w', 'n', 'r', 't']:
            for session in sessions:
                indices+=([i + '_'+session+ '_'+state for i in groups['test']+ groups['control']])
        df_spectrum = pd.DataFrame(index = indices, columns = ['group', 'session', 'state', 'ref_values'] + freqs)



        for group in groups:
            for mouse in groups[group]:
                ds = xr.open_dataset(precompute_dir + '/spectrums/spectrum_scoring_{}.nc'.format(mouse))
                freqs = ds.coords['freqs_{}'.format(spectrum_method)].values
                freqs = np.round(freqs, 2)

                ZT = ds.coords['time_epochs'].values
                print('collecting relative spectrums ', mouse)
                states = ['w', 'n', 'r', 't']
                for session in sessions:
                    if period == 'dark' and date_ref.at['MTA-{}'.format(mouse), 'short_end'] == 'yes' and session == 'r1':
                        for state in states :
                            df_spectrum.at['{}_{}_{}'.format(mouse,session, state), 'state'] = state
                            df_spectrum.at['{}_{}_{}'.format(mouse,session, state), 'session'] = session
                            df_spectrum.at['{}_{}_{}'.format(mouse,session, state), 'group'] = group
                            df_spectrum.at['{}_{}_{}'.format(mouse,session, state), 'ref_values'] = np.nan
                            empty = np.ones(freqs.size)
                            empty[:] = np.nan
                            df_spectrum.at['{}_{}_{}'.format(mouse,session, state), 4:] = empty
                    else :
                        print(session, mouse, period)
                        ref_values, relative_spectrums = get_relative_spectrums_one_mouse_one_condition_day_night_Sha(mouse, session, spectrum_method, period, auc, agg_power )
                        # print(relative_spectrums[state].size)
                        for state in states :
                            df_spectrum.at['{}_{}_{}'.format(mouse,session, state), 'state'] = state
                            df_spectrum.at['{}_{}_{}'.format(mouse,session, state), 'session'] = session
                            df_spectrum.at['{}_{}_{}'.format(mouse,session, state), 'group'] = group
                            df_spectrum.at['{}_{}_{}'.format(mouse,session, state), 'ref_values'] = ref_values
                            df_spectrum.at['{}_{}_{}'.format(mouse,session, state), 4:] = relative_spectrums[state]
        # print(df_spectrum)
            # df_spectrum = df_spectrum.dropna()
        if period == 'light':
            Nrows =4
        elif period == 'dark':
            Nrows=5
        fig,ax = plt.subplots(nrows = Nrows, ncols= 4, sharex = True, sharey = True)
        fig_diff, ax_diff = plt.subplots(nrows = Nrows, ncols= 4, sharex = True, sharey = True)

        fig.suptitle(miRNA +' -- Spectrum per state -- ' + period + ' -exclude mice ='+str(exclude_mice))
        fig_diff.suptitle(miRNA + ' -- Diff DCR - ctrl')

        for i in range(Nrows):
            ax[i,0].set_ylabel(sessions[i])
            ax_diff[i,0].set_ylabel(sessions[i])

        for row, session in enumerate(sessions):
            for col, state in enumerate(states):
                # print(df_spectrum[(df_spectrum.session == session)& (df_spectrum.state == state)])
                data_dcr = df_spectrum[(df_spectrum.session == session) & (df_spectrum.group == 'test') & (df_spectrum.state == state)]
                data_ctrl = df_spectrum[(df_spectrum.session == session) & (df_spectrum.group == 'control') & (df_spectrum.state == state)]
                # print(data_dcr)
                # print(data_dcr.columns)
                # print(freqs)
                data_dcr = data_dcr.dropna()
                data_ctrl = data_ctrl.dropna()
                data_dcr = (data_dcr[freqs]).values
                data_ctrl = (data_ctrl[freqs]).values
                # print()
                # for i, maxi in enumerate(np.max(data_dcr, axis = 1)):('')
                #     data_dcr[i] = data_dcr[i]/maxi
                # for i, maxi in enumerate(np.max(data_ctrl, axis = 1)):
                #     data_ctrl[i] = data_ctrl[i]/maxi
                diff = np.mean(data_dcr, axis=0) - np.mean(data_ctrl, axis=0)
                # ax[row, col].plot(data.T, color = 'black', alpha = .1)
                # ax[row, col].plot(data.mean(axis = 0), color = 'red')
                ax_diff[row, col].plot(freqs, diff )
        for group in groups:
            if exclude_mice:
                dirname = excel_dir + '/RT_PCR_excluded_mice/power_spectrum_baseline_auc_{}_agg_power_{}/{}/{}/'.format(auc, agg_power,miRNA,group)
            else :
                dirname = excel_dir + '/power_spectrum_baseline_auc_{}_agg_power_{}/{}/{}/'.format(auc, agg_power,miRNA,group)

            if not os.path.exists(dirname):
                os.makedirs(dirname)

            for session in sessions:
                for state in states:
                    name = '{}_spectrum_{}_{}_{}_{}_exclude_mice_{}.xlsx'.format(miRNA,spectrum_method,state, session,period, exclude_mice)
                    df = df_spectrum[(df_spectrum.group == group) & (df_spectrum.session == session) & (df_spectrum.state==state)]
                    df.to_excel(dirname+name)
        for group in ['control', 'test']:
            for row, session in enumerate(sessions):
                for col, state in enumerate(states):
                    print(group, session, state)
                    data = df_spectrum[(df_spectrum.session == session) & (df_spectrum.group == group) & (df_spectrum.state == state)]
                    data = data.dropna()
                    data = (data[freqs]).values

                    print(data.shape)
                    plot = np.mean(data, axis=0)
                    # inf, med, sup = np.percentile(data, [25,50,75], axis=0)
                    # inf = np.array(inf, dtype = 'float64')
                    # sup = np.array(sup, dtype = 'float64')
                    # med = np.array(med, dtype = 'float64')

                    ax[row, col].plot(freqs,plot)

                    # ax[row, col].plot(freqs, med )
                    # ax[row, col].fill_between(x =freqs, y1 = inf, y2 = sup, alpha = .5 )


                # plt.show()



def get_clear_spectral_score_one_mouse_one_condition_day_night(mouse, condition, state, period = 'light'):
    debug = 0
    # print(condition, period)
    ds = xr.open_dataset(precompute_dir + '/spectrums/spectrum_scoring_{}.nc'.format(mouse))
    times = ds.coords['time_epochs'].values/3600
    # score = ds['score'].values
    dirname = precompute_dir + '/tdw_score/'
    ds_tdw = xr.open_dataset(dirname + 'tdw_score_{}.nc'.format(mouse))
    score = ds_tdw['new_score'].values

    if condition == 'bl0':
        if period == 'dark':
            t1 = 4
            mask = (times<t1)

    elif condition == 'bl1':
        if period == 'light':
            t1, t2 = 4, 16
            mask = (times>t1)&(times<t2)
        elif period == 'dark':
            t1, t2 = 16, 28
            mask = (times>=t1) & (times<t2)
    elif condition == 'bl2':
        if period == 'light':
            t1, t2 = 28, 40
            mask = (times>=t1) & (times<t2)
        elif period == 'dark':
            t1, t2 = 40, 52
            mask = (times>=t1) & (times<t2)
    elif condition == 'sd':
        if period == 'light':
            t1, t2 = 58, 64    #####caution 6 first hours is SD 48 to 54
            mask = (times>=t1) & (times<t2)
        elif period == 'dark':
            t1, t2 = 64, 76
            mask = (times>=t1) & (times<t2)
    elif condition == 'r1':
        if period == 'light':
            t1, t2 = 76, 88
            mask = (times>=t1) & (times<t2)
        elif period == 'dark':
            t1 = 88
            mask = (times>=t1)
    else :
        print('Condition does not exist !')

    # mask = (times>=t1) & (times<t2)
    score = score[mask]

    if state == 'w':
            score_no_artifact = (score == 'w' )|(score == 't')
    else :
        score_no_artifact = score == state
    score_no_transition = score_no_artifact.copy()
    count = 0
    event_length = []
    for num, value in enumerate(score_no_artifact) :
        if value:
            count +=1
        elif value == False:
            event_length.append(count)
            if count == 1 :
                score_no_transition[num-1] = False
            elif count == 2 :
                score_no_transition[num-1] = False
                score_no_transition[num-2] = False
                count = 0
            elif count > 2 :
                score_no_transition[num-1] = False
                score_no_transition[num-count] = False
            count = 0
    if count == 1 :
        score_no_transition[num-1] = False
    elif count == 2 :
        score_no_transition[num-1] = False
        score_no_transition[num-2] = False
        count = 0
    elif count > 2 :
        score_no_transition[num-1] = False
        score_no_transition[num-count] = False


    if debug :
        fig, ax = plt.subplots()
        s = ds.coords['time_epochs'].values.size
        for i in range(5):
            ax.axvline(i*s/(3600))
        h = ds.coords['time_epochs'].values/3600
        ax.plot(h, np.ones(s))
        ax.plot(h[mask], np.ones(s)[mask])
        ax.plot(h[mask], score == state)

        fig, ax = plt.subplots()
        ax.scatter(h[mask], score_no_artifact)
        ax.plot(h[mask], score_no_transition, color = 'green')
        # print(count)
        plt.show()
    if sum(score_no_transition) == 0:
        print('_________ No {} for {} during {} {}________'.format(state, mouse, condition, period))
    return score_no_transition, mask

def get_relative_spectrums_one_mouse_one_condition_day_night(mouse, condition, spectrum_method = 'somno', period = 'light'):
    debug = 0
    cleared_score_w, _= get_clear_spectral_score_one_mouse_one_condition_day_night(mouse, condition, 'w', period = period)
    cleared_score_r, _ = get_clear_spectral_score_one_mouse_one_condition_day_night(mouse, condition, 'r', period = period)
    cleared_score_t, _ = get_clear_spectral_score_one_mouse_one_condition_day_night(mouse, condition, 't', period = period)
    cleared_score_n, mask_condition = get_clear_spectral_score_one_mouse_one_condition_day_night(mouse, condition, 'n', period = period)

    ds = xr.open_dataset(precompute_dir + '/spectrums/spectrum_scoring_{}.nc'.format(mouse))

    spectrum = ds['{}_spectrum'.format(spectrum_method)].values
    freqs = ds.coords['freqs_{}'.format(spectrum_method)].values
    spectrum = spectrum[mask_condition, :]

    total_epochs = np.sum(cleared_score_w) + np.sum(cleared_score_r) + np.sum(cleared_score_n)+ np.sum(cleared_score_t)

    proportion_w = np.sum(cleared_score_w)/total_epochs
    proportion_r = np.sum(cleared_score_r)/total_epochs
    proportion_n = np.sum(cleared_score_n)/total_epochs
    proportion_t = np.sum(cleared_score_t)/total_epochs

    absolute_spectrum_w = np.nanmedian(spectrum[cleared_score_w, :], axis = 0)
    absolute_spectrum_r = np.nanmedian(spectrum[cleared_score_r, :], axis = 0)
    absolute_spectrum_n = np.nanmedian(spectrum[cleared_score_n, :], axis = 0)
    absolute_spectrum_t = np.nanmedian(spectrum[cleared_score_t, :], axis = 0)

    df = (freqs[11]-freqs[1])/10


    if np.sum(cleared_score_w) !=0:
        sum_w = np.sum(absolute_spectrum_w)*df
    else:
        sum_w =0

    if np.sum(cleared_score_t) !=0:
        sum_t = np.sum(absolute_spectrum_t)*df
    else:
        sum_t =0

    if np.sum(cleared_score_r) !=0:
        sum_r = np.sum(absolute_spectrum_r)*df
    else:
        sum_r =0

    if np.sum(cleared_score_n) !=0:
        sum_n = np.sum(absolute_spectrum_n)*df
    else:
        sum_n = 0


    #########CAUTION power spectrum is area under curv#########
    # auc_w = auc(freqs, absolute_spectrum_w)
    # auc_r = auc(freqs, absolute_spectrum_r)
    # auc_n = auc(freqs, absolute_spectrum_n)
    # ref_values = proportion_w*auc_w + proportion_r*auc_r + proportion_n*auc_n


    ####FALSE
    # mean_w = np.mean(absolute_spectrum_w)
    # mean_r = np.mean(absolute_spectrum_r)
    # mean_n = np.mean(absolute_spectrum_n)
    # ref_values = proportion_w*mean_w + proportion_r*mean_r + proportion_n*mean_n


    #####ALi's
    # sum_w = np.sum(absolute_spectrum_w)*df
    # sum_r = np.sum(absolute_spectrum_r)*df
    # sum_n = np.sum(absolute_spectrum_n)*df
    ref_values = (proportion_w*sum_w + proportion_r*sum_r + proportion_n*sum_n +  proportion_t*sum_t)
    # print(ref_values)

    ref_values = (proportion_w*sum_w + proportion_r*sum_r + proportion_n*sum_n)
    # print(ref_values)
    # print(ref_values/df)
    # print(1/df)
    # exit()

    relative_spectrum_w = 100*absolute_spectrum_w / ref_values
    relative_spectrum_r = 100*absolute_spectrum_r / ref_values
    relative_spectrum_n = 100*absolute_spectrum_n / ref_values
    relative_spectrum_t = 100*absolute_spectrum_t / ref_values

    # fig, ax = plt.subplots()
    # ax.plot(relative_spectrum_n)
    # ax.plot(relative_spectrum_n/df)
    # plt.show()
    # exit()
    # print(ref_values)
    # fig, ax = plt.subplots()
    # global_spec = np.mean(spectrum[cleared_score_w + cleared_score_r + cleared_score_n], axis = 0)
    # print(auc(freqs,global_spec))
    # ax.plot(freqs,global_spec )
    # ax.plot(freqs,global_spec/ref_values )
    #
    #
    # fig, ax = plt.subplots()
    # ax.plot(freqs, relative_spectrum_w)
    # ax.plot(freqs, relative_spectrum_r)
    # ax.plot(freqs, relative_spectrum_n)
    # plt.show()

    if debug == True:
        fig, ax = plt.subplots()
        ax.plot(absolute_spectrum_w, color = 'blue')
        ax.plot(absolute_spectrum_r, color = 'blue')
        ax.plot(absolute_spectrum_n, color = 'blue')
        ax.plot(absolute_spectrum_t, color = 'blue')

        ax.plot(relative_spectrum_w, color = 'green')
        ax.plot(relative_spectrum_r, color = 'green', alpha = .75)
        ax.plot(relative_spectrum_n, color = 'green', alpha = .50)
        ax.plot(relative_spectrum_t, color = 'green', alpha = .25)
        plt.show()
    return ref_values, {'w' : relative_spectrum_w, 'r':relative_spectrum_r, 'n' :relative_spectrum_n, 't' :relative_spectrum_t}



def plot_spectrum_compare_day_night(spectrum_method = 'somno', period = 'light', exclude_mice = False):
    date_ref = pd.read_excel(work_dir + 'datetime_reference_miRNA.xlsx', index_col = 0)
    if period == 'light':
        sessions = ['bl1', 'bl2', 'sd', 'r1' ]
    elif period == 'dark':
        sessions = ['bl0', 'bl1', 'bl2', 'sd', 'r1' ]

    for miRNA in ['128', '137', '665']:
        mice_control = get_mice(miRNA, 'control', exclude_mice)
        mice_test = get_mice(miRNA, 'test', exclude_mice)

        groups = {'control' :mice_control , 'test' : mice_test}
        ds = xr.open_dataset(precompute_dir + '/spectrums/spectrum_scoring_B2884.nc')   #B2576
        freqs = ds.coords['freqs_{}'.format(spectrum_method)].values
        freqs = np.round(freqs, 2).tolist()

        indices = []
        for state in ['w', 'n', 'r', 't']:
            for session in sessions:
                indices+=([i + '_'+session+ '_'+state for i in groups['test']+ groups['control']])
        df_spectrum = pd.DataFrame(index = indices, columns = ['group', 'session', 'state', 'ref_values'] + freqs)



        for group in groups:
            for mouse in groups[group]:
                ds = xr.open_dataset(precompute_dir + '/spectrums/spectrum_scoring_{}.nc'.format(mouse))
                freqs = ds.coords['freqs_{}'.format(spectrum_method)].values
                freqs = np.round(freqs, 2)

                ZT = ds.coords['time_epochs'].values
                print('collecting relative spectrums ', mouse)
                states = ['w', 'n', 'r', 't']
                for session in sessions:
                    if period == 'dark' and date_ref.at['MTA-{}'.format(mouse), 'short_end'] == 'yes' and session == 'r1':
                        for state in states :
                            df_spectrum.at['{}_{}_{}'.format(mouse,session, state), 'state'] = state
                            df_spectrum.at['{}_{}_{}'.format(mouse,session, state), 'session'] = session
                            df_spectrum.at['{}_{}_{}'.format(mouse,session, state), 'group'] = group
                            df_spectrum.at['{}_{}_{}'.format(mouse,session, state), 'ref_values'] = np.nan
                            empty = np.ones(freqs.size)
                            empty[:] = np.nan
                            df_spectrum.at['{}_{}_{}'.format(mouse,session, state), 4:] = empty
                    else :
                        print(session, mouse, period)
                        ref_values, relative_spectrums = get_relative_spectrums_one_mouse_one_condition_day_night(mouse, session, spectrum_method, period )
                        # print(relative_spectrums[state].size)
                        for state in states :
                            df_spectrum.at['{}_{}_{}'.format(mouse,session, state), 'state'] = state
                            df_spectrum.at['{}_{}_{}'.format(mouse,session, state), 'session'] = session
                            df_spectrum.at['{}_{}_{}'.format(mouse,session, state), 'group'] = group
                            df_spectrum.at['{}_{}_{}'.format(mouse,session, state), 'ref_values'] = ref_values
                            df_spectrum.at['{}_{}_{}'.format(mouse,session, state), 4:] = relative_spectrums[state]
        # print(df_spectrum)
            # df_spectrum = df_spectrum.dropna()
        if period == 'light':
            Nrows =4
        elif period == 'dark':
            Nrows=5
        fig,ax = plt.subplots(nrows = Nrows, ncols= 4, sharex = True, sharey = True)
        fig_diff, ax_diff = plt.subplots(nrows = Nrows, ncols= 4, sharex = True, sharey = True)

        fig.suptitle(miRNA +' -- Spectrum per state -- ' + period + ' -exclude mice ='+str(exclude_mice))
        fig_diff.suptitle(miRNA + ' -- Diff DCR - ctrl')

        for i in range(Nrows):
            ax[i,0].set_ylabel(sessions[i])
            ax_diff[i,0].set_ylabel(sessions[i])

        for row, session in enumerate(sessions):
            for col, state in enumerate(states):
                # print(df_spectrum[(df_spectrum.session == session)& (df_spectrum.state == state)])
                data_dcr = df_spectrum[(df_spectrum.session == session) & (df_spectrum.group == 'test') & (df_spectrum.state == state)]
                data_ctrl = df_spectrum[(df_spectrum.session == session) & (df_spectrum.group == 'control') & (df_spectrum.state == state)]
                # print(data_dcr)
                # print(data_dcr.columns)
                # print(freqs)
                data_dcr = data_dcr.dropna()
                data_ctrl = data_ctrl.dropna()
                data_dcr = (data_dcr[freqs]).values
                data_ctrl = (data_ctrl[freqs]).values
                # print()
                # for i, maxi in enumerate(np.max(data_dcr, axis = 1)):
                #     data_dcr[i] = data_dcr[i]/maxi
                # for i, maxi in enumerate(np.max(data_ctrl, axis = 1)):
                #     data_ctrl[i] = data_ctrl[i]/maxi
                diff = np.mean(data_dcr, axis=0) - np.mean(data_ctrl, axis=0)
                # ax[row, col].plot(data.T, color = 'black', alpha = .1)
                # ax[row, col].plot(data.mean(axis = 0), color = 'red')
                ax_diff[row, col].plot(freqs, diff )
        for group in groups:
            if exclude_mice:
                dirname = excel_dir + '/RT_PCR_excluded_mice/power_spectrum/{}/{}/'.format(miRNA,group)
            else :
                dirname = excel_dir + '/power_spectrum/{}/{}/'.format(miRNA,group)

            if not os.path.exists(dirname):
                os.makedirs(dirname)

            for session in sessions:
                for state in states:
                    name = '{}_spectrum_{}_{}_{}_{}_exclude_mice_{}.xlsx'.format(miRNA,spectrum_method,state, session,period, exclude_mice)
                    df = df_spectrum[(df_spectrum.group == group) & (df_spectrum.session == session) & (df_spectrum.state==state)]
                    df.to_excel(dirname+name)
        for group in ['control', 'test']:
            for row, session in enumerate(sessions):
                for col, state in enumerate(states):
                    print(group, session, state)
                    data = df_spectrum[(df_spectrum.session == session) & (df_spectrum.group == group) & (df_spectrum.state == state)]
                    data = data.dropna()
                    data = (data[freqs]).values

                    print(data.shape)
                    plot = np.mean(data, axis=0)
                    # inf, med, sup = np.percentile(data, [25,50,75], axis=0)
                    # inf = np.array(inf, dtype = 'float64')
                    # sup = np.array(sup, dtype = 'float64')
                    # med = np.array(med, dtype = 'float64')

                    ax[row, col].plot(freqs,plot)

                    # ax[row, col].plot(freqs, med )
                    # ax[row, col].fill_between(x =freqs, y1 = inf, y2 = sup, alpha = .5 )


                # plt.show()
    # plt.show()
def plot_spectrum_by_mouse(spectrum_method = 'somno'):
    date_ref = pd.read_excel(work_dir + 'datetime_reference_miRNA.xlsx', index_col = 0)
    for miRNA in ['128', '137', '665']:
        groups = {'control' : get_mice(miRNA, 'control'), 'test' : get_mice(miRNA, 'test')}
        for group in groups:
            for mouse in groups[group]:
                ds = xr.open_dataset(precompute_dir + '/spectrums/spectrum_scoring_{}.nc'.format(mouse))
                freqs = ds.coords['freqs_{}'.format(spectrum_method)].values
                print(mouse)
                period_color = {'dark' : 'black', 'light': 'grey'}

                states = ['w', 't', 'n', 'r']
                sessions = ['bl1', 'bl2', 'sd', 'r1']
                fig, ax = plt.subplots(nrows = 4, ncols= 4, sharex = True)
                ax[0,0].set_title('wake')
                ax[0,1].set_title('tdw')
                ax[0,2].set_title('nrem')
                ax[0,3].set_title('rem')
                ax[0,0].set_ylabel('bl1')
                ax[1,0].set_ylabel('bl2')
                ax[2,0].set_ylabel('sd')
                ax[3,0].set_ylabel('r1')

                fig.suptitle('{} - {} - {} - {} - light:{}, dark:{}'.format(mouse, spectrum_method, miRNA, group,period_color['light'],period_color['dark']))
                for period in ['dark', 'light']:
                    for col, session in enumerate(sessions):
                        if period == 'dark' and date_ref.at['MTA-{}'.format(mouse), 'short_end'] == 'yes' and session == 'r1':
                            continue
                        else :
                            ref_values, relative_spectrums = get_relative_spectrums_one_mouse_one_condition_day_night(mouse, session, spectrum_method, period )
                            for row, state in enumerate(states) :
                                ax[col,row].plot(freqs, relative_spectrums[state], color = period_color[period])
                                ax[col,row].set_xlim(0,15)
                dirname = work_dir+'/pyFig/power_spectrum/{}/{}/'.format(spectrum_method, miRNA, group)
                if not os.path.exists(dirname):
                    os.makedirs(dirname)
                plt.savefig(dirname + mouse+'.png')
                plt.close()
def plot_spectrum_one_mouse(mouse, spectrum_method = 'welch'):
    date_ref = pd.read_excel(work_dir + 'datetime_reference_miRNA.xlsx', index_col = 0)
    miRNA, genotype = get_mouse_info(mouse)
    ds = xr.open_dataset(precompute_dir + '/spectrums/spectrum_scoring_{}.nc'.format(mouse))
    freqs = ds.coords['freqs_{}'.format(spectrum_method)].values
    print(mouse)
    period_color = {'dark' : 'black', 'light': 'grey'}

    states = ['w','t', 'n', 'r']
    sessions = ['bl1', 'bl2', 'sd', 'r1']
    fig, ax = plt.subplots(nrows = 4, ncols= 3, sharex = True)
    ax[0,0].set_title('wake')
    ax[0,1].set_title('tdw')
    ax[0,2].set_title('nrem')
    ax[0,3].set_title('rem')
    ax[0,0].set_ylabel('bl1')
    ax[1,0].set_ylabel('bl2')
    ax[2,0].set_ylabel('sd')
    ax[3,0].set_ylabel('r1')

    fig.suptitle('{} - {} - {} - {} - light:{}, dark:{}'.format(mouse, spectrum_method, miRNA, genotype,period_color['light'],period_color['dark']))
    for period in ['dark', 'light']:
        for col, session in enumerate(sessions):
            if period == 'dark' and date_ref.at['MTA-{}'.format(mouse), 'short_end'] == 'yes' and session == 'r1':
                continue
            else :
                ref_values, relative_spectrums = get_relative_spectrums_one_mouse_one_condition_day_night(mouse, session, spectrum_method, period )
                for row, state in enumerate(states) :
                    ax[col,row].plot(freqs, relative_spectrums[state], color = period_color[period])
                    ax[col,row].set_xlim(0,15)
    dirname = work_dir+'/pyFig/power_spectrum/{}/{}/'.format(spectrum_method, miRNA, genotype)
    if not os.path.exists(dirname):
        os.makedirs(dirname)
    plt.savefig(dirname + mouse+'.png')
    plt.close()


def get_clear_4firstdarkhour_spectral_score_one_mouse_one_condition_day_night(mouse, condition, state):
    debug = 0
    # print(condition, period)
    ds = xr.open_dataset(precompute_dir + '/spectrums/spectrum_scoring_{}.nc'.format(mouse))
    times = ds.coords['time_epochs'].values/3600
    # score = ds['score'].values
    dirname = precompute_dir + '/tdw_score/'
    ds_tdw = xr.open_dataset(dirname + 'tdw_score_{}.nc'.format(mouse))
    score = ds_tdw['new_score'].values


    if condition == 'dark1':
        t1, t2 = 16, 20
        mask = (times>=t1) & (times<t2)
    elif condition == 'dark2':
        t1, t2 = 40, 44
        mask = (times>=t1) & (times<t2)
    elif condition == 'post_sd':
        t1, t2 = 58, 62
        mask = (times>=t1) & (times<t2)
    elif condition == 'dark3':
        t1, t2 = 64, 68
        mask = (times>=t1) & (times<t2)
    elif condition == 'dark4':
        t1, t2 = 88, 92
        mask = (times>=t1) & (times<t2)
    else :
        print('Condition does not exist !')

    # mask = (times>=t1) & (times<t2)
    score = score[mask]

    if state == 'w':
            score_no_artifact = (score == 'w' )|(score == 't')
    else :
        score_no_artifact = score == state
    score_no_transition = score_no_artifact.copy()
    count = 0
    event_length = []
    for num, value in enumerate(score_no_artifact) :
        if value:
            count +=1
        elif value == False:
            event_length.append(count)
            if count == 1 :
                score_no_transition[num-1] = False
            elif count == 2 :
                score_no_transition[num-1] = False
                score_no_transition[num-2] = False
                count = 0
            elif count > 2 :
                score_no_transition[num-1] = False
                score_no_transition[num-count] = False
            count = 0
    if count == 1 :
        score_no_transition[num-1] = False
    elif count == 2 :
        score_no_transition[num-1] = False
        score_no_transition[num-2] = False
        count = 0
    elif count > 2 :
        score_no_transition[num-1] = False
        score_no_transition[num-count] = False


    if debug :
        fig, ax = plt.subplots()
        s = ds.coords['time_epochs'].values.size
        for i in range(5):
            ax.axvline(i*s/(3600))
        h = ds.coords['time_epochs'].values/3600
        ax.plot(h, np.ones(s))
        ax.plot(h[mask], np.ones(s)[mask])
        ax.plot(h[mask], score == state)

        fig, ax = plt.subplots()
        ax.scatter(h[mask], score_no_artifact)
        ax.plot(h[mask], score_no_transition, color = 'green')
        # print(count)
        plt.show()
    if sum(score_no_transition) == 0:
        print('_________ No {} for {} during {}________'.format(state, mouse, condition))
    return score_no_transition, mask

def get_relative_4firstdarkhour_spectrums_one_mouse_one_condition_day_night_Sha(mouse, condition, spectrum_method = 'somno', auc = False, agg_power = 'mean'):
    cleared_score_w, _= get_clear_4firstdarkhour_spectral_score_one_mouse_one_condition_day_night(mouse, condition, 'w')
    cleared_score_r, _ = get_clear_4firstdarkhour_spectral_score_one_mouse_one_condition_day_night(mouse, condition, 'r')
    cleared_score_t, _ = get_clear_4firstdarkhour_spectral_score_one_mouse_one_condition_day_night(mouse, condition, 't')
    cleared_score_n, mask_condition = get_clear_4firstdarkhour_spectral_score_one_mouse_one_condition_day_night(mouse, condition, 'n')

    ds = xr.open_dataset(precompute_dir + '/spectrums/spectrum_scoring_{}.nc'.format(mouse))

    spectrum = ds['{}_spectrum'.format(spectrum_method)].values
    freqs = ds.coords['freqs_{}'.format(spectrum_method)].values
    spectrum = spectrum[mask_condition, :]

    # total_epochs = np.sum(cleared_score_w) + np.sum(cleared_score_r) + np.sum(cleared_score_n)+ np.sum(cleared_score_t)
    #
    # proportion_w = np.sum(cleared_score_w)/total_epochs
    # proportion_r = np.sum(cleared_score_r)/total_epochs
    # proportion_n = np.sum(cleared_score_n)/total_epochs
    # proportion_t = np.sum(cleared_score_t)/total_epochs


    if agg_power == 'median':
        absolute_spectrum_w = np.nanmedian(spectrum[cleared_score_w, :], axis = 0)
        absolute_spectrum_r = np.nanmedian(spectrum[cleared_score_r, :], axis = 0)
        absolute_spectrum_n = np.nanmedian(spectrum[cleared_score_n, :], axis = 0)
        absolute_spectrum_t = np.nanmedian(spectrum[cleared_score_t, :], axis = 0)
    elif agg_power=='mean':
        absolute_spectrum_w = np.nanmean(spectrum[cleared_score_w, :], axis = 0)
        absolute_spectrum_r = np.nanmean(spectrum[cleared_score_r, :], axis = 0)
        absolute_spectrum_n = np.nanmean(spectrum[cleared_score_n, :], axis = 0)
        absolute_spectrum_t = np.nanmean(spectrum[cleared_score_t, :], axis = 0)

    ref_values = get_ref_values_Sha(mouse,spectrum_method, period = 'light', auc = auc, agg_power = agg_power)

    relative_spectrum_w = 100*absolute_spectrum_w / ref_values
    relative_spectrum_r = 100*absolute_spectrum_r / ref_values
    relative_spectrum_n = 100*absolute_spectrum_n / ref_values
    relative_spectrum_t = 100*absolute_spectrum_t / ref_values

    return ref_values, {'w' : relative_spectrum_w, 'r':relative_spectrum_r, 'n' :relative_spectrum_n, 't' :relative_spectrum_t}




def plot_four_first_dark_hour_psd(spectrum_method = 'somno', exclude_mice = False, auc = False, agg_power = 'mean'):
    sessions = ['dark1', 'dark2', 'post_sd', 'dark3', 'dark4']
    date_ref = pd.read_excel(work_dir + 'datetime_reference_miRNA.xlsx', index_col = 0)

    for miRNA in ['128', '137', '665']:
        mice_control = get_mice(miRNA, 'control', exclude_mice)
        mice_test = get_mice(miRNA, 'test', exclude_mice)

        groups = {'control' :mice_control , 'test' : mice_test}
        ds = xr.open_dataset(precompute_dir + '/spectrums/spectrum_scoring_B2884.nc')   #B2576
        freqs = ds.coords['freqs_{}'.format(spectrum_method)].values
        freqs = np.round(freqs, 2).tolist()

        indices = []
        for state in ['w', 'n', 'r', 't']:
            for session in sessions:
                indices+=([i + '_'+session+ '_'+state for i in groups['test']+ groups['control']])
        df_spectrum = pd.DataFrame(index = indices, columns = ['group', 'session', 'state', 'ref_values'] + freqs)



        for group in groups:
            for mouse in groups[group]:
                ds = xr.open_dataset(precompute_dir + '/spectrums/spectrum_scoring_{}.nc'.format(mouse))
                freqs = ds.coords['freqs_{}'.format(spectrum_method)].values
                freqs = np.round(freqs, 2)

                ZT = ds.coords['time_epochs'].values
                print('collecting relative spectrums ', mouse)
                states = ['w', 'n', 'r', 't']
                for session in sessions:
                    if session == 'dark4' and date_ref.at['MTA-{}'.format(mouse), 'short_end'] == 'yes' :
                        for state in states :
                            df_spectrum.at['{}_{}_{}'.format(mouse,session, state), 'state'] = state
                            df_spectrum.at['{}_{}_{}'.format(mouse,session, state), 'session'] = session
                            df_spectrum.at['{}_{}_{}'.format(mouse,session, state), 'group'] = group
                            df_spectrum.at['{}_{}_{}'.format(mouse,session, state), 'ref_values'] = np.nan
                            empty = np.ones(freqs.size)
                            empty[:] = np.nan
                            df_spectrum.at['{}_{}_{}'.format(mouse,session, state), 4:] = empty
                    else :
                        print(session, mouse)
                        ref_values, relative_spectrums = get_relative_4firstdarkhour_spectrums_one_mouse_one_condition_day_night_Sha(mouse, session, spectrum_method, auc, agg_power )
                        # print(relative_spectrums[state].size)
                        for state in states :
                            df_spectrum.at['{}_{}_{}'.format(mouse,session, state), 'state'] = state
                            df_spectrum.at['{}_{}_{}'.format(mouse,session, state), 'session'] = session
                            df_spectrum.at['{}_{}_{}'.format(mouse,session, state), 'group'] = group
                            df_spectrum.at['{}_{}_{}'.format(mouse,session, state), 'ref_values'] = ref_values
                            df_spectrum.at['{}_{}_{}'.format(mouse,session, state), 4:] = relative_spectrums[state]
        # print(df_spectrum)
            # df_spectrum = df_spectrum.dropna()

        Nrows=len(sessions)
        fig,ax = plt.subplots(nrows = Nrows, ncols= 4, sharex = True, sharey = True)
        fig_diff, ax_diff = plt.subplots(nrows = Nrows, ncols= 4, sharex = True, sharey = True)

        fig.suptitle(miRNA +' -- Spectrum per state -- exclude mice ='+str(exclude_mice))
        fig_diff.suptitle(miRNA + ' -- Diff DCR - ctrl')

        for i in range(Nrows):
            ax[i,0].set_ylabel(sessions[i])
            ax_diff[i,0].set_ylabel(sessions[i])

        for row, session in enumerate(sessions):
            for col, state in enumerate(states):
                # print(df_spectrum[(df_spectrum.session == session)& (df_spectrum.state == state)])
                data_dcr = df_spectrum[(df_spectrum.session == session) & (df_spectrum.group == 'test') & (df_spectrum.state == state)]
                data_ctrl = df_spectrum[(df_spectrum.session == session) & (df_spectrum.group == 'control') & (df_spectrum.state == state)]
                # print(data_dcr)
                # print(data_dcr.columns)
                # print(freqs)
                data_dcr = data_dcr.dropna()
                data_ctrl = data_ctrl.dropna()
                data_dcr = (data_dcr[freqs]).values
                data_ctrl = (data_ctrl[freqs]).values
                # print()
                # for i, maxi in enumerate(np.max(data_dcr, axis = 1)):
                #     data_dcr[i] = data_dcr[i]/maxi
                # for i, maxi in enumerate(np.max(data_ctrl, axis = 1)):
                #     data_ctrl[i] = data_ctrl[i]/maxi
                diff = np.mean(data_dcr, axis=0) - np.mean(data_ctrl, axis=0)
                # ax[row, col].plot(data.T, color = 'black', alpha = .1)
                # ax[row, col].plot(data.mean(axis = 0), color = 'red')
                ax_diff[row, col].plot(freqs, diff )
        for group in groups:
            if exclude_mice:
                dirname = excel_dir + '/RT_PCR_excluded_mice/fourfirstdarkhour_power_spectrum_baseline_auc_{}_agg_power_{}/{}/{}/'.format(auc, agg_power,miRNA,group)
            else :
                dirname = excel_dir + '/fourfirstdarkhour_power_spectrum_baseline_auc_{}_agg_power_{}/{}/{}/'.format(auc, agg_power,miRNA,group)

            if not os.path.exists(dirname):
                os.makedirs(dirname)

            for session in sessions:
                for state in states:
                    name = '{}_spectrum_{}_{}_{}_exclude_mice_{}.xlsx'.format(miRNA,spectrum_method,state, session, exclude_mice)
                    df = df_spectrum[(df_spectrum.group == group) & (df_spectrum.session == session) & (df_spectrum.state==state)]
                    df.to_excel(dirname+name)
        for group in ['control', 'test']:
            for row, session in enumerate(sessions):
                for col, state in enumerate(states):
                    print(group, session, state)
                    data = df_spectrum[(df_spectrum.session == session) & (df_spectrum.group == group) & (df_spectrum.state == state)]
                    data = data.dropna()
                    data = (data[freqs]).values

                    print(data.shape)
                    plot = np.mean(data, axis=0)
                    # inf, med, sup = np.percentile(data, [25,50,75], axis=0)
                    # inf = np.array(inf, dtype = 'float64')
                    # sup = np.array(sup, dtype = 'float64')
                    # med = np.array(med, dtype = 'float64')

                    ax[row, col].plot(freqs,plot)

                    # ax[row, col].plot(freqs, med )
                    # ax[row, col].fill_between(x =freqs, y1 = inf, y2 = sup, alpha = .5 )


                # plt.show()




if __name__ == '__main__':
    # mouse = 'B3512'
    # mouse = 'B2534'
    # mouse = 'B2533'
    # mouse = 'B2767'
    # mouse = 'B2761'
    # mouse = 'B2700'
    # mouse = 'B2762'
    # mouse = 'B2763'
    # mouse = 'B3072'
    # mouse = 'B3140'
    # mouse = 'B3513'
    # mouse = 'B3512'
    # mouse = 'B3766'
    # mouse = 'B4112'
    # mouse = 'B4113'
    # mouse = 'B4975'
    # mouse = 'B4976'
    # mouse = 'B4907'


    # mouse = 'B4904'
    # mouse  ='B2762'
    # mouse  ='B2763'
    # mouse  ='B3072'
    # mouse  ='B3140'
    # mouse  ='B3513'
    # mouse = 'B4977'
    # mouse = 'B2884'
    mouse = 'B2576'
    # get_ref_values_Sha(mouse, 'welch', 'light')
    # compute_ref_values_baseline_Sha('665', 'control')
    plot_four_first_dark_hour_psd(spectrum_method = 'welch', exclude_mice = False, auc = False, agg_power = 'mean')
    # plt.show()
    plot_four_first_dark_hour_psd(spectrum_method = 'welch', exclude_mice = True, auc = False, agg_power = 'mean')
    # plt.show()
    # plot_spectrum_compare_day_night_Sha(spectrum_method = 'welch', period = 'dark', auc = False, agg_power = 'mean' )
    # plot_spectrum_compare_day_night_Sha(spectrum_method = 'welch', period = 'light', auc = False, agg_power = 'mean' )
    #
    # plot_spectrum_compare_day_night_Sha(spectrum_method = 'welch', period = 'dark', exclude_mice = True, auc = False, agg_power = 'mean' )
    # plot_spectrum_compare_day_night_Sha(spectrum_method = 'welch', period = 'light',exclude_mice = True,  auc = False, agg_power = 'mean' )

    # plot_spectrum_compare_day_night(spectrum_method = 'welch', period = 'dark')
    # plt.show()
    # exit()

        # group = 'Control'
    # rec = 'b2'
    # compute_all()
    # store_scoring_and_spectrums_one_mouse_one_session_day_(group, mouse)
    # get_clear_spectral_score_one_mouse_one_condition(mouse, 'bl1', 'r')
    # get_relative_spectrums_one_mouse_one_condition_day_night(mouse, 'bl1', 'welch')

    # plot_spectrum_compare_global(spectrum_method = 'somno')
    # plot_spectrum_compare_day_night(spectrum_method = 'somno', period = 'dark')
    # plot_spectrum_compare_day_night(spectrum_method = 'somno', period = 'dark')
    # plot_spectrum_compare_day_night(spectrum_method = 'welch', period = 'light')
    # plot_spectrum_compare_day_night(spectrum_method = 'welch', period = 'dark')
    # plot_spectrum_compare_day_night(spectrum_method = 'welch', period = 'light', exclude_mice = True)
    # plot_spectrum_compare_day_night(spectrum_method = 'welch', period = 'dark', exclude_mice = True)


    # get_clear_spectral_score_one_mouse_one_condition_day_night(mouse, 'bl1', 'r', 'dark')
    # get_relative_spectrums_one_mouse_one_condition_day_night(mouse = mouse, condition ='bl1', spectrum_method = 'somno', period = 'dark')
    plot_spectrum_by_mouse('welch')
    # plot_spectrum_one_mouse(mouse)
    plt.show()
