from configuration import *
from select_mice_cata_Malo import get_mice
import multiprocessing
from joblib import Parallel, delayed
from sklearn.metrics import auc

local_path = os.path.dirname(os.path.realpath(__file__))
print(local_path)


def plot_ratio_compare(spectrum_method = 'somno'):
    control_list = get_mice('Control')
    DCR_list = get_mice('DCR-HCRT')
    groups = {'Control' : control_list, 'DCR-HCRT' : DCR_list}
    ds = xr.open_dataset(precompute_dir + '/spectrums/spectrum_scoring_B2533.nc')
    freqs = ds.coords['freqs_{}'.format(spectrum_method)].values.tolist()
    df_ratio = pd.DataFrame(index = control_list+ DCR_list, columns = np.arange(96))
    indices = []
    for state in ['w', 'n', 'r', 'a']:
        for session in ['bl1', 'bl2', 'sd', 'r1' ]:
            indices+=([i + '_'+session+ '_'+state for i in control_list+ DCR_list])
    df_spectrum = pd.DataFrame(index = indices, columns = ['group', 'session', 'state'] + freqs)

    fig,ax = plt.subplots(nrows = 4, ncols= 4, sharex = True, sharey = True)
    fig. suptitle('Spectrum per state')
    for group in groups:
        for mouse in groups[group]:
            ds = xr.open_dataset(precompute_dir + '/spectrums/spectrum_scoring_{}.nc'.format(mouse))
            print(ds)
            # mask = [0:int(12*3600*sr)]
            score = ds['score'].values
            score_behavior = score.copy()
            for f, n in zip(['1', '2', '3'], ['w', 'n', 'r']) :
                score_behavior = np.where(score_behavior == f, n, score_behavior)

            freqs = ds.coords['freqs_{}'.format(spectrum_method)].values
            times = ds.coords['epochs'].values
            t_start = 8.
            real_times = times + t_start

            spectrums = ds['{}_spectrum'.format(spectrum_method)].values
            # print(spectrums.values[np.isnan(spectrums.values)])
            # print(spectrums.max(dims = 'freqs_{}'.format(spectrum_method)))
            # spectrums = (spectrums/spectrums.max(dims = 'freqs_{}'.format(spectrum_method))).values
            # print(np.max(spectrums, axis = 1).shape)

            # spectrums = spectrums/np.max(spectrums, axis = 1)
            # spectrums = spectrums*freqs
            sr =200

            delta_power = ds['{}_delta_power'.format(spectrum_method)].values
            states = ['w', 'n', 'r', 'a']

            hours = np.arange(int(times.size*4/3600))
            delta_power_by_hour = np.zeros(hours.size)
            std_delta_power_by_hour = delta_power_by_hour.copy()
            for h in hours[:-1]:
                subdata = delta_power[int(3600/4*h):int(3600/4*(h+1))]
                subscore = score[int(3600/4*h):int(3600/4*(h+1))]
                # artifact_free = (subscore != '1') & (subscore != '2') & (subscore != '3')
                nrem = subscore == 'n'
                subdata = subdata[nrem]
                m = np.mean(subdata)
                std = np.std(subdata)
                delta_power_by_hour[h] = m
                std_delta_power_by_hour[h] = std
            norm_delta_power_by_hour = delta_power_by_hour/max(delta_power_by_hour) ####A REPRENDRE
            df_ratio.at[mouse, :] = norm_delta_power_by_hour

            mydict = {'bl1' : (0, 24), 'bl2' : (24,48), 'sd':(48,72), 'r1':(78, 102)}
            for s, state in enumerate(states) :
                for row, phase in enumerate(mydict):
                    ind1, ind2 = mydict[phase][0]*900, mydict[phase][1]*900
                    score_phase = score[ind1 : ind2]
                    spectrum_phase = spectrums[ind1 :ind2,:]
                    score_no_transition = score_phase == state
                    count = 0
                    event_length = []
                    for num, value in enumerate(score_no_transition) :
                        if value:
                            count +=1
                        elif value == False:
                            event_length.append(count)
                            if count == 1:
                                score_no_transition[num-1] = False
                            count = 0
                    state_spectrum = spectrum_phase[score_no_transition]
                    if state_spectrum.shape[0] <= 1:
                        continue

                    m = state_spectrum.mean(axis=0)
                    df_spectrum.at['{}_{}_{}'.format(mouse,phase, state), 'state'] = state
                    df_spectrum.at['{}_{}_{}'.format(mouse,phase, state), 'session'] = phase
                    df_spectrum.at['{}_{}_{}'.format(mouse,phase, state), 'group'] = group
                    df_spectrum.at['{}_{}_{}'.format(mouse,phase, state), 3:] = m
        # df_spectrum = df_spectrum.dropna()
    fig_diff, ax_diff = plt.subplots(nrows = 4, ncols= 4, sharex = True, sharey = True)
    fig_diff.suptitle('Diff DCR - ctrl')
    for row, session in enumerate(mydict):
        for col, state in enumerate(states):
            if state == 'a':
                continue
            # print(df_spectrum[(df_spectrum.session == session)& (df_spectrum.state == state)])
            data_dcr = df_spectrum[(df_spectrum.session == session) & (df_spectrum.group == 'DCR-HCRT') & (df_spectrum.state == state)]
            data_ctrl = df_spectrum[(df_spectrum.session == session) & (df_spectrum.group == 'Control') & (df_spectrum.state == state)]
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

    for group in ['Control', 'DCR-HCRT']:
        for row, session in enumerate(mydict):
            for col, state in enumerate(states):
                if group == 'Control' and state == 'a':
                    continue
                print(group, session, state)
                # print(df_spectrum[(df_spectrum.session == session)& (df_spectrum.state == state)])
                data = df_spectrum[(df_spectrum.session == session) & (df_spectrum.group == group) & (df_spectrum.state == state)]
                data = (data[freqs]).values
                # for i, maxi in enumerate(np.max(data, axis = 1)):
                #     data[i] = data[i]/maxi
                plot = np.mean(data, axis=0)
                inf, med, sup = np.percentile(data, [25,50,75], axis=0)
                inf = np.array(inf, dtype = 'float64')
                sup = np.array(sup, dtype = 'float64')
                med = np.array(med, dtype = 'float64')
                # print(np.isfinite(inf))
                # ax[row, col].plot(data.T, color = 'black', alpha = .1)
                # ax[row, col].plot(data.mean(axis = 0), color = 'red')
                ax[row, col].plot(freqs, med )
                ax[row, col].fill_between(x =freqs, y1 = inf, y2 = sup, alpha = .5 )


                # plt.show()
    plt.show()
    exit()


    fig, ax = plt.subplots()
    ax.plot(hours+t_start, delta_power_by_hour, color = 'darkgreen')
    ax.fill_between(hours+t_start, delta_power_by_hour-std_delta_power_by_hour, delta_power_by_hour+std_delta_power_by_hour, color = 'darkgreen', alpha = .5)
    fig.suptitle('Delta power for {}'.format(mouse))

    fig, ax = plt.subplots()
    ax.plot(times, delta_power, color = 'black')

    fig,ax = plt.subplots(nrows = 4, ncols= 4, sharex = True, sharey = True)
    # fig,ax = plt.subplots(nrows = 2, ncols= 4, sharex = True)
    mydict = {'bl1' : (0, 24), 'bl2' : (24,48), 'sd':(48,72), 'r1':(78, 102)}
    for s, state in enumerate(states) :
        for row, phase in enumerate(mydict):
            ind1, ind2 = mydict[phase][0]*900, mydict[phase][1]*900
            print(ind1,ind2)
            score_phase = score[ind1 : ind2]
            spectrum_phase = spectrums[ind1 :ind2,:]
            score_no_transition = score_phase == state
            count = 0
            event_length = []
            for num, value in enumerate(score_no_transition) :
                if value:
                    count +=1
                elif value == False:
                    event_length.append(count)
                    if count == 1:
                        score_no_transition[num-1] = False
                    count = 0

            state_spectrum = spectrum_phase[score_no_transition]
            if state_spectrum.shape[0] <= 1:
                continue

            print(state_spectrum.shape)
            print(state)
            # ax[1,s].plot(freqs,state_spectrum.T, color = 'black', alpha = .01)

            m = state_spectrum.mean(axis=0)
            std = state_spectrum.std(axis=0)
            p90, median, p30 = np.percentile(state_spectrum, q = (90, 50, 30), axis=0)
            # print(m)
            ax[row,s].plot(freqs,median, color = 'darkgreen')
            # ax[row,s].fill_between(freqs, p30, p90, color = 'darkgreen', alpha = .3)
            # ax[s].fill_between(freqs, m-std, m+std, color = 'darkgreen', alpha = .3)
            # ax[1,s].plot(freqs, m+std, color = 'darkgreen', alpha = .6, ls = '--')
            lim = max(median+p90)
            ax[row,s].set_ylim(0, lim*1.25)
            # ax[0,s].set_ylim(0, 30*10**-5)
            # ax[1,s].set_ylim(0, lim*4)
            # ax[1,s].set_ylim(0, 30*10**-5)
            ax[0,s].set_title(state)
    plt.show()


def get_clean_score_one_condition(mouse, condition, state):
    debug = 0
    ds = xr.open_dataset(precompute_dir + '/spectrums/spectrum_scoring_{}.nc'.format(mouse))
    times = ds.coords['times_somno'].values/3600
    score = ds['score'].values
    if condition == 'bl1':
        t1 = 12
        mask_light = (times<t1)
        t1, t2 = 12, 24
        mask_dark = (times>=t1) & (times<t2)
        mask = (times<t2)
    elif condition == 'bl2':
        t1, t2 = 24, 36
        mask_light = (times>=t1) & (times<t2)
        t1, t2 = 36, 48
        mask_dark = (times>=t1) & (times<t2)
        t1, t2 = 24, 48
        mask = (times>=t1) & (times<t2)
    elif condition == 'sd':
        t1, t2 = 54, 60    #####caution 6 first hours is SD 48 to 54
        mask_light = (times>=t1) & (times<t2)
        t1, t2 = 60, 72
        mask_dark = (times>=t1) & (times<t2)
        t1, t2 = 54, 72
        mask = (times>=t1) & (times<t2)
    elif condition == 'r1':
        t1, t2 = 72, 84
        mask_light = (times>=t1) & (times<t2)
        t1 = 84
        mask_dark = (times>=t1)
        t1 = 72
        mask = (times>=t1)
    else :
        print('Condition does not exist !')

    mask_dark = mask_dark[mask]
    mask_light = mask_light[mask]
    score = score[mask]

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
        s = ds.coords['times_somno'].values.size
        for i in range(5):
            ax.axvline(i*s/(3600))
        h = ds.coords['times_somno'].values/3600
        ax.plot(h, np.ones(s))
        ax.plot(h[mask], np.ones(s)[mask])
        ax.plot(h[mask], score == state)

        fig, ax = plt.subplots()
        ax.scatter(h[mask], score_no_artifact)
        ax.plot(h[mask], score_no_transition, color = 'green')
        print(count)
        plt.show()
    if sum(score_no_transition) == 0:
        print('_________ No {} for {} during {} {}________'.format(state, mouse, condition, period))
    Nd= sum(score_no_transition[mask_dark])
    Nl= sum(score_no_transition[mask_light])
    return {'clean_score' :score_no_transition,
            'condition_mask': mask,
            'mask_dark': mask_dark,
            'mask_light': mask_light,
            'Nd': Nd,
            'Nl':Nl }


def get_delta_ref_one_mouse(mouse, state, spectrum_method = 'somno'):
    debug = 1
    ds = xr.open_dataset(precompute_dir + '/spectrums/spectrum_scoring_{}.nc'.format(mouse))
    times = ds.coords['times_somno'].values/3600
    absolute_delta_power = ds['{}_delta_power'.format(spectrum_method)].values
    delta_ref = np.array([])
    for bl in ['bl1', 'bl2']:
        mydict = get_clean_score_one_condition(mouse, bl, state)
        mask = mydict['condition_mask']
        clean_score = mydict['clean_score']
        times_bl = times[mask]-times[mask][0]
        absolute = absolute_delta_power[mask]
        time_of_ref = (times_bl>=8) & (times_bl<=12)
        score_ref = clean_score[time_of_ref]
        delta_ref = np.concatenate((delta_ref, absolute[time_of_ref]))

    ref = np.mean(delta_ref)

    return ref


def get_delta_ratio_one_condition(mouse, condition, state, spectrum_method = 'somno'):
    debug = 1
    print('get delta ratio for {} during {}'.format(mouse, condition))
    ds = xr.open_dataset(precompute_dir + '/spectrums/spectrum_scoring_{}.nc'.format(mouse))
    times = ds.coords['times_somno'].values/3600

    mydict = get_clean_score_one_condition(mouse, condition, state)

    mask = mydict['condition_mask']
    mask_dark = mydict['mask_dark']
    mask_light = mydict['mask_light']
    clean_score = mydict['clean_score']

    times = times[mask]-times[mask][0]
    if condition == 'sd':
        times+=6

    absolute_delta_power = ds['{}_delta_power'.format(spectrum_method)].values
    absolute_delta_power = absolute_delta_power[mask]
    ref = get_delta_ref_one_mouse(mouse, state)
    relative_delta_power = 100*absolute_delta_power/ref

    dark_intervals = 6

    ZT = []
    delta = []

    if condition == 'sd':
        light_intervals = 8
        light_intervals_borders = np.linspace(6,12,light_intervals+1)
        for interval in np.arange(light_intervals):
            ZT1 = light_intervals_borders[interval]
            ZT2 = light_intervals_borders[interval+1]
            sliding_mask = (times>=ZT1) & (times<ZT2)
            sliding_score = clean_score[sliding_mask]
            sliding_delta = relative_delta_power[sliding_mask]
            sliding_delta = sliding_delta[sliding_score]
            ZT.append((ZT2+ZT1)/2)
            delta.append(np.mean(sliding_delta))

    else :
        light_intervals = 12
        light_intervals_borders = np.linspace(0,12,light_intervals+1)
        for interval in np.arange(light_intervals):
            ZT1 = light_intervals_borders[interval]
            ZT2 = light_intervals_borders[interval+1]
            sliding_mask = (times>=ZT1) & (times<ZT2)
            sliding_score = clean_score[sliding_mask]
            sliding_delta = relative_delta_power[sliding_mask]
            sliding_delta = sliding_delta[sliding_score]
            ZT.append((ZT2+ZT1)/2)
            delta.append(np.mean(sliding_delta))

    dark_intervals_borders = np.linspace(12,24,dark_intervals+1)
    for interval in np.arange(dark_intervals):
        ZT1 = dark_intervals_borders[interval]
        ZT2 = dark_intervals_borders[interval+1]
        sliding_mask = (times>=ZT1) & (times<ZT2)
        sliding_score = clean_score[sliding_mask]
        sliding_delta = relative_delta_power[sliding_mask]
        sliding_delta = sliding_delta[sliding_score]
        ZT.append((ZT2+ZT1)/2)
        delta.append(np.mean(sliding_delta))

    if debug :
        fig,ax = plt.subplots()
        ax.scatter(ZT, delta)
        ax.plot(ZT, delta)
        # plt.show()

    return {'ZT':np.array(ZT), 'delta':np.array(delta) }
def get_time_course():
    print('get_time_course')
    mouse = 'B4977'
    time_course = []
    for i, session in enumerate(['bl1', 'bl2', 'sd', 'r1' ]):
        dict = get_delta_ratio_one_condition(mouse, session, 'n')
        ZT = (dict['ZT']+i*24).tolist()
        time_course += ZT

    return time_course


def plot_compare_ratio(state, spectrum_method = 'somno'):
    debug = 1
    control_list = get_mice('Control')
    DCR_list = get_mice('DCR-HCRT')
    groups = {'Control' : control_list, 'DCR-HCRT' : DCR_list}
    time_course = get_time_course()
    indices = []
    df_ratio = pd.DataFrame(index = control_list+DCR_list, columns = ['group'] + time_course )
    for group in groups:
        for mouse in groups[group]:
            print(mouse)
            df_ratio.loc[mouse, 'group'] = group
            delta = np.array([])
            for i, session in enumerate(['bl1', 'bl2', 'sd', 'r1' ]):
                dict = get_delta_ratio_one_condition(mouse, session, state)
                delta = np.append(delta, dict['delta'])
            df_ratio.loc[mouse, time_course] = delta
    for group in groups :
        dirname = excel_dir + '/{}/delta_ratio/'.format(group)
        if not os.path.exists(dirname):
            os.makedirs(dirname)
        name = 'delta_ratio_{}_spectrum_{}_{}.xlsx'.format(state,spectrum_method,group)
        df = df_ratio[df_ratio.group == group]
        df.to_excel(dirname+name)
    if debug :
        fig, ax = plt.subplots()
        fig.suptitle('Delta ratio for ' + state)
        for group in groups:
            subdata = df_ratio[df_ratio.group == group]
            subdata = subdata[time_course].values

            # inf, med, sup = np.nanpercentile(np.array(subdata, dtype = 'float64'), q=[25, 50,75], axis = 0)
            # inf = np.array(inf, dtype = 'float64')
            # sup = np.array(sup, dtype = 'float64')
            # med = np.array(med, dtype = 'float64')
            # ax.plot(time_course, med, label = group)
            # ax.fill_between(time_course, inf, sup, alpha = .2)
            #
            ax.plot(time_course, np.nanmean(subdata, axis = 0), label = group)
            ax.plot(time_course, subdata.T, color = 'black', alpha = .1)

            # ax.scatter(time_course, np.nanmean(subdata, axis = 0))
        ax.legend()
    # plt.show()

def get_event_number_day_light(state):
    for group in ['Control', 'DCR-HCRT']:
        print(group)
        d =  0
        l =  0
        distrib = { 'bl1' : np.zeros(21600),
                    'bl2' : np.zeros(21600),
                    'sd' : np.zeros(16200),
                    'r1' : np.zeros(21600)
                    }
        fig, ax = plt.subplots(nrows = 4, sharex = True)
        for mouse in get_mice(group):
            for condition in ['bl1', 'bl2', 'sd', 'r1']:
                dict = get_clean_score_one_condition(mouse, condition, state)
                d+=dict['Nd']
                l+=dict['Nl']
                clean_score = dict['clean_score']
                dist = distrib[condition]
                dist +=clean_score*1
        plot = distrib['bl1']
        ax[0].plot((4*np.arange(plot.size)/3600), plot)
        plot = distrib['bl2']
        ax[1].plot((4*np.arange(plot.size)/3600), plot)
        plot = distrib['sd']
        ax[2].plot((4*np.arange(plot.size)/3600)+6, plot)
        plot = distrib['r1']
        ax[3].plot((4*np.arange(plot.size)/3600), plot)
        plt.show()
        print(l)
        print(d)
        print(d/l)





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

    group = 'DCR-HCRT'
    # rec = 'b1'

    # mouse = 'B4904'
    # mouse  ='B2762'
    # mouse  ='B2763'
    # mouse  ='B3072'
    # mouse  ='B3140'
    mouse  ='B3513'
    # mouse = 'B4977'
    # mouse = 'B2767'
    # group = 'Control'
    # rec = 'b2'
    # get_event_number_day_light('n')
    # get_delta_ref_one_mouse(mouse, 'n')
    get_delta_ratio_one_condition(mouse, condition = 'bl1', state = 'n')
    get_delta_ratio_one_condition(mouse, condition = 'bl2', state = 'n')
    # get_time_course()
    # plot_compare_ratio('n')
    # plot_compare_ratio('w')
    plt.show()
