from configuration_Malo import *
from select_mice_cata_Malo import get_mice_for_spectrum

local_path = os.path.dirname(os.path.realpath(__file__))
print(local_path)

def compute_all():
    dcr_mice = get_mice_for_spectrum(group = 'DCR-HCRT')
    control_mice = get_mice_for_spectrum(group = 'Control')
    animals_by_group = {'DCR-HCRT' : dcr_mice, 'Control' : control_mice}
    for group in animals_by_group :
        mice = animals_by_group[group]
        for mouse in mice :
            print(mouse)
            for rec in ['b1', 'b2', 'sd', 'r1']:
                try:
                    store_scoring_and_spectrums_one_mouse_one_session(group, mouse, rec)
                except :
                    print('******* ERROR ******', mouse,rec, '*********')

def compute_one_group(group =  'Control'):
    mice = get_mice_for_spectrum(group)
    for mouse in mice :
        print(mouse)
        for rec in ['b1', 'b2', 'sd', 'r1']:
            try:
                store_scoring_and_spectrums_one_mouse_one_session(group, mouse, rec)
            except :
                print('******* ERROR ******', mouse,rec, '*********')


def store_scoring_and_spectrums_one_mouse_one_session(group, mouse,rec):
    print( 'compute {} session'.format(rec))
    #######         Extract scoring         #######
    score = np.loadtxt(data_dir + '/' + group + '/' + mouse + 'DCR' + rec + '.txt', dtype = str)
    print(score.shape)
    ######      Caution ! Do not remove 1, 2, 3. It correspond to invalid EEG
    # for f, n in zip(['1', '2', '3'], ['w', 'n', 'r']) :
    #     score = np.where(score == f, n, score)

    #######         Extract spectrums from txt        #######
    file = open(data_spectrum_dir + '/{}/MTA-{}/{}DCR{}s.txt'.format(group, mouse, mouse, rec))
    real_time_lines = []
    data_lines = []
    real_times = []
    data = []
    delta_power = []


    for i, line in enumerate(file.readlines()) :
        if line == 'Power Spectrum\n':
            real_time_lines.append(i+1)
            data_lines.append(i+2)
        if i in real_time_lines:
            real_times.append(line[:-1])
        if i in data_lines :
            spec = np.array(line.split('\t')[:-1], dtype = float)[4:]    ####Remove .5 Hz     #######
            data.append(spec)
            delta_power.append(np.mean(spec[:15])) #####delta = .75 to 4Hz

    #######         Storing          #######
    freqs = np.arange(4, spec.size+4, 1) * .25
    times = np.arange(len(real_times))*int(4)

    coords = {'real_times' : real_times, 'times' : times, 'freqs' : freqs}
    ds = xr.Dataset(coords = coords)
    data = np.array(data)
    ds['spectrums'] = xr.DataArray(data, dims = ['times', 'freqs'])
    ds['delta_power'] = xr.DataArray(np.array(delta_power), dims = 'times')
    ds['score'] = xr.DataArray(score, dims = 'times')


    dirname = precompute_dir + '/{}/spectrum/'.format(group)
    if not os.path.exists(dirname):
        os.makedirs(dirname)
    ds.to_netcdf(dirname + mouse + 'DCR'+ rec  + '.nc', mode='w')


def plot_anomalies(group, mouse, rec, hours):
    ds = xr.open_dataset(precompute_dir + '/' + group +'/spectrum/' + mouse + 'DCR' + rec + '.nc')
    score = ds['score'].values
    freqs = ds.coords['freqs'].values[:60]
    times = ds.coords['times'].values
    delta_power = ds['delta_power']
    spectrums = ds['spectrums'].values[:,:60]
    mask = (times>hours[0]*3600) & (times<hours[1]*3600)
    fig, ax = plt.subplots()
    ax.plot(delta_power[mask])
    plt.show()
    for i in np.where(mask)[0]:
        print(i/3600*4)
        print(delta_power[i].values)
        fig, ax = plt.subplots()
        ax.plot(freqs, spectrums[i])
        plt.show()
        # exit()



def plot_spectrum_by_state_one_mouse_one_session(group, mouse, rec):
    ds = xr.open_dataset(precompute_dir + '/' + group +'/spectrum/' + mouse + 'DCR' + rec + '.nc')
    score = ds['score'].values
    score_behavior = score.copy()
    for f, n in zip(['1', '2', '3'], ['w', 'n', 'r']) :
        score_behavior = np.where(score_behavior == f, n, score_behavior)

    freqs = ds.coords['freqs'].values[:60]
    times = ds.coords['times']
    t_start = ds.coords['real_times'].values[0]
    t_start = t_start.split('.')[0]
    t_start = t_start.split(':')
    t_start = float(t_start[0]) +float(t_start[1])/60 + float(t_start[2])/3600



    delta_power = ds['delta_power']
    states = ['w', 'n', 'r', 'a']

    hours = np.arange(int(times.size*4/3600))
    delta_power_by_hour = np.zeros(hours.size)
    std_delta_power_by_hour = delta_power_by_hour.copy()
    for h in hours[:-1]:
        # print(h)
        subdata = delta_power[int(3600/4*h):int(3600/4*(h+1))]
        subscore = score[int(3600/4*h):int(3600/4*(h+1))]
        artifact_free = (subscore != '1') & (subscore != '2') & (subscore != '3')
        subdata = subdata[artifact_free]
        m = np.mean(subdata)
        std = np.std(subdata)
        delta_power_by_hour[h] = m
        std_delta_power_by_hour[h] = std

    fig, ax = plt.subplots()
    ax.plot(hours+t_start, delta_power_by_hour, color = 'darkgreen')
    ax.fill_between(hours+t_start, delta_power_by_hour-std_delta_power_by_hour, delta_power_by_hour+std_delta_power_by_hour, color = 'darkgreen', alpha = .5)
    fig.suptitle('Delta power for {}{}'.format(mouse, rec))

    fig, ax = plt.subplots()
    ax.plot(times, delta_power, color = 'black')

    spectrums = ds['spectrums'].values[:,:60]
    print(spectrums.shape)
    fig,ax = plt.subplots(nrows = 2, ncols= 4, sharex = True)

    for s, state in enumerate(states) :
        state_spectrum = spectrums[score == state]
        print(state_spectrum.shape)
        ax[1,s].plot(freqs,state_spectrum.T, color = 'black', alpha = .01)

        m = state_spectrum.mean(axis=0)
        std = state_spectrum.std(axis=0)
        ax[1,s].plot(freqs,m, color = 'darkgreen')
        ax[0,s].plot(freqs,m, color = 'darkgreen')
        ax[0,s].fill_between(freqs, m-std, m+std, color = 'darkgreen', alpha = .3)
        ax[1,s].plot(freqs, m+std, color = 'darkgreen', alpha = .6, ls = '--')
        lim = max(m+std)
        ax[0,s].set_ylim(0, lim*1.25)
        # ax[0,s].set_ylim(0, 30*10**-5)
        # ax[1,s].set_ylim(0, lim*4)
        ax[1,s].set_ylim(0, 30*10**-5)
        ax[0,s].set_title(state)
    plt.show()

def plot_spectrum_by_state_one_mouse(group, mouse):
    ds_b1 = xr.open_dataset(precompute_dir + '/' + group +'/spectrum/' + mouse  + 'DCRb1.nc')
    ds_b2 = xr.open_dataset(precompute_dir + '/' + group +'/spectrum/' + mouse + 'DCRb2.nc')
    ds_sd = xr.open_dataset(precompute_dir + '/' + group +'/spectrum/' + mouse + 'DCRsd.nc')
    ds_r = xr.open_dataset(precompute_dir + '/' + group +'/spectrum/' + mouse + 'DCRr1.nc')
    score = np.concatenate((ds_b1['score'].values,ds_b2['score'].values,ds_sd['score'].values,ds_r['score'].values))
    score_behavior = score.copy()
    for f, n in zip(['1', '2', '3'], ['w', 'n', 'r']) :
        score_behavior = np.where(score_behavior == f, n, score_behavior)

    freqs = ds_b1.coords['freqs'].values[:60]
    t_start = ds_b1.coords['real_times'].values[0]
    t_start = t_start.split('.')[0]
    t_start = t_start.split(':')
    t_start = float(t_start[0]) +float(t_start[1])/60 + float(t_start[2])/3600


    fig, ax = plt.subplots()
    ax.plot(ds_b1.coords['times'].values, ds_b1['delta_power'].values)


    delta_power = np.concatenate((ds_b1['delta_power'].values,ds_b2['delta_power'].values,ds_sd['delta_power'].values,ds_r['delta_power'].values))

    times = np.arange(delta_power.size)*4
    delta_power = np.concatenate((ds_b1['delta_power'].values,ds_b2['delta_power'].values,ds_sd['delta_power'].values,ds_r['delta_power'].values))

    times = np.arange(delta_power.size)*4
    print(delta_power.shape)
    states = ['w', 'n', 'r', 'a']


    fig, ax = plt.subplots()
    ax.plot(times, delta_power)
    fig.suptitle('Delta power every 4s')


    step = .5
    hours = np.arange(int(times.size*(1/step)*4/3600))
    delta_power_by_hour = np.zeros(hours.size)
    std_delta_power_by_hour = delta_power_by_hour.copy()
    for h in hours[:-1]:
        if h <24 :
            delta_ref = np.mean(delta_power[int(8*3600):int(12*3600)])
        elif h >=24 and h<48 :
            delta_ref = np.mean(delta_power[int((8+24)*3600/4):int((12+24)*3600/4)])
        elif h >=48 and h<72 :
            delta_ref = np.mean(delta_power[int((8+48)*3600/4):int((12+48)*3600/4)])
        elif h >=72 :
            delta_ref = np.mean(delta_power[int((8+72)*3600/4):int((12+72)*3600/4)])

        subdata = delta_power[int(step*3600/4*h):int(step*3600/4*(h+1))]/delta_ref
        subscore =score[int(step*3600/4*h):int(step*3600/4*(h+1))]
        mask = (subscore != '1') & (subscore != '2') & (subscore != '3')
        # mask = (subscore != '1') | (subscore != '2') | (subscore != '3')
        print(sum(mask))
        m = np.mean(subdata[mask])
        std = np.std(subdata[mask])
        delta_power_by_hour[h] = m
        std_delta_power_by_hour[h] = std

    fig, ax = plt.subplots()
    ax.plot(hours+t_start, delta_power_by_hour, color = 'darkgreen')
    ax.fill_between(hours+t_start, delta_power_by_hour-std_delta_power_by_hour, delta_power_by_hour+std_delta_power_by_hour, color = 'darkgreen', alpha = .5)
    fig.suptitle('Delta power for {}'.format(mouse))

    plt.show()
    exit()
    spectrums = np.concatenate((ds_b1['spectrums'].values[:,:60],ds_b2['spectrums'].values[:,:60],ds_sd['spectrums'].values[:,:60],ds_r['spectrums'].values[:,:60]))

    print(spectrums.shape)
    fig,ax = plt.subplots(ncols= 4, sharex = True)

    for s, state in enumerate(states) :
        state_spectrum = spectrums[score == state]
        print(state_spectrum.shape)
        m = state_spectrum.mean(axis=0)
        std = state_spectrum.std(axis=0)
        ax[s].plot(freqs,m, color = 'darkgreen')
        ax[s].fill_between(freqs, m-std, m+std, color = 'darkgreen', alpha = .3)
        lim = max(m+std)
        # ax[0,s].set_ylim(0, lim*1.25)
        # ax[1,s].set_ylim(0, lim*4)
        ax[s].set_ylim(0, 30*10**-5)
        ax[s].set_title(state)
        fig.suptitle('Power spectrum by state -- ' + mouse)
    plt.show()


def plot_compare_spectrum_DCR():
    groups = ['DCR-HCRT', 'Control']
    fig_delta, ax_delta = plt.subplots()
    fig_spec, ax_spec =plt.subplots(nrows = 4)
    ds_b1 = xr.open_dataset(precompute_dir + '/Control/spectrum/B4904DCRb1.nc')
    ds_b2 = xr.open_dataset(precompute_dir + '/Control/spectrum/B4904DCRb2.nc')
    ds_sd = xr.open_dataset(precompute_dir + '/Control/spectrum/B4904DCRsd.nc')
    ds_r = xr.open_dataset(precompute_dir + '/Control/spectrum/B4904DCRr1.nc')
    freqs = ds_b2.coords['freqs'].values[:60]
    delta_power = np.concatenate((ds_b1['delta_power'].values,ds_b2['delta_power'].values,ds_sd['delta_power'].values,ds_r['delta_power'].values))
    times = np.arange(delta_power.size)*4
    step = .5
    hours = np.arange(int(times.size*(1/step)*4/3600))

    for group in groups:
        mice = get_mice_for_spectrum(group)
        df_deltas = pd.DataFrame(np.zeros((len(mice), hours.size)),index = mice, columns = hours)
        df_spectrums_a = pd.DataFrame(np.zeros((len(mice), freqs.size)),index = mice, columns = freqs)
        df_spectrums_w = pd.DataFrame(np.zeros((len(mice), freqs.size)),index = mice, columns = freqs)
        df_spectrums_n = pd.DataFrame(np.zeros((len(mice), freqs.size)),index = mice, columns = freqs)
        df_spectrums_r = pd.DataFrame(np.zeros((len(mice), freqs.size)),index = mice, columns = freqs)
        for mouse in mice :
            if mouse == 'B3766' or mouse == 'B4112' or mouse == 'B4113':
                continue
            # try :
            if 1 :
                ds_b1 = xr.open_dataset(precompute_dir + '/' + group +'/spectrum/' + mouse  + 'DCRb1.nc')
                ds_b2 = xr.open_dataset(precompute_dir + '/' + group +'/spectrum/' + mouse + 'DCRb2.nc')
                ds_sd = xr.open_dataset(precompute_dir + '/' + group +'/spectrum/' + mouse + 'DCRsd.nc')
                ds_r = xr.open_dataset(precompute_dir + '/' + group +'/spectrum/' + mouse + 'DCRr1.nc')
                score = np.concatenate((ds_b1['score'].values,ds_b2['score'].values,ds_sd['score'].values,ds_r['score'].values))
                score_behavior = score.copy()
                for f, n in zip(['1', '2', '3'], ['w', 'n', 'r']) :
                    score_behavior = np.where(score_behavior == f, n, score_behavior)

                freqs = ds_b1.coords['freqs'].values[:60]
                t_start = ds_b1.coords['real_times'].values[0]
                t_start = t_start.split('.')[0]
                t_start = t_start.split(':')
                t_start = float(t_start[0]) +float(t_start[1])/60 + float(t_start[2])/3600


                delta_power = np.concatenate((ds_b1['delta_power'].values,ds_b2['delta_power'].values,ds_sd['delta_power'].values,ds_r['delta_power'].values))
                times = np.arange(delta_power.size)*4
                states = ['w', 'n', 'r', 'a']
                step = .5
                hours = np.arange(int(times.size*(1/step)*4/3600))
                delta_power_by_hour = np.zeros(hours.size)
                std_delta_power_by_hour = delta_power_by_hour.copy()
                for h in hours:
                    if h <24 :
                        delta_ref = np.mean(delta_power[int(8*3600):int(12*3600)])
                    elif h >=24 and h<48 :
                        delta_ref = np.mean(delta_power[int((8+24)*3600/4):int((12+24)*3600/4)])
                    elif h >=48 and h<72 :
                        delta_ref = np.mean(delta_power[int((8+48)*3600/4):int((12+48)*3600/4)])
                    elif h >=72 :
                        delta_ref = np.mean(delta_power[int((8+72)*3600/4):int((12+72)*3600/4)])

                    subdata = delta_power[int(step*3600/4*h):int(step*3600/4*(h+1))]/delta_ref
                    subscore =score[int(step*3600/4*h):int(step*3600/4*(h+1))]
                    mask = (subscore != '1') & (subscore != '2') & (subscore != '3')
                    m = np.mean(subdata[mask])
                    std = np.std(subdata[mask])
                    delta_power_by_hour[h] = m
                    std_delta_power_by_hour[h] = std

                df_deltas.loc[mouse] = delta_power_by_hour

                spectrums = np.concatenate((ds_b1['spectrums'].values[:,:60],ds_b2['spectrums'].values[:,:60],ds_sd['spectrums'].values[:,:60],ds_r['spectrums'].values[:,:60]))

                state_spectrum = spectrums[score == 'a']
                m = state_spectrum.mean(axis=0)
                df_spectrums_a.loc[mouse] = m
                state_spectrum = spectrums[score == 'w']
                m = state_spectrum.mean(axis=0)
                df_spectrums_w.loc[mouse] = m
                state_spectrum = spectrums[score == 'n']
                m = state_spectrum.mean(axis=0)
                df_spectrums_n.loc[mouse] = m
                state_spectrum = spectrums[score == 'r']
                m = state_spectrum.mean(axis=0)
                df_spectrums_r.loc[mouse] = m
        df_deltas.to_excel(precompute_dir + '/{}/Spectrum/delta_power.xlsx'.format(group))
        df_spectrums_a.to_excel(precompute_dir + '/{}/Spectrum/spectrums_a.xlsx'.format(group))
        df_spectrums_w.to_excel(precompute_dir + '/{}/Spectrum/spectrums_w.xlsx'.format(group))
        df_spectrums_n.to_excel(precompute_dir + '/{}/Spectrum/spectrums_n.xlsx'.format(group))
        df_spectrums_r.to_excel(precompute_dir + '/{}/Spectrum/spectrums_r.xlsx'.format(group))

        ax_delta.plot(hours+t_start, df_deltas.mean(axis = 0), label = group)
        ax_spec[0].plot(freqs, df_spectrums_w.mean(axis = 0), label = group)
        ax_spec[1].plot(freqs, df_spectrums_n.mean(axis = 0), label = group)
        ax_spec[2].plot(freqs, df_spectrums_r.mean(axis = 0), label = group)
        ax_spec[3].plot(freqs, df_spectrums_a.mean(axis = 0), label = group)

        ax_spec[0].set_ylabel('wake')
        ax_spec[1].set_ylabel('nrem')
        ax_spec[2].set_ylabel('rem')
        ax_spec[3].set_ylabel('cata')

    ax_spec[0].legend()
    ax_spec[1].legend()
    ax_spec[2].legend()
    ax_spec[3].legend()

    ax_delta.legend()
    plt.show()
            # except :
            #     print('******** ERROR during precompute collecting ******* Mouse : ', mouse)
    #
    #
    #
    # fig,ax = plt.subplots(ncols= 4, sharex = True)
    # for s, state in enumerate(states) :
    #     state_spectrum = spectrums[score == state]
    #     print(state_spectrum.shape)
    #     m = state_spectrum.mean(axis=0)
    #     std = state_spectrum.std(axis=0)
    #     ax[s].plot(freqs,m, color = 'darkgreen')
    #     ax[s].fill_between(freqs, m-std, m+std, color = 'darkgreen', alpha = .3)
    #     lim = max(m+std)
    #     # ax[0,s].set_ylim(0, lim*1.25)
    #     # ax[1,s].set_ylim(0, lim*4)
    #     ax[s].set_ylim(0, 30*10**-5)
    #     ax[s].set_title(state)
    #     fig.suptitle('Power spectrum by state -- ' + mouse)
    # plt.show()

# def plot_compare_spectrum(KO = 'Mir_137'):
    # 1- get mean spectrums by mouse and  state
        # control mice  = {
                            # mice :{
                                    # 'w': 0,
                                    # 'r': 0,
                                    # 'n': 0,
                                    # 'a': 0
                                    # }
                        # }
        # KO_mice...
    # for mouse in controle mice :
    #     ds = xr.open_dataset(local_path + '/precompute/' + file_name + '.nc')
    #     score = ds['score'].values
    #     freqs = ds.coords['freqs'].values[:60]
    #     spectrums = ds['spectrums'].values[:,:60]
        # for state in states :
        # control_mice['mouse'] = np.mean(spectrums[score == state], axis = 0)
    # 2 plot mean by group ie mean of means and DIFF
    # fig, ax = plt.suplots(nrows =2, ncols = 4)
    # KO,control = [],[]
    # for s, state in enumerate(states) :
    #     for mouse in mice
    #         control.append(control_mice[mouse][state])
    #
    #     ax.plot[0, s].plot(freqs, np.mean(KO,  axis = 0), color = 'red')
    #     ax.plot[0, s].fill_between(freqs, np.mean(KO,  axis = 0) + np.std(KO,   axis = 0), color = 'red', alpha = .5)
    #     ax.plot[0, s].plot(freqs, np.mean(control, axis = 0), color = 'black')
    #     ax.plot[0, s].plot(freqs, np.mean(control, axis = 0), color = 'black')
    #     ax.plot[0, s].fill_between(freqs, + np.std(control, axis = 0), color = 'black', alpha = .5)
    #     diff = np.mean(KO,  axis = 0) - np.mean(control, axis = 0)
    #     ax.plot[1, s].plot(freqs, diff, color = 'black')
    #     ax.plot[1, s].axhline(0, color = 'red')
    #     t = .5
    #     signif_mask = (diff < -t) | (diff > t)
    #     ax.plot(freqs[signif_mask], diff[signif_mask], color = 'red')



if __name__ == '__main__':
    mouse = 'B3512'
    group = 'DCR-HCRT'
    rec = 'b1'
    # mouse = 'B4904'
    # group = 'Control'
    # rec = 'b2'
    # compute_all()
    hours = [11424/3600,12000/3600]
    # compute_one_group(group)
    # store_scoring_and_spectrums_one_mouse_one_session(group, mouse,rec)
    # plot_spectrum_by_state_one_mouse_one_session(group, mouse, rec)
    plot_anomalies(group, mouse, rec, hours)

    # plot_spectrum_by_state_one_mouse(group, mouse)
    # plot_compare_spectrum_DCR()
    # plot_compare_spectrum(KO = 'Mir_137')
    # plot_compare_spectrum(KO = 'Mir_665')
    # plot_compare_spectrum(KO = 'Mir_128')
