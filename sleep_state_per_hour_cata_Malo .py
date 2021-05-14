from configuration import *
import os
os.system('pwd')
os.system('which python')
os.system('python --version')
import numpy as np
import pandas as pd
import xarray as xr
import matplotlib.pyplot as plt
import seaborn as sns
from select_mice_cata_Malo import get_mice,get_mouse_info

#
def compute_all(exclude_mice = False):
    for miRNA in ['128', '137', '665']:
        for genotype in ['test', 'control']:
            print(miRNA, genotype)
            precompute_sleep_state_by_epoch(miRNA, genotype, exclude_mice)
            sleep_state_statistics(miRNA, genotype, exclude_mice)
            sleep_bouts(miRNA, genotype, exclude_mice)
            fourfirsthours_sleep_bouts(miRNA, genotype,exclude_mice=exclude_mice)
            inter_REM_interval(miRNA, genotype, exclude_mice)
            tdw_wake_ratio(miRNA, genotype, exclude_mice)
            REM_sleep_latency_by_hour(miRNA, genotype, exclude_mice)
            for mouse in get_mice(miRNA, genotype):
                if exclude_mice and mouse in RT_PCR_execption:
                    continue
                else:
                    REM_sleep_latency_one_mouse(mouse, exclude_mice)


def precompute_sleep_state_by_epoch(miRNA, genotype, exclude_mice=False):
    mice = get_mice(miRNA, genotype, exclude_mice)

    days = [ 'b1', 'b2', 'sd', 'r1' ]

    all_mice_all_days = {}
    numbers = ['1', '2', '3', '4' ,'5' ,'6']
    letters = ['w', 'n', 'r', 'w', 'n', 'r']
    dirname = precompute_dir + '/tdw_score/'

    for mouse in mice :
        ds = xr.open_dataset(dirname + '/tdw_score_{}.nc'.format(mouse))

        # ds = xr.open_dataset(precompute_dir + '/spectrums/spectrum_scoring_{}.nc'.format(mouse))
        one_mouse = ds['new_score'].values
        for f, n in zip(['1', '2', '3'], ['w', 'n', 'r']) :
            one_mouse = np.where(one_mouse == f, n, one_mouse)
        date_ref = pd.read_excel(work_dir + 'datetime_reference_miRNA.xlsx', index_col = 0)
        if date_ref.at['MTA-{}'.format(mouse), 'short_end'] == 'yes':
            end = np.ones(10800)
            end[:]=np.nan
            one_mouse = np.append(one_mouse, end)
        all_mice_all_days[mouse] = one_mouse
    all_mice_all_days = pd.DataFrame.from_dict(all_mice_all_days)
    print(all_mice_all_days)
    # all_mice_all_days = all_mice_all_days.T

    print('find w, n, r')
    all_mice_w_epoch = (all_mice_all_days == 'w') | (all_mice_all_days == 't')
    all_mice_n_epoch = all_mice_all_days == 'n'
    all_mice_r_epoch = all_mice_all_days == 'r'
    all_mice_tdw_epoch = all_mice_all_days == 't'


    print('create Dataset')
    if exclude_mice :
        path = precompute_dir + '/RT_PCR_excluded_mice/sleep_by_epoch/'
    else :
        path = precompute_dir + '/sleep_by_epoch/'
    if not os.path.exists(path):
        os.makedirs(path)
    # ds.to_netcdf(path + 'summed_sleep_by_epoch.nc', mode='w')
    epochs = np.arange(all_mice_all_days.shape[0])
    # print(epochs)
    # print(type(mice))
    coords = {'mice':mice, 'epochs':epochs }
    ds = xr.Dataset(coords = coords)
    # print(all_mice_w_epoch)

    ds['wake_all_mice'] = xr.DataArray(all_mice_w_epoch*1, dims = ['epochs','mice'])
    ds['nrem_all_mice'] = xr.DataArray(all_mice_n_epoch*1, dims = ['epochs','mice'])
    ds['rem_all_mice'] = xr.DataArray(all_mice_r_epoch*1, dims = ['epochs','mice'])
    ds['tdw_all_mice'] = xr.DataArray(all_mice_tdw_epoch*1, dims = ['epochs','mice'])
    # print(ds)
    # exit()
    print('save')
    if exclude_mice :
        ds.to_netcdf(path + '{}_{}_exclude_mice_sleep_by_epoch.nc'.format(miRNA, genotype), mode='w')
    else :
        ds.to_netcdf(path + '{}_{}_sleep_by_epoch.nc'.format(miRNA, genotype), mode='w')

    #
def sleep_state_statistics( miRNA, genotype,exclude_mice=False):
    if exclude_mice :
        path = precompute_dir+'/RT_PCR_excluded_mice/sleep_by_epoch/{}_{}_exclude_mice_sleep_by_epoch.nc'.format(miRNA, genotype)
    else :
        path = precompute_dir+'/sleep_by_epoch/{}_{}_sleep_by_epoch.nc'.format(miRNA, genotype)
    ds = xr.open_dataset(path)
    level = [0, 1]
    # print(ds)
    all_mice_w_epoch = ds['wake_all_mice'].to_pandas()
    all_mice_n_epoch = ds['nrem_all_mice'].to_pandas()
    all_mice_r_epoch = ds['rem_all_mice'].to_pandas()
    all_mice_t_epoch = ds['tdw_all_mice'].to_pandas()

    number_of_epochs = all_mice_w_epoch.shape[0]
    epoch_duration = 4
    duration_in_hours = int(number_of_epochs*epoch_duration/3600)
    time = int(1) # window to look at in HOUR. CAUTION must be a multiple of 100
    hours = np.arange(0, int(duration_in_hours),time)
    window = int(time * 3600 / epoch_duration)

    fake = np.zeros(number_of_epochs)

    all_mice_w_by_time = []
    all_mice_r_by_time = []
    all_mice_n_by_time = []
    all_mice_t_by_time = []
    for hour in hours :
        i1, i2 = int(hour*window), int((hour+1)*window)
        # tarace[hour] = all_mice_w_epoch[i1:i2].sum(axis =0)
        all_mice_w_by_time.append(all_mice_w_epoch[i1:i2].sum(axis =0)*4/60)
        all_mice_r_by_time.append(all_mice_r_epoch[i1:i2].sum(axis =0)*4/60)
        all_mice_n_by_time.append(all_mice_n_epoch[i1:i2].sum(axis =0)*4/60)
        all_mice_t_by_time.append(all_mice_t_epoch[i1:i2].sum(axis =0)*4/60)

    all_mice_w_by_time = pd.concat(all_mice_w_by_time, axis = 1)
    all_mice_r_by_time = pd.concat(all_mice_r_by_time, axis = 1)
    all_mice_n_by_time = pd.concat(all_mice_n_by_time, axis = 1)
    all_mice_t_by_time = pd.concat(all_mice_t_by_time, axis = 1)

    if exclude_mice:
        dirname = excel_dir+  '/RT_PCR_excluded_mice/time_dynamics/{}/{}/'.format(miRNA, genotype)
    else :
        dirname = excel_dir+  '/time_dynamics/{}/{}/'.format(miRNA, genotype)
    if not os.path.exists(dirname):
        os.makedirs(dirname)
    all_mice_w_by_time.to_excel(dirname+   '/wake_event_by_mouse_by_hour.xlsx')
    all_mice_n_by_time.to_excel(dirname+   '/NREM_event_by_mouse_by_hour.xlsx')
    all_mice_r_by_time.to_excel(dirname+   '/REM_event_by_mouse_by_hour.xlsx')
    all_mice_t_by_time.to_excel(dirname+   '/tdw_event_by_mouse_by_hour.xlsx')



#
# def plot_sleep_state_accross_time(miRNA, genotype):
#     # data_dir = 'C:/Users/maxime.juventin/Desktop/scripts_ML/data/'
#     path = precompute_dir+'/sleep_by_epoch/{}_{}_sleep_by_epoch.nc'.format(miRNA, genotype)
#     ds = xr.open_dataset(path)
#
#     print(ds)
#     number_of_epochs = ds['wake_all_mice'].shape[0]
#     epoch_duration = 4
#     duration_in_hours = int(number_of_epochs*epoch_duration/3600)
#     wake = ds['wake_all_mice'].values
#     nrem = ds['nrem_all_mice'].values
#     rem = ds['rem_all_mice'].values
#
#     time = int(1) # window to look at in HOUR.
#     window = int(time * 3600 / epoch_duration)
#     wake_by_time = np.zeros(int(duration_in_hours/time))
#     rem_by_time = np.zeros(int(duration_in_hours/time))
#     cata_by_time = np.zeros(int(duration_in_hours/time))
#     nrem_by_time = np.zeros(int(duration_in_hours/time))
#     for h, i in enumerate(np.arange(0, number_of_epochs, window)):
#         wake_by_time[h] = wake[i : i+window].sum()
#         rem_by_time[h] = rem[i : i+window].sum()
#         nrem_by_time[h] = nrem[i : i+window].sum()
#     print(wake_by_time.shape)
#     fig, ax = plt.subplots(nrows=3)
#     times = np.arange(0, duration_in_hours,time)
#     ax[0].plot(times, wake_by_time*epoch_duration/60, color = 'black', label = 'wake')
#     ax[1].plot(times, rem_by_time*epoch_duration/60, color = 'red', label='rem')
#     ax[2].plot(times, nrem_by_time*epoch_duration/60, color = 'blue', label = 'nrem')
#     for i,j in zip([0,28-time,76-time], [4-time,52-time,100-time]):
#         ax[0].axvspan(i,j , color = 'black', alpha = .3)
#         ax[1].axvspan(i,j , color = 'black', alpha = .3)
#         ax[2].axvspan(i,j , color = 'black', alpha = .3)
#     ax[0].set_title('Sleep state amount for '+ genotype + ' , sum per ' + str(time) + ' hours')
#     ax[0].legend()
#     ax[1].legend()
#     ax[2].legend()
#
#     plt.show()
#
# #
# def plot_compare_sleep_state_accross_time(control ='Control', test ='DCR_HCRT'):
#     # data_dir = 'C:/Users/maxime.juventin/Desktop/scripts_ML/data/'
#     path_control = precompute_dir +  '/' +control+'/sleep_by_epoch.nc'
#     path_cre = precompute_dir + '/'+ test+'/sleep_by_epoch.nc'
#     ds_control = xr.open_dataset(path_control)
#     ds_cre = xr.open_dataset(path_cre)
#
#     number_of_epochs = ds_cre['wake_all_mice'].shape[0]
#     epoch_duration = 4
#     duration_in_hours = int(number_of_epochs*epoch_duration/3600)
#     wake_control = ds_control['wake_all_mice'].to_pandas()
#     nrem_control = ds_control['nrem_all_mice'].to_pandas()
#     rem_control = ds_control['rem_all_mice'].to_pandas()
#     cata_control = ds_control['cata_all_mice'].to_pandas()
#
#     wake_cre = ds_cre['wake_all_mice'].to_pandas()
#     nrem_cre = ds_cre['nrem_all_mice'].to_pandas()
#     rem_cre = ds_cre['rem_all_mice'].to_pandas()
#     cata_cre = ds_cre['cata_all_mice'].to_pandas()
#
#     time = int(1) # window to look at in HOUR. CAUTION must be a multiple of 100
#     window = int(time * 3600 / epoch_duration)
#     wake_by_time_control = []
#     rem_by_time_control = []
#     nrem_by_time_control = []
#     cata_by_time_control = []
#
#     wake_by_time_cre = []
#     rem_by_time_cre = []
#     nrem_by_time_cre = []
#     cata_by_time_cre = []
#
#     for h, i in enumerate(np.arange(0, number_of_epochs, window)):
#
#         wake_by_time_control.append(wake_control[i : i+window].sum())
#         rem_by_time_control.append(rem_control[i : i+window].sum())
#         nrem_by_time_control.append(nrem_control[i : i+window].sum())
#         cata_by_time_control.append(cata_control[i : i+window].sum())
#
#         wake_by_time_cre.append(wake_cre[i : i+window].sum())
#         rem_by_time_cre.append(rem_cre[i : i+window].sum())
#         nrem_by_time_cre.append(nrem_cre[i : i+window].sum())
#         cata_by_time_cre.append(cata_cre[i : i+window].sum())
#     wake_by_time_control = pd.concat(wake_by_time_control, axis = 1)
#     cata_by_time_control = pd.concat(cata_by_time_control, axis = 1)
#     rem_by_time_control = pd.concat(rem_by_time_control, axis = 1)
#     nrem_by_time_control = pd.concat(nrem_by_time_control, axis = 1)
#     wake_by_time_cre = pd.concat(wake_by_time_cre, axis = 1)
#     rem_by_time_cre = pd.concat(rem_by_time_cre, axis = 1)
#     nrem_by_time_cre = pd.concat(nrem_by_time_cre, axis = 1)
#     cata_by_time_cre = pd.concat(cata_by_time_cre, axis = 1)
#
#     # print(wake_by_time.shape)
#     fig, ax = plt.subplots(nrows=3)
#     times = np.arange(0, duration_in_hours)
#     plot  = np.arange(3)
#     controls = [wake_by_time_control, rem_by_time_control, nrem_by_time_control, cata_by_time_control]
#     cres = [wake_by_time_cre, rem_by_time_cre, nrem_by_time_cre, cata_by_time_cre]
#     # colors = [ 'black', 'red', 'blue']
#     # colors = [ 'black', 'green', 'black', 'green', 'black', 'green']
#     labels = ['control', 'control' ,'control', 'control', 'DCR', 'DCR', 'DCR' 'DCR']
#
#     for i, data_control, data_cre, label in zip(plot, controls, cres,labels) :
#         m = data_control.mean(axis=0)*epoch_duration/60
#         s =data_control.std(axis=0)*epoch_duration/60
#         ax[i].plot(times, m, color = 'black', label = label)
#         ax[i].fill_between(times, m-s, m+s , color = 'black', alpha = .2,)
#         m = data_cre.mean(axis=0)*epoch_duration/60
#         s =data_cre.std(axis=0)*epoch_duration/60
#         ax[i].plot(times, m, color = 'green', label = label)
#         ax[i].fill_between(times, m-s, m+s , color = 'green', alpha = .2)
#
#
#     nights = { 'dark0' : [0, 4],
#                  'dark1' : [16, 28],
#                  'dark2' : [40, 52],
#                  'dark3':[64, 76],
#                  'dark4':[88, 100]}
#
#     for night in nights :
#         h1, h2 = nights[night][0], nights[night][1]
#         ax[0].axvspan(h1, h2 , color = 'black', alpha = .3)
#         ax[1].axvspan(h1, h2 , color = 'black', alpha = .3)
#         ax[2].axvspan(h1, h2 , color = 'black', alpha = .3)
#         ax[3].axvspan(h1, h2 , color = 'black', alpha = .3)
#     height0_1, height0_2 =ax[0].get_ylim()[1]*.95, ax[0].get_ylim()[1]*1.05
#     height1_1, height1_2 = ax[1].get_ylim()[1]*.95, ax[1].get_ylim()[1]*1.05
#     height2_1, height2_2 = ax[2].get_ylim()[1]*.95, ax[2].get_ylim()[1]*1.05
#     height3_1, height3_2 = ax[3].get_ylim()[1]*.95, ax[3].get_ylim()[1]*1.05
#     for x1, x2 in zip([0,28,76], [4,52,100]):
#         lim = 110
#         x1 += 5
#         x2 +=5
#         ax[0].axhspan(height0_1, height0_2, xmin =x1/lim, xmax = x2/lim, fill = False, edgecolor = 'black', hatch = '////')
#         ax[1].axhspan(height1_1, height1_2, xmin =x1/lim, xmax = x2/lim, fill = False, edgecolor = 'black', hatch = '////')
#         ax[2].axhspan(height2_1, height2_2, xmin =x1/lim, xmax = x2/lim, fill = False, edgecolor = 'black', hatch = '////')
#         ax[2].axhspan(height3_1, height3_2, xmin =x1/lim, xmax = x2/lim, fill = False, edgecolor = 'black', hatch = '////')
#     x1 = 52 + 5
#     x2 = 58 +5
#     ax[0].axhspan(height0_1, height0_2, xmin =x1/lim, xmax = x2/lim, color = 'red')
#     ax[1].axhspan(height1_1, height1_2, xmin =x1/lim, xmax = x2/lim, color = 'red')
#     ax[2].axhspan(height2_1, height2_2, xmin =x1/lim, xmax = x2/lim, color = 'red')
#     ax[3].axhspan(height3_1, height3_2, xmin =x1/lim, xmax = x2/lim, color = 'red')
#         # ax[0].axhspan(i,j , color = 'black', alpha = .3)
#         # ax[1].axhspan(i,j , color = 'black', alpha = .3)
#         # ax[2].axhspan(i,j , color = 'black', alpha = .3)
#     ax[0].set_title('Sleep state amount DCR vs control for , sum per ' + str(time) + ' hours')
#     ax[0].legend()
#     ax[1].legend()
#     ax[2].legend()
#     ax[3].legend()
#     ax[0].set_ylabel('wake')
#     ax[1].set_ylabel('rem')
#     ax[2].set_ylabel('nrem')
#     ax[3].set_ylabel('cata')
#
#     plt.show()

def sleep_bouts(miRNA, genotype,exclude_mice=False):
    # B0 begins at 4h, B1 etc begins at 8

    if exclude_mice :
        path = precompute_dir+'/RT_PCR_excluded_mice/sleep_by_epoch/{}_{}_exclude_mice_sleep_by_epoch.nc'.format(miRNA, genotype)
    else :
        path = precompute_dir+'/sleep_by_epoch/{}_{}_sleep_by_epoch.nc'.format(miRNA, genotype)
    print(precompute_dir)
    print(path)
    ds = xr.open_dataset(path)


    all_mice_w_epoch = ds['wake_all_mice'].to_pandas()
    all_mice_n_epoch = ds['nrem_all_mice'].to_pandas()
    all_mice_t_epoch = ds['tdw_all_mice'].to_pandas()
    all_mice_r_epoch = ds['rem_all_mice'].to_pandas()    #.unstack().T*1
    data_by_state = {'wake':all_mice_w_epoch,
                    'tdw' : all_mice_t_epoch,
                    'nrem' : all_mice_n_epoch,
                    'rem' :all_mice_r_epoch}
    number_of_epochs = ds['wake_all_mice'].shape[0]
    # print(number_of_epochs)
    epoch_duration = 4
    duration_in_hours = int(number_of_epochs*epoch_duration/3600)
    day_cycles = {'dark0': [0,4],
                 'light1' : [4, 16],
                 'dark1' : [16, 28],
                 'light2' : [28, 40],
                 'dark2' : [40, 52],
                 'sd' : [52, 58],
                 'light3' : [58,64],
                 'dark3':[64, 76],
                 'light4':[76,88 ],
                 'dark4':[88, 100]}
    # time = int(1) # window to look at in HOUR. CAUTION must be a multiple of 100
    # window = int(time * 3600 / epoch_duration)
    # hours_in_epochs = hours*3600/epoch_duration
    # print(hours_in_epochs)
    mice = ds.coords['mice'].values

    for s, state in enumerate(data_by_state):
        data = data_by_state[state]
        for c, cycle in enumerate(day_cycles):
            cycle_hours = day_cycles[cycle]
            # print(cycle_hours)
            i1, i2 = int(cycle_hours[0]), int(cycle_hours[1])
            i1 = int(i1 *3600/epoch_duration)
            i2 = int(i2 *3600/epoch_duration)
            selected_data = data[i1: i2]
            # print(diff.where(diff == -1))
            all_mice_count = {}
            all_mice_percent_frag = {}
            all_mice_mean_bout_duration = {}
            all_mice_bout_amount = {}
            all_mice_state_time = {}
            # print('je suis la')
            for mouse in mice :
                one_and_zeros_one_mouse = selected_data[mouse].values
                counter = 0
                bouts = []
                state_time = np.sum(one_and_zeros_one_mouse)*4/60  #min
                all_mice_state_time[mouse] = state_time

                for i in one_and_zeros_one_mouse:
                    if i ==0:
                        bouts.append(counter)
                        counter = 0
                    if i ==1:
                        counter +=1
                #### True duration last bout unknown
                bouts.append(counter)
                bouts=np.array(bouts)
                mask = bouts!=0
                bouts = bouts[mask]
                all_mice_bout_amount[mouse] = [len(bouts), len(bouts)/12]

                bouts = bouts*4
                mean_bouts = np.mean(bouts)/60 #min
                all_mice_mean_bout_duration[mouse] = mean_bouts

                if bouts.size >0 :
                    if np.max(bouts) > 2**12:
                        bins_sleep = np.append(2**np.arange(2,12), np.max(bouts))
                    else :
                        bins_sleep = 2**np.arange(2,13)
                    bins_sleep= np.append(np.array([0]), bins_sleep)

                    #bins_sleep += [max(bouts)]
                    # count, bins= np.histogram(bouts, bins = bins_sleep)
                    state_percentage_frag = []
                    count = []
                    # print(bins_sleep)
                    for b in np.arange(1,bins_sleep.size):
                        # print(bins_sleep[b])
                        # print(bins_sleep[b+1])
                        # print()
                        mask = (bouts>bins_sleep[b-1]) & (bouts<=bins_sleep[b])
                        # print(np.unique(bouts[mask]))
                        # print(bins_sleep[b])
                        state_percentage_frag.append(np.sum(bouts[mask])/(60*state_time))
                        count.append(np.sum(mask))
                    state_percentage_frag = np.array(state_percentage_frag)
                    count = np.array(count)
                    # exit()

                    all_mice_count[mouse] = count
                    all_mice_percent_frag[mouse] = state_percentage_frag
                    # fig, ax = plt.subplots()
                    # ax.plot(count)
                    # axb=ax.twinx()
                    # axb.plot(state_percentage_frag, color ='k')
                    # plt.show()
                    # exit()

                else :
                    bins_sleep = 2**np.arange(2,13)
                    bins_sleep= np.append(np.array([0]), bins_sleep)
                    count = np.zeros(bins_sleep.size-1)
                    state_percentage_frag = np.zeros(bins_sleep.size-1)
                    # count_vs_total_time = np.zeros(bins_sleep.size-1)

                    all_mice_count[mouse] = count
                    all_mice_percent_frag[mouse] = state_percentage_frag

            bins_sleep = 2**np.arange(2,13)
            df_percent_frag = pd.DataFrame.from_dict(all_mice_percent_frag, orient = 'index', columns = bins_sleep )
            if exclude_mice :
                output_path = excel_dir + '/RT_PCR_excluded_mice/sleep_percent_fragmentation/{}/{}/{}/'.format(miRNA, genotype, state)
            else :
                output_path = excel_dir + '/sleep_percent_fragmentation/{}/{}/{}/'.format(miRNA, genotype, state)
            if not os.path.exists(output_path):
                os.makedirs(output_path)
            df_percent_frag.to_excel(output_path +'{}_{}_{}_{}.xlsx'.format(miRNA, genotype, state,cycle))

            df = pd.DataFrame.from_dict(all_mice_count, orient = 'index', columns = bins_sleep )
            if exclude_mice :
                output_path = excel_dir + '/RT_PCR_excluded_mice/sleep_fragmentation/{}/{}/{}/'.format(miRNA, genotype, state)
            else :
                output_path = excel_dir + '/sleep_fragmentation/{}/{}/{}/'.format(miRNA, genotype, state)
            if not os.path.exists(output_path):
                os.makedirs(output_path)
            df.to_excel(output_path +'{}_{}_{}_{}.xlsx'.format(miRNA, genotype, state,cycle))

            df_duration = pd.DataFrame.from_dict(all_mice_mean_bout_duration, orient = 'index', columns = ['mean_duration'] )
            if exclude_mice :
                output_path = excel_dir + '/RT_PCR_excluded_mice/state_duration/{}/{}/{}/'.format(miRNA, genotype, state)
            else:
                output_path = excel_dir + '/state_duration/{}/{}/{}/'.format(miRNA, genotype, state)
            if not os.path.exists(output_path):
                os.makedirs(output_path)
            df_duration.to_excel(output_path +'{}_{}_{}_{}.xlsx'.format(miRNA, genotype, state,cycle))

            df_time = pd.DataFrame.from_dict(all_mice_state_time, orient = 'index', columns = ['time'] )
            if exclude_mice:
                output_path = excel_dir + '/RT_PCR_excluded_mice/state_time/{}/{}/{}/'.format(miRNA, genotype, state)
            else :
                output_path = excel_dir + '/state_time/{}/{}/{}/'.format(miRNA, genotype, state)
            if not os.path.exists(output_path):
                os.makedirs(output_path)
            df_time.to_excel(output_path +'{}_{}_{}_{}.xlsx'.format(miRNA, genotype, state,cycle))

            df_period_amount = pd.DataFrame.from_dict(all_mice_bout_amount, orient = 'index', columns = ['period', 'freq'] )
            if exclude_mice :
                output_path = excel_dir + '/RT_PCR_excluded_mice/period_amount/{}/{}/{}/'.format(miRNA, genotype, state)
            else :
                output_path = excel_dir + '/period_amount/{}/{}/{}/'.format(miRNA, genotype, state)
            if not os.path.exists(output_path):
                os.makedirs(output_path)
            df_period_amount.to_excel(output_path +'{}_{}_{}_{}.xlsx'.format(miRNA, genotype, state,cycle))

def fourfirsthours_sleep_bouts(miRNA, genotype,exclude_mice=False):
    # B0 begins at 4h, B1 etc begins at 8

    if exclude_mice :
        path = precompute_dir+'/RT_PCR_excluded_mice/sleep_by_epoch/{}_{}_exclude_mice_sleep_by_epoch.nc'.format(miRNA, genotype)
    else :
        path = precompute_dir+'/sleep_by_epoch/{}_{}_sleep_by_epoch.nc'.format(miRNA, genotype)
    print(precompute_dir)
    print(path)
    ds = xr.open_dataset(path)


    all_mice_w_epoch = ds['wake_all_mice'].to_pandas()
    all_mice_n_epoch = ds['nrem_all_mice'].to_pandas()
    all_mice_t_epoch = ds['tdw_all_mice'].to_pandas()
    all_mice_r_epoch = ds['rem_all_mice'].to_pandas()    #.unstack().T*1
    data_by_state = {'wake':all_mice_w_epoch,
                    'tdw' : all_mice_t_epoch,
                    'nrem' : all_mice_n_epoch,
                    'rem' :all_mice_r_epoch}
    number_of_epochs = ds['wake_all_mice'].shape[0]
    # print(number_of_epochs)
    epoch_duration = 4
    duration_in_hours = int(number_of_epochs*epoch_duration/3600)
    day_cycles = {
                 'dark1' : [16, 20],
                 'dark2' : [40, 44],
                 'post_sd' : [58, 62],
                 'dark3':[64, 68],
                 'dark4':[88, 92]}
    # time = int(1) # window to look at in HOUR. CAUTION must be a multiple of 100
    # window = int(time * 3600 / epoch_duration)
    # hours_in_epochs = hours*3600/epoch_duration
    # print(hours_in_epochs)
    mice = ds.coords['mice'].values

    for s, state in enumerate(data_by_state):
        data = data_by_state[state]
        for c, cycle in enumerate(day_cycles):
            cycle_hours = day_cycles[cycle]
            # print(cycle_hours)
            i1, i2 = int(cycle_hours[0]), int(cycle_hours[1])
            i1 = int(i1 *3600/epoch_duration)
            i2 = int(i2 *3600/epoch_duration)
            selected_data = data[i1: i2]
            # print(diff.where(diff == -1))
            all_mice_count = {}
            all_mice_percent_frag = {}
            all_mice_mean_bout_duration = {}
            all_mice_bout_amount = {}
            all_mice_state_time = {}
            # print('je suis la')
            for mouse in mice :
                one_and_zeros_one_mouse = selected_data[mouse].values
                counter = 0
                bouts = []
                state_time = np.sum(one_and_zeros_one_mouse)*4/60  #min
                all_mice_state_time[mouse] = state_time

                # first = 0
                for i in one_and_zeros_one_mouse:
                    if i ==0:
                        bouts.append(counter)
                        counter = 0
                    if i ==1:
                        # if first ==0:
                        #     remove_first = True
                        counter +=1
                    # first+=1
                ### True duration last bout unknown
                bouts.append(counter)

                bouts=np.array(bouts)
                mask = bouts!=0
                bouts = bouts[mask]
                all_mice_bout_amount[mouse] = [len(bouts), len(bouts)/12]

                bouts = bouts*4
                mean_bouts = np.mean(bouts)/60 #min
                all_mice_mean_bout_duration[mouse] = mean_bouts

                if bouts.size >0 :
                    if np.max(bouts) > 2**12:
                        bins_sleep = np.append(2**np.arange(2,12), np.max(bouts))
                    else :
                        bins_sleep = 2**np.arange(2,13)
                    bins_sleep= np.append(np.array([0]), bins_sleep)

                    #bins_sleep += [max(bouts)]
                    # count, bins= np.histogram(bouts, bins = bins_sleep)
                    state_percentage_frag = []
                    count = []
                    # print(bins_sleep)
                    for b in np.arange(1,bins_sleep.size):
                        # print(bins_sleep[b])
                        # print(bins_sleep[b+1])
                        # print()
                        mask = (bouts>bins_sleep[b-1]) & (bouts<=bins_sleep[b])
                        # print(np.unique(bouts[mask]))
                        # print(bins_sleep[b])
                        state_percentage_frag.append(np.sum(bouts[mask])/(60*state_time))
                        count.append(np.sum(mask))
                    state_percentage_frag = np.array(state_percentage_frag)
                    count = np.array(count)

                    # exit()

                    all_mice_count[mouse] = count
                    all_mice_percent_frag[mouse] = state_percentage_frag
                    # fig, ax = plt.subplots()
                    # ax.plot(count)
                    # axb=ax.twinx()
                    # axb.plot(state_percentage_frag, color ='k')
                    # plt.show()
                    # exit()

                else :
                    bins_sleep = 2**np.arange(2,13)
                    bins_sleep= np.append(np.array([0]), bins_sleep)
                    count = np.zeros(bins_sleep.size-1)
                    state_percentage_frag = np.zeros(bins_sleep.size-1)
                    # count_vs_total_time = np.zeros(bins_sleep.size-1)

                    all_mice_count[mouse] = count
                    all_mice_percent_frag[mouse] = state_percentage_frag

            bins_sleep = 2**np.arange(2,13)
            df_percent_frag = pd.DataFrame.from_dict(all_mice_percent_frag, orient = 'index', columns = bins_sleep )
            if exclude_mice :
                output_path = excel_dir + '/RT_PCR_excluded_mice/fourfirsthours/sleep_percent_fragmentation/{}/{}/{}/'.format(miRNA, genotype, state)
            else :
                output_path = excel_dir + '/fourfirsthours/sleep_percent_fragmentation/{}/{}/{}/'.format(miRNA, genotype, state)
            if not os.path.exists(output_path):
                os.makedirs(output_path)
            df_percent_frag.to_excel(output_path +'{}_{}_{}_{}.xlsx'.format(miRNA, genotype, state,cycle))

            df = pd.DataFrame.from_dict(all_mice_count, orient = 'index', columns = bins_sleep )
            if exclude_mice :
                output_path = excel_dir + '/RT_PCR_excluded_mice/fourfirsthours/sleep_fragmentation/{}/{}/{}/'.format(miRNA, genotype, state)
            else :
                output_path = excel_dir + '/fourfirsthours/sleep_fragmentation/{}/{}/{}/'.format(miRNA, genotype, state)
            if not os.path.exists(output_path):
                os.makedirs(output_path)
            df.to_excel(output_path +'{}_{}_{}_{}.xlsx'.format(miRNA, genotype, state,cycle))

            df_duration = pd.DataFrame.from_dict(all_mice_mean_bout_duration, orient = 'index', columns = ['mean_duration'] )
            if exclude_mice :
                output_path = excel_dir + '/RT_PCR_excluded_mice/fourfirsthours/fourfirsthours/state_duration/{}/{}/{}/'.format(miRNA, genotype, state)
            else:
                output_path = excel_dir + '/fourfirsthours/state_duration/{}/{}/{}/'.format(miRNA, genotype, state)
            if not os.path.exists(output_path):
                os.makedirs(output_path)
            df_duration.to_excel(output_path +'{}_{}_{}_{}.xlsx'.format(miRNA, genotype, state,cycle))

            df_time = pd.DataFrame.from_dict(all_mice_state_time, orient = 'index', columns = ['time'] )
            if exclude_mice:
                output_path = excel_dir + '/RT_PCR_excluded_mice/fourfirsthours/state_time/{}/{}/{}/'.format(miRNA, genotype, state)
            else :
                output_path = excel_dir + '/fourfirsthours/state_time/{}/{}/{}/'.format(miRNA, genotype, state)
            if not os.path.exists(output_path):
                os.makedirs(output_path)
            df_time.to_excel(output_path +'{}_{}_{}_{}.xlsx'.format(miRNA, genotype, state,cycle))

            df_period_amount = pd.DataFrame.from_dict(all_mice_bout_amount, orient = 'index', columns = ['period', 'freq'] )
            if exclude_mice :
                output_path = excel_dir + '/RT_PCR_excluded_mice/fourfirsthours/period_amount/{}/{}/{}/'.format(miRNA, genotype, state)
            else :
                output_path = excel_dir + '/fourfirsthours/period_amount/{}/{}/{}/'.format(miRNA, genotype, state)
            if not os.path.exists(output_path):
                os.makedirs(output_path)
            df_period_amount.to_excel(output_path +'{}_{}_{}_{}.xlsx'.format(miRNA, genotype, state,cycle))



def tdw_wake_ratio(miRNA, genotype, exclude_mice = False):
    if exclude_mice :
        path = precompute_dir+'/RT_PCR_excluded_mice/sleep_by_epoch/{}_{}_exclude_mice_sleep_by_epoch.nc'.format(miRNA, genotype)
    else :
        path = precompute_dir+'/sleep_by_epoch/{}_{}_sleep_by_epoch.nc'.format(miRNA, genotype)
    ds = xr.open_dataset(path)
    level = [0, 1]
    # print(ds)
    all_mice_w_epoch = ds['wake_all_mice'].to_pandas()
    all_mice_t_epoch = ds['tdw_all_mice'].to_pandas()

    number_of_epochs = all_mice_w_epoch.shape[0]
    epoch_duration = 4
    duration_in_hours = int(number_of_epochs*epoch_duration/3600)
    time = int(1) # window to look at in HOUR. CAUTION must be a multiple of 100
    hours = np.arange(0, int(duration_in_hours),time)
    window = int(time * 3600 / epoch_duration)

    fake = np.zeros(number_of_epochs)

    # all_mice_w_by_time = []
    # all_mice_t_by_time = []
    all_mice_t_w_ratio_by_time = []
    for hour in hours :
        i1, i2 = int(hour*window), int((hour+1)*window)
        # tarace[hour] = all_mice_w_epoch[i1:i2].sum(axis =0)
        w = all_mice_w_epoch[i1:i2].sum(axis =0)*4/60
        t = all_mice_t_epoch[i1:i2].sum(axis =0)*4/60
        # all_mice_w_by_time.append(w)
        # all_mice_t_by_time.append(t)
        all_mice_t_w_ratio_by_time.append(t/w)
    halfdays = np.arange(8)
    one_hour_to_window =  int(3600 / epoch_duration)
    all_mice_t_w_ratio_by_halfday = []
    for halfday in halfdays:
        i1 = int((12*halfday+4)*one_hour_to_window)
        i2 = int((12*(halfday+1)+4)*one_hour_to_window)
        # tarace[hour] = all_mice_w_epoch[i1:i2].sum(axis =0)
        w = all_mice_w_epoch[i1:i2].sum(axis =0)*4/60
        t = all_mice_t_epoch[i1:i2].sum(axis =0)*4/60
        # all_mice_w_by_time.append(w)
        # all_mice_t_by_time.append(t)
        all_mice_t_w_ratio_by_halfday.append(t/w)


    all_mice_t_w_ratio_by_time = pd.concat(all_mice_t_w_ratio_by_time, axis = 1)
    all_mice_t_w_ratio_by_halfday = pd.concat(all_mice_t_w_ratio_by_halfday, axis = 1)
    if exclude_mice :
        dirname = excel_dir+  '/RT_PCR_excluded_mice/tdw_w_ratio/{}/{}/'.format(miRNA, genotype)
    else:
        dirname = excel_dir+  '/tdw_w_ratio/{}/{}/'.format(miRNA, genotype)

    if not os.path.exists(dirname):
        os.makedirs(dirname)
    all_mice_t_w_ratio_by_time.to_excel(dirname+ '/tdw_w_ratio_by_mouse_by_hour.xlsx')
    all_mice_t_w_ratio_by_halfday.to_excel(dirname+ '/tdw_w_ratio_by_mouse_by_12hours.xlsx')



def inter_REM_interval(miRNA, genotype, exclude_mice = False):
    # B0 begins at 4h, B1 etc begins at 8
    if exclude_mice :
        path = precompute_dir+'/RT_PCR_excluded_mice/sleep_by_epoch/{}_{}_exclude_mice_sleep_by_epoch.nc'.format(miRNA, genotype)
    else :
        path = precompute_dir+'/sleep_by_epoch/{}_{}_sleep_by_epoch.nc'.format(miRNA, genotype)
    # print(precompute_dir)
    # print(path)
    ds = xr.open_dataset(path)

    all_mice_r_epoch = ds['rem_all_mice'].to_pandas()    #.unstack().T*1
    number_of_epochs = ds['rem_all_mice'].shape[0]

    epoch_duration = 4
    duration_in_hours = int(number_of_epochs*epoch_duration/3600)
    day_cycles = {'dark0': [0,4],
                 'light1' : [4, 16],
                 'dark1' : [16, 28],
                 'light2' : [28, 40],
                 'dark2' : [40, 52],
                 'sd' : [52, 58],
                 'light3' : [58,64],
                 'dark3':[64, 76],
                 'light4':[76,88 ],
                 'dark4':[88, 100]}
    mice = ds.coords['mice'].values

    data = all_mice_r_epoch
    for c, cycle in enumerate(day_cycles):
        cycle_hours = day_cycles[cycle]
        # print(cycle_hours)
        i1, i2 = int(cycle_hours[0]), int(cycle_hours[1])
        i1 = int(i1 *3600/epoch_duration)
        i2 = int(i2 *3600/epoch_duration)
        selected_data = data[i1: i2]
        # print(diff.where(diff == -1))
        all_mice_IRI = {}

        # print('je suis la')
        for mouse in mice :
            one_and_zeros_one_mouse = selected_data[mouse].values
            # print(np.sum(one_and_zeros_one_mouse))
            diff = np.diff(one_and_zeros_one_mouse)
            ini = np.where(diff ==1)[0]+1
            if one_and_zeros_one_mouse[0]==1:
                ini = np.append(np.array([0]), ini)
            end = np.where(diff ==-1)[0]
            if one_and_zeros_one_mouse[-1]==1:
                end = np.append(end, np.array([one_and_zeros_one_mouse.size]))
            IRI = np.mean(ini[1:]-end[:-1])
            IRI *= (4/60)
            all_mice_IRI[mouse]=IRI
            # fig, ax = plt.subplots()
            # ax.plot(diff)
            # ax.plot(one_and_zeros_one_mouse)
            # ax.scatter(ini, np.ones(ini.size), color ='green')
            # ax.scatter(end, np.ones(end.size), color ='red')
            # for i,iri in enumerate(IRI):
            #     ax.plot([end[:-1][i], end[:-1][i]+iri],[1,1])
            # plt.show()
            # exit()

        df_IRI = pd.DataFrame.from_dict(all_mice_IRI, orient = 'index', columns = ['IRI'] )
        if exclude_mice:
            output_path = excel_dir + '/RT_PCR_excluded_mice/IRI/{}/{}/'.format(miRNA, genotype)
        else :
            output_path = excel_dir + '/IRI/{}/{}/'.format(miRNA, genotype)
        if not os.path.exists(output_path):
            os.makedirs(output_path)
        df_IRI.to_excel(output_path +'{}_{}_{}.xlsx'.format(miRNA, genotype,cycle))



def REM_sleep_latency_one_mouse(mouse,exclude_mice):
    min_wake_duration = 6
    ds = xr.open_dataset(precompute_dir + '/spectrums/spectrum_scoring_{}.nc'.format(mouse))
    times = ds['time_epochs'].values/3600
    score = ds['score'].values
    score_behavior = score.copy()
    miRNA, genotype = get_mouse_info(mouse)
    mice = get_mice(miRNA, genotype,exclude_mice)



    halfday_times= {'dark0': [0,4],
                 'light1' : [4, 16],
                 'dark1' : [16, 28],
                 'light2' : [28, 40],
                 'dark2' : [40, 52],
                 'sd' : [52, 58],
                 'light3' : [58,64],
                 'dark3':[64, 76],
                 'light4':[76,88 ],
                 'dark4':[88, 100]}

    for f, n in zip(['1', '2', '3'], ['w', 'n', 'r']) :
        score_behavior = np.where(score_behavior == f, n, score_behavior)

    for halfday in halfday_times:
        if exclude_mice:
            dirname = excel_dir + '/RT_PCR_excluded_mice/REM_latency/{}/{}/'.format(miRNA, genotype)
        else :
            dirname = excel_dir + '/REM_latency/{}/{}/'.format(miRNA, genotype)
        if not os.path.exists(dirname):
            os.makedirs(dirname)
        filename =dirname+'REM_latency_{}_{}_{}.xlsx'.format(miRNA, genotype, halfday)
        if not os.path.exists(filename):
            df = pd.DataFrame(index = mice, columns = np.arange(5))
        else :
            df = pd.read_excel(filename, index_col = 0)


        t1 = halfday_times[halfday][0]
        t2 = halfday_times[halfday][1]
        time_mask = (times>t1) & (times<t2)

        masked_score = score_behavior[time_mask]
        rem = score_behavior[time_mask] == 'r'
        nrem = score_behavior[time_mask] == 'n'
        wake = score_behavior[time_mask] == 'w'

        # fig, ax = plt.subplots()
        latencies = []
        ini_rem = np.diff(rem*1)
        pos_ini_rem = np.where(ini_rem == 1)[0]+1
        for pos in pos_ini_rem:
            wake_duration = 0
            for latency, ind in enumerate(np.arange(pos-1)[::-1]):
                s = masked_score[ind]
                if s == 'w':
                    wake_duration += 1
                if s == 'w' and wake_duration >=min_wake_duration:
                    if latency == 5:
                        pos_ini_rem = pos_ini_rem[pos_ini_rem!=pos]
                        break
                    else :
                        latency -= 5
                        latencies.append(latency)
                        break
                elif s == 'n' and wake_duration < min_wake_duration:
                    wake_duration = 0
                elif s =='r' and wake_duration == latency:
                    pos_ini_rem = pos_ini_rem[pos_ini_rem!=pos]
                    break
                elif s == 'r' and wake_duration != latency:
                    if masked_score[ind+1] == 'w':
                        latency -= wake_duration
                        latencies.append(latency)
                    else :
                        latencies.append(latency)
                    break
                    # break

            # ax.plot(np.arange(pos-latency, pos), np.ones(latency), color = 'black')
        # print(df.columns.to_list())
        # print(len(latencies))
        # print(latencies)
        if len(df.columns.to_list())<len(latencies):
            df = df.reindex(columns = np.arange(len(latencies)))
        df.loc[mouse, np.arange(len(latencies))] = np.array(latencies)
        df.to_excel(filename)


def REM_sleep_latency_by_hour(miRNA, genotype, exclude_mice = False):
    min_wake_duration = 6
    epoch_duration = 4

    mice = get_mice(miRNA, genotype, exclude_mice)
    df = pd.DataFrame(index = mice, columns = np.arange(100))
    for mouse in mice :
        ds = xr.open_dataset(precompute_dir + '/spectrums/spectrum_scoring_{}.nc'.format(mouse))
        times = ds['time_epochs'].values/3600
        score = ds['score'].values

        score_behavior = score.copy()
        miRNA, genotype = get_mouse_info(mouse)

        for f, n in zip(['1', '2', '3'], ['w', 'n', 'r']) :
            score_behavior = np.where(score_behavior == f, n, score_behavior)
        hours = np.arange(0, int(round(times[-1])))

        window = int(3600 / epoch_duration)

        for hour in hours :
            i1, i2 = int(hour*window), int((hour+1)*window)

            t1 = hour
            t2 = hour+1
            time_mask = (times>t1) & (times<t2)

            masked_score = score_behavior[time_mask]
            rem = score_behavior[time_mask] == 'r'
            nrem = score_behavior[time_mask] == 'n'
            wake = score_behavior[time_mask] == 'w'

            # fig, ax = plt.subplots()
            latencies = []
            ini_rem = np.diff(rem*1)
            pos_ini_rem = np.where(ini_rem == 1)[0]+1
            for pos in pos_ini_rem:
                wake_duration = 0
                for latency, ind in enumerate(np.arange(pos-1)[::-1]):
                    s = masked_score[ind]
                    if s == 'w':
                        wake_duration += 1
                    if s == 'w' and wake_duration >=min_wake_duration:
                        if latency == 5:
                            pos_ini_rem = pos_ini_rem[pos_ini_rem!=pos]
                            break
                        else :
                            latency -= 5
                            latencies.append(latency)
                            break
                    elif s == 'n' and wake_duration < min_wake_duration:
                        wake_duration = 0
                    elif s =='r' and wake_duration == latency:
                        pos_ini_rem = pos_ini_rem[pos_ini_rem!=pos]
                        break
                    elif s == 'r' and wake_duration != latency:
                        if masked_score[ind+1] == 'w':
                            latency -= wake_duration
                            latencies.append(latency)
                        else :
                            latencies.append(latency)
                        break
                        # break
            latencies = 4*np.array(latencies)/60
            if len(latencies)>0:
                latency = np.mean(latencies)
            else :
                latency = np.nan
            df.loc[mouse, hour] = latency

    if exclude_mice:
        dirname = excel_dir + '/RT_PCR_excluded_mice/REM_latency_by_hour/{}/{}/'.format(miRNA, genotype)
    else :
        dirname = excel_dir + '/REM_latency_by_hour/{}/{}/'.format(miRNA, genotype)
    if not os.path.exists(dirname):
        os.makedirs(dirname)

    df.to_excel(dirname + '/REM_latency_by_hour.xlsx')


def compute_all_REM_latency():
        for miRNA in ['128', '137', '665']:
            for genotype in ['test', 'control']:
                for mouse in get_mice(miRNA, genotype):
                    REM_sleep_latency_one_mouse(mouse)


if __name__ == '__main__' :
    # precompute_sleep_state_by_epoch('128', 'control')
    # REM_sleep_latency_one_mouse('B2884')
    # compute_all_REM_latency()
    compute_all()
    compute_all(exclude_mice =True)
    # sleep_bouts('128', 'control')
    # sleep_state_statistics('128', 'control', True)
    # tdw_wake_ratio('128', 'control')
    # inter_REM_interval('128', 'control')
    # sleep_state_statistics_all()
    # plt.show()

    # REM_sleep_latency_by_hour('128', 'control', exclude_mice = False)
