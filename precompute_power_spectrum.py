
from configuration import *
import scipy.signal
from select_mice_cata_Malo import get_mice
import mne.time_frequency
import multiprocessing
from joblib import Parallel, delayed


windows = 4   #### in s

#
def vectorial_power_spectrum(mouse):
    print(mouse)
    ds = xr.open_dataset(precompute_dir + '/raw/raw_{}.nc'.format(mouse))
    raw = ds['signal'].values
    times = ds.coords['times_second'].values
    t_start = times[0]
    sr = ds['sampling_rate'].values
    # n_epochs = ((times[-1] - times[0]))//windows
    n_epochs = int(raw.size//(windows*sr))

    point_per_epochs = int(4*sr)#800    ##### 4 sec at about 200 HZ   199,6

    ###### highpass filter 0.3 Hz
    N =3
    f_cut = .3
    nyq = sr/2
    W = f_cut/nyq
    b, a = scipy.signal.butter(N, W, btype = 'highpass', output = 'ba')
    raw = scipy.signal.filtfilt(b,a, raw)



    bandwidth = 1.01/windows
    epochs = np.arange(n_epochs)
    sample_per_epoch = int(windows*sr)


    point_per_epochs = 800


    sample_per_2epoch = int(100*windows*sr)
    n_2epochs = int(raw.size//(100*windows*sr))
    rest = int(raw.size - sample_per_2epoch*n_2epochs)
    # rest = int(raw.size - sample_per_epoch*n_epochs)
    print(rest)
    raw = raw[:-rest]
    stacked_sigs = raw.copy()
    # stacked_sigs = stacked_sigs.reshape((int(n_epochs), sample_per_epoch))
    print(stacked_sigs.shape)
    stacked_sigs = stacked_sigs.reshape((int(n_2epochs), sample_per_2epoch))
    # filling = (stacked_sigs[:,-1])
    # filling = filling.reshape((filling.size,1))
    filling = stacked_sigs[:,-6:]
    stacked_sigs = np.concatenate((stacked_sigs,filling), axis = 1)
    stacked_sigs = stacked_sigs.reshape((int(n_epochs), point_per_epochs))
    print(stacked_sigs.shape)
    # exit()
    ###### 6hz oscillation trial
    # stacked_sigs = np.arange(8000)*6*2*np.pi/sr
    # stacked_sigs = stacked_sigs.reshape(10,800)
    # stacked_sigs = np.sin(stacked_sigs)

    freqs_welch, welch = scipy.signal.welch(stacked_sigs, fs = sr, nperseg = int(3.99*sr) )
    welch = welch[:,(freqs_welch>=.75) & (freqs_welch<=47.5)]
    freqs_welch = freqs_welch[(freqs_welch>.75) & (freqs_welch<=47.5)]
    print('here')
    fft = np.fft.fft(stacked_sigs)
    fft = np.abs(fft)**2
    # fft = fft[:,0:int(fft.size/2)]
    # freqs_fft = np.linspace(0, sr/2, point_per_epochs/2)
    # print(fft.shape, freqs_fft.shape)
    freqs_fft = np.fft.fftfreq(point_per_epochs, 1/sr)
    idx = np.argsort(freqs_fft)
    freqs_fft = freqs_fft[idx]
    fft = fft[:,idx]
    fft = fft[:,(freqs_fft>=.75)&(freqs_fft<=47.5)]
    freqs_fft = freqs_fft[(freqs_fft>=.75)&(freqs_fft<=47.5)]

    multitaper, freqs_multitaper = mne.time_frequency.psd_array_multitaper(stacked_sigs, sfreq = sr, fmax = 47.5, bandwidth = 4*bandwidth, n_jobs = 2)
    multitaper = multitaper[:,freqs_multitaper>=.75]
    freqs_multitaper = freqs_multitaper[freqs_multitaper>=.75]
    # times = np.arange(point_per_epochs)/sr
    # for i in np.arange(2000):
    #     fig, ax = plt.subplots(nrows = 4)
    #     ax[0].plot(times, stacked_sigs[i])
    #     ax[1].plot(freqs_multitaper, multitaper[i])
    #     ax[1].set_ylabel('multitaper')
    #     ax[2].plot(freqs_welch, welch[i])
    #     ax[2].set_ylabel('welch')
    #     ax[3].plot(freqs_fft, fft[i])
    #     ax[3].set_ylabel('fft')
    #     plt.show()

    coords = {  'freqs_welch': freqs_welch,
                'freqs_fft': freqs_fft,
                'freqs_multitaper': freqs_multitaper,
                'epochs' : epochs}
    ds_spectrum = xr.Dataset(coords = coords)
    ds_spectrum['welch_spectrum'] = xr.DataArray(welch, dims = ['epochs', 'freqs_welch'])
    ds_spectrum['fft_spectrum'] = xr.DataArray(fft, dims = ['epochs', 'freqs_fft'])
    ds_spectrum['multitaper_spectrum'] = xr.DataArray(multitaper, dims = ['epochs', 'freqs_multitaper'])
    ds_spectrum['sampling_rate'] = sr
    saving_path = precompute_dir + '/spectrums/'
    if not os.path.exists(saving_path):
        os.makedirs(saving_path)
    print(saving_path)
    ds_spectrum.to_netcdf(saving_path + 'spectrums_{}.nc'.format(mouse))

def compute_one_welch(raw,times, w, t_start, sr):
    mask = (times>t_start + windows*w) & (times<=t_start + windows*(w+1))
    sig = raw[mask]
    sig = sig - np.mean(sig)
    freqs_welch, pxx = scipy.signal.welch(sig, fs = sr, nperseg = int(3.99*sr) )
    pxx = pxx[(freqs_welch>.4) & (freqs_welch<=10)]
    freqs_welch = freqs_welch[(freqs_welch>.4) & (freqs_welch<=10)]
    return (int(w), np.array(pxx))

def read_one_mouse_spectrums(mouse):
    ds = xr.open_dataset(precompute_dir + '/spectrums/spectrums_{}.nc'.format(mouse))
    print(ds)
    freqs = ds.coords['freqs'].values
    spectrums = ds['spectrum_by_epoch'].values
    fig, ax = plt.subplots()
    ax.plot(freqs, spectrums[15])
    plt.show()

def compute_all_spectrums_data():
    date_ref = pd.read_excel(work_dir + 'datetime_reference_DICER.xls', index_col = 0)
    mice = date_ref.index.to_list()
    for mouse in mice:
        mouse = mouse[4:]
        vectorial_power_spectrum(mouse)

def hack_dask():
    run_key = sys.argv[1]
    print(run_key, type(run_key))
    sliding_power_spectrum(run_key)

def store_all_score_and_spectrum():
    dcr_mice = get_mice(group = 'DCR-HCRT')
    control_mice = get_mice(group = 'Control')
    animals_by_group = {'DCR-HCRT' : dcr_mice, 'Control' : control_mice}
    for group in animals_by_group :
        mice = animals_by_group[group]
        results = Parallel(n_jobs=2)(delayed(store_scoring_and_spectrums_one_mouse_one_session)(group,mouse) for mouse in mice)
        # for mouse in mice :
        #     print(mouse)
        #     # for rec in ['b1', 'b2', 'sd', 'r1']:
        #         # try:
        #     store_scoring_and_spectrums_one_mouse_one_session(group, mouse )
                # except :
                    # print('******* ERROR ******', mouse,rec, '*********')

def compute_one_group(group =  'Control'):
    mice = get_mice_for_spectrum(group)
    for mouse in mice :
        print(mouse)
        for rec in ['b1', 'b2', 'sd', 'r1']:
            try:
                store_scoring_and_spectrums_one_mouse_one_session(group, mouse, rec)
            except :
                print('******* ERROR ******', mouse,rec, '*********')


def store_scoring_and_spectrums_one_mouse_one_session(group, mouse):
    # print( 'compute {} session'.format(rec))
    #######         Extract scoring         #######

    score_b1 = np.loadtxt(data_dir + '/Scoring/' + group + '/' + mouse + 'DCRb1.txt', dtype = str)
    score_b2 = np.loadtxt(data_dir + '/Scoring/' + group + '/' + mouse + 'DCRb2.txt', dtype = str)
    score_sd = np.loadtxt(data_dir + '/Scoring/' + group + '/' + mouse + 'DCRsd.txt', dtype = str)
    score_r1 = np.loadtxt(data_dir + '/Scoring/' + group + '/' + mouse + 'DCRr1.txt', dtype = str)
    all_score = np.zeros(score_b1.size + score_b2.size + score_sd.size + score_r1.size, dtype = str)
    one_day = score_b1.size
    all_score[:one_day] = score_b1
    all_score[one_day:int(one_day*2)] = score_b2
    all_score[int(2*one_day):int(3*one_day)] = score_sd
    all_score[int(3*one_day):int(4*one_day)] = score_r1
    file_b1 = open('/Users/maximejuventin/Desktop/power_spectrum_somnologica/Power_spectrum/{}/MTA-{}/{}DCRb1s.txt'.format(group, mouse, mouse))
    file_b2 = open('/Users/maximejuventin/Desktop/power_spectrum_somnologica/Power_spectrum/{}/MTA-{}/{}DCRb2s.txt'.format(group, mouse, mouse))
    file_sd = open('/Users/maximejuventin/Desktop/power_spectrum_somnologica/Power_spectrum/{}/MTA-{}/{}DCRsds.txt'.format(group, mouse, mouse))
    file_r1 = open('/Users/maximejuventin/Desktop/power_spectrum_somnologica/Power_spectrum/{}/MTA-{}/{}DCRr1s.txt'.format(group, mouse, mouse))
    # file_r1 = open(data_spectrum_dir + '/{}/MTA-{}/{}DCRr1s.txt'.format(group, mouse, mouse))
    real_time_lines = []
    data_lines = []
    real_times_somno = []
    data = []
    delta_power_somno = []
    for file in [file_b1, file_b2, file_sd, file_r1] :
        for i, line in enumerate(file.readlines()) :
            if line == 'Power Spectrum\n':
                real_time_lines.append(i+1)
                data_lines.append(i+2)
            if i in real_time_lines:
                real_times_somno.append(line[:-1])
            if i in data_lines :
                spec = np.array(line.split('\t')[:-1], dtype = float)[4:]    ####Remove .5 Hz     #######
                data.append(spec)
                delta_power_somno.append(np.sum(spec[:15])) #####delta = .75 to 4Hz

    print('finish')

    freqs_somno = np.arange(4, spec.size+4, 1) * .25
    times_somno = np.arange(len(real_times_somno))*int(4)

    somno_spectrums = np.array(data)

    somno_spectrums = somno_spectrums[:,(freqs_somno >=.75) & (freqs_somno <= 47.5)]
    freqs_somno = freqs_somno[(freqs_somno >=.75) & (freqs_somno <= 47.5)]
    print(somno_spectrums.shape)



    ######      Caution ! Do not remove 1, 2, 3. It correspond to invalid EEG
    # for f, n in zip(['1', '2', '3'], ['w', 'n', 'r']) :
    #     score = np.where(score == f, n, score)

    #######         Extract spectrums from txt        #######
    ds_spectrum = xr.open_dataset(precompute_dir + '/spectrums/spectrums_{}.nc'.format(mouse))
    print(ds_spectrum)
    welch_spectrums = ds_spectrum['welch_spectrum']
    freqs_welch = ds_spectrum.coords['freqs_welch'].values

    fft_spectrums = ds_spectrum['fft_spectrum']
    freqs_fft = ds_spectrum.coords['freqs_fft'].values

    multitaper_spectrums = ds_spectrum['multitaper_spectrum']
    freqs_multitaper = ds_spectrum.coords['freqs_multitaper'].values

    welch_mask = (freqs_welch>1.) & (freqs_welch<=4)
    welch_delta_power = welch_spectrums[:,welch_mask]
    welch_delta_power = welch_delta_power.sum(axis =1)

    fft_mask = (freqs_fft>1.) & (freqs_fft<=4)
    fft_delta_power = fft_spectrums[:,fft_mask]
    fft_delta_power = fft_delta_power.sum(axis =1)

    multitaper_mask = (freqs_multitaper>1.) & (freqs_multitaper<=4)
    multitaper_delta_power = multitaper_spectrums[:,multitaper_mask]
    multitaper_delta_power = multitaper_delta_power.sum(axis =1)

    #######         Storing          #######
    # times = np.arange(len(real_times))*int(4)

    # coords = {'real_times' : real_times, 'times' : times, 'freqs' : freqs}
    coords = ds_spectrum.coords
    coords['real_times_somno'] = real_times_somno
    coords['freqs_somno'] = freqs_somno
    coords['times_somno'] = times_somno


    ds = xr.Dataset(coords = ds_spectrum.coords)
    # ds = xr.Dataset(coords = coords)
    ds['somno_spectrum'] = xr.DataArray(somno_spectrums, dims = ['epochs', 'freqs_somno'])
    ds['somno_delta_power'] = xr.DataArray(np.array(delta_power_somno), dims = 'epochs')

    ds['welch_spectrum'] = welch_spectrums
    ds['welch_delta_power'] = welch_delta_power

    ds['fft_spectrum'] = fft_spectrums
    ds['fft_delta_power'] = fft_delta_power

    ds['multitaper_spectrum'] = multitaper_spectrums
    ds['multitaper_delta_power'] = multitaper_delta_power

    ds['score'] = xr.DataArray(all_score, dims = 'epochs')

    dirname = precompute_dir + '/spectrums/'
    if not os.path.exists(dirname):
        os.makedirs(dirname)
    print(dirname)
    ds.to_netcdf(dirname + 'spectrum_scoring_{}.nc'.format(mouse))

if __name__ == '__main__':
    mouse = 'B2533'
    # mouse = 'B4907'
    # mouse  ='B2763'


    # hack_dask()
    # vectorial_power_spectrum(mouse)
    # read_one_mouse_spectrums(mouse)
    # compute_all_spectrums_data()
    store_all_score_and_spectrum()
