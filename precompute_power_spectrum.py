
from configuration import *
import scipy.signal
from select_mice_cata_Malo import get_mice, get_mouse_info
import mne.time_frequency
import multiprocessing
from joblib import Parallel, delayed


windows = 4   #### in s


def compute_own_spectrum(mouse):
    print('compute own spectrum mouse : ',mouse )
    ds = xr.open_dataset(precompute_dir + '/raw/raw_{}.nc'.format(mouse))
    raw = ds['signal'].values.astype('float32')
    sr = ds['sampling_rate'].values
    #
    # times = ds.coords['times_second'].values
    times = np.arange(raw.size)/(sr*3600)

    t_start = times[0]
    n_epochs = int(raw.size//(windows*sr))

    point_per_epochs = int(4*sr)#800    ##### 4 sec at about 200 HZ   199,6
    ###### highpass filter 0.5 Hz
    N =3
    f_cut = .5
    nyq = sr/2
    W = f_cut/nyq
    b, a = scipy.signal.butter(N, W, btype = 'highpass', output = 'ba')
    raw = scipy.signal.filtfilt(b,a, raw)
    #
    bandwidth = 1.01/windows
    epochs = np.arange(n_epochs)
    sample_per_epoch = int(windows*sr)
    mylist = []
    index = []
    for i in range(n_epochs):
        mylist.append(4*i*sr)
    real_sample_by_epoch= np.diff(np.array(mylist, dtype='int'))
    ind = 0
    for fr in real_sample_by_epoch:
        fr = int(fr)
        ind += fr
        if fr == 799:
            index.append(ind)
    index = np.array(index, dtype ='int')+1
    ref = 69120000
    if raw.size + index.size != ref:
        index = index[:-(raw.size + index.size-ref)]
    raw = np.insert(raw, index, raw[index])
    point_per_epochs = 800
    stacked_sigs = raw.reshape((-1, point_per_epochs)).astype('float32')

    print(stacked_sigs.dtype)
    #
    freqs_welch, welch = scipy.signal.welch(stacked_sigs, fs = sr, nperseg = int(3.99*sr) )
    welch = welch[:,(freqs_welch>=.75) & (freqs_welch<=47.5)]
    freqs_welch = freqs_welch[(freqs_welch>.75) & (freqs_welch<=47.5)]
    print(welch.shape)
    print('here')

    # fft = scipy.fft.fft(stacked_sigs)
    # Fourier = np.fft.fft(stacked_sigs)
    # Fourier = np.fft.fft(np.random.rand(welch.shape))
    # Fourier = np.fft.fft(np.random.rand(10,10))
    # Fourier = F(stacked_sigs)
    # print(Fourier.shape)
    # Fourier = np.abs(Fourier)**2
    # Fourier =  np.abs(np.fft.fft(stacked_sigs))**2
    # freqs_Fourier = np.fft.fftfreq(point_per_epochs, 1/sr)
    # idx = np.argsort(freqs_Fourier)
    # freqs_Fourier = freqs_Fourier[idx]
    # Fourier = Fourier[:,idx]
    # Fourier = Fourier[:,(freqs_Fourier>=.75)&(freqs_Fourier<=47.5)]
    # freqs_Fourier = freqs_Fourier[(freqs_Fourier>=.75)&(freqs_Fourier<=47.5)]

    multitaper, freqs_multitaper = mne.time_frequency.psd_array_multitaper(stacked_sigs, sfreq = sr, fmax = 47.5, bandwidth = 4*bandwidth, n_jobs = 2)
    multitaper = multitaper[:,freqs_multitaper>=.75]
    freqs_multitaper = freqs_multitaper[freqs_multitaper>=.75]

    freqs_welch  = freqs_welch.astype('float32')
    # freqs_Fourier  = freqs_Fourier.astype('float32')
    freqs_multitaper  = freqs_multitaper.astype('float32')
    epochs = epochs.astype('int32')
    welch = welch.astype('float32')
    # Fourier = Fourier.astype('float32')
    multitaper = multitaper.astype('float32')

    mydict = {  'freqs_welch' : freqs_welch,
                # 'freqs_Fourier' : freqs_Fourier,
                'freqs_multitaper' : freqs_multitaper,
                'epochs' : epochs,
                'welch' : welch,
                # 'Fourier' : Fourier,
                'multitaper' : multitaper}
    return mydict

def store_scoring_and_spectrums_one_mouse(mouse):
    # print( 'compute {} session'.format(rec))
    # date_ref = pd.read_excel(work_dir + 'datetime_reference_miRNA.xlsx', index_col = 0)
    # if date_ref.at['MTA-{}'.format(mouse), 'short_end'] == 'yes':
    #     print('END exception')
    #     n = int((3*24 +12+4)*3600*sr) +1

    miRNA, genotype = get_mouse_info(mouse)
    #######         Extract scoring         #######
    print('get scoring mouse : ',mouse )
    # filename = scoring_data_dir + '/{}/{}{}/MTA-{}/'.format(miRNA, genotype, miRNA,mouse)
    filename = work_dir + '/data/Scoring/{}/{}{}/MTA-{}/'.format(miRNA, genotype, miRNA,mouse)

    score_b0 = np.loadtxt(filename + mouse + 'b0.txt', dtype = str)
    score_b1 = np.loadtxt(filename + mouse + 'b1.txt', dtype = str)
    score_b2 = np.loadtxt(filename + mouse + 'b2.txt', dtype = str)
    score_sd = np.loadtxt(filename + mouse + 'sd.txt', dtype = str)
    score_r1 = np.loadtxt(filename + mouse + 'r.txt', dtype = str)
    print(score_b0.size, score_b1.size, score_b2.size, score_sd.size, score_r1.size)
    # all_score = np.zeros(score_b0.size + score_b1.size + score_b2.size + score_sd.size + score_r1.size, dtype = str)
    # one_day = score_b1.size
    all_score = np.concatenate((score_b0, score_b1, score_b2, score_sd, score_r1)).astype('str')
    # all_score[:one_day] = score_b1
    # all_score[:one_day] = score_b1
    # all_score[one_day:int(one_day*2)] = score_b2
    # all_score[int(2*one_day):int(3*one_day)] = score_sd
    # all_score[int(3*one_day):int(4*one_day)] = score_r1
    score_b0 = 0
    score_b1 = 0
    score_b2 = 0
    score_sd = 0
    score_r1 = 0
    # time_epochs = np.arange(all_score.size)*4

    # exit()

    # exit()
    # print('get scoring somnologica spectrum : ',mouse )
    # file_b1 = open(data_dir + '/Power_spectrum/{}/MTA-{}/{}DCRb1s.txt'.format(group, mouse, mouse))
    # file_b2 = open(data_dir + '/Power_spectrum/{}/MTA-{}/{}DCRb2s.txt'.format(group, mouse, mouse))
    # file_sd = open(data_dir + '/Power_spectrum/{}/MTA-{}/{}DCRsds.txt'.format(group, mouse, mouse))
    # file_r1 = open(data_dir + '/Power_spectrum/{}/MTA-{}/{}DCRr1s.txt'.format(group, mouse, mouse))
    # file_r1 = open(data_dir + '/Power_spectrum/{}/MTA-{}/{}DCRr1s.txt'.format(group, mouse, mouse))
    # real_time_lines = []
    # data_lines = []
    # real_times_somno = []
    # data = []
    # for file in [file_b1, file_b2, file_sd, file_r1] :
    #     for i, line in enumerate(file.readlines()) :
    #         if line == 'Power Spectrum\n':
    #             real_time_lines.append(i+1)
    #             data_lines.append(i+2)
    #         if i in real_time_lines:
    #             real_times_somno.append(line[:-1])
    #         if i in data_lines :
    #             spec = np.array(line.split('\t')[:-1], dtype = float)[4:]    ####Remove .5 Hz     #######
    #             data.append(spec)
    # score_b1 = 0
    # score_b2 = 0
    # score_sd = 0
    # score_r1 = 0
    # #
    # freqs_somno = np.arange(4, spec.size+4, 1, dtype = 'float32') * .25
    # times_somno = np.arange(len(real_times_somno), dtype='int32')*int(4)
    #
    # somno_spectrums = np.array(data, dtype = 'float32')
    #
    # somno_spectrums = somno_spectrums[:,(freqs_somno >=.75) & (freqs_somno <= 47.5)]
    # freqs_somno = freqs_somno[(freqs_somno >=.75) & (freqs_somno <= 47.5)]
    #

    print('compute own spectrum mouse : ',mouse )
    # ds = xr.open_dataset(edf_data_dir + '/precompute/raw/raw_{}.nc'.format(mouse))
    ds = xr.open_dataset(precompute_dir + '/raw/raw_{}.nc'.format(mouse))
    raw = ds['signal'].values.astype('float32')
    sr = ds['sampling_rate'].values
    #
    # times = ds.coords['times_second'].values
    times = np.arange(raw.size)/(sr*3600)

    t_start = times[0]
    n_epochs = int(raw.size//(windows*sr))
    # time_epochs = np.arange(n_epochs)*4
    print(n_epochs)
    all_score = all_score[:n_epochs]
    point_per_epochs = int(4*sr)#800    ##### 4 sec at about 200 HZ   199,6
    ###### highpass filter 0.5 Hz
    N =3
    f_cut = .5
    nyq = sr/2
    W = f_cut/nyq
    b, a = scipy.signal.butter(N, W, btype = 'highpass', output = 'ba')
    raw = scipy.signal.filtfilt(b,a, raw)
    # print(raw.size)
    bandwidth = 1.01/windows
    epochs = np.arange(n_epochs)
    sample_per_epoch = int(windows*sr)
    mylist = []
    index = []
    for i in range(n_epochs):
        mylist.append(4*i*sr)
    real_sample_by_epoch= np.diff(np.array(mylist, dtype='int'))
    ind = 0
    for fr in real_sample_by_epoch:
        fr = int(fr)
        ind += fr
        if fr == 799:
            index.append(ind)
    index = np.array(index, dtype ='int')+1
    # ref = 69120000
    ref = int(n_epochs*800)
    # print(time_epochs.size)
    # print(ref, '  ref')


    if raw.size + index.size != ref:
        index = index[:-(raw.size + index.size-ref)]

    raw = np.insert(raw, index, raw[index])

    # print(raw.size)

    point_per_epochs = 800
    stacked_sigs = raw.reshape((-1, point_per_epochs)).astype('float32')

    # print(stacked_sigs.dtype)
    #
    freqs_welch, welch = scipy.signal.welch(stacked_sigs, fs = sr, nperseg = int(3.99*sr) )
    welch = welch[:,(freqs_welch>=.75) & (freqs_welch<=47.5)]
    freqs_welch = freqs_welch[(freqs_welch>.75) & (freqs_welch<=47.5)]

    # fft = scipy.fft.fft(stacked_sigs)
    # Fourier = np.fft.fft(stacked_sigs)
    # Fourier = np.fft.fft(np.random.rand(welch.shape))
    # Fourier = np.fft.fft(np.random.rand(10,10))
    # Fourier = F(stacked_sigs)
    # print(Fourier.shape)
    # Fourier = np.abs(Fourier)**2
    # # Fourier =  np.abs(np.fft.fft(stacked_sigs))**2
    # freqs_Fourier = np.fft.fftfreq(point_per_epochs, 1/sr)
    # idx = np.argsort(freqs_Fourier)
    # freqs_Fourier = freqs_Fourier[idx]
    # Fourier = Fourier[:,idx]
    # Fourier = Fourier[:,(freqs_Fourier>=.75)&(freqs_Fourier<=47.5)]
    # freqs_Fourier = freqs_Fourier[(freqs_Fourier>=.75)&(freqs_Fourier<=47.5)]

    # multitaper, freqs_multitaper = mne.time_frequency.psd_array_multitaper(stacked_sigs, sfreq = sr, fmax = 47.5, bandwidth = 4*bandwidth, n_jobs = 2)
    # multitaper = multitaper[:,freqs_multitaper>=.75]
    # freqs_multitaper = freqs_multitaper[freqs_multitaper>=.75]

    freqs_welch  = freqs_welch.astype('float32')
    # freqs_Fourier  = freqs_Fourier.astype('float32')
    # freqs_multitaper  = freqs_multitaper.astype('float32')
    epochs = epochs.astype('int32')
    welch = welch.astype('float32')
    # Fourier = Fourier.astype('float32')
    # multitaper = multitaper.astype('float32')
    # exit()

    ######      Caution ! Do not remove 1, 2, 3. It correspond to invalid EEG
    # for f, n in zip(['1', '2', '3'], ['w', 'n', 'r']) :
    #     score = np.where(score == f, n, score)


    coords = {  'freqs_welch': freqs_welch,
                # 'freqs_fft': freqs_fft,
                # 'freqs_multitaper': freqs_multitaper,
                'epochs' : epochs,
                'time_epochs' : epochs*4,
                # 'real_times_somno' : real_times_somno,
                # 'freqs_somno' : freqs_somno
                }

    ds = xr.Dataset(coords = coords)
    # ds['somno_spectrum'] = xr.DataArray(somno_spectrums, dims = ['epochs', 'freqs_somno'])
    ds['welch_spectrum'] = xr.DataArray(welch, dims = ['epochs', 'freqs_welch'])
    # ds['fft_spectrum'] = xr.DataArray(fft, dims = ['epochs', 'freqs_fft'])
    # ds['multitaper_spectrum'] = xr.DataArray(multitaper, dims = ['epochs', 'freqs_multitaper'])
    ds['sampling_rate'] = sr

    ds['score'] = xr.DataArray(all_score, dims = 'epochs')
    dirname = precompute_dir + '/spectrums/'
    if not os.path.exists(dirname):
        os.makedirs(dirname)
    print(dirname)
    ds.to_netcdf(dirname + 'spectrum_scoring_{}.nc'.format(mouse))

# def F(x):
#     fft = np.fft.fft(x)
#     return fft
# def store_all_score_and_spectrum():
#     dcr_mice = get_mice(group = 'DCR-HCRT')
#     control_mice = get_mice(group = 'Control')
#     mice = dcr_mice+control_mice
#     results = Parallel(n_jobs=2)(delayed(store_scoring_and_spectrums_one_mouse_one_session)(mouse) for mouse in mice)

def store_all_score_and_spectrum(parallel = True):
    for miRNA in ['128', '137', '665']:
        for genotype in ['test', 'control']:
            mice = get_mice(miRNA, genotype)
            if parallel:
                results = Parallel(n_jobs=3)(delayed(store_scoring_and_spectrums_one_mouse)(mouse) for mouse in mice)
            else :
                for mouse in mice :
                    try :
                        store_scoring_and_spectrums_one_mouse(mouse)
                    except :
                        print('-----------   Trouble with {}'.format(mouse))


def theta_dominated_wake_encoding_one_mouse(mouse, spectrum_method = 'welch'):
    dirname = precompute_dir + '/spectrums/'
    ds = xr.open_dataset(dirname + 'spectrum_scoring_{}.nc'.format(mouse))
    freqs = ds['freqs_{}'.format(spectrum_method)].values
    epoch_spectrums = ds['{}_spectrum'.format(spectrum_method)].values
    wake = ds['score'].values == 'w'
    new_score = ds['score'].values.copy()
    wake_spectrums = epoch_spectrums[wake,:]
    matrix_freqs = np.tile(freqs[:,np.newaxis], wake_spectrums.shape[0]).T
    matrix_indice = np.tile(np.arange(freqs.size)[:,np.newaxis], wake_spectrums.shape[0]).T
    peak = np.max(wake_spectrums[:,(freqs>3.5) & (freqs<15)], axis = 1)
    ####Error spectrum computation
    if np.sum(peak==0) !=0 :
        print('{} spectrums are KO'.format(np.sum(peak==0)))
        correction_mask = peak != 0
        wake_index = np.where(wake)[0]
        wake[wake_index[peak==0]] = np.zeros(np.sum(peak==0))
        new_score[wake_index[peak==0]] = 'ko_spectrum'
        peak = peak[correction_mask]
        wake_spectrums = wake_spectrums[correction_mask]
        matrix_freqs = matrix_freqs[np.arange(correction_mask.size)[correction_mask]]
        matrix_indice = matrix_indice[np.arange(correction_mask.size)[correction_mask]]

    peak_freq = matrix_freqs[:,(freqs>3.5) & (freqs<15)][wake_spectrums[:,(freqs>3.5) & (freqs<15)] == peak[:,np.newaxis]]
    peak_freq = peak_freq[:, np.newaxis]
    peak_in_theta = (peak_freq>6.5) & (peak_freq<12)
    peak_freq_index = np.where(matrix_freqs == peak_freq)[1][:,np.newaxis]
    mask_theta_band = (matrix_indice>=peak_freq_index-4) & (matrix_indice<=peak_freq_index+4)
    theta_power = np.sum(wake_spectrums[mask_theta_band].reshape(-1,9), axis = 1)[:,np.newaxis]
    fullband_no_delta = np.sum(wake_spectrums[:,(freqs>3.5)&(freqs<45)], axis = 1)[:,np.newaxis]
    theta_majo = theta_power/fullband_no_delta > .228
    tdw_amongst_wake = (theta_majo & peak_in_theta).reshape(wake_spectrums.shape[0])


    #####DO the exact same thing but in loop, slower if long loop, but here it isn't that slow
    # a = np.zeros(wake_spectrums.shape[0])
    # for w in np.arange(wake_spectrums.shape[0]):
    #     mask = (freqs>3.5) & (freqs<15)
    #     f = freqs[mask][wake_spectrums[w][mask]== max(wake_spectrums[w][mask])][0]
    #     if f>6.5 and f<12 :
    #         f_ind = np.where(freqs==f)[0][0]
    #         theta_power = np.sum(wake_spectrums[w][f_ind -4:f_ind+4+1])
    #         fullband_no_delta = np.sum(wake_spectrums[w][(freqs>3.5)&(freqs<45)])
    #         majo = theta_power/fullband_no_delta
    #         majo = majo > .228 ####22,8%
    #         if majo:
    #             a[w]=1
    # tdw_amongst_wake = a.astype('bool')

    wake_epochs = np.where(wake)[0]
    tdw_epochs = wake_epochs[np.arange(wake_epochs.size)[tdw_amongst_wake]]
    new_score[tdw_epochs] = 't'

    ######### Anne's exclusion process for tdw - signel and borders#####################
    # true_tdw = 0
    # last_pos = new_score.size-1
    # # print(last_pos)
    # for pos,s in enumerate(new_score):
    #     if s == 't':
    #         validation = 0
    #         if pos !=0 and (new_score[pos-1] == 'w' or  new_score[pos-1] == 't'):
    #             validation +=1
    #         if pos != last_pos and new_score[pos+1] != 'n':
    #             validation +=1
    #         if pos not in [0,1,2, last_pos, last_pos-1,last_pos-2] and (new_score[pos-3] == 'w' and new_score[pos-2] == 'w' and new_score[pos-1] == 'w' and new_score[pos+1] == 'w' and new_score[pos+2] == 'w' and new_score[pos+3] == 'w'):
    #             validation = 0
    #         if validation == 2:
    #             true_tdw +=1
    # print(true_tdw/np.sum(new_score=='t')*100)
    ###########################################################################

    coords = {'epochs' : np.arange(new_score.size)}
    ds_tdw = xr.Dataset(coords = coords)
    ds_tdw['new_score'] = xr.DataArray(new_score, dims = 'epochs')
    dirname = precompute_dir + '/tdw_score/'
    if not os.path.exists(dirname):
        os.makedirs(dirname)
    ds_tdw.to_netcdf(dirname + 'tdw_score_{}.nc'.format(mouse))

def theta_dominated_wake_encoding_all_mouse(parallel = False):
    for miRNA in ['128', '137', '665']:
        for genotype in ['test', 'control']:
            mice = get_mice(miRNA, genotype)
            print(mice)
            if parallel:
                results = Parallel(n_jobs=3)(delayed(theta_dominated_wake_encoding_one_mouse)(mouse) for mouse in mice)
                # print('heloo')
            else :
                for mouse in mice :
                    theta_dominated_wake_encoding_one_mouse(mouse)

                    # try :
                    #     theta_dominated_wake_encoding_one_mouse(mouse)
                    # except :
                    #     print('-----------   Trouble with {}'.format(mouse))


def HF_store_scoring_and_spectrums_one_mouse(mouse):
    # print( 'compute {} session'.format(rec))
    # date_ref = pd.read_excel(work_dir + 'datetime_reference_miRNA.xlsx', index_col = 0)
    # if date_ref.at['MTA-{}'.format(mouse), 'short_end'] == 'yes':
    #     print('END exception')
    #     n = int((3*24 +12+4)*3600*sr) +1

    miRNA, genotype = get_mouse_info(mouse)
    #######         Extract scoring         #######
    print('get scoring mouse : ',mouse )
    filename = scoring_data_dir + '/{}/{}{}/MTA-{}/'.format(miRNA, genotype, miRNA,mouse)

    score_b0 = np.loadtxt(filename + mouse + 'b0.txt', dtype = str)
    score_b1 = np.loadtxt(filename + mouse + 'b1.txt', dtype = str)
    score_b2 = np.loadtxt(filename + mouse + 'b2.txt', dtype = str)
    score_sd = np.loadtxt(filename + mouse + 'sd.txt', dtype = str)
    score_r1 = np.loadtxt(filename + mouse + 'r.txt', dtype = str)
    print(score_b0.size, score_b1.size, score_b2.size, score_sd.size, score_r1.size)
    # all_score = np.zeros(score_b0.size + score_b1.size + score_b2.size + score_sd.size + score_r1.size, dtype = str)
    # one_day = score_b1.size
    all_score = np.concatenate((score_b0, score_b1, score_b2, score_sd, score_r1)).astype('str')
    # all_score[:one_day] = score_b1
    # all_score[:one_day] = score_b1
    # all_score[one_day:int(one_day*2)] = score_b2
    # all_score[int(2*one_day):int(3*one_day)] = score_sd
    # all_score[int(3*one_day):int(4*one_day)] = score_r1
    score_b0 = 0
    score_b1 = 0
    score_b2 = 0
    score_sd = 0
    score_r1 = 0
    time_epochs = np.arange(all_score.size)*4

    # exit()

    # exit()
    # print('get scoring somnologica spectrum : ',mouse )
    # file_b1 = open(data_dir + '/Power_spectrum/{}/MTA-{}/{}DCRb1s.txt'.format(group, mouse, mouse))
    # file_b2 = open(data_dir + '/Power_spectrum/{}/MTA-{}/{}DCRb2s.txt'.format(group, mouse, mouse))
    # file_sd = open(data_dir + '/Power_spectrum/{}/MTA-{}/{}DCRsds.txt'.format(group, mouse, mouse))
    # file_r1 = open(data_dir + '/Power_spectrum/{}/MTA-{}/{}DCRr1s.txt'.format(group, mouse, mouse))
    # file_r1 = open(data_dir + '/Power_spectrum/{}/MTA-{}/{}DCRr1s.txt'.format(group, mouse, mouse))
    # real_time_lines = []
    # data_lines = []
    # real_times_somno = []
    # data = []
    # for file in [file_b1, file_b2, file_sd, file_r1] :
    #     for i, line in enumerate(file.readlines()) :
    #         if line == 'Power Spectrum\n':
    #             real_time_lines.append(i+1)
    #             data_lines.append(i+2)
    #         if i in real_time_lines:
    #             real_times_somno.append(line[:-1])
    #         if i in data_lines :
    #             spec = np.array(line.split('\t')[:-1], dtype = float)[4:]    ####Remove .5 Hz     #######
    #             data.append(spec)
    # score_b1 = 0
    # score_b2 = 0
    # score_sd = 0
    # score_r1 = 0
    # #
    # freqs_somno = np.arange(4, spec.size+4, 1, dtype = 'float32') * .25
    # times_somno = np.arange(len(real_times_somno), dtype='int32')*int(4)
    #
    # somno_spectrums = np.array(data, dtype = 'float32')
    #
    # somno_spectrums = somno_spectrums[:,(freqs_somno >=.75) & (freqs_somno <= 47.5)]
    # freqs_somno = freqs_somno[(freqs_somno >=.75) & (freqs_somno <= 47.5)]
    #

    print('compute own spectrum mouse : ',mouse )
    ds = xr.open_dataset(edf_data_dir + '/precompute/raw/raw_{}.nc'.format(mouse))
    raw = ds['signal'].values.astype('float32')
    sr = ds['sampling_rate'].values
    #
    # times = ds.coords['times_second'].values
    times = np.arange(raw.size)/(sr*3600)

    t_start = times[0]
    n_epochs = int(raw.size//(windows*sr))

    point_per_epochs = int(4*sr)#800    ##### 4 sec at about 200 HZ   199,6
    ###### highpass filter 0.5 Hz
    N =3
    f_cut = .5
    nyq = sr/2
    W = f_cut/nyq
    b, a = scipy.signal.butter(N, W, btype = 'highpass', output = 'ba')
    raw = scipy.signal.filtfilt(b,a, raw)
    # print(raw.size)
    bandwidth = 1.01/windows
    epochs = np.arange(n_epochs)
    sample_per_epoch = int(windows*sr)
    mylist = []
    index = []
    for i in range(n_epochs):
        mylist.append(4*i*sr)
    real_sample_by_epoch= np.diff(np.array(mylist, dtype='int'))
    ind = 0
    for fr in real_sample_by_epoch:
        fr = int(fr)
        ind += fr
        if fr == 799:
            index.append(ind)
    index = np.array(index, dtype ='int')+1
    # ref = 69120000
    ref = int(time_epochs.size*800)
    # print(time_epochs.size)
    # print(ref, '  ref')
    if raw.size + index.size != ref:
        index = index[:-(raw.size + index.size-ref)]
    raw = np.insert(raw, index, raw[index])
    point_per_epochs = 800
    stacked_sigs = raw.reshape((-1, point_per_epochs)).astype('float32')

    # print(stacked_sigs.dtype)
    #
    freqs_welch, welch = scipy.signal.welch(stacked_sigs, fs = sr, nperseg = int(3.99*sr) )
    welch = welch[:,(freqs_welch>=55) & (freqs_welch<=80)]
    freqs_welch = freqs_welch[(freqs_welch>55) & (freqs_welch<=80)]

    # fft = scipy.fft.fft(stacked_sigs)
    # Fourier = np.fft.fft(stacked_sigs)
    # Fourier = np.fft.fft(np.random.rand(welch.shape))
    # Fourier = np.fft.fft(np.random.rand(10,10))
    # Fourier = F(stacked_sigs)
    # print(Fourier.shape)
    # Fourier = np.abs(Fourier)**2
    # # Fourier =  np.abs(np.fft.fft(stacked_sigs))**2
    # freqs_Fourier = np.fft.fftfreq(point_per_epochs, 1/sr)
    # idx = np.argsort(freqs_Fourier)
    # freqs_Fourier = freqs_Fourier[idx]
    # Fourier = Fourier[:,idx]
    # Fourier = Fourier[:,(freqs_Fourier>=.75)&(freqs_Fourier<=47.5)]
    # freqs_Fourier = freqs_Fourier[(freqs_Fourier>=.75)&(freqs_Fourier<=47.5)]

    # multitaper, freqs_multitaper = mne.time_frequency.psd_array_multitaper(stacked_sigs, sfreq = sr, fmax = 47.5, bandwidth = 4*bandwidth, n_jobs = 2)
    # multitaper = multitaper[:,freqs_multitaper>=.75]
    # freqs_multitaper = freqs_multitaper[freqs_multitaper>=.75]

    freqs_welch  = freqs_welch.astype('float32')
    # freqs_Fourier  = freqs_Fourier.astype('float32')
    # freqs_multitaper  = freqs_multitaper.astype('float32')
    epochs = epochs.astype('int32')
    welch = welch.astype('float32')
    # Fourier = Fourier.astype('float32')
    # multitaper = multitaper.astype('float32')
    # exit()

    ######      Caution ! Do not remove 1, 2, 3. It correspond to invalid EEG
    # for f, n in zip(['1', '2', '3'], ['w', 'n', 'r']) :
    #     score = np.where(score == f, n, score)


    coords = {  'freqs_welch': freqs_welch,
                # 'freqs_fft': freqs_fft,
                # 'freqs_multitaper': freqs_multitaper,
                'epochs' : epochs,
                'time_epochs' : time_epochs,
                # 'real_times_somno' : real_times_somno,
                # 'freqs_somno' : freqs_somno
                }

    ds = xr.Dataset(coords = coords)
    # ds['somno_spectrum'] = xr.DataArray(somno_spectrums, dims = ['epochs', 'freqs_somno'])
    ds['welch_spectrum'] = xr.DataArray(welch, dims = ['epochs', 'freqs_welch'])
    # ds['fft_spectrum'] = xr.DataArray(fft, dims = ['epochs', 'freqs_fft'])
    # ds['multitaper_spectrum'] = xr.DataArray(multitaper, dims = ['epochs', 'freqs_multitaper'])
    ds['sampling_rate'] = sr

    ds['score'] = xr.DataArray(all_score, dims = 'epochs')
    dirname = precompute_dir + '/spectrums/'
    if not os.path.exists(dirname):
        os.makedirs(dirname)
    print(dirname)
    ds.to_netcdf(dirname + 'HF_spectrum_scoring_{}.nc'.format(mouse))

# def F(x):
#     fft = np.fft.fft(x)
#     return fft
# def store_all_score_and_spectrum():
#     dcr_mice = get_mice(group = 'DCR-HCRT')
#     control_mice = get_mice(group = 'Control')
#     mice = dcr_mice+control_mice
#     results = Parallel(n_jobs=2)(delayed(store_scoring_and_spectrums_one_mouse_one_session)(mouse) for mouse in mice)

def HF_store_all_score_and_spectrum(parallel = True):
    for miRNA in ['128', '137', '665']:
        for genotype in ['test', 'control']:
            mice = get_mice(miRNA, genotype)
            if parallel:
                results = Parallel(n_jobs=3)(delayed(HF_store_scoring_and_spectrums_one_mouse)(mouse) for mouse in mice)
            else :
                for mouse in mice :
                    try :
                        HF_store_scoring_and_spectrums_one_mouse(mouse)
                    except :
                        print('-----------   Trouble with {}'.format(mouse))




if __name__ == '__main__':
    # mouse = 'B5626'  #format errors
    # mouse = 'B2884'
    mouse = 'B4627'
    # mouse = 'B5619'
    # HF_store_scoring_and_spectrums_one_mouse(mouse)
    # exit()
    # compute_own_spectrum(mouse)
    store_scoring_and_spectrums_one_mouse(mouse)
    # store_scoring_and_spectrums_one_mouse_one_session('Control', mouse)
    # get_all_raw_datas()
    # store_all_score_and_spectrum(False)
    theta_dominated_wake_encoding_one_mouse(mouse)
    # theta_dominated_wake_encoding_one_mouse('B2576')
    # theta_dominated_wake_encoding_all_mouse(parallel = False)
