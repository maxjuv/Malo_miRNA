
from configuration import *
import scipy.signal
import mne.time_frequency
import multiprocessing
from joblib import Parallel, delayed

windows = 4   #### in s


def sliding_power_spectrum(mouse):
    ds = xr.open_dataset(precompute_dir + '/raw/raw_{}.nc'.format(mouse))
    raw = ds['signal'].values
    times = ds.coords['times_second'].values
    t_start = times[0]
    sr = ds['sampling_rate'].values
    n_epochs = ((times[-1] - times[0]))//windows

    ###### highpass filter 0.3 Hz
    N =3
    f_cut = .3
    nyq = sr/2
    W = f_cut/nyq
    b, a = scipy.signal.butter(N, W, btype = 'highpass', output = 'ba')
    raw = scipy.signal.filtfilt(b,a, raw)



    bandwidth = 1.01/windows
    epochs = np.arange(n_epochs)
    for w in epochs:
        mask = (times>t_start + windows*w) & (times<=t_start + windows*(w+1))
        sig = raw[mask]
        sig = sig - np.mean(sig)

        #
        # ##### Mutltitaper
        # # psds, freqs_mnt = mne.time_frequency.psd_array_multitaper(sig, sfreq = sr, fmax = 10, bandwidth = bandwidth)
        # # psds = psds[freqs_mnt>.4]
        # # freqs_mnt = freqs_mnt[freqs_mnt>.4]
        #
        #
        freqs_welch, pxx = scipy.signal.welch(sig, fs = sr, nperseg = int(3.99*sr) )
        pxx = pxx[(freqs_welch>.4) & (freqs_welch<=10)]
        freqs_welch = freqs_welch[(freqs_welch>.4) & (freqs_welch<=10)]
        if w == 0:
            spectrums = np.array(pxx)
        else :
            spectrums = np.vstack((spectrums, np.array(pxx)))
        #
        # if w%1000 == 0:
        #     print('{:.1f} % achieved'.format(w/n_epochs))
        #     print('{} on {}'.format(w, n_epochs))

        # fig, ax = plt.subplots(nrows =3)
        # ax[0].plot(freqs_mnt, psds)
        # ax[0].set_ylabel('multitaper')
        # ax[1].plot(freqs_welch, pxx)
        # ax[1].set_ylabel('welch')
        # ax[2].plot(np.arange(sig.size)/sr, sig)
        # ax[2].set_ylabel('raw')
        # plt.show()
        # exit()

    # print(freqs_welch.shape)
    # print(spectrums.shape)
    coords = {'freqs': freqs_welch, 'epochs' : epochs}
    ds_spectrum = xr.Dataset(coords = coords)
    ds_spectrum['spectrum_by_epoch'] = xr.DataArray(spectrums, dims = ['epochs', 'freqs'])

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
        sliding_power_spectrum(mouse)

def hack_dask():
    run_key = sys.argv[1]
    print(run_key, type(run_key))
    sliding_power_spectrum(run_key)


if __name__ == '__main__':
    # mouse = 'B2533'
    hack_dask()
    # sliding_power_spectrum(mouse)
    # read_one_mouse_spectrums(mouse)
    # compute_all_spectrums_data()
