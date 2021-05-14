# !/Users/maximejuventin/anaconda3/bin/python

from configuration import *
from select_mice_cata_Malo import get_mice, get_mouse_info
import multiprocessing
from joblib import Parallel, delayed
import pyedflib


def edf_to_dataset(mouse):
    print('Extract raw data from mouse {}'.format(mouse))

    ##### Read excel baseline 1 begin
    date_ref = pd.read_excel(work_dir + 'datetime_reference_miRNA.xlsx', index_col = 0)
    d = date_ref.at['MTA-{}'.format(mouse), 'baseline1']
    true_EEG = date_ref.at['MTA-{}'.format(mouse), 'channel']
    begin_date = datetime.datetime(d.year, d.month, d.day, 4,0 ,0) ####4h before light 8-4
    sr = date_ref.at['MTA-{}'.format(mouse), 'sampling_rate']
    miRNA, genotype = get_mouse_info(mouse)


    ##### Read EDF
    # filename = edf_data_dir + '/{}/{}{}/{}.edf'.format(miRNA, genotype, miRNA,mouse)
    filename = work_dir + '/data/MiRNA_edf/{}/{}{}/{}.edf'.format(miRNA, genotype, miRNA,mouse)
    print(filename)
    # exit()
    f = pyedflib.EdfReader(filename)
    n = f.signals_in_file
    if n >1:
        print('ERROR, more than 1 EEG in EDF file')
        exit()
    channel_labels = f.getSignalLabels()

    if channel_labels[0] != true_EEG:
        print('ERROR EEG channel does not correspond between excel and EDF')
        print('EDF : {} vs. Excel : {}'.format(channel_labels[0],true_EEG))
        exit()
    channel_num = 0


    # sr = f.getSampleFrequency(channel_num)
    # print('real sr = {:.5f}'.format(f.getNSamples()[0]/f.getFileDuration()))
    # sr = f.getNSamples()[0]/f.getFileDuration()
    # print(sr)
    # print(f.getSampleFrequencies())
    # print(np.float(f.samplefrequency(channel_num)))
    # sr = 199.985

    # exit()
    # sr = f.samplefrequency(channel_num)
    # DO NOT USE f.samplefrequency(channel_num))
    print(channel_labels)

    ####### Slice signal from Begin to the end
    start_rec_date = f.getStartdatetime()
    pre_rec_duration = begin_date - start_rec_date
    index_start = int(pre_rec_duration.total_seconds()*sr)
    d, h, m, s = pre_rec_duration.days, pre_rec_duration.seconds//3600, (pre_rec_duration.seconds//60)%60, (pre_rec_duration.seconds%60)
    print('Recording started {}\nBaseline1 started {}\nDuration in between is {}'.format(start_rec_date, begin_date, pre_rec_duration))
    print(date_ref.at['MTA-{}'.format(mouse), 'short_end'])
    # exit()
    if date_ref.at['MTA-{}'.format(mouse), 'short_end'] == 'yes':
        print('END exception')
        n = int((3*24 +12+4)*3600*sr) +1
    else :
        n = int((4*24+4)*3600*sr) +1
    times = np.arange(f.getNSamples()[0])/sr
    # times = times.astype('float32')
    ###########sr somnologica = 199.985
    start_sleep_lab = start_rec_date + datetime.timedelta(seconds=int(pre_rec_duration.total_seconds()*sr)/sr)
    stop_sleep_lab = start_sleep_lab + datetime.timedelta(seconds = n/sr)
    print(start_sleep_lab, '----', stop_sleep_lab)
    # exit()
    # stop_rec_date = start_rec_date + datetime.timedelta(seconds = f.getNSamples()[0]/sr)
    # real_times = np.arange(start_rec_date, stop_rec_date, datetime.timedelta(seconds = 1/sr))
    # print(real_times[index_start])
    # print(real_times.shape, times.shape)
    #
    # print(1/sr, (real_times[20]-real_times[10])/10)
    # exit()
    sliced_signal = f.readSignal(chn = channel_num, start =  index_start, n = n)
    sliced_time_in_hours =times[index_start:index_start+n]/3600 -(24*d + h + m/60 + s/3600) +8
    sliced_time_in_seconds = sliced_time_in_hours*3600
    sliced_signal = sliced_signal.astype('float32')

    # print(real_times[index_start], real_times[index_start+n])
    # print(times.size)
    # print(sliced_signal.size)
    # print(3600*sliced_signal.size/sr)
    ##### Store to dataset and save
    # coords = {'times_second' : sliced_time_in_seconds, 'times_hour' : sliced_time_in_hours}

    # coords = {'times_second' : sliced_time_in_seconds}
    # ds = xr.Dataset(coords = coords)
    ds = xr.Dataset()

    # ds['signal'] = xr.DataArray(sliced_signal, dims = 'times_second')
    ds['signal'] = xr.DataArray(sliced_signal)
    ds['sampling_rate'] = sr

    # exit()
    # saving_path = precompute_dir +'/raw/'
    saving_path = precompute_dir +'/raw/'
    # saving_path = edf_data_dir +'/precompute/raw/'
    if not os.path.exists(saving_path):
        os.makedirs(saving_path)
    print(saving_path)
    ds.to_netcdf(saving_path + 'raw_{}.nc'.format(mouse))


def read_data_one_mouse(mouse):
    ds = xr.open_dataset(precompute_dir + 'raw_{}.nc'.format(mouse))
    print(ds)
    print(ds['signal'].values)

def get_all_raw_data(parallel = True):
    for miRNA in ['128', '137', '665']:
        for genotype in ['test', 'control']:
            mice = get_mice(miRNA, genotype)
            for mouse in mice :
                if parallel :
                    results = Parallel(n_jobs=2)(delayed(edf_to_dataset)(mouse) for mouse in mice)
                try :
                    edf_to_dataset(mouse)
                except :
                    print('-----------   Trouble with {}'.format(mouse))

def hack_dask():
    run_key = sys.argv[1]
    print(run_key, type(run_key))
    edf_to_dataset(run_key)

if __name__ == '__main__':
    mouse = 'B2884'  #format errors
    # mouse = 'B2882' #format errors
    # mouse = 'B2880' #format errors
    # mouse = 'B2878' #format errors
    # mouse = 'B2879' #format errors
    # mouse = 'B2877' #format errors


    mouse = 'B4627'
    # mouse = 'B2574'
    # mouse = 'B5627'
    edf_to_dataset(mouse)
    # get_all_raw_data(False)
