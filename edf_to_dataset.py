# !/Users/maximejuventin/anaconda3/bin/python

from configuration import *

import pyedflib


def edf_to_dataset(mouse):
    print('Extract raw data from mouse {}'.format(mouse))

    ##### Read excel baseline 1 begin
    date_ref = pd.read_excel(work_dir + 'datetime_reference_DICER.xls', index_col = 0)
    d = date_ref.at['MTA-{}'.format(mouse), 'baseline1']
    true_EEG = date_ref.at['MTA-{}'.format(mouse), 'channel']
    begin_date = datetime.datetime(d.year, d.month, d.day, 8,0 ,0)

    ##### Read EDF
    filename = data_dir + '{}.edf'.format(mouse)
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
    sr = f.samplefrequency(channel_num)
    print(channel_labels)

    ####### Slice signal from Begin to the end
    start_rec_date = f.getStartdatetime()
    pre_rec_duration = begin_date - start_rec_date
    index_start = int(pre_rec_duration.total_seconds()*sr)
    d, h, m, s = pre_rec_duration.days, pre_rec_duration.seconds//3600, (pre_rec_duration.seconds//60)%60, (pre_rec_duration.seconds%60)
    print('Recording started {}\nBaseline1 started {}\nDuration in between is {}'.format(start_rec_date, begin_date, pre_rec_duration))
    n = int(4*24*3600*sr)
    times = np.arange(f.getNSamples())/sr
    sliced_signal = f.readSignal(chn = channel_num, start =  index_start, n = n)
    sliced_time_in_hours =times[index_start:index_start+n]/3600 -(24*d + h + m/60 + s/3600) +8
    sliced_time_in_seconds = sliced_time_in_hours*3600

    print(times.size)
    print(sliced_signal.size)
    print(3600*sliced_signal.size/sr)
    ##### Store to dataset and save
    # coords = {'times_second' : sliced_time_in_seconds, 'times_hour' : sliced_time_in_hours}
    coords = {'times_second' : sliced_time_in_seconds}
    ds = xr.Dataset(coords = coords)
    # ds = xr.Dataset()

    ds['signal'] = xr.DataArray(sliced_signal, dims = 'times_second')
    # ds['signal'] = xr.DataArray(sliced_signal)
    ds['sampling_rate'] = sr

    saving_path = precompute_dir +'/raw/'
    if not os.path.exists(saving_path):
        os.makedirs(saving_path)
    print(saving_path)
    ds.to_netcdf(saving_path + 'raw_{}.nc'.format(mouse))


def read_data_one_mouse(mouse):
    ds = xr.open_dataset(precompute_dir + 'raw_{}.nc'.format(mouse))
    print(ds)
    print(ds['signal'].values)

def get_all_raw_data():
    date_ref = pd.read_excel(work_dir + 'datetime_reference_DICER.xls', index_col = 0)
    mice = date_ref.index.to_list()
    for mouse in mice:
        mouse = mouse[4:]
        edf_to_dataset(mouse)

def hack_dask():
    run_key = sys.argv[1]
    print(run_key, type(run_key))
    edf_to_dataset(run_key)

if __name__ == '__main__':
    # mouse = 'B2533'
    mouse = 'B4907'
    # mouse = 'B4977'
    # edf_to_dataset(mouse)
    hack_dask()
    # read_data_one_mouse(mouse)
    # get_all_raw_data()
