from configuration import *
from select_mice_cata_Malo import get_mice

def plot_hypnogram_one_mouse(mouse, time_by_line = 12):
    ds = xr.open_dataset(precompute_dir + '/spectrums/spectrum_scoring_{}.nc'.format(mouse))
    times = ds.coords['times_somno'].values/3600
    score = ds['score'].values
    score_behavior = score.copy()
    for f, n in zip(['1', '2', '3'], ['w', 'n', 'r']) :
        score_behavior = np.where(score_behavior == f, n, score_behavior)
    print('Remaining possible score : ', np.unique(score_behavior))
    height_a = 3
    height_r = 1
    height_n = 2
    height_w = 3

    a = (score_behavior == 'a')*height_a
    r = (score_behavior == 'r')*height_r
    n = (score_behavior == 'n')*height_n
    w = (score_behavior == 'w')*height_w

    hypno = a + r +n + w
    fig, ax = plt.subplots(figsize=(20,10))
    h = time_by_line
    loop = int(96/h)
    for i in range(loop+1):
        t1, t2 = h*i, h*(i+1)
        mask = (times > t1) & (times<t2)
        ax.plot(times[mask] - h*i, 3.5*i + hypno[mask], color = 'black')
    # ax.set_xlim(0,4)
    fig.suptitle('{} -- a = {}, r = {}, n = {}, w = {}'.format(mouse, height_a, height_r, height_n, height_w) )
    return fig
def save_all_hypnogram(time_by_line = 12):
    control_list = get_mice('Control')
    DCR_list = get_mice('DCR-HCRT')
    groups = {'Control' : control_list, 'DCR-HCRT' : DCR_list}
    for group in groups:
        dirname = work_dir+'/pyFig/{}/hypnogram/'.format(group)
        if not os.path.exists(dirname):
            os.makedirs(dirname)
        for mouse in groups[group]:
            fig = plot_hypnogram_one_mouse(mouse, time_by_line = time_by_line)
            plt.savefig(dirname + mouse+'.png')
if __name__ == '__main__':

    # mouse = 'B2700'
    # mouse = 'B2761'
    # plot_hypnogram_one_mouse(time_by_line = 12)
    save_all_hypnogram()
