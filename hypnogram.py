from configuration import *
from select_mice_cata_Malo import get_mice,get_mouse_info

def plot_hypnogram_one_mouse(mouse, time_by_line = 12):
    miRNA, genotype = get_mouse_info(mouse)
    print(mouse)
    ds = xr.open_dataset(precompute_dir + '/spectrums/spectrum_scoring_{}.nc'.format(mouse))
    times = ds.coords['time_epochs'].values/3600
    score = ds['score'].values
    score_behavior = score.copy()
    # group_color = {'Control':'black', 'DCR-HCRT':'seagreen'}
    group_color = {'control':'black', 'test':'black'}
    for f, n in zip(['1', '2', '3'], ['w', 'n', 'r']) :
        score_behavior = np.where(score_behavior == f, n, score_behavior)
    print('Remaining possible score : ', np.unique(score_behavior))
    height_r = 1
    height_n = 2
    height_w = 3
    height_a = 4

    a = (score_behavior == 's')*height_a
    r = (score_behavior == 'r')*height_r
    n = (score_behavior == 'n')*height_n
    w = (score_behavior == 'w')*height_w

    hypno = r +n + w + a
    fig, ax = plt.subplots(figsize=(20,10))
    h = time_by_line
    loop = int(96/h)
    ax.plot(times[times<4]+8, hypno[times<4], color = 'k')
    hypno = hypno[times>=4]
    times = times[times>=4]
    for i in range(loop+1):
        t1, t2 = h*i+4, h*(i+1)+4
        mask = (times > t1) & (times<t2)
        ax.plot(times[mask] - h*i -4, 3.5*(i+1) + hypno[mask], color = 'k')
    # ax.set_xlim(0,4)
    fig.suptitle('{} -- r = {}, n = {}, w = {}'.format(mouse, height_r, height_n, height_w) )
    return fig
def save_all_hypnogram(time_by_line = 12):
    for miRNA in ['128', '137', '665']:
        for genotype in ['test', 'control']:
            mice = get_mice(miRNA, genotype)
            dirname = work_dir+'/pyFig/hypnogram/{}/{}{}/'.format(miRNA, genotype, miRNA)
            if not os.path.exists(dirname):
                os.makedirs(dirname)
            for mouse in mice :
                try :
                    fig = plot_hypnogram_one_mouse(mouse, time_by_line = time_by_line)
                    plt.savefig(dirname + mouse+'.png')
                except:
                    print(mouse)
if __name__ == '__main__':

    mouse = 'B3868'
    # mouse = 'B2761'
    # plot_hypnogram_one_mouse(mouse,time_by_line = 12)
    # plt.show()
    save_all_hypnogram()
