from configuration import *
from select_mice_cata_Malo import get_mice

def average_baseline(spectrum_method = 'somno'):
    for group in ['Control', 'DCR-HCRT']:
        dirname = excel_dir + '{}/power_spectrum'.format(group)
        if group == 'Control':
            states = ['w', 'n', 'r']
        elif group == 'DCR-HCRT':
            states = ['w', 'n', 'r', 'a']
        for state in states :
            for period in ['light', 'dark']:
                df1 = pd.read_excel(dirname+'/spectrum_{}_{}_bl1_{}.xlsx'.format(spectrum_method, state, period), index_col =0)
                df2 = pd.read_excel(dirname+'/spectrum_{}_{}_bl2_{}.xlsx'.format(spectrum_method, state, period), index_col =0)
                mice1 = df1.index
                mice2 = df2.index
                freqs = df1.columns[4:].to_list()
                df = pd.DataFrame(index = get_mice(group), columns = ['ref_values'] + freqs )
                for ind in range(mice1.size) :
                    mouse1 = mice1[ind]
                    mouse2 = mice2[ind]
                    mouse = mouse1[:5]
                    df.at[mouse, 'ref_values'] = (df1.at[mouse1, 'ref_values']+ df2.at[mouse2, 'ref_values'])/2
                    transient = np.vstack((df1.loc[mouse1, freqs].values, df2.loc[mouse2, freqs].values))
                    df.loc[mouse,freqs] = np.nanmean(transient, axis =0)
                df.to_excel(dirname + '/spectrum_{}_{}_{}_baseline_{}.xlsx'.format(group,spectrum_method, state, period))

def average_recovery(spectrum_method = 'somno'):
    for group in ['Control', 'DCR-HCRT']:
        dirname = excel_dir + '{}/power_spectrum'.format(group)
        if group == 'Control':
            states = ['w', 'n', 'r']
        elif group == 'DCR-HCRT':
            states = ['w', 'n', 'r', 'a']
        for state in states :
            for period in ['light', 'dark']:
                df1 = pd.read_excel(dirname+'/spectrum_{}_{}_sd_{}.xlsx'.format(spectrum_method, state, period), index_col =0)
                df2 = pd.read_excel(dirname+'/spectrum_{}_{}_r1_{}.xlsx'.format(spectrum_method, state, period), index_col =0)
                mice1 = df1.index
                mice2 = df2.index
                freqs = df1.columns[4:].to_list()
                df = pd.DataFrame(index = get_mice(group), columns = ['ref_values'] + freqs )
                for ind in range(mice1.size) :
                    mouse1 = mice1[ind]
                    mouse2 = mice2[ind]
                    mouse = mouse1[:5]
                    df.at[mouse, 'ref_values'] = (df1.at[mouse1, 'ref_values']+ df2.at[mouse2, 'ref_values'])/2
                    transient = np.vstack((df1.loc[mouse1, freqs].values, df2.loc[mouse2, freqs].values))
                    df.loc[mouse,freqs] = np.nanmean(transient, axis =0)
                df.to_excel(dirname + '/spectrum_{}_{}_{}_recovery_{}.xlsx'.format(group, spectrum_method, state, period))





if __name__ == '__main__':
    average_baseline()
    average_recovery()
