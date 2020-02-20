
from configuration import*


def clusterize():
    date_ref = pd.read_excel(work_dir + 'datetime_reference_DICER.xls', index_col = 0)
    mice = date_ref.index.to_list()

    nb_core = 1
    jobname = 'malo'
    memory = 30
    # function = 'edf_to_dataset.py'
    function = 'power_spectrum.py'

    error_path = precompute_dir + '/errors/'

    for mouse in mice:
        mouse = mouse[4:]
        print('srun --cpus-per-task='+str(nb_core)+
                    ' --job-name=' + str(jobname)+
                    ' --mem=' + str(memory) +
                    'G --error=' + str(error_path) + str(mouse)+ '.txt'+
                    'G python ' + function +
                    ' '+ str(mouse)+'&')
        os.system('srun --cpus-per-task='+str(nb_core)+
                    ' --job-name=' + str(jobname)+
                    ' --mem=' + str(memory) +
                    'G --error=' + str(error_path) + str(mouse)+ '.txt'+
                    'G python ' + function +
                    ' '+ str(mouse)+'&')



if __name__ == '__main__':
    clusterize()
