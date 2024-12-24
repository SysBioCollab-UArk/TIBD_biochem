import pandas as pd
import numpy as np
import os
import csv
from io import StringIO
import re
import itertools
import matplotlib.pyplot as plt
from matplotlib.lines import Line2D
from warnings import warn

def dixon_q_test(data, alpha=0.05):
    """
    Perform Dixon's Q Test for outliers in a small dataset.

    Parameters:
        data (list): List of numerical data points (3 <= len(data) <= 10).
        alpha (float): Significance level (default: 0.05).

    Returns:
        tuple: (is_outlier, q_stat, q_critical, suspected_value)
            is_outlier (bool): Whether an outlier was detected.
            q_stat (float): The calculated Q statistic.
            q_critical (float): The critical Q value for the test.
            suspected_value (float): The suspected outlier value.
    """
    # Predefined critical values for n=3 to n=10
    critical_values = {
        3: {0.10: 0.941, 0.05: 0.970, 0.01: 0.994},
        4: {0.10: 0.679, 0.05: 0.765, 0.01: 0.880},
        5: {0.10: 0.564, 0.05: 0.624, 0.01: 0.698},
        6: {0.10: 0.486, 0.05: 0.560, 0.01: 0.625},
        7: {0.10: 0.434, 0.05: 0.507, 0.01: 0.568},
        8: {0.10: 0.391, 0.05: 0.468, 0.01: 0.526},
        9: {0.10: 0.359, 0.05: 0.437, 0.01: 0.493},
        10: {0.10: 0.336, 0.05: 0.412, 0.01: 0.466},
    }

    # Validate input
    n = len(data)
    if n < 3 or n > 10:
        raise ValueError("Dixon's Q Test is only valid for 3 <= n <= 10.")

    if alpha not in [0.10, 0.05, 0.01]:
        raise ValueError("Supported alpha levels are 0.10, 0.05, and 0.01.")

    # Sort data while keeping track of the original indices
    data_sorted_with_indices = sorted(enumerate(data), key=lambda x: x[1])
    data_sorted = [x[1] for x in data_sorted_with_indices]
    original_indices = [x[0] for x in data_sorted_with_indices]

    # Calculate Q statistic for smallest and largest values
    q_stat_low = (data_sorted[1] - data_sorted[0]) / (data_sorted[-1] - data_sorted[0])
    q_stat_high = (data_sorted[-1] - data_sorted[-2]) / (data_sorted[-1] - data_sorted[0])

    # Determine the most extreme value
    q_stat = max(q_stat_low, q_stat_high)
    outlier_value = data_sorted[0] if q_stat == q_stat_low else data_sorted[-1]
    outlier_idx = original_indices[0] if q_stat == q_stat_low else original_indices[-1]

    # Get critical Q value
    q_critical = critical_values[n][alpha]

    # Determine if outlier
    is_outlier = q_stat > q_critical

    print('==========')
    print('data:', data)
    print('suspected_value:', outlier_value)
    print('idx:', outlier_idx)
    print('==========')

    return is_outlier, outlier_idx, outlier_value, q_stat, q_critical

# # Example usage
# data = [25.1, 24.9, 30.5, 25.3]  # Example dataset
# alpha = 0.05
# is_outlier, q_stat, q_critical, suspected_value = dixon_q_test(data, alpha=alpha)
#
# print(f"Outlier detected: {is_outlier}")
# print(f"Q Statistic: {q_stat:.3f}, Q Critical: {q_critical:.3f}")
# print(f"Suspected Outlier Value: {suspected_value}")
# quit()

def display_platemap(pmap):
    # get platemap rows
    rows = np.unique([key[0] for key in pmap.keys()])
    cols = np.unique([key[1:] for key in pmap.keys()])
    for row in rows:
        wells = ['%s%s' % (row, col) for col in cols]
        print(['%s: (%s, %s, %s, %s, %s)' %
               (well, pmap[well][0], pmap[well][1], pmap[well][2], pmap[well][3], pmap[well][4]) for well in wells])


def read_platemap_data(dirname):

    # STEP 1: Read in the platemap info
    file_path = os.path.join(dirname, 'platemap.xlsx')
    platemap_xlsx = pd.read_excel(file_path, header=None)
    platemap_csv = platemap_xlsx.to_csv(header=False, index=False)

    # make the platemap a dictionary with keys corresponding to the different platemap files
    platemap = {}

    reader = csv.reader(StringIO(platemap_csv), delimiter=',')
    for line in reader:
        '''print(line)'''
        if re.search('^Platemap:', line[0]):
            pmap_name = re.search('^Platemap:\s*(.+)$', line[0]).group(1).strip()
            platemap[pmap_name] = {}  # make each individual platemap a dict with keys corresponding to wells
            while line[0] != 'Col_Row':
                line = next(reader)
            cols = [n for n in line]  # column labels (integers)
            line = next(reader)
            print(line)
            while len(line[0]) == 1:  # row label is a letter, e.g., 'A', 'B', etc.
                row = line[0]
                for i in range(1, len(line)):
                    well_id = '%s%s' % (row, cols[i] if len(cols[i]) == 2 else '0%s' % cols[i])
                    sample = line[i] if line[i] == '' else int(line[i])
                    platemap[pmap_name][well_id] = [sample, '', '', '', '']  # sample, target, Cq, cell_line, substrate
                try:
                    line = next(reader)
                except StopIteration:  # in case we reach the end of the file
                    break
                print(line)
                for i in range(1, len(line)):
                    well_id = '%s%s' % (row, cols[i] if len(cols[i]) == 2 else '0%s' % cols[i])
                    target = line[i]
                    print(well_id, target)
                    platemap[pmap_name][well_id][1] = target
                try:
                    line = next(reader)
                except StopIteration:  # in case we reach the end of the file
                    break
                print(line)

    # STEP 2: Read in the Cq values
    for pmap_name in platemap.keys():

        file_path = os.path.join(dirname, 'RawData', '%s -  Quantification Summary.xlsx' % pmap_name)
        qpcr_data_xlsx = pd.read_excel(file_path, sheet_name=0, header=None)
        qpcr_data_csv = qpcr_data_xlsx.to_csv(header=False, index=False)

        reader = csv.reader(StringIO(qpcr_data_csv), delimiter=',')
        next(reader)  # read the first line
        for line in reader:
            well_id = line[1]
            Cq = line[6] if line[6] == '' else float(line[6])
            platemap[pmap_name][well_id][2] = Cq

    # STEP 3: Read in cell line and substrate information for each sample
    file_path = os.path.join(dirname, 'cell_culture_samples.xlsx')
    samples_xlsx = pd.read_excel(file_path, header=None)
    samples_csv = samples_xlsx.to_csv(header=False, index=False)

    samples_data = np.genfromtxt(StringIO(samples_csv), dtype=None, delimiter=',', names=True, encoding="utf_8_sig")

    for pmap_name in platemap.keys():
        # make a dict mapping sample # to (cell line, substrate) tuple
        samples_dict = dict(
            zip(
                samples_data['Sample'],
                [(cell, substr) for cell, substr in zip(samples_data['CellLine'], samples_data['Substrate'])]
            )
        )
        # print(samples_dict)
        for well in platemap[pmap_name].keys():
            print(pmap_name)
            print(well)
            print(platemap[pmap_name][well])
            sample = platemap[pmap_name][well][0]
            print(sample)
            print('----------')
            if sample != '':  # make sure this is not an empty well
                cell_line = samples_dict[sample][0]
                substrate = samples_dict[sample][1]
                platemap[pmap_name][well][3:4] = [cell_line, substrate]

    return platemap


def extract_cell_substrate_data(platemap, control):

    # extract lists of targets, cell lines, and substrates from the platemaps
    targets = []
    cell_lines = []
    substrates = []
    for pmap_name in platemap.keys():
        pmap = platemap[pmap_name]
        targets += [pmap[well][1] for well in pmap.keys() if pmap[well][1] != control and pmap[well][1] != '']
        cell_lines += [pmap[well][3] for well in pmap.keys() if pmap[well][3] != '']
        substrates += [pmap[well][4] for well in pmap.keys() if pmap[well][4] != '']
    targets = np.unique(targets)
    cell_lines = np.unique(cell_lines)
    substrates = np.unique(substrates)
    print(targets)
    print(cell_lines)
    print(substrates)

    Cq_data = {}  # store qPCR measurements here

    # loop over (target, cell line, substrate) groups
    for group in itertools.product(targets, cell_lines, substrates):
        print(group)
        target = group[0]
        cell_line = group[1]
        substrate = group[2]
        ctrl_key = '%s_%s' % (control, target)
        qpcr_data_temp = {target: [], ctrl_key: []}
        # loop over platemaps
        for pmap_name in platemap.keys():
            pmap = platemap[pmap_name]
            # save values of (sample #, Cq) for the target gene
            samples_Cq = [(pmap[well][0], pmap[well][2]) for well in pmap.keys() if
                          pmap[well][2] != '' and  # missing data
                          pmap[well][1] == target and pmap[well][3] == cell_line and pmap[well][4] == substrate]
            # if found values, get the Cq values for the control gene
            if len(samples_Cq) > 0:
                qpcr_data_temp[target].append(samples_Cq)
                samples = np.unique([well[0] for well in samples_Cq])
                qpcr_data_temp[ctrl_key].append([(pmap[well][0], pmap[well][2]) for well in pmap.keys() if
                                                 pmap[well][2] != '' and  # missing data
                                                 pmap[well][0] in samples and pmap[well][1] == control and
                                                 pmap[well][3] == cell_line and pmap[well][4] == substrate])

        # now consolidate Cq values for the target and control genes based on the sample number
        samples = np.unique([well[0] for q in qpcr_data_temp[target] for well in q])
        Cq_target = [[] for sample in samples]
        Cq_ctrl = [[] for sample in samples]
        for i in range(len(qpcr_data_temp[target])):
            print('q_t:', qpcr_data_temp[target][i])
            print('q_c:', qpcr_data_temp[ctrl_key][i])
            print(Cq_target)
            print(Cq_ctrl)
            for j, sample in enumerate(samples):
                print('j: %d, sample: %d' % (j, sample))
                print('  %s_t:' % sample, [well[1] for well in qpcr_data_temp[target][i] if well[0] == sample])
                print('  %s_c:' % sample, [well[1] for well in qpcr_data_temp[ctrl_key][i] if well[0] == sample])
                print(Cq_target)
                print(Cq_ctrl)
                Cq_target[j].append([well[1] for well in qpcr_data_temp[target][i] if well[0] == sample])
                Cq_ctrl[j].append([well[1] for well in qpcr_data_temp[ctrl_key][i] if well[0] == sample])
                print(Cq_target)
                print(Cq_ctrl)
        for i in range(len(Cq_target)):
            print('pmap %d' % i)
            print('  %s:' % target, Cq_target[i])
            print('  %s:' % ctrl_key, Cq_ctrl[i])
        # if haven't seen this (cell line, substrate) pair before, add it as a key to qpcr_vals dict
        if (cell_line, substrate) not in Cq_data.keys():
            Cq_data[(cell_line, substrate)] = {}
        # loop over the samples identified for this (cell line, substrate) pair
        for i in range(len(samples)):
            # if haven't seen this sample number before, add it as a key to qpcr_vals[(cell_line, substrate)] dict
            if samples[i] not in Cq_data[(cell_line, substrate)].keys():
                Cq_data[(cell_line, substrate)][samples[i]] = {}
            # now, add the Cq values for the target and control genes to the existing dict
            Cq_data[(cell_line, substrate)][samples[i]][target] = [cq for cq in Cq_target[i] if len(cq) > 0]
            Cq_data[(cell_line, substrate)][samples[i]][ctrl_key] = [cq for cq in Cq_ctrl[i] if len(cq) > 0]

        print((cell_line, substrate))
        for key in Cq_data[(cell_line, substrate)].keys():
            print('%s:' % key, Cq_data[(cell_line, substrate)][key])

    return Cq_data


def plot_relative_mRNA(Cq_data, ref_gene, ctrl_sample, add_subplot=None):
    print()
    dCt = {}  # Delta Ct values
    # loop over (cell line, substrate) pairs
    for key in Cq_data.keys():
        print(key)
        print(Cq_data[key].keys())
        dCt[key] = {}
        # loop over sample numbers
        for key2 in Cq_data[key].keys():
            print('  %s' % key2)
            avg = {}
            targets = []
            # loop over targets + controls
            for key3 in Cq_data[key][key2].keys():
                print('    %s:' % key3, Cq_data[key][key2][key3], end=' ')
                avg[key3] = np.array([np.mean(cq) for cq in Cq_data[key][key2][key3]])
                print(avg[key3])
                # save targets
                if ref_gene not in key3:
                    targets.append(key3)
            print(targets)
            # loop over targets
            for target in targets:
                # if this is the first time we've seen this target, add it to the dct[key] dict
                if target not in dCt[key].keys():
                    dCt[key][target] = []
                # calculate Delta_Ct
                dCt[key][target].append(list(avg[target] - avg['%s_%s' % (ref_gene, target)]))

    print()
    for key in dCt.keys():
        print('%s:' % str(key), dCt[key])
    print()

    genes = []
    ddCt = {}  # Delta Delta Ct values
    print('dCt[ctrl_sample].keys():', list(dCt[ctrl_sample].keys()))
    for cell_substr in dCt.keys():
        print('%s:' % str(cell_substr), list(dCt[cell_substr].keys()))
        #####
        for gene in dCt[cell_substr].keys():
            print('  %s:' % gene, dCt[cell_substr][gene])
            print('  %s' % gene, [np.array(d) for d in dCt[cell_substr][gene]])
            # print('  %s_ctrl' % gene, dCt[ctrl_sample][gene])
            print('  %s_ctrl' % gene, np.array([dd for d in dCt[ctrl_sample][gene] for dd in d]))
            print('  %s_ctrl' % gene, np.mean([dd for d in dCt[ctrl_sample][gene] for dd in d]))
            print('  diff', [np.array(d) - np.mean([xx for x in dCt[ctrl_sample][gene] for xx in x])
                             for d in dCt[cell_substr][gene]])
            print()
        #####
        ddCt[cell_substr] = dict(
            zip(
                list(dCt[cell_substr].keys()),
                [[np.array(d) - np.mean([xx for x in dCt[ctrl_sample][gene] for xx in x])
                  for d in dCt[cell_substr][gene]] for gene in dCt[cell_substr].keys()
                 if gene in dCt[ctrl_sample].keys()]
            )
        )
        genes += [gene for gene in dCt[cell_substr].keys() if gene not in genes]
        print(ddCt[cell_substr])
    genes.sort(key=str.lower)  # case-insensitive sorting
    print(genes)
    print()

    if add_subplot is None:
        fig, ax = plt.subplots(constrained_layout=True)
    else:
        fig = add_subplot[0]
        fig_axes = fig.axes
        sharey = None
        if add_subplot[2] is True and len(fig_axes) > 0:
            sharey = fig_axes[-1] # share y-axes between all plots
        ax = fig.add_subplot(add_subplot[1], sharey=sharey)
    xtick_subs = {'Bone Clone': 'Bone', 'Parental': 'Par', 'Plastic': 'TC', 'LC Aligned': 'AC', 'LC Random': 'RC'}
    cycle = plt.rcParams['axes.prop_cycle'].by_key()['color']  # standard colors
    width = 0.75 / len(genes)
    markers = []
    labels = []
    for i, key in enumerate(ddCt.keys()):
        print(key)
        for key2 in ddCt[key].keys():
            print('  %s:' % key2, ddCt[key][key2])
            j = genes.index(key2)
            print('  j: %d' % j)
            data = [dd for d in ddCt[key][key2] for dd in d]
            print('    data:', data)
            print('    2^-data:', 2 ** -np.array(data))
            print('    2^-mean(data):', 2 ** -np.mean(data))
            # #######################################
            # check for outliers using Dixon's Q Test
            if len(data) >= 3:
                alpha = 0.05
                is_outlier, outlier_idx, outlier_value, q_stat, q_critical = \
                    dixon_q_test(2 ** -np.array(data), alpha=alpha)
                print(f"Outlier detected: {is_outlier}")
                print(f"Q Statistic: {q_stat:.3f}, Q Critical: {q_critical:.3f}")
                print(f"Suspected Outlier Value: {outlier_value}")
                if is_outlier:
                    warn(f"Outlier detected: {key} : {key2}\nValue: {outlier_value}\n"
                         f"Q Statistic: {q_stat:.3f}\nQ Critical: {q_critical:.3f}")
                    # remove outlier from the flattened data array
                    data.pop(outlier_idx)
                    print('    data:', data)
                    print('    2^-data:', 2 ** -np.array(data))
                    print('    2^-mean(data):', 2 ** -np.mean(data))
                    # set outlier value to NaN in the structured ddCt array
                    for d in ddCt[key][key2]:
                        if outlier_idx < len(d):
                            d[outlier_idx] = None
                            break
                        else:
                            outlier_idx -= len(d)
                    print('   ', ddCt[key][key2])
            # #######################################
            label = None
            if key2 not in labels:  # need to do this to make sure we get all the gene labels
                label = key2
                labels.append(label)
            offset = width * (j - 1)
            ax.bar(i + offset, 2 ** -np.mean(data), width=0.8 * width, color=cycle[j % 10], label=label)
            for k, d in enumerate(ddCt[key][key2]):
                markers.append(Line2D.filled_markers[1:][k])
                ax.plot([i + offset] * len(d), 2 ** -np.array(d), ls='', marker=Line2D.filled_markers[1:][k],
                         mfc=cycle[j % 10], ms=8, mew=1, color='k')
    ax.set_xticks(np.arange(len(ddCt.keys())),
               ['%s (%s)' % (xtick_subs[key[0]], xtick_subs[key[1]]) for key in ddCt.keys()])
    ax.set_ylabel('relative mRNA')
    # add legend for targets
    handles, labels = ax.get_legend_handles_labels()
    sorted_handles_labels = sorted(zip(handles, labels), key=lambda x: x[1])
    handles, labels = zip(*sorted_handles_labels)
    leg1 = ax.legend(handles, labels, ncols=len(handles), bbox_to_anchor=(0.5, -0.12), loc='center', columnspacing=1,
                     frameon=False)
    # add legend for BioRep markers
    markers = [Line2D.filled_markers[1:][k] for k in range(len(np.unique(markers)))]
    handles = [
        Line2D([0], [0], marker=m, color='black', ls='None', ms=8, mfc='None') for m in markers
    ]
    labels = ['TechRep %d' % (k+1) for k in range(len(markers))]
    ax.legend(handles, labels, ncols=len(handles), bbox_to_anchor=(0.5, -0.18), loc='center', columnspacing=1,
              frameon=False)
    # add the first legend back
    ax.add_artist(leg1)

    return ddCt, fig


if __name__ == '__main__':

    basedir = '/Users/leonardharris/Library/CloudStorage/Box-Box/UArk Sys Bio Collab/Projects/TIBD/qPCR/OCT_NOV_2024/'
    dirnames =  ['BioRep1', 'BioRep2', 'BioRep3']  # 'TEST'
    figname = 'qPCR_relative_mRNA_ALL.pdf'
    fig = plt.figure(figsize=(6.4 * 2, 4.8 * len(dirnames)), constrained_layout=True)

    for i, dirname in enumerate([os.path.join(basedir, d) for d in dirnames]):

        platemap = read_platemap_data(dirname=dirname)
        for pmap_name in platemap.keys():
            print('Platemap:', pmap_name)
            display_platemap(platemap[pmap_name])
            print()

        Cq_data = extract_cell_substrate_data(platemap, control='18S')

        # print Cq values to the screen
        for key in Cq_data.keys():
            print()
            print(key)
            for key2 in Cq_data[key].keys():
                print('   %d:' % key2)
                for key3 in Cq_data[key][key2].keys():
                    print('      %s:' % key3, Cq_data[key][key2][key3])

        cell_line_groups = [['Bone Clone'], ['Parental']]
        ctrl_sample = [('Bone Clone', 'Plastic'), ('Parental', 'Plastic')]
        substr_groups = [['LC Random', 'LC Aligned', 'Plastic']]
        subs = {'Bone Clone': 'Bone', 'Parental': 'Par', 'Plastic': 'TC', 'LC Aligned': 'AC', 'LC Random': 'RC'}
        for j, group in enumerate(itertools.product(cell_line_groups, substr_groups)):
            Cq_data_subset = dict((key, Cq_data[key]) for key in Cq_data if key[0] in group[0] and key[1] in group[1])
            ddCt, fig = plot_relative_mRNA(Cq_data_subset, ref_gene='18S', ctrl_sample=ctrl_sample[j],
                                           add_subplot=(fig, int('%d2%d' % (len(dirnames), 2*i+j+1)), False))
                                         # add_subplot(figure, which figure=(row,col,idx), sharey=True|False)
    fig.savefig(os.path.join(basedir, figname), format='pdf')

    plt.show()
