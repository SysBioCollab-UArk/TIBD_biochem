import pandas as pd
import numpy as np
import os
import csv
from io import StringIO
import re
import itertools
import sys
import matplotlib.pyplot as plt
from matplotlib.lines import Line2D
from scipy.stats import linregress
import warnings
warnings.simplefilter("always")  # Ensures warning is shown


def warn_no_traceback(msg):
    """
    Emits a warning message without printing traceback or metadata.
    The message is printed to stderr.

    Parameters:
        msg (str): The warning message to display
    """
    with warnings.catch_warnings(record=True) as w:
        # Trigger the warning (caught internally)
        warnings.warn(msg)
        # Print just the message if warning was captured
        if w:
            print(f"{w[0].message}", file=sys.stderr)


def dixon_q_test(data, alpha=0.05):
    """
    Perform Dixon's Q Test for outliers in a small dataset.

    Parameters:
        data (list): List of numerical data points (3 <= len(data) <= 10).
        alpha (float): Significance level (default: 0.05).

    Returns:
        tuple: (is_outlier, outlier_idx, outlier_value, q_stat, q_critical)
            is_outlier (bool): Whether an outlier was detected.
            outlier_idx (int): Index of the suspected outlier.
            outlier_value (float): The suspected outlier value.
            q_stat (float): The calculated Q statistic.
            q_critical (float): The critical Q value for the test.
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
    if len(data) not in critical_values.keys():
        raise ValueError("Dixon's Q Test is only valid for n = %s. Your n = %d." %
                         (str(list(critical_values.keys())), len(data)))

    if alpha not in critical_values[len(data)].keys():
        raise ValueError("Supported alpha levels are %s. Your alpha = %g." %
                         (str(list(critical_values[len(data)].keys())), alpha))

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
    q_critical = critical_values[len(data)][alpha]

    # Determine if outlier
    is_outlier = q_stat > q_critical

    return is_outlier, outlier_idx, outlier_value, q_stat, q_critical


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
            '''print(line)'''
            while len(line[0]) == 1:  # row label is a letter, e.g., 'A', 'B', etc.
                row = line[0]
                for i in range(1, len(line)):
                    well_id = '%s%s' % (row, cols[i] if len(cols[i]) == 2 else '0%s' % cols[i])
                    sample = line[i].strip()  # if not line[i].isdigit() else int(line[i])
                    platemap[pmap_name][well_id] = [sample, '', '', '', '']  # sample, target, Cq, cell_line, substrate
                try:
                    line = next(reader)
                except StopIteration:  # in case we reach the end of the file
                    break
                '''print(line)'''
                for i in range(1, len(line)):
                    well_id = '%s%s' % (row, cols[i] if len(cols[i]) == 2 else '0%s' % cols[i])
                    target = line[i]
                    platemap[pmap_name][well_id][1] = target
                try:
                    line = next(reader)
                except StopIteration:  # in case we reach the end of the file
                    break
                '''print(line)'''

    # STEP 2: Read in the Cq values
    for pmap_name in platemap.keys():

        file_path = os.path.join(dirname, 'RawData', '%s - Quantification Summary.xlsx' % pmap_name)
        print("Reading '%s'" % os.path.basename(file_path))
        qpcr_data_xlsx = pd.read_excel(file_path, sheet_name=0, header=None)
        qpcr_data_csv = qpcr_data_xlsx.to_csv(header=False, index=False)

        reader = csv.reader(StringIO(qpcr_data_csv), delimiter=',')
        next(reader)  # read the first line
        for line in reader:
            well_id = line[1]
            Cq = line[6] if line[6] == '' else float(line[6])
            platemap[pmap_name][well_id][2] = Cq

    # STEP 3: Read in additional information (cell line/substrate or dilution) for each sample
    # '''
    file_path = os.path.join(dirname, 'cell_culture_samples.xlsx')
    samples_xlsx = pd.read_excel(file_path).astype(str)  #, header=None)
    samples_data = samples_xlsx.to_records(index=False)
    colnames = samples_data.dtype.names

    for pmap_name in platemap.keys():
        # make a dict mapping sample # to (cell line, substrate) tuple
        samples_dict = dict(
            zip(
                samples_data['Sample'],
                # [(cell, substr) for cell, substr in zip(samples_data['CellLine'], samples_data['Substrate'])]
                [(col1, col2) for col1, col2 in zip(samples_data[colnames[1]], samples_data[colnames[2]])]
            )
        )
        for well in platemap[pmap_name].keys():
            sample = platemap[pmap_name][well][0]
            if sample != '':  # make sure this is not an empty well
                platemap[pmap_name][well][3:4] = [samples_dict[sample][0], samples_dict[sample][1]]
    # '''
    return platemap


def extract_cell_substrate_data(platemap, control=None):

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
    '''print(targets)
    print(cell_lines)
    print(substrates)'''

    Cq_data = {}  # store qPCR measurements here

    # loop over (target, cell line, substrate) groups
    for group in itertools.product(targets, cell_lines, substrates):
        '''print(group)'''
        target = group[0]
        cell_line = group[1]
        substrate = group[2]
        qpcr_data_temp = {target: []}
        if control is not None:
            ctrl_key = '%s_%s' % (control, target)
            qpcr_data_temp[ctrl_key] = []
        # loop over platemaps
        for pmap_name in platemap.keys():
            pmap = platemap[pmap_name]
            # save values of (sample #, Cq) for the target gene
            samples_Cq = [(pmap[well][0], pmap[well][2]) for well in pmap.keys() if
                          pmap[well][2] != '' and  # missing data
                          pmap[well][1] == target and pmap[well][3] == cell_line and
                          pmap[well][4] == substrate]
            # if found values, get the Cq values for the control gene
            if len(samples_Cq) > 0:
                qpcr_data_temp[target].append(samples_Cq)
                if control is not None:
                    samples = np.unique([well[0] for well in samples_Cq]).tolist()
                    qpcr_data_temp[ctrl_key].append([(pmap[well][0], pmap[well][2]) for well in pmap.keys() if
                                                     pmap[well][2] != '' and  # missing data
                                                     pmap[well][0] in samples and
                                                     pmap[well][1] == control and
                                                     pmap[well][3] == cell_line and pmap[well][4] == substrate])

        # now consolidate Cq values for the target and control genes based on the sample number
        samples = np.unique([well[0] for q in qpcr_data_temp[target] for well in q])
        '''print('samples:', samples)'''
        Cq_target = [[] for _ in samples]
        if control is not None:
            Cq_ctrl = [[] for _ in samples]
        for i in range(len(qpcr_data_temp[target])):
            for j, sample in enumerate(samples):
                Cq_target[j].append([well[1] for well in qpcr_data_temp[target][i] if well[0] == sample])
                if control is not None:
                    Cq_ctrl[j].append([well[1] for well in qpcr_data_temp[ctrl_key][i] if well[0] == sample])
        #####
        # if control is not None:
        #     for i in range(len(qpcr_data_temp[ctrl_key])):
        #         for j, sample in enumerate(samples):
        #             Cq_ctrl[j].append([well[1] for well in qpcr_data_temp[ctrl_key][i] if well[0] == sample])
        #####
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
            if control is not None:
                Cq_data[(cell_line, substrate)][samples[i]][ctrl_key] = [cq for cq in Cq_ctrl[i] if len(cq) > 0]

        '''print(' ', (cell_line, substrate))
        for key in Cq_data[(cell_line, substrate)].keys():
            print('    %s:' % key, Cq_data[(cell_line, substrate)][key])'''

    return Cq_data


def outlier_warning_msg(sample, gene, idx, value, q_stat, q_critical):
    msg = f"Outlier detected:\n  Sample: {sample}\n  Gene: {gene}\n  Index: {idx}\n  Value: {value}\n  " + \
          f"Q Statistic: {q_stat:.3f}\n  Q Critical: {q_critical:.3f}"
    return msg


def calc_relative_mRNA(Cq_data, ref_gene, ctrl_sample, alpha=0.05, add_subplot=None, **kwargs):

    # process kwargs
    fig_title = kwargs.get('fig_title', None)
    fontsizes = kwargs.get('fontsizes', {})
    fs_title = fontsizes.get('title', None)
    fs_axis_labels = fontsizes.get('axis_labels', None)
    fs_axis_ticks = fontsizes.get('axis_ticks', None)
    fs_legend = fontsizes.get('legend', None)

    dCt = {}  # Delta Ct values
    # loop over (cell line, substrate) pairs
    for key in Cq_data.keys():
        '''print('%s:' % str(key), list(Cq_data[key].keys()))'''
        dCt[key] = {}
        # loop over sample numbers
        for key2 in Cq_data[key].keys():
            '''print('  %s' % key2)'''
            avg = {}
            targets = []
            # loop over targets + controls
            for key3 in Cq_data[key][key2].keys():
                '''print('    %s:' % key3, Cq_data[key][key2][key3], end=', ')'''
                avg[key3] = np.array([np.mean(cq) for cq in Cq_data[key][key2][key3]])
                '''print('avg:', avg[key3])'''
                # save targets
                if ref_gene not in key3:
                    targets.append(key3)
            # loop over targets
            for target in targets:
                # if this is the first time we've seen this target, add it to the dct[key] dict
                if target not in dCt[key].keys():
                    dCt[key][target] = []
                # calculate Delta_Ct
                dCt[key][target].append(list(avg[target] - avg['%s_%s' % (ref_gene, target)]))

    # get the names of the genes we have data for
    genes = list(np.unique([gene for cell_substr in dCt.keys() for gene in dCt[cell_substr].keys()]))

    ddCt = {}  # Delta Delta Ct values
    # First, calculate ddCt values for the control (cell line, substrate) pair and look for outliers
    DONE = False
    while not DONE:
        DONE = True
        # Calculate ddCt values
        ddCt[ctrl_sample] = dict(
            zip(
                list(dCt[ctrl_sample].keys()),
                [[np.array(d) - np.mean([xx for x in dCt[ctrl_sample][gene] for xx in x])
                  for d in dCt[ctrl_sample][gene]] for gene in dCt[ctrl_sample].keys()]
            )
        )
        '''print('== ddCt ==')
        for key in ddCt.keys():
            print(key)
            for key2 in ddCt[key].keys():
                print('\t', key2, ddCt[key][key2])
        print()'''

        # Identify and remove outliers using Dixon's Q Test
        for key2 in ddCt[ctrl_sample].keys():
            data = [dd for d in ddCt[ctrl_sample][key2] for dd in d]
            # check for outliers if there are at least 3 data points
            if len(data) >= 3:
                is_outlier, outlier_idx, outlier_value, q_stat, q_critical = \
                    dixon_q_test(2 ** -np.array(data), alpha=alpha)
                if is_outlier:
                    DONE = False
                    warn_no_traceback(
                        outlier_warning_msg(ctrl_sample, key2, outlier_idx, outlier_value, q_stat, q_critical)
                    )
                    # remove outlier from dCt array (not ddCt)
                    for d in dCt[ctrl_sample][key2]:
                        if outlier_idx < len(d):
                            d.pop(outlier_idx)
                            break
                        else:
                            outlier_idx -= len(d)
        '''print('== dCt ==')
        print(ctrl_sample)
        for key2 in dCt[ctrl_sample].keys():
            print('\t', key2, dCt[ctrl_sample][key2])
        print()
    print('NO MORE OUTLIERS!\n')'''

    # Now, calculate ddCt values for all other (cell line, substrate) pairs ...
    for cell_substr in [key for key in dCt.keys() if key != ctrl_sample]:
        ddCt[cell_substr] = dict(
            zip(
                list(dCt[cell_substr].keys()),
                [[np.array(d) - np.mean([xx for x in dCt[ctrl_sample][gene] for xx in x])
                  for d in dCt[cell_substr][gene]] for gene in dCt[cell_substr].keys()
                 if gene in dCt[ctrl_sample].keys()]
            )
        )
    '''print('== ddCt ==')
    for cell_substr in [key for key in dCt.keys() if key != ctrl_sample]:
        print(cell_substr)
        for key2 in ddCt[cell_substr].keys():
            print('\t', key2, ddCt[cell_substr][key2])
    print()'''

    # ...and identify and remove outliers using Dixon's Q Test
    DONE = False
    while not DONE:
        DONE = True
        for key in ddCt.keys():
            for key2 in ddCt[key].keys():
                data = [dd for d in ddCt[key][key2] for dd in d if not np.isnan(dd)]
                # check for outliers if there are at least 3 data points
                if len(data) >= 3:
                    is_outlier, outlier_idx, outlier_value, q_stat, q_critical = \
                        dixon_q_test(2 ** -np.array(data), alpha=alpha)
                    if is_outlier:
                        DONE = False
                        warn_no_traceback(
                            outlier_warning_msg(key, key2, outlier_idx, outlier_value, q_stat, q_critical)
                        )
                        # remove outlier from ddCt array
                        for d in ddCt[key][key2]:
                            if outlier_idx < len(d):
                                d[outlier_idx] = None
                                break
                            else:
                                outlier_idx -= len(d)
        '''print('== ddCt ==')
        for key in ddCt.keys():
            print(key)
            for key2 in ddCt[key].keys():
                print('\t', key2, ddCt[key][key2])
        print()
    print('NO MORE OUTLIERS!\n')'''

    # Plot relative mRNA levels
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
    '''print('Plotting relative mRNA expression')'''
    for i, key in enumerate(ddCt.keys()):
        '''print(key)'''
        for key2 in ddCt[key].keys():
            '''print('  ', key2)'''
            j = genes.index(key2)
            data = [dd for d in ddCt[key][key2] for dd in d if not np.isnan(dd)]
            '''print('    ', 2 ** -np.array(data))'''
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
    ax.set_title(label=fig_title, fontweight='bold', fontsize=fs_title)
    ax.set_xticks(np.arange(len(ddCt.keys())),
               ['%s (%s)' % (xtick_subs[key[0]], xtick_subs[key[1]]) for key in ddCt.keys()])
    ax.set_ylabel('relative mRNA', fontsize=fs_axis_labels)
    ax.tick_params(axis='both', which='major', labelsize=fs_axis_ticks)
    # add legend for targets
    handles, labels = ax.get_legend_handles_labels()
    sorted_handles_labels = sorted(zip(handles, labels), key=lambda x: x[1])
    handles, labels = zip(*sorted_handles_labels)
    leg1 = ax.legend(handles, labels, ncols=len(handles), bbox_to_anchor=(0.5, -0.12), loc='center', columnspacing=1,
                     frameon=False, fontsize=fs_legend)
    # add legend for BioRep markers
    markers = [Line2D.filled_markers[1:][k] for k in range(len(np.unique(markers)))]
    handles = [
        Line2D([0], [0], marker=m, color='black', ls='None', ms=8, mfc='None') for m in markers
    ]
    labels = ['TechRep %d' % (k+1) for k in range(len(markers))]
    ax.legend(handles, labels, ncols=min(5, len(handles)), bbox_to_anchor=(0.5, -0.18), loc='center', columnspacing=1,
              handletextpad=0.1, frameon=False, fontsize=fs_legend)
    # add the first legend back
    ax.add_artist(leg1)

    return ddCt, fig


def calc_standard_curves(dirname, r2threshold=0.98):
    platemap_sc = read_platemap_data(dirname=os.path.join(dirname, 'RawData', 'StandardCurve'))
    Cq_data_sc = extract_cell_substrate_data(platemap_sc, control=None)
    # print platemap to the screen
    '''for pmap_name in platemap_sc.keys():
        print('Platemap:', pmap_name)
        display_platemap(platemap_sc[pmap_name])
        print()'''
    # print Cq values to the screen
    genes_dilution = {}
    for key in Cq_data_sc.keys():
        '''print(key)'''
        for key2 in Cq_data_sc[key].keys():
            '''print('   %s:' % key2)'''
            for key3 in Cq_data_sc[key][key2].keys():
                if key3 not in list(genes_dilution.keys()):
                    genes_dilution[key3] = key[1]
                elif genes_dilution[key3] != key[1]:  # make sure this is the same dilution factor
                    raise Exception('Dilution factor for gene %s already recorded as %d. This value is %d.' %
                                    (key3, genes_dilution[key3], key[1]))

                '''print('      %s:' % key3, Cq_data_sc[key][key2][key3], np.mean(Cq_data_sc[key][key2][key3]))
        print()'''

    # Get log10(conc) vs. Cq values for each gene and calculate standard curves
    log10conc_Cq = {}
    std_curve = {}
    for gene in genes_dilution.keys():
        log10conc_Cq[gene] = np.vstack([
            (
                np.log10(key[0]), np.mean(data),
                np.std(data, ddof=1) / np.sqrt(len(data)) if len(data) > 1 else np.std(data, ddof=0)
            )
            for key in Cq_data_sc.keys() for key2 in Cq_data_sc[key].keys()
            if gene in Cq_data_sc[key][key2].keys() and (data := Cq_data_sc[key][key2][gene])
        ])
        # sort the array based on the x-axis values from least to most negative
        log10conc_Cq[gene] = log10conc_Cq[gene][log10conc_Cq[gene][:, 0].argsort()[::-1]]
        # fit a line to the data
        slope = 1
        rsquared = 0
        idxs = [i for i in range(len(log10conc_Cq[gene]))]
        # iterate and remove data points from the end of the data set until an acceptable fit is achieved
        while slope > 0 or rsquared < r2threshold:
            fit_line = linregress(log10conc_Cq[gene][idxs, 0], log10conc_Cq[gene][idxs, 1])
            min_max = (min(log10conc_Cq[gene][idxs, 1]), max(log10conc_Cq[gene][idxs, 1]))
            slope = fit_line.slope  # make sure the slope is negative
            rsquared = fit_line.rvalue ** 2
            '''print(f'{gene}: y = {slope:.1f}*x + {fit_line.intercept:.1f}, R^2 = {rsquared:.3f}')
            print(log10conc_Cq[gene][idxs])'''
            idxs.pop(-1)  # remove the last point
        '''print()'''
        # std_curve[gene] = (fit_line, genes_dilution[gene])  # second value is the dilution factor
        std_curve[gene] = {'fit_line': fit_line, 'dilution': genes_dilution[gene], 'range': min_max}
        '''print(std_curve[gene])'''

    # plot standard curves for all genes
    fig_sc = plt.figure(constrained_layout=True)
    cycle = plt.rcParams['axes.prop_cycle'].by_key()['color']  # standard colors
    genes = genes_dilution.keys()
    ncols = int(np.ceil(np.sqrt(len(genes))))
    nrows = int(np.ceil(len(genes) / ncols))
    for i, gene in enumerate([ctrl_gene] + sorted([g for g in genes if g != ctrl_gene])):
        ax = fig_sc.add_subplot(nrows, ncols, i + 1)
        color = 'k' if i == 0 else cycle[(i - 1) % 10]
        ax.set_title(gene, fontweight='bold', color=color)
        ax.errorbar(log10conc_Cq[gene][:, 0], log10conc_Cq[gene][:, 1], yerr=log10conc_Cq[gene][:, 2],
                    fmt='o', ms=8, color=color, capsize=6)
        # add fit line to the plot
        fit_line = std_curve[gene]['fit_line']
        yvals = fit_line.slope * log10conc_Cq[gene][:, 0] + fit_line.intercept
        ax.plot(log10conc_Cq[gene][:, 0], yvals, color='0.5', ls='--', lw=3,
                label=f'y = {fit_line.slope:.1f}*x+{fit_line.intercept:.1f}\n$R^2$ = {fit_line.rvalue ** 2:.3g}')
        ax.set_xlabel(r'log$_{10}$[concentration ($\mu$g/$\mu$l)]')
        ax.set_ylabel(r'$C_t$')
        # Retrieve handles and labels, reverse their order, and then create the legend
        handles, labels = plt.gca().get_legend_handles_labels()
        handles, labels = handles[::-1], labels[::-1]
        ax.legend(handles, labels, loc='best', frameon=False)

    return std_curve, fig_sc


if __name__ == '__main__':

    kwargs = {'fontsizes': {'title': 18, 'axis_labels': 16, 'axis_ticks': 16}}  #, 'legend': 12}}

    # basedir = '.'
    # dirnames = ['TEST']
    # figname =  'qPCR_relative_mRNA_TEST.pdf'

    basedir = '/Users/leonardharris/Library/CloudStorage/Box-Box/UArk Sys Bio Collab/Projects/TIBD/qPCR/MAY_AUG_2025' # '.'
    dirnames = ['BioRep1']  #, 'BioRep1_OLD'] #, 'BioRep2', 'BioRep3']  # ['BioRep3_OLD']
    figname = 'qPCR_relative_mRNA_BioRep1.pdf'

    fig_rel_mRNA = plt.figure(figsize=(6.4 * 2, 4.8 * len(dirnames)), constrained_layout=True)

    for i, dirname in enumerate([os.path.join(basedir, d) for d in dirnames]):

        print('=== %s ===\n' % os.path.split(dirname)[1])

        # Read platemap info
        platemap = read_platemap_data(dirname=dirname)
        # print platemap to the screen
        '''for pmap_name in platemap.keys():
            print('Platemap:', pmap_name)
            display_platemap(platemap[pmap_name])
            print()
        quit()'''

        # Get Cq data from platemap
        ctrl_gene = '18S'
        Cq_data_rel = extract_cell_substrate_data(platemap, control=ctrl_gene)
        '''
        Cq_data_abs = extract_cell_substrate_data(platemap, control=None)
        '''
        # print Cq values to the screen
        '''Cq_data = Cq_data_rel  # Cq_data_abs
        for key in Cq_data.keys():
            print(key)
            for key2 in Cq_data[key].keys():
                print('   %s:' % key2)
                for key3 in Cq_data[key][key2].keys():
                    print('      %s:' % key3, Cq_data[key][key2][key3])
            print()
        quit()'''

        # Get standard curves
        '''std_curve, fig_sc = calc_standard_curves(dirname, r2threshold=0.88)
        # print standard curves to the screen and plot data points on top of the curves
        for n, gene in enumerate([ctrl_gene] + sorted([g for g in std_curve.keys() if g != ctrl_gene])):
            fit_line = std_curve[gene]['fit_line']
            print('%s:' % gene, f'Ct = {fit_line.slope:.2f}*log10(conc) + {fit_line.intercept:.2f}, '
                                f'R^2 = {fit_line.rvalue ** 2:.3g}')
            Cq_gene = []  # collect data for this gene to plot on top of the standard curve
            for key in Cq_data_abs.keys():
                # print(key)
                for key2 in Cq_data_abs[key].keys():
                    # print('   %s:' % key2)
                    for key3 in Cq_data_abs[key][key2].keys():
                        if re.search('^%s' % gene, key3):
                            Cq_gene += [item for sublist in Cq_data_abs[key][key2][key3] for item in sublist]
                            # print('      (%s, %s):' % (gene, key3), Cq_gene)
            Cq_gene = np.array(Cq_gene)
            # print(Cq_gene)
            log10conc = (Cq_gene - fit_line.intercept) / fit_line.slope
            fig_sc.axes[n].plot(log10conc, Cq_gene, 'r^', mfc='None')
        print()
        # save standard curve figure
        fig_sc.savefig(os.path.join(dirname, '%s_standard_curves.pdf' % os.path.split(dirname)[1]), format='pdf')'''

        # Calculate relative and absolute mRNA levels
        cell_line_groups = [['Bone Clone'], ['Parental']]
        ctrl_sample = [('Bone Clone', 'Plastic'), ('Parental', 'Plastic')]
        substr_groups = [['LC Random', 'LC Aligned', 'Plastic']]
        subs = {'Bone Clone': 'Bone', 'Parental': 'Par', 'Plastic': 'TC', 'LC Aligned': 'AC', 'LC Random': 'RC'}
        for j, group in enumerate(itertools.product(cell_line_groups, substr_groups)):
            # Calculate relative mRNA levels
            fig_title = '%s: %s' % (dirnames[i], ' '.join(cell_line_groups[j]))
            Cq_data_subset = dict((key, Cq_data_rel[key]) for key in Cq_data_rel
                                  if key[0] in group[0] and key[1] in group[1])
            ddCt, fig_rel_mRNA = \
                calc_relative_mRNA(Cq_data_subset, ref_gene=ctrl_gene, ctrl_sample=ctrl_sample[j], alpha=0.1,
                                   add_subplot=(fig_rel_mRNA, int('%d2%d' % (len(dirnames), 2*i+j+1)), False),
                                   fig_title=fig_title, **kwargs)
                                 # add_subplot(figure, which figure=(row,col,idx), sharey=True|False)
            # print ddCt values to the screen
            # for key in ddCt.keys():
            #     print(key)
            #     for key2 in ddCt[key].keys():
            #         print('  %s:' % key2, ddCt[key][key2])
            #         data = [dd for d in ddCt[key][key2] for dd in d]
            #         print('    data:', data)
            #         print('    2^-data:', 2 ** -np.array(data))
            #         print('    2^-mean(data):', 2 ** -np.mean(data))

            # Calculate absolute mRNA levels
            '''plot_idx = {'18S': 0, 'Gli2': 1, 'ITGB3': 2, 'PTHLH': 3}
            Cq_data_subset = dict((key, Cq_data_abs[key]) for key in Cq_data_abs if key[0] in group[0] and
                                  key[1] in group[1])
            ##########
            print(Cq_data_abs)
            # calculate average 18S Ct value
            # check for systemic differences across samples # todo
            Cq_ctrl_gene = [cq[i] for key in Cq_data_abs.keys() for key2 in Cq_data_abs[key].keys()
                            for key3 in Cq_data_abs[key][key2].keys()
                            for cq in Cq_data_abs[key][key2][key3] for i in range(len(cq)) if key3 == ctrl_gene]
            print(Cq_ctrl_gene)
            print(len(Cq_ctrl_gene))
            plt.figure(constrained_layout=True)
            plt.plot([i for i in range(len(Cq_ctrl_gene))], Cq_ctrl_gene, 'o')
            plt.ylabel('Ct (%s)' % ctrl_gene)'''
            # plt.show()
            # quit()
            ##########
            '''abs_mRNA = {}
            for key in Cq_data_subset.keys():
                print(key)
                abs_mRNA[key] = {}
                for key2 in Cq_data_subset[key].keys():
                    print('  ', key2)
                    abs_mRNA[key][key2] = {}
                    for key3 in Cq_data_subset[key][key2].keys():
                        # gene = ctrl_gene if ctrl_gene in key3 else key3
                        print('    ', key3)
                        fit_line = std_curve[key3]['fit_line']
                        dilution_factor = std_curve[key3]['dilution']
                        min_max = std_curve[key3]['range']
                        print('         min_max:', min_max)
                        abs_mRNA[key][key2][key3] = []
                        ct_temp = []
                        for Ct_vals in Cq_data_subset[key][key2][key3]:
                            print('         rel:', Ct_vals)
            #                 print('       %s:' % key3,
            #                       f'Ct = {fit_line.slope:.2f}*log10(conc) + {fit_line.intercept:.2f}, '
            #                       f'R^2 = {fit_line.rvalue ** 2:.3g}')
                            abs_mRNA[key][key2][key3] += \
                                [10 ** ((ct - fit_line.intercept) / fit_line.slope) * dilution_factor
                                 for ct in Ct_vals if min_max[0] <= ct <= min_max[1]]
                            ct_temp += [ct for ct in Ct_vals if min_max[0] <= ct <= min_max[1]]
                        print('         abs:', abs_mRNA[key][key2][key3])
                        print('          ct:', ct_temp)
                        fig_sc.axes[plot_idx[key3]].plot(np.log10(abs_mRNA[key][key2][key3] / dilution_factor), ct_temp,
                                                         'bv', mfc='None')'''

            # print(Cq_data_subset)
            # print(abs_mRNA)
            # ctrl_gene_abs_mRNA = [m[i] for key in abs_mRNA.keys() for key2 in abs_mRNA[key].keys()
            #                       for key3 in abs_mRNA[key][key2].keys() if key3 == ctrl_gene
            #                       for m in abs_mRNA[key][key2][key3] for i in range(len(m))]
            # print(ctrl_gene_abs_mRNA)
            # print(len(ctrl_gene_abs_mRNA))
            # print(np.mean(ctrl_gene_abs_mRNA))
            # quit()
            #####

    # save relative mRNA figure
    fig_rel_mRNA.savefig(os.path.join(basedir, figname), format='pdf')

    # fig_rel_mRNA.axes[0].set_ylim(top=3)

    plt.show()
