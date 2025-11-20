import math
import pandas as pd
import numpy as np
import os
import csv
from io import StringIO
import re
import itertools
import sys
import matplotlib.pyplot as plt
import matplotlib.lines as mlines
from matplotlib.lines import Line2D
import matplotlib.legend
from scipy.stats import linregress
import seaborn as sns
from itertools import combinations
from statannotations.Annotator import Annotator
import warnings, traceback
warnings.simplefilter("always")  # Ensures warning is shown
set_warnings_traceback = False  # set to True if want full traceback for warnings


# --- Monkey patch to enable per-annotation line colors ---
# Save original functions
_orig_annotate_pair = Annotator._annotate_pair
_orig_plot_line = Annotator._plot_line

def _annotate_pair_patched(self, annotation, ax_to_data, ann_list, orig_value_lim):
    # Store reference to current annotation
    self._current_annotation = annotation
    # Call original method
    return _orig_annotate_pair(self, annotation, ax_to_data, ann_list, orig_value_lim)

def _plot_line_patched(self, line_x, line_y):
    ann = getattr(self, "_current_annotation", None)
    color = getattr(ann, "line_color", self.color)
    lw = getattr(ann, "line_width", self.line_width)
    self.ax.plot(line_x, line_y, lw=lw, c=color, clip_on=False)

# Apply patches
Annotator._annotate_pair = _annotate_pair_patched
Annotator._plot_line = _plot_line_patched


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


def warn_with_traceback(message, category, filename, lineno, file=None, line=None):
    log = file if hasattr(file, "write") else sys.stderr
    traceback.print_stack(file=log)
    log.write(warnings.formatwarning(message, category, filename, lineno, line))


# display full traceback for warnings, if requested
if set_warnings_traceback:
    warnings.showwarning = warn_with_traceback


def outlier_warning_msg(sample, gene, idx, value, q_stat, q_critical):
    msg = f"Outlier detected:\n  Sample: {sample}\n  Gene: {gene}\n  Index: {idx}\n  Value: {value}\n  " + \
          f"Q Statistic: {q_stat:.3f}\n  Q Critical: {q_critical:.3f}"
    return msg


def strip_whitespace_from_df(df):
    # Trim column names
    df.columns = df.columns.str.strip()
    # Trim leading/trailing whitespace in all string columns
    obj = df.select_dtypes(include="object").columns
    df[obj] = df[obj].apply(lambda s: s.str.strip())
    # return modified df
    return df


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


def get_filepath_ignore_spaces(dirpath, prefix, suffix):
    files = os.listdir(dirpath)
    pattern = re.compile(rf"^{re.escape(prefix)}\s*-\s*{re.escape(suffix)}$")
    matching_files = [f for f in files if pattern.match(f)]
    if len(matching_files) == 1:
        file_path = os.path.join(dirpath, matching_files[0])
        return file_path
    else:
        raise FileNotFoundError(f"Expected one match, found {len(matching_files)}: {matching_files}")


def read_platemap_info(dirname):

    # STEP 1: Read in the platemap info
    file_path = os.path.join(dirname, 'platemap.xlsx')
    platemap_xlsx = pd.read_excel(file_path, header=None)
    platemap_csv = platemap_xlsx.to_csv(header=False, index=False)

    # make the platemap a dictionary with keys corresponding to the different platemap files
    platemap = {}

    reader = csv.reader(StringIO(platemap_csv), delimiter=',')
    reader = ([cell.strip() for cell in line] for line in reader)  # strip all whitespace
    for line in reader:
        '''print(line)'''
        if re.search(r'^Platemap:', line[0]):
            pmap_name = re.search(r'^Platemap:\s*(.+)$', line[0]).group(1).strip()
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

        file_path = get_filepath_ignore_spaces(
            os.path.join(dirname, 'RawData'), pmap_name, "Quantification Summary.xlsx")
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
    file_path = os.path.join(dirname, 'cell_culture_samples.xlsx')
    samples_xlsx = pd.read_excel(file_path).astype(str).map(str.strip)  #, header=None)
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

    # check for concentration units (standard curve)
    mass_units = colnames[1].split('_')
    vol_units = colnames[2].split('_')
    if len(mass_units) == len(vol_units) == 2:
        conc_units = '%s/%s' % (mass_units[1].replace('u', r'$\mu$'), vol_units[1].replace('u', r'$\mu$'))
        '''print('conc_units:', conc_units)'''
        mass_amount = samples_data[colnames[1]]
        if len(np.unique(mass_amount)) > 1:
            raise Exception("cDNA amounts in 2nd column of 'cell_culture_samples.xlsx' must all be equal: %s" % dirname)
        mass_amount = int(mass_amount[0])
        '''print('mass_amount:', mass_amount)'''
        return platemap, conc_units, mass_amount
    else:
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


def plot_mRNA_expression(ax, mRNA_expr_dict, genes, ref_sample=None, **kwargs):

    # process kwargs
    plot_title = kwargs.get('plot_title', None)
    fontsizes = kwargs.get('fontsizes', {})
    fs_title = fontsizes.get('title', None)
    fs_axis_labels = fontsizes.get('axis_labels', None)
    fs_axis_ticks = fontsizes.get('axis_ticks', None)
    fs_legend = fontsizes.get('legend', None)
    ylabel = kwargs.get('ylabel', 'mRNA expression')

    xtick_subs = {'Bone Clone': 'Bone', 'Parental': 'Par', 'Plastic': 'TC', 'LC Aligned': 'AC', 'LC Random': 'RC'}

    # Create mRNA data dictionary with 'condition', 'expression', and 'techRep' keys
    mRNA_data = dict([(key2, {'condition': [], 'expression': [], 'techRep': []}) for key in mRNA_expr_dict.keys()
                      for key2 in mRNA_expr_dict[key].keys()])
    for key in mRNA_expr_dict.keys():
        '''print('(%s, %s)' % (key[0], key[1]))'''
        for key2 in mRNA_expr_dict[key].keys():
            '''print('  %s:' % key2, mRNA_expr_dict[key][key2])'''
            data = []
            techRep = []
            rep = 0
            for d in mRNA_expr_dict[key][key2]:
                if len(d) > 0:
                    techRep.append(rep)
                    data += list(d)
                    rep += 1
            mRNA_data[key2]['condition'] += \
                ['%s (%s)' % (xtick_subs.get(key[0], key[0]), xtick_subs.get(key[1], key[1]))] * len(data)
            mRNA_data[key2]['expression'] += data
            mRNA_data[key2]['techRep'] += techRep

    # Create mRNA dataframe
    rows = []
    for gene, inner in mRNA_data.items():
        for cond, expr, tech in zip(inner['condition'], inner['expression'], inner['techRep']):
            rows.append({
                'condition': cond,
                'gene': str(gene),  # or gene.item() if you want plain str
                'expression': float(expr),  # converts np.float64 → Python float
                'techRep': tech
            })
    mRNA_expr_df = pd.DataFrame(rows)

    # Plot barplot with seaborn
    sample_order = np.sort(mRNA_expr_df['condition'].unique())  # enforce order of samples
    if ref_sample is not None:
        ref_sample = ('%s (%s)' %
                      (xtick_subs.get(ref_sample[0], ref_sample[0]), xtick_subs.get(ref_sample[1], ref_sample[1])))
        sample_order = [ref_sample] + [sample for sample in sample_order if sample != ref_sample]
    sns.barplot(ax=ax, data=mRNA_expr_df, x='condition', y='expression', hue='gene', order=sample_order,
                hue_order=genes, estimator='mean', errorbar='se', capsize=0.1, palette='bright')
    ax.legend(title=None)
    # shrink bar widths to create real gaps between bars
    shrink_factor = 0.75  # 0.75 of original width; tweak as needed
    for bar in ax.patches:
        orig_width = bar.get_width()
        center = bar.get_x() + orig_width / 2.0
        new_width = orig_width * shrink_factor
        bar.set_width(new_width)
        bar.set_x(center - new_width / 2.0)

    # Add stripplot with data points to each bar
    markers = Line2D.filled_markers[1:mRNA_expr_df['techRep'].nunique() + 1]
    for rep, marker in enumerate(markers):
        subset = mRNA_expr_df[mRNA_expr_df['techRep'] == rep]
        sns.stripplot(data=subset, x='condition', y='expression', hue='gene', order=sample_order, hue_order=genes,
                      dodge=True, jitter=True, marker=marker, size=10, linewidth=1, edgecolor='k', alpha=1,
                      palette='muted', ax=ax)

    # Create legends
    # add a legend for gene targets
    handles, labels = ax.get_legend_handles_labels()
    gene_leg = ax.legend(handles[:3], labels[:3], ncols=3, loc='center', bbox_to_anchor=(0.5, -0.12),
                         columnspacing=1, handletextpad=0.5, frameon=False, fontsize=None)  # fs_legend)
    ax.add_artist(gene_leg)
    # add a second legend for tech rep symbols
    techRep_handles = [
        mlines.Line2D([], [], color='black', marker=marker, linestyle='None', markersize=8,
                      markeredgewidth=1.5, markerfacecolor='none', label='TechRep %d' % rep)
        for rep, marker in enumerate(markers)
    ]
    techRep_leg = ax.legend(handles=techRep_handles, ncols=min(5, len(techRep_handles)), loc='center',
                            bbox_to_anchor=(0.5, -0.18), columnspacing=1, handletextpad=0.1,
                            frameon=False, fontsize=None)  # fs_legend)

    # Add statannotations *after* tweaking bar widths
    pairs = [] # define comparison pairs (compare conditions within each gene)
    # (('Bone (TC)', 'Gli2'), ('Bone (AC)', 'Gli2'))
    for gene in mRNA_expr_df['gene'].unique():
        df_gene = mRNA_expr_df[mRNA_expr_df['gene'] == gene]
        pairs.extend(list(combinations(
            list(set(zip(df_gene['condition'].tolist(), df_gene['gene'].tolist()))), 2))
        )
    pairs = sorted([tuple(sorted(p)) for p in pairs])  # sort pairs in a deterministic order

    annotator = Annotator(ax, pairs, data=mRNA_expr_df, x='condition', y='expression', hue='gene',
                          order=sample_order, hue_order=genes)
    alpha = 0.05  # p-value significance level
    annotator.configure(test='t-test_ind', text_format='star', comparisons_correction='bonferroni', loc='inside',
                        text_offset=0, line_height=0.04, line_width=2, verbose=1, hide_non_significant=True,
                        alpha=alpha,  # Only show brackets for p ≤ alpha
                        pvalue_thresholds=[
                            (alpha / 1000, r"${\ast}{\ast}{\ast}{\ast}$"),
                            (alpha / 100, r"${\ast}{\ast}{\ast}$"),
                            (alpha / 10, r"${\ast}{\ast}$"),
                            (alpha, r"$\ast$"),
                            (1, "ns")])
    # color the brackets and stars the same colors as the bars
    annotator.apply_test()  # gives list of Annotation objects
    colors = sns.color_palette('bright')
    colors = [colors[genes.index(pair[0][1])] for pair in pairs]
    pair_color_dict = dict(zip(pairs, colors))
    for ann in annotator.annotations:
        pair = tuple(sorted([ann.structs[0]['group'], ann.structs[1]['group']]))
        ann.line_color = pair_color_dict[pair]  # New feature enabled with patch
    annotator.annotate()

    # increase fontsize for stars over brackets
    for text in ax.texts:
        if text.get_text() in \
                {r"$\ast$", r"${\ast}{\ast}$", r"${\ast}{\ast}{\ast}$", r"${\ast}{\ast}{\ast}{\ast}$", "ns"}:
            text.set_fontsize(14)

    # add a third legend for significance testing levels
    if any(r.text != "ns" for r in annotator._get_results(num_comparisons=len(annotator.pairs))):
        ax.add_artist(techRep_leg)  # add the tech rep legend back in
        sig_legend_handles = [
            mlines.Line2D([], [], color='black', marker='', linestyle='None',
                          label=rf'{label} p $\leq$ {thresh:.0e}')
            for thresh, label in annotator.pvalue_thresholds[:-1]  # skip 'ns' if desired
        ]
        sig_leg = ax.legend(handles=sig_legend_handles, ncols=len(annotator.pvalue_thresholds), loc='center',
                            bbox_to_anchor=(0.5, -0.24), columnspacing=0, handletextpad=0,
                            frameon=False, fontsize=None)  # fs_legend)

    # Finish up plot
    ax.set_xlabel(None)
    ax.set_ylabel(ylabel, fontsize=fs_axis_labels)
    ax.tick_params(axis='both', which='major', labelsize=fs_axis_ticks)
    # Switch to scientific notation if >= 1000 (use abs() in case values are negative)
    if max(abs(ax.get_yticks())) >= 1000:
        ax.ticklabel_format(axis='y', style='sci', scilimits=(0, 0))
        ax.yaxis.get_offset_text().set_fontsize(fs_axis_ticks)
    ax.set_title(plot_title, fontsize=fs_title, fontweight='bold')
    sns.despine()


    def add_legend_to_fig(legend):
        handles = legend.legend_handles
        labels = [text.get_text() for text in legend.get_texts()]
        ax.figure.legend(handles, labels, ncols=legend._ncols, loc='center', columnspacing=legend.columnspacing,
                         handletextpad=legend.handletextpad, frameon=False, fontsize=fs_legend,
                         bbox_to_anchor=(0.5, offset))


    # Put the legend at the bottom of the figure
    legends = [child for child in ax.get_children() if isinstance(child, matplotlib.legend.Legend)]
    fig_legends = [child for child in ax.figure.get_children() if isinstance(child, matplotlib.legend.Legend)]
    offset = 0.02 + 0.03 * (max(len(legends), len(fig_legends)) - 1)  # 0.08
    # adjust figure area to make space for the figure legend
    bottom = offset + 0.02  # 0.1
    ax.figure.get_layout_engine().set(
        rect=(0, bottom, 1, 1 - bottom))  # (left, bottom, width, height)
    # loop over local legends
    for i, legend in enumerate(legends):
        # legends that already exist in the figure
        if i < len(fig_legends):
            # if there are more elements in the local legend than the figure legend, add the local legend to the figure
            if len(legend.legend_handles) > len(fig_legends[i].legend_handles):
                add_legend_to_fig(legend)
        # new legends to add to the figure
        else:
            add_legend_to_fig(legend)
        offset -= 0.03
        # delete local legend
        legend.remove()


def find_and_remove_outliers(data_dict, key, key2, remove_dict=None, pop=True, alpha=0.05):
    if remove_dict is None:
        remove_dict = data_dict
    is_outlier = False
    data = [dd for d in data_dict[key][key2] for dd in d if not np.isnan(dd)]
    # make sure there are at least 3 data points
    if len(data) >= 3:
        is_outlier, outlier_idx, outlier_value, q_stat, q_critical = dixon_q_test(data, alpha=alpha)
        if is_outlier:
            warn_no_traceback(outlier_warning_msg(key, key2, outlier_idx, outlier_value, q_stat, q_critical))
            # remove outlier from 'remove_dict'
            for d in remove_dict[key][key2]:
                if outlier_idx < len(d):
                    if pop:
                        d.pop(outlier_idx)
                    else:
                        d[outlier_idx] = None
                    break
                else:
                    outlier_idx -= len(d)
    return is_outlier


def calc_relative_mRNA(Cq_data, ref_gene, ref_sample, alpha=0.05, add_subplot=None, **kwargs):

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
    rel_mRNA = {}  # 2^-ddCt values
    # First, calculate ddCt values for the control (cell line, substrate) pair and look for outliers
    # NOTE the control sample values must equal 1, so if an outlier is detected the ddCt values must be recalculated
    print('Checking for outliers in data for control sample:', ref_sample)
    DONE = False
    while not DONE:
        DONE = True
        # Calculate ddCt values
        '''print('Calculate ddCt values for control sample')'''
        ddCt[ref_sample] = dict(
            zip(
                list(dCt[ref_sample].keys()),
                [[np.array(d) - np.mean([xx for x in dCt[ref_sample][gene] for xx in x])
                  for d in dCt[ref_sample][gene]] for gene in dCt[ref_sample].keys()]
            )
        )
        '''
        print('== ddCt ==')
        for key in ddCt.keys():
            print(key)
            for key2 in ddCt[key].keys():
                print('\t', key2, ddCt[key][key2])
        print()
        '''

        # Identify and remove outliers using Dixon's Q Test
        rel_mRNA[ref_sample] = {}
        for key2 in ddCt[ref_sample].keys():
            '''print('  ', key2)'''
            rel_mRNA[ref_sample][key2] = [2 ** -np.array(d) for d in ddCt[ref_sample][key2]]
            # remove outliers from dCt dict (not ddCt)
            found_outlier = find_and_remove_outliers(rel_mRNA, ref_sample, key2, remove_dict=dCt, pop=True,
                                                     alpha=alpha)
            if found_outlier:
                DONE = False
            '''
            print('     Outlier?', found_outlier)
            print('     DONE?', DONE)
            '''
        '''
        print('== dCt ==')
        print(ctrl_sample)
        for key2 in dCt[ctrl_sample].keys():
            print('\t', key2, dCt[ctrl_sample][key2])
        print()

    print('NO MORE OUTLIERS!\n')
    '''

    # Now, calculate ddCt values for all other (cell line, substrate) pairs ...
    for cell_substr in [key for key in dCt.keys() if key != ref_sample]:
        ddCt[cell_substr] = dict(
            zip(
                list(dCt[cell_substr].keys()),
                [[np.array(d) - np.mean([xx for x in dCt[ref_sample][gene] for xx in x])
                  for d in dCt[cell_substr][gene]] for gene in dCt[cell_substr].keys()
                 if gene in dCt[ref_sample].keys()]
            )
        )
        # ...and identify and remove outliers using Dixon's Q Test
        print('Checking for outliers in data for sample:', str(cell_substr))
        for key in ddCt.keys():
            '''print(key)'''
            rel_mRNA[key] = {}
            for key2 in ddCt[key].keys():
                '''print('  ', key2)'''
                rel_mRNA[key][key2] = [2 ** -np.array(d) for d in ddCt[key][key2]]
                DONE = False
                while not DONE:
                    # remove outliers from rel_mRNA dict
                    DONE = not find_and_remove_outliers(rel_mRNA, key, key2, pop=False, alpha=alpha)
                    '''print('     DONE?', DONE)'''
        '''
        print('== rel_mRNA ==')
        for key in rel_mRNA.keys():
            print(key)
            for key2 in ddCt[key].keys():
                print('\t', key2, rel_mRNA[key][key2])
        print()
        
    print('NO MORE OUTLIERS!\n')
    '''

    # Plot relative mRNA levels
    print('Plotting relative mRNA expression levels')
    if add_subplot is None:
        fig, ax = plt.subplots(constrained_layout=True)
    else:
        fig = add_subplot[0]
        fig_axes = fig.axes
        sharey = None
        if add_subplot[2] is True and len(fig_axes) > 0:
            sharey = fig_axes[-1] # share y-axes between all plots
        ax = fig.add_subplot(add_subplot[1], sharey=sharey)

    plot_mRNA_expression(ax, rel_mRNA, list(genes), ref_sample, **kwargs)

    return rel_mRNA, fig


def plot_Cq_data(ax, Cq_data, **kwargs):

    # process kwargs
    xlabel = kwargs.get('xlabel', 'index')
    ylabel = kwargs.get('ylabel', 'Cq values')
    title = kwargs.get('title', None)
    ylim = kwargs.get('ylim', (None, None))

    # plot data
    avg = np.mean(Cq_data)
    sdev = np.std(Cq_data)
    ax.plot(np.arange(len(Cq_data)), Cq_data, 'o', label=r'$C_t = %g \pm %g$' % (avg, sdev))
    ax.set_xlabel(xlabel)
    ax.set_ylabel(ylabel)
    ax.set_title(title, fontweight='bold')
    handles, labels = ax.get_legend_handles_labels()
    empty_handles = [Line2D([], [], linestyle='None', marker=None) for _ in labels]
    ax.legend(empty_handles, labels, loc='best', handlelength=0, handletextpad=0, frameon=False)
    ax.set_ylim(ylim)


def calc_standard_curves(dirname, r2threshold=0.98):
    dirpath = os.path.join(dirname, 'RawData', 'StandardCurve')
    platemap_std_curve, conc_units, cDNA_amount = read_platemap_info(dirname=dirpath)
    Cq_data_std_curve = extract_cell_substrate_data(platemap_std_curve, control=None)
    dilution_xlsx = pd.read_excel(os.path.join(dirpath, 'dilution.xlsx'), header=0)
    dilution_dict = dict(zip(dilution_xlsx['Gene'], dilution_xlsx['Dilution']))
    '''
    # print platemap to the screen
    for pmap_name in platemap_std_curve.keys():
        print('Platemap:', pmap_name)
        display_platemap(platemap_std_curve[pmap_name])
        print()
    '''

    # Get list of genes (and print Cq values to the screen)
    genes = []
    for key in Cq_data_std_curve.keys():
        '''print(key)'''
        for key2 in Cq_data_std_curve[key].keys():
            '''print('   %s:' % key2)'''
            for key3 in Cq_data_std_curve[key][key2].keys():
                '''print('      %s:' % key3, Cq_data_std_curve[key][key2][key3],
                      '<%g>' % np.mean(Cq_data_std_curve[key][key2][key3]))'''
                if key3 not in genes:
                    genes.append(key3)

    # Get log10(conc) vs. Cq values for each gene and calculate standard curves
    log10conc_Cq = {}
    std_curve = {'conc_units': conc_units, 'cDNA_amount': cDNA_amount}
    for gene in genes:
        log10conc_Cq[gene] = np.vstack([
            (
                np.log10(float(key[0]) / float(key[1])), np.mean(data),
                np.std(data, ddof=1) / np.sqrt(len(data)) if len(data) > 1 else np.std(data, ddof=0)
            )
            for key in Cq_data_std_curve.keys() for key2 in Cq_data_std_curve[key].keys()
            if gene in Cq_data_std_curve[key][key2].keys() and (data := Cq_data_std_curve[key][key2][gene])
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
        std_curve[gene] = {'fit_line': fit_line, 'Cq_range': min_max, 'dilution': dilution_dict[gene]}
        '''print(std_curve[gene])'''

    # plot standard curves for all genes
    fig_sc = plt.figure(constrained_layout=True)
    cycle = plt.rcParams['axes.prop_cycle'].by_key()['color']  # standard colors
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
        ax.plot(log10conc_Cq[gene][:, 0], yvals, color='0.5', ls='--', lw=2,
                label=f' y = {fit_line.slope:.1f}*x+{fit_line.intercept:.1f} \n $R^2$ = {fit_line.rvalue ** 2:.3g} ')
        ax.set_xlabel(r'log$_{10}$ concentration (%s)' % conc_units)
        ax.set_ylabel(r'$C_t$')
        # plot horizontal lines indicating the range of Cq values
        ax.axhline(std_curve[gene]['Cq_range'][0], ls=':', color='r')
        ax.axhline(std_curve[gene]['Cq_range'][1], ls=':', color='r')
        # Hide the handle in the legend
        ax.legend(loc='best', handlelength=0, handletextpad=0, facecolor='white', framealpha=1, edgecolor='white',
                  borderpad=0)

    return std_curve, fig_sc


# TODO: Modify this function to fit data to total # of cells vs. RNA yield to estimate cell counts
def get_cell_counts_hemocytometer(dirname):

    counts_df = pd.read_excel(os.path.join(dirname, 'cell_counts_hemocytometer.xlsx'), header=0)

    # solution volume
    vol_column = [col for col in counts_df.columns if 'Soln_Vol' in col][0]
    vol_units = vol_column.replace('Soln_Vol', '').replace('(', '').replace(')', '').strip()
    if counts_df[vol_column].dropna().nunique() == 1:
        soln_vol = counts_df[vol_column].dropna().unique()[0]
    else:
        raise Exception("Only one value of the cells solution volume in 'cell_counts.xlsx' is allowed. The "
                        "following values were read:", counts_df[vol_column].dropna().unique())

    # dilution factor
    if counts_df['Dilution'].dropna().nunique() == 1:
        dilution = counts_df['Dilution'].dropna().unique()[0]
    else:
        raise Exception("Only one value of the dilution factor in 'cell_counts.xlsx' is allowed. The "
                        "following values were read:", counts_df['Dilution'].dropna().unique())

    # hemocytometer cell counts
    counts_cols = [col for col in counts_df.columns if col not in
                   ['Sample', 'Dilution', 'Soln_Vol (%s)' % vol_units]]
    num_squares = len(np.unique([re.sub(r'(\w)\d+', r'\1', col) for col in counts_cols]))
    tot_columns = []
    for i in range(0, len(counts_cols), num_squares):
        columns = counts_cols[i:i + num_squares]
        tot_columns.append('tot_' + ''.join(columns))
        counts_df[tot_columns[-1]] = counts_df[columns].sum(axis=1)

    # average cell counts over all hemocytometer squares
    counts_df['avg_count'] = counts_df[tot_columns].mean(axis=1)

    # cell concentrations
    counts_df['cell_conc (cells/mL)'] = counts_df['avg_count'] / num_squares * dilution * 1e4

    # total cell counts
    if vol_units.lower() != 'ml':
        raise Exception("Units for cell concentrations (cells/mL) and solution volume (%s) must be consistent." %
                        vol_units)
    counts_df['tot_cells'] = counts_df['cell_conc (cells/mL)'] * soln_vol

    return counts_df.rename(columns={'tot_cells': 'NumCells'})


def get_cell_counts(dirname):
    if os.path.exists(os.path.join(dirname, 'cell_counts.xlsx')):
        return pd.read_excel(os.path.join(dirname, 'cell_counts.xlsx'), header=0)
    elif os.path.exists(os.path.join(dirname, 'cell_counts_hemocytometer.xlsx')):
        return get_cell_counts_hemocytometer(dirname)
    else:
        raise Exception("Can't find 'cell_counts.xlsx' or 'cell_counts_hemocytometer.xlsx' in %s" % dirname)


def get_rna_soln_and_cell_volumes(dirname):

    # === Get mRNA yields ===
    rna_df = pd.read_excel(os.path.join(dirname, 'rna_extraction.xlsx'), header=0)
    # solution volume
    vol_column = [col for col in rna_df.columns if 'Soln_Vol' in col][0]
    vol_units = vol_column.replace('Soln_Vol', '').replace('(', '').replace(')', '').strip()
    if rna_df[vol_column].dropna().nunique() == 1:
        rna_vol = rna_df[vol_column].dropna().unique()[0]
    else:
        raise Exception("Only one value of the RNA solution volume in 'rna_extraction.xlsx' is allowed. The "
                        "following values were read:", rna_df[vol_column].dropna().unique())
    rna_df = pd.DataFrame({'Sample': rna_df['Sample'].tolist(),
                           'Soln_Vol': [rna_vol] * len(rna_df),
                           'Soln_Vol_Units:': [vol_units] * len(rna_df)})
    '''print(rna_df.columns)
    print('RNA extraction fluid:\n', rna_df)'''

    # === Get cell counts ===
    counts_df = get_cell_counts(dirname)

    # === Get cell volumes ===
    cell_vols_df = pd.read_excel(os.path.join(dirname, 'cell_vols.xlsx'))

    # === Read in sample info ===
    samples_df = pd.read_excel(os.path.join(dirname, 'cell_culture_samples.xlsx'), header=0)
    samples_df = strip_whitespace_from_df(samples_df)

    # === Merge dataframes ===
    merged_df = rna_df.merge(counts_df, on='Sample', how='outer').rename(columns={vol_column: 'RNA_' + vol_column})
    merged_df = merged_df.merge(samples_df[['Sample', 'CellLine']], on='Sample', how='outer')
    merged_df = merged_df.merge(cell_vols_df[['CellLine', 'Volume', 'Units']], on='CellLine', how='left').rename(
        columns={'Volume': 'CellVolume', 'Units': "CellVolume_Units"})
    '''print(merged_df.columns)
    print(merged_df)'''

    return merged_df


def calc_absolute_mRNA(dirname, Cq_data, ref_gene, Cq_ref_gene, std_curve, ref_sample=None, alpha=0.05,
                       add_subplot=None, **kwargs):

    rna_and_cell_vols_df = get_rna_soln_and_cell_volumes(dirname)
    rna_and_cell_vols_df = rna_and_cell_vols_df.set_index('Sample')  # set index to sample number to speed up lookups

    abs_mRNA_all = {}
    for key in Cq_data.keys():
        '''print('(%s, %s)' % (key[0], key[1]))'''
        abs_mRNA_all[key] = {}
        for key2 in Cq_data[key].keys():
            '''print('  ', key2)'''
            # get volume of RNA extraction fluid, number of cells in the sample, and estimated cell volume
            soln_vol, num_cells, cell_vol = rna_and_cell_vols_df.loc[int(key2), ['Soln_Vol', 'NumCells', 'CellVolume']]
            sv_DIV_nc_DIV_cv = np.inf if num_cells == 0 else soln_vol / num_cells / cell_vol
            abs_mRNA_all[key][key2] = {}
            # loop over genes, with the ctrl_gene first
            if ctrl_gene in list(Cq_data[key][key2].keys()):
                for key3 in list(dict.fromkeys([ctrl_gene] + list(Cq_data[key][key2].keys()))):
                    '''print('    ', key3)'''
                    fit_line = std_curve[key3]['fit_line']
                    min_max = std_curve[key3]['Cq_range']
                    dilution_factor = std_curve[key3]['dilution']
                    abs_mRNA_all[key][key2][key3] = []
                    for Ct_vals in Cq_data[key][key2][key3]:
                        '''print('         Ct:', Ct_vals)'''
                        # scale the absolute mRNA concentrations by the ratio of the local and global 18S Ct values
                        # if local mean(Ct_18S) < global mean(Ct_18S), increase Ct -> decrease absolute conc
                        scale_factor = np.mean([ct for rep in Cq_data[key][key2][ctrl_gene] for ct in rep]) / \
                                       np.mean(Cq_ref_gene)
                        '''print('            scale_factor:', scale_factor)'''
                        abs_mRNA_all[key][key2][key3].append(
                            [10 ** ((ct - fit_line.intercept) / fit_line.slope) * dilution_factor * scale_factor
                             * sv_DIV_nc_DIV_cv  # soln_vol / num_cells / cell_vol
                             for ct in Ct_vals if ct <= min_max[1]])  # min_max[0] <= ct <= min_max[1] # TODO
                    '''print('         abs_mRNA:', abs_mRNA_all[key][key2][key3])'''
            else:
                print('No data for the internal control gene %s: (%s, %s), sample %s...Skipping.' %
                      (ctrl_gene, key[0], key[1], key2))

    # Get average absolute mRNA levels for each technical replicate
    abs_mRNA_techRep = {}
    genes = []
    for cell_substr in abs_mRNA_all.keys():
        '''print(cell_substr)'''
        abs_mRNA_techRep[cell_substr] = {}
        genes = np.unique(list(genes) + [key3 for key2 in abs_mRNA_all[cell_substr].keys()
                                         for key3 in abs_mRNA_all[cell_substr][key2].keys() if key3 != ctrl_gene])
        for gene in genes:
            '''print('\t%s' % gene)'''
            abs_mRNA_techRep[cell_substr][gene] = []
            for key2 in abs_mRNA_all[cell_substr].keys():
                '''print('\t\t%s' % key2)'''
                # use np.ma.masked_invalid() to ignore all NaN and +/-inf
                abs_mRNA_techRep[cell_substr][gene].append(
                    [np.ma.masked_invalid(x).mean() if len(x) > 0 and np.any(np.isfinite(x)) else np.nan
                     for x in abs_mRNA_all[cell_substr][key2].get(gene, [[np.nan]])]
                )
                '''
                print('\t\t\t', abs_mRNA_all[cell_substr][key2].get(gene, [[np.nan]]))
            print('\tabs_mRNA_techRep[%s][%s]:' % (str(cell_substr), gene), abs_mRNA_techRep[cell_substr][gene])
            '''

    # Identify and remove outliers using Dixon's Q Test
    for cell_substr in abs_mRNA_techRep.keys():
        print('Checking for outliers in data for sample:', str(cell_substr))
        '''print(cell_substr)'''
        for key2 in abs_mRNA_techRep[cell_substr].keys():
            '''print('\t%s:' % key2, abs_mRNA_techRep[cell_substr][key2])'''
            DONE = False
            while not DONE:
                DONE = not find_and_remove_outliers(abs_mRNA_techRep, cell_substr, key2, pop=True, alpha=alpha)
                '''print('\t\tDONE?', DONE)'''
    '''print()'''

    # Plot absolute mRNA levels
    print('Plotting absolute mRNA expression levels')
    if add_subplot is None:
        fig, ax = plt.subplots(constrained_layout=True)
    else:
        fig = add_subplot[0]
        fig_axes = fig.axes
        sharey = None
        if add_subplot[2] is True and len(fig_axes) > 0:
            sharey = fig_axes[-1]  # share y-axes between all plots
        ax = fig.add_subplot(add_subplot[1], sharey=sharey)

    plot_mRNA_expression(ax, abs_mRNA_techRep, list(genes), ref_sample, **kwargs)

    return abs_mRNA_techRep, fig


if __name__ == '__main__':

    # ##################################
    # TODO: code snippet for reading data from hemocytometer data file
    # dirname = '/Users/leonardharris/Downloads'
    # df = get_cell_counts_hemocytometer(dirname)
    # print(df)
    # seeded = [None] * 18
    # seeded[0:3] = [50000] * 3
    # seeded[3:6] = [100000] * 3
    # seeded[6:9] = [175000] * 3
    # seeded[9:12] = [250000] * 3
    # seeded[12:15] = [350000] * 3
    # seeded[15:18] = [400000] * 3
    # seeded += seeded
    #
    # plt.figure(constrained_layout=True)
    # plt.plot(seeded, df['NumCells'], 'o')
    # plt.xlabel('Seeded number of cells')
    # plt.ylabel('Observed number of cells')
    # plt.ticklabel_format(style='sci', axis='both', scilimits=(0, 0))
    #
    # # plot diagonal dashed line
    # lims = np.array([plt.gca().get_xlim(), plt.gca().get_ylim()])
    # min_val = lims.min()
    # max_val = lims.max()
    # plt.plot([min_val, max_val], [min_val, max_val], 'k--')
    # plt.xlim(min_val, max_val)
    # plt.ylim(min_val, max_val)
    #
    # plt.show()
    # quit()
    # ##################################

    # basedir = '.'
    # dirnames = ['TEST']
    # fig_prefix =  'TEST_qPCR'

    basedir = '/Users/leonardharris/Library/CloudStorage/Box-Box/UArk Sys Bio Collab/Projects/TIBD/qPCR/MAY_AUG_2025'
    dirnames = ['BioRep2', 'BioRep3']  #, 'BioRep1', 'BioRep1_OLD']
    fig_prefix = 'TIBD_qPCR'

    mRNA_kwargs = {'fontsizes': {'title': 18, 'axis_labels': 16, 'axis_ticks': 16, 'legend': 14}, 'sharey': True}

    fig_rel_mRNA = plt.figure(figsize=(6.4 * 2, 4.8 * len(dirnames)), constrained_layout=True)
    fig_abs_mRNA = plt.figure(figsize=(6.4 * 2, 4.8 * len(dirnames)), constrained_layout=True)
    fig_ctrl_gene = plt.figure(figsize=(6.4, 4.8 * 0.5 * len(dirnames)), constrained_layout=True)

    # Internal control gene
    ctrl_gene = '18S'

    # Loop over directories
    for i, dirname in enumerate([os.path.join(basedir, d) for d in dirnames]):

        print('\n=== %s ===\n' % os.path.split(dirname)[1])

        # Read platemap info
        platemap = read_platemap_info(dirname=dirname)
        '''
        # print platemap to the screen
        for pmap_name in platemap.keys():
            print('Platemap:', pmap_name)
            display_platemap(platemap[pmap_name])
            print()
        '''

        # Get Ct data based on platemap info
        Cq_data_rel = extract_cell_substrate_data(platemap, control=ctrl_gene)
        Cq_data_abs = extract_cell_substrate_data(platemap, control=None)
        '''
        # print Cq values to the screen
        Cq_data = Cq_data_abs  # Cq_data_abs Cq_data_rel
        for key in Cq_data.keys():
            print(key)
            for key2 in Cq_data[key].keys():
                print('   %s:' % key2)
                for key3 in Cq_data[key][key2].keys():
                    print('      %s:' % key3, Cq_data[key][key2][key3])
            print()
        '''

        # Get Ct values for internal control gene
        ax_ctrl_gene_ref = None if i == 0 else ax_ctrl_gene
        ax_ctrl_gene = fig_ctrl_gene.add_subplot(len(dirnames), 1, i + 1, sharey=ax_ctrl_gene_ref)
        Cq_ctrl_gene = [cq[j] for key in Cq_data_abs.keys() for key2 in Cq_data_abs[key].keys()
                        for key3 in Cq_data_abs[key][key2].keys()
                        for cq in Cq_data_abs[key][key2][key3] for j in range(len(cq)) if key3 == ctrl_gene]
        bottom = min(math.floor(min(Cq_ctrl_gene) / 10) * 10, np.inf if i == 0 else bottom) - 1
        top = max(math.ceil(max(Cq_ctrl_gene) / 10) * 10, -np.inf if i == 0 else top) + 1
        plot_Cq_data(ax_ctrl_gene, Cq_ctrl_gene, xlabel='sample index', ylabel='Ct (%s)' % ctrl_gene, title=dirnames[i],
                     ylim=(bottom, top))

        # Get standard curves
        print('\nCalculating standard curves')
        std_curve, fig_std_curve = calc_standard_curves(dirname, r2threshold=0.95)
        fig_std_curve.suptitle('%s' % os.path.split(dirname)[1], fontweight='bold')

        # Print standard curves to the screen and plot data points on top of the curves
        for n, gene in enumerate([ctrl_gene] + sorted([g for g in std_curve.keys() if
                                                       g not in [ctrl_gene, 'conc_units', 'cDNA_amount']])):
            fit_line = std_curve[gene]['fit_line']
            print('%s:' % gene, f'Ct = {fit_line.slope:.2f} * log10(conc) + {fit_line.intercept:.2f}, '
                                f'R^2 = {fit_line.rvalue ** 2:.3g}')
            Cq_gene = []  # collect data for this gene to plot on top of the standard curve
            for key in Cq_data_abs.keys():
                '''print(key)'''
                for key2 in Cq_data_abs[key].keys():
                    '''print('   %s:' % key2)'''
                    for key3 in Cq_data_abs[key][key2].keys():
                        if re.search('^%s' % gene, key3):
                            Cq_gene += [item for sublist in Cq_data_abs[key][key2][key3] for item in sublist]
                            '''print('      (%s, %s):' % (gene, key3), Cq_gene)'''
            Cq_gene_range = std_curve[gene]['Cq_range']
            Cq_gene_within_range = [cq for cq in Cq_gene if Cq_gene_range[0] <= cq <= Cq_gene_range[1]]
            Cq_gene_out_of_range = [cq for cq in Cq_gene if cq < Cq_gene_range[0] or cq > Cq_gene_range[1]]
            for Cq_vals, marker, color in zip([Cq_gene_within_range, Cq_gene_out_of_range], ['^', 'x'], ['r', 'gold']):
                log10conc = (np.array(Cq_vals) - fit_line.intercept) / fit_line.slope
                fig_std_curve.axes[n].plot(log10conc, Cq_vals, ls='', marker=marker, mfc='None', color=color,
                                           zorder=0)
        # Save standard curve figure
        outfile = '%s_standard_curves.pdf' % os.path.split(dirname)[1]
        fig_std_curve.savefig(os.path.join(dirname, outfile), format='pdf')
        print('\nSaving %s to...' % outfile)
        print('...%s' % os.path.abspath(dirname))

        # Calculate relative and absolute mRNA levels
        cell_line_groups = [['Bone Clone'], ['Parental']]
        ctrl_sample = [('Bone Clone', 'Plastic'), ('Parental', 'Plastic')]
        substr_groups = [['LC Random', 'LC Aligned', 'Plastic']]
        subs = {'Bone Clone': 'Bone', 'Parental': 'Par', 'Plastic': 'TC', 'LC Aligned': 'AC', 'LC Random': 'RC'}
        for j, group in enumerate(itertools.product(cell_line_groups, substr_groups)):

            print('\n*** Sample group: %s ***' % str(group))

            plot_title = '%s: %s' % (dirnames[i], ' '.join(cell_line_groups[j]))

            # Calculate relative mRNA levels
            print('\nCalculating relative mRNA expression levels')
            Cq_data_subset = dict((key, Cq_data_rel[key]) for key in Cq_data_rel if key[0] in group[0] and
                                  key[1] in group[1])
            rel_mRNA, fig_rel_mRNA = \
                calc_relative_mRNA(
                    Cq_data_subset, ctrl_gene, ref_sample=ctrl_sample[j], alpha=0.1,
                    add_subplot=(fig_rel_mRNA, int('%d2%d' % (len(dirnames), 2*i+j+1)),
                                 mRNA_kwargs.get('sharey', False)),
                    # add_subplot(figure, which_figure = (row, col, idx), sharey = True | False)
                    plot_title=plot_title, ylabel='relative mRNA', **mRNA_kwargs
                )
            '''
            # print relative mRNA values to the screen
            for key in rel_mRNA.keys():
                print(key)
                for key2 in rel_mRNA[key].keys():
                    print('  %s:' % key2, rel_mRNA[key][key2])
                    data = [dd for d in rel_mRNA[key][key2] for dd in d]
                    print('    data:', data)
                    print('    mean(data):', np.mean(data))
            '''

            # Calculate absolute mRNA levels
            print('\nCalculating absolute mRNA expression levels')
            Cq_data_subset = dict((key, Cq_data_abs[key]) for key in Cq_data_abs if key[0] in group[0] and
                                  key[1] in group[1])
            abs_mRNA, fig_abs_mRNA = \
                calc_absolute_mRNA(
                    dirname, Cq_data_subset, ctrl_gene, Cq_ctrl_gene, std_curve, ref_sample=ctrl_sample[j], alpha=0.1,
                    add_subplot = (fig_abs_mRNA, int('%d2%d' % (len(dirnames), 2*i+j+1)),
                                   mRNA_kwargs.get('sharey', False)),
                    # add_subplot(figure, which_figure = (row, col, idx), sharey = True | False)
                    plot_title = plot_title, ylabel = 'absolute mRNA (%s)' % std_curve['conc_units'],
                    **mRNA_kwargs
                )
            '''
            # print absolute mRNA values to the screen
            for key in abs_mRNA.keys():
                print('(%s, %s)' % (key[0], key[1]))
                for key2 in abs_mRNA[key].keys():
                    print('  %s:' % key2, abs_mRNA[key][key2])
                    data = [dd for d in abs_mRNA[key][key2] for dd in d]
                    print('    data:', np.array(data))
                    print('    mean(data):', np.mean(data))
            '''

    # Save figures
    print('\nSaving...')
    outfile = fig_prefix + '_rel_mRNA.pdf'
    print('   %s' % outfile)
    fig_rel_mRNA.savefig(os.path.join(basedir, outfile), format='pdf')
    outfile = fig_prefix + '_abs_mRNA.pdf'
    print('   %s' % outfile)
    fig_abs_mRNA.savefig(os.path.join(basedir, outfile), format='pdf')
    outfile = fig_prefix + '_ctrl_%s.pdf' % ctrl_gene
    print('   %s' % outfile)
    fig_ctrl_gene.savefig(os.path.join(basedir, outfile), format='pdf')
    print('...to %s' % os.path.abspath(basedir))

    plt.show()
