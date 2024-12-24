import numpy as np
import matplotlib.pyplot as plt
from matplotlib import container
import re

infile = 'PCRData_CSV_011422.csv'
raw_data = np.genfromtxt(infile, names=True, dtype=None, delimiter=',', encoding='UTF-8')

# print(raw_data)
print(raw_data.dtype.names)

drugs = np.unique(raw_data['Drug'])
print(drugs)

drugs_map = {'AVB3InhibitorCilengitide': 'Cilengitide',
             'RecombinantTGFb': r'TGF-$\beta$',
             'TGFbInhibitor1D11': '1D11',
             'Gli2InhibitorGANT58': 'GANT58',
             'p38MAPKInhibitorSB202190': 'SB202190',
             'SMAD3InhibitorSIS3': 'SIS3'}


def get_drug_target_dose_data(drug, target, dose):
    print('%s: %s (dose = %g)' % (drug, target, dose))
    xy_data = [(d['Hours'], d['Count'], d['Units']) for d in raw_data
               if d['Drug'] == drug and d['Target'] == target and d['Dose'] == dose]
    xdata = [xy_data[i][0] for i in range(len(xy_data))]
    ydata = [xy_data[i][1] for i in range(len(xy_data))]
    xdata, counts = np.unique(xdata, return_counts=True)  # time
    ydata_avg = []
    ydata_se = []
    start = 0
    for n in range(len(counts)):
        ydata_avg.append(np.mean(ydata[start:start + counts[n]]))
        ydata_se.append(np.std(ydata[start:start + counts[n]], ddof=1) / np.sqrt(counts[n]))
        start += counts[n]
    ydata_avg = np.array(ydata_avg)
    ydata_se = np.array(ydata_se)

    return xdata, ydata_avg, ydata_se

##########
'''
import csv

drugs = ['TGFbInhibitor1D11', 'TGFbInhibitor1D11']  # , 'AVB3InhibitorCilengitide', 'AVB3InhibitorCilengitide']
targets = [['PTHrP', 'Gli2'], ['PTHrP', 'Gli2']]  # , ['ITGB3', 'PTHrP', 'Gli2'], ['ITGB3', 'PTHrP', 'Gli2']]
doses = [0, 10]  # , 0, 20]

# drugs = ['TGFbInhibitor1D11', 'TGFbInhibitor1D11']
# targets = [['PTHrP'], ['PTHrP']]
# doses = [0, 10]

drug_dict = {'TGFbInhibitor1D11': 'x1D11',
             'AVB3InhibitorCilengitide': 'Cilengitide'}

target_dict = {'PTHrP': 'mPTHrP_Obs',
               'Gli2': 'mGli2_Obs',
               'ITGB3': 'mITGB3_Obs'}

# 'TGFbInhibitor1D11', ['PTHrP', 'Gli2'], 0
# 'TGFbInhibitor1D11', ['PTHrP', 'Gli2'], 10
# 'AVB3InhibitorCilengitide', ['ITGB3'], 0
# 'AVB3InhibitorCilengitide', ['ITGB3'], 20

for i, (drug, target, dose) in enumerate(zip(drugs, targets, doses)):
    print(i)
    # print('===', drug, target, dose, '===')
    headers = [target_dict[tgt] for tgt in target]
    exp_data = np.array([get_drug_target_dose_data(drug, tgt, dose) for tgt in target])
    exp_time = exp_data[:, 0, :].T
    exp_mean = exp_data[:, 1, :].T
    exp_serr = exp_data[:, 2, :].T
    # print(headers)
    # print(exp_time)
    # print(exp_mean)
    # print(exp_serr)
    for suffix, data in zip(['time', 'avg', 'se'], [exp_time, exp_mean, exp_serr]):
        file = open('../tumorCell_exp_data_%s_%d.csv' % (suffix, i), 'w')
        csvwriter = csv.writer(file)
        csvwriter.writerow(headers)
        csvwriter.writerows(data)
        file.close()

    # plot data for each experiment
    plt.figure()
    plt.title('%s: %d %s' % (drug_dict[drug], dose, r'$\mu$g/ml'), fontsize=20)
    for j, tgt in enumerate(target):
        plt.errorbar(exp_time[:, j], exp_mean[:, j], exp_serr[:, j], ls='-', marker='o', ms=10, capsize=10, lw=3,
                     label='%s' % tgt)
        plt.xlim(left=0)
        plt.xlabel('time (h)', fontsize=20)
        plt.ylabel('mRNA', fontsize=20)
        plt.xticks(fontsize=20)
        plt.yticks(fontsize=20)
        plt.legend(fontsize=16)
        plt.tight_layout()
plt.show()
quit()
'''
##########

# Force
force_dict = {}
markers = ['o', '^', 's']
for target, marker in zip(np.unique([d['Target'] for d in raw_data if d['Drug'] == 'Force']), markers):
    plt.figure('stiffness')
    xdata = [d['Dose'] for d in raw_data if d['Drug'] == 'Force' and d['Target'] == target]
    ydata = np.array([d['Count'] for d in raw_data if d['Drug'] == 'Force' and d['Target'] == target])
    stiffness, counts = np.unique(xdata, return_counts=True)
    start_idx = 0
    mRNA_avg = []
    mRNA_se = []
    for i in range(len(stiffness)):
        amount = []
        for j in range(len(ydata)):
            if xdata[j] == stiffness[i]:
                amount.append(ydata[j])
        mRNA_avg.append(np.mean(amount))
        mRNA_se.append(np.std(amount, ddof=1) / np.sqrt(counts[i]))

    # plots
    plt.errorbar(np.log10(stiffness), mRNA_avg, mRNA_se, ls='-', marker=marker, ms=10, capsize=10, lw=3, label=target)
    plt.xlim(left=0)
    plt.xlabel(r'log$_{10}$ stiffness (kPa)', fontsize=20)
    plt.ylabel('mRNA', fontsize=20)
    plt.xticks(fontsize=20)
    plt.yticks(fontsize=20)
    plt.legend()
    # Remove error bars from legend
    handles, labels = plt.gca().get_legend_handles_labels()
    handles = [h[0] if isinstance(h, container.ErrorbarContainer) else h for h in handles]
    plt.legend(handles, labels, loc=0, fontsize=16)
    plt.tight_layout()

# Other drugs
for drug in drugs[drugs != 'Force']:
    for target in np.unique([d['Target'] for d in raw_data if d['Drug'] == drug]):
        plt.figure('%s: %s' % (drug, target))
        print('%s: %s' % (drug, target))
        doses = np.unique([d['Dose'] for d in raw_data if d['Drug'] == drug and d['Target'] == target])
        # TODO: Get units here
        ncol = 2 if len(doses) > 3 else 1
        for dose in doses:
            xy_data = [(d['Hours'], d['Count'], d['Units']) for d in raw_data if d['Drug'] == drug
                       and d['Target'] == target and d['Dose'] == dose]
            xdata = [xy_data[i][0] for i in range(len(xy_data))]
            ydata = [xy_data[i][1] for i in range(len(xy_data))]
            units = np.unique([xy_data[i][2] for i in range(len(xy_data))])
            if len(units) > 1:
                print('Error: Please ensure that units for all doses are the same (%s, %s)' % (drug, target))
            else:
                units = re.sub(r'^u', r'$\\mu$', units[0])
            xdata, counts = np.unique(xdata, return_counts=True)  # time
            ydata_avg = []
            ydata_se = []
            start = 0
            for n in range(len(counts)):
                ydata_avg.append(np.mean(ydata[start:start+counts[n]]))
                ydata_se.append(np.std(ydata[start:start+counts[n]], ddof=1) / np.sqrt(counts[n]))
                start += counts[n]
            ydata_avg = np.array(ydata_avg)

            # plots
            plt.errorbar(xdata, ydata_avg, ydata_se, ls='-', marker='o', ms=10, capsize=10, lw=3, label='%d' % dose)
            plt.xlim(left=0)
            plt.xlabel('time (h)', fontsize=20)
            plt.ylabel('%s mRNA' % target, fontsize=20)
            plt.xticks(fontsize=20)
            plt.yticks(fontsize=20)
            plt.legend()
            # Remove error bars from legend
            handles, labels = plt.gca().get_legend_handles_labels()
            handles = [h[0] if isinstance(h, container.ErrorbarContainer) else h for h in handles]
            plt.legend(handles, labels, loc=0, ncol=ncol, columnspacing=1, fontsize=16,
                       title='%s (%s)' % (drugs_map.get(drug, drug), units), title_fontsize=16)
            plt.tight_layout()
            # plt.savefig()

            #####
            if target in ['Gli2', 'PTHrP'] and dose == 0:
                plt.figure(target)
                plt.errorbar(xdata, ydata_avg, ydata_se, ls='-', marker='o', ms=10, capsize=10, lw=3,
                             label='%s' % drugs_map.get(drug, drug))
                plt.xlim(left=0, right=80)  ###
                plt.xlabel('time (h)', fontsize=20)
                plt.ylabel('%s mRNA (control)' % target, fontsize=20)
                plt.xticks(fontsize=20)
                plt.yticks(fontsize=20)
                plt.legend()
                # Remove error bars from legend
                handles, labels = plt.gca().get_legend_handles_labels()
                handles = [h[0] if isinstance(h, container.ErrorbarContainer) else h for h in handles]
                plt.legend(handles, labels, loc=0, ncol=ncol, columnspacing=1, fontsize=16)
                plt.tight_layout()
                # switch back to original plot
                plt.figure('%s: %s' % (drug, target))
            #####

plt.show()
