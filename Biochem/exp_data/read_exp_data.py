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

# for drug in drugs:
#     print(drugs_map.get(drug, drug))
# quit()

# Force
force_dict = {}
markers = ['o', '^', 's']
for target, marker in zip(np.unique([d['Target'] for d in raw_data if d['Drug'] == 'Force']), markers):
    # plt.figure()
    xdata = [d['Dose'] for d in raw_data if d['Drug'] == 'Force' and d['Target'] == target]
    ydata = np.array([d['Count'] for d in raw_data if d['Drug'] == 'Force' and d['Target'] == target])
    # print(xdata)
    # print(ydata)
    stiffness, counts = np.unique(xdata, return_counts=True)
    # print(stiffness, counts)
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
    # plt.plot(np.log10(xdata), ydata/ydata[0], 'o', lw=3, label=target)
    # plt.errorbar(np.log10(stiffness), mRNA_avg/mRNA_avg[0], mRNA_se, ls='-', marker=marker, ms=10, capsize=10, lw=3, label=target)
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
        # plt.title('drug: %s, target: %s' % (drug, target))
        doses = np.unique([d['Dose'] for d in raw_data if d['Drug'] == drug and d['Target'] == target])
        ncol = 2 if len(doses) > 3 else 1
        for dose in doses:
            # print('dose: %g' % dose)
            xy_data = [(d['Hours'], d['Count'], d['Units']) for d in raw_data if d['Drug'] == drug
                       and d['Target'] == target and d['Dose'] == dose]
            xdata = [xy_data[i][0] for i in range(len(xy_data))]
            ydata = [xy_data[i][1] for i in range(len(xy_data))]
            units = np.unique([xy_data[i][2] for i in range(len(xy_data))])
            if len(units) > 1:
                print('Error: Please ensure that units for all doses are the same (%s, %s)' % (drug, target))
            else:
                units = re.sub(r'^u', r'$\\mu$', units[0])
            xdata, counts = np.unique(xdata, return_counts=True)
            ydata_avg = []
            ydata_se = []
            start = 0
            for n in range(len(counts)):
                ydata_avg.append(np.mean(ydata[start:start+counts[n]]))
                ydata_se.append(np.std(ydata[start:start+counts[n]], ddof=1) / np.sqrt(counts[n]))
                start += counts[n]
            ydata_avg = np.array(ydata_avg)

            # plots
            # plt.plot(xdata, ydata_avg/ydata_avg[0], 'o-', lw=3, ms=10, label='%d %s' % (dose, units))
            # plt.plot(xdata, ydata_avg, 'o-', lw=3, ms=10, label='%d %s' % (dose, units))
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
                plt.xlim(left=0)
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
