import matplotlib.pyplot as plt
import numpy as np

# Example data
categories = ['Category 1', 'Category 2', 'Category 3', 'Category 4']
averages = [10, 15, 12, 18]
std_devs = [2, 3, 1.5, 2.5]
sample_counts = [30, 25, 40, 20]

# Plotting the bar chart
fig, ax = plt.subplots()
bars = ax.bar(categories, averages, yerr=std_devs, capsize=5, color='blue', alpha=0.7)

# Adding labels and title
ax.set_xlabel('Categories')
ax.set_ylabel('Average Value')
ax.set_title('Bar Chart with Average, Standard Deviation, and Sample Count')

# Adding sample count annotations on each bar
for bar, count in zip(bars, sample_counts):
    height = bar.get_height()
    ax.text(bar.get_x() + bar.get_width() / 2, height + 0.1, f'n={count}', ha='center', va='bottom')

# Adding legend
ax.legend(['Average with Std Dev'])

# Display the plot
plt.show()
