How to run `read_qpcr_data.py` to get plots of relative mRNA expression 
for all cell lines/substrate pairs:

1. Get the platemap data by calling:


    platemap = read_platemap_data(dirname=dirname)
   
   * `dirname` is the name of the directory where all the required input files
   are located. For the example data included here, `dirname = 'TEST'`.

   * The platemap file must be an Excel spreadsheet with the name 
   `platemap.xlsx`. The spreadsheet should have all platemaps used in the 
   experiment. 

   * Each platemap in `platemap.xlsx` must be named based on the output 
   filenames in the `RawData` subfolder within the `dirname` target directory.

   * All data files in the `RawData` directory must be named 
   `[prefix] -  Quantification Summary.xlsx`, where `[prefix]` is the platemap 
   name in `platemap.xlsx`.

   * To connect the sample numbers in the platemap to experimental conditions 
   (cell line and substrate), a file called `cell_culture_samples.xlsx` must 
   also be included in the `dirname` target directory.


2. Extract the Cq (i.e., Ct) data by calling:


    Cq_data = extract_cell_substrate_data(platemap, control=ctrl_gene)

   * `platemap` is the platemap data extracted in Step 1, `ctrl_gene` is the
   control gene used for these experiments (e.g., 18S).


3. Calculate and plot relative mRNA expression using:


    ddCt, fig = calc_relative_mRNA(Cq_data, ref_gene=ctrl_gene, ctrl_sample=ctrl_sample, add_subplot=(fig, which_plot, sharey))

   * `Cq_data` could be all or a subset of the data returned in Step 2, 
   `ctrl_gene` is the same as in Step 2 (e.g., 18S), `ctrl_sample` is the 
   cell line/substrate pair used for comparison (e.g., bone clone cells on 
   tissue plastic).

   * `add_subplot` indicates that a new plot should be added to the current 
   figure. Useful when running over conditions in a loop.
     * `fig` is the current figure.
     * `which_plot` is a 3-digit integer indicating the number of rows, 
     columns, and the plot index to add. For example, `322` means the figure 
     has 3 rows, 2 columns, and we're adding a new plot at position `2`, 
     which corresponds to the 1st row, 2nd column in this case).
     * `sharey` is a Boolean variable (`True` or `False`) indicating whether 
     the y-axis values should be shared among all plots in the figure or not.

   * `ddCt` are the &Delta;&Delta;Ct values calculated from the data and `fig` 
   is figure the plots are being added to. Note that the `fig` returned is the 
   same as the `fig` passed to the `add_subplot` option, i.e., the figure to 
   be modified is passed in and then the modified figure is returned to the 
   user.
