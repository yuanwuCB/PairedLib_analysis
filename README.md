# PairedLib_analysis
Python scripts to analyze sgRNA-target paired library data.

Dependency:

(1). Pandas v2.0.3.
(2). Numpt v1.24.3.
(3). Seaborn v0.12.2.
(4). Matplotlib v3.7.3.

For every batch analysis of NGS data, run mapping.py first to map reads to the reference file (So-called whitelist). And then run Combine_Mapped_To_One_Table.py to make one excel table to identify sites that pass the coverage filter in all editors. Finally, use PairedLib_analyze.py to do analysis or plotting. 
