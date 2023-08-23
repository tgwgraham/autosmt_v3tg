# autosmt_v3tg
NIS Elements macro and associated Python scripts for running automated SMT experiments.

The core NIS Elements macro is v3tg.mac, which rasters over a grid on the surface of a coverslip.

The script masterscript.py runs continuously to locate suitable cells and write out mini-macros to re-center the stage on the cell of interest and resize the ROI.
Communication with v3tg.mac involves reading and writing temporary files, which is clunky, but it works.

You can optionally run realtime_analysis.py in the background to track movies and make plots of the analysis output. This requires the quot package
https://github.com/alecheckert/quot
