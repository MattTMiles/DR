import numpy as np
import os

pulsar_list = "/fred/oz002/users/mmiles/MSP_DR/pulsar_list.txt"

for pulsar in [line.rstrip() for line in open(pulsar_list)]:
    print(pulsar)
    with open("/fred/oz002/users/mmiles/MSP_DR/paper_plots/latex_figure_appendix_file.txt","a") as f:
        f.write(r"\begin{figure*}"+" \n")
        f.write("\t "+r"\centerline{\includegraphics[width=0.9\paperwidth]{Appendix_figs/"+pulsar+"_grand.dly_4_notebook_comparison.pdf}} \n")
        f.write("\t "+r"\vspace{-0.5cm}"+" \n")
        #f.write("\t \centering \n")
        f.write("\t "+r"\centerline{\includegraphics[width=\paperwidth]{Appendix_figs/"+pulsar+"_noise_reduce_steps.pdf}} \n")
        f.write("\t "+r"\caption[]{Same as Figure D1, for PSR "+pulsar+"} \n")
        f.write("\t \label{fig: app_"+pulsar+"} \n")
        f.write("\t "+r"\vspace{-15pt}"+" \n")
        f.write(r"\end{figure*}"+" \n")
        f.write("\n")
    f.close()
    

