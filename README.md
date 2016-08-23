# BETA
BETA: Bayesian Estimation of Transcription Factor Activity


Comparative analysis of genome-wide expression profiles can provide critical biological insights, and the ability to compare these profiles across studies is becoming increasingly important with the growing number of data-sets in public repositories. In order to gain mechanistic insights, gene expression profiles can be combined with information on DNA-binding sites of transcription factors (TFs) to detect transcription factor activity (by analysis of target gene sets). Bayesian Estimation of Transcription Factor Activity (BETA) is a method to quantify TF activity (rather than simply detecting it), and allows for direct comparison between multiple studies. Specifically, BETA estimates the probability density function for transcription factor activity, which is defined as the log-odds ratio between the observed frequency of TF targets among differentially-expressed genes compared with the expected TF target frequency among a set of background genes. Individual BETA activities can be compared to 0 in order to detect significant TF activity, and two BETA activities can be compared to each other in order to detect quantitative differences in TF activities.

![image](https://cloud.githubusercontent.com/assets/21067499/17899839/b8e75278-692a-11e6-801f-e3d6d316c7d0.png)


When including the results of the BETA package in a publication, please cite the following paper:


Thakar, Juilee, Boris M. Hartmann, Nada Marjanovic, Stuart C. Sealfon, and Steven H. Kleinstein. "Comparative analysis of anti-viral transcriptomics reveals novel effects of influenza immune antagonism." BMC immunology 16, no. 1 (2015): 46.



For questions, comments, or requests, please contact Juilee at Juilee_Thakar@URMC.Rochester.edu


