Welcome to the code for reproducing the analysis of origins in the context of 3D chromatin structure from the following paper:
<p><b>3D chromatin connectivity underlies replication origin efficiency in mouse embryonic stem cells</b>

by Karolina Jodkowska*, Vera Pancaldi*, Maria Rigau, Ricardo Almeida, José M. Fernández-Justel, Osvaldo Graña-Castro, Sara Rodríguez-Acebes, Miriam Rubio-Camarillo, Enrique Carrillo-de Santa Pau, David Pisano, Fátima Al-Shahrour, Alfonso Valencia, María Gómez, Juan Méndez
Currently available at https://www.biorxiv.org/content/10.1101/644971v2

Code for the integration of Replication Origin datasets and chromatin structure can be found in 3 main R noteboooks:
RepOri_revision_VOCHiC_resubAugust2022.Rmd
RepOri_revision_VPCHiC_resubAugust2022.Rmd
RepOri_revision_PCHiC_resubAugust2022.Rmd

The analysis of replication origin datasets projected onto chromatin interaction networks is performed using the ChAseR package available here
https://bitbucket.org/eraineri/chaser/src/master/

Another necessary package to process the chromatin network is igraph.

Data to run these scripts or results of some of the most intensive computations (randomized networks) can be found at the following Figshare repositories :
<p>General origin positions and efficiency with data for enrichment in specific regions or chromatin states:https://doi.org/10.6084/m9.figshare.20425917.v1


<p>Data to perform the VOCHiC analysis (Virtual Origin Capture HiC using HiC data from Bonev et al. 2017) https://doi.org/10.6084/m9.figshare.20442777.v1

  <p>Data to perform the origin analysis on Promoter Capture HiC networks from Schoenfelder et al. 2015 https://doi.org/10.6084/m9.figshare.20443077.v1
