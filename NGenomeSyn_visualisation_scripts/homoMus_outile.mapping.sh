/shared/ifbstor1/software/miniconda/envs/mummer4-4.0.0rc1/bin/nucmer  --mum     --mincluster 500  -t 6   homoMus_outile.A.fa  homoMus_outile.B.fa    -p homoMus_outile 
/shared/ifbstor1/software/miniconda/envs/mummer4-4.0.0rc1/bin/delta-filter   -1 -i 90 -l 2000   homoMus_outile.delta  >  homoMus_outile.filter.delta 
/shared/ifbstor1/software/miniconda/envs/mummer4-4.0.0rc1/bin/show-coords -c -r   homoMus_outile.filter.delta     >  homoMus_outile.filter.coords
perl  ../../NGenomeSyn/bin/GetTwoGenomeSyn.pl  Coords2Link   homoMus_outile.filter.coords 2000 homoMus_outile.link  
/shared/ifbstor1/projects/pten_enh_phylo/project_pten_nazim/synteny/NGenomeSyn/bin/NGenomeSyn  -InConf   homoMus_outile.conf   -OutPut    homoMus_outile.svg 
