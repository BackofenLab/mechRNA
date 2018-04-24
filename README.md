**MechRNA**: Mechanism Inference for ncRNA
===================
### What is MechRNA?
MechRNA is a computational tool for integration RNA-RNA interactions, RNA-Protein interactions and correlation data in order to infer potential mechanisms for a lncRNA of interest.

### How do I install MechRNA?

1. Ensure you have the required dependancies:
 
- [IntaRNA2](https://github.com/BackofenLab/IntaRNA) >= 2.2.0
- Python (tested with version 2.7 and 3.5) with packages:
	- numpy
	- scipy
	- statsmodels
	- configparser (Python2 only)

2. Clone or download MechRNA:

		git clone https://bitbucket.org/compbio/mechrna.git
		

3. Download the required data into /mechrna/data/ which can found [here](https://zenodo.org/record/1115534/files/mechrna.data.grch38.tar.gz). The file structure should be as shown below:

		/mechrna/data/
		/mechrna/data/correlation/
		/mechrna/data/peaks/
		/mechrna/data/refs/
	
	For screening mode, please also download the [correlation file](https://zenodo.org/record/1115534) for your cancer type of interest. This file should be placed in the /mechrna/data/correlations/ folder.

### How do I run MechRNA?

You can use 
```
python mechrna.py -h

```

to get a description of each parameter. For more details, please check doc/MechRNA_manual.pdf (coming soon).

To create a new hypothesis-driven mode project: specify (1) project name and (2) lncRNA transcript ID (3) Files listing the targets, rbps, and mechanisms you are intested in (all are used if not specified):
```
python mechrna.py -p my_project -l ENSTXXXXXXXXXXX -a [-T my_targets] [-B my_rbps] [-M my_mechanisms]

```

To create a new screening mode project: specify (1) project name and (2) lncRNA transcript ID (3) correlation file:
```
python mechrna.py -p my_project -l ENSTXXXXXXXXXXX -c ./data/correlation/correl.genenet.[cancer_type].tcga	

```

The hypothesis-driven mode parameters can also be specified for screening mode. Screening mode is computationally expensive and should be executed with multiple workers [--num-worker] and ideally on a cluster [--mode sge/pbs]. Slurm support coming soon!

Feel free to change **my_project** to any name you want. The TSV file with the prediction results will be generated in **my_project**.

#### Example run with lncRNA 7SL

This is the validation case from Gawronski et. al (DOI: )

To run MechRNA on this lncRNA, type the following when in the mechrna root directory:

```
python mechrna.py -p my_project -l ENST00000618786 -a -T ./example/7SL.ids -M ./example/7SL.mechs

```

---


### Contact & Support

Feel free to drop any inquiry to [agawrons at sfu dot ca](mailto:).
