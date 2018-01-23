**MechRNA**: Mechanism Inference for ncRNA
===================
### What is MechRNA?
MechRNA is a computational tool for integration RNA-RNA interactions, RNA-Protein interactions and correlation data in order to infer potential mechanisms for a lncRNA of interest.

### How do I install MechRNA?

1. Ensure you have the required dependancies:
 
- [ViennaRNA](https://www.tbi.univie.ac.at/RNA/) >= 2.4.3
- Python version >= 2.7 with packages:
	- numpy
	- scipy
	- statsmodels

2. Clone or download MechRNA:

```
git clone --recursive https://bitbucket.org/compbio/mechrna.git

```

3. Build [IntaRNA2](https://github.com/BackofenLab/IntaRNA)

```
cd ./mechrna/IntaRNA
bash ./autotools-init.sh
./configure
make

```
Installation is not necessary, MechRNA will use the local executable.

4. Download the required data into /mechrna/data/ which can found [here](). The file structure should be as shown below:

```
/mechrna/data/
/mechrna/data/correlations/
/mechrna/data/peaks/
/mechrna/data/refs/

```


### How do I run MechRNA?
You can use 
```
python mechrna.py -h
```
to get a description of each parameter. For more details, please check doc/MechRNA_manual.pdf (coming soon).


#### Example run with lncRNA 7SL

This is the validation case from Gawronski et. al (DOI: )

To run MechRNA on this lncRNA, type the following when in the mechrna root directory:

```
python mechrna.py -p my_project -l ./sequences/7SL.fa -a -T ./example/7SL.ids -M ./example/7SL.mechs

```

Feel free to change **my_project** to any name you want. The TSV file with the prediction results will be generated in **my_project**.


---


### Contact & Support

Feel free to drop any inquiry to [agawrons at sfu dot ca](mailto:).
