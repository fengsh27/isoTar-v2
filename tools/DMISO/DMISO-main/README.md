# DMISO

Author: Amlan Talukder

Date: August 1, 2021

DMISO is a software used to predict interaction between the provided miRNA and target sequences. 
It was also trained to predict interaction between isomiR and their target sequences. 
It was developed by the computational System biology group at University of Central Florida.


INSTALLATION
--------------------------------------------------------------------------------------------
   1. Install Python 3
   2. Install NumPy, itertools, Keras (2.3.1) packages for Python 3.

EXECUTION 
--------------------------------------------------------------------------------------------------------------------------------------

   You can run the software by running the script "dmiso.py" with the following command:
   
   ----------------------------------------------------------------------------------------
   
	usage: dmiso.py [-h] -p PAIR -m MIRNA/ISOMIR -t MRNA [-o OUTPUT]

	Predicts interactions between miRNA and target pair. It takes an input file of
	query miRNA/isomiR and target pairs and predict the interactions between the
	pairs. Alternatively, it can take separate files for miRNA/isomiR and target
	sequences and outputs all the possible interactions.

	optional arguments:
	  -h, --help       show this help message and exit
	  -o OUTPUT        Path for DMISO outputs

	required arguments:
	  -p PAIR          Path for miRNA (or isomiR) and mRNA pair file in a tsv
		           format. The file must contain data in the following order:
		           miRNA id, target id, miRNA sequence, target sequence
	  -m MIRNA/ISOMIR  Path for miRNA or isomiR sequence in fasta format.
	  -t MRNA          Path for mRNA sequence in fasta format.

	Example: "python3 dmiso.py -p examples/test_pairs.txt -o
	examples/test_output.txt" OR "python3 dmiso.py -m examples/test_miRNAs.fa -t
	examples/test_mRNAs.fa -o examples/test_output.txt"


Required inputs
---------------------------------------------------------------------------------------------
The tool takes a file path as mandatory parameter. The file must contain miRNA id, target id, miRNA sequence and target sequence in the tab delimited format.
Alternatively, it can take separate files for miRNA/isomiR and target sequences and outputs all the possible interactions.

MODEL
----------------------------------------------------------------------------------------------------------------------------------
The model file is stored under "models" directory.


RESULTS
----------------------------------------------------------------------------------------------------------------------------------
The result file contains the miRNA and target information provided in the input file followed by a prediction score and prediction value.

```
miRNA ID	Target ID	miRNA Sequence	Target Sequence	Prediction Score	Prediction
hsa-miR-7111-3p|MIMAT0028120	ENST00000569083.1_3|ENSG00000134419.15_7|RPS15A|protein_coding	UCCUCUUCUCCCUCCUCCCAG	CUUAAUUAAAAGAAGUUAAUGCUAAGAAUUUCUGUGGUGCAGUUUGACUUAAG	0.073714145	0
hsa-miR-4286|MIMAT0016916	ENST00000338086.9_5|ENSG00000055917.16_9|PUM2|protein_coding	ACCCCACUCCUGGUACCA	GUAGAUUAUUGGAAGAUUUCAGAAACAACCGCUUCCCAAACCUUCAGCUUAGA	0.07340459	0
hsa-miR-27a-3p|MIMAT0000084	ENST00000369158.1_2|ENSG00000203811.1_4|H3C14|protein_coding	UUCACAGUGGCUAAGUUCUUCU	CGGCAAGGCCCCGAGGAAGCAGCUGGCCACCAAGGCGGCCCGCAAGAGCGCGCCGG	0.99942267	1
hsa-miR-8089|MIMAT0031016	ENST00000264658.11_7|ENSG00000108306.13_13|FBXL20|protein_coding	UGAGGUCAGGGGAUUGGGAUUC	GGCUCACGCCUGUAAUCCCAGCACUUUGGGAGGCCGAGGCGGGCGGAUCACAAAGUGUCAGGAGUUUGAGAACAG	0.9987809	1
```

LICENSE & CREDITS
-------------------------------------------------------------------------------------------------
The software is freely available for academic use.
plase contact Xiaoman (Shawn) Li (xiaoman@mail.ucf.edu) for further information. 


CONTACT INFO
-------------------------------------------------------------------------------------------------
If you are encountering any problem regarding to DMISO, please refer the manual first.
If problem still can not be solved, please feel free to contact us:
Amlan Talukder (amlan@knights.ucf.edu)
Haiyan (Nancy) Hu (haihu@cs.ucf.edu)
Xiaoman (Shawn) Li (xiaoman@mail.ucf.edu)
