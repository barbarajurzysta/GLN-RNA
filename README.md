# GLN-RNA

## Usage examples
All RNAs:  
`./countGLN.sh -i downloaded.txt -o out.txt -m -1 -t df -j`

Single RNA:  
`python3 countGLN.py -x pliki/PDBID.xyz -b pliki/PDBID_bonds.csv --out df --maxbonds -1 --name PDBID --notjoined`

## Other files:
* results_no_maxbonds.txt - GLN results with parameter maxbonds=-1 (we don't break the tails because of the bonds)
* analysis_new.ipynb - some plots and statistics of the results
* nrlist_3.235_4.0A.csv - non-redundant list of RNAs from the BGSU database (http://rna.bgsu.edu/rna3dhub/nrlist)
* downloaded.txt - list of downloaded and analysed RNAs (RNAs from the non-redundant list except for the ones with no hydrogen bonds)
