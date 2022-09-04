#!/bin/bash
  
function Usage()
{
cat <<EOF
Usage: $0 [-h] -i INPUT_FILE [-o OUTPUT_FILE] [-m MAXBONDS] [-t OUTPUT_TYPE] [-j]
 -h help
 -i INPUT_FILE - path to file with pdb ids and chain names (format <PDBid chain> in each line)
 -o OUTPUT_FILE
 -m MAXBONDS - maxbonds parameter for defining tails, default: 0
 -t OUTPUT_TYPE - min/middle/total, default: total
 -j - count GLN only on continuous fragments of tails
EOF
  
} 
  
# Jesli nie podano argumentow wypisz Usage
# 
if [ "$*" = "" ]
then
    Usage 
    # Zwroc kod bledu
    exit 20
fi
  
# Przygotowanie do pobierania opcji 
set -- `getopt hi:o:m:t:j $*`

while [ "$1" != -- ]
do
    case $1 in
        -h)    HFLG=1;;
        -i)    INPUT_FILE=$2; shift;; 
        -o)    OUTPUT_FILE=$2; shift;; 
        -m)    MAXBONDS=$2; shift;;
        -t)    OUT=$2; shift;;
	-j)    NJ=1;;
    esac
    shift   # nastepna opcja
done
shift  



if [ "$HFLG" = 1 ] 
then
    Usage 
    exit 20
fi


if [ ! "$OUTPUT_FILE" ]; then 
  OUTPUT_FILE="gln_results_no_maxbonds.txt"
fi


if [ ! "$MAXBONDS" ]; then 
  MAXBONDS=0
fi


if [ ! "$OUT" ]; then 
  OUT="total"
fi

if [ "$OUT" == "df" ]; then 
echo "PDB;Loop;BType;Tails;Nmmax;Nclass;Nwhole;Nmax;NmaxTail;NmaxTails;Nmin;NminTail;NminTails;Cmmax;Cclass;Cwhole;Cmax;CmaxTail;CmaxTails;Cmin;CminTail;CminTails;WithGaps" >> $OUTPUT_FILE
fi


IFS=$'\n'       # make newlines the only separator
set -f          # disable globbing
for LINE in $(cat < "$INPUT_FILE"); do
 #RNA=${LINE// /_} - zamiana " " na "_", teraz niepotrzebna
 RNA=$LINE
 echo "$RNA"
 if [ "$NJ" == 1 ]; then
  echo `python3 countGLN.py -x pliki/"$RNA".xyz -b pliki/"$RNA"_bonds.csv --out "$OUT" --maxbonds "$MAXBONDS" --name "$RNA" --notjoined` >> $OUTPUT_FILE
 else
  echo `python3 countGLN.py -x pliki/"$RNA".xyz -b pliki/"$RNA"_bonds.csv --out "$OUT" --maxbonds "$MAXBONDS" --name "$RNA"` >> $OUTPUT_FILE
 fi
done
