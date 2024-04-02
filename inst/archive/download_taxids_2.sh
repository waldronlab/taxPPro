#! /bin/bash

# I think this is the version (2) that I need rather than download_taxids.sh

# rm -rf *txt
# RUn in terminal, not using Rstudio
exec > >(tee log.txt) 2>&1
DAT=$(date +%m%d%Y)
getTaxID() (
	echo -e "Query (Bacteria/Archaea): \"$2\""
	esearch -db taxonomy -query "Bacteria$2" | efetch -format uid > bac_"$1"_"$DAT".txt
	esearch -db taxonomy -query "Archaea$2" | efetch -format uid > arc_"$1"_"$DAT".txt
	cat bac_"$1"_$DAT.txt arc_"$1"_$DAT.txt | sort | uniq > proc_"$1"_$DAT.txt
	wc -l *$1*.txt | grep -e "$1"
)
echo -e "Downloading taxids..."


echo
fnamekey=species # needs update
query='[SubTree] AND species[Rank] NOT unclassified[prop] NOT uncultured[prop] AND ("above species level"[prop] OR specified[prop])' # needs update
echo -e "Getting species (*$fnamekey*txt files)." # needs update
getTaxID "$fnamekey" "$query"

echo
fnamekey=strains # needs update
query='[SubTree] AND strain[Rank]' # needs update
echo -e "Getting strains (*$fnamekey*txt files)." # needs update
getTaxID "$fnamekey" "$query"






## mv bac*txt arc*txt proc*txt log.txt ../extdata
