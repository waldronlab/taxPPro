#! /bin/bash
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
fnamekey=exclude_all # needs update
query='[SubTree] NOT unclassified[prop] NOT uncultured[prop] AND ("above species level"[prop] OR specified[prop])' # needs update
echo -e "Exclude unclassified, uncultured, and informal names (*$fnamekey*txt files)." # needs update
getTaxID "$fnamekey" "$query"

echo
fnamekey=include_all # needs update
query='[SubTree]' # needs update
echo -e "Include unclassified, uncultured, and informal names (*$fnamekey*txt files)." # needs update
getTaxID "$fnamekey" "$query"

echo
fnamekey=exclude_unclassified # needs update
query='[SubTree] NOT unclassified[prop]' # needs update
echo -e "Exclude unclassified (*$fnamekey*txt files)." # needs update
getTaxID "$fnamekey" "$query"

echo
fnamekey=exclude_uncultured # needs update
query='[SubTree] NOT uncultured[prop]' # needs update
echo -e "Exclude uncultured (*$fnamekey*txt files)." # needs update
getTaxID "$fnamekey" "$query"

echo
fnamekey=exclude_informal # needs update
query='[SubTree] AND ("above species level"[prop] OR specified[prop])' # needs update
echo -e "Exclude informal names (*$fnamekey*txt files)." # needs update
getTaxID "$fnamekey" "$query"

echo
fnamekey=exclude_unclassified_informal # needs update
query='[SubTree] NOT unclassified[prop] AND ("above species level"[prop] OR specified[prop])' # needs update
echo -e "Exclude unclassified and informal names (*$fnamekey*txt files)." # needs update
getTaxID "$fnamekey" "$query"

echo
fnamekey=exclude_unclassified_uncultured # needs update
query='[SubTree] NOT unclassified[prop] NOT uncultured[prop]' # needs update
echo -e "Exclude unclassified and uncultured (*$fnamekey*txt files)." # needs update
getTaxID "$fnamekey" "$query"

echo
fnamekey=exclude_uncultured_informal # needs update
query='[SubTree] NOT uncultured[prop] AND ("above species level"[prop] OR specified[prop])' # needs update
echo -e "Exclude uncultured and informal names (*$fnamekey*txt files)." # needs update
getTaxID "$fnamekey" "$query"

mv bac*txt arc*txt proc*txt log.txt ../extdata
