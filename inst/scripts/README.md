These datasets were downloaded on Jan 23, 2023

The source link is: https://www.ncbi.nlm.nih.gov/Taxonomy/taxonomyhome.html/index.cgi?chapter=statistics&period=&from=&to=

Some help: https://www.ncbi.nlm.nih.gov/books/NBK21100/

+ Bacteria[subtree] will retrieve a nonhierarchical list of all of the taxa listed within the Bacteria.
+ The ‘specified’ property flags formal names at and below the species level.

################################################################################
################################################################################
## The searchers were:

## Include all
Bacteria[SubTree]
Archaea[SubTree]

## Exclude all
Bacteria[SubTree] NOT unclassified[prop] NOT uncultured[prop] AND ("above species level"[prop] OR specified[prop])
Archaea[SubTree] NOT unclassified[prop] NOT uncultured[prop] AND ("above species level"[prop] OR specified[prop])

## exclude unclassified
Bacteria[SubTree] NOT unclassified[prop]
Archaea[SubTree] NOT unclassified[prop]

## exclude uncultured
Bacteria[SubTree] NOT uncultured[prop]
Archaea[SubTree] NOT uncultured[prop]

## exclude informal names
Bacteria[SubTree] AND ("above species level"[prop] OR specified[prop])
Archaea[SubTree] AND ("above species level"[prop] OR specified[prop])

## Exclude unclassified and informal names
Bacteria[SubTree] NOT unclassified[prop] AND ("above species level"[prop] OR specified[prop])
Archaea[SubTree] NOT unclassified[prop] AND ("above species level"[prop] OR specified[prop])

## Exclude unclassified and uncultured
Bacteria[SubTree] NOT unclassified[prop] NOT uncultured[prop]
Archaea[SubTree] NOT unclassified[prop] NOT uncultured[prop]

## Exclude uncultured and informal names
Bacteria[SubTree] NOT uncultured[prop] AND ("above species level"[prop] OR specified[prop])
Archaea[SubTree] NOT uncultured[prop] AND ("above species level"[prop] OR specified[prop])

################################################################################
################################################################################



echo 
fnamekey=include_all # needs update
query='[SubTree]' # needs update
echo -e "Include unclassified, uncultured, and informal names (*$fnamekey*txt files)." # needs update
################################################################################
echo -e "Query (Bacteria/Archaea): \"$query\""
esearch -db taxonomy -query "Bacteria$query" | efetch -format uid > bac_"$fnamekey"_"$DAT".txt
esearch -db taxonomy -query "Archaea$query" | efetch -format uid > arc_"$fnamekey"_"$DAT".txt
cat bac_"$fnamekey"_$DAT.txt arc_"$fnamekey"_$DAT.txt | sort | uniq > proc_"$fnamekey"_$DAT.txt
wc -l *$fnamekey*.txt | grep -e "$fnamekey"

echo 
fnamekey=exclude_all # needs update
query='[SubTree] NOT unclassified[prop] NOT uncultured[prop] AND ("above species level"[prop] OR specified[prop])' # needs update
echo -e "Exlude unclassified, uncultured, and informal names (*$fnamekey*txt files)" # needs update
################################################################################
echo -e "Query (Bacteria/Archaea): \"$query\"."
esearch -db taxonomy -query "Bacteria$query" | efetch -format uid > bac_"$fnamekey"_"$DAT".txt
esearch -db taxonomy -query "Archaea$query" | efetch -format uid > arc_"$fnamekey"_"$DAT".txt
cat bac_"$fnamekey"_$DAT.txt arc_"$fnamekey"_$DAT.txt | sort | uniq > proc_"$fnamekey"_$DAT.txt
wc -l *$fnamekey*.txt | grep -e "$fnamekey"

echo 
fnamekey=exclude_unclassified # needs update
query='[SubTree] NOT unclassified[prop]' # needs update
echo -e "Exclude unclassified (*$fnamekey*txt files)." # needs update
################################################################################
echo -e "Query (Bacteria/Archaea): \"$query\""
esearch -db taxonomy -query "Bacteria$query" | efetch -format uid > bac_"$fnamekey"_"$DAT".txt
esearch -db taxonomy -query "Archaea$query" | efetch -format uid > arc_"$fnamekey"_"$DAT".txt
cat bac_"$fnamekey"_$DAT.txt arc_"$fnamekey"_$DAT.txt | sort | uniq > proc_"$fnamekey"_$DAT.txt
wc -l *$fnamekey*.txt | grep -e "$fnamekey"

