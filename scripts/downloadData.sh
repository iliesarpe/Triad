!/bin/bash

urls=(
	"https://nrvis.com/download/data/labeled/fb-CMU-Carnegie49.zip"
	"https://snap.stanford.edu/data/soc-pokec-relationships.txt.gz"
	"https://snap.stanford.edu/data/bigdata/communities/com-friendster.ungraph.txt.gz"
	"https://snap.stanford.edu/data/bigdata/communities/com-orkut.ungraph.txt.gz"
	"https://snap.stanford.edu/data/soc-LiveJournal1.txt.gz"
	"https://nrvis.com/download/data/bio/bio-mouse-gene.zip"
	"https://nrvis.com/download/data/graph500/graph500-scale23-ef16_adj.zip"
	"https://snap.stanford.edu/data/gplus_combined.txt.gz"
	"https://nrvis.com/download/data/labeled/proteins-all.zip"
	"https://nrvis.com/download/data/ca/ca-hollywood-2009.zip"
	"https://snap.stanford.edu/data/higgs-social_network.edgelist.gz"
	"https://nrvis.com/download/data/bn/bn-human-BNU_1_0025919_session_1.zip"
	"https://snap.stanford.edu/data/twitch_gamers.zip"
)

output_dir="../data/raw"
mkdir -p "$output_dir"

for url in "${urls[@]}"; do
  echo "Downloading $url..."
  filename=$(basename "$url")
  wget -P "$output_dir" "$url"
done
