1. The script downDblp is used to download and build the raw data needed to construct the temporal DBLP snapshots used for the case study

2. The considered datasets are:
BM: bio-mouse-gene.txt (https://networkrepository.com/bio-mouse-gene.php)
TW: large_twitch_edges.txt(https://snap.stanford.edu/data/twitch_gamers.html)        
SP: soc-pokec-relationships.txt (https://snap.stanford.edu/data/soc-Pokec.html)
BNH: bn-human-BNU_1_0025919_S1.txt (https://networkrepository.com/bn-human-BNU-1-0025919-session-1.php)
fb-CMU: fb-CMU-Carnegie49.txt (https://networkrepository.com/fb-CMU-Carnegie49.php) 
PT: proteins-all.txt (https://networkrepository.com/proteins-all.php) 
YT: soc-YouTube-ASU.txt (https://networkrepository.com/soc-YouTube-ASU.php)
HW: ca-hollywood-2009.txt (https://networkrepository.com/hollywood-2009.php) 
GP: gplus_combined.txt (https://snap.stanford.edu/data/ego-Gplus.html) 
FR: soc-friendster.txt (https://snap.stanford.edu/data/com-Friendster.html)
OR: com-orkut.txt (https://snap.stanford.edu/data/com-Orkut.html) 
HG: higgs-social-followers.txt (https://snap.stanford.edu/data/higgs-twitter.html) 
LJ: soc-LiveJournal1.txt (https://snap.stanford.edu/data/soc-LiveJournal1.html)
G500: graph500-scale23-ef16-adj.txt (https://networkrepository.com/graph500-scale23-ef16-adj.php)

Make sure that they are in txt format such that each line consists of
<src dst>

Then use the driver.py in the main folder to build an HDF5 file for each of such txt file (we will provide a tool in the final version of the current library)

