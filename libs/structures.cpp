#include "structures.h"


std::vector<std::string> Utilities::splitString(std::string s, std::string delimiter)
{
	size_t pos_start = 0, pos_end, delim_len = delimiter.length();
	std::string token;
	std::vector<std::string> res;

	while ((pos_end = s.find(delimiter, pos_start)) != std::string::npos) 
	{
		token = s.substr (pos_start, pos_end - pos_start);
		pos_start = pos_end + delim_len;
		res.push_back (token);
	}
	res.push_back (s.substr (pos_start));
	return res;
}


std::vector<long long int> Utilities::sublinearBinomialSample(long long int trials, double probability, std::mt19937_64& generator)
{
    long long int position  = 0;
    std::vector<long long int> sampled;
    std::uniform_real_distribution<double> unif(0.0, 1.0);
    double x = unif(generator);
    double c = std::log(1-probability);
    if(x < probability)
    {
        sampled.push_back(position);
    }
    while(position < trials)
    {
        x = unif(generator);
        double y = std::log(x)/c;
        position = position + std::floor(y) + 1;
        if(position < trials)
        {
            sampled.push_back(position);
        }
    }
    return sampled;
}

void Structures::csvToHdf5(std::string fname, std::string outname, bool isColored, bool remap, bool huge){
	std::vector<std::vector<Node_t>> edges;
	std::cout << "in: " << fname << " out: "<< outname << "remap: " << remap << '\n';

	long long int totNodes = 0;
	long long int totEdges = 0;

	std::unordered_map<Node_t, Node_t> mapIDs;
	std::unordered_map<Node_t, Node_t> revmapIDs;
	std::unordered_map<Node_t, std::unordered_set<Node_t>> adjs;
	Node_t ID = 0;
	std::string line{};
    std::ifstream myfile(fname);
    int i{0};
    while(std::getline(myfile, line)){
        std::istringstream iss{line};
		if(!remap)
		{
			if(i==0)
			{
				long long int N, M;
				iss >> N >> M;
				totNodes = N;
				totEdges = M;
			}
			else if(i>0 && i <= totEdges)
			{
				long long int src, dst;
				iss >> src >> dst;
				std::vector<long long int> edge{src, dst};
				edges.push_back(edge);
			}
			else if(isColored && (i>totEdges))
			{
				std::cout << "Do something as graph is colored" << std::endl;
			}
		}
		else
		{
			Node_t src, dst;
			iss >> src >> dst;
			if(src == dst)
				continue;
			if(mapIDs.find(src) == mapIDs.end())
			{
				mapIDs[src] = ID;
				revmapIDs[ID] = src;
				ID++;
			}
			if(mapIDs.find(dst) == mapIDs.end())
			{
				mapIDs[dst] = ID;
				revmapIDs[ID] = dst;
				ID++;
			}
			Node_t srcID = (mapIDs[src] < mapIDs[dst]) ? mapIDs[src] : mapIDs[dst];
			Node_t dstID = (mapIDs[src] < mapIDs[dst]) ? mapIDs[dst] : mapIDs[src];
			if(adjs.find(srcID) != adjs.end())
			{
				if(adjs[srcID].find(dstID) == adjs[srcID].end())
				{
					adjs[srcID].insert(dstID);
					std::vector<long long int> edge{srcID, dstID};
					edges.push_back(edge);
				}
			}
			else
			{
				std::unordered_set<Node_t> neig{dstID};
				adjs[srcID] = neig;
			}
		}
		i++;
    }
	H5Easy::File hdf5File(outname, H5Easy::File::OpenOrCreate);
	H5Easy::DumpOptions options{H5Easy::Compression(), H5Easy::DumpMode::Overwrite};
	if(remap)
	{
		totEdges = edges.size();
		totNodes = ID;
		std::vector<Node_t> origIDs(ID, 0);
		for(int j=0; j< totNodes; j++)
		{
			origIDs[j] = revmapIDs[j];
		}
		std::cout << "dumping original-IDs\n";
		H5Easy::dump(hdf5File, "/original-node-IDs", origIDs, options);
	}
	std::string pathWrite = "/edges/batch-0";
	std::cout << "dumping edges\n";
	size_t MEC{100000000};
	if(edges.size() < MEC)
	{
		H5Easy::dump(hdf5File, pathWrite, edges, options);
		H5Easy::dumpAttribute(hdf5File, pathWrite, "batches", 0, options);
	}
	else
	{
		size_t tot_seen = 0;
		int it{0};
		size_t tot_edges = edges.size();
		while(tot_seen < tot_edges)
		{
			size_t starting_pos{tot_seen};
			size_t ending_pos{std::min(tot_seen+MEC, edges.size())};
			std::vector<std::vector<Node_t>> to_insert;
			to_insert.insert(to_insert.begin(), edges.begin()+starting_pos, edges.begin()+ending_pos);
			std::string currname;
			currname = "/edges/batch-"+std::to_string(it);
				
			H5Easy::dump(hdf5File, currname, to_insert, options);
			tot_seen += to_insert.size();
			it++;
		}
		H5Easy::dumpAttribute(hdf5File, "/edges/batch-0", "batches", it);
	}
	std::cout << "dumping number of nodes\n";
	H5Easy::dumpAttribute(hdf5File, pathWrite, "nodes", totNodes, options);
	std::cout << "dumping number of edges\n";
	H5Easy::dumpAttribute(hdf5File, pathWrite, "edges", totEdges, options);
    return;
}

void Structures::generateDBLPTemporal(std::string outname, json jsonlabs)
{
	std::vector<std::vector<Node_t>> edges;

	H5Easy::File hdf5File(outname, H5Easy::File::OpenOrCreate);
	H5Easy::DumpOptions options{H5Easy::Compression(), H5Easy::DumpMode::Overwrite};


	std::unordered_map<std::string, Node_t> authorMap = jsonlabs["authormapping"].get<std::unordered_map<std::string, Node_t>>();
	std::unordered_map<std::string, std::unordered_map<std::string, std::unordered_map<std::string, long long int>>> authorVenuesYears = jsonlabs["authorsVenuesYears"].get<std::unordered_map<std::string, std::unordered_map<std::string, std::unordered_map<std::string, long long int>>>>();
	std::vector<std::vector<int>> yearRanges = jsonlabs["yearLegend"].get<std::vector<std::vector<int>>>();
	std::unordered_map<std::string, std::unordered_map<std::string, std::vector<Node_t>>> edgesYears = jsonlabs["graphEdgesYear"].get<std::unordered_map<std::string, std::unordered_map<std::string, std::vector<Node_t>>>>();
	std::vector<json> confs;
	for(int y=0; y< yearRanges.size(); y++)
	{
		long long int totNodes = 0;
		long long int totEdges = 0;
		std::cout << "year: " << y << std::endl;
		std::unordered_map<std::string, std::unordered_map<std::string, long long int>> authorVenues = authorVenuesYears[std::to_string(y)];
		std::unordered_map<std::string, std::vector<Node_t>> currentEdges = edgesYears[std::to_string(y)];
		
		std::unordered_map<Node_t, Node_t> mapIDs;
		std::unordered_map<Node_t, Node_t> revmapIDs;
		std::unordered_map<Node_t, std::unordered_set<Node_t>> adjs;
		std::vector<std::vector<long long int>> edges;

		Node_t maxID = 0;
		Node_t currID = 0;
		for(auto& auth : currentEdges)
		{
			auto srcO = stoll(auth.first);
			auto src = 0;
			if(mapIDs.find(srcO) == mapIDs.end())
			{
				mapIDs[srcO] = currID;
				revmapIDs[currID] = srcO;
				currID++;
			}
			src = mapIDs[srcO];
			for(auto dst : auth.second)
			{
				if(mapIDs.find(dst) == mapIDs.end())
				{
					mapIDs[dst] = currID;
					revmapIDs[currID] = dst;
					currID++;
				}
				dst = mapIDs[dst];
				std::vector<long long int> edge{src, dst};
				edges.push_back(edge);
			}
		}
		totNodes = currID;
		totEdges = edges.size();

		std::unordered_map<Node_t, std::string> authRevMap;
		std::vector<std::string> authNames(totNodes);
		for(auto& pair : authorMap)
		{
			if(mapIDs.find(pair.second) != mapIDs.end()) // Author active in current graph
			{
				authRevMap[mapIDs[pair.second]] = pair.first;
				authNames[mapIDs[pair.second]] = pair.first;
			}
		}

		std::unordered_set<std::string> ML_VENUES{"NIPS", "NeurIPS", "ICML", "ICLR", "AAAI", "IJCAI", "AAMAS"};
		std::unordered_set<std::string> THEORY_VENUES{"STOC", "FOCS", "SODA","ITCS", "ICALP", "ISAAC", "ESA", "FOSSACS", "STACS", "COLT", "MFCS", "STACS"};
		std::unordered_set<std::string> DM_CONFERENCES{"KDD", "WWW","CIKM", "SDM", "DAMI", "TKDD", "TKDE", "WSDM", "SIGKDD"};
		std::unordered_set<std::string> DB_CONFERENCES{"VLDB", "PVLDB", "SIGMOD", "PODS", "ICDE", "PODC", "SIGIR"};
		std::unordered_set<std::string> COMP_BIO{"RECOMB", "CIBB", "ISMB", "WABI", "PSB"};
		std::unordered_set<std::string> SOFTWARE{"ASE", "DEBS", "ESEC", "FSE", "ESEM", "ICPE", "ICSE", "ISSTA", "MOBILESoft", "MoDELS", "CBSE", "PASTE", "SIGSOFT"};
		const int DB = 0;
		const int DM = 1;
		const int ML = 2;
		const int TH = 3;
		const int BIO = 4;
		const int SW = 5;
		const int OTH = 6;
		std::vector<std::string> legend{"Data-management","Data-mining", "Machine-learning", "TCS", "Comput. Bio.", "Soft. Eng.", "Other"};
		std::vector<Node_t> authorLabels(totNodes);
		std::vector<int> totLabs(legend.size());
		std::vector<std::vector<Node_t>> parts(legend.size());
		int i = 0;
		for(int i=0; i < totNodes; i++)
		{
			std::string author = authNames[i];
			int ml = 0, theo = 0, dm = 0, db = 0, other = 0, bio=0, sw=0;
			for(auto& venues : authorVenues[author])
			{
				std::string venue = venues.first;
				long long int publications = venues.second;
				std::vector<std::string> chunks = Utilities::splitString(venue, " ");
				for(auto chunk : chunks)
				{
					if(ML_VENUES.find(chunk) != ML_VENUES.end())
					{
						ml += publications;
					}
					else if(THEORY_VENUES.find(chunk) != THEORY_VENUES.end())
					{
						theo += publications;
					}
					else if(DM_CONFERENCES.find(chunk) != DM_CONFERENCES.end())
					{
						dm += publications;
					}
					else if(DB_CONFERENCES.find(chunk) != DB_CONFERENCES.end())
					{
						db += publications;
					}
					else if(COMP_BIO.find(chunk) != COMP_BIO.end())
					{
						bio += publications;
					}
					else if(SOFTWARE.find(chunk) != SOFTWARE.end())
					{
						sw += publications;
					}
					else
					{
						other += publications;
					}
				}
			}
			if((ml+theo+db+dm+bio+sw) > 0) // at least one paper in the considered confs
			{
				int max = std::max({ml, theo, dm, db, bio, sw});
				if(max == db)
					authorLabels[i] = DB;
				else if(max == dm)
					authorLabels[i] = DM;
				else if(max == theo)
					authorLabels[i] = TH;
				else if(max == ml)
					authorLabels[i] = ML;
				else if(max == bio)
					authorLabels[i] = BIO;
				else if(max == sw)
					authorLabels[i] = SW;
			}
			else
			{
				authorLabels[i] = OTH;
			}
			totLabs[authorLabels[i]]++;
			parts[authorLabels[i]].push_back(i);
		}
		int j = 0;
		/*
		for(auto& labcnt : totLabs)
		{
			std::cout << "cathegory: " << j++ << " has count: " << labcnt << '\n';
		}
		*/

		std::cout << "dumping author names\n";
		H5Easy::dump(hdf5File, "authors/names/year"+std::to_string(y), authNames, options);
		std::cout << "dumping author classification\n";
		H5Easy::dump(hdf5File, "authors/labels/year"+std::to_string(y), authorLabels, options);
		H5Easy::dumpAttribute(hdf5File, "authors/labels/year"+std::to_string(y), "totParts", 5,  options);
		H5Easy::dump(hdf5File, "authors/legend", legend, options);
		
		H5Easy::File hdf5File(outname, H5Easy::File::OpenOrCreate);
		H5Easy::DumpOptions options{H5Easy::Compression(), H5Easy::DumpMode::Overwrite};
		std::string pathWrite = "/edges/year" + std::to_string(y);
		std::cout << "dumping edges\n";
		size_t MEC{100000000};
		if(edges.size() < MEC)
		{
			H5Easy::dump(hdf5File, pathWrite, edges, options);
			H5Easy::dumpAttribute(hdf5File, pathWrite, "batches", 0, options);
		}
		else
		{
			size_t tot_seen = 0;
			int it{0};
			size_t tot_edges = edges.size();
			while(tot_seen < tot_edges)
			{
				size_t starting_pos{tot_seen};
				size_t ending_pos{std::min(tot_seen+MEC, edges.size())};
				std::vector<std::vector<Node_t>> to_insert;
				to_insert.insert(to_insert.begin(), edges.begin()+starting_pos, edges.begin()+ending_pos);
				std::string currname;
				currname = "/edges/batch-"+std::to_string(it);
					
				H5Easy::dump(hdf5File, currname, to_insert, options);
				tot_seen += to_insert.size();
				it++;
			}
			H5Easy::dumpAttribute(hdf5File, "/edges/batch-0", "batches", it);
		}
		std::cout << "dumping number of nodes\n";
		H5Easy::dumpAttribute(hdf5File, pathWrite, "nodes", totNodes, options);
		std::cout << "dumping number of edges\n";
		H5Easy::dumpAttribute(hdf5File, pathWrite, "edges", totEdges, options);

		std::cout << "tot nodes: " << totNodes << " tot edges: " << totEdges << std::endl;
   	 	Structures::Graph G{0,0}; //Create emptyGraph
		G.setNodes(totNodes);
		G.setEdges(totEdges);
		for(auto edge : edges)
		{
			G.addEdge(edge[0], edge[1]);
		}
		G.fixDegreeDistrib();
		std::cout << "graph built!\n";

		auto adj = G.GetAdjsFast();
		auto N = G.GetN();
		auto degs = G.getNodeDegrees();
		std::vector<double> normClust(N, -1);
		std::vector<double> normClosure(N, -1);

    	std::vector<double> localTriangleCounts;
		auto mygraphCopy = G;

		for(long long int vertex=0; vertex < N; vertex++)
		{
			auto mydeg = degs[vertex].second;
			// Processing clustCoeff
			if(mydeg > 1)
			{
				double currDen = static_cast<double>(mydeg * (mydeg-1)/2.0);
				normClust[vertex] = currDen;
			}
			// Processing closureCoeff
			if(mydeg > 0)
			{
				long long int pathCount = 0;
				for(auto& nei : adj[vertex])
				{
					pathCount += (degs[nei].second-1);
				}
				if(pathCount > 0)
				{
					normClosure[vertex] = (pathCount/2.);
				}
			}
		}
		std::cout << "computed norm vertices\n";
		Structures::InputParams paramsAppx;
		paramsAppx.epsilon = 0.05;
		paramsAppx.maxwork = 2;
		paramsAppx.eta = 0.01;
		paramsAppx.iterations = 100;
		paramsAppx.steps = 100;
		paramsAppx.qval = 0.01;
		paramsAppx.mode = 100;

		Structures::Statistics edgeSampleStats;
		edgeSampleStats.buckets = parts;
		edgeSampleStats.nodeToBucket = authorLabels;
		edgeSampleStats.seed = 254*(y+1);
		paramsAppx.normMetricVertices = normClust;
		std::vector<double> clusts = CountingAlgs::SimpleAdapt(mygraphCopy, edgeSampleStats, std::ref(paramsAppx));
		mygraphCopy = G;
		paramsAppx.normMetricVertices = normClosure;
		std::vector<double> closures = CountingAlgs::SimpleAdapt(mygraphCopy, edgeSampleStats, std::ref(paramsAppx));
		json currJ = {
			{"legend", legend},
			{"year", yearRanges[y]},
			{"labels", totLabs},
			{"nodes", totNodes},
			{"edges", totEdges},
			{"clustering", clusts},
			{"closure", closures}
		};
		confs.push_back(currJ);
	}
	json outJson = {
		{"years", yearRanges},
		{"conf-years", confs}
	};
	std::ofstream outj1;
	outj1.open("results/casestudy/dblptemporal.json");
	outj1 << outJson << '\n';
	outj1.close();
	return;
}

void Structures::DumpExactCoeffs(H5Easy::File& datafile, Structures::Graph& G, std::vector<double>& triCounts){
	auto adj = G.GetAdjsFast();
	auto N = G.GetN();
    auto degs = G.getNodeDegrees();
	std::vector<double> normClust(N, -1);
	std::vector<double> clustVals(N, 0);
	std::vector<double> normClosure(N, -1);
	std::vector<double> closureVals(N, 0);

	for(long long int i=0; i < N; i++)
	{
		auto mydeg = degs[i].second;
		// Processing clustCoeff
        if(mydeg > 1)
        {
            double currDen = static_cast<double>(mydeg * (mydeg-1)/2.0);
            clustVals[i] = triCounts[i]/currDen;
			normClust[i] = currDen;
        }
		// Processing closureCoeff
        if(mydeg > 0)
		{
			long long int pathCount = 0;
			for(auto& nei : adj[i])
			{
				pathCount += (degs[nei].second-1);
			}
			if(pathCount > 0)
			{
				closureVals[i] = (2*triCounts[i])/(1.* pathCount);
				normClosure[i] = (pathCount/2.);
			}
		}
	}

	H5Easy::DumpOptions options{H5Easy::Compression(), H5Easy::DumpMode::Overwrite};
	std::string pathWrite = "/clusteringCoeffs";
	std::cout << "dumping clust vals\n";
	H5Easy::dump(datafile, pathWrite + "/values", clustVals, options);
	H5Easy::dump(datafile, pathWrite + "/normNodes", normClust, options);
	std::cout << "dumping closure vals\n";
	pathWrite = "/closureCoeffs";
	H5Easy::dump(datafile, pathWrite + "/values", closureVals, options);
	H5Easy::dump(datafile, pathWrite + "/normNodes", normClosure, options);
    return;
}

void Structures::DumpPartitions(H5Easy::File& datafile, Structures::Graph& G, Structures::PartitionParams& params){
	auto adj = G.GetAdjsFast();
	auto N = G.GetN();
    auto degs = G.getNodeDegrees();
	auto maxDeg = G.getMaxDegree();
	auto totBucketsDegs = params.numberDegreeBuckets;
	H5Easy::DumpOptions options{H5Easy::Compression(), H5Easy::DumpMode::Overwrite};
	std::string pathWrite = "/vertexSetPartitions";
	if(totBucketsDegs > 0) // Compute partition based on bucketed Degree
	{
		int stepDeg = ceil(1.*maxDeg/totBucketsDegs);
		std::vector<int> partition(N, 0);
		std::vector<std::vector<long int>> degreeBucketNodes(totBucketsDegs); // Materialize buckets
		
		for(int n=0; n < N; n++)
		{
			auto mydeg = degs[n].second;
			if(mydeg > 1)
			{
				auto mydegIdx = ceil(mydeg/stepDeg);
				if(mydegIdx > 0) // indices go from 0,...,totBucketsDegs-1
					mydegIdx--;
				partition[n] = mydegIdx;
			}
		}
		std::cout << "dumping partition degree-length\n";
		H5Easy::dump(datafile, pathWrite + "/degreeLengthBased", partition, options);
		H5Easy::dumpAttribute(datafile, pathWrite + "/degreeLengthBased", "totParts", totBucketsDegs, options);
	}
	
	int totCoresPart = params.numberCoreBuckets;
	if(totCoresPart > 0)
	{
		std::vector<std::vector<NodeType>> cores(maxDeg);
		CountingAlgs::GetCoreDecomposition(G, std::ref(cores));
		int cv = 0; // compute max-core number
		for(auto core : cores)
		{
			if(core.size() > 0)
			{
				cv++;
			}
		}
		std::vector<int> partition(N, 0);
		int stepCore = ceil(1.*cv/totCoresPart);
		for(int c=0; c< cv; c++)
		{
			auto myCoreIdx = ceil(1.*c/stepCore);
			if(myCoreIdx > 0) // indices go from 0,...,totBucketsDegs-1
				myCoreIdx--;
			for(auto& node: cores[c])
			{
				partition[node] = myCoreIdx;
			}
		}
		std::cout << "dumping partition core-based\n";
		H5Easy::dump(datafile, pathWrite + "/coreBased", partition, options);
		H5Easy::dumpAttribute(datafile, pathWrite + "/coreBased", "totParts", totCoresPart, options);
	}
	
	std::cout << "computing partition based on log bins " << '\n';
	int tau = 2; //as form the original paper
	int omega = 2;
	std::vector<int> logPartition(N, 0);
	int maxVal=0;
	for(int n=0; n < N; n++)
	{
		auto mydeg = degs[n].second;
		if(mydeg > 1)
		{
			int mypart = 1;
			if(mydeg <= tau)
				mypart = mydeg;
			else
			{
				mypart = std::floor(std::log(1+((omega-1)*(mydeg-tau)))/std::log(omega))+ tau;
			}
			logPartition[n] = mypart-1;
			maxVal = std::max(maxVal, mypart);
		}
	}
	std::cout << "dumping partition log-deg-based\n";
	H5Easy::dump(datafile, pathWrite + "/logDegBased", logPartition, options);
	H5Easy::dumpAttribute(datafile, pathWrite + "/logDegBased", "totParts", maxVal, options);

	if(params.labels)
	{
		std::string pathFileLabels = params.pathLabels;
		std::vector<long long int> nodeLabs(N, 0);
		std::vector<Node_t> nodeIds = H5Easy::load<std::vector<Node_t>>(datafile, "/original-node-IDs");
		std::unordered_map<Node_t, Node_t> mapIds;
		for(Node_t i=0; i < nodeIds.size(); i++)
		{
			mapIds[nodeIds[i]] = i;
		}

		std::string line{};
		std::cout << "path: " << pathFileLabels << '\n';
		std::ifstream mylabfile(pathFileLabels);
		std::unordered_map<int, int> remapLabels;
		int SLAB = 0;
		long long int maxLab = 0;
		while(std::getline(mylabfile, line))
		{
			std::istringstream iss{line};
			long long int ID, label;
			iss >> ID >> label;
			if(remapLabels.find(label) == remapLabels.end())
				remapLabels[label] = SLAB++;
			//maxLab = std::max(maxLab, label);
			nodeLabs[mapIds[ID]] = remapLabels[label];
		}
		std::vector<int> maps(SLAB);
		for(auto& el : remapLabels)
			maps[el.second] = el.first;

		std::cout << "dumping original labels\n";
		H5Easy::dump(datafile, pathWrite + "/originalLabels", nodeLabs, options);
		H5Easy::dumpAttribute(datafile, pathWrite + "/originalLabels", "totParts", SLAB, options);
		H5Easy::dump(datafile, "/mappings/labels", maps, options);
	}

    return;
}

void Structures::BuildGraph(H5Easy::File& datafile, Structures::Graph& G, bool oldread){
	std::vector<std::vector<Node_t>> edges; //= H5Easy::load<std::vector<std::vector<long long int>>>(datafile, "edges");
	long long int M;
	long long int N;
	if(oldread)
	{
		edges = H5Easy::load<std::vector<std::vector<Node_t>>>(datafile, "edges");
		M = H5Easy::loadAttribute<long long int>(datafile, "edges", "edges");
		N = H5Easy::loadAttribute<long long int>(datafile, "edges", "nodes");
	}
	else
	{
		edges = H5Easy::load<std::vector<std::vector<Node_t>>>(datafile, "edges/batch-0");
		M = H5Easy::loadAttribute<long long int>(datafile, "edges/batch-0", "edges");
		if(edges.size() != M)	
		{
			int batches = H5Easy::loadAttribute<int>(datafile, "edges/batch-0", "batches");
			std::cout << "tot_batches: " << batches << '\n';
			for(int b=1; b < batches; b++)
			{
				std::string currData ="edges/batch-" + std::to_string(b);
				std::vector<std::vector<Node_t>> currEdges = H5Easy::load<std::vector<std::vector<Node_t>>>(datafile, currData);
				edges.insert(edges.end(), currEdges.begin(), currEdges.end());
				currEdges.clear();
			}
		}
		N = H5Easy::loadAttribute<long long int>(datafile, "edges/batch-0", "nodes");
		if(edges.size() !=M)
		{
			std::cout << "Some error occurred, edges size do not match!" << std::endl;
			exit(0);
		}
	}
	std::cout << "edges has size: " << edges.size() << " while M is: " << M << std::endl;
	

	G.setNodes(N);
	G.setEdges(M);
	for(auto& edge : edges)
	{
		G.addEdge(edge[0], edge[1]);
	}
    G.fixDegreeDistrib();
    return;
}

void Utilities::dumpTxt(H5Easy::File& datafile, std::string outfile, bool oldread){
	std::vector<std::vector<Node_t>> edges; //= H5Easy::load<std::vector<std::vector<long long int>>>(datafile, "edges");
	long long int M;
	long long int N;
	if(oldread)
	{
		edges = H5Easy::load<std::vector<std::vector<Node_t>>>(datafile, "edges");
		M = H5Easy::loadAttribute<long long int>(datafile, "edges", "edges");
		N = H5Easy::loadAttribute<long long int>(datafile, "edges", "nodes");
	}
	else
	{
		edges = H5Easy::load<std::vector<std::vector<Node_t>>>(datafile, "edges/batch-0");
		M = H5Easy::loadAttribute<long long int>(datafile, "edges/batch-0", "edges");
		if(edges.size() != M)	
		{
			int batches = H5Easy::loadAttribute<int>(datafile, "edges/batch-0", "batches");
			std::cout << "tot_batches: " << batches << '\n';
			for(int b=1; b < batches; b++)
			{
				std::string currData ="edges/batch-" + std::to_string(b);
				std::vector<std::vector<Node_t>> currEdges = H5Easy::load<std::vector<std::vector<Node_t>>>(datafile, currData);
				edges.insert(edges.end(), currEdges.begin(), currEdges.end());
				currEdges.clear();
			}
		}
		N = H5Easy::loadAttribute<long long int>(datafile, "edges/batch-0", "nodes");
		if(edges.size() !=M)
		{
			std::cout << "Some error occurred, edges size do not match!" << std::endl;
			exit(0);
		}
	}
	std::cout << "edges has size: " << edges.size() << " while M is: " << M << std::endl;
	std::ofstream outf;
	outf.open(outfile);
	outf << N << " " << M << '\n';
	for(auto& edge : edges)
	{
		outf << edge[0] << " " << edge[1] << '\n';
	}
	outf.close();
}

std::vector<std::vector<Node_t>> Structures::LoadPartition(H5Easy::File& datafile, std::string path)
{
	std::vector<int> assignment = H5Easy::load<std::vector<int>>(datafile, path);
	int numParts = H5Easy::loadAttribute<int>(datafile, path, "totParts");
	std::vector<std::vector<Node_t>> parts(numParts);
	for(int i=0; i< assignment.size(); i++)
	{
		parts[assignment[i]].push_back(i);
	}
	return parts;
}


std::vector<double> Structures::LoadExactSolPartition(H5Easy::File& datafile, std::string path, std::vector<double>& ccVals){
	std::vector<int> assignment = H5Easy::load<std::vector<int>>(datafile, path);
	int numParts = H5Easy::loadAttribute<int>(datafile, path, "totParts");
	std::vector<double> partVals(numParts, 0);
	std::vector<long long int> partSizes(numParts, 0);
	for(int i=0; i< assignment.size(); i++)
	{
		partVals[assignment[i]] += ccVals[i];
		partSizes[assignment[i]] += 1;

	}
	for(int j=0; j < numParts; j++)
	{
		if(partSizes[j] >0)
			partVals[j] /= partSizes[j];
		//std::cout << "Partition i: " << j << " has nodes: " << partSizes[j] << " and value: " << partVals[j] << std::endl;
	}
	return partVals;
}

std::vector<long long int> Structures::LoadPartitionSize(H5Easy::File& datafile, std::string path){
	std::vector<int> assignment = H5Easy::load<std::vector<int>>(datafile, path);
	int numParts = H5Easy::loadAttribute<int>(datafile, path, "totParts");
	std::vector<long long int> partSizes(numParts, 0);
	for(int i=0; i< assignment.size(); i++)
	{
		partSizes[assignment[i]] += 1;

	}
	return partSizes;
}


void Structures::BuildGraph(std::string fname, Structures::Graph& G){
    std::string line{};
    std::ifstream myfile(fname);
    int i{0};
    int edges{0};
    //std::cout << "Begin init: " << std::endl;
    while(std::getline(myfile, line)){
        std::istringstream iss{line};
        if(i==0)
        {
            int N, M;
            iss >> N >> M;
            G.setNodes(N);
            G.setEdges(M);
            edges = M;
        }
        else if(i>0 && i <= edges)
        {
            int src, dst;
            iss >> src >> dst;
            G.addEdge(src, dst);
        }
        else if(i>edges)
        {
            int id, col;
            iss >> id >> col;
            Node curnode{id, col};
            G.addNode(curnode);
        }
        i++;
    }
    G.fixDegreeDistrib();
    return;
}

std::vector<std::pair<long long int, long long int>> CountingAlgs::GetNodeDegrees(Structures::Graph G)
{
        std::vector<std::pair<long long int, long long int>> degrees(G.GetN()); 
        auto adjlist = G.GetAdjs();
        std::vector<Structures::Node> allNodes= G.GetNodes();
        int i{0};
        for(auto adj : adjlist)
        {
            for(auto node : adj)
            {
                if(allNodes[node].getColor() == 0)
                    degrees[i].first++;
                else
                    degrees[i].second++;
            }
            i++;
        }
        return degrees;
}

std::vector<double> CountingAlgs::clustCoeffsWithChiba(Structures::Graph G, Structures::Statistics& mystats)
{
    // Simpmle implementation of the Chiba and Nishizeki algorithm 
    auto tStart = std::chrono::steady_clock::now();
    auto n = G.GetN();
    auto adjlist = G.GetAdjsFast();
    std::vector<double> clust(n);
    std::vector<bool> alive(n, true);
    std::vector<bool> marked(n, false);
    auto sortedDegs = G.getSortedUncoloredDegrees(); 
    for(long long int i = 0; i < n; i++)
    {
        for(auto nei : adjlist[sortedDegs[i].first])
        {
            if(alive[nei])
			{
                marked[nei] = true;
			}
        }
        for(auto nei : adjlist[sortedDegs[i].first])
        {
            if((marked[nei]) && (alive[nei]))
            {
                for(auto nei1 : adjlist[nei])
                {
                    if((marked[nei1]) && (alive[nei1]) && (nei1 != sortedDegs[i].first))
                    {
                        clust[sortedDegs[i].first]++;
                        clust[nei]++;
                        clust[nei1]++;
                    }
                }
            }
            marked[nei] = false;
        }
        alive[sortedDegs[i].first] = false;
    }

    auto tFinish = std::chrono::steady_clock::now();
    std::chrono::duration<double> totRT = tFinish-tStart;
	mystats.runtime.push_back(totRT.count());

    return clust;
}

std::tuple<int, std::vector<double>> CountingAlgs::SmallTriangleCount(Structures::Graph& G, int threshold, std::vector<bool>& isSmall, int ub)
{

    Node_t nodes = G.GetCurrN();
    std::vector<double> triCountsSmall(nodes, 0);
    long long int removed{0};
    auto bucketDegs =  G.MaterializeDegreeBuckets();
    long long int totRemoved = 0;
    for(auto node : bucketDegs[1]) // Nodes with degree 1
    {
            G.removeNode(node);
            isSmall[node] = true;
            removed++;
    }
    
    bool firstIt = true;
    // TODO: fix this dor when 4*numer of nodes with deg 2 > maxAllowedWork
    int boundDegree = 2;
    int cntIts = 1;
    double maxAllowedWork = static_cast<double>(threshold*nodes);
    // Iteratively remove all nodes with degree less than beta
    while(true)
    {
        auto minDegAlive = G.getMinDegree();
        long long int maximumDegree = G.getMaxDegree();
                
        // Find threshold on the degrees for small nodes
        double totalWork = 0;

        // degreeDist has maxDeg+1 entries, from 0,1,...,maxDeg. Nodes with deg 0 and 1 at this point have CC = 0
        // The only work we need to do is for nodes with deg >= 2
        auto degreeDist = G.getDegreeDistrib();
        double currentBoundDegree = 1;
        if(firstIt)
        {
            if(ub==-1)
            {
                while((totalWork < maxAllowedWork) && (currentBoundDegree <= maximumDegree) )
                {
                    currentBoundDegree++;
                    totalWork += std::pow(currentBoundDegree, 2)*degreeDist[currentBoundDegree];
                }
                currentBoundDegree--;
            }
            else
            {
                    currentBoundDegree =ub;
            }
            if(firstIt)
            {
                boundDegree = currentBoundDegree;
                firstIt = false;
            }
        }
        else
            break;
        if((boundDegree < 2) || (minDegAlive > boundDegree))
            break;
        auto nodeDegrees = G.getNodeDegrees();
        auto adjlistfast = G.GetAdjsFast();
        auto bucketDegs =  G.MaterializeDegreeBuckets();
        for(int i=2; i <= boundDegree; i++)
        {
            for(auto& currNodeID : bucketDegs[i])
            {
                double currTri=0;
                for(auto& nei1 : G.GetNodeNeighborhood(currNodeID))
                {
                    if(!G.isNodeAlive(nei1))
                        continue;
                    auto mindegNode = (nodeDegrees[nei1].second >= nodeDegrees[currNodeID].second)? currNodeID : nei1;
                    auto maxdegNode = (nodeDegrees[nei1].second < nodeDegrees[currNodeID].second)? currNodeID : nei1;
                    for(auto& nei2 : G.GetNodeNeighborhood(mindegNode))
                    {
                        if(!G.isNodeAlive(nei2))
                            continue;
                        if(nei1 >= nei2 || (currNodeID==nei2)) // Avoid counting each triangle twice
                            continue;
                        if(adjlistfast[maxdegNode].find(nei2) != adjlistfast[maxdegNode].end())
                        {
                            triCountsSmall[nei1]+=1;
                            triCountsSmall[nei2]+=1;
                            currTri +=1;
                        }
                    }
                }
                triCountsSmall[currNodeID] += currTri;
                G.removeNode(currNodeID);
                isSmall[currNodeID] = true;
                removed++;
            }
        }
        // Remove also nodes with degree 1
        for(auto& currNodeID : bucketDegs[1])
        {
            G.removeNode(currNodeID);
            isSmall[currNodeID] = true;
            removed++;
        }
		for(auto& currNodeID : bucketDegs[0])
        {
            G.removeNode(currNodeID);
            isSmall[currNodeID] = true;
            removed++;
        }
    }
    return std::make_tuple(boundDegree, triCountsSmall);
}

int CountingAlgs::GetCoreDecomposition(Structures::Graph G, std::vector<std::vector<NodeType>>& cores, long long int N)
{
    // Very naive implementation not O(n+m) !!
    auto degreeBuckets = G.MaterializeDegreeBucketsSets();
    auto nodeDegrees = G.getNodeDegrees();
    long long int nodes = G.GetN();
    if(N != -1)
        nodes = N;
    auto currAdjs = G.GetAdjsFast();
    std::vector<bool> decompAlive(nodes, true);
    int k = 0;
    for(int i=0; i < nodes; i++)
    {
        if(!G.isNodeAlive(i))
            continue;
        int firstNonEmpty=0;
        while(degreeBuckets[firstNonEmpty].size() == 0)
        {
            firstNonEmpty++;
        }
        auto it = degreeBuckets[firstNonEmpty].begin();
        NodeType toRemove = (*it);
        degreeBuckets[firstNonEmpty].erase(it);

        for(auto& neigh : currAdjs[toRemove])
        {
            if(!decompAlive[neigh])
                continue;
            degreeBuckets[nodeDegrees[neigh].second].erase(neigh);
            nodeDegrees[neigh].second-=1;
            degreeBuckets[nodeDegrees[neigh].second].insert(neigh);
        }
        k = std::max(k, firstNonEmpty);
        cores[k].push_back(toRemove);
        decompAlive[toRemove] = false;
    }
    return k;
}

std::vector<double> CountingAlgs::SimpleAdapt(Structures::Graph G, Structures::Statistics& mystats, Structures::InputParams& params)
{

	double epsilon = params.epsilon; 
	double eta = params.eta; 
	int iterations = params.iterations; 
	double qval = params.qval;
	int mode = params.mode;
	long long int samplesI = params.samples;
	long long int maxwork = params.maxwork;
	std::vector<double> normVertices = params.normMetricVertices;
		
	std::vector<double> epsilons(mystats.buckets.size(), epsilon);
	if(params.nonUni)
	{
		epsilons = params.nonUniEps;
		epsilon = *std::max_element(epsilons.begin(), epsilons.end());
	}
    
    auto MD = G.getMaxDegree();
    auto n0 = G.GetCurrN();
    auto m0 = G.GetM();
    std::vector<double> degreeSquared(n0, 0.0); // Denominator in local clustering coefficient cc(v), d_v * (d_v - 1)/2, only with d_v > 1
	degreeSquared = normVertices;
    auto tSTART = std::chrono::steady_clock::now();

    long long int totn = 0;
    auto t0 = std::chrono::steady_clock::now();
    auto t1 = std::chrono::steady_clock::now();
    std::chrono::duration<double> tot = t1-t0;

    auto t0Mid = std::chrono::steady_clock::now();
    std::vector<std::pair<NodeType, long long int>> originalDegrees = G.getNodeDegrees();

    std::vector<bool> smallBool(n0, false);

    const int MAX_WORK_SMALL = maxwork; // constant for removing small degree nodes, O(MAX_WORK_SMALL * |V|)

    auto t0A = std::chrono::steady_clock::now();
    auto [boundDeg ,smallCountsTriangles] = SmallTriangleCount(std::ref(G), MAX_WORK_SMALL, std::ref(smallBool), -1);
    auto t1A = std::chrono::steady_clock::now();
    std::chrono::duration<double> totA = t1A-t0A;
	mystats.timepeel = totA.count();

    auto n1 = G.GetCurrN();
    auto m1 = G.GetM();
	mystats.n1= n1;
	mystats.m1= m1;
    auto nodeDegrees = G.getNodeDegrees();
    std::vector<double> retCC(n0,0.);
	std::vector<double> smallBucketCountsNodes(mystats.buckets.size(),0);
	// The following updates the counts over each partition according to the "small-degree" nodes that are removed
	std::vector<long long int> nodesActiveInParts(mystats.buckets.size(),0);
    t0A = std::chrono::steady_clock::now();
    for(int i=0; i< smallCountsTriangles.size(); i++)
    {
        if(degreeSquared[i]>0)
		{
			if(!smallBool[i])
			{
				nodesActiveInParts[mystats.nodeToBucket[i]]++;
			}
            retCC[i] = smallCountsTriangles[i]/degreeSquared[i];
			double mysmallcc = (1.*smallCountsTriangles[i])/degreeSquared[i];
			smallBucketCountsNodes[mystats.nodeToBucket[i]] += mysmallcc/(1.*mystats.buckets[mystats.nodeToBucket[i]].size());
		}
        
    }
    auto t1Mid = std::chrono::steady_clock::now();
    std::chrono::duration<double> totMid = t1Mid-t0Mid;
	mystats.singleTimes.push_back(totMid.count());


    t0Mid = std::chrono::steady_clock::now();
    std::mt19937_64 gen64(mystats.seed);
    auto currAdjs = G.GetAdjsFast();
    auto edgeList = G.getEdgeList();
    double samplingProb{1.0/edgeList.size()};
    double Q{qval};
    int idx = 0;
	int bucketidx=0;
	int samplesOptVar;
	if(!(mode == MODEFIXEDSS))
	{
		samplesOptVar = 500;
	}
	else
	{
		samplesOptVar = static_cast<int>(samplesI/10.);
	}

    std::uniform_int_distribution<> uniformEdges(0, edgeList.size()-1);

    std::vector<double> edgeWeights(edgeList.size(), 1);
    std::discrete_distribution<> imp(edgeWeights.begin(), edgeWeights.end());
    std::vector<double> probImp = imp.probabilities();

	GRBEnv env = GRBEnv(true);
	env.set("LogFile", "grblog.log");
	env.start();
	GRBModel model = GRBModel(env);
	// Model minimax problem
	GRBVar qModel = model.addVar(0.0, 0.5, 0.0, GRB_CONTINUOUS, "q");
	GRBVar tMax = model.addVar(0.0, m1/samplesOptVar, 0.0, GRB_CONTINUOUS, "t");

	model.setObjective(1*tMax, GRB_MINIMIZE);



	std::vector<long long int> edgesOptVar;
	for(int j=0; j<samplesOptVar; j++)
	{
		edgesOptVar.push_back(imp(gen64));
	}
	std::vector<std::vector<double>> atermEdge(mystats.buckets.size(), std::vector<double>(samplesOptVar,0.0));
	std::vector<std::vector<double>> btermEdge(mystats.buckets.size(), std::vector<double>(samplesOptVar,0.0));
	std::vector<std::vector<double>> ahatEdge(mystats.buckets.size(), std::vector<double>(samplesOptVar,0.0));
	std::vector<std::vector<double>> bhatEdge(mystats.buckets.size(), std::vector<double>(samplesOptVar,0.0));
	std::vector<double> aBar(mystats.buckets.size(), 0);
	std::vector<double> bBar(mystats.buckets.size(), 0);
	int idxEdgeVar =0;
	for(auto& edgePos : edgesOptVar)
	{
		auto currEdge = edgeList[edgePos];
		auto v1 = currEdge.getSRC();
		auto v2 = currEdge.getDST();

		auto degv1 = nodeDegrees[v1].second;
		auto degv2 = nodeDegrees[v2].second;

		auto mindegnode = (degv1 >= degv2) ? v2 : v1;
		auto maxdegnode = (degv1 < degv2) ? v2 : v1;

		auto edgeProb = probImp[edgePos];

		long long int tris=0;
		for(auto& v3 : currAdjs[mindegnode])
		{
			if(currAdjs[v3].find(maxdegnode) !=currAdjs[v3].end())
			{
				tris++;
				atermEdge[mystats.nodeToBucket[v3]][idxEdgeVar] += (1./degreeSquared[v3]);
				btermEdge[mystats.nodeToBucket[v3]][idxEdgeVar] -= (2./degreeSquared[v3]);
			}
		}
		if(tris > 0)
		{
			btermEdge[mystats.nodeToBucket[v1]][idxEdgeVar] -= (1.*tris/degreeSquared[v1]);
			btermEdge[mystats.nodeToBucket[v2]][idxEdgeVar] -= (1.*tris/degreeSquared[v2]);
		}
		idxEdgeVar++;
	}
	for(int i=0; i< atermEdge.size(); i++)
	{
		for(int j=0; j < atermEdge[i].size(); j++)
		{
			aBar[i] += atermEdge[i][j];
			bBar[i] += btermEdge[i][j];
		}
	}
	for(int i=0; i < atermEdge.size(); i++)
	{
		double meanA = aBar[i]/(1.*samplesOptVar);
		double meanB = bBar[i]/(1.*samplesOptVar);
		for(int j=0; j < atermEdge[i].size(); j++)
		{
			ahatEdge[i][j] = atermEdge[i][j] - meanA;
			bhatEdge[i][j] = btermEdge[i][j] - meanB;
		}
	}
	for(int i=0; i< atermEdge.size(); i++)
	{
		if((mystats.buckets[i].size()) > 0 && (nodesActiveInParts[i] > 0))
		{
			double A = 0;
			double B = 0;
			double C = 0;
			for(int j=0; j < atermEdge[i].size(); j++)
			{
				A += std::pow(ahatEdge[i][j], 2);
				B += std::pow(bhatEdge[i][j], 2);
				C += ahatEdge[i][j] * bhatEdge[i][j];
			}
			C *= 2;
			double norm = (1./samplingProb) * (1./mystats.buckets[i].size());
			norm = std::pow(norm, 2);
			norm /= (samplesOptVar-1.0);
			A = A*norm;
			B = B*norm;
			C = C*norm;
			if((std::fpclassify(A) == FP_ZERO) && (std::fpclassify(B) == FP_ZERO) && (std::fpclassify(C) == FP_ZERO))
				continue;
			model.addQConstr(1*tMax >= A + B*qModel*qModel + C*qModel);
		}
	}
	model.optimize();
	double minVar = model.get(GRB_DoubleAttr_ObjVal);
	double bestQMinVar = qModel.get(GRB_DoubleAttr_X);
	mystats.optimizedQ = bestQMinVar;
	params.qopt = bestQMinVar;
		
	Q = bestQMinVar;
	if(mode==MODEQVAL)
	{
		Q =qval;
	}
    double CQ{(1.0-(2.0*Q))};
    t1A = std::chrono::steady_clock::now();
    totA = t1A-t0A;
	mystats.timevar =  totA.count();
		
    t0A = std::chrono::steady_clock::now();
	std::vector<double> upperBounditEdges(mystats.buckets.size(), 0);
	double maxRangeB = 0;
	long long int hoeffdingSamp =0;
	double pseudoDimBound = 0;
	long long int pdimSamp =0;
	double curr_eta = eta;
	if(!(mode == MODEFIXEDSS))
	{
		auto t0VectComp = std::chrono::steady_clock::now();
		long long int boundPseudoDimension = 0;
		std::vector<std::vector<double>> nodeDistribBuckets(n0);
		std::vector<long long int> distinctBucketsNode(n0, 0);

		for(int i=0; i < n0; i++)
		{
			std::vector<bool> distinctBuckets(mystats.buckets.size(), false);
			if(smallBool[i]) // Node is removed cannot contribute further to cc value
				continue;
			nodeDistribBuckets[i].resize(mystats.buckets.size(), 0);
			for(auto& nei : currAdjs[i])
			{
				if(degreeSquared[nei] > 0)
				{
					double degSqNei = degreeSquared[nei];
					nodeDistribBuckets[i][mystats.nodeToBucket[nei]] += CQ/degSqNei;
					distinctBuckets[mystats.nodeToBucket[nei]] = true;
				}
			}
			distinctBuckets[mystats.nodeToBucket[i]] = true;
			distinctBucketsNode[i] = std::count(distinctBuckets.begin(), distinctBuckets.end(), true);
		}
		auto t1VectComp = std::chrono::steady_clock::now();
		std::chrono::duration<double> totVectComp = t1VectComp-t0VectComp;
		mystats.singleTimes.push_back(totVectComp.count());
		t0VectComp = std::chrono::steady_clock::now();

		std::vector<long long int> maxMinDegs;
		std::map<long long int, long long int> degreemaxDistrib;
		std::vector<double> degreeDistribMax(mystats.buckets.size(), 0);

		double minThresDeg = (boundDeg+1) * boundDeg/2.0;

		int zerocounter = 0;
		bucketidx=0;

		long long int totEdgesHeavy0 = 0;
		long long int avgUb = 0;
		double maxRangeInner = 0;

		long long int edgecnt = 0;
		for(auto& curredge : edgeList)
		{
			auto srcNode = curredge.getSRC();
			auto dstNode = curredge.getDST();

			auto srcDeg = nodeDegrees[srcNode].second;
			auto dstDeg = nodeDegrees[dstNode].second;

			auto minDegEdge = std::min(srcDeg, dstDeg);
			auto minDegNode = (srcDeg < dstDeg) ? srcNode : dstNode;

			long long int minDeg = std::min(srcDeg, dstDeg);

			if(minDeg < 2) // current edge cannot close any triangles
				continue;
			double srcDegSq = degreeSquared[srcNode];
			double dstDegSq = degreeSquared[dstNode];

			double qWeight = Q*minDeg;
			nodeDistribBuckets[srcNode][mystats.nodeToBucket[srcNode]] += qWeight/srcDegSq;
			nodeDistribBuckets[dstNode][mystats.nodeToBucket[dstNode]] += qWeight/dstDegSq;
			boundPseudoDimension = std::max(distinctBucketsNode[minDegNode], boundPseudoDimension);

			for(int j=0; j < mystats.buckets.size(); j++)
			{
				if(mystats.buckets[j].size() > 0) // Partition is non-empty
				{

					double bucketUb = nodeDistribBuckets[minDegNode][j];
					bucketUb *= (1.0/(samplingProb));
					bucketUb *= (1.0/mystats.buckets[j].size());
					upperBounditEdges[j] = std::max(bucketUb, upperBounditEdges[j]); // Non-uniform bound on the range
				}
			}
			nodeDistribBuckets[srcNode][mystats.nodeToBucket[srcNode]] -= qWeight/srcDegSq;
			nodeDistribBuckets[dstNode][mystats.nodeToBucket[dstNode]] -= qWeight/dstDegSq;
			edgecnt++;
		}
		for(int i=0; i< upperBounditEdges.size(); i++)
		{
			maxRangeB = std::max(maxRangeB, upperBounditEdges[i]);
		}
		mystats.upperBoundRange = maxRangeB;
		mystats.allBoundsRanges = upperBounditEdges;
		t1VectComp = std::chrono::steady_clock::now();
		totVectComp = t1VectComp-t0VectComp;
		mystats.singleTimes.push_back(totVectComp.count());
		hoeffdingSamp = std::ceil(std::pow(maxRangeB, 2)/(2*std::pow(epsilon, 2)) * std::log(2*mystats.buckets.size()/eta));
		pseudoDimBound = std::floor(std::log2(boundPseudoDimension)) + 1;

		pdimSamp = std::ceil(std::pow(1.*maxRangeB, 2)/(2.*std::pow(epsilon, 2)) * (pseudoDimBound + std::log(1./eta)));
		double epsilonBennett = (epsilon*maxRangeB)/(m1-1);
		double logBennett = (1+epsilonBennett)*log(1+epsilonBennett) - epsilonBennett;

		long long int bennettBound = std::ceil(std::pow(1.*maxRangeB, 2)/(m1-1.) * 1./logBennett * log(2.*mystats.buckets.size()/eta));
		t1Mid = std::chrono::steady_clock::now();
		totMid = t1Mid-t0Mid;
		mystats.singleTimes.push_back(totMid.count());

		// Start sampling procedure
		t0Mid = std::chrono::steady_clock::now();
		
		double timeSpentCheking = 0;
		iterations = (std::ceil(1.3*maxRangeB/epsilon * log(4*mystats.buckets.size()/curr_eta))+1);
	}
	else if(mode == MODEFIXEDSS)
	{
		iterations = samplesI;
	}
    t1A = std::chrono::steady_clock::now();
    totA = t1A-t0A;
	mystats.timeub =  totA.count();

    t0A = std::chrono::steady_clock::now();
    std::vector<long long int> sampledEdgePositions;

    std::vector<double> repTris(n0, 0);
	std::vector<std::vector<double>> toVar(mystats.buckets.size(), std::vector<double>(sampledEdgePositions.size(), 0));
	int itEdge = 0;

    auto t0Var = std::chrono::steady_clock::now();
    auto t1Var = std::chrono::steady_clock::now();
    std::chrono::duration<double> totVar = t1Var - t0Var;
	
	std::vector<double> output(mystats.buckets.size(), 0);
	double timeVariance = 0;

	bool stoppingCondition = false;
	long long int totSamples = 0;
	long long int currIteration = 1;
	if(mode==MODEQVAL)
	{
		iterations = params.iterations;
	}
	std::vector<double> ubsSD(mystats.buckets.size(), 0);
	while(!stoppingCondition)
	{

		sampledEdgePositions.clear();
		iterations *= std::pow(1.4, currIteration-1); // Sample Schedule

		totSamples += iterations;
		for(int j=0; j<iterations; j++)
		{
			sampledEdgePositions.push_back(imp(gen64));
		}

	
		for(int k=0; k < toVar.size(); k++)
		{
			toVar[k].resize(totSamples);
		}

		for(auto& edgePos : sampledEdgePositions)
		{
			auto currEdge = edgeList[edgePos];
			auto v1 = currEdge.getSRC();
			auto v2 = currEdge.getDST();

			auto degv1 = nodeDegrees[v1].second;
			auto degv2 = nodeDegrees[v2].second;

			auto mindegnode = (degv1 >= degv2) ? v2 : v1;
			auto maxdegnode = (degv1 < degv2) ? v2 : v1;

			auto edgeProb = probImp[edgePos];

			long long int tris=0;
			for(auto& v3 : currAdjs[mindegnode])
			{
				if(currAdjs[v3].find(maxdegnode) !=currAdjs[v3].end())
				{
					repTris[v3] += CQ/(edgeProb);
					tris++;
					double functionValue = (CQ)/((edgeProb) * degreeSquared[v3]);
					functionValue /= mystats.buckets[mystats.nodeToBucket[v3]].size();
					toVar[mystats.nodeToBucket[v3]][itEdge] += functionValue;
				}
			}
			if(tris > 0)
			{
				repTris[v1] += (1.*tris*Q)/(edgeProb);
				repTris[v2] += (1.*tris*Q)/(edgeProb);
				double functionValue = (1.*tris*Q)/((edgeProb) * degreeSquared[v1]);
				functionValue /= mystats.buckets[mystats.nodeToBucket[v1]].size();
				toVar[mystats.nodeToBucket[v1]][itEdge] += functionValue;

				functionValue = (1.*tris*Q)/((edgeProb) * degreeSquared[v2]);
				functionValue /= mystats.buckets[mystats.nodeToBucket[v2]].size();
				toVar[mystats.nodeToBucket[v2]][itEdge] += functionValue;
			}
			itEdge++;
		}

		int itBuck = 0;

    	t0Var = std::chrono::steady_clock::now();
		bool stoppingmet = true;
		double currsup = 0;
		for(auto& sampB : toVar)
		{
			double mean = 0;
			double var = 0;
			if(mystats.buckets[itBuck].size() > 0)
			{
				if((mode == MODEFIXEDSS) || (mode == MODEQVAL))
				{
					for(auto samp : sampB)
					{
						mean += samp;
					}
					mean /= totSamples;
					output[itBuck] =  mean+smallBucketCountsNodes[itBuck];
				}
				else
				{
					//PrPl-EB structs:
					std::vector<double> mu_hat_t{0.}; // mu_hat_t[0]=0
					double sum_squared = std::pow(upperBounditEdges[itBuck], 2)/4;
					std::vector<double> sigma_t{sum_squared}; // mu_hat_t[0]=R_j^2/4 by Popovicious
					std::vector<double> lambda_vector;
					long long int timeT = 1;

					for(auto samp : sampB)
					{
						mean += samp;
						mu_hat_t.push_back(mean/timeT);
						sum_squared += std::pow(samp-mu_hat_t[timeT-1], 2);
						sigma_t.push_back(sum_squared/timeT);
						double lambda_val = std::min(1./(2.*upperBounditEdges[itBuck]), std::pow(2*std::log(2*mystats.buckets.size()/eta)/(totSamples*std::log(totSamples+1)*sigma_t[timeT-1]), 0.5));
						lambda_vector.push_back(lambda_val);
						timeT++;
					}
					double numer_tilde = 0;
					double denum_tilde = 0;
					double psi_term = 0;
					for(int sidx=0; sidx < sampB.size(); sidx++)
					{
						numer_tilde += sampB[sidx] * lambda_vector[sidx];
						denum_tilde += lambda_vector[sidx];
						if(sidx==0)
							continue;
						double log_lambda = (-std::log(1-(upperBounditEdges[itBuck]*lambda_vector[sidx-1])) - (upperBounditEdges[itBuck]*lambda_vector[sidx-1]))/4;
						psi_term += log_lambda*std::pow(sampB[sidx]-mu_hat_t[sidx-1], 2);
					}
					double c_var = std::pow(2/upperBounditEdges[itBuck], 2);
					double mean_tilde = numer_tilde/denum_tilde;
					double CI_lambda = (std::log(2*mystats.buckets.size()/eta) + (c_var*psi_term))/denum_tilde;
					
					if(nodesActiveInParts[itBuck] == 0)
					{
						mean_tilde = 0;
						CI_lambda = 0;
						output[itBuck] = smallBucketCountsNodes[itBuck];
					}
					else
					{
						output[itBuck] =  mean_tilde+smallBucketCountsNodes[itBuck];
					}

					mean /= totSamples;
					for(auto samp : sampB)
					{
						var += std::pow((mean-samp), 2);
					}
					var /= (totSamples-1.0);
					double logterm = std::log(3.0*mystats.buckets.size()/curr_eta);
					double varianceterm = std::pow((2.*var*logterm/(totSamples)), 0.5);
					double currBoundRange = upperBounditEdges[itBuck];
					double epshat = varianceterm + (currBoundRange *((3*logterm)/(totSamples-1.0)));


					double nbf = (totSamples/m1 < (m1/2.)) ? (1-(totSamples/m1)) : (1-(totSamples/m1))*(1+1/totSamples);
					double kappa = 7./3. + 3./std::pow(2.,0.5);
					double epshatBernfling = std::pow((2.*var*nbf*logterm/(totSamples)), 0.5) + (kappa *currBoundRange * logterm/(totSamples));//varianceterm + (maxRangeB *((3*logterm)/(totSamples-1.0)));
					epshat = std::min(epshat, CI_lambda);
					ubsSD[itBuck] = CI_lambda;



					currsup = std::max(currsup, CI_lambda);
					if((params.nonUni) && (epshat > epsilons[itBuck]))
					{
						std::cout << "Using eps bucket: " << itBuck << " " << epsilons[itBuck] << '\n';
						stoppingmet = false;
					}
					else if(epshat > epsilon)
						stoppingmet = false;
				}
			}
			itBuck++;
		}
		if(stoppingmet)
		{
			stoppingCondition = true;
			mystats.upperBoundSD = currsup;
			mystats.nonUniformUbs = ubsSD;
		}
		if((mode == MODEQVAL) || (mode == MODEFIXEDSS))
		{
			stoppingCondition = true;
			mystats.upperBoundSD = currsup;
			mystats.nonUniformUbs = ubsSD;
		}


    	t1Var = std::chrono::steady_clock::now();
    	totVar = t1Var - t0Var;
		timeVariance += totVar.count();
		currIteration++;
		curr_eta /= 2.0;
	}
    t1A = std::chrono::steady_clock::now();
    totA = t1A-t0A;
	mystats.timeadapt =  totA.count();
	mystats.singleTimes.push_back(timeVariance);
    t1Mid = std::chrono::steady_clock::now();
    totMid = t1Mid-t0Mid;
	mystats.singleTimes.push_back(totMid.count());

	mystats.ss = itEdge;
    t1 = std::chrono::steady_clock::now();
    tot = t1-t0;
    auto tFINISH = std::chrono::steady_clock::now();
    std::chrono::duration<double> totRT = tFINISH-tSTART;
	mystats.runtime.push_back(totRT.count());
    return output;
}

std::vector<double> CountingAlgs::TriadPar(Structures::Graph G, Structures::Statistics& mystats, Structures::InputParams& params)
{

	double epsilon = params.epsilon; 
	double eta = params.eta; 
	int iterations = params.iterations; 
	double qval = params.qval;
	int mode = params.mode;
	long long int samplesI = params.samples;
	long long int maxwork = params.maxwork;
	std::vector<double> normVertices = params.normMetricVertices;
	int totThreads = params.threads;
	std::cout << "Threads: " << totThreads << '\n';
		
	std::vector<double> epsilons(mystats.buckets.size(), epsilon);
	if(params.nonUni)
	{
		epsilons = params.nonUniEps;
		epsilon = *std::max_element(epsilons.begin(), epsilons.end());
	}
    
    auto MD = G.getMaxDegree();
    auto n0 = G.GetCurrN();
    auto m0 = G.GetM();
    std::vector<double> degreeSquared(n0, 0.0); // Denominator in local clustering coefficient cc(v), d_v * (d_v - 1)/2, only with d_v > 1
	degreeSquared = normVertices;
    auto tSTART = std::chrono::steady_clock::now();

    long long int totn = 0;
    auto t0 = std::chrono::steady_clock::now();
    auto t1 = std::chrono::steady_clock::now();
    std::chrono::duration<double> tot = t1-t0;

    auto t0Mid = std::chrono::steady_clock::now();
    std::vector<std::pair<NodeType, long long int>> originalDegrees = G.getNodeDegrees();

    std::vector<bool> smallBool(n0, false);
    const int MAX_WORK_SMALL = maxwork; // constant for removing small degree nodes, O(MAX_WORK_SMALL * |V|)

    auto t0A = std::chrono::steady_clock::now();
    auto [boundDeg ,smallCountsTriangles] = SmallTriangleCount(std::ref(G), MAX_WORK_SMALL, std::ref(smallBool), -1);
    auto t1A = std::chrono::steady_clock::now();
    std::chrono::duration<double> totA = t1A-t0A;
	mystats.timepeel = totA.count();

    auto n1 = G.GetCurrN();
    auto m1 = G.GetM();
	mystats.n1= n1;
	mystats.m1= m1;
    auto nodeDegrees = G.getNodeDegrees();
    std::vector<double> retCC(n0,0.);
	std::vector<double> smallBucketCountsNodes(mystats.buckets.size(),0);
	std::vector<long long int> nodesActiveInParts(mystats.buckets.size(),0);
    t0A = std::chrono::steady_clock::now();
    for(int i=0; i< smallCountsTriangles.size(); i++)
    {
        if(degreeSquared[i]>0)
		{
			if(!smallBool[i])
			{
				nodesActiveInParts[mystats.nodeToBucket[i]]++;
			}
            retCC[i] = smallCountsTriangles[i]/degreeSquared[i];
			double mysmallcc = (1.*smallCountsTriangles[i])/degreeSquared[i];
			smallBucketCountsNodes[mystats.nodeToBucket[i]] += mysmallcc/(1.*mystats.buckets[mystats.nodeToBucket[i]].size());
		}
        
    }
    auto t1Mid = std::chrono::steady_clock::now();
    std::chrono::duration<double> totMid = t1Mid-t0Mid;
	mystats.singleTimes.push_back(totMid.count());


    t0Mid = std::chrono::steady_clock::now();
    std::mt19937_64 gen64(mystats.seed);
    auto currAdjs = G.GetAdjsFast();
    auto edgeList = G.getEdgeList();
    double samplingProb{1.0/edgeList.size()};
    double Q{qval};
    int idx = 0;
	int bucketidx=0;
	int samplesOptVar=500;

    std::uniform_int_distribution<> uniformEdges(0, edgeList.size()-1);

    std::vector<double> edgeWeights(edgeList.size(), 1);
    std::discrete_distribution<> imp(edgeWeights.begin(), edgeWeights.end());
    std::vector<double> probImp = imp.probabilities();

	GRBEnv env = GRBEnv(true);
	env.set("LogFile", "grblog.log");
	env.start();
	// Create an empty model
	GRBModel model = GRBModel(env);
	// Model minimax problem
	GRBVar qModel = model.addVar(0.0, 0.5, 0.0, GRB_CONTINUOUS, "q");
	GRBVar tMax = model.addVar(0.0, m1/samplesOptVar, 0.0, GRB_CONTINUOUS, "t");

	model.setObjective(1*tMax, GRB_MINIMIZE);



	std::vector<long long int> edgesOptVar;
	for(int j=0; j<samplesOptVar; j++)
	{
		edgesOptVar.push_back(imp(gen64));
	}
	std::vector<std::vector<double>> atermEdge(mystats.buckets.size(), std::vector<double>(samplesOptVar,0.0));
	std::vector<std::vector<double>> btermEdge(mystats.buckets.size(), std::vector<double>(samplesOptVar,0.0));
	std::vector<std::vector<double>> ahatEdge(mystats.buckets.size(), std::vector<double>(samplesOptVar,0.0));
	std::vector<std::vector<double>> bhatEdge(mystats.buckets.size(), std::vector<double>(samplesOptVar,0.0));
	std::vector<double> aBar(mystats.buckets.size(), 0);
	std::vector<double> bBar(mystats.buckets.size(), 0);
	int idxEdgeVar =0;
	for(auto& edgePos : edgesOptVar)
	{
		auto currEdge = edgeList[edgePos];
		auto v1 = currEdge.getSRC();
		auto v2 = currEdge.getDST();

		auto degv1 = nodeDegrees[v1].second;
		auto degv2 = nodeDegrees[v2].second;

		auto mindegnode = (degv1 >= degv2) ? v2 : v1;
		auto maxdegnode = (degv1 < degv2) ? v2 : v1;

		auto edgeProb = probImp[edgePos];

		long long int tris=0;
		for(auto& v3 : currAdjs[mindegnode])
		{
			if(currAdjs[v3].find(maxdegnode) !=currAdjs[v3].end())
			{
				tris++;
				atermEdge[mystats.nodeToBucket[v3]][idxEdgeVar] += (1./degreeSquared[v3]);
				btermEdge[mystats.nodeToBucket[v3]][idxEdgeVar] -= (2./degreeSquared[v3]);
			}
		}
		if(tris > 0)
		{
			btermEdge[mystats.nodeToBucket[v1]][idxEdgeVar] -= (1.*tris/degreeSquared[v1]);
			btermEdge[mystats.nodeToBucket[v2]][idxEdgeVar] -= (1.*tris/degreeSquared[v2]);
		}
		idxEdgeVar++;
	}
	for(int i=0; i< atermEdge.size(); i++)
	{
		for(int j=0; j < atermEdge[i].size(); j++)
		{
			aBar[i] += atermEdge[i][j];
			bBar[i] += btermEdge[i][j];
		}
	}
	for(int i=0; i < atermEdge.size(); i++)
	{
		double meanA = aBar[i]/(1.*samplesOptVar);
		double meanB = bBar[i]/(1.*samplesOptVar);
		for(int j=0; j < atermEdge[i].size(); j++)
		{
			ahatEdge[i][j] = atermEdge[i][j] - meanA;
			bhatEdge[i][j] = btermEdge[i][j] - meanB;
		}
	}
	for(int i=0; i< atermEdge.size(); i++)
	{
		if((mystats.buckets[i].size()) > 0 && (nodesActiveInParts[i] > 0))
		{
			double A = 0;
			double B = 0;
			double C = 0;
			for(int j=0; j < atermEdge[i].size(); j++)
			{
				A += std::pow(ahatEdge[i][j], 2);
				B += std::pow(bhatEdge[i][j], 2);
				C += ahatEdge[i][j] * bhatEdge[i][j];
			}
			C *= 2;
			double norm = (1./samplingProb) * (1./mystats.buckets[i].size());
			norm = std::pow(norm, 2);
			norm /= (samplesOptVar-1.0);
			A = A*norm;
			B = B*norm;
			C = C*norm;
			if((std::fpclassify(A) == FP_ZERO) && (std::fpclassify(B) == FP_ZERO) && (std::fpclassify(C) == FP_ZERO))
				continue;
			model.addQConstr(1*tMax >= A + B*qModel*qModel + C*qModel);
		}
	}
	model.optimize();
	double minVar = model.get(GRB_DoubleAttr_ObjVal);
	double bestQMinVar = qModel.get(GRB_DoubleAttr_X);
	mystats.optimizedQ = bestQMinVar;
	params.qopt = bestQMinVar;
		
	Q = bestQMinVar;
	if(mode==MODEQVAL)
	{
		Q =qval;
	}
    double CQ{(1.0-(2.0*Q))};
    t1A = std::chrono::steady_clock::now();
    totA = t1A-t0A;
	mystats.timevar =  totA.count();
		
    t0A = std::chrono::steady_clock::now();
	std::vector<double> upperBounditEdges(mystats.buckets.size(), 0);
	double maxRangeB = 0;
	long long int hoeffdingSamp =0;
	double pseudoDimBound = 0;
	long long int pdimSamp =0;
	double curr_eta = eta;
	if(mode == MODEFIXEDSS)
	{
		iterations = samplesI;
	}
    t1A = std::chrono::steady_clock::now();
    totA = t1A-t0A;
	mystats.timeub =  totA.count();

    t0A = std::chrono::steady_clock::now();
    std::vector<long long int> sampledEdgePositions;

	std::vector<std::vector<double>> toVar(mystats.buckets.size(), std::vector<double>(sampledEdgePositions.size(), 0));
	int itEdge = 0;

    auto t0Var = std::chrono::steady_clock::now();
    auto t1Var = std::chrono::steady_clock::now();
    std::chrono::duration<double> totVar = t1Var - t0Var;
	
	std::vector<double> output(mystats.buckets.size(), 0);
	double timeVariance = 0;

	bool stoppingCondition = false;
	long long int totSamples = 0;
	long long int currIteration = 1;
	std::vector<double> ubsSD(mystats.buckets.size(), 0);
    std::vector<double> repTris(n0, 0);
	while(!stoppingCondition)
	{
		sampledEdgePositions.clear();
		iterations *= std::pow(1.4, currIteration-1); // Sample Schedule

		totSamples += iterations;
		for(int j=0; j<iterations; j++)
		{
			sampledEdgePositions.push_back(imp(gen64));
		}

	
		for(int k=0; k < toVar.size(); k++)
		{
			toVar[k].resize(totSamples);
		}

		#pragma omp parallel for num_threads(totThreads)
		for(int ep=0; ep < sampledEdgePositions.size(); ep++)
		{
			auto edgePos = sampledEdgePositions[ep];
			auto currEdge = edgeList[edgePos];
			auto v1 = currEdge.getSRC();
			auto v2 = currEdge.getDST();

			auto degv1 = nodeDegrees[v1].second;
			auto degv2 = nodeDegrees[v2].second;

			auto mindegnode = (degv1 >= degv2) ? v2 : v1;
			auto maxdegnode = (degv1 < degv2) ? v2 : v1;

			auto edgeProb = probImp[edgePos];

			long long int tris=0;
			std::unordered_map<Node_t, double> lazycnt;
			for(auto& v3 : currAdjs[mindegnode])
			{
				if(currAdjs[v3].find(maxdegnode) !=currAdjs[v3].end())
				{
					if(lazycnt.find(v3) != lazycnt.end())
						lazycnt[v3] += CQ/(edgeProb);
					else
						lazycnt[v3] = CQ/(edgeProb);
					tris++;
				}
			}
			if(tris > 0)
			{
				
				double functionValue = (1.*tris*Q)/((edgeProb) * degreeSquared[v1]);
				functionValue /= mystats.buckets[mystats.nodeToBucket[v1]].size();
				toVar[mystats.nodeToBucket[v1]][ep] += functionValue;

				functionValue = (1.*tris*Q)/((edgeProb) * degreeSquared[v2]);
				functionValue /= mystats.buckets[mystats.nodeToBucket[v2]].size();
				toVar[mystats.nodeToBucket[v2]][ep] += functionValue;

				for(auto& pair : lazycnt)
				{
					auto v3 = pair.first;
					double functionValue = pair.second / degreeSquared[v3];
					functionValue /= mystats.buckets[mystats.nodeToBucket[v3]].size();
					toVar[mystats.nodeToBucket[v3]][ep] += functionValue;
				}
			}
		}

		int itBuck = 0;

    	t0Var = std::chrono::steady_clock::now();
		bool stoppingmet = true;
		double currsup = 0;
		for(auto& sampB : toVar)
		{
			double mean = 0;
			double var = 0;
			if(mystats.buckets[itBuck].size() > 0)
			{
				if((mode == MODEFIXEDSS) || (mode == MODEQVAL))
				{
					for(auto samp : sampB)
					{
						mean += samp;
					}
					mean /= totSamples;
					output[itBuck] =  mean+smallBucketCountsNodes[itBuck];
				}
			}
			itBuck++;
		}
		if((mode == MODEQVAL) || (mode == MODEFIXEDSS))
		{
			stoppingCondition = true;
			mystats.upperBoundSD = currsup;
			mystats.nonUniformUbs = ubsSD;
		}


    	t1Var = std::chrono::steady_clock::now();
    	totVar = t1Var - t0Var;
		timeVariance += totVar.count();
		currIteration++;
		curr_eta /= 2.0;
	}
    t1A = std::chrono::steady_clock::now();
    totA = t1A-t0A;
	mystats.timeadapt =  totA.count();
	mystats.singleTimes.push_back(timeVariance);
    t1Mid = std::chrono::steady_clock::now();
    totMid = t1Mid-t0Mid;
	mystats.singleTimes.push_back(totMid.count());

	mystats.ss = itEdge;
    t1 = std::chrono::steady_clock::now();
    tot = t1-t0;
    auto tFINISH = std::chrono::steady_clock::now();
    std::chrono::duration<double> totRT = tFINISH-tSTART;
	mystats.runtime.push_back(totRT.count());
    return output;
}

std::vector<double> CountingAlgs::ApproximateCoeffsSota(Structures::Graph& G, Structures::Statistics& mystats, Structures::InputParams& params)
{
	auto tStart = std::chrono::steady_clock::now();
    auto n = G.GetN();
	auto tLimit = params.timeLimit;
    auto adjlistfast = G.GetAdjsFast();
	double epsilon = params.epsilon;
	double eta = params.eta;
    std::vector<double> approx(n, 0);
    std::mt19937_64 gen;
    gen.seed(mystats.seed);
	long long int nBound = n;
	std::vector<double> normVertices = params.normMetricVertices;
	int clustOrClos = params.clustOrClosure;
	auto buckets = mystats.nodeToBucket;
	std::vector<double> result(mystats.buckets.size(), 0);

    long long sampleSize = std::ceil(1./(2.*(std::pow(epsilon,2))) * std::log((2.*mystats.buckets.size())/eta)); // eta/2
	long long int totNodes = sampleSize * mystats.buckets.size();
    long long sampleSizeWedge = 1;
	long long int computedExactly = 0;
	long long int totWedgesExplored=0;
    for(long long int i=0; i< mystats.buckets.size(); i++)
    {
		if(mystats.buckets[i].size() < 1)
		{
			continue;
		}
		if(params.nonUni)
		{
			double boundEps = params.nonUniformErrors[i];
			sampleSize = std::ceil(1./(2.*(std::pow(boundEps, 2))) * std::log((2.*mystats.buckets.size())/eta)); // eta/2
    		sampleSizeWedge = 1;
			//std::cout << "bucket: " << i << " samp1: " << sampleSize << " samp2: " << sampleSizeWedge << " epsilon: " << boundEps << std::endl;
		}
		std::vector<Node_t> nodesToProc;
		long long int norm = sampleSize;//= mystats.buckets[i].size();
		
		nodesToProc.clear();
		std::uniform_int_distribution<> distr(0, mystats.buckets[i].size()-1);
		for(int j=0; j < sampleSize; j++)
		{
			auto sampledNode = mystats.buckets[i][distr(gen)];
			nodesToProc.push_back(sampledNode);	
		}

		double value = 0;
		for(auto& node : nodesToProc)
		{
			auto tNow = std::chrono::steady_clock::now();
			std::chrono::duration<double> totRTNow = tNow-tStart;
			if(totRTNow.count() > tLimit)
			{
				mystats.runtime.push_back(tLimit);
				return result;
			}

        	auto tot_degree = adjlistfast[node].size();
			int metricnode = 0;
			double tot_wedges = 0.0;
			if(tot_degree > 1)
			{
				tot_wedges = normVertices[node];
				if(clustOrClos == 1)
				{
					tot_wedges *= 2;
				}
			}
			else
			{
				continue;
			}
			if(clustOrClos == 0) // clusteringCoeff
			{
				std::vector<NodeType> neighbors(adjlistfast[node].size());
				std::copy(adjlistfast[node].begin(), adjlistfast[node].end(), neighbors.begin());
				std::uniform_int_distribution<> distr(0, neighbors.size()-1);
				auto first = distr(gen);
				auto second = distr(gen);
				if(first == second) 
				{
					while(first == second)
						second = distr(gen);
				}
				//Check edge between first second
				if(adjlistfast[neighbors[first]].find(neighbors[second]) != adjlistfast[neighbors[first]].end())
				{
					value +=1;
				}
			}
			else // closure coefficient
			{
				std::vector<NodeType> neighbors;
				std::unordered_map<long long int, std::vector<long long int>> samples;
				std::vector<long long int> cumulative;
				long long int cumulativeDegree = 0;
				for(auto& neigh : adjlistfast[node])
				{
					if(adjlistfast[neigh].size() > 1)
					{
						neighbors.push_back(neigh);
						cumulativeDegree += adjlistfast[neigh].size()-1;
						samples[cumulative.size()] = std::vector<long long int>();
						cumulative.push_back(cumulativeDegree);
					}
				}

				if(neighbors.size() > 1) 
				{
					std::uniform_int_distribution<> distr(0, cumulativeDegree-1);
					for(long long j=0; j < sampleSizeWedge; j++)
					{
						auto wedgeIdx = distr(gen);
						auto firstpos = std::upper_bound(cumulative.begin(), cumulative.end(), wedgeIdx) - cumulative.begin();
						samples[firstpos].push_back(wedgeIdx);
					}
					for(int idx=0; idx < cumulative.size(); idx++)
					{
						auto firstNeigh = neighbors[idx];
						std::vector<NodeType> neighbors2(adjlistfast[firstNeigh].size());
						std::copy(adjlistfast[firstNeigh].begin(), adjlistfast[firstNeigh].end(), neighbors2.begin());
						for(auto& wedgeIdx : samples[idx])
						{
							long long int second = (idx > 0) ? (wedgeIdx - cumulative[idx-1]) : wedgeIdx;
							if(node == neighbors2[second]) 
							{
								if(second == (neighbors2.size()-1))
									second--;
								else
									second++;
							}

							if(adjlistfast[neighbors2[second]].find(node) != adjlistfast[neighbors2[second]].end())
							{
								value +=1;
							}
						}
					}
				}
			}
		}
		result[i] = value/(1.*norm);
	}
    auto tFinish = std::chrono::steady_clock::now();
    std::chrono::duration<double> totRT = tFinish-tStart;
	mystats.runtime.push_back(totRT.count());
    return result;
}
