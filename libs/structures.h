#ifndef STRUCTURES_H
#define STRUCTURES_H

#include <vector>
#include <iostream>
#include <functional>
#include <string>
#include <sstream>
#include <fstream>
#include <cmath>
#include <random>
#include <set>
#include <unordered_map>
#include <unordered_set>
#include <nlohmann/json.hpp>
#include <highfive/H5Easy.hpp>
#include <chrono>
#include <queue>
#include "gurobi_c++.h"

using json = nlohmann::json;


#define MODEEX 0
#define MODEQVAL 1
#define MODETOHDF5 2
#define MODEEXCOUNTS 3
#define MODEPARTITION 4
#define MODEGETSTATS 5
#define MODESINGLERUN 6
#define MODETOTXT 10
#define MODEFIXEDSS 13
#define MODEPROCESSDBLPTEMPORAL 15
#define MODERUNWITHADAPTEPS 16
#define MODEPARALLEL 17

typedef long NodeType;
typedef long ColorType;
typedef long long int Node_t;
typedef long Edge_t;

namespace Structures {
    class Node;
    class Edge;
    class Graph;
	struct Statistics;
	struct InputParams;
	struct PartitionParams;
    void BuildGraph(std::string fname, Structures::Graph& G);
	void generateDBLPTemporal(std::string outfile, json jsonlabs);
    void csvToHdf5(std::string fname, std::string outname, bool isColored=false, bool remap=false, bool huge=false);
    void BuildGraph(H5Easy::File& datafile, Structures::Graph& G, bool oldread=false);
	void DumpExactCoeffs(H5Easy::File& datafile, Structures::Graph& G, std::vector<double>& triCounts);
	void DumpPartitions(H5Easy::File& datafile, Structures::Graph& G, Structures::PartitionParams& params);
	std::vector<std::vector<Node_t>> LoadPartition(H5Easy::File& datafile, std::string path);
	std::vector<double> LoadExactSolPartition(H5Easy::File& datafile, std::string path, std::vector<double>& ccVals);
	std::vector<long long int> LoadPartitionSize(H5Easy::File& datafile, std::string path);
}

struct Structures::PartitionParams{
	long long int seed=100;
	int numberDegreeBuckets=-1;
	int numberCoreBuckets=-1;
	bool labels=false;
	bool huge=false;
	std::string pathLabels;
};

struct Structures::Statistics{
	std::vector<double> runtime;
	std::vector<double> singleTimes;
	std::vector<double> sampleSize;
	std::vector<double> approxAvg;
	std::vector<double> approxSup;
	std::vector<long long int> wedgesExplored;
	std::vector<std::vector<Node_t>> buckets;
	std::vector<long long int> nodeToBucket;
	std::vector<double> variancePartitions;
	long long int nodesCountedExactly;
	long long int ss;
	double upperBoundRange;
	double upperBoundSD;
	std::vector<double> allBoundsRanges;
	std::vector<double> nonUniformUbs;
	long long int seed=100;
	double optimizedQ;
	long long int n1;
	long long int m1;
	double timepeel;
	double timevar;
	double timeub;
	double timeadapt;
};

struct Structures::InputParams{
	long long int seed=100;
	bool nonUni=false;
	double epsilon; 
	double timeLimit; 
	double eta; 
	int clustOrClosure=-1; // 0 for clustering and 1 for closure
	int iterations; 
	int steps; 
	double qval;
	long long int maxwork=10;
	long long int samples=0; // Samples to be taken if mode is  MODEFIXEDSSS
	int mode; // Used for testing a fixed sample size inside our adaptive algorigm (e.g., experiment with the value of q)
	std::vector<double> normMetricVertices; // This containes either (d(v) 2) or #of two paths starting at v, v\in V
	std::vector<double> debugExactVals;
	std::vector<double> nonUniformErrors;
	double qopt;
	int threads;
	std::vector<double> nonUniEps;
};

// Structures implementation
class Structures::Node{
    private: 
        NodeType id;
        ColorType c;
    public:
        Node() : id{-1}, c{-1}{}
        Node(NodeType ID, ColorType col) : id{ID}, c{col} {}
        NodeType getID() const {return id;}
        ColorType getColor() const {return c;}

        friend bool operator<(const Node& n1, const Node& n2)
        {
            if(n1.getID() == n2.getID())
                return false;
            else
                return n1.getID() < n2.getID();
        }
};

class Structures::Edge{
    private: 
        NodeType src;
        NodeType dst;
    public:
        Edge() : src{-1}, dst{-1}{}
        Edge(NodeType src, NodeType dst) : src{src}, dst{dst} {}
        NodeType getSRC() const {return src;}
        NodeType getDST() const {return dst;}

};

namespace CountingAlgs {
    std::vector<std::pair<long long int, long long int>> GetNodeDegrees(Structures::Graph G);
}

namespace Utilities {
    // Function to perform binomial sampling Bin(trials, probability) in expected time O(trials * probability)
    std::vector<long long int> sublinearBinomialSample(long long int trials, double probability, std::mt19937_64& generator);
	std::vector<std::string> splitString(std::string s, std::string delimiter);
	void dumpTxt(H5Easy::File& datafile, std::string outfile, bool oldread);
}

class Structures::Graph
{
    private:
        Node_t N;
        Edge_t M;
        long long int OriginalMaxDegree; // MaxDegree in the original graph before pruning
        std::vector<std::vector<NodeType>> AdjList; // for each node i=0,..., n-1 keep N_i (set of nodes adjacent to i) QUERYTIME: is j in N_i -> O(|N_i|) inefficient
        std::vector<std::unordered_set<NodeType>> AdjListFast;// for each node i=0,..., n-1 keep N_i (set of nodes adjacent to i) QUERYTIME: is j in N_i -> O(log |N_i|) a little bit better
        std::vector<std::unordered_set<NodeType>> origAdjs;// for each node i=0,..., n-1 keep N_i (set of nodes adjacent to i) QUERYTIME: is j in N_i -> O(log |N_i|) a little bit better
        //std::vector<std::unordered_map<NodeType>> AdjListFast;// for each node i=0,..., n-1 keep N_i (set of nodes adjacent to i) QUERYTIME: is j in N_i -> O(log |N_i|) a little bit better
        std::vector<Structures::Node> NodeList; // i=0,...,n-1 pos gives a possibly attributed node with id i
        std::vector<double> DegreeDistrib; // position i=0,...,d_max keeps |D_i| in the current subgraph D_i = {v \in V : d_v = i} (set of nodes with degree i)
        std::vector<bool> NodeAlive; // position i=0,...,n-1 contains a boolean value, true if node i was not removed, false otherwise
        Node_t NumAlive; // |{i=0,...,n-1}| such that "NodeAlive[i] == true"
        std::vector<std::pair<NodeType, long long int>> uncoloredNodeDegrees; //position i keeps a pair <node_id, degree_id> where degree_id is the degree in the current subgraph
        
        std::unordered_set<Node_t> SetAlive;

    public:
        Graph(Node_t N, Edge_t M) : N{N}, M{M} {
        }
        void setNodes(Node_t n)
        {
            N = n;
            AdjListFast.resize(n);
            NodeList.resize(n);
            NodeAlive.resize(n, true);
            NumAlive = n;
            DegreeDistrib.resize(n-1, 0);
            uncoloredNodeDegrees.resize(n);
            for(NodeType i=0; i < n; i++)
            {
                uncoloredNodeDegrees[i].first = i;
                SetAlive.insert(i);
                uncoloredNodeDegrees[i].second = 0;
            }
            return;
        }

        const Node_t GetN()
        {
            return N;
        }
        const Edge_t GetM()
        {
            return M;
        }
        void setEdges(Edge_t m)
        {
            M = m;
        }
        const Node_t GetCurrN()
        {
            return NumAlive;
        }
        const std::vector<bool>& getAlive()
        {
            return std::cref(NodeAlive); }
        void removeNode(Node_t id)
        {
            if(NumAlive < 1 or (!NodeAlive[id]))
            {
                std::cerr << "ERROR when removing: " << id << std::endl;
                return;
            }
            auto myDegree = AdjListFast[id].size(); // current degree before removing anything
            for(auto& neighbor : AdjListFast[id])
            {
                if(NodeAlive[neighbor])
                {
                    AdjListFast[neighbor].erase(id);
                    // Decreasing the degree of the 'neighbor' by 1
                    DegreeDistrib[uncoloredNodeDegrees[neighbor].second]--;
                    uncoloredNodeDegrees[neighbor].second--;
                    DegreeDistrib[uncoloredNodeDegrees[neighbor].second]++;
                }
            }
            DegreeDistrib[myDegree] -= 1.0;
            AdjListFast[id].clear();
            uncoloredNodeDegrees[id].second = 0;
            NodeAlive[id] = false;
            SetAlive.erase(id); // Remove it from the set of alive nodes
            NumAlive--;
            M-= myDegree;
            return;
        }

        long long int getMinDegree() // return the minimum degree in the current graph O(d_max), less efficient when gap in min degree is far from constant. otherwise constant
        {
            long long int minDeg = 1;
            for(int i=2; i < DegreeDistrib.size(); i++)
            {
                minDeg++;
                if(DegreeDistrib[i]>0.5) // There exists a node with the current degree
                {
                    break;
                }
            }
            return minDeg;
        }
        std::vector<Structures::Edge> getEdgeList()
        {
            std::vector<Structures::Edge> edgeList;
            for(auto& node : SetAlive)
            {
                for(auto& neighbor : AdjListFast[node])
                {
                    if(!NodeAlive[neighbor])
                        std::cout << "ERROR in ALIVES" << std::endl;
                    if(node < neighbor) // Do not enumerate each edge twice
                    {
                        Structures::Edge e = Structures::Edge(node, neighbor);
                        edgeList.push_back(e);
                    }
                }
            }
            return edgeList;
        }

        long long int getMaxDegree() // return the maximum degree in the current graph O(d_max), less efficient when max degree shrinks much, otherwise almost constant
        {
            long long int maxDeg = DegreeDistrib.size();
            for(auto it = DegreeDistrib.rbegin(); it!= DegreeDistrib.rend(); ++it)
            {
                maxDeg--;
                if((*it)>0.5) // There exists a node with the current degree
                {
                    break;
                }
            }
            return maxDeg;
        }
        std::vector<std::vector<NodeType>> MaterializeDegreeBuckets()
        {
            //auto maxDeg = getMaxDegree();
            auto maxDeg = OriginalMaxDegree;
            std::vector<std::vector<NodeType>> buckets(maxDeg+1);
            for(int node = 0; node < NodeAlive.size(); node++) // O(n)
            {
                if(NodeAlive[node])
                    buckets[uncoloredNodeDegrees[node].second].push_back(node);
            }
            //std::cout << "Max degree is: " << maxDeg << std::endl;
            return buckets;
        }
        std::vector<std::unordered_set<NodeType>> MaterializeDegreeBucketsSets()
        {
            //auto maxDeg = getMaxDegree();
            auto maxDeg = OriginalMaxDegree;
            std::vector<std::unordered_set<NodeType>> buckets(maxDeg+1);
            for(int node = 0; node < NodeAlive.size(); node++) // O(n)
            {
                if(NodeAlive[node])
                    buckets[uncoloredNodeDegrees[node].second].insert(node);
            }
            //std::cout << "Max degree is: " << maxDeg << std::endl;
            return buckets;
        }
        
        void addEdge(const NodeType& src, const NodeType& dst)
        {
			if(src == dst)
				return;
            uncoloredNodeDegrees[dst].second++;
            uncoloredNodeDegrees[src].second++;
            AdjListFast[src].insert(dst);
            AdjListFast[dst].insert(src);
        }
        const std::vector<std::pair<NodeType, long long int>>& getSortedUncoloredDegrees()
        {
            std::sort(uncoloredNodeDegrees.begin(), uncoloredNodeDegrees.end(), [](auto& left, auto& right){ return left.second > right.second; });
            return std::cref(uncoloredNodeDegrees);

        }
        const std::vector<std::pair<NodeType, long long int>>& getNodeDegrees()
        {
            return std::cref(uncoloredNodeDegrees);
        }
        void fixDegreeDistrib()
        {
            auto max_degree = std::max_element(uncoloredNodeDegrees.begin(), uncoloredNodeDegrees.end(), [](auto&left, auto&right) {return left.second < right.second;})->second;
            DegreeDistrib.resize(static_cast<unsigned long long int>(max_degree+1), 0); // a degree of a node cannot exceed max
            for(auto& nodeDeg : uncoloredNodeDegrees)
            {
                DegreeDistrib[nodeDeg.second]++;
            }
            OriginalMaxDegree = max_degree;
            return;
        }
        const std::vector<double>& getDegreeDistrib() // O(1)
        {
            return std::cref(DegreeDistrib);
        }
        void addNode(const Structures::Node v)
        {
            NodeList[v.getID()] = v;
        }
        const std::vector<Structures::Node>& GetNodes()
        {
            return std::cref(NodeList);
        }
        const std::vector<std::vector<NodeType>>& GetAdjs()
        {
            return std::cref(AdjList);
        }
        const std::unordered_set<NodeType>& GetNodeNeighborhood(NodeType id)
        {
            if(!isNodeAlive(id))
                std::cerr << "ERROR NODE SHOULD BE DEAD: " << std::endl;
            return std::cref(AdjListFast[id]);
        }
        std::vector<std::unordered_set<NodeType>>& GetAdjsFast()
        {
            return std::ref(AdjListFast);
        }
        std::vector<std::unordered_set<NodeType>>& GetOrigAdjs()
        {
            return std::ref(origAdjs);
        }
        bool isNodeAlive(NodeType id)
        {
            return NodeAlive[id];
        }
		bool hasEdge(Node_t src, Node_t dst)
		{
			bool answer = (AdjListFast[src].find(dst) != AdjListFast[src].end());
			return answer;
		}
		void swap(Node_t src1, Node_t dst1, Node_t src2, Node_t dst2)
		{
			// Remove first edge
			AdjListFast[src1].erase(dst1);
			AdjListFast[dst1].erase(src1);
			// Remove second edge
			AdjListFast[src2].erase(dst2);
			AdjListFast[dst2].erase(src2);

			// Add new swapped edges
			AdjListFast[src1].insert(dst2);
			AdjListFast[dst2].insert(src1);

			AdjListFast[src2].insert(dst1);
			AdjListFast[dst1].insert(src2);
			return;
		}
};

namespace CountingAlgs {
    std::vector<std::pair<long long int, long long int>> GetNodeDegrees(Structures::Graph G);
    std::vector<double> clustCoeffsWithChiba(Structures::Graph G, Structures::Statistics& mystats);

    std::tuple<int, std::vector<double>> SmallTriangleCount(Structures::Graph& G, int threshold, std::vector<bool>& isSmall, int ub=-1);
    std::vector<double> SimpleAdapt(Structures::Graph G, Structures::Statistics& mystats, Structures::InputParams& params);
	std::vector<double> TriadPar(Structures::Graph G, Structures::Statistics& mystats, Structures::InputParams& params);

    int GetCoreDecomposition(Structures::Graph G, std::vector<std::vector<NodeType>>& cores, long long int N=-1);

	std::vector<double> ApproximateCoeffsSota(Structures::Graph& G, Structures::Statistics& mystats, Structures::InputParams& params);


}

#endif //STRUCTURES_H
