#include "structures.h"
#include <iostream>
#include <ostream>
#include <functional>
#include <fstream>
#include <iomanip>
#include <nlohmann/json.hpp>
#include <highfive/H5Easy.hpp>

using json = nlohmann::json;

int main(int argc, char *argv[])
{
    std::string jsonfile(argv[1]);
    std::ifstream f(jsonfile);
    json exe_spec = json::parse(f);
    f.close();

    int mode = exe_spec["mode"].get<int>();
	if(mode == MODETOHDF5) // Take file from csv format and make it hdf5
	{
		std::string fedges = exe_spec["dataset"].get<std::string>();
		std::string outfile = exe_spec["outFile"].get<std::string>();
		bool remap = exe_spec["remap"].get<bool>();
		bool huge = exe_spec["huge"].get<bool>();

		H5Easy::DumpOptions options{H5Easy::Compression()};
		Structures::csvToHdf5(fedges, outfile, false, remap, huge);
		return 0;
	}
	else if(mode == MODEPROCESSDBLPTEMPORAL) // Generate DBLP data from csv and dicts
	{
		std::string jsonlabs = exe_spec["fileData"].get<std::string>();
		std::string outfile = exe_spec["outFile"].get<std::string>();

		std::ifstream flab(jsonlabs);
		json labelAndIdMap = json::parse(flab);
		flab.close();

		std::cout << "Running genBDLP\n";
		H5Easy::DumpOptions options{H5Easy::Compression()};
		Structures::generateDBLPTemporal(outfile, labelAndIdMap);

		return 0;
	}
	else if(mode == MODEEXCOUNTS) // Compute and store exact triangle counts, and clustering and closure coefficient values for each node v\in V
	{
		std::string fedges = exe_spec["dataset"].get<std::string>();
		H5Easy::File hdf5Data(fedges, H5Easy::File::ReadWrite);
   	 	Structures::Graph mygraph{0,0}; //Create emptyGraph
		bool oldread = exe_spec["oldread"].get<bool>();
		bool huge = exe_spec["huge"].get<bool>();
		Structures::BuildGraph(std::ref(hdf5Data), std::ref(mygraph), oldread);
		std::cout << "loaded graph has: " << mygraph.GetN() << " and " << mygraph.GetM() << '\n';
		Structures::Statistics chibaStats;
    	std::vector<double> localTriangleCounts;
		int totIters{5};
		if(huge) totIters = 1;
		for(int i=0; i < totIters; i++)
		{
			localTriangleCounts.clear();
			localTriangleCounts = CountingAlgs::clustCoeffsWithChiba(mygraph, std::ref(chibaStats));
		}
		H5Easy::DumpOptions options{H5Easy::Compression(), H5Easy::DumpMode::Overwrite};
		std::string pathWrite = "/exact-counts";
		std::cout << "dumping exactCounts\n";
		H5Easy::dump(hdf5Data, pathWrite + "/solution", localTriangleCounts, options);
		std::cout << "dumping exactRuntime\n";
		H5Easy::dump(hdf5Data, pathWrite + "/runtime", chibaStats.runtime, options);

		// Compute clustCoeff(v) and closureCoeff(v), v \in V and dump them
		Structures::DumpExactCoeffs(std::ref(hdf5Data), std::ref(mygraph), std::ref(localTriangleCounts));
		return 0;	
	}
	else if(mode == MODEPARTITION) // Compute and store exact triangle counts, and clustering and closure coefficient values for each node v\in V
	{
		std::string fedges = exe_spec["dataset"].get<std::string>();
		int partsBucks = exe_spec["bucketsDegs"].get<int>();
		int partsCores = exe_spec["partsCores"].get<int>();
		bool labs = exe_spec["hasLabels"].get<bool>();
		std::string labsPath = exe_spec["labsPath"].get<std::string>();
		bool huge = exe_spec["huge"].get<bool>();
		bool oldread = exe_spec["oldread"].get<bool>();
		Structures::PartitionParams params = {.seed=150, .numberDegreeBuckets=partsBucks, .numberCoreBuckets=partsCores, 
			.labels=labs, .huge=huge, .pathLabels=labsPath};
		H5Easy::File hdf5Data(fedges, H5Easy::File::ReadWrite);

   	 	Structures::Graph mygraph{0,0}; //Create emptyGraph
		Structures::BuildGraph(std::ref(hdf5Data), std::ref(mygraph), oldread);
		std::cout << "loaded graph has: " << mygraph.GetN() << " and " << mygraph.GetM() << '\n';

		DumpPartitions(std::ref(hdf5Data), std::ref(mygraph), std::ref(params));
		return 0;	
	}
	else if(mode == MODETOTXT) // Compute and store exact triangle counts, and clustering and closure coefficient values for each node v\in V
	{
		std::string fedges = exe_spec["dataset"].get<std::string>();
		std::string outfile = exe_spec["outfile"].get<std::string>();
		bool oldread = exe_spec["oldread"].get<bool>();
		H5Easy::File hdf5Data(fedges, H5Easy::File::ReadOnly);

		Utilities::dumpTxt(std::ref(hdf5Data), outfile, oldread);

		return 0;	
	}
	else if(mode == MODEGETSTATS) // Compute and store exact triangle counts, and clustering and closure coefficient values for each node v\in V
	{
		std::string fedges = exe_spec["dataset"].get<std::string>();
		H5Easy::File hdf5Data(fedges, H5Easy::File::ReadOnly);

    	std::string outJson = exe_spec["outJson"].get<std::string>();
		Structures::Graph mygraph{0,0}; //Create emptyGraph
		bool oldread = exe_spec["oldread"].get<bool>();
		Structures::BuildGraph(std::ref(hdf5Data), std::ref(mygraph), oldread);
		std::cout << "loaded graph has: " << mygraph.GetN() << " and " << mygraph.GetM() << '\n';
		auto MD = mygraph.getMaxDegree();
    	auto N = mygraph.GetN();
    	auto M = mygraph.GetM();


		std::vector<std::string> metrics{"/clusteringCoeffs","/closureCoeffs"};
		std::vector<std::string> pathsPartitions = exe_spec["paths"].get<std::vector<std::string>>();
		std::unordered_map<std::string, std::vector<json>> confM;
		std::unordered_map<std::string, double> avges;
		for(auto& metric : metrics) // Either clustcoeff or closure coeff
		{
			std::vector<json> confs;
			std::vector<double> exactValues = H5Easy::load<std::vector<double>>(hdf5Data, metric + "/values");// Load exact metric values
			double avgMetric = 0;
			for(auto & val : exactValues)
			{
				avgMetric += val;
			}
			avgMetric /= exactValues.size();
			avges[metric] = avgMetric;

			for(auto& path : pathsPartitions)
			{
				auto exactSol = Structures::LoadExactSolPartition(std::ref(hdf5Data), path, std::ref(exactValues));
				auto partSizes =  Structures::LoadPartitionSize(std::ref(hdf5Data), path);
				json conf = {
				{"Dataset", fedges},
				{"Partition", path},
				{"PartSizes", partSizes},
				{"Sol-CN", exactSol}
				};
				confs.push_back(conf);
			}
			confM[metric] = confs;
		}
		json outA = {
			{"Dataset", fedges},
			{"N-nodes", N},
			{"M-edges", M},
			{"Max-Deg", MD},
			{"Averages", avges},
			{"Buckets", confM}
			};
		std::ofstream outj1;
		outj1.open(outJson);
		outj1 << outA << '\n';
		outj1.close();
		return 0;	
	}
	else if(mode == MODEQVAL) // Compute and store exact triangle counts, and clustering and closure coefficient values for each node v\in V
	{
		std::string fedges = exe_spec["dataset"].get<std::string>();
    	std::string outJson = exe_spec["outJson"].get<std::string>();
		H5Easy::File hdf5Data(fedges, H5Easy::File::ReadOnly);

   	 	Structures::Graph mygraph{0,0}; //Create emptyGraph
		bool oldread = exe_spec["oldread"].get<bool>();
		Structures::BuildGraph(std::ref(hdf5Data), std::ref(mygraph), oldread);
		std::cout << "loaded graph has: " << mygraph.GetN() << " and " << mygraph.GetM() << '\n';
    	auto N = mygraph.GetN();
    	auto M = mygraph.GetM();

		Structures::InputParams paramsAppx;
		paramsAppx.epsilon = 0.05;
		paramsAppx.maxwork = 2;
		paramsAppx.eta = 0.1;
		paramsAppx.iterations = 100;
		paramsAppx.steps = 100;
		paramsAppx.qval = 0.01;
		paramsAppx.mode = MODEQVAL;

		int thMain = 1;//exe_spec["mainThreads"].get<int>();
		std::vector<double> itersPerc = exe_spec["itersPercVect"].get<std::vector<double>>();
    	std::vector<double> qvect = exe_spec["qvalVect"].get<std::vector<double>>();
		std::vector<std::string> pathsPartitions = exe_spec["paths"].get<std::vector<std::string>>();

		long long int Cseeds[10] {12,320,456,32345,2321,858,9048934,38475984,83128391,42};
		std::vector<std::string> metrics{"/clusteringCoeffs","/closureCoeffs"};
		std::unordered_map<std::string, std::vector<json>> exps;
		for(auto& metric: metrics) // Either clustcoeff or closure coeff
		{
			std::vector<double> exactValues = H5Easy::load<std::vector<double>>(hdf5Data, metric + "/values");// Load exact metric values
			paramsAppx.normMetricVertices = H5Easy::load<std::vector<double>>(hdf5Data, metric + "/normNodes");// Load denumerator for each vertex
			
			std::vector<json> confs;
			for(auto& path : pathsPartitions)
			{
				Structures::Statistics edgeSampleStats;
				auto buckets = Structures::LoadPartition(std::ref(hdf5Data), path);
				auto exactSol = Structures::LoadExactSolPartition(std::ref(hdf5Data), path, std::ref(exactValues));
				edgeSampleStats.buckets = buckets;
				edgeSampleStats.nodeToBucket = H5Easy::load<std::vector<long long int>>(hdf5Data, path);
				std::unordered_map<std::string, std::unordered_map<std::string, std::vector<std::vector<double>> >> resq;
				std::unordered_map<std::string, std::vector<std::vector<double>> > optQ;
				#pragma omp parallel for num_threads(thMain)
				for(auto& q : qvect)
				{
					auto myparams = paramsAppx;
					myparams.qval = q;
					std::unordered_map<std::string, std::vector<std::vector<double>> > fixedq;
					std::vector<std::vector<double>> qalgPerc;
					for(auto& perc : itersPerc)
					{
						std::string percString = std::to_string(perc);
						long long int iterations = std::ceil(M*perc);
						myparams.iterations = iterations;
						std::vector<double> qalg;
						for(int it=0; it<10; it++)
						{
							auto mygraphCopy = mygraph;
							edgeSampleStats.seed = Cseeds[it];
							std::vector<double> adaptiveCoeffs = CountingAlgs::SimpleAdapt(mygraphCopy, edgeSampleStats, std::ref(myparams));
							qalg.push_back(myparams.qopt);
							if(it == 0)
							{
								std::vector<std::vector<double>> currRes{adaptiveCoeffs};
								fixedq[percString] = currRes;
							}
							else
								fixedq[percString].push_back(adaptiveCoeffs);
						}
						qalgPerc.push_back(qalg);
					}
					std::string percQ = std::to_string(q);
					#pragma omp critical
					resq[percQ] = fixedq;
					#pragma omp critical
					optQ[percQ] = qalgPerc;
				}
				json conf = {
				{"Dataset", fedges},
				{"Partition", path},
				{"Sol-CN", exactSol},
				{"VaryingQ", resq},
				{"OptQ", optQ},
				{"RT-ES", edgeSampleStats.runtime},
				{"samplePerc", itersPerc},
				{"QValues", qvect}
				};
				confs.push_back(conf);
			}
			exps[metric] = confs;
		}
		json outA = {
			{"Dataset", fedges},
			{"N-nodes", N},
			{"M-edges", M},
			{"Results", exps}
			};
		std::ofstream outj1;
		outj1.open(outJson);
		outj1 << outA << '\n';
		outj1.close();
		
		return 0;	
	}
	else if(mode == MODEFIXEDSS) // 
	{
		std::string fedges = exe_spec["dataset"].get<std::string>();
    	std::string outJson = exe_spec["outJson"].get<std::string>();
		H5Easy::File hdf5Data(fedges, H5Easy::File::ReadOnly);

   	 	Structures::Graph mygraph{0,0}; //Create emptyGraph
		bool oldread = exe_spec["oldread"].get<bool>();
		long long int maxwork=exe_spec["worksmall"].get<long long int>();
		long long int samples=exe_spec["samples"].get<long long int>();
		Structures::BuildGraph(std::ref(hdf5Data), std::ref(mygraph), oldread);
		std::cout << "loaded graph has: " << mygraph.GetN() << " and " << mygraph.GetM() << '\n';
    	auto N = mygraph.GetN();
    	auto M = mygraph.GetM();

		Structures::InputParams paramsAppx;
		paramsAppx.epsilon = 0.05;
		paramsAppx.maxwork = maxwork;
		paramsAppx.eta = 0.1;
		paramsAppx.iterations = samples;
		paramsAppx.samples = samples;
		paramsAppx.steps = 100;
		paramsAppx.qval = 0.01;
		paramsAppx.mode = MODEFIXEDSS;

		std::vector<std::string> pathsPartitions = exe_spec["paths"].get<std::vector<std::string>>();

		long long int Cseeds[10] {12,320,456,32345,2321,858,9048934,38475984,83128391,42};
		std::vector<std::string> metrics{"/clusteringCoeffs","/closureCoeffs"};
		std::unordered_map<std::string, std::vector<json>> exps;
		for(auto& metric: metrics) // Either clustcoeff or closure coeff
		{
			std::vector<double> exactValues = H5Easy::load<std::vector<double>>(hdf5Data, metric + "/values");// Load exact metric values
			paramsAppx.normMetricVertices = H5Easy::load<std::vector<double>>(hdf5Data, metric + "/normNodes");// Load denumerator for each vertex
			
			std::vector<json> confs;
			for(auto& path : pathsPartitions)
			{
				Structures::Statistics edgeSampleStats;
				auto buckets = Structures::LoadPartition(std::ref(hdf5Data), path);
				auto exactSol = Structures::LoadExactSolPartition(std::ref(hdf5Data), path, std::ref(exactValues));
				edgeSampleStats.buckets = buckets;
				edgeSampleStats.nodeToBucket = H5Easy::load<std::vector<long long int>>(hdf5Data, path);
				std::vector<std::vector<double>> approxes;
				std::vector<double> qvect;

				for(int i=0; i < 5; i++)
				{
					auto myparams = paramsAppx;
					edgeSampleStats.seed = Cseeds[i];
					std::unordered_map<std::string, std::vector<std::vector<double>> > fixedq;
					auto mygraphCopy = mygraph;
					std::vector<double> adaptiveCoeffs = CountingAlgs::SimpleAdapt(mygraphCopy, edgeSampleStats, std::ref(myparams));
					approxes.push_back(adaptiveCoeffs);
					qvect.push_back(edgeSampleStats.optimizedQ);
				}
				std::vector<double> runtimes = edgeSampleStats.runtime;
				json conf = {
				{"Partition", path},
				{"Sol-CN", exactSol},
				{"Approx", approxes},
				{"Samples", samples},
				{"RT-ours", runtimes},
				{"QValues", qvect}
				};
				confs.push_back(conf);
			}
			exps[metric] = confs;
		}
		json outA = {
			{"Dataset", fedges},
			{"N-nodes", N},
			{"M-edges", M},
			{"Results", exps}
			};
		std::ofstream outj1;
		outj1.open(outJson);
		outj1 << outA << '\n';
		outj1.close();
		
		return 0;	
	}
	else if(mode == MODERUNWITHADAPTEPS) // Compute and store exact triangle counts, and clustering and closure coefficient values for each node v\in V
	{
		std::string fedges = exe_spec["dataset"].get<std::string>();
    	std::string outJson = exe_spec["outJson"].get<std::string>();
		H5Easy::File hdf5Data(fedges, H5Easy::File::ReadOnly);

   	 	Structures::Graph mygraph{0,0}; //Create emptyGraph
		bool oldread = exe_spec["oldread"].get<bool>();
		Structures::BuildGraph(std::ref(hdf5Data), std::ref(mygraph), oldread);
		std::cout << "loaded graph has: " << mygraph.GetN() << " and " << mygraph.GetM() << '\n';
    	auto N = mygraph.GetN();
    	auto M = mygraph.GetM();
    	std::vector<double> epsVect = exe_spec["epsilon"].get<std::vector<double>>();
		long long int maxwork=exe_spec["worksmall"].get<long long int>();
		double timeLimit=exe_spec["timeLimit"].get<double>();

		Structures::InputParams paramsAppx;
		paramsAppx.epsilon = epsVect[0];
		paramsAppx.eta = 0.01;
		paramsAppx.maxwork = maxwork;
		paramsAppx.iterations = 100;
		paramsAppx.steps = 100;
		paramsAppx.qval = 0.01;
		paramsAppx.mode = MODESINGLERUN;

		long long int Cseeds[10] {12,320,456,32345,2321,858,9048934,38475984,83128391,42};
		std::vector<std::string> pathsPartitions = exe_spec["paths"].get<std::vector<std::string>>();

		std::vector<std::string> metrics{"/clusteringCoeffs", "/closureCoeffs"};
		std::unordered_map<std::string, std::vector<json>> exps;
		int meth = 0;
		for(auto& metric: metrics) // Either clustcoeff or closure coeff
		{
			std::cout << " metric: " << metric << '\n';
			paramsAppx.normMetricVertices = H5Easy::load<std::vector<double>>(hdf5Data, metric + "/normNodes");// Load denumerator for each vertex

			std::vector<json> confs;
			if(metric == "/closureCoeffs")
				meth = 1;
			else
				meth = 0;
			for(auto& path : pathsPartitions)
			{
				auto wsParams = paramsAppx;
				wsParams.timeLimit = timeLimit;
				Structures::Statistics edgeSampleStats;
				std::vector<double> exactValues = H5Easy::load<std::vector<double>>(hdf5Data, metric + "/values");// Load exact metric values

				auto exactSol = Structures::LoadExactSolPartition(std::ref(hdf5Data), path, std::ref(exactValues));
				auto buckets = Structures::LoadPartition(std::ref(hdf5Data), path);

				edgeSampleStats.buckets = buckets;
				edgeSampleStats.nodeToBucket = H5Easy::load<std::vector<long long int>>(hdf5Data, path);

				std::vector<std::vector<double>> approxes;
				std::vector<long long int> samplesS;
				std::vector<std::vector<double>> nonUniBounds;
				std::vector<std::vector<double>> appxBL;
				std::vector<double> qvect;
				Structures::Statistics wsStats;
				std::vector<std::string> legendTimes {"timepeel", "timevar", "timeub", "timeadapt"};
				std::vector<std::vector<double>> breakTimes;
				for(int i=0; i < 5; i++)
				{
					edgeSampleStats.seed = Cseeds[i];
					std::unordered_map<std::string, std::vector<std::vector<double>> > fixedq;
					auto mygraphCopy = mygraph;

					auto myparams = paramsAppx;
					wsParams.clustOrClosure = meth;
					wsStats.buckets = buckets;
					wsParams.seed =Cseeds[i];
					wsStats.nodeToBucket = H5Easy::load<std::vector<long long int>>(hdf5Data, path);
					wsParams.debugExactVals = exactValues;

					std::vector<double> adaptiveCoeffs = CountingAlgs::SimpleAdapt(mygraphCopy, edgeSampleStats, std::ref(myparams));
					samplesS.push_back(edgeSampleStats.ss);
					breakTimes.push_back(std::vector<double>{edgeSampleStats.timepeel, edgeSampleStats.timevar, edgeSampleStats.timeub, edgeSampleStats.timeadapt});
					wsParams.timeLimit = edgeSampleStats.runtime.back()*10;

					approxes.push_back(adaptiveCoeffs);
					qvect.push_back(edgeSampleStats.optimizedQ);
					nonUniBounds.push_back(edgeSampleStats.nonUniformUbs);

					double supBound = edgeSampleStats.upperBoundSD;
					std::vector<double> ubsToSet;
					for(auto& ub : edgeSampleStats.nonUniformUbs)
					{
						if(ub == 0)
						{
							ubsToSet.push_back(paramsAppx.epsilon);
						}
						else if(ub < std::pow(10, -6))
						{
							ubsToSet.push_back(epsVect[0]);
						}
						else if(std::isnan(ub))
						{
							ubsToSet.push_back(paramsAppx.epsilon);
						}
						else
							ubsToSet.push_back(ub);
					}
					wsParams.nonUni =true;
					wsParams.nonUniformErrors = ubsToSet;
					wsParams.epsilon = supBound;
					std::vector<double> partitionAndWedgeSampleCoeffs = CountingAlgs::ApproximateCoeffsSota(mygraph, std::ref(wsStats), std::ref(wsParams));
					appxBL.push_back(partitionAndWedgeSampleCoeffs);
				}
				std::vector<double> runtimes = edgeSampleStats.runtime;
				std::vector<double> runtimesBL = wsStats.runtime;
				json conf = {
				{"Partition", path},
				{"Sol-CN", exactSol},
				{"Approx", approxes},
				{"Samples", samplesS},
				{"RT-ours", runtimes},
				{"fineTimes", breakTimes},
				{"TimeBL", runtimesBL},
				{"Non-uni-bounds", nonUniBounds},
				{"Approx-BL", appxBL},
				{"QValues", qvect}
				};
				confs.push_back(conf);
			}
			exps[metric] = confs;
		}
		json outA = {
			{"Dataset", fedges},
			{"N-nodes", N},
			{"M-edges", M},
			{"Results", exps}
			};
		std::ofstream outj1;
		outj1.open(outJson);
		outj1 << outA << '\n';
		outj1.close();
		return 0;
	}
	else if(mode == MODEPARALLEL) // 
	{
		std::cout << "Parallel mode\n";
		std::string fedges = exe_spec["dataset"].get<std::string>();
    	std::string outJson = exe_spec["outJson"].get<std::string>();
		H5Easy::File hdf5Data(fedges, H5Easy::File::ReadOnly);

   	 	Structures::Graph mygraph{0,0}; //Create emptyGraph
		bool oldread = exe_spec["oldread"].get<bool>();
		int threads = exe_spec["threads"].get<int>();
		Structures::BuildGraph(std::ref(hdf5Data), std::ref(mygraph), oldread);
		std::cout << "loaded graph has: " << mygraph.GetN() << " and " << mygraph.GetM() << '\n';
    	auto N = mygraph.GetN();
    	auto M = mygraph.GetM();
		long long int maxwork=exe_spec["worksmall"].get<long long int>();
		long long int samples=exe_spec["samples"].get<long long int>();

		Structures::InputParams paramsAppx;
		paramsAppx.epsilon = 0.05;
		paramsAppx.eta = 0.01;
		paramsAppx.maxwork = maxwork;
		paramsAppx.iterations = samples;
		paramsAppx.samples = samples;
		paramsAppx.steps = 100;
		paramsAppx.qval = 0.01;
		paramsAppx.mode = MODEFIXEDSS;

		std::vector<std::string> pathsPartitions = exe_spec["paths"].get<std::vector<std::string>>();
		std::vector<int> th_vec{1,2,4,8,16};
		int execs = 3;

		long long int Cseeds[10] {12,320,456,32345,2321,858,9048934,38475984,83128391,42};
		std::vector<std::string> metrics{"/closureCoeffs", "/clusteringCoeffs"};
		std::unordered_map<std::string, std::vector<json>> exps;
		int meth = 0;
		for(auto& metric: metrics) // Either clustcoeff or closure coeff
		{
			std::cout << " metric: " << metric << '\n';
			paramsAppx.normMetricVertices = H5Easy::load<std::vector<double>>(hdf5Data, metric + "/normNodes");// Load denumerator for each vertex
			Structures::Statistics edgeSampleStats;
			std::vector<json> confs;
			for(auto& path : pathsPartitions)
			{
				std::vector<double> exactValues = H5Easy::load<std::vector<double>>(hdf5Data, metric + "/values");// Load exact metric values

				auto exactSol = Structures::LoadExactSolPartition(std::ref(hdf5Data), path, std::ref(exactValues));
				auto buckets = Structures::LoadPartition(std::ref(hdf5Data), path);

				edgeSampleStats.buckets = buckets;
				edgeSampleStats.nodeToBucket = H5Easy::load<std::vector<long long int>>(hdf5Data, path);
				for(auto& th : th_vec)
				{
					std::vector<double> peelT;
					std::vector<double> optQT;
					std::vector<double> sampleT;
					std::vector<double> totT;
					paramsAppx.threads = th;
					for(int iter=0; iter<execs; iter++)
					{
						edgeSampleStats.seed = Cseeds[iter];
						auto mygraphCopy = mygraph;
						auto wsParams = paramsAppx;
						if(metric == "/closureCoeffs")
							meth = 1;
						else
							meth = 0;
						wsParams.clustOrClosure = meth;
						std::vector<double> adaptiveCoeffs = CountingAlgs::TriadPar(mygraphCopy, edgeSampleStats, std::ref(paramsAppx));
						peelT.push_back(edgeSampleStats.timepeel);
						optQT.push_back(edgeSampleStats.timevar);
						sampleT.push_back(edgeSampleStats.timeadapt);
						totT.push_back(edgeSampleStats.runtime.back());
					}

					json currconf = {
						{"Partition", path},
						{"threads", th},
						{"timePeel", peelT},
						{"timeOpt", optQT},
						{"timeSample", sampleT},
						{"totRuntime", totT},
					};
					confs.push_back(currconf);
				}
			}
			exps[metric] = confs;
		}
		json outA = {
			{"Dataset", fedges},
			{"N-nodes", N},
			{"M-edges", M},
			{"Results", exps}
			};
		std::ofstream outj1;
		outj1.open(outJson);
		outj1 << outA << '\n';
		outj1.close();
		return 0;	
	}
}
