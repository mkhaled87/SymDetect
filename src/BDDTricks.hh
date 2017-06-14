/*
 * BDDTricks.hh
 *
 *  created on: 03.05.2017
 *      author: M.Khaled
 */

#ifndef BDDDTRICKS_HH_
#define BDDDTRICKS_HH_

#include <array>
#include <vector>

#include "cuddObj.hh"
#include "VectorTricks.hh"


class BddTricks {
public:

	/* returns an empty permutation map ready for manipulation */
	static
	std::vector<int> GetReadyPermuteMap(const Cudd& cuddManager){
		size_t nAllBddVars = cuddManager.ReadSize();
		std::vector<int> permute(nAllBddVars);
		for(size_t i=0; i<nAllBddVars; i++) permute[i]=i;
		return permute;
	}

	/* Permutes a given BDD from original transition relation to the new order of transition relation*/
	static
	BDD PermutreBdd(const Cudd& cuddManager, BDD& orgBdd, std::vector<size_t> oldBddVars, std::vector<size_t> newBddVars){

		if(oldBddVars.size() != newBddVars.size())
			throw std::runtime_error("Error::PermutreBddToNewVars:: Mismatch in sizes of bdd-vectors passed !");

		// Permute f_org_rel to reflect the new BddVars
		std::vector<int> permute = BddTricks::GetReadyPermuteMap(cuddManager);

		for(size_t i=0; i<oldBddVars.size(); i++){
			permute[oldBddVars[i]] = newBddVars[i];
			permute[newBddVars[i]] = oldBddVars[i];
		}

		return orgBdd.Permute(permute.data());
	}

	/* returns a BDD that is a projection of the supplied BDD to the supplied BddVars,
	 * all other BDD vars are treated as dont-cares*/
	static
	BDD ProjectBDD(const Cudd& cuddManager, const BDD& srcBDD, std::vector<size_t> projVars){
		size_t nAllBddVars = cuddManager.ReadSize();
		std::vector<BDD> otherBDDs;

		for(size_t i=0; i<nAllBddVars; i++)
			if (!(find(projVars.begin(),projVars.end(),i) != projVars.end()))
				otherBDDs.push_back(cuddManager.bddVar(i));

		BDD otherBdDDsCube = cuddManager.bddComputeCube(otherBDDs.data(), NULL, otherBDDs.size());
		return srcBDD.ExistAbstract(otherBdDDsCube);
	}


	static
	void PrintBDD(std::string NameIt, const BDD& bddObj){
		std::cout << "Info for BDD: " << NameIt << std::endl;
		//cout << "\t    f : (" << bddObj << ")" << endl;
		std::cout << "       cover: " << std::endl;
		bddObj.PrintCover();
		std::cout << std::endl;
	}

	static
    	bool ckechAgainstSuportVars(BDD toTest, std::vector<size_t> vars){
		std::vector<unsigned int> sup = toTest.SupportIndices();
		std::vector<size_t> sup2;
		for(size_t i=0; i<sup.size(); i++) sup2.push_back(sup[i]);
		return VectorTricks::isVectorIncluded(sup2, vars);
	}

};
#endif
