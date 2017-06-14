/*
 * SymbolicDetector.hh
 *
 *  created on: 26.04.2017
 *      author: M. Khaled
 */

#ifndef SYMDETECTOR_HH_
#define SYMDETECTOR_HH_

#include <vector>
#include <algorithm>
#include <sstream>
#include <cassert>
#include <cmath>
#include <iostream>
#include <fstream>
#include <stdexcept>

#include "dddmp.h"
#include "cuddObj.hh"

#include "nfts.hh"
#include "nfa.hh"
#include "BDDTricks.hh"
#include "VectorTricks.hh"
#include "ExtraUtils.hh"

/* class: SymbolicDetector
 *
 * A symbolic detector that uses a symbolic model to construct an NFA
 * The NFA is checking the decttability of the symbolic model and also
 * works as a detector.
 *
 * Properties: 
 * - 
 *
 *
 */
namespace scots {

template<class inputType_, class outType_>
class ioTuple {
private:
	bool blankInput;
	bool blankOutput;
	inputType_ input;
	outType_ output;

public:
	ioTuple() :
			blankInput(true), blankOutput(true) {
	}
	ioTuple(inputType_ in, outType_ out) :
			blankInput(false), blankOutput(false), input(in), output(out) {
	}
	//ioTuple(outType_ out):blankInput(true), blankOutput(false), output(out){}
	//ioTuple(inputType_ in):blankInput(false), blankOutput(true), input(in){}

	void setAsInputOnlyTuple(inputType_ in) {
		input = in;
		blankInput = false;
		blankOutput = true;
	}

	void setAsOutputOnlyTuple(outType_ out) {
		output = out;
		blankInput = true;
		blankOutput = false;
	}

	void setTuple(inputType_ in, outType_ out) {
		input = in;
		output = out;
		blankInput = false;
		blankOutput = false;
	}

	bool isBlankInput() {
		return blankInput;
	}
	bool isBlankOutput() {
		return blankOutput;
	}

	inline inputType_ getInput() {
		return input;
	}
	inline outType_ getOutput() {
		return output;
	}

	bool operator==(ioTuple rhs) const {

		if (blankInput == true && rhs.isBlankInput() && blankOutput==true && rhs.isBlankOutput())
			return true;

		if (blankInput == true && rhs.isBlankInput() && blankOutput==false && !rhs.isBlankOutput() && ExtraUtils::doubleEquals(output.size(), output.data(), rhs.getOutput().data()))
			return true;

		if (blankInput == false && !rhs.isBlankInput() && blankOutput==true && rhs.isBlankOutput() && ExtraUtils::doubleEquals(input.size() ,input.data(), rhs.getInput().data()))
			return true;

		if (blankInput == false && !rhs.isBlankInput() && blankOutput==false && !rhs.isBlankOutput() && ExtraUtils::doubleEquals(input.size(), input.data(), rhs.getInput().data()) && ExtraUtils::doubleEquals(output.size(), output.data(), rhs.getOutput().data()))
			return true;

		return false;
	}

};

template<class stateType_, class inputType_, class outType_>
class SymbolicDetector {
protected:
	Cudd* ddmgr_;
	NFTS<stateType_, inputType_, outType_>* ntfs;
	NFA<BDD, ioTuple<inputType_, outType_>>* nfa;

public:
	/* constructor: NFA :: default
	 */
	SymbolicDetector(NFTS<stateType_, inputType_, outType_>* outSystem_) :
			ntfs(outSystem_) {
		ddmgr_ = outSystem_->getCuddManager();

		size_t nfaMaxStates = (size_t)outSystem_->getStateSpace()->getSize() + (size_t)pow(outSystem_->getStateSpace()->getSize(),2) + 1;
		size_t nfaMaxInputs = ((size_t)outSystem_->getOutSpace()->getSize()+1)*((size_t)outSystem_->getInputSpace()->getSize());
		nfa = new NFA<BDD, ioTuple<inputType_, outType_>>(nfaMaxStates, nfaMaxInputs);
	}

	SymbolicDetector(NFTS<stateType_, inputType_, outType_>* outSystem_, size_t nfaMaxStates, size_t nfaMaxInputs) :
			ntfs(outSystem_) {
		ddmgr_ = outSystem_->getCuddManager();
		nfa = new NFA<BDD, ioTuple<inputType_, outType_>>(nfaMaxStates, nfaMaxInputs);
	}

	~SymbolicDetector() {
		delete nfa;
	}

	const NFA<BDD, ioTuple<inputType_, outType_>>* getNFA(){
		return nfa;
	}

	const NFTS<stateType_, inputType_, outType_>* getNFTS(){
		return ntfs;
	}

	std::vector<BDD> constructTwoStateSubsets(BDD statesSet,
			std::vector<size_t> stateVars) {

		std::vector<BDD> stateVarsBdd;
		for (size_t i = 0; i < stateVars.size(); i++)
			stateVarsBdd.push_back(ddmgr_->bddVar(stateVars[i]));

		std::vector<BDD> subsets;

		int I1_minterm[stateVars.size()];
		int I2_minterm[stateVars.size()];

		CuddMintermIterator I1(statesSet, stateVars, stateVars.size());
		size_t i1 = 0;
		while (!I1.done()) {
			CuddMintermIterator I2(statesSet, stateVars, stateVars.size());
			size_t i2 = 0;
			while (!I2.done()) {

				if (i1 == i2 || i2 < i1) {
					++I2;
					++i2;
					continue;
				}

				I1.copyMinterm(I1_minterm);
				I2.copyMinterm(I2_minterm);

				BDD I1_bdd = ddmgr_->bddComputeCube(stateVarsBdd.data(),
						I1_minterm, stateVars.size());
				BDD I2_bdd = ddmgr_->bddComputeCube(stateVarsBdd.data(),
						I2_minterm, stateVars.size());

				subsets.push_back(I1_bdd + I2_bdd);

				++I2;
				++i2;
			}
			++I1;
			++i1;
		}
		return subsets;
	}

	bool detectByNFTS(const std::vector<BDD> appliedInputs, const std::vector<BDD> observedOutputs, stateType_* outState){
		SymbolicSet stateContainer(*ntfs->getStateSpace());
		BDD post = ntfs->getPost(ntfs->getStateSpace()->getSymbolicSet(), appliedInputs, observedOutputs);
		stateContainer.setSymbolicSet(post);
		if(stateContainer.getSize() == 1){
			// getting state bddVars ready
			std::vector<size_t> vctStateVars;
			for (size_t i = 0; i < ntfs->getStateSpace()->getDimension(); i++)
				for (size_t j = 0; j < ntfs->getStateSpace()->getNofBddVars()[i]; j++)
					vctStateVars.push_back(ntfs->getStateSpace()->getIndBddVars()[i][j]);

			int minterm[ntfs->getNumStateVars()];
			CuddMintermIterator Ix(post, vctStateVars, vctStateVars.size());
			Ix.copyMinterm(minterm);
			stateContainer.mintermToElement(minterm, outState->data());
			return true;
		}else{
			return false;
		}
	}

	std::vector<stateType_> getSingletonNFAStates(){
		// getting state bddVars ready
		std::vector<size_t> vctStateVars;
		for (size_t i = 0; i < ntfs->getStateSpace()->getDimension(); i++)
			for (size_t j = 0; j < ntfs->getStateSpace()->getNofBddVars()[i]; j++)
				vctStateVars.push_back(ntfs->getStateSpace()->getIndBddVars()[i][j]);

		std::vector<BDD> states = nfa->getStateValues();
		std::vector<stateType_> ret;

		SymbolicSet stateContainer(*ntfs->getStateSpace());
		for(size_t i=0; i<states.size(); i++){
			stateContainer.setSymbolicSet(states[i]);			
			if(stateContainer.getSize() == 1){
				int minterm[ntfs->getNumStateVars()];
				CuddMintermIterator Ix(states[i], vctStateVars, vctStateVars.size());
				Ix.copyMinterm(minterm);
				stateType_ outState;
				stateContainer.mintermToElement(minterm, outState.data());
	
				ret.push_back(outState);
			}
		}
		return ret;
	}

	bool detectByNFA(const std::vector<inputType_> appliedInputs, const std::vector<outType_> observedOutputs, stateType_* outState, int verbosity =0){

		if (observedOutputs.size() != (appliedInputs.size()+1)) {
			std::cout
					<< "SymbolicDetector::dtetctByNFA (error): mismatch of size btw input/output arrays!"
					<< std::endl;
			return false;
		}
		if(verbosity > 0 ) std::cout << "Detecting by NFA: " << std::endl;
		std::vector<ioTuple<inputType_, outType_>> ioPairs;

		ioTuple<inputType_, outType_> firstIoPair;
		firstIoPair.setAsOutputOnlyTuple(observedOutputs[0]);
		ioPairs.push_back(firstIoPair);
		for(size_t i=0; i<appliedInputs.size(); i++){
			ioTuple<inputType_, outType_> tmp(appliedInputs[i], observedOutputs[i+1]);
			ioPairs.push_back(tmp);
		}

		SymbolicSet stateContainer(*ntfs->getStateSpace());
		if(verbosity > 0 ) std::cout << "\tsimulating transitions from empty state: " << std::endl;
		std::vector<BDD> post = nfa->getPost(ddmgr_->bddZero(), ioPairs);
		if(verbosity > 0 ) std::cout << "\tNFA post size: " << post.size() << std::endl;
		if(post.size() != 1)
			return false;

		stateContainer.setSymbolicSet(post[0]);
		if(verbosity > 0 ) std::cout << "\tNFA post size: " << post.size() << std::endl;
		if(stateContainer.getSize() == 1){
			// getting state bddVars ready
			std::vector<size_t> vctStateVars;
			for (size_t i = 0; i < ntfs->getStateSpace()->getDimension(); i++)
				for (size_t j = 0; j < ntfs->getStateSpace()->getNofBddVars()[i]; j++)
					vctStateVars.push_back(ntfs->getStateSpace()->getIndBddVars()[i][j]);

			int minterm[ntfs->getNumStateVars()];
			CuddMintermIterator Ix(post[0], vctStateVars, vctStateVars.size());
			Ix.copyMinterm(minterm);
			stateContainer.mintermToElement(minterm, outState->data());
			return true;
		}else{
			return false;
		}
	}

	void constructExactDetectorNFA(int verbosity=0) {

		/* init */
		/******************************************/

		// getting output map and output vars ready
		BDD H = ntfs->getOutputMap().getSymbolicSet();
		SymbolicSet outContainer(*ntfs->getOutSpace());
		std::vector<size_t> vctOutVars(ntfs->getNumOutVars());
		std::vector<BDD> vctOutVarsBdd;
		for (size_t i = 0; i < ntfs->getNumOutVars(); i++) {
			vctOutVars[i] = ntfs->getOutVars()[i];
			vctOutVarsBdd.push_back(ddmgr_->bddVar(vctOutVars[i]));
		}
		BDD Y = BddTricks::ProjectBDD(*ddmgr_, H, vctOutVars);

		// getting state bddVars ready
		std::vector<size_t> vctStateVars;
		std::vector<BDD> vctStateVarsBdd;
		SymbolicSet stateContainer(*ntfs->getStateSpace());
		for (size_t i = 0; i < ntfs->getStateSpace()->getDimension(); i++) {
			for (size_t j = 0; j < ntfs->getStateSpace()->getNofBddVars()[i];
					j++) {
				vctStateVars.push_back(
						ntfs->getStateSpace()->getIndBddVars()[i][j]);
				vctStateVarsBdd.push_back(
						ddmgr_->bddVar(
								ntfs->getStateSpace()->getIndBddVars()[i][j]));
			}
		}

		// getting input bddVars ready
		std::vector<size_t> vctInputVars;
		std::vector<BDD> vctInputVarsBdd;
		SymbolicSet inputContainer(*ntfs->getInputSpace());
		for (size_t i = 0; i < ntfs->getInputSpace()->getDimension(); i++) {
			for (size_t j = 0; j < ntfs->getInputSpace()->getNofBddVars()[i];
					j++) {
				vctInputVars.push_back(
						ntfs->getInputSpace()->getIndBddVars()[i][j]);
				vctInputVarsBdd.push_back(
						ddmgr_->bddVar(
								ntfs->getInputSpace()->getIndBddVars()[i][j]));
			}
		}

		BDD U = BddTricks::ProjectBDD(*ddmgr_, H, vctInputVars);

		// vars and init of the NFA detector
		std::vector<BDD> Q1, Q2;
		size_t q0 = nfa->addState(ddmgr_->bddZero());	// the init state

		// Minterm tricks
		int global_minterm[ntfs->getNumStateVars() + ntfs->getNumInputVars()
				+ ntfs->getNumStateVars() + ntfs->getNumOutVars()];
		int* pStateMinterm = global_minterm;
		int* pInpMinterm = global_minterm + ntfs->getNumStateVars();
		int* pOutMinterm = global_minterm
				+ (ntfs->getNumStateVars() + ntfs->getNumInputVars()
						+ ntfs->getNumStateVars());

		/* (1) for-all y in Y */
		/***********************************************/
		if(verbosity > 0) std::cout << "Constructing an exact detector: " << std::endl;
		if(verbosity > 0) std::cout << "1] Main iteration over all outputs " << std::endl;
		if(verbosity > 1) {std::cout << "\tOutput Vars: "; VectorTricks::PrintVector(vctOutVars);}
		CuddMintermIterator Iy(Y, vctOutVars, vctOutVars.size());
		while (!Iy.done()) {

			// y as bdd, minterm and as value
			outType_ yValue;

			Iy.copyMinterm(pOutMinterm);
			BDD y = ddmgr_->bddComputeCube(vctOutVarsBdd.data(), pOutMinterm,
					ntfs->getNumOutVars());
			outContainer.mintermToElement(global_minterm, yValue.data());

			if(verbosity > 1){
				std::cout << "\tcurrent output: ";
				VectorTricks::PrintArray(yValue.data(), yValue.size());
			}
				
			BDD X_y = BddTricks::ProjectBDD(*ddmgr_, H * y, vctStateVars);
			stateContainer.setSymbolicSet(X_y);

			size_t inpIdx, stateIdx;
			if (stateContainer.getSize() == 1) {
				if (!VectorTricks::isVectorElement(Q1, X_y))
					Q1.push_back(X_y);				// Q1     = Q1      U {X_y}

				if (!nfa->getStateByValue(X_y, &stateIdx))
					stateIdx = nfa->addState(X_y);

				ioTuple<inputType_, outType_> tmp;
				tmp.setAsOutputOnlyTuple(yValue);
				if (!nfa->getInputByValue(tmp, &inpIdx))
					inpIdx = nfa->addInput(tmp);// SIGMA  = SIGMA   U {(-,y)}

				nfa->addTransition(q0, inpIdx, stateIdx);// DELTA  = DELTA   U {(q0,(-,y),X_y)}
			} else if (stateContainer.getSize() > 1) {
				std::vector<BDD> twoStateSubsets = constructTwoStateSubsets(X_y,
						vctStateVars);

				for (size_t i = 0; i < twoStateSubsets.size(); i++)
					if (!VectorTricks::isVectorElement(Q1, twoStateSubsets[i]))
						Q1.push_back(twoStateSubsets[i]);	// Q1     = Q1 U Zs

				ioTuple<inputType_, outType_> tmp;
				tmp.setAsOutputOnlyTuple(yValue);
				if (!nfa->getInputByValue(tmp, &inpIdx))
					inpIdx = nfa->addInput(tmp);// SIGMA  = SIGMA   U {(-,y)}

				for (size_t i = 0; i < twoStateSubsets.size(); i++) {
					if (!nfa->getStateByValue(twoStateSubsets[i], &stateIdx))
						stateIdx = nfa->addState(twoStateSubsets[i]);

					nfa->addTransition(q0, inpIdx, stateIdx); // DELTA  = DELTA   U {(q0,(-,y),X)}
				}
			}

			if(verbosity > 0) Iy.printProgress();
			++Iy;			
		}
		VectorTricks::AppendVector(Q2, Q1);
		Q1.clear();

		if(verbosity > 0) std::cout << std::endl << "2] Second iteration ";
		/* (2)-(3)*/
		/***********************************************/
		while (Q2.size() != 0) {
			if(verbosity > 1) std::cout << std::endl << "Q2:" << Q2.size() << ": ";

			// for each q2 \in Q2
			for (size_t i = 0; i < Q2.size(); i++) {
				if(verbosity > 1) {
					static size_t last_progress=0;
					size_t prog = ((double)i/(double)Q2.size())*80;
					if(prog>last_progress || Q2.size()<80){
						std::cout << ".";
						std::flush(std::cout);
						last_progress = prog;
					}
				}
				BDD q2 = Q2[i];
				size_t q2Idx;
				if(!nfa->getStateByValue(q2, &q2Idx))
					throw std::runtime_error("Invalid element in the NFA !");

				// for each state x \in q2
				CuddMintermIterator Ix(q2, vctStateVars, vctStateVars.size());
				while (!Ix.done()) {
					Ix.copyMinterm(pStateMinterm);
					BDD x = ddmgr_->bddComputeCube(vctStateVarsBdd.data(),
							pStateMinterm, vctStateVarsBdd.size());
					BDD y0 = BddTricks::ProjectBDD(*ddmgr_, H * x, vctOutVars);

					// for each y \in Y
					CuddMintermIterator Iy(Y, vctOutVars, vctOutVars.size());
					while (!Iy.done()) {
						Iy.copyMinterm(pOutMinterm);
						BDD y = ddmgr_->bddComputeCube(vctOutVarsBdd.data(),
								pOutMinterm, vctOutVarsBdd.size());
						outType_ yValue;
						outContainer.mintermToElement(global_minterm,
								yValue.data());

						// foreach u \in U
						CuddMintermIterator Iu(U, vctInputVars,
								vctInputVars.size());
						while (!Iu.done()) {
							Iu.copyMinterm(pInpMinterm);
							BDD u = ddmgr_->bddComputeCube(
									vctInputVarsBdd.data(), pInpMinterm,
									vctInputVarsBdd.size());
							inputType_ uValue;
							inputContainer.mintermToElement(global_minterm,
									uValue.data());

							std::vector<BDD> appliedInputs;
							appliedInputs.push_back(u);

							std::vector<BDD> observedOutputs;
							observedOutputs.push_back(y0);
							observedOutputs.push_back(y);

							BDD postStates = ntfs->getPost(q2, appliedInputs,
									observedOutputs);
							stateContainer.setSymbolicSet(postStates);

							size_t inpIdx, stateIdx;
							if (stateContainer.getSize() == 1) {
								ioTuple<inputType_, outType_> tmp(uValue,
										yValue);
								if (!nfa->getInputByValue(tmp, &inpIdx))
									inpIdx = nfa->addInput(tmp); // SIGMA  = SIGMA   U {(-,y)}

								if (!nfa->getStateByValue(postStates,
										&stateIdx)) {
									stateIdx = nfa->addState(postStates);
									Q1.push_back(postStates);
								}

								nfa->addTransition(q2Idx, inpIdx, stateIdx); // DELTA  = DELTA   U {(q0,(-,y),X_y)}
							} else if (stateContainer.getSize() > 1) {
								std::vector<BDD> twoStateSubsets =
										constructTwoStateSubsets(postStates,
												vctStateVars);

								ioTuple<inputType_, outType_> tmp(uValue,
										yValue);
								if (!nfa->getInputByValue(tmp, &inpIdx))
									inpIdx = nfa->addInput(tmp); // SIGMA  = SIGMA   U {(-,y)}

								for (size_t i = 0; i < twoStateSubsets.size();
										i++) {
									if (!nfa->getStateByValue(
											twoStateSubsets[i], &stateIdx)) {
										stateIdx = nfa->addState(
												twoStateSubsets[i]);
										Q1.push_back(twoStateSubsets[i]); // Q1     = Q1 U Zs
									}



									nfa->addTransition(q2Idx, inpIdx, stateIdx); // DELTA  = DELTA   U {(q0,(-,y),X)}
								}

							}

							++Iu;
						}
						++Iy;
					}
					++Ix;
				}
			}
			Q2 = Q1;
			Q1.clear();
		}
		if(verbosity > 0) std::cout << std::endl << "Finished constructing the NFA detector !" << std::endl;
	}

};
/* close class def */
} /* close namespace */

#endif /* SYMDETECTOR_HH_ */
