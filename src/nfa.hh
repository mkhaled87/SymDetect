/*
 * nfa.hh
 *
 *  created on: 26.04.2017
 *      author: M. Khaled
 */

#ifndef SCOTSNFA_HH_
#define SCOTSNFA_HH_

#include <vector>
#include <algorithm>
#include <sstream>
#include <cassert>
#include <cmath>
#include <iostream>
#include <fstream>
#include <stdexcept>
#include <math.h>

#include "dddmp.h"
#include "cuddObj.hh"
#include "BDDTricks.hh"

/* class: NFA
 *
 * bdd-based non-deterministic finite automata implementation
 *
 *
 * Properties: 
 * - 
 *
 *
 */
namespace scots {

template<class stateType_, class inputType_>
class NFA_transition{
private:
	stateType_ sourceState;
	inputType_ input;
	stateType_ postState;
public:
	NFA_transition(stateType_ sourceState_, inputType_ input_, stateType_ postState_):
		sourceState(sourceState_), input(input_), postState(postState_)
	{}

	stateType_ getSourceState(){return sourceState;}
	stateType_ getPostState(){return postState;}
	inputType_ getInput(){return input;}
};

template<class stateType_, class inputType_>
class NFA {
private:
	Cudd* pMgr;			// the BDD manager
	bool ownManager = true;		// whether we have our own manager;

	std::vector<size_t> statesIdxBddVars;
	std::vector<size_t> inputsIdxBddVars;
	std::vector<size_t> postStatesIdxBddVars;

	std::vector<BDD> statesBddVars;
	std::vector<BDD> inputsBddVars;
	std::vector<BDD> postStatesBddVars;

	std::vector<stateType_> stateValues;
	std::vector<inputType_> inputValues;

	size_t maxStates;
	size_t maxInputs;

	size_t idxNextState = 0;
	size_t idxNextInput = 0;

	BDD states;			// to hold states
	BDD inputs;			// to hold the inputs
	BDD postStates;		// to hold post-states
	BDD transitions;		// to hold

public:

	/* constructor: NFA
	 */
	NFA(size_t max_num_states = ((size_t) -1), size_t max_num_inputs =
			((size_t) -1)) {
		if (ownManager) {
			pMgr = new Cudd();
		}

		if (max_num_states <= 0 || max_num_inputs <= 0) {
			throw std::invalid_argument(
					"Error: scots::NFA: number of states/input should be > 0.");
		}

		maxStates = max_num_states;
		maxInputs = max_num_inputs;

		for (size_t i = 0; i < log2(max_num_states); i++) {
			BDD var = pMgr->bddVar();
			statesBddVars.push_back(var);
			statesIdxBddVars.push_back(var.NodeReadIndex());
		}
		for (size_t i = 0; i < log2(max_num_inputs); i++) {
			BDD var = pMgr->bddVar();
			inputsBddVars.push_back(var);
			inputsIdxBddVars.push_back(var.NodeReadIndex());
		}
		for (size_t i = 0; i < log2(max_num_states); i++) {
			BDD var = pMgr->bddVar();
			postStatesBddVars.push_back(var);
			postStatesIdxBddVars.push_back(var.NodeReadIndex());
		}

		states = pMgr->bddZero();
		postStates = pMgr->bddZero();
		inputs = pMgr->bddZero();
		transitions = pMgr->bddZero();

	}

	/* constructor: NFA :: assigning a CUDD manager
	 */
	NFA(Cudd mgr, size_t max_num_states = ((size_t) -1), size_t max_num_inputs =
			((size_t) -1)) :
			NFA(max_num_states, max_num_inputs) {
		pMgr = &mgr;
		ownManager = false;
	}

	~NFA() {
		if (ownManager) {
			delete pMgr;
		}
	}

	/* function: stateIdxToPhase :: produces a BDD from the index of state
	 */
	BDD stateIdxToBDD(size_t idx) {

		if (idx >= maxStates) {
			throw std::invalid_argument(
					"Error: scots::NFA: state idx must be < maximum number of states.");
		}

		std::vector<int> phase(statesBddVars.size());
		for (size_t i = 0; i < statesBddVars.size(); i++)
			phase[i] = 0;
		for (size_t i = 0; idx; idx /= 2, i++)
			phase[i] = (idx % 2);

		return pMgr->bddComputeCube(statesBddVars.data(), phase.data(),
				statesBddVars.size());
	}

	/* function: postStateIdxToPhase :: produces a BDD from the index of post state
	 */
	BDD postStateIdxToBDD(size_t idx) {

		if (idx >= maxStates) {
			std::cout << "idx= " << idx << std::endl;
			throw std::invalid_argument(
					"Error: scots::NFA: post state idx must be < maximum number of states.");
		}

		std::vector<int> phase(postStatesBddVars.size());
		for (size_t i = 0; i < postStatesBddVars.size(); i++)
			phase[i] = 0;
		for (size_t i = 0; idx; idx /= 2, i++)
			phase[i] = (idx % 2);

		return pMgr->bddComputeCube(postStatesBddVars.data(), phase.data(),
				postStatesBddVars.size());
	}

	/* function: inputIdxToPhase :: produces a BDD from the index of input
	 */
	BDD inputIdxToBDD(size_t idx) {

		if (idx >= maxInputs) {
			throw std::invalid_argument(
					"Error: scots::NFA: input idx must be < maximum number of inputs.");
		}

		std::vector<int> phase(inputsBddVars.size());
		for (size_t i = 0; i < inputsBddVars.size(); i++)
			phase[i] = 0;
		for (size_t i = 0; idx; idx /= 2, i++)
			phase[i] = (idx % 2);

		return pMgr->bddComputeCube(inputsBddVars.data(), phase.data(),
				inputsBddVars.size());
	}

	/* function: addState :: adds a new state to the NFA
	 */
	size_t addState(stateType_ value) {
		BDD state = stateIdxToBDD(idxNextState);
		states += state;

		BDD postState = postStateIdxToBDD(idxNextState);
		postStates += postState;

		stateValues.push_back(value);

		idxNextState++;
		return (idxNextState - 1);
	}

	/* function: addInput :: adds a new input to the NFA
	 */
	size_t addInput(inputType_ value) {
		BDD input = inputIdxToBDD(idxNextInput++);
		inputs += input;

		inputValues.push_back(value);

		return (idxNextInput - 1);
	}

	inline size_t countStates(BDD states) const {
		if (!BddTricks::ckechAgainstSuportVars(states, statesIdxBddVars)) {
			throw std::invalid_argument(
					"Error: scots::NFA: BDD has support one or more vars out of expected !");
		}
		return states.CountMinterm(statesBddVars.size());
	}

	inline size_t countPostStates(BDD states)  const {
		if (!BddTricks::ckechAgainstSuportVars(states, postStatesIdxBddVars)) {
			throw std::invalid_argument(
					"Error: scots::NFA: BDD has support one or more vars out of expected !");
		}
		return states.CountMinterm(postStatesBddVars.size());
	}

	inline size_t countInputs(BDD inputs)  const {
		if (!BddTricks::ckechAgainstSuportVars(inputs, inputsIdxBddVars)) {
			throw std::invalid_argument(
					"Error: scots::NFA: BDD has support one or more vars out of expected !");
		}
		return inputs.CountMinterm(inputsBddVars.size());
	}

	BDD addTransition(size_t stateIdx, size_t inputIdx, size_t postStateIdx) {

		BDD state = stateIdxToBDD(stateIdx);
		BDD input = inputIdxToBDD(inputIdx);
		BDD postState = postStateIdxToBDD(postStateIdx);

		if (state == pMgr->bddZero() || countStates(state) > 1) {
			throw std::invalid_argument(
					"Error: scots::NFA: BDD zero state or BDD has more than one state !");
		}

		if (input == pMgr->bddZero() || countInputs(input) > 1) {
			throw std::invalid_argument(
					"Error: scots::NFA: BDD zero inputs or BDD has more than one input !");
		}

		if (postState == pMgr->bddZero() || countPostStates(postState) > 1) {
			throw std::invalid_argument(
					"Error: scots::NFA: BDD zero post state or BDD has more than one post state !");
		}

		BDD transition = (state * input * postState);

		transitions += transition;
		return transition;
	}

	size_t numStates() const {
		return countStates(states);
	}
	size_t numInputs() const {
		return countInputs(inputs);
	}
	size_t numTransitions() const {
		return transitions.CountMinterm(
				statesBddVars.size() + inputsBddVars.size()
						+ postStatesBddVars.size());
	}

	stateType_ getStateValue(size_t idx) const{
		return stateValues[idx];
	}
	std::vector<stateType_> getStateValues() const{
		return stateValues;
	}

	inputType_ getInputValue(size_t idx) const{
		return inputValues[idx];
	}
	std::vector<inputType_> getInputValues() const{
		return inputValues;
	}

	BDD getTransitionsBDD() const{
		return transitions;
	}

	std::vector<NFA_transition<stateType_, inputType_>> dumpTransitions() const{
		int minterm[statesBddVars.size() + inputsBddVars.size() + postStatesBddVars.size()];

		std::vector<size_t> allvars;
		VectorTricks::AppendVector(allvars, statesIdxBddVars);
		VectorTricks::AppendVector(allvars, inputsIdxBddVars);
		VectorTricks::AppendVector(allvars, postStatesIdxBddVars);

		std::vector<NFA_transition<stateType_, inputType_>> ret;
		CuddMintermIterator Ielements(transitions, allvars, allvars.size());
		while(!Ielements.done()){
			Ielements.copyMinterm(minterm);

			size_t stateIdx=0;
		    for(size_t c=1, i=0; i<statesIdxBddVars.size(); c*=2, i++)
		      stateIdx+=minterm[statesIdxBddVars[i]]*c;

			size_t postStateIdx=0;
		    for(size_t c=1, i=0; i<postStatesIdxBddVars.size(); c*=2, i++)
		    	postStateIdx+=minterm[postStatesIdxBddVars[i]]*c;

			size_t inputIdx=0;
		    for(size_t c=1, i=0; i<inputsIdxBddVars.size(); c*=2, i++)
		    	inputIdx+=minterm[inputsIdxBddVars[i]]*c;


		    NFA_transition<stateType_, inputType_> trans(getStateValue(stateIdx), getInputValue(inputIdx), getStateValue(postStateIdx));
		    ret.push_back(trans);

			++Ielements;
		}

		return ret;
	}
	std::vector<stateType_> getPost(stateType_ s0, std::vector<inputType_> appliedInputs) {
		std::vector<stateType_> ret;
		
		size_t s0Idx;
		if(!getStateByValue(s0, &s0Idx))
			throw std::runtime_error("NFA::getPost: Invalid source state !");
		
		BDD srcState = stateIdxToBDD(s0Idx);
		for (size_t i = 0; i < appliedInputs.size(); i++) {

			size_t inpIdx;
			if(!getInputByValue(appliedInputs[i], &inpIdx))
				throw std::runtime_error("NFA::getPost: Invalid input !");
			
			BDD input = inputIdxToBDD(inpIdx);
			BDD post = transitions*srcState*input;

			if(post == pMgr->bddZero()){
				std::cout << "NFA::getPost: (warninig) input sequence leads to no posts !";
				srcState=pMgr->bddZero();
				break;
			}

			post = BddTricks::ProjectBDD(*pMgr, post, postStatesIdxBddVars);
			srcState = BddTricks::PermutreBdd(*pMgr, post, postStatesIdxBddVars, statesIdxBddVars);
		}

		srcState = BddTricks::ProjectBDD(*pMgr, srcState, statesIdxBddVars);

		int minterm[statesIdxBddVars.size()];
		CuddMintermIterator Istates(srcState, statesIdxBddVars, statesIdxBddVars.size());
		while(!Istates.done()){
			Istates.copyMinterm(minterm);

			size_t stateIdx=0;
		    	for(size_t c=1, i=0; i<statesIdxBddVars.size(); c*=2, i++)
		      		stateIdx+=minterm[statesIdxBddVars[i]]*c;

		        BDD state = getStateValue(stateIdx);
		        ret.push_back(state);

			++Istates;
		}

		return ret;
	}

	bool getStateByValue(const stateType_ stateValue, size_t* outIdx) {
		for (size_t i = 0; i < stateValues.size(); i++) {
			if (stateValues[i] == stateValue) {
				*outIdx = i;
				return true;
			}
		}
		return false;
	}
	bool getInputByValue(const inputType_ inputValue, size_t* outIdx) {
		for (size_t i = 0; i < inputValues.size(); i++) {
			if (inputValues[i] == inputValue) {
				*outIdx = i;
				return true;
			}
		}
		return false;
	}

};
/* close class def */
} /* close namespace */

#endif /* SCOTSNFA_ */
