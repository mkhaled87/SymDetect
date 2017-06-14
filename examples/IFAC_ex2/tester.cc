/*
 * tester.cc
 *
 *  created on: 11.05.2017
 *      author: M. Khaled
 */


#include <array>
#include <iostream>

#include "cuddObj.hh"

#include "nfts.hh"
#include "SymbolicDetector.hh"

using namespace std;
using namespace scots;

#define SS_DIM 1
#define IS_DIM 1
#define OS_DIM 1

typedef array<double,SS_DIM> state_type;
typedef array<double,IS_DIM> input_type;
typedef array<double,OS_DIM> output_type;

std::vector<state_type> sys_post(const state_type &x, const input_type &u)
{
	std::vector<state_type> ret;

	if(x[0] == 10 && u[0] == 0) {
		state_type x_post1,x_post2;
		x_post1[0] = 20;
		x_post2[0] = 30;
		ret.push_back(x_post1);
		ret.push_back(x_post2);
	}
	if(x[0] == 10 && u[0] == 1) {
		state_type x_post1;
		x_post1[0] = 10;
		ret.push_back(x_post1);
	}

	if(x[0] == 20 && u[0] == 0) {
		state_type x_post1;
		x_post1[0] = 20;
		ret.push_back(x_post1);
	}
	if(x[0] == 20 && u[0] == 1) {
		state_type x_post1;
		x_post1[0] = 20;
		ret.push_back(x_post1);
	}


	if(x[0] == 30 && u[0] == 0) {
		state_type x_post1;
		x_post1[0] = 30;
		ret.push_back(x_post1);
	}
	
	if(x[0] == 30 && u[0] == 1) {
		state_type x_post1;
		x_post1[0] = 20;
		ret.push_back(x_post1);
	}
	return ret;
}

/* computation of the growth bound (the result is stored in r)  */
void r_post (state_type &r, const input_type &u)
{
	r[0] = 0*u[0];
}

void out_of_state(output_type &y, const state_type &x, const input_type &u)
{
	if(x[0] == 10 && u[0] == 0) { y[0] = 0; return; }
	if(x[0] == 10 && u[0] == 1) { y[0] = 0; return; }

	if(x[0] == 20 && u[0] == 0) { y[0] = 1; return; }
	if(x[0] == 20 && u[0] == 1) { y[0] = 1; return; }
	
	if(x[0] == 30 && u[0] == 0) { y[0] = 1; return; }
	if(x[0] == 30 && u[0] == 1) { y[0] = 1; return; }
}

int main() {
  /* there is one unique manager to organize the bdd variables */
  Cudd mgr;

  double sslb[SS_DIM] ={10};
  double ssub[SS_DIM] ={30};
  double sseta[SS_DIM]={10};

  double islb[IS_DIM] ={0};
  double isub[IS_DIM] ={1};
  double iseta[IS_DIM]={1};

  double oslb[OS_DIM] ={0};
  double osub[OS_DIM] ={1};
  double oseta[OS_DIM]={1};

  SymbolicSet ss(mgr,SS_DIM,sslb,ssub,sseta);
  ss.addGridPoints();
  ss.writeToFile("test_ss.bdd");

  SymbolicSet is(mgr,IS_DIM,islb,isub,iseta);
  is.addGridPoints();
  is.writeToFile("test_is.bdd");

  SymbolicSet sspost(mgr,SS_DIM,sslb,ssub,sseta);
  sspost.addGridPoints();
  sspost.writeToFile("test_sspost.bdd");

  SymbolicSet os(mgr,OS_DIM,oslb,osub,oseta);
  os.addGridPoints();
  os.writeToFile("test_os.bdd");

  NFTS<state_type, input_type, output_type> outSystem(&ss,&is,&sspost,&os);
  outSystem.buildNFTS(sys_post, r_post, out_of_state);

  std:: cout << "writing the constructed relations: " << std::endl;
  outSystem.getTransitionRelation().writeToFile("test_rel.bdd");
  outSystem.getOutputMap().writeToFile("test_outmap.bdd");

  SymbolicDetector<state_type, input_type, output_type> detector(&outSystem, (int)pow(ss.getSize(),2) + 1, (int)((os.getSize()+1)*is.getSize()));

  std:: cout << "constructing the detector: ";
  detector.constructExactDetectorNFA(1);
  std:: cout << " done !!" << std::endl << std::endl;

  const NFA<BDD, ioTuple<input_type, output_type>>* nfa = detector.getNFA();

  std:: cout << "Priniting states in the detector's NFA: " << std::endl;
  scots::SymbolicSet stateContainer(ss);
  std::vector<BDD> nfaStateValues = nfa->getStateValues();
  for(size_t i=0; i<nfaStateValues.size(); i++){
	  BDD nfaStates =  nfaStateValues[i];
	  stateContainer.setSymbolicSet(nfaStates);
	  std::cout << stateContainer.elementsToString();
  }
  std:: cout << std::endl;
  std:: cout << std::endl;


  std:: cout << "Priniting inputs in the detector's NFA: " << std::endl;
  std::vector<ioTuple<input_type, output_type>> nfaInputValues = nfa->getInputValues();
  for(size_t i=0; i<nfaInputValues.size(); i++){
	  std::cout << "(";
	  ioTuple<input_type, output_type> nfaIo = nfaInputValues[i];
	  std::cout << "(";
	  if(!nfaIo.isBlankInput()){
		  for(size_t j=0; j<IS_DIM; j++){
			  std::cout << nfaIo.getInput()[j];
			  if(j != (IS_DIM-1))
				  std::cout << ", ";
		  }
	  }
	  else{
		  std::cout << "-";
	  }
      std::cout << "), ";

	  std::cout << "(";
	  for(size_t j=0; j<OS_DIM; j++){
		  std::cout << nfaIo.getOutput()[j];
		  if(j != (OS_DIM-1))
			  std::cout << ", ";
	  }
      std::cout << ")";
      std::cout << ")";
      if(i != (nfaInputValues.size()-1))
    	  std::cout << ", ";
  }
  std::cout << std::endl;
  std::cout << std::endl;

  std:: cout << "Priniting transitions in the detector's NFA: " << std::endl;
  std::vector<NFA_transition<BDD, ioTuple<input_type, output_type>>>  nfaTransitions = nfa->dumpTransitions();

  for(size_t i=0; i<nfaTransitions.size(); i++){
	  BDD state     = nfaTransitions[i].getSourceState();
	  BDD postState = nfaTransitions[i].getPostState();
	  ioTuple<input_type, output_type> input = nfaTransitions[i].getInput();

	  stateContainer.setSymbolicSet(state);
	  std::string strState = stateContainer.elementsToString();

	  stateContainer.setSymbolicSet(postState);
	  std::string strPostState = stateContainer.elementsToString();

	  std::stringstream ss;
	  ss << "(";
	  if(!input.isBlankInput()){
		  for(size_t j=0; j<IS_DIM; j++){
			  ss << input.getInput()[j];
			  if(j != (IS_DIM-1))
				  ss << ", ";
		  }
	  }
	  else{
		  ss << "-";
	  }
	  ss << "),(";
	  for(size_t j=0; j<OS_DIM; j++){
		  ss << input.getOutput()[j];
		  if(j != (OS_DIM-1))
			  ss << ", ";
	  }
	  ss << "))";
	  std::string strInput = ss.str();

	  std::cout << strState << " |--------- \t" << strInput << "\t --------> " << strPostState << std::endl;

  }


  return 0;
}
