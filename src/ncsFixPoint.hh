/*
 * FixedPoint.hh
 *
 *  created on: 09.10.2015
 *      author: rungger
 */

#ifndef ncsFixPoint_HH_
#define ncsFixPoint_HH_

#include <array>
#include <iostream>
#include <stdexcept>

#include "cuddObj.hh"
#include "../asncs/Misc.hh"


/*
 * class: ncsFixPoint
 *
 * it provides fixed point computations to
 * synthesize controllers for NCS Models provided as BDD objects
 *
 */

class ncsFixPoint {
protected:
  /* var: ddmgr_ */
  Cudd* ddmgr_;

  /* var: permute 
   * stores the permutation array used to swap pre with post variables */
  int* permute_;

  /*vars: nisVars_, nssVars_*/
  /* store number of vars in the ss and is */
  size_t nisVars_, nssVars_;

  /*vars: postVars_, preVars_, inpVars_*/
  /* store bddVars for pre, post and input*/
  std::vector<size_t> postVars_;
  std::vector<size_t> preVars_;
  std::vector<size_t> inpVars_;

  /* helper BDDs */
  /* transition relation */
  BDD R_;

  /* transition relation with cubePost_ abstracted */
  BDD RR_;  

   /* cubes with input and post variables; used in the existential abstraction  */
  BDD cubePost_;
  BDD cubeInput_;
  
public:

  /* function: ncsFixPoint
   *
   * initialize the ncsFixPoint object with a <SymbolicModel> containing the
   * transition relation
   */
  ncsFixPoint(Cudd* cuddManager, BDD& transRelation, std::vector<size_t>& preVars, std::vector<size_t>& inpVars, std::vector<size_t>& postVars) {

	  if(postVars.size() != preVars.size()){
		  std::ostringstream os;
		  os << "Error: ncsFixPoint: preVars and postVars should have same number of bddVars.";
		  throw std::invalid_argument(os.str().c_str());
	  }


    ddmgr_ = cuddManager;
    nssVars_ = preVars.size();
    nisVars_ = inpVars.size();

    postVars_ = postVars;
    preVars_ = preVars;
    inpVars_ = inpVars;

     /* the permutation array */
    size_t n=ddmgr_->ReadSize();
    permute_ = new int[n];
    for(size_t i=0; i<n; i++)
      permute_[i]=i;
    for(size_t i=0; i<nssVars_; i++)
      permute_[preVars_[i]]=postVars_[i];
    /* create a cube with the input Vars */
    BDD* vars = new BDD[nisVars_];
    for (size_t i=0; i<nisVars_; i++)
      vars[i]=ddmgr_->bddVar(inpVars_[i]);
    cubeInput_ = ddmgr_->bddComputeCube(vars,NULL,nisVars_);
    delete[] vars;
    /* create a cube with the post Vars */
    vars = new BDD[nssVars_];
    for (size_t i=0; i<nssVars_; i++)
      vars[i]=ddmgr_->bddVar(postVars_[i]);
    cubePost_ = ddmgr_->bddComputeCube(vars,NULL,nssVars_);
    delete[] vars;

    /* copy the transition relation */
    R_ = transRelation;
    RR_ = R_.ExistAbstract(cubePost_);
  }

  ~ncsFixPoint() {
    delete[] permute_;
  }

  /* function: pre 
   *
   * computes the enforcable predecessor 
   *  
   * pre(Z) = { (x,u) | exists x': (x,u,x') in transitionRelation 
   *                    and (x,u,x') in transitionRelation  => x' in Z } 
   *
   */
  BDD pre(BDD Z)  {
    /* project onto state alphabet */
    Z=Z.ExistAbstract(cubePost_*cubeInput_);
    /* swap variables */
    Z=Z.Permute(permute_);
    /* find the (state, inputs) pairs with a post outside the safe set */
    BDD nZ = !Z;
    BDD F = R_.AndAbstract(nZ,cubePost_); 
    /* the remaining (state, input) pairs make up the pre */
    BDD nF = !F;
    BDD preZ= RR_.AndAbstract(nF,cubePost_);
    return preZ;
  }

  
  /* function: reach 
   *  
   * computation of the minimal fixed point mu Z.pre(Z) | T
   *
   */
  BDD reach(const BDD &T, int verbose=0)  {

	  std::cout << "Started to compute the reach controller ! " << std::endl;
	  std::cout << "The target set is ";
	  PrintBDD("T", T);

    /* check if target is a subset of the state space */
    std::vector<unsigned int> sup = T.SupportIndices();
    for(size_t i=0; i<sup.size();i++) {
      int marker=0;
      for(size_t j=0; j<nssVars_; j++) {
        if (sup[i]==preVars_[j])
          marker=1;
      }
      if(!marker) {
        std::ostringstream os;
        os << "Error: reach: the target set depends on variables outside of the state space.";
        throw std::invalid_argument(os.str().c_str());
      }
    }
    if(verbose) 
      std::cout << "Iterations: ";

    BDD Z = ddmgr_->bddOne();
    BDD ZZ = ddmgr_->bddZero();


    /* the controller */
    BDD C = ddmgr_->bddZero();

    std::cout << std::endl;
    /* as long as not converged */
    size_t i;
    for(i=1; ZZ != Z; i++ ) {

    	std::cout << i << "]:" << std::endl;

      Z=ZZ;

      BDD Pre_Z = pre(Z);
      PrintBDD("PreZ", Pre_Z);
      ZZ= Pre_Z | T;
      PrintBDD("PreZ_or_T", ZZ);

      /* new (state/input) pairs */
      BDD N = ZZ & (!(C.ExistAbstract(cubeInput_)));
      PrintBDD("N", N);

      /* add new (state/input) pairs to the controller */
      C=C | N;
      PrintBDD("C", C);

      /* print progress */
      if(verbose) {
        std::cout << ".";
        std::flush(std::cout);
        if(!(i%80))
          std::cout << std::endl;
      }
    }
    if(verbose) 
      std::cout << " number: " << i << std::endl;
    return C;
  }


}; /* close class def */

#endif /* ncsFixPoint_HH_ */
