/*
 * NTFTS.hh
 *
 *  created on: 09.10.2015
 *      author: rungger + m.khaled
 *
 *  modified version of the the symbolic model to accomodate for output
 *
 */

#ifndef NFTS_HH_
#define NFTS_HH_

#include <iostream>
#include <stdexcept>


#include "TicToc.hh"

#include "cuddObj.hh"
#include "cudd.h"
#include "dddmp.h"


#include "SymbolicSet.hh"
#include "SymbolicModel.hh"

#include "BDDTricks.hh"


namespace scots {
/*
 * class: NFTS
 *
 * Derived class from <SymbolicModel>
 * 
 * Constructs an output-based symbolic model according to
 *
 */

template<class stateType_, class inputType_, class outType_>
class NFTS: public SymbolicModel {

protected:

  /* var: outSpace_ */
  SymbolicSet *outSpace_;

  /* var: nosVars_ 
   * number of output space bdd variables */
  size_t nosVars_;

  /* var: outVars_ 
   * array of indices of the output space bdd variables  */
  size_t* outVars_;

  /* var: outMap_ 
   * the bdd representation of the state/output relation  X x O
   * the bdd variables of the map are given by
   *  preVars_ x inpVars_ x outVars_  */
  BDD outMap_;

public:

  /* constructor: NFTS 
   *
   * Representation of the transition relation and output-map as BDDs 
   *   transitions:
   *   preX  x U x postX
   *
   *   output-map:
   *   preX  x U x output
   *  
   * provide SymbolicSet for preX 
   * provide SymbolicSet for U
   * provide SymbolicSet for postX 
   * provide SymbolicSet for output 
   *
   * the SymbolicSet for preX and postX need to be identical, except the
   * BDD variable IDs need to differ
   * 
   */
  NFTS(SymbolicSet *stateSpace, 
       SymbolicSet *inputSpace, 
       SymbolicSet *stateSpacePost,
       SymbolicSet *outSpace): 
       SymbolicModel(stateSpace,inputSpace,stateSpacePost), outSpace_(outSpace){
    
    /* prrepare the nosVars */
    nosVars_=0;
    for(size_t i=0; i<outSpace_->getDimension(); i++) 
      for(size_t j=0; j<outSpace_->getNofBddVars()[i]; j++) 
       nosVars_++;

    /* prrepare the outVars_ */
    outVars_ = new size_t[nosVars_];
    for(size_t k=0, i=0; i<outSpace_->getDimension(); i++) {
      for(size_t j=0; j<outSpace_->getNofBddVars()[i]; k++, j++) {
       outVars_[k]=outSpace_->getIndBddVars()[i][j];
      }
    }

    /* initialize the transition relation */
    outMap_=stateSpace_->getSymbolicSet()*inputSpace_->getSymbolicSet();

  }


  ~NFTS(void){
    delete[] outVars_;
  }

  
  Cudd* getCuddManager(){
    return ddmgr_;
  }


  /* function:  getNumOutVars 
   * get the number of output bddVars */
  inline size_t getNumOutVars(void) {
    return nosVars_;
  };

  /* function:  getNumStateVars
     * get the number of state bddVars */
  inline size_t getNumStateVars(void){
	  return stateSpace_->getNVars();
  }

  /* function:  getNumInputVars
       * get the number of state bddVars */
  inline size_t getNumInputVars(void){
	  return inputSpace_->getNVars();
  }

  /* function:  getOutVars 
   * get the array of indicies of bddVars of outputs */
  inline const size_t* getOutVars(void) const {
    return outVars_;
  };

  /* function:  getTransSize 
   * get the number of elements in the output map */
  inline double getInputSize(void) {
    return getInputSpace()->getSize();
  };

  /* function:  getTransSize 
   * get the number of elements in the output map */
  inline double getStateSize(void) {
    return getStateSpace()->getSize();
  };

  /* function:  getOutSize 
   * get the number of elements in the output */
  inline double getOutSize(void) {
    return getOutSpace()->getSize();
  };

  /* function:  getOutMapSize 
   * get the number of elements in the output map */
  inline double getOutMapSize(void) {
    return outMap_.CountMinterm(nssVars_+ nisVars_ + nosVars_);
  };

  /* function:  getTransSize 
   * get the number of elements in the output map */
  inline double getTransSize(void) {
    return getSize();
  };

  /* function:  getStateSpace */
  inline const SymbolicSet* getStateSpace(void) {
      return stateSpace_;
  };

  /* function:  getInputSpace */
  inline const SymbolicSet* getInputSpace(void) {
      return inputSpace_;
  };

  /* function:  getOutSpace */
  inline const SymbolicSet* getOutSpace(void) {
      return outSpace_;
  };

  /* function:  getPostStatSpace */
  inline const SymbolicSet* getPostStatSpace(void) {
      return stateSpacePost_;
  };


  /* function:  getOutputMap 
   * get the SymbolicSet which represents output-map in X x U x O */
  inline SymbolicSet getOutputMap(void) const {

    SymbolicSet si(*stateSpace_,*inputSpace_);
    SymbolicSet outmap(si,*outSpace_);
    /* fill SymbolicSet with transition relation */
    outmap.setSymbolicSet(outMap_);
    return outmap;
  }; 


  /* function:  buildNFTS
   *
   * - provide the solution of the system at sampling time and
   * - provide the solution of the linear system associated with the growth bound 
   * at sampling time 
   * - assigns the states to their corresponding output points
   *
   * see the example directory for the specific format
   *
   */
  template<class F1, class F2, class F3>
  void buildNFTS(F1 &system_post, F2 &radius_post, F3 &out_of_state, int verbosity = 1) {

    /* create the BDD's with numbers 0,1,2,.., #gridPoints */
    size_t dim=stateSpace_->getDimension();
    const size_t* nvars= stateSpace_->getNofBddVars();
    BDD **bddVars = new BDD*[dim];
    for(size_t n=0, i=0; i<dim; i++) {
      bddVars[i]= new BDD[nvars[i]];
      for(size_t j=0; j<nvars[i]; j++)  {
        bddVars[i][j]=ddmgr_->bddVar(postVars_[n+j]);
      }
      n+=nvars[i];
    }
    const size_t* ngp= stateSpace_->getNofGridPoints();
    BDD **num = new BDD*[dim];
    for(size_t i=0; i<dim; i++) {
      num[i] = new BDD[ngp[i]];
      int *phase = new int[nvars[i]];
      for(size_t j=0;j<nvars[i];j++)
        phase[j]=0;
      for(size_t j=0;j<ngp[i];j++) {
        int *p=phase;
        int x=j;
        for (; x; x/=2) *(p++)=0+x%2;
        num[i][j]= ddmgr_->bddComputeCube(bddVars[i],(int*)phase,nvars[i]);
      }
      delete[] phase;
      delete[] bddVars[i];
    }
    delete[] bddVars;

    /* bdd nodes in pre and input variables */
    DdManager *mgr = ddmgr_->getManager();
    size_t ndom=nssVars_+nisVars_;
    int*  phase = new int[ndom];
    DdNode**  dvars = new DdNode*[ndom];
    for(size_t i=0;i<nssVars_; i++)
      dvars[i]=Cudd_bddIthVar(mgr,preVars_[i]);
    for(size_t i=0;i<nisVars_; i++)
      dvars[nssVars_+i]=Cudd_bddIthVar(mgr,inpVars_[i]);

    /* initialize cell radius
     * used to compute the growth bound */
    stateType_ eta;
    stateSpace_->copyEta(&eta[0]);
    stateType_ z;
    stateSpace_->copyZ(&z[0]);
    stateType_ r;

    stateType_ first;
    stateSpace_->copyFirstGridPoint(&first[0]);

    transitionRelation_=ddmgr_->bddZero();
    outMap_ = ddmgr_->bddZero();
    const int* minterm;

    /* compute constraint set against the post is checked */
    size_t n=ddmgr_->ReadSize();
    int* permute = new int[n];
    for(size_t i=0; i<nssVars_; i++)
      permute[preVars_[i]]=postVars_[i];
    BDD ss = stateSpace_->getSymbolicSet();
    BDD constraints=ss.Permute(permute);
    delete[] permute;

    int* phase_y = new int[nosVars_];
    BDD* dvars_y = new BDD[nosVars_];
    for(size_t i=0; i<nosVars_; i++)
	dvars_y[i] = ddmgr_->bddVar(outVars_[i]);

    /** big loop over all state elements and input elements **/
    for(begin(); !done(); next()) {
      if(verbosity > 0 )
    	  progress();

      minterm=currentMinterm();

      /* current state */
      stateType_ x;
      stateSpace_->mintermToElement(minterm,&x[0]);
      /* current input */
      inputType_ u;
      inputSpace_->mintermToElement(minterm,&u[0]);
      /* cell radius (including measurement errors) */
      for(size_t i=0; i<dim; i++)
        r[i]=eta[i]/2.0+z[i];

      /* integrate system and radius growth bound */
      /* the result is stored in x and r */
      std::vector<stateType_> x_posts = system_post(x,u);
      radius_post(r,u);

      outType_ y;
      out_of_state(y,x,u);

      for(size_t k=0; k<x_posts.size(); k++){
		  stateType_ x_post = x_posts[k];
		  /* determine the cells which intersect with the attainable set*/
		  /* start with the computation of the indices */
		  BDD post=ddmgr_->bddOne();
		  for(size_t i=0; i<dim; i++) {
			int lb = std::lround(((x_post[i]-r[i]-z[i]-first[i])/eta[i]));
			int ub = std::lround(((x_post[i]+r[i]+z[i]-first[i])/eta[i]));
			if(lb<0 || ub>=(int)ngp[i]) {
			  post=ddmgr_->bddZero();
			  break;
			}
			BDD zz=ddmgr_->bddZero();
			for(int j=lb; j<=ub; j++) {
			  zz|=num[i][j];
			}
			post &= zz;
		  }
		  if(!(post==ddmgr_->bddZero()) && post<= constraints) {
			/* compute bdd for the current x and u element and add x' */
			for(size_t i=0;i<nssVars_; i++)
			  phase[i]=minterm[preVars_[i]];
			for(size_t i=0;i<nisVars_; i++)
			  phase[nssVars_+i]=minterm[inpVars_[i]];
			BDD current(*ddmgr_,Cudd_bddComputeCube(mgr,dvars,phase,ndom));

			outSpace_->elementToMinterm(y.data(), phase_y);
			BDD cube_y = ddmgr_->bddComputeCube(dvars_y,phase_y,nosVars_);

			BDD current_for_y = current;
	
			outMap_ += current_for_y&cube_y;
	
			current&=post;
			transitionRelation_ +=current;
		  }
      }
    }

    for(size_t i=0; i<dim; i++) 
      delete[] num[i];

    delete[] num;
    delete[] dvars;
    delete[] dvars_y;
    delete[] phase;
    delete[] phase_y;
  }

  /* function:  getPost
   *
   * simulates the NFTS by BDD operations to compute the post of a sequence of inputs
   * while ensuring that the observed outputs reduces non-deterministity
   *
   */
  BDD getPost(BDD x0, std::vector<BDD> appliedInputs, std::vector<BDD> observedOutputs) {

	BDD H = getOutputMap().getSymbolicSet();			// the out map
	BDD T = getTransitionRelation().getSymbolicSet();	// the transition relation

	if (H * x0 * observedOutputs[0] == ddmgr_->bddZero()) {
		std::cout
				<< "NFTS::getPost (Warning): initial set has non-correct output !"
				<< std::endl;
		return ddmgr_->bddZero();
	}

	if (observedOutputs.size() != (appliedInputs.size()+1)) {
		std::cout
				<< "NFTS::getPost (error): mismatch in sizes btw input/output arrays!"
				<< std::endl;
		return ddmgr_->bddZero();
	}

	// getting state/post bddVars ready
	std::vector<size_t> stateBddVars, postBddVars;
	for (size_t i = 0; i < getStateSpace()->getDimension(); i++)
		for (size_t j = 0; j < getStateSpace()->getNofBddVars()[i];
				j++)
			stateBddVars.push_back(
					getStateSpace()->getIndBddVars()[i][j]);

	for (size_t i = 0; i < getPostStatSpace()->getDimension(); i++)
		for (size_t j = 0; j < getPostStatSpace()->getNofBddVars()[i];
				j++)
			postBddVars.push_back(
					getPostStatSpace()->getIndBddVars()[i][j]);

	BDD post;
	BDD pre = x0;
	for (size_t i = 0; i < appliedInputs.size(); i++) {

		// finiding post states
		post = T * pre * appliedInputs[i];
		post = BddTricks::ProjectBDD(*ddmgr_, post, postBddVars);// project as pre-states

		pre = BddTricks::PermutreBdd(*ddmgr_, post, postBddVars,
				stateBddVars);		// permute to states

		pre = pre * H * observedOutputs[i+1];
		pre = BddTricks::ProjectBDD(*ddmgr_, pre, stateBddVars);// project as pre-states

	}

	return pre;
  }  
}; /* close class def */
} /* close namespace */

#endif /* SYMBOLICMODELGROWTHBOUND_HH_ */
