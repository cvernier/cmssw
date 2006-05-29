/**
   \file
   Implementation of class ScheduleValidator

   \author Stefano ARGIRO
   \version $Id: ScheduleValidator.cc,v 1.9 2006/03/02 22:58:11 paterno Exp $
   \date 10 Jun 2005
*/

static const char CVSId[] = "$Id: ScheduleValidator.cc,v 1.9 2006/03/02 22:58:11 paterno Exp $";

#include "FWCore/ParameterSet/src/ScheduleValidator.h"
#include "FWCore/Utilities/interface/EDMException.h"

#include <sstream>
#include <iterator>

using namespace edm;
using namespace edm::pset;
using namespace std;

 
ScheduleValidator::ScheduleValidator(const ScheduleValidator::PathContainer& 
				     pathFragments,
				     const ParameterSet& processPSet): 
  nodes_(pathFragments),
  processPSet_(processPSet){
 
  vector<string> paths = processPSet.getParameter<vector<string> >("@paths");

  for(vector<string>::const_iterator pathIt= paths.begin();
      pathIt!=paths.end(); 
      ++pathIt){
    NodePtr head = findPathHead(*pathIt);
    gatherLeafNodes(head);
  }//for

}


NodePtr ScheduleValidator::findPathHead(string pathName){
  // cout << "in findPathHead" << endl;
  for (PathContainer::iterator pathIt= nodes_.begin();
       pathIt!=nodes_.end();++pathIt){
//cout << "  looking at " << (*pathIt)->type() << " " << (*pathIt)->name << endl;
    if ((*pathIt)->type()!="path" &&
        (*pathIt)->type()!="endpath") continue;
    if ((*pathIt)->name==pathName) return ((*pathIt)->wrapped_);

  }// for
  //cout << "did not find " << pathName << endl;
  // This can only cause a problem and should probably throw
  // an exception - jbk
  NodePtr ret;
  return ret;
}

void ScheduleValidator::gatherLeafNodes(NodePtr& basenode){

//cout << "gatherLeafNodes " << basenode.get() << endl;
//basenode->print(cout);

  if (basenode->type()=="," || basenode->type() =="&"){
    OperatorNode* onode = dynamic_cast<OperatorNode*>(basenode.get());
    gatherLeafNodes(onode->left_);
    gatherLeafNodes(onode->right_);

  } else {
    leaves_.push_back(basenode);
  }

}// gatherLeafNodes




void ScheduleValidator::validate(){

  dependencies_.clear();

  // iterate on leaf nodes
  std::list<NodePtr>::iterator leafIt;
  for (leafIt=leaves_.begin(); leafIt!=leaves_.end(); ++leafIt){
    
    DependencyList dep;

    // follow the tree up thru parent nodes  
    Node* p  = (*leafIt)->getParent();
    Node* son= (*leafIt).get();
    
    // make sure we don't redescend right of our parent
    
    while  (p){
      // if we got operator '&' continue
      if (p->type()=="&") {
	son =p;
	p= p->getParent(); 
	continue;
      }

      OperatorNode* node = dynamic_cast<OperatorNode*>(p);
      // make sure we don't redescend where we came from
      if (son->name != node->left_->name) findDeps(node->left_, dep);
      
      son=p;
      p= p->getParent();
    } // while
        
    dep.sort();   
    dep.unique(); // removes duplicates

    // insert the list of deps
    Dependencies::iterator depIt= dependencies_.find((*leafIt)->name);
    if (depIt!= dependencies_.end()){
      DependencyList& old_deplist = (*depIt).second;
    
      // if the list is different from an existing one
      // then we have an inconsitency
      if (old_deplist !=  dep) {

	ostringstream olddepstr,newdepstr;
	copy(old_deplist.begin(), old_deplist.end(), 
	      ostream_iterator<string>(olddepstr,","));
	copy(dep.begin(), dep.end(), 
	      ostream_iterator<string>(newdepstr,","));

	throw edm::Exception(errors::Configuration,"InconsistentSchedule")
	  << "Inconsistent schedule for module "
	  << (*leafIt)->name
	  << "\n"
	  << "Depends on " << olddepstr.str() 
	  << " but also on " << newdepstr.str()
	  << "\n";
      }
    }
    else {
      dependencies_[(*leafIt)->name] = dep;
    }

  }//for leaf
     
}// validate

void ScheduleValidator::findDeps(NodePtr& node, DependencyList& dep){

  // if we have an operand, add it to the list of dependencies
  if (node->type() =="operand") dep.push_back(node->name);
  // else follow the tree, unless the leaf is contained in the node
  else{

    OperatorNode* opnode = dynamic_cast<OperatorNode*>(node.get());
    
    findDeps(opnode->left_,dep);
    findDeps(opnode->right_,dep);
	   
  }
}// findDeps


std::string 
ScheduleValidator::dependencies(const std::string& modulename) const{


  Dependencies::const_iterator depIt = dependencies_.find(modulename);
  if (depIt== dependencies_.end()){
    ostringstream err;
    throw edm::Exception(errors::Configuration,"ScheduleDependencies")
      << "Error : Dependecies for "
      << modulename << " were not calculated";
  }

  ostringstream deplist;
  copy((*depIt).second.begin(), (*depIt).second.end(), 
	      ostream_iterator<string>(deplist,","));
  return deplist.str();

}
