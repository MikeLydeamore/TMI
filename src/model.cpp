#include <Rcpp.h>
#include "MarkovChainSimulator/MarkovChain/MarkovChain.hpp"
#include "MarkovChainSimulator/MarkovChain/MarkovChain.cpp"
#include "MarkovChainSimulator/MarkovChain/Serialiser.hpp"
#include "MarkovChainSimulator/MarkovChain/Serialiser.cpp"
using namespace Rcpp;

class SerialiserFullR : public Serialiser {
  
private:
  std::map<std::string, std::vector<double>> mResults;

public:
  
  SerialiserFullR() {}
  
 
  
  virtual void serialise(double t, state_values states)
  {
    
    mResults["t"].push_back(t);
    for (auto& state : states)
    {
      mResults[state.first].push_back(state.second);
    }
  }
  
  virtual void serialiseHeader(state_values states)
  {
  }
  
  virtual void serialiseFinally(double t, state_values states)
  {
    serialise(t, states);
  }
  
  List getResults()
  {
    List list(mResults.size());
    CharacterVector namevec;
    int i = 0;
    for (auto& result : mResults)
    {
      namevec.push_back(result.first);
      list[i] = result.second;
      i++;
    }
    list.attr("names") = namevec;
    return list;
  }
  
};

std::map<std::string, double> convertListToMap(List list)
{
  std::map<std::string, double> ret;
  CharacterVector names = list.names();
  for (int i = 0 ; i < names.size() ; i++)
  {
    ret[as<std::string>(names[i])] = as<double>(list[i]);
  }
  
  return (ret);
}

class LinearisedSISModel
{
private:
  double mLambda;
  double mGamma;
  
public:
  LinearisedSISModel(double lambda, double gamma) : mLambda(lambda), mGamma(gamma) {}
  
  void setupModel(MarkovChain &rChain)
  {
    rChain.addState("S", 1);
    rChain.addState("I", 0);
    
    rChain.addTransition(new TransitionIndividual("S", "I", mLambda));
    rChain.addTransition(new TransitionIndividual("I", "S", mGamma));
  }
};


// [[Rcpp::export]]
List runLinearisedSISModel(double lambda, double gamma, double max_time) {
  SerialiserFullR serialiser;
  
  LinearisedSISModel model(lambda, gamma);
  
  MarkovChain chain;
  model.setupModel(chain);
  
  chain.setSerialiser(&serialiser);
  chain.setMaxTime(max_time);
  chain.solve(MarkovChain::SOLVER_TYPE_GILLESPIE);
  
  return (serialiser.getResults());
}


