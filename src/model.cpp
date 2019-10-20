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


class IndividualSISModel
{
private:
  double mBeta;
  double mGamma;
  int mNumIndividuals;
  int mInitialInfected;
  
public:
  IndividualSISModel(double beta, double gamma, int numIndividuals, int initialInfected) : 
  mBeta(beta), mGamma(gamma), mNumIndividuals(numIndividuals), mInitialInfected(initialInfected) {}
  
  void setupModel(MarkovChain &rChain)
  {
    std::cout << "Setting up" << std::endl;
    std::vector<std::string> infected_states(mNumIndividuals, "");
    for (int i = 1 ; i <= mNumIndividuals; i++)
    {
      std::string indiv_string = std::to_string(i);
      if (i <= mInitialInfected)
      {
        rChain.addState("S"+indiv_string, 0);
        rChain.addState("I"+indiv_string, 1);
      }
      else
      {
        rChain.addState("S"+indiv_string, 1);
        rChain.addState("I"+indiv_string, 0);
      }
      
      infected_states[i-1] = "I"+indiv_string;
    }
    
    //Transitions!
    for (int i = 1; i <= mNumIndividuals; i++)
    {
      std::string indiv_string = std::to_string(i);
      rChain.addTransition(new TransitionMassAction("S"+indiv_string, "I"+indiv_string, mBeta/(mNumIndividuals-1), infected_states));
      rChain.addTransition(new TransitionIndividual("I"+indiv_string, "S"+indiv_string, mGamma));
    }
  }
};

// [[Rcpp::export]]
List runIndividualSISModel(double beta, double gamma, int num_individuals, int num_infected, double max_time)
{
  SerialiserFullR serialiser;
  
  IndividualSISModel model(beta, gamma, num_individuals, num_infected);
  
  MarkovChain chain;
  model.setupModel(chain);
  
  chain.setSerialiser(&serialiser);
  chain.setMaxTime(max_time);
  chain.solve(MarkovChain::SOLVER_TYPE_GILLESPIE);
  
  return (serialiser.getResults());
}

class IndividualSIIModel
{
private:
  double mBeta;
  double mGamma;
  int mNumIndividuals;
  int mInitialInfected;
  
public:
  IndividualSIIModel(double beta, double gamma, int numIndividuals, int initialInfected) : 
  mBeta(beta), mGamma(gamma), mNumIndividuals(numIndividuals), mInitialInfected(initialInfected) {}
  
  void setupModel(MarkovChain &rChain)
  {
    std::cout << "Setting up" << std::endl;
    std::vector<std::string> infected_states(2*mNumIndividuals, "");
    for (int i = 1 ; i <= mNumIndividuals; i++)
    {
      std::string indiv_string = std::to_string(i);
      if (i <= mInitialInfected)
      {
        rChain.addState("S"+indiv_string, 0);
        rChain.addState("I"+indiv_string, 1);
        rChain.addState("U"+indiv_string, 1);
      }
      else
      {
        rChain.addState("S"+indiv_string, 1);
        rChain.addState("I"+indiv_string, 0);
        rChain.addState("U"+indiv_string, 0);
      }
      
      infected_states[i-1] = "I"+indiv_string;
      infected_states[mNumIndividuals+i-1] = "U"+indiv_string;
    }
    
    //Transitions!
    for (int i = 1; i <= mNumIndividuals; i++)
    {
      std::string indiv_string = std::to_string(i);
      rChain.addTransition(new TransitionMassAction("S"+indiv_string, "I"+indiv_string, mBeta/(mNumIndividuals-1), infected_states));
      rChain.addTransition(new TransitionIndividual("I"+indiv_string, "U"+indiv_string, 2*mGamma));
      rChain.addTransition(new TransitionIndividual("U"+indiv_string, "S"+indiv_string, 2*mGamma));
    }
  }
};

// [[Rcpp::export]]
List runIndividualSIIModel(double beta, double gamma, int num_individuals, int num_infected, double max_time)
{
  SerialiserFullR serialiser;
  
  IndividualSIIModel model(beta, gamma, num_individuals, num_infected);
  
  MarkovChain chain;
  model.setupModel(chain);
  
  chain.setSerialiser(&serialiser);
  chain.setMaxTime(max_time);
  chain.solve(MarkovChain::SOLVER_TYPE_GILLESPIE);
  
  return (serialiser.getResults());
}
