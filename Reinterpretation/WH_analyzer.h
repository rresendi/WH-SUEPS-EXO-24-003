#ifndef analysis_user_h
#define analysis_user_h
#include <algorithm> // for std::sort
#include <fastjet/ClusterSequence.hh>
#include <fastjet/PseudoJet.hh>
#include "SampleAnalyzer/Process/Analyzer/AnalyzerBase.h"
#include "SampleAnalyzer/Interfaces/root/RootMainHeaders.h"

namespace MA5
{
class user : public AnalyzerBase
{
  INIT_ANALYSIS(user, "user")

  public : 
      MAbool Initialize(const MA5::Configuration& cfg, const std::map<std::string,std::string>& parameters);
      void Finalize(const SampleFormat& summary, const std::vector<SampleFormat>& files);
      MAbool Execute(SampleFormat& sample, const EventFormat& event);
      void FillHistoAllRegions(const std::string &baseName, double value, const std::vector<std::string> &regions);

  private : 
};
}

#endif