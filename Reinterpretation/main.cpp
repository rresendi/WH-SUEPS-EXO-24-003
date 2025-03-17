// SampleHeader header
#include "SampleAnalyzer/Process/Core/SampleAnalyzer.h"
#include "SampleAnalyzer/User/Analyzer/analysisList.h"
using namespace MA5;

// -----------------------------------------------------------------------
// Info function
// -----------------------------------------------------------------------
int Info(SampleAnalyzer& manager)
{
  INFO << "BEGIN " << __FILE__ << endmsg;
  manager.AnalyzerList().Print();
  INFO << "END " << __FILE__ << endmsg;
  return 0;
}
// -----------------------------------------------------------------------
// main program
// -----------------------------------------------------------------------
int main(int argc, char *argv[])
{
  // Creating a manager
  SampleAnalyzer manager;
  BuildUserTable(manager.AnalyzerList());

  // Identifying --info argument
  if (argc==2)
  {
    std::string arg=argv[1];
    if (arg=="--info") return Info(manager);
  }

  // ---------------------------------------------------
  //                    INITIALIZATION
  // ---------------------------------------------------
  INFO << "    * Initializing all components" << endmsg;

  // Initializing the manager
  if (!manager.Initialize(argc,argv,"pdg.ma5")) return 1;

  // Creating data format for storing data
  EventFormat myEvent;
  std::vector<SampleFormat> mySamples;

  // Getting pointer to the analyzer
  std::map<std::string, std::string> parametersA1;
  AnalyzerBase* analyzer1 = 
      manager.InitializeAnalyzer("MadAnalysis5job","MadAnalysis5job.saf",parametersA1);
  if (analyzer1==0) return 1;

  // Post initialization (creates the new output directory structure)
  if(!manager.PostInitialize()) return 1;

  // Initializing PhysicsService for MC
  PHYSICS->mcConfig().Reset();
  // definition of the multiparticle "hadronic"
  manager.AddDefaultHadronic();
  // definition of the multiparticle "invisible"
  manager.AddDefaultInvisible();

  // ---------------------------------------------------
  //                      EXECUTION
  // ---------------------------------------------------
  INFO << "    * Running over files ..." << endmsg;

  // Loop over files
  while(1)
  {
    // Opening input file
    mySamples.push_back(SampleFormat());
    SampleFormat& mySample=mySamples.back();
    StatusCode::Type result1 = manager.NextFile(mySample);
    if (result1!=StatusCode::KEEP)
    {
      if (result1==StatusCode::SKIP) continue;
      else if (result1==StatusCode::FAILURE) {mySamples.pop_back(); break;}
    }
    
    int event_counter = 0;

    // Loop over events
    while (1)
    {
        std::cout << "Fetching next event (event " << event_counter << ")..." << std::endl;
        StatusCode::Type result2 = manager.NextEvent(mySample, myEvent);
    
        std::cout << "NextEvent returned: " << result2 << " at event " << event_counter << std::endl;
    
        if (result2 != StatusCode::KEEP)
        {
            if (result2 == StatusCode::SKIP) {
                std::cout << "Event skipped." << std::endl;
                continue;
            }
            else if (result2 == StatusCode::FAILURE) {
                std::cerr << "ERROR: Event processing failed. Breaking loop." << std::endl;
                break;
            }
        }
    
        std::cout << "Updating progress bar..." << std::endl;
        manager.UpdateProgressBar();
    
        std::cout << "Executing analyzer1..." << std::endl;
        if (!analyzer1->Execute(mySample, myEvent)) {
            std::cout << "Event did not pass analyzer1. Skipping to next event." << std::endl;
            continue;
        }
    
        std::cout << "Event processed successfully. Total events processed: " << ++event_counter << std::endl;
    
        // Stop after processing 1000 events
        if (event_counter >= 1000) {
            std::cout << "Reached event limit, stopping loop." << std::endl;
            break;
        }
    }    
  }

  // ---------------------------------------------------
  //                     FINALIZATION
  // ---------------------------------------------------
  INFO << "    * Finalizing all components ..." << endmsg;

  // Finalizing all components
  manager.Finalize(mySamples,myEvent);
  return 0;
}
