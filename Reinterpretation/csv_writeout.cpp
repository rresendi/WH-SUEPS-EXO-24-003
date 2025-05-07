  //////////////////////////////////////////////
  //           Writing out the CSV            //
  //////////////////////////////////////////////

  std::ofstream outFile("kinematics_output.csv", std::ios::app);

  // Write the header one time

  if (outFile.tellp() == 0) { 
    outFile << "EventNumber,ak15NTracks,wTransverseMass,wPt,ak4Pt,boostedSphericity,metPt,looseLep,tightLep,NJets,ak15Pt,wPhi\n"; 
  }

  // Event information to be stored
  // Event number
  static int eventCounter = 0;
  eventCounter++;

  outFile
      << eventCounter << ","
      << Ak15jetConstituents.at(0).size() << ","                  // size_t or unassigned long
      << wTransverseMass << ","              // double
      << recoW.pt() << ","                          // assumed float
      << Ak4jets.at(0).pt() << ","                        // assumed float
      << boostedSphericity << ","            // double
      << smearedMET.pt() << ","                        // double
      << loose_leptons.size() << ","                     // int
      << tight_leptons.size() << ","                     // int
      << Ak4jets.size() << ","                        // int
      << Ak15Jets.at(0).pt() << ","                       // float
      << recoW.phi() << "\n";                       // assumed float

  outFile.close();
