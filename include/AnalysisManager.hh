#ifndef ANALYSISMANAGER_HH
#define ANALYSISMANAGER_HH

#include "globals.hh"

class AnalysisManager
{
public:
  AnalysisManager();
  ~AnalysisManager();

  void Book();
  void Finish();
  void SetNeutronDetectors(G4int neutronDetectorCount);
  void ConfigureDetectorOutput(G4bool dssd, G4bool ic, G4bool siMon);

  void Detector_DSSD(G4double edep, G4double energy, G4double A, G4double Z,
                     G4double theta, G4double phi, G4double posx, G4double posy);
  void Detector_IC(G4double edep, G4double energy, G4double A, G4double Z);
  void Detector_SiMon(G4double edep, G4double energy, G4double A, G4double Z);

  void TargetTrack(G4int eventID, G4int trackID,
                   G4double ekin, G4double theta, G4double phi);

  void LSHitTrack(G4int eventID, G4int detID,
                  G4int trackID, G4int parentID,
                  G4double edep, G4int Z, G4int A);

  void LSNeutronCrossing(G4int eventID, G4int detID,
                         G4int trackID, G4int parentID,
                         G4double ekin, G4double tof,
                         G4double posx, G4double posy, G4double posz,
                         G4double theta, G4double phi);

private:
  G4bool factoryOn = false;
  G4int LScinNum = 0;
  G4bool writeDSSD = false;
  G4bool writeIC = false;
  G4bool writeSiMon = false;
  G4int dssdNtupleId = -1;
  G4int icNtupleId = -1;
  G4int siMonNtupleId = -1;
};

#endif
