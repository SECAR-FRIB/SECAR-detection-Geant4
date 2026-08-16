#include "AnalysisManager.hh"

#include "G4AnalysisManager.hh"

AnalysisManager::AnalysisManager() : factoryOn(false) {}

AnalysisManager::~AnalysisManager()
{
  delete G4AnalysisManager::Instance();
}

void AnalysisManager::SetNeutronDetectors(G4int neutronDetectorCount)
{
  LScinNum = neutronDetectorCount;
}

void AnalysisManager::ConfigureDetectorOutput(G4bool dssd, G4bool ic, G4bool siMon)
{
  writeDSSD = dssd;
  writeIC = ic;
  writeSiMon = siMon;
}

void AnalysisManager::Book()
{
  auto* manager = G4AnalysisManager::Instance();
  manager->SetDefaultFileType("root");
  manager->SetVerboseLevel(1);
  manager->SetNtupleMerging(true);

  if (!manager->OpenFile("SECAR_sim_output.root")) {
    G4cout << "\n---> AnalysisManager::Book(): cannot open output file" << G4endl;
    return;
  }

  manager->SetFirstNtupleId(1);

  manager->CreateNtuple("TargetTrack", "Neutrons emitted at target");
  manager->CreateNtupleDColumn(1, "eventID");
  manager->CreateNtupleDColumn(1, "trackID");
  manager->CreateNtupleDColumn(1, "isNeutron");
  manager->CreateNtupleDColumn(1, "ekin_tar");
  manager->CreateNtupleDColumn(1, "theta_tar");
  manager->CreateNtupleDColumn(1, "phi_tar");
  manager->FinishNtuple(1);

  manager->CreateNtuple("LSHitTrack", "Detected liquid-scintillator hits");
  manager->CreateNtupleDColumn(2, "eventID");
  manager->CreateNtupleDColumn(2, "detID");
  manager->CreateNtupleDColumn(2, "trackID");
  manager->CreateNtupleDColumn(2, "parentID");
  manager->CreateNtupleDColumn(2, "edep_sum");
  manager->CreateNtupleDColumn(2, "Z_dep");
  manager->CreateNtupleDColumn(2, "A_dep");
  manager->FinishNtuple(2);

  manager->CreateNtuple("LSNeutronCrossing", "Neutrons entering liquid scintillators");
  manager->CreateNtupleDColumn(3, "eventID");
  manager->CreateNtupleDColumn(3, "detID");
  manager->CreateNtupleDColumn(3, "trackID");
  manager->CreateNtupleDColumn(3, "parentID");
  manager->CreateNtupleDColumn(3, "ekin_entry");
  manager->CreateNtupleDColumn(3, "ToF_entry");
  manager->CreateNtupleDColumn(3, "posX_entry");
  manager->CreateNtupleDColumn(3, "posY_entry");
  manager->CreateNtupleDColumn(3, "posZ_entry");
  manager->CreateNtupleDColumn(3, "theta_entry");
  manager->CreateNtupleDColumn(3, "phi_entry");
  manager->FinishNtuple(3);

  if (writeDSSD) {
    dssdNtupleId = manager->CreateNtuple("DSSD", "DSSD variables");
    for (const auto* name : {"totalEdep", "ekin", "A", "Z", "theta", "phi", "posx", "posy"})
      manager->CreateNtupleDColumn(dssdNtupleId, name);
    manager->FinishNtuple(dssdNtupleId);
  }

  if (writeIC) {
    icNtupleId = manager->CreateNtuple("IC", "IC variables");
    for (const auto* name : {"totalEdep", "ekin", "A", "Z"})
      manager->CreateNtupleDColumn(icNtupleId, name);
    manager->FinishNtuple(icNtupleId);
  }

  if (writeSiMon) {
    siMonNtupleId = manager->CreateNtuple("SiMon", "SiMon variables");
    for (const auto* name : {"totalEdep", "ekin", "A", "Z"})
      manager->CreateNtupleDColumn(siMonNtupleId, name);
    manager->FinishNtuple(siMonNtupleId);
  }

  factoryOn = true;
}

void AnalysisManager::Detector_DSSD(G4double edep, G4double energy,
                                    G4double A, G4double Z,
                                    G4double theta, G4double phi,
                                    G4double posx, G4double posy)
{
  if (!writeDSSD) return;
  auto* manager = G4AnalysisManager::Instance();
  const G4double values[] = {edep, energy, A, Z, theta, phi, posx, posy};
  for (G4int column = 0; column < 8; ++column)
    manager->FillNtupleDColumn(dssdNtupleId, column, values[column]);
  manager->AddNtupleRow(dssdNtupleId);
}

void AnalysisManager::Detector_IC(G4double edep, G4double energy,
                                  G4double A, G4double Z)
{
  if (!writeIC) return;
  auto* manager = G4AnalysisManager::Instance();
  const G4double values[] = {edep, energy, A, Z};
  for (G4int column = 0; column < 4; ++column)
    manager->FillNtupleDColumn(icNtupleId, column, values[column]);
  manager->AddNtupleRow(icNtupleId);
}

void AnalysisManager::Detector_SiMon(G4double edep, G4double energy,
                                     G4double A, G4double Z)
{
  if (!writeSiMon) return;
  auto* manager = G4AnalysisManager::Instance();
  const G4double values[] = {edep, energy, A, Z};
  for (G4int column = 0; column < 4; ++column)
    manager->FillNtupleDColumn(siMonNtupleId, column, values[column]);
  manager->AddNtupleRow(siMonNtupleId);
}

void AnalysisManager::LSNeutronCrossing(G4int eventID, G4int detID,
                                        G4int trackID, G4int parentID,
                                        G4double ekin, G4double tof,
                                        G4double posx, G4double posy, G4double posz,
                                        G4double theta, G4double phi)
{
  auto* manager = G4AnalysisManager::Instance();
  const G4double values[] = {G4double(eventID), G4double(detID), G4double(trackID),
                             G4double(parentID), ekin, tof, posx, posy, posz, theta, phi};
  for (G4int column = 0; column < 11; ++column)
    manager->FillNtupleDColumn(3, column, values[column]);
  manager->AddNtupleRow(3);
}

void AnalysisManager::TargetTrack(G4int eventID, G4int trackID,
                                  G4double ekin, G4double theta, G4double phi)
{
  auto* manager = G4AnalysisManager::Instance();
  manager->FillNtupleDColumn(1, 0, eventID);
  manager->FillNtupleDColumn(1, 1, trackID);
  manager->FillNtupleDColumn(1, 2, 1.0);
  manager->FillNtupleDColumn(1, 3, ekin);
  manager->FillNtupleDColumn(1, 4, theta);
  manager->FillNtupleDColumn(1, 5, phi);
  manager->AddNtupleRow(1);
}

void AnalysisManager::LSHitTrack(G4int eventID, G4int detID,
                                 G4int trackID, G4int parentID,
                                 G4double edep, G4int Z, G4int A)
{
  auto* manager = G4AnalysisManager::Instance();
  manager->FillNtupleDColumn(2, 0, eventID);
  manager->FillNtupleDColumn(2, 1, detID);
  manager->FillNtupleDColumn(2, 2, trackID);
  manager->FillNtupleDColumn(2, 3, parentID);
  manager->FillNtupleDColumn(2, 4, edep);
  manager->FillNtupleDColumn(2, 5, Z);
  manager->FillNtupleDColumn(2, 6, A);
  manager->AddNtupleRow(2);
}

void AnalysisManager::Finish()
{
  auto* manager = G4AnalysisManager::Instance();
  manager->Write();
  manager->CloseFile();
}
