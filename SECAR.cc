//
// ********************************************************************
// * License and Disclaimer                                           *
// *                                                                  *
// * The  Geant4 software  is  copyright of the Copyright Holders  of *
// * the Geant4 Collaboration.  It is provided  under  the terms  and *
// * conditions of the Geant4 Software License,  included in the file *
// * LICENSE and available at  http://cern.ch/geant4/license .  These *
// * include a list of copyright holders.                             *
// *                                                                  *
// * Neither the authors of this software system, nor their employing *
// * institutes,nor the agencies providing financial support for this *
// * work  make  any representation or  warranty, express or implied, *
// * regarding  this  software system or assume any liability for its *
// * use.  Please see the license in the file  LICENSE  and URL above *
// * for the full disclaimer and the limitation of liability.         *
// *                                                                  *
// * This  code  implementation is the result of  the  scientific and *
// * technical work of the GEANT4 collaboration.                      *
// * By using,  copying,  modifying or  distributing the software (or *
// * any work based  on the software)  you  agree  to acknowledge its *
// * use  in  resulting  scientific  publications,  and indicate your *
// * acceptance of all terms of the Geant4 Software license.          *
// ********************************************************************
//
// Author: Pelagia Tsintari, pelagia.tsin@gmail.com
// 
// Code based on the advanced example radioprotection

#include "DetectorConstruction.hh"
#include "AnalysisManager.hh"
#include "ActionInitialization.hh"

#include "G4RunManagerFactory.hh"

#include "G4UImanager.hh"
#include "G4VisExecutive.hh"
#include "G4UIExecutive.hh"

#include "G4Timer.hh"
#include "G4Types.hh"
#include "Randomize.hh"
#include "time.h"

#include "G4Scintillation.hh"
#include "G4OpBoundaryProcess.hh"
#include "G4EmStandardPhysics.hh"
#include "QGSP_BERT_HP.hh"
//#include "FTFP_BERT_HP.hh"
#include "G4PhysListFactory.hh"
#include "G4OpticalPhysics.hh"
#include "G4DecayPhysics.hh"
#include "G4RadioactiveDecayPhysics.hh"
#include "G4ParticleHPManager.hh"

int main(int argc, char** argv)
{
  
  //Instantiate the G4Timer object, to monitor the CPU time spent for the entire execution
  G4Timer* theTimer = new G4Timer();
  theTimer->Start();

  G4UIExecutive* ui = nullptr;
  if (argc == 1) {
    ui = new G4UIExecutive(argc, argv);
  }
  //Set the Random engine
  //The following guarantees random generation also for different runs in multithread
  CLHEP::RanluxEngine defaultEngine( 1234567, 4 );
  G4Random::setTheEngine( &defaultEngine );
  G4int seed = time( NULL );
  G4Random::setTheSeed( seed );

  G4RunManager* pRunManager = G4RunManagerFactory::CreateRunManager(G4RunManagerType::Default);
  // pRunManager->SetNumberOfThreads(1);

  // User action initialization
  AnalysisManager* analysisMan = new AnalysisManager();

  DetectorConstruction* detector = new DetectorConstruction(analysisMan);  
  pRunManager -> SetUserInitialization(detector);

  // Use NRESP model 
  G4ParticleHPManager::GetInstance()->SetUseNRESP71Model( true );
  G4PhysListFactory physicsFactory;
  G4VModularPhysicsList* physics = physicsFactory.GetReferencePhysList("QGSP_BERT_HP"); // icludes both Neutron and Particle HP
  G4OpticalPhysics* opticalPhysics = new G4OpticalPhysics();
  physics->RegisterPhysics(opticalPhysics);
  pRunManager -> SetUserInitialization(physics); 
   
  ActionInitialization* actions = new ActionInitialization(analysisMan, detector);
  pRunManager->SetUserInitialization(actions);

  // Initialize visualization with the default graphics system
  auto visManager = new G4VisExecutive(argc, argv);
  visManager->Initialize();
  // Get the pointer to the User Interface manager
  auto UImanager = G4UImanager::GetUIpointer();

  if (!ui) {
    // batch mode
    G4String command = "/control/execute";
    G4String fileName = argv[1];
    UImanager->ApplyCommand(command + " " + fileName);
  }
  else {
    // interactive mode
    UImanager->ApplyCommand("/control/execute macros/vis.mac");
    ui->SessionStart();
    delete ui;
  }    
  
  //Stop the time benchmark here
  theTimer->Stop();
  G4cout << "The simulation took: " << theTimer->GetRealElapsed()/60 << " m to run (real time)"<< G4endl;

  
  delete visManager; 
  delete pRunManager;
  return 0;
}
