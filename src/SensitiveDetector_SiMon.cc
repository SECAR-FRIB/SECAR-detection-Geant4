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

#include "SensitiveDetector_SiMon.hh"
#include "G4HCofThisEvent.hh"
#include "G4Step.hh"
#include "G4ThreeVector.hh"
#include "G4SDManager.hh"
#include "G4ios.hh"
#include "G4SystemOfUnits.hh"
#include "G4DynamicParticle.hh"
#include "G4INCLRandom.hh"

SensitiveDetector_SiMon::SensitiveDetector_SiMon(const G4String& name, const G4String& hitsCollectionName, AnalysisManager* analysis_manager) 
 : G4VSensitiveDetector(name),
   fHitsCollection(NULL)
{
  collectionName.insert(hitsCollectionName);
  analysis = analysis_manager;
  count = 0;
}

SensitiveDetector_SiMon::~SensitiveDetector_SiMon() 
{}

void SensitiveDetector_SiMon::Initialize(G4HCofThisEvent* hce)
{
  // Create hits collection
  fHitsCollection = new SensitiveDetector_SiMonHitsCollection(SensitiveDetectorName, collectionName[0]); 

  // Add this collection in hce
  G4int hcID = G4SDManager::GetSDMpointer()->GetCollectionID(collectionName[0]);
  hce->AddHitsCollection( hcID, fHitsCollection ); 
}

G4bool SensitiveDetector_SiMon::ProcessHits(G4Step* aStep, G4TouchableHistory*)
{  
  G4Track* aTrack = aStep->GetTrack();
  
  G4double edep = aStep->GetTotalEnergyDeposit();
  G4double ekin = aStep->GetPreStepPoint()->GetKineticEnergy();
  // G4double ekin_p = aStep->GetPostStepPoint()->GetKineticEnergy();
  G4int charge = aTrack->GetDefinition()->GetPDGCharge();  
  G4int baryon = aTrack->GetDefinition()->GetBaryonNumber();
  G4String particleName = aTrack->GetDefinition()->GetParticleName();

  if (ekin==0.) return false;
  G4String volumeName = aStep -> GetTrack() -> GetVolume() -> GetName();
  if(volumeName != "IC_gas_dE_phys") return false; 
  
  
  SensitiveDetectorHit_SiMon* newHit = new SensitiveDetectorHit_SiMon();

  newHit -> SetEdep(edep);
  newHit -> SetEkin(ekin);
  // newHit -> SetEkin_p(ekin_p);
  newHit -> SetParticleName(particleName);
  newHit -> SetParticleCharge(charge);
  newHit -> SetBaryonNum(baryon);

  fHitsCollection -> insert( newHit );

  return true;
}

void SensitiveDetector_SiMon::EndOfEvent(G4HCofThisEvent*)
{

 // Initialisation of total energy deposition per event to zero
 G4double totalEdep = 0.0*MeV;
 G4double Ekin = 0.0*MeV;
 G4int Z = 0;
 G4int A = 0;
 
 G4int NbHits = fHitsCollection->entries();
  
  for (G4int i=0;i<NbHits;i++)
  {
    G4double edep = (*fHitsCollection)[i]->GetEdep();
    totalEdep     = totalEdep + edep;
    Ekin  = (*fHitsCollection)[0]->GetEkin();
    A     = (*fHitsCollection)[0]->GetBaryonNumber();
    Z     = (*fHitsCollection)[0]->GetParticleCharge();

  } 
  G4double sigmaE = 0.01;
  totalEdep = G4RandGauss::shoot(totalEdep, sigmaE); 
  if(totalEdep>0.1) analysis->Detector_SiMon(totalEdep, Ekin, A, Z);
}

