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
//

#include "SteppingAction.hh"
#include "AnalysisManager.hh"
#include "EventAction.hh"
#include "DetectorConstruction.hh" 
#include "TrackInfo.hh"

//#include "G4ios.hh"
//#include "G4SteppingManager.hh"
#include "G4Step.hh"
#include "G4Track.hh"
//#include "G4ParticleTypes.hh"
#include "G4ParticleDefinition.hh"
#include "G4StepPoint.hh"
#include "G4StepStatus.hh"
#include "G4VTouchable.hh"
//#include "G4ThreeVector.hh"
#include "G4VPhysicalVolume.hh"
#include "G4SystemOfUnits.hh"
#include "G4VProcess.hh"
#include "G4Neutron.hh"

namespace {

  inline TrackInfo* GetOrMakeTrackInfo(G4Track* t)
  {
    auto* info = static_cast<TrackInfo*>(t->GetUserInformation());
    if (!info) {
      info = new TrackInfo();
      t->SetUserInformation(info);
    }
    return info;
  }

  inline bool IsNeutron(const G4Track* t) {
    return t->GetDefinition() == G4Neutron::Definition();
  }

  inline bool IsIon(const G4Track* t) {
    // For ions, atomic number/mass are meaningful.
    // For neutrons/protons/gammas they are 0.
    return t->GetDefinition()->GetParticleType() == "nucleus";
  }

inline G4int GetZ(const G4Track* t)
{
  const auto* def = t->GetDefinition();

  // For ions
  if (def->GetParticleType() == "nucleus")
    return def->GetAtomicNumber();

  // For baryons (proton, neutron)
  const G4double charge = def->GetPDGCharge();   // in units of eplus
  return (G4int)std::lround(charge / eplus);
}

inline G4int GetA(const G4Track* t)
{
  const auto* def = t->GetDefinition();

  // For ions
  if (def->GetParticleType() == "nucleus")
    return def->GetAtomicMass();

  // For baryons (proton, neutron)
  return def->GetBaryonNumber();
}

  inline bool IsTargetVolumeName(const G4String& name)
  {
    // support both names
    return (name == "Sphere_phys" || name == "Tar_phys");
  }
}

SteppingAction::SteppingAction(EventAction* eventAction, DetectorConstruction* det, AnalysisManager* pAnalysis)
: G4UserSteppingAction(),
  fEventAction(eventAction),
  detector(det),
  analysis(pAnalysis)  
{}

//SteppingAction::~SteppingAction(){}

void SteppingAction::UserSteppingAction(const G4Step* aStep)
{ 
  if (!aStep) return;
  
  G4Track* aTrack = aStep -> GetTrack();
  if (!aTrack) return;
  G4StepPoint* initialPoint = aStep->GetPreStepPoint();
  if (!initialPoint) return;
  G4VPhysicalVolume* volume = initialPoint->GetTouchableHandle()->GetVolume();
  if (!volume) return;

  const G4String volumeName = volume->GetName();

  // Track metadata / kinematics at this step
  auto* info = GetOrMakeTrackInfo(aTrack);

  const G4int trackID  = aTrack->GetTrackID();
  const G4int parentID = aTrack->GetParentID();

  const G4ThreeVector pos = initialPoint->GetPosition();
  const G4ThreeVector dir = aTrack->GetMomentumDirection();
  const G4double ekin      = aTrack->GetKineticEnergy();
  const G4double tof       = aTrack->GetGlobalTime();

  const bool neutron = IsNeutron(aTrack);
  const bool ion     = IsIon(aTrack);

  const G4int Z = GetZ(aTrack);
  const G4int A = GetA(aTrack);
  
  // ------------------------------------------------------------
  // 1) TARGET: "what is produced at target" (source or reaction)
  // Record ONCE per track when it is in target volume.
  // This works for:
  //   - primaries spawned inside target
  //   - particles entering target
  //   - secondaries born inside target
  // ------------------------------------------------------------
  if (IsTargetVolumeName(volumeName))
  {
    if (!info->recordedTarget)
    {
      info->recordedTarget = true;
      // Store a "target record" per track
      if (neutron)
      {
        fEventAction->AddTargetNeutron(trackID, parentID, ekin, dir, pos, tof);
      }
    }
  }  
  // ------------------------------------------------------------
  // 3) LIQUID SCINTILLATORS
  // Record entry kinematics at boundary; record edep per step
  // Use trackID for later correlation (no single “neutron at target” scalars).
  // ------------------------------------------------------------
  // Flag for turning on/off the Liquid Scintillators. 
  // It is controlled from the header file of the DetectorConstruction 
  if (detector && detector->LScinOn())
  { 
    const G4int LScinNum = detector->GetLScinNum();
    for(int i=0; i<LScinNum; i++)
    {
      if (volume != detector->GetLScin(i)) continue;
      if (neutron && initialPoint->GetStepStatus() == fGeomBoundary && ekin > 0.0001)
        fEventAction->AddLSNeutronCrossing(i, trackID, parentID, ekin, tof, pos, dir);

      // Accumulate energy deposition for THIS track in THIS LS
      const G4double edep = aStep->GetTotalEnergyDeposit();
      if (edep <= 0.0) continue;

      // Detect charged particles (includes ions, protons, electrons, etc.)
      //const G4double q = aTrack->GetDefinition()->GetPDGCharge();
      if (Z == 0) continue;

      // IMPORTANT: no per-step threshold here; accumulate everything
      fEventAction->AddLSDetected(i, trackID, parentID, Z, A, edep, tof);
    }
  }
}
