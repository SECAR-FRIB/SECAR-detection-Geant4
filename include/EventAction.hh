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
// 

#ifndef EventAction_h
#define EventAction_h 1

#include "G4UserEventAction.hh"
#include "globals.hh"
#include "G4ThreeVector.hh"

#include <unordered_map>
#include <vector>

class AnalysisManager;
class DetectorConstruction;

class EventAction : public G4UserEventAction
{
public:
  EventAction(AnalysisManager*, DetectorConstruction*);
  ~EventAction() override = default;

  void BeginOfEventAction(const G4Event* evt) override;
  void EndOfEventAction(const G4Event* evt) override;

  // ---- called from SteppingAction ----
  void AddTargetNeutron(G4int trackID, G4int parentID,
                        G4double ekin,
                        const G4ThreeVector& dir,
                        const G4ThreeVector& pos,
                        G4double tof);

  void AddLSDetected(G4int detID, G4int trackID, G4int parentID,
                     G4int Z, G4int A,
                     G4double edep, G4double tof);

  void AddLSNeutronCrossing(G4int detID, G4int trackID, G4int parentID,
                            G4double ekin, G4double tof,
                            const G4ThreeVector& pos,
                            const G4ThreeVector& dir);

private:
  AnalysisManager* analysis = nullptr;
  DetectorConstruction* detector = nullptr;

  // ------------ Target neutron per-track record ------------
  struct TargetRecord {
    bool hasTarget = false;

    G4int trackID  = -1;
    G4int parentID = -1;

    G4double ekin_tar = 0.0;
    G4ThreeVector dir_tar;
  };

  std::unordered_map<G4int, TargetRecord> fTargetByTrack;
  std::vector<G4int> fTargetOrder;

  // ------------ LS per (detID, trackID) record ------------
  struct LSHitRecord {
    G4int detID   = -1;
    G4int trackID = -1;
    G4int parentID = -1;

    // detected energy (sum over steps)
    G4double edep_sum = 0.0;

    G4int Z_dep = 0;
    G4int A_dep = 0;

    bool hasDep   = false;
  };

  // encode (detID, trackID) into one key for unordered_map
  static inline std::uint64_t LSKey(G4int detID, G4int trackID) {
    return (std::uint64_t(std::uint32_t(detID)) << 32) | std::uint32_t(trackID);
  }

  std::unordered_map<std::uint64_t, LSHitRecord> fLSByKey;
  std::vector<std::uint64_t> fLSOrder; // stable iteration

  struct NeutronCrossingRecord {
    G4int detID = -1;
    G4int trackID = -1;
    G4int parentID = -1;
    G4double ekin = 0.0;
    G4double tof = 0.0;
    G4ThreeVector pos;
    G4ThreeVector dir;
  };
  std::unordered_map<std::uint64_t, NeutronCrossingRecord> fCrossings;
  std::vector<std::uint64_t> fCrossingOrder;
};

#endif

    
