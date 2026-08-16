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

#include "EventAction.hh"
#include "AnalysisManager.hh"
#include "DetectorConstruction.hh"

#include "G4Event.hh"
#include <algorithm>

EventAction::EventAction(AnalysisManager* analysisMan, DetectorConstruction* det)
: G4UserEventAction(),
  analysis(analysisMan),
  detector(det)
{}

void EventAction::BeginOfEventAction(const G4Event*)
{
  fTargetByTrack.clear();
  fTargetOrder.clear();

  fLSByKey.clear();
  fLSOrder.clear();
  fCrossings.clear();
  fCrossingOrder.clear();
}

void EventAction::AddLSNeutronCrossing(G4int detID, G4int trackID, G4int parentID,
                                       G4double ekin, G4double tof,
                                       const G4ThreeVector& pos,
                                       const G4ThreeVector& dir)
{
  const auto key = LSKey(detID, trackID);
  if (fCrossings.find(key) != fCrossings.end()) return;

  NeutronCrossingRecord record;
  record.detID = detID;
  record.trackID = trackID;
  record.parentID = parentID;
  record.ekin = ekin;
  record.tof = tof;
  record.pos = pos;
  record.dir = dir;
  fCrossings.emplace(key, record);
  fCrossingOrder.push_back(key);
}

// ---------- Target/strip recorders ----------
void EventAction::AddTargetNeutron(G4int trackID, G4int parentID, G4double ekin, const G4ThreeVector& dir, const G4ThreeVector& pos, G4double tof)
{
  auto it = fTargetByTrack.find(trackID);
  if (it == fTargetByTrack.end())
  {
    TargetRecord rec;
    rec.trackID = trackID;
    rec.parentID = parentID;
    fTargetOrder.push_back(trackID);
    it = fTargetByTrack.emplace(trackID, rec).first;
  }

  auto& rec = it->second;
  if (rec.hasTarget) return;

  rec.hasTarget = true;
  rec.ekin_tar = ekin;
  rec.dir_tar  = dir;
}

void EventAction::AddLSDetected(G4int detID, G4int trackID, G4int parentID, G4int Z_, G4int A_, G4double edep, G4double /*tof*/)
{
  const auto key = LSKey(detID, trackID);

  auto it = fLSByKey.find(key);
  if (it == fLSByKey.end())
  {
    LSHitRecord r;
    r.detID = detID;
    r.trackID = trackID;
    r.parentID = parentID;
    fLSOrder.push_back(key);
    it = fLSByKey.emplace(key, r).first;
  }

  auto& r = it->second;
  r.hasDep = true;
  r.edep_sum += edep;
  r.Z_dep = Z_;
  r.A_dep = A_;
}

// ---------- End of event: write per-track TargetTrack + per-hit LSHitTrack ----------
void EventAction::EndOfEventAction(const G4Event* evt)
{
  const G4int eventID = evt->GetEventID();

  // 1) per-track target/strip record
  for (const auto trackID : fTargetOrder)
  {
    const auto it = fTargetByTrack.find(trackID);
    if (it == fTargetByTrack.end()) continue;

    const TargetRecord& r = it->second;
    if (!r.hasTarget) continue;
    analysis->TargetTrack(eventID, r.trackID, r.ekin_tar,
                          r.dir_tar.theta(), r.dir_tar.phi());
  }

  for (const auto key : fCrossingOrder)
  {
    const auto it = fCrossings.find(key);
    if (it == fCrossings.end()) continue;
    const auto& n = it->second;
    analysis->LSNeutronCrossing(eventID, n.detID, n.trackID, n.parentID,
                                n.ekin, n.tof,
                                n.pos.x(), n.pos.y(), n.pos.z(),
                                n.dir.theta(), n.dir.phi());
  }

  // 2) per (detID, trackID) LS hit record
  for (const auto key : fLSOrder)
  {
    const auto it = fLSByKey.find(key);
    if (it == fLSByKey.end()) continue;

    const LSHitRecord& h = it->second;

    // Optional: only write if there was some deposit
    // if (!h.hasDep) continue;

    // Write ONLY if total deposited energy for that (detID,trackID) is above threshold
    const G4double LS_TOTAL_EDEP_THR = 0.001; // choose your threshold in MeV (or whatever units you're using)
    if (!h.hasDep) continue;
    if (h.edep_sum < LS_TOTAL_EDEP_THR) continue;
    
    analysis->LSHitTrack(eventID,
                         h.detID,
                         h.trackID,
                         h.parentID,
                         h.edep_sum,
                         h.Z_dep, h.A_dep);
  }
}
