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
// Authors: Pelagia Tsintar, pelagia.tsin@gmail.com
// 
//


#ifndef PrimaryGeneratorAction_hh
#define PrimaryGeneratorAction_hh 1

#include "G4VUserPrimaryGeneratorAction.hh"
#include "globals.hh"
#include "AnalysisManager.hh"
#include "DetectorConstruction.hh"

#include <fstream>

using namespace std; 


class G4GeneralParticleSource;
class G4ParticleGun;
class G4Event;
class PrimaryGeneratorAction : public G4VUserPrimaryGeneratorAction
{
public:
PrimaryGeneratorAction(AnalysisManager*);
~PrimaryGeneratorAction();

virtual void GeneratePrimaries(G4Event* anEvent);
void OpenExternalFile();
private:
    void SmearIfRewound(G4double& energy_n, G4double& theta_n);
    bool ReadNeutronRow(G4double& energy_n, G4double& theta_n, bool& rewound);

// method to access particle gun
G4ParticleGun* GetParticleGun() {return fParticleGun;};

private:
AnalysisManager* analysis;
G4GeneralParticleSource* gps;
G4ParticleGun*  fParticleGun;
G4bool source = false;
G4bool neutron = true;
G4bool two_neutrons = false;
G4bool ions = false;
ifstream input_file;
G4int n;

};
#endif
