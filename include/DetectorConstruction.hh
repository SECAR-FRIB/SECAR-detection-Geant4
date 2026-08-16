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

#ifndef DetectorConstruction_H 
#define DetectorConstruction_H 1

#include "G4VUserDetectorConstruction.hh"
#include "globals.hh"
#include "AnalysisManager.hh"

#include <fstream>

class G4Box;
class G4Tubs;
class G4LogicalVolume;
class G4VPhysicalVolume;
class G4Material;
class DetectorMessenger;

class DetectorConstruction : public G4VUserDetectorConstruction
{
public:
    DetectorConstruction(AnalysisManager* analysis);
    ~DetectorConstruction();

    G4VPhysicalVolume* Construct();

    void ConstructSDandField();
    void NeutronDetectorsMap();
	
    G4VPhysicalVolume* GetLScin(int i) const {return fLScin_phys[i];};
    G4int GetLScinNum() const {return LScinNum;};
    G4bool LScinOn() const {return LScin;};

private:
  
	AnalysisManager* analysis;

    //Enable your detectors
    G4bool chamber = false;
    G4bool source = false; // !! To use 4Pi emmision change also in PrimaryGeneratorAction.hh !!
    // When source is false all of the covering material is set to vaccum, thus one can still get what was 
    // emmitted from the target location.
    G4bool target = false;
    G4bool Si_monitor = false;
    G4bool stripper = false;
    G4bool LScin = true;
    G4bool coverLS = true;
    G4bool MCP = false;
    G4bool IC = false;
        G4double gas_thick_one = 12.7; //12.7*mm;
        G4double gas_thick_dE = 50.8; //50.8*mm
        G4double gas_thick_two = 3.2; //3.2*mm 
    G4bool DSSD = false;
    //G4bool LENDA = false;
    
    //Reading the neutron detectors' map
    G4String name;
    G4double xpos, ypos, zpos, xcovfpos, ycovfpos, xcovbpos, ycovbpos;
    G4double xangle, yangle, zangle;
    std::string line;

    G4int LScinNum;
    std::ifstream LScinMap;
    std::vector<G4String>  LScinName;
    std::vector<G4double>  LScinXpos;
    std::vector<G4double>  LScinYpos;
    std::vector<G4double>  LScinCovFXpos;
    std::vector<G4double>  LScinCovFYpos;
    std::vector<G4double>  LScinCovBXpos;
    std::vector<G4double>  LScinCovBYpos;
    std::vector<G4double>  LScinZpos;
    std::vector<G4double>  LScinXangle;
    std::vector<G4double>  LScinYangle;
    std::vector<G4double>  LScinZangle;
    std::vector<G4VPhysicalVolume*> fLScin_phys;
    std::vector<G4VPhysicalVolume*> coverLScin_phys;
    std::vector<G4VPhysicalVolume*> coverbLScin_phys;
    std::vector<G4VPhysicalVolume*> coverfLScin_phys;
    
};
#endif


