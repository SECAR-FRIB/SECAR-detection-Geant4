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

#include "PrimaryGeneratorAction.hh"
#include "AnalysisManager.hh"
#include "DetectorConstruction.hh"
#include "G4Event.hh"
#include "G4Step.hh"
#include "G4GeneralParticleSource.hh"
#include "G4ParticleGun.hh"
#include "G4ParticleTable.hh"
#include "G4ParticleDefinition.hh"
#include "G4SystemOfUnits.hh"
#include "G4IonTable.hh"
#include "G4Geantino.hh"
#include "G4PhysicalConstants.hh"
#include "Randomize.hh"

#include <sstream>
#include <vector>
#include <cstdlib>

PrimaryGeneratorAction::PrimaryGeneratorAction(AnalysisManager* pAnalysis)
{
    analysis = pAnalysis;  

    if(source)
    {
        fParticleGun  = new G4ParticleGun(1);
        fParticleGun->SetParticleEnergy(0*eV);
        fParticleGun->SetParticlePosition(G4ThreeVector(0.,0.,-31.*mm));
        fParticleGun->SetParticleMomentumDirection(G4ThreeVector(0.,0.,1.)); 
    }
    if(neutron || two_neutrons)
    {
        n=1;
        fParticleGun  = new G4ParticleGun(1);
        fParticleGun->SetParticlePosition(G4ThreeVector(0.,0.,0.));
        fParticleGun->SetParticleMomentumDirection(G4ThreeVector(0.,0.,1.));
        G4ParticleTable* particleTable = G4ParticleTable::GetParticleTable();
        G4ParticleDefinition* particle = particleTable->FindParticle("neutron");
        fParticleGun->SetParticleDefinition(particle);
    }
    if(ions)
    {
        n=1;
        fParticleGun  = new G4ParticleGun(1);
        fParticleGun->SetParticlePosition(G4ThreeVector(0.,0.,0.));
        fParticleGun->SetParticleMomentumDirection(G4ThreeVector(0.,0.,1.));
        fParticleGun->SetParticleEnergy(0*MeV);
    }
    else gps = new G4GeneralParticleSource();
}

PrimaryGeneratorAction::~PrimaryGeneratorAction()
{
    if(source) delete fParticleGun; 
    else delete gps;
}	

void PrimaryGeneratorAction::GeneratePrimaries(G4Event* anEvent)
{
    if(source)
    {
        fParticleGun->SetParticlePosition(G4ThreeVector(0,0,-30.*mm));

        // Distribution uniform in solid angle
        G4double cosTheta = 2*G4UniformRand() - 1.; 
        G4double sinTheta = std::sqrt(1. - cosTheta*cosTheta);
        G4double phi = twopi*G4UniformRand(); 
        G4double vx = sinTheta*std::cos(phi),
                 vy = sinTheta*std::sin(phi),
                 vz = cosTheta;

        fParticleGun->SetParticleMomentumDirection(G4ThreeVector(vx,vy,vz));

        fParticleGun->GeneratePrimaryVertex(anEvent);
    }
    else if(neutron)
    {  
        bool rewind = false;
        G4double energy_n, theta_n;
        if (!ReadNeutronRow(energy_n, theta_n, rewind)) return;
        // This here is so if need to rewind the file (more events than the file)
        // the are not exactly the same introducing a maximum of 1 keV difference in energy
        // and a maximum of 0.1 mrad difference in angle
        G4double energy = energy_n;
        G4double theta  = theta_n;
        if (rewind)
        {
            // Energy smearing ±1 keV
            G4double dE = 0.001; // MeV = 1 keV
            G4double Emin = std::max(0.0, energy - dE);
            G4double Emax = energy + dE;
            energy_n = CLHEP::RandFlat::shoot(Emin, Emax);

            // Angle smearing ±0.1 mrad
            G4double dth = 0.0001; // rad = 0.1 mrad
            const G4double eps = 1e-6;    // rad (keeps away from 0)
            G4double thmin = std::max(eps, theta - dth);
            G4double thmax = std::min(pi- eps, theta + dth);
            if (thmin >= thmax) {
                theta_n = std::min(std::max(theta, eps), pi - eps); // just clamp
            } else {
                theta_n = CLHEP::RandFlat::shoot(thmin, thmax);
            }
        }

        fParticleGun->SetParticlePosition(G4ThreeVector(0,0,0));
        // Uniform spherical distribution
        // Converting the spherical coordinates (polar and azimuthal angles) to Cartesian coordinates. 
        G4double phi = CLHEP::RandFlat::shoot(0.0, 2.0*pi);
        G4double vx = sin(theta_n) * cos(phi);
        G4double vy = sin(theta_n) * sin(phi);
        G4double vz = cos(theta_n);

        fParticleGun->SetParticleMomentumDirection(G4ThreeVector(vx,vy,vz));
        fParticleGun->SetParticleEnergy(energy_n*MeV);  
        fParticleGun->GeneratePrimaryVertex(anEvent);
    }
    else if(ions)
    {  
        if (fParticleGun->GetParticleDefinition() == G4Geantino::Geantino()) {
            G4int Z = 27 , A = 58;
            G4double ionCharge = 0.*eplus;
            G4double excitEnergy = 0.*keV;
            G4ParticleDefinition* ion = G4IonTable::GetIonTable()->GetIon(Z,A,excitEnergy);
            fParticleGun->SetParticleCharge(ionCharge);
            fParticleGun->SetParticleDefinition(ion);
        }

        G4double energy, theta;
        if(n==1) OpenExternalFile();
        
        input_file >> theta >> energy;

        G4double theta_r = 0.0;
        G4double energy_r = 0.0;
        if(energy>0) energy_r = CLHEP::RandFlat::shoot(energy-0.002, energy+0.002);
        else energy_r = CLHEP::RandFlat::shoot(energy, energy+0.002);
        if(theta>0) theta_r = CLHEP::RandFlat::shoot(theta-0.002, theta+0.002);
        else theta_r = theta+0.000001;


        fParticleGun->SetParticlePosition(G4ThreeVector(0,0,-20*cm));
        // Uniform spherical distribution
        // Converting the spherical coordinates (polar and azimuthal angles) to Cartesian coordinates. 
        G4double phi = CLHEP::RandFlat::shoot(0.0,360.);
        G4double vx = sin(theta_r*deg)*cos(phi*deg);
        G4double vy = sin(theta_r*deg)*sin(phi*deg);
        G4double vz = cos(theta_r*deg);

        fParticleGun->SetParticleMomentumDirection(G4ThreeVector(vx,vy,vz));
        fParticleGun->SetParticleEnergy(energy_r*MeV);  
        fParticleGun->GeneratePrimaryVertex(anEvent);
        n++;
        //if(n==nmax) {n=0; input_file.close();} //This is the line where the file ends
        if(input_file.eof()){ input_file.close(); n=1; } 
    }
    else if(two_neutrons)
    {
        if (!input_file.is_open()) OpenExternalFile();

        char comma;
        bool rewound = false;

        G4double E1, th1, E2, th2;

        // Try to read BOTH lines as one event
        if (!(input_file >> E1 >> comma >> th1 >> E2 >> comma >> th2))
        {
            // rewind only here (at event boundary)
            input_file.clear();
            input_file.seekg(0, std::ios::beg);

            std::string header;
            std::getline(input_file, header);

            rewound = true;

            if (!(input_file >> E1 >> comma >> th1 >> E2 >> comma >> th2))
            {
                G4Exception("PrimaryGeneratorAction", "BadInput", FatalException,
                            "Could not read a full 2-neutron event (two lines) after rewind.");
                return;
            }
        }

        if (rewound)
        {
            SmearIfRewound(E1, th1);
            SmearIfRewound(E2, th2);
        }

        const G4ThreeVector vtx(0,0,0);
        fParticleGun->SetParticlePosition(vtx);

        auto shoot = [&](G4double E, G4double theta)
        {
            G4double phi = CLHEP::RandFlat::shoot(0.0, 2.0*pi);
            G4double vx = std::sin(theta) * std::cos(phi);
            G4double vy = std::sin(theta) * std::sin(phi);
            G4double vz = std::cos(theta);

            fParticleGun->SetParticleMomentumDirection(G4ThreeVector(vx,vy,vz));
            fParticleGun->SetParticleEnergy(E*MeV);
            fParticleGun->GeneratePrimaryVertex(anEvent);
        };
        shoot(E1, th1);
        shoot(E2, th2);

        /* Diagnostics print
        G4cout << "Evt " << anEvent->GetEventID()
            << " : (E1,th1)=(" << E1 << "," << th1 << ")"
            << " (E2,th2)=(" << E2 << "," << th2 << ")"
            << (rewound ? " [REWOUND]" : "")
            << G4endl;*/
    }
    else gps->GeneratePrimaryVertex(anEvent);
}
void PrimaryGeneratorAction::OpenExternalFile()
{   
    // Open neutron distribution file
    if(neutron || two_neutrons)
    {
        const char* configuredPath = std::getenv("SECAR_NEUTRON_INPUT");
        const std::string inputPath = configuredPath
            ? configuredPath
            : "../neutronInput/neutron_input.csv";
        input_file.open(inputPath, std::ios::in);
        
        if (!input_file.is_open()) {
            G4Exception("PrimaryGeneratorAction", "MissingInput", FatalException,
                        "Could not open the neutron CSV. Set SECAR_NEUTRON_INPUT "
                        "or provide ../neutronInput/neutron_input.csv.");
            return;
        }

        std::string header;
        std::getline(input_file, header);  // skip CSV header
    }
    if(ions)
    input_file.open("input/input_rec_3_211.txt",std::ios::in);
    
    // !!!!! Important !!!!! Make sure the file doesn't have an empty line at the end
    // If there is remove it, it will create warnings while running !!!!

    // if (input_file.is_open() && n==1) G4cout << "Neutrons are genereted using external file!! "<< G4endl;
    // else G4cout << "File not found!! "<< G4endl;
    
}

bool PrimaryGeneratorAction::ReadNeutronRow(G4double& energy_n,
                                            G4double& theta_n,
                                            bool& rewound)
{
    if (!input_file.is_open()) OpenExternalFile();
    if (!input_file.is_open()) return false;

    for (G4int attempt = 0; attempt < 2; ++attempt)
    {
        std::string line;
        while (std::getline(input_file, line))
        {
            if (line.empty()) continue;

            std::vector<G4double> fields;
            std::stringstream row(line);
            std::string field;
            try {
                while (std::getline(row, field, ',')) fields.push_back(std::stod(field));
            } catch (const std::exception&) {
                G4Exception("PrimaryGeneratorAction", "BadInput", FatalException,
                            "A neutron CSV row contains a non-numeric field.");
                return false;
            }

            if (fields.size() >= 4) {
                energy_n = fields[2];
                theta_n = fields[3];
                return true;
            }
            if (fields.size() >= 2) {
                energy_n = fields[0];
                theta_n = fields[1];
                return true;
            }

            G4Exception("PrimaryGeneratorAction", "BadInput", FatalException,
                        "A neutron CSV row must contain either 2 or at least 4 columns.");
            return false;
        }

        input_file.clear();
        input_file.seekg(0, std::ios::beg);
        std::string header;
        std::getline(input_file, header);
        rewound = true;
    }

    G4Exception("PrimaryGeneratorAction", "BadInput", FatalException,
                "Could not read a neutron data row after rewinding the CSV file.");
    return false;
}
void PrimaryGeneratorAction::SmearIfRewound(G4double& energy_n, G4double& theta_n)
{
    // Energy smearing ±1 keV
    const G4double dE = 0.001; // MeV
    const G4double Emin = std::max(0.0, energy_n - dE);
    const G4double Emax = energy_n + dE;
    energy_n = CLHEP::RandFlat::shoot(Emin, Emax);

    // Angle smearing ±0.1 mrad
    const G4double dth = 0.0001; // rad
    const G4double eps = 1e-6;
    const G4double thmin = std::max(eps, theta_n - dth);
    const G4double thmax = std::min(pi - eps, theta_n + dth);

    if (thmin >= thmax) theta_n = std::min(std::max(theta_n, eps), pi - eps);
    else                theta_n = CLHEP::RandFlat::shoot(thmin, thmax);
}
