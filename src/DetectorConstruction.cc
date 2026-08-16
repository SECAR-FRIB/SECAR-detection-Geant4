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

#include "DetectorConstruction.hh"
#include "SensitiveDetector.hh"
#include "SensitiveDetector_IC.hh"

#include "globals.hh"
#include "G4NistManager.hh"
#include "G4Element.hh"
#include "G4Material.hh"

#include "G4Box.hh"
#include "G4Tubs.hh"
#include "G4Sphere.hh"
#include "G4Cons.hh"
#include "G4PVPlacement.hh"
#include "G4LogicalVolume.hh"
#include "G4LogicalVolumeStore.hh"
#include "G4SubtractionSolid.hh"
#include "G4Color.hh"
#include "G4VisAttributes.hh"
#include "G4SDManager.hh"

#include "Randomize.hh"
#include "G4ThreeVector.hh"
#include "G4RotationMatrix.hh"

#include "G4SystemOfUnits.hh"

//#include "CADMesh.hh" //uncomment if you want to use the CAD drawing of JENSA chamber - comment out lines 127 - 187

DetectorConstruction::DetectorConstruction(AnalysisManager* analysis_manager)
:G4VUserDetectorConstruction(),
  fLScin_phys(0)
{
  analysis = analysis_manager;
  LScinNum = -1;
}

DetectorConstruction::~DetectorConstruction(){}

G4VPhysicalVolume* DetectorConstruction::Construct()
{
  analysis->ConfigureDetectorOutput(DSSD, IC, Si_monitor);

  G4NistManager* nist = G4NistManager::Instance();

  G4Material* vacuum = nist->FindOrBuildMaterial("G4_Galactic");
  //G4Material* air = nist->FindOrBuildMaterial("G4_AIR");
  G4Material* aluminum = nist->FindOrBuildMaterial("G4_Al");
  G4Material* mylar = nist->FindOrBuildMaterial("G4_MYLAR");
  //G4Material* kapton = nist->FindOrBuildMaterial("G4_KAPTON");
  G4Material* silicon = nist->FindOrBuildMaterial("G4_Si");
  G4Material* gold = nist->FindOrBuildMaterial("G4_Au");
  G4Material* nickel = nist->FindOrBuildMaterial("G4_Ni");
  G4Material* carbon = nist->FindOrBuildMaterial("G4_C");
  G4Material* stainless = nist->FindOrBuildMaterial("G4_STAINLESS-STEEL");
  G4Material* polyethylene = nist->FindOrBuildMaterial("G4_POLYETHYLENE");
  G4Material* pyrex = nist->FindOrBuildMaterial("G4_Pyrex_Glass");
  G4Material* HDPE = new G4Material("hdpe",0.93*g/cm3,1);
    HDPE->AddMaterial(polyethylene,1);
  //G4Material* carbon = nist->FindOrBuildMaterial("G4_C");
  
  G4Element* C = nist->FindOrBuildElement(6);
  G4Element* H = nist->FindOrBuildElement(1);
  G4Element* He = nist->FindOrBuildElement(2);
  G4Element* F = nist->FindOrBuildElement(9);

  // (!!!) Note: Geant needs the mass ratio which is different than the atomic ratio!!
  // Atomic H/C = 1.212 => (0.092/1.008u) / (0.908/12.011u)
  G4Element* elHydrogen = new G4Element("Hydrogen","H",1.,1.00794*g/mole);
  G4Element* elCarbon = new G4Element("Carbon","C",6.,12.0107*g/mole);
  G4Material* EJ301 = new G4Material("EJ301", 0.874*g/cm3, 2, kStateLiquid);
    EJ301->AddElement(elHydrogen, 9.2*perCent); 
    EJ301->AddElement(elCarbon, 90.8*perCent); 

  // test with all hydrogen target
  G4Material* h_tar = new G4Material("h_tar", 3.5841e-7*g/cm3, 1, kStateSolid); // 6.5 torr = 3.5841e-07 g/cm3
    h_tar->AddElement(elHydrogen,1); 

  G4Material* c_tar = new G4Material("c_tar", 1.15193*g/cm3, 1, kStateSolid);
    c_tar->AddElement(elCarbon,1); 

  G4Material* isobutane = new G4Material("isobutane",0.000168*g/cm3,2,kStateGas,293.15*kelvin); 
    isobutane->AddElement(C,82.64*perCent); 
    isobutane->AddElement(H,17.36*perCent); 
    //0.00015=45.9, 0.0001615=50.8Torr 0.000159=50Torr, 0.00013=40Torr, 0.00017=55Torr
    // multiply all lise densities with 0.00249343/0.0024164 = 1.031878
    //0.00014446292=45, 0.00016510048=50

  G4Material* cf4 = new G4Material("CF4",9.6279e-5*g/cm3,2,kStateGas,293.15*kelvin); // 20 torr = 9.6279e-05 g/cm3
    cf4->AddElement(C,1); //);
    cf4->AddElement(F,4); //);
    //0.00015=45.9, 0.0001615=50.8Torr 0.000159=50Torr, 0.00013=40Torr, 0.00017=55Torr
    // multiply all lise densities with 0.00249343/0.0024164 = 1.031878
    //0.00014446292=45, 0.00016510048=50
  
  /*G4Material* stainless = new G4Material("stainless-steel",0.2*g/cm3,2); 
    stainless->AddElement(Fe,89*perCent);
    stainless->AddElement(Ni,11*perCent);
    //stainless->AddElement(Ni,8*perCent);*/

  G4Material* helium = new G4Material("HeliumGas",0.00272*g/cm3,1,kStateGas,293.15*kelvin); //,16.331*atmosphere
  helium->AddElement(He, 1.0);

  G4RotationMatrix* rotX90 = new G4RotationMatrix(); rotX90->rotateX(90*deg);
  G4RotationMatrix* rotY90 = new G4RotationMatrix(); rotY90->rotateY(90*deg);
  G4RotationMatrix* rotZ90 = new G4RotationMatrix(); rotZ90->rotateZ(90*deg);
  G4RotationMatrix* rotY50 = new G4RotationMatrix(); rotY50->rotateY(50*deg);
  G4RotationMatrix* rotY310 = new G4RotationMatrix(); rotY310->rotateY(310*deg);

  //Define volumes
  G4double worldx = 300.0*cm; 
  G4double worldy = 300.0*cm;
  G4double worldz = 300.0*cm;
  G4Material *world_mat;
  //if(source) world_mat = air;
  world_mat = vacuum;

  // World volume, containing all geometry
  G4Box* world = new G4Box("world_box", worldx, worldy, worldz);
  G4LogicalVolume* logical_world = new G4LogicalVolume(world, world_mat, "world_log", 0,0,0);
  G4VPhysicalVolume* physical_world = new G4PVPlacement(0,G4ThreeVector(),logical_world,"world_phys",0,false, true);
  //set the logical world volume invisible
  logical_world -> SetVisAttributes(G4VisAttributes::GetInvisible());

  if(target)
  {
    //Target
    // G4double Tar_x = 6.00*mm;
    // G4double Tar_y = 6.00*mm;
    // G4double Tar_z = (0.5*3.914)*um; //4.2*um ;//4.2021505*um = 0.3908 mg/cm2; //3.7=HDPE=6.8, 1.82 C
    //4.62*um; used during the experiment (it was wrong)
    // G4double Tar_z = 137.0*mm; //Extended gas target effective thickness He
    // G4Box* Target = new G4Box("Target",Tar_x,Tar_y,Tar_z/2);

    // Jet gas target
    G4double Tar_d = 3.30*mm; // nozzle E
    G4double Tar_z = 10*mm; // distance nozzle to catcher - not sure

    G4Tubs* Target = new G4Tubs("Target", 0, Tar_d/2, Tar_z/2, 0, 360*deg);
    G4LogicalVolume* logical_Tar = new G4LogicalVolume(Target, helium, "Tar_log", 0,0,0); 

    //Visualisation attributes
    new G4PVPlacement(rotX90, G4ThreeVector(0,0,0), logical_Tar, "Tar_phys", logical_world, false, 0, true);
    G4VisAttributes vis_Tar(G4Colour(250, 0, 127)); //numbers stand for colors of red, green, and blue respectively
    vis_Tar.SetForceSolid(true);
    logical_Tar -> SetVisAttributes(vis_Tar);
    
    if(stripper)
    {
      // center target to stripper foil = 481.6 mm
      G4double foil_thick = 0.1864*um;
      G4double dist_to_foil = 413.1*mm;
      G4Box* Stripper_foil = new G4Box("Stripper_foil",5.0*cm,5.0*cm,foil_thick/2);
      G4LogicalVolume* logical_StripFoil = new G4LogicalVolume(Stripper_foil, carbon, "Strip_log", 0,0,0);
      new G4PVPlacement(0, G4ThreeVector(0,0,dist_to_foil), logical_StripFoil, "Strip_phys", logical_world, false, 0, true);
      // G4VisAttributes vis_Stripper(G4Colour(23, 255, 255)); 
      // vis_Stripper.SetForceCloud(true);
      // logical_Tar -> SetVisAttributes(vis_Stripper);
    }
  }

  //Source
  G4Material *source_cov_mat;
  G4Material *source_front_mat;
  if(source)
  {
    source_cov_mat = nickel;
    source_front_mat = gold;
  }
  else 
  {
    source_cov_mat = vacuum;
    source_front_mat = vacuum;
  }

  G4double source_rmin = (15.87/2)*mm;
  G4double source_rmax = (31.75/2)*mm;
  G4double source_z = 2.5*mm;

  G4Tubs* source_vol = new G4Tubs("source", source_rmin, source_rmax, source_z, 0, 360*deg);
  G4LogicalVolume* log_source = new G4LogicalVolume(source_vol, source_cov_mat, "log_source", 0,0,0); //nickel
  //new G4PVPlacement(0, G4ThreeVector(0,0,-30*mm), log_source,"source_phys", logical_world, false, 0, true);
  log_source->SetVisAttributes(G4Color::Red());
  G4Tubs* source_b = new G4Tubs("source", 0, source_rmax, source_z, 0, 360*deg);
  G4LogicalVolume* log_source_b = new G4LogicalVolume(source_b, source_cov_mat, "log_source_b", 0,0,0); //nickel
  //new G4PVPlacement(0, G4ThreeVector(0,0,-35*mm), log_source_b,"source_back_phys", logical_world, false, 0, true);
  // log_source_b->SetVisAttributes(G4Color::Red());
  G4Tubs* source_cov = new G4Tubs("source", 0, source_rmin, 0.013*um, 0, 360*deg);
  G4LogicalVolume* log_source_cov = new G4LogicalVolume(source_cov, source_front_mat, "log_source_cov", 0,0,0); //gold
  //new G4PVPlacement(0, G4ThreeVector(0,0,-27.474*mm), log_source_cov,"source_cov_phys", logical_world, false, 0, true);
  // log_source_cov->SetVisAttributes(G4Color::Yellow());

  /*G4Tubs* source_mnt = new G4Tubs("source_mnt", 5.0*mm, source_rmax, 70*mm, 0, 360*deg);
  G4LogicalVolume* log_source_mnt = new G4LogicalVolume(source_mnt, stainless, "log_source_mnt", 0,0,0);
  new G4PVPlacement(rotX90, G4ThreeVector(0,100*mm,0), log_source_mnt,"source_mnt_phys", logical_world, false, 0, true);
  log_source_mnt->SetVisAttributes(G4Color::White()); */

  G4double innerRadiusOfTheSphere = 3.*mm;
  G4double outerRadiusOfTheSphere = 3.5*mm;
  G4Sphere* Sphere = new G4Sphere("Sphere_solid", innerRadiusOfTheSphere, outerRadiusOfTheSphere, 0.0, 360*deg, 0.0, 180*deg);
  G4LogicalVolume* logical_Sphere = new G4LogicalVolume(Sphere, vacuum, "Sphere_log", 0, 0, 0);
  if(source) new G4PVPlacement(0, G4ThreeVector(0,0,-29*mm), logical_Sphere, "Sphere_phys", logical_world, false, 0, true);
  else new G4PVPlacement(0, G4ThreeVector(0,0,0), logical_Sphere, "Sphere_phys", logical_world, false, 0, true);
  logical_Sphere->SetVisAttributes(G4Color::Yellow());

  if(chamber)
  {
    //JENSA Chamber
    G4double chmbInnerR = (315.3/2)*mm;
    G4double chmbOuterR = (323.8/2)*mm;
    G4double chmbLength = (150*mm); //(534.0/2)*mm;

    G4double beampipeInnerR = (94.0/2)*mm;
    G4double beampipeOuterR = (102.0/2)*mm;
    G4double beampipeLength = (40.23/2)*mm;
    
    //Create a horizontal and vertical tube
    G4Tubs* chamber_h = new G4Tubs("chamber_h", chmbInnerR , chmbOuterR, chmbLength, 0, 360*deg);
    G4Tubs* chamber_v = new G4Tubs("chamber_v", chmbInnerR , chmbOuterR, chmbLength*1.3, 0, 360*deg);

    //Create th perpenticular cut in shape of another tube that goes through it
    G4Tubs* chamberCut_v = new G4Tubs("chamberCut_v", 0, chmbOuterR, chmbLength*1.1, 0, 360*deg);
    
    //Create the front and back cut for the beam pipe
    G4Tubs* beapPipeCut = new G4Tubs("beapPipeCut", 0, beampipeOuterR, beampipeLength, 0, 360*deg);
    
    //Do the cuts on the horizontal chamber
    G4SubtractionSolid* chamberCut1 = new G4SubtractionSolid("Subtraction1", chamber_h, chamberCut_v, rotX90, G4ThreeVector());
    G4SubtractionSolid* chamberCut3 = new G4SubtractionSolid("Subtraction3", chamberCut1, beapPipeCut, rotY90, G4ThreeVector(chmbLength,0,0));
    G4SubtractionSolid* chamberCut5 = new G4SubtractionSolid("Subtraction5", chamberCut3, beapPipeCut, rotY90, G4ThreeVector(-chmbLength,0,0));
    //Do the cuts on the vertical chamber
    G4SubtractionSolid* chamberCut2 = new G4SubtractionSolid("Subtraction2", chamber_v, chamberCut_v, rotY90, G4ThreeVector());
    G4SubtractionSolid* chamberCut4 = new G4SubtractionSolid("Subtraction4", chamberCut2, beapPipeCut, rotX90, G4ThreeVector(0,-chmbLength,0));
    G4SubtractionSolid* chamberCut6 = new G4SubtractionSolid("Subtraction6", chamberCut4, beapPipeCut, rotX90, G4ThreeVector(0,chmbLength,0));
    
    //Create and place the cutted horizontal chamber
    G4LogicalVolume* log_chamber_h = new G4LogicalVolume(chamberCut5, stainless, "chamber_h_log");
    new G4PVPlacement(rotY90, G4ThreeVector(0.0,0.0,0.0), log_chamber_h, "chamber_h_phys", logical_world, false, 0, true);
    log_chamber_h->SetVisAttributes(G4Color::White());
    
    //Create and place the cutted vertical chamber
    G4LogicalVolume* log_chamber_v = new G4LogicalVolume(chamberCut6, stainless, "chamber_v_log");
    new G4PVPlacement(rotX90, G4ThreeVector(0.0,0.0,0.0), log_chamber_v, "chamber_v_phys", logical_world, false, 0, true);
    log_chamber_v->SetVisAttributes(G4Color::White());

    G4Cons* receiverCone = new G4Cons("receiverCone", 281.6/2, 309.1/2, 50.5/2, 59.7/2, 201.3/2, 0, 360*deg);
    G4LogicalVolume* log_receiverCone = new G4LogicalVolume(receiverCone, stainless, "receiverCone_log");
    new G4PVPlacement(rotX90, G4ThreeVector(0.0,-145.7,0.0), log_receiverCone, "receiverCone_phys", logical_world, false, 0, true);
    log_receiverCone->SetVisAttributes(G4Color::White());

    //Create and place the beam pipe downstream of the chamber
    G4Tubs* beamPipe_f = new G4Tubs("chamber_f", beampipeInnerR, beampipeOuterR, beampipeLength, 0, 360*deg);
    G4LogicalVolume* log_beamPipe_f = new G4LogicalVolume(beamPipe_f, stainless, "beamPipe_f_log");
    new G4PVPlacement(0, G4ThreeVector(0.0,0.0,chmbLength+8.0*mm+beampipeLength), log_beamPipe_f, "beamPipe_f_phys", logical_world, false, 0, true);
    log_beamPipe_f->SetVisAttributes(G4Color::White());
    
    //Create and place the beam pipe flange downstream of the chamber
    G4Tubs* beampipeFlange = new G4Tubs("beampipeFlange", beampipeOuterR*mm , (130/2)*mm, (11.3/2)*mm, 0, 360*deg);
    G4LogicalVolume* log_beampipeFlange = new G4LogicalVolume(beampipeFlange, stainless, "beampipeFlange_log");
    new G4PVPlacement(0, G4ThreeVector(0,0,chmbLength+8.0*mm+(beampipeLength*2)+(11.3/2)), log_beampipeFlange, "beampipeFlange_phys", logical_world, false, 0, true);
    log_beampipeFlange->SetVisAttributes(G4Color::White());

    //Create and place the flages on the left and right of the chamber
    G4Tubs* sideFlange = new G4Tubs("sideFlange", 0, chmbOuterR, chmbLength/20, 0, 360*deg);
    G4LogicalVolume* log_sideFlange = new G4LogicalVolume(sideFlange, stainless, "sideFlange_log");
    new G4PVPlacement(rotY90, G4ThreeVector(chmbLength+chmbLength/20,0.0,0.0), log_sideFlange, "sideFlange_L_phys", logical_world, false, 0, true);
    new G4PVPlacement(rotY90, G4ThreeVector(-chmbLength-chmbLength/20,0.0,0.0), log_sideFlange, "sideFlange_R_phys", logical_world, false, 1, true);
    log_sideFlange->SetVisAttributes(G4Color::White());
  }

  //PIPs detectors inside JENSA
  if(Si_monitor){

    G4Tubs* solidSiMon = new G4Tubs("solid_SiMon", 0, 6.91*mm, 0.5*mm, 0, 360*deg); 
    G4LogicalVolume* log_SiMon = new G4LogicalVolume(solidSiMon, silicon, "SiMon_log");
    new G4PVPlacement(rotY310, G4ThreeVector(68.101*mm,0.0,57.14*mm), log_SiMon, "SiMon_L_phys", logical_world, false, 0, true);
    new G4PVPlacement(rotY50, G4ThreeVector(-68.101*mm,0.0,57.14*mm), log_SiMon, "SiMon_R_phys", logical_world, false, 1, true);
    log_SiMon->SetVisAttributes(G4Color::Yellow());

  }

  // CAD geometry of the JENSA chamber
  // To use uncomment the CadMesh.hh at the top of the file
  // When used it creates a warning that some particles get stuck in the volume
  /*
  auto chamber = CADMesh::TessellatedMesh::FromSTL("../cadFiles/chamber.stl");
  auto logical_chamber = new G4LogicalVolume( chamber->GetSolid(), stainless, "chamber_log", 0, 0, 0 );
  new G4PVPlacement(0, G4ThreeVector(-6.2*cm,-9.8*cm,11.0*cm), logical_chamber, "chamber_phys", logical_world, false, 0, true);
  */

  // Detectors
  // To activate a detector change the respective flag from the DetectorConstruction.hh
  // Liquid Scintillators
  if(LScin)
  {
    //Get the positions of the neutron detectors
    NeutronDetectorsMap();

    G4Tubs* solidLScin = new G4Tubs("solid_LScin", 0, 25.4*mm, 25.4*mm, 0, 360*deg); 
    G4LogicalVolume* logical_LScin = new G4LogicalVolume(solidLScin,EJ301,"logic_LScin");
    G4VisAttributes vis_LScin(G4Color::Blue());
    vis_LScin.SetForceSolid(true);
    logical_LScin->SetVisAttributes(vis_LScin);

    G4Tubs* coverLScin = new G4Tubs("cover_LScin", 25.4*mm, 26.95*mm, 26.95*mm, 0, 360*deg); 
    G4LogicalVolume* logical_coverLScin = new G4LogicalVolume(coverLScin,aluminum,"logic_coverLScin");
    G4Tubs* coverfLScin = new G4Tubs("coverf_LScin", 0, 26.95*mm, 0.76*mm, 0, 360*deg);
    G4LogicalVolume* logical_coverfLScin = new G4LogicalVolume(coverfLScin,aluminum,"logic_coverfLScin");
    G4Tubs* coverbLScin = new G4Tubs("coverb_LScin", 0, 26.95*mm, 3.15*mm, 0, 360*deg);
    G4LogicalVolume* logical_coverbLScin = new G4LogicalVolume(coverbLScin,pyrex,"logic_coverbLScin");
    G4VisAttributes vis_coverLScin(G4Color::Grey());
    vis_coverLScin.SetForceSolid(true);
    logical_coverLScin->SetVisAttributes(vis_coverLScin);
    logical_coverfLScin->SetVisAttributes(vis_coverLScin);
    G4VisAttributes vis_coverbLScin(G4Color::Yellow());
    vis_coverbLScin.SetForceSolid(true);
    logical_coverbLScin->SetVisAttributes(vis_coverbLScin);
    
    for(G4int j=0; j<LScinNum; j++)
    {
      //double d=pi()/180;
      G4ThreeVector position = G4ThreeVector(LScinXpos[j]*mm,LScinYpos[j]*mm,LScinZpos[j]*mm);
      G4RotationMatrix* rotation = new G4RotationMatrix();
      rotation->rotateX(LScinXangle[j]*deg); rotation->rotateY(LScinYangle[j]*deg); rotation->rotateZ(LScinZangle[j]*deg);
      fLScin_phys.push_back(new G4PVPlacement(rotation,position,logical_LScin,LScinName[j],logical_world,false,j,true));
      if(coverLS){
        coverLScin_phys.push_back(new G4PVPlacement(rotation,position,logical_coverLScin,"coverLScin_phys",logical_world,false,j,true));
        G4ThreeVector position_covb = G4ThreeVector(LScinCovBXpos[j]*mm,LScinCovBYpos[j]*mm,LScinZpos[j]*mm);
        //G4ThreeVector position_covb = G4ThreeVector(LScinCovBXpos[j]*mm,LScinCovBYpos[j]*mm,-30.1*mm);
        coverbLScin_phys.push_back(new G4PVPlacement(rotation,position_covb,logical_coverbLScin,"coverbLScin_phys",logical_world,false,j,true));
        G4ThreeVector position_covf = G4ThreeVector(LScinCovFXpos[j]*mm,LScinCovFYpos[j]*mm,LScinZpos[j]*mm);
        //G4ThreeVector position_covf = G4ThreeVector(LScinCovFXpos[j]*mm,LScinCovFYpos[j]*mm,27.71*mm);
        coverfLScin_phys.push_back(new G4PVPlacement(rotation,position_covf,logical_coverfLScin,"coverfLScin_phys",logical_world,false,j,true));
      }
    }
  }
  
  //MCP - Upstream & Downstream
  if(MCP)
  {
    G4double dist_to_MCP = 600.0*mm;
    G4Tubs* MCPs = new G4Tubs("MCPs", 0, 3.5*cm, (0.50/2.0)*um, 0, CLHEP::twopi);
    G4LogicalVolume* logical_MCPs = new G4LogicalVolume(MCPs, mylar, "MCPs_log");
    new G4PVPlacement(0, G4ThreeVector(0,0,dist_to_MCP), logical_MCPs ,"MCPs_phys", logical_world, false, 0, true);

    G4VisAttributes vis_MCPs(G4Colour(0, 200, 255));
    vis_MCPs.SetForceSolid(true);
    logical_MCPs->SetVisAttributes(vis_MCPs);
  }

  //IC - dE Ionization Chamber
  G4double dist_to_IC = 700.0*mm;
  if(IC)
  {
    G4Tubs* IC_chamber = new G4Tubs("IC_chamber", 5.2*cm, 5.5*cm, (67.0/2)*mm, 0, CLHEP::twopi); 
    G4LogicalVolume* logical_IC_chamber = new G4LogicalVolume(IC_chamber, aluminum, "IC_chamber_log");
    new G4PVPlacement(0, G4ThreeVector(0,0,dist_to_IC+33.5), logical_IC_chamber ,"IC_chamber_phys", logical_world, false, 0, true);
  
    //IC window - Mylar foil
    G4double window_thick = 4.0*um;
    G4Tubs* IC_window = new G4Tubs("IC window", 0, 5.0*cm, (window_thick/2), 0, CLHEP::twopi);
    G4LogicalVolume* logical_IC_window = new G4LogicalVolume(IC_window, mylar, "IC_window_log");
    new G4PVPlacement(0, G4ThreeVector(0,0,dist_to_IC-0.1), logical_IC_window ,"IC_window_phys", logical_world, false, 0, true);

    //IC gas_1 - CF4
    G4Tubs* IC_gas_one = new G4Tubs("IC_gas_1", 0, 5.0*cm, (gas_thick_one/2)*mm, 0, CLHEP::twopi); 
    G4LogicalVolume* logical_IC_gas_one = new G4LogicalVolume(IC_gas_one, cf4, "IC_gas_1_log");
    new G4PVPlacement(0, G4ThreeVector(0,0,(dist_to_IC+gas_thick_one/2)*mm), logical_IC_gas_one ,"IC_gas_1_phys", logical_world, false, 0, true);

    //IC gas_dE - CF4
    G4Tubs* IC_gas_dE = new G4Tubs("IC_gas_dE", 0, 5.0*cm, (gas_thick_dE/2)*mm, 0, CLHEP::twopi); 
    G4LogicalVolume* logical_IC_gas_dE = new G4LogicalVolume(IC_gas_dE, cf4, "IC_gas_dE_log");
    new G4PVPlacement(0, G4ThreeVector(0,0,(dist_to_IC+gas_thick_one+gas_thick_dE/2)*mm), logical_IC_gas_dE ,"IC_gas_dE_phys", logical_world, false, 0, true);

    //IC gas_2 - CF4
    G4Tubs* IC_gas_two = new G4Tubs("IC_gas_2", 0, 5.0*cm, (gas_thick_two/2)*mm, 0, CLHEP::twopi); 
    G4LogicalVolume* logical_IC_gas_two = new G4LogicalVolume(IC_gas_two, cf4, "IC_gas_2_log");
    new G4PVPlacement(0, G4ThreeVector(0,0,(dist_to_IC+gas_thick_one+gas_thick_dE+gas_thick_two/2)*mm), logical_IC_gas_two ,"IC_gas_2_phys", logical_world, false, 0, true);

    G4VisAttributes vis_IC_chamber(G4Colour(192, 192, 192));
    vis_IC_chamber.SetForceSolid(true);
    logical_IC_chamber-> SetVisAttributes(vis_IC_chamber);

    G4VisAttributes vis_ICw(G4Colour(200, 0, 0));
    vis_ICw.SetForceSolid(true);
    logical_IC_window-> SetVisAttributes(vis_ICw);

    G4VisAttributes vis_IC_gas_one(G4Colour(100, 10, 0));
    vis_IC_gas_one.SetForceSolid(true);
    logical_IC_gas_one -> SetVisAttributes(vis_IC_gas_one);

    G4VisAttributes vis_IC_dE(G4Colour(250, 0, 127));
    vis_IC_dE.SetForceSolid(true);
    logical_IC_gas_dE -> SetVisAttributes(vis_IC_dE);

    G4VisAttributes vis_IC_gas_two(G4Colour(100, 10, 0));
    vis_IC_gas_two.SetForceSolid(true);
    logical_IC_gas_two -> SetVisAttributes(vis_IC_gas_two);
  }

  //DSSD Charged particle E detector
  if(DSSD)
  {
    G4double Dlayer_x = 3.20*cm; //3.2*cm; 
    G4double Dlayer_y = 3.20*cm; //3.2*cm;
    G4double Dlayer_z = 0.5*um; 

    G4Box* Dead_Layer = new G4Box("Dead_Layer",Dlayer_x,Dlayer_y,(Dlayer_z/2));
    G4LogicalVolume* logical_Dead_Layer = new G4LogicalVolume(Dead_Layer, aluminum, "Dead_Layer_log", 0,0,0);
    new G4PVPlacement(0, G4ThreeVector(0,0,(dist_to_IC+gas_thick_one+gas_thick_dE+gas_thick_two+Dlayer_z/2)*mm),logical_Dead_Layer,"Dead_Layer_phys",logical_world,false,0,true);

    G4double Det_x = 3.20*cm; //3.2*cm; 
    G4double Det_y = 3.20*cm; //3.2*cm; 
    G4double Det_z = 2.0*mm; 

    G4Box* Detector = new G4Box("Detector",Det_x,Det_y,(Det_z/2.0));
    G4LogicalVolume* logical_Detector = new G4LogicalVolume(Detector, silicon, "Detector_log", 0,0,0);
    new G4PVPlacement(0, G4ThreeVector(0,0,(dist_to_IC+gas_thick_one+gas_thick_dE+gas_thick_two+Dlayer_z+Det_z/2)*mm),logical_Detector,"Detector_phys",logical_world,false,0,true);

    G4VisAttributes vis_Dead_Layer(G4Colour(100, 0, 150));
    vis_Dead_Layer.SetForceSolid(true);
    logical_Dead_Layer -> SetVisAttributes(vis_Dead_Layer);

    G4VisAttributes vis_Detector(G4Colour(0, 0, 255));
    vis_Detector.SetForceSolid(true);
    logical_Detector -> SetVisAttributes(vis_Detector);
  }
  
  return physical_world; 
}

void DetectorConstruction::ConstructSDandField()
{
  G4SDManager::GetSDMpointer()->SetVerboseLevel(2);

  G4SDManager* SDman = G4SDManager::GetSDMpointer();
  
  if(DSSD)
  {
   auto SD1 = new SensitiveDetector("DSSD", "DetectorHitsCollection_1", analysis);
   SDman->AddNewDetector(SD1);
   SetSensitiveDetector("Detector_log", SD1);
  }
  
  if(IC)
  {
    auto SD2 = new SensitiveDetector_IC("IC", "DetectorHitsCollection_2", analysis);
    SDman->AddNewDetector(SD2);
    SetSensitiveDetector("IC_gas_dE_log", SD2);
  }
}

void DetectorConstruction::NeutronDetectorsMap()
{
  if(LScin)
  {
    LScinMap.open("../detectorMap/LScin_12.dat",std::ios::in);
    if (LScinMap.is_open()) G4cout<<"The LScinMap external file is open!"<<G4endl;
    while (true) 
    {
      ++LScinNum;
      //G4cout<<LScinNum<<G4endl;
      LScinMap >> name >> xpos >> ypos >> zpos >> xangle >> yangle >> zangle >> xcovfpos >> ycovfpos >> xcovbpos >> ycovbpos;
      LScinName.push_back(name);
      LScinXpos.push_back(xpos);
      LScinYpos.push_back(ypos);
      LScinZpos.push_back(zpos);
      LScinXangle.push_back(xangle);
      LScinYangle.push_back(yangle);
      LScinZangle.push_back(zangle);
      LScinCovFXpos.push_back(xcovfpos);
      LScinCovFYpos.push_back(ycovfpos);
      LScinCovBXpos.push_back(xcovbpos);
      LScinCovBYpos.push_back(ycovbpos);
      std::getline(LScinMap, line);

      if (LScinMap.eof()){break;}
    }
    LScinMap.close();
    analysis->SetNeutronDetectors(LScinNum);
  }
}

