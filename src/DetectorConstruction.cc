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
/// \file DetectorConstruction.cc
/// \brief Implementation of the DetectorConstruction class
//
//

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......
//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

#include "DetectorConstruction.hh"
#include "DetectorMessenger.hh"
#include "G4Material.hh"
#include "G4NistManager.hh"

#include "G4Box.hh"
#include "G4Tubs.hh"
#include "G4LogicalVolume.hh"
#include "G4PVPlacement.hh"

#include "G4GeometryManager.hh"
#include "G4PhysicalVolumeStore.hh"
#include "G4LogicalVolumeStore.hh"
#include "G4SolidStore.hh"
#include "G4RunManager.hh"

#include "G4UnitsTable.hh"
#include "G4SystemOfUnits.hh"

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

DetectorConstruction::DetectorConstruction()
:G4VUserDetectorConstruction(),
 fPBox(0), fLBox(0), fMaterial(0), fDetectorMessenger(0)
{
  fBoxSize = 1*m;
  DefineMaterials();
  SetMaterial("air");  
  fDetectorMessenger = new DetectorMessenger(this);
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

DetectorConstruction::~DetectorConstruction()
{ delete fDetectorMessenger;}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

G4VPhysicalVolume* DetectorConstruction::Construct()
{
  return ConstructVolumes();
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

void DetectorConstruction::DefineMaterials()
{

  // specific element name for thermal neutronHP
  // (see G4ParticleHPThermalScatteringNames.cc)

  G4int ncomponents, natoms;
  G4double massfraction;

  // pressurized water
  G4Element* H  = new G4Element("TS_H_of_Water" ,"H" , 1., 1.0079*g/mole);
  G4Element* O  = new G4Element("Oxygen"        ,"O" , 8., 16.00*g/mole);
  G4Material* H2O = 
  new G4Material("Water_ts", 1.000*g/cm3, ncomponents=2,
                         kStateLiquid, 593*kelvin, 150*bar);
  H2O->AddElement(H, natoms=2);
  H2O->AddElement(O, natoms=1);
  H2O->GetIonisation()->SetMeanExcitationEnergy(78.0*eV);
  
  // heavy water
  G4Isotope* H2 = new G4Isotope("H2",1,2);
  G4Element* D  = new G4Element("TS_D_of_Heavy_Water", "D", 1);
  D->AddIsotope(H2, 100*perCent);  
  G4Material* D2O = new G4Material("HeavyWater", 1.11*g/cm3, ncomponents=2,
                        kStateLiquid, 293.15*kelvin, 1*atmosphere);
  D2O->AddElement(D, natoms=2);
  D2O->AddElement(O, natoms=1);
  
  // graphite
  G4Isotope* C12 = new G4Isotope("C12", 6, 12);  
  G4Element* C   = new G4Element("TS_C_of_Graphite","C", ncomponents=1);
  C->AddIsotope(C12, 100.*perCent);
  G4Material* graphite = 
  new G4Material("graphite", 2.27*g/cm3, ncomponents=1,
                         kStateSolid, 293*kelvin, 1*atmosphere);
  graphite->AddElement(C, natoms=1);    
  
  // air
  G4Element* N = new G4Element("Nitrogen", "N", 7., 14.01*g/mole);
  G4Material* Air = new G4Material("air", 1.290*mg/cm3, ncomponents=2, kStateGas, 293*kelvin, 1*atmosphere);
  Air->AddElement(N, massfraction=70.*perCent);
  Air->AddElement(O, massfraction=30.*perCent);

  // iron
  G4Isotope* Fe56 = new G4Isotope("Fe56", 26, 56);
  G4Element* Fe = new G4Element("TS_Iron_Metal", "Fe", ncomponents=1);
  Fe->AddIsotope(Fe56, 100.*perCent);
  G4Material* iron = new G4Material("iron", 7.874*g/cm3, ncomponents=1, kStateSolid, 293*kelvin, 1*atmosphere);
  iron->AddElement(Fe, natoms=1);

  // polyethilene
  G4Element* Hpe = new G4Element("TS_H_of_Polyethylene", "H", 1, 1.0079*g/mole);
  G4Element* Cpe = new G4Element("Carbon", "C", 6, 12.01*g/mole);
  G4Material* polyethylene = new G4Material("polyethylene", 0.93*g/cm3, ncomponents=2, kStateSolid, 293*kelvin, 1*atmosphere);
  polyethylene->AddElement(Hpe, natoms=4);
  polyethylene->AddElement(Cpe, natoms=2);

 ///G4cout << *(G4Material::GetMaterialTable()) << G4endl;
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

G4Material* DetectorConstruction::MaterialWithSingleIsotope( G4String name,
                           G4String symbol, G4double density, G4int Z, G4int A)
{
 // define a material from an isotope
 //
 G4int ncomponents;
 G4double abundance, massfraction;

 G4Isotope* isotope = new G4Isotope(symbol, Z, A);
 
 G4Element* element  = new G4Element(name, symbol, ncomponents=1);
 element->AddIsotope(isotope, abundance= 100.*perCent);
 
 G4Material* material = new G4Material(name, density, ncomponents=1);
 material->AddElement(element, massfraction=100.*perCent);

 return material;
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

G4VPhysicalVolume* DetectorConstruction::ConstructVolumes()
{
  // Cleanup old geometry
  G4GeometryManager::GetInstance()->OpenGeometry();
  G4PhysicalVolumeStore::GetInstance()->Clean();
  G4LogicalVolumeStore::GetInstance()->Clean();
  G4SolidStore::GetInstance()->Clean();

  // Get materials
  auto iron = G4Material::GetMaterial("iron");
  auto polyethylene = G4Material::GetMaterial("polyethylene");
  auto D2O = G4Material::GetMaterial("HeavyWater");
  
  // world
  G4Box*
  sBox = new G4Box("Container",                         //its name
                   fBoxSize/2,fBoxSize/2,fBoxSize/2);   //its dimensions

  fLBox = new G4LogicalVolume(sBox,                     //its shape
                             fMaterial,                 //its material
                             fMaterial->GetName());     //its name

  fPBox = new G4PVPlacement(0,                          //no rotation
                            G4ThreeVector(),            //at (0,0,0)
                            fLBox,                      //its logical volume
                            fMaterial->GetName(),       //its name
                            0,                          //its mother  volume
                            false,                      //no boolean operation
                            0);                         //copy number
  

  // target
  G4double RMinO = 1.25*cm;	//Outter min radius
  G4double RMaxO = 2.5*cm;	//Outter max radius
  G4double Lz = 15*cm;		//Length

  G4double RMinI = 0*cm;	//Inner min radius
  G4double RMaxI = 1.25*cm;	//Inner max radius

  G4Tubs* sTargOut = new G4Tubs("shell",
		  RMinO,RMaxO,Lz/2,0.,2*M_PI);

  G4LogicalVolume* LTargOut = new G4LogicalVolume(sTargOut,
		  iron,
		  "Shell");

  auto PTargOut = new G4PVPlacement(0,
		  G4ThreeVector(14.*cm,0*cm,0*cm),
		  LTargOut,
		  "Shell",
		  fLBox,
		  false,
		  0);

  G4Tubs* sTargIn = new G4Tubs("core",
		  RMinI,RMaxI,Lz/2,0.,2*M_PI);

  G4LogicalVolume* LTargIn = new G4LogicalVolume(sTargIn,
		  D2O,
		  "Core");

  auto PTargIn = new G4PVPlacement(0,
		  G4ThreeVector(0.*cm,0.*cm,0.*cm),
		  LTargIn,
		  "Core",
		  LTargOut,
		  false,
		  0);

  // detector
  G4double Dx = 1*cm;		//Detector depth
  G4double Dy = 10*cm;		//Detector lateral length
  G4double Dz = 10*cm;

  G4Box* sDet = new G4Box("det",
		  Dx/2,Dy/2,Dz/2);

  G4LogicalVolume* LDet = new G4LogicalVolume(sDet,
		  polyethylene,
		  "Det");

  auto PDet = new G4PVPlacement(0,
		  G4ThreeVector(23.5*cm,0.*cm,0.*cm),
		  LDet,
		  "Det",
		  fLBox,
		  false,
		  0);

  // shielding-colimator
  G4double B1 = 23*cm;
  G4double B2 = 15*cm;
  G4double Bh = 53*cm;
  G4double Bp = 7*cm;
  G4double Bp2 = 38*cm;

  G4RotationMatrix* zRot = new G4RotationMatrix;
  zRot->rotateZ(M_PI/2*rad);

  G4Box* sShL = new G4Box("sh_long",
		  B1/2,B2/2,Bh/2);
  G4Box* sShS = new G4Box("sh_short",
		  B2/2,B2/2,Bh/2);
  G4Box* sShT = new G4Box("sh_top",
		  B1/2,B1/2,B2/2);
  G4Box* sShP = new G4Box("sh_pinhole",
      B2/2,Bp/2,B1/2);

  G4LogicalVolume* LShL = new G4LogicalVolume(sShL,
		  polyethylene,
		  "ShLong");
  G4LogicalVolume* LShS = new G4LogicalVolume(sShS,
		  polyethylene,
		  "ShShort");
  G4LogicalVolume* LShT = new G4LogicalVolume(sShT,
		  polyethylene,
		  "ShTop");
  G4LogicalVolume* LShP = new G4LogicalVolume(sShP,
      polyethylene,
      "ShPinh");

  auto PShL1 = new G4PVPlacement(0,G4ThreeVector(16.5*cm,19*cm,0.*cm),LShL,"ShLong1",fLBox,false,0);
  auto PShL2 = new G4PVPlacement(0,G4ThreeVector(16.5*cm,-19*cm,0.*cm),LShL,"ShLong2",fLBox,false,0);
  auto PShL3 = new G4PVPlacement(zRot,G4ThreeVector(-2.5*cm,15*cm,0.*cm),LShL,"ShLong3",fLBox,false,0);
  auto PShL4 = new G4PVPlacement(zRot,G4ThreeVector(-2.5*cm,-15*cm,0.*cm),LShL,"ShLong4",fLBox,false,0);
  auto PShL5 = new G4PVPlacement(zRot,G4ThreeVector(35.5*cm,0.*cm,0.*cm),LShL,"ShLong5",fLBox,false,0);

  auto PShS1 = new G4PVPlacement(0,G4ThreeVector(35.5*cm,19*cm,0.*cm),LShS,"ShShort1",fLBox,false,0);
  auto PShS2 = new G4PVPlacement(0,G4ThreeVector(35.5*cm,-19*cm,0.*cm),LShS,"ShShort2",fLBox,false,0);

  auto PShT1 = new G4PVPlacement(0,G4ThreeVector(16.5*cm,0.*cm,19*cm),LShT,"ShTop1",fLBox,false,0);     // (Bh/2-7.5)
  auto PShT2 = new G4PVPlacement(0,G4ThreeVector(16.5*cm,0.*cm,-19*cm),LShT,"ShTop2",fLBox,false,0);    // -(Bh/2-7.5)

  auto PShP1 = new G4PVPlacement(0,G4ThreeVector(-2.5*cm,0.*cm,15*cm),LShP,"ShPinh1",fLBox,false,0);
  auto PShP2 = new G4PVPlacement(0,G4ThreeVector(-2.5*cm,0.*cm,-15*cm),LShP,"ShPinh2",fLBox,false,0);


  PrintParameters();
  
  //always return the root volume
  //
  return fPBox;
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

void DetectorConstruction::PrintParameters()
{
  G4cout << "\n The Box is " << G4BestUnit(fBoxSize,"Length")
         << " of " << fMaterial->GetName() 
         << "\n \n" << fMaterial << G4endl;
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

void DetectorConstruction::SetMaterial(G4String materialChoice)
{
  // search the material by its name
  G4Material* pttoMaterial =
     G4NistManager::Instance()->FindOrBuildMaterial(materialChoice);   
  
  if (pttoMaterial) { 
    if(fMaterial != pttoMaterial) {
      fMaterial = pttoMaterial;
      if(fLBox) { fLBox->SetMaterial(pttoMaterial); }
      G4RunManager::GetRunManager()->PhysicsHasBeenModified();
    }
  } else {
    G4cout << "\n--> warning from DetectorConstruction::SetMaterial : "
           << materialChoice << " not found" << G4endl;
  }              
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

void DetectorConstruction::SetSize(G4double value)
{
  fBoxSize = value;
  G4RunManager::GetRunManager()->ReinitializeGeometry();
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

