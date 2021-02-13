#include "CDEXPrimaryGeneratorAction.hh"

#include "G4Event.hh"
#include "G4GeneralParticleSource.hh"
#include "G4ParticleTable.hh"
#include "G4ParticleDefinition.hh"
#include "globals.hh"
#include "G4SystemOfUnits.hh"
#include "G4PhysicalVolumeStore.hh"
#include "G4RunManager.hh"

//#include "TMath.h"
#include "Randomize.hh"

CDEXPrimaryGeneratorAction::CDEXPrimaryGeneratorAction(CDEXDetectorConstruction* det):
	G4VUserPrimaryGeneratorAction(),
	PrimaryE(0),
	InitialE(1 * keV),
	PrimaryName(""),
	SrcType("Bucket"),
	ImprintID(1)
{
    fCDEXGPS = new G4GeneralParticleSource();
	fPrimaryMessenger = new CDEXPrimaryGeneratorMessenger(this);
	fDetCons = det;
	G4ParticleTable* particleTable = G4ParticleTable::GetParticleTable();
	G4String particleName = "e-";
	fCDEXGPS->SetParticleDefinition(particleTable->FindParticle(particleName));

	/*
	G4PhysicalVolumeStore* PVStore = G4PhysicalVolumeStore::GetInstance();
	G4int i = 0;
	G4VPhysicalVolume* tempPV = NULL;
	while (i < G4int(PVStore->size())) {
		tempPV = (*PVStore)[i];
		G4cout << i << " " << " " << tempPV->GetName() << " " << G4endl;
		i++;
	}
	*/
}

CDEXPrimaryGeneratorAction::~CDEXPrimaryGeneratorAction()
{
    delete fCDEXGPS;
}

void CDEXPrimaryGeneratorAction::GeneratePrimaries(G4Event* anEvent)
{
	
	G4ParticleTable* particleTable = G4ParticleTable::GetParticleTable();
	//G4String particleName = "e-";
	//G4double particleEnergy = 0 * MeV;
	G4String mode = fDetCons->GetMode();
	G4int EvtID = anEvent->GetEventID();
	if (mode == "Unit") {
		if (SrcType == "Wire") {
			G4double Radius = fDetCons->GetWireRadius();
			G4double Length = fDetCons->GetWireLength();
			G4ThreeVector WirePos = fDetCons->GetWirePos();

			//G4cout << "==========================Primary Info==========================" << G4endl;
			//G4cout << "Wire Position: " << WirePos << G4endl;
			//G4cout << "Wire Radius: " << Radius << G4endl;
			//G4cout << "Wire Length: " << Length << G4endl;
			//G4cout << "================================================================" << G4endl;

			fCDEXGPS->GetCurrentSource()->GetEneDist()->SetEnergyDisType("Mono");
			fCDEXGPS->GetCurrentSource()->GetAngDist()->SetAngDistType("iso");
			fCDEXGPS->GetCurrentSource()->GetPosDist()->SetPosDisType("Volume");
			fCDEXGPS->GetCurrentSource()->GetPosDist()->SetPosDisShape("Cylinder");
			fCDEXGPS->GetCurrentSource()->GetPosDist()->SetCentreCoords(WirePos);
			fCDEXGPS->GetCurrentSource()->GetPosDist()->SetRadius(Radius);
			fCDEXGPS->GetCurrentSource()->GetPosDist()->SetHalfZ(Length / 2);

			fCDEXGPS->GetCurrentSource()->GetPosDist()->ConfineSourceToVolume("Wire");
		}
		else if (SrcType == "PENShell") {
			G4double Radius = fDetCons->GetPENShellRadius();
			G4double Length = fDetCons->GetPENShellLength();

			if (G4RunManager::GetRunManager()->GetRunManagerType() == 1) {
				G4cout << "==========================Primary Info==========================" << G4endl;
				G4cout << "Sample Region Radius: " << Radius << G4endl;
				G4cout << "Sample Region Length: " << Length << G4endl;
				G4cout << "================================================================" << G4endl;
			}

			fCDEXGPS->GetCurrentSource()->GetEneDist()->SetEnergyDisType("Mono");
			fCDEXGPS->GetCurrentSource()->GetAngDist()->SetAngDistType("iso");
			fCDEXGPS->GetCurrentSource()->GetPosDist()->SetPosDisType("Volume");
			fCDEXGPS->GetCurrentSource()->GetPosDist()->SetPosDisShape("Cylinder");
			fCDEXGPS->GetCurrentSource()->GetPosDist()->SetRadius(Radius * 2);
			fCDEXGPS->GetCurrentSource()->GetPosDist()->SetHalfZ(Length / 2);

			fCDEXGPS->GetCurrentSource()->GetPosDist()->ConfineSourceToVolume("PENShell");
		}
		else if (SrcType == "InnerShell") {
			G4double Radius = fDetCons->GetPENShellRadius();
			G4double Length = fDetCons->GetPENShellLength();

			fCDEXGPS->GetCurrentSource()->GetEneDist()->SetEnergyDisType("Mono");
			fCDEXGPS->GetCurrentSource()->GetAngDist()->SetAngDistType("iso");
			fCDEXGPS->GetCurrentSource()->GetPosDist()->SetPosDisType("Volume");
			fCDEXGPS->GetCurrentSource()->GetPosDist()->SetPosDisShape("Cylinder");
			fCDEXGPS->GetCurrentSource()->GetPosDist()->SetRadius(Radius * 2);
			fCDEXGPS->GetCurrentSource()->GetPosDist()->SetHalfZ(Length / 2);

			fCDEXGPS->GetCurrentSource()->GetPosDist()->ConfineSourceToVolume("InnerShell");
		}
		else if (SrcType == "OuterShell") {
			G4double Radius = fDetCons->GetPENShellRadius();
			G4double Length = fDetCons->GetPENShellLength();

			fCDEXGPS->GetCurrentSource()->GetEneDist()->SetEnergyDisType("Mono");
			fCDEXGPS->GetCurrentSource()->GetAngDist()->SetAngDistType("iso");
			fCDEXGPS->GetCurrentSource()->GetPosDist()->SetPosDisType("Volume");
			fCDEXGPS->GetCurrentSource()->GetPosDist()->SetPosDisShape("Cylinder");
			fCDEXGPS->GetCurrentSource()->GetPosDist()->SetRadius(Radius * 2);
			fCDEXGPS->GetCurrentSource()->GetPosDist()->SetHalfZ(Length / 2);

			fCDEXGPS->GetCurrentSource()->GetPosDist()->ConfineSourceToVolume("OuterShell");
		}
		else if (SrcType == "Point") {
			fCDEXGPS->GetCurrentSource()->GetPosDist()->SetPosDisType("Point");
			fCDEXGPS->GetCurrentSource()->GetPosDist()->SetCentreCoords(G4ThreeVector(0, 0, 25 * mm));
			fCDEXGPS->GetCurrentSource()->GetEneDist()->SetEnergyDisType("Mono");
			fCDEXGPS->GetCurrentSource()->GetAngDist()->SetAngDistType("iso");
			fCDEXGPS->GetCurrentSource()->GetEneDist()->SetEmax(InitialE);
			fCDEXGPS->GetCurrentSource()->GetEneDist()->SetEmin(InitialE);
		}
		else if (SrcType == "InnerShellSurface") {
			fCDEXGPS->GetCurrentSource()->GetPosDist()->SetPosDisType("Surface");
			fCDEXGPS->GetCurrentSource()->GetEneDist()->SetEnergyDisType("Mono");
			fCDEXGPS->GetCurrentSource()->GetAngDist()->SetAngDistType("iso");
			fCDEXGPS->GetCurrentSource()->GetPosDist()->SetPosDisShape("Cylinder");
			fCDEXGPS->GetCurrentSource()->GetPosDist()->SetRadius(51 * mm);
			fCDEXGPS->GetCurrentSource()->GetPosDist()->SetHalfZ(51 * mm / 2);
		}
		else {
			G4cout << "Error: Src type not found! Using Geant4 default settings." << G4endl;
		}
	}
	else if (mode == "SArUnit") {
		if (SrcType == "Crystal") {
			G4double Radius = fDetCons->GetSArBrickRadius();
			G4double Length = fDetCons->GetSArBrickHeight();
			if (G4RunManager::GetRunManager()->GetRunManagerType() == 1) {
				G4cout << "==========================Primary Info==========================" << G4endl;
				G4cout << "Sample Region Radius: " << Radius << G4endl;
				G4cout << "Sample Region Length: " << Length << G4endl;
				G4cout << "================================================================" << G4endl;
			}

			fCDEXGPS->GetCurrentSource()->GetEneDist()->SetEnergyDisType("Mono");
			fCDEXGPS->GetCurrentSource()->GetAngDist()->SetAngDistType("iso");
			fCDEXGPS->GetCurrentSource()->GetPosDist()->SetPosDisType("Volume");
			fCDEXGPS->GetCurrentSource()->GetPosDist()->SetPosDisShape("Cylinder");
			fCDEXGPS->GetCurrentSource()->GetPosDist()->SetCentreCoords(G4ThreeVector(0, 0, Length / 2));
			fCDEXGPS->GetCurrentSource()->GetPosDist()->SetRadius(Radius);
			fCDEXGPS->GetCurrentSource()->GetPosDist()->SetHalfZ(Length / 2);

			fCDEXGPS->GetCurrentSource()->GetPosDist()->ConfineSourceToVolume("SArCrystal");
		}

		//Needs to be fixed
		
		else if (SrcType == "Container") {
			G4double Radius = fDetCons->GetSArBrickRadius();
			G4double Length = fDetCons->GetSArBrickHeight();
			if (G4RunManager::GetRunManager()->GetRunManagerType() == 1) {
				G4cout << "==========================Primary Info==========================" << G4endl;
				G4cout << "Sample Region Radius: " << Radius << G4endl;
				G4cout << "Sample Region Length: " << Length << G4endl;
				G4cout << "================================================================" << G4endl;
			}
			fCDEXGPS->GetCurrentSource()->GetEneDist()->SetEnergyDisType("Mono");
			fCDEXGPS->GetCurrentSource()->GetAngDist()->SetAngDistType("iso");
			fCDEXGPS->GetCurrentSource()->GetPosDist()->SetPosDisType("Volume");
			fCDEXGPS->GetCurrentSource()->GetPosDist()->SetPosDisShape("Cylinder");
			fCDEXGPS->GetCurrentSource()->GetPosDist()->SetCentreCoords(G4ThreeVector(0, 0, Length / 2));
			fCDEXGPS->GetCurrentSource()->GetPosDist()->SetRadius(Radius);
			fCDEXGPS->GetCurrentSource()->GetPosDist()->SetHalfZ(Length / 2);

			fCDEXGPS->GetCurrentSource()->GetPosDist()->ConfineSourceToVolume("SArBrick");
		}
		
		else {
			G4cout << "Error: Src type not found! Using Geant4 default settings." << G4endl;
		}
	}
	else if (mode == "Array-1") {
		if (SrcType == "SingleOuterShell") {
			G4ThreeVector CentCoord;
			G4double zCoord;
			if (ImprintID % 2 == 0) {
				zCoord = -(ImprintID / 2 - 1) * 65 - 32.5;
				CentCoord = G4ThreeVector(0, 0, zCoord * mm);
			}
			if (ImprintID % 2 == 1) {
				zCoord = (ImprintID - 1) / 2 * 65 + 32.5;
				CentCoord = G4ThreeVector(0, 0, zCoord * mm);
			}
			G4double Radius = fDetCons->GetPENShellRadius();
			G4double Length = fDetCons->GetPENShellLength();

			if (G4RunManager::GetRunManager()->GetRunManagerType() == 1) {
				G4cout << "==========================Primary Info==========================" << G4endl;
				G4cout << "Sample Region Radius: " << Radius << G4endl;
				G4cout << "Sample Region Length: " << Length << G4endl;
				G4cout << "================================================================" << G4endl;
			}

			fCDEXGPS->GetCurrentSource()->GetEneDist()->SetEnergyDisType("Mono");
			fCDEXGPS->GetCurrentSource()->GetAngDist()->SetAngDistType("iso");
			fCDEXGPS->GetCurrentSource()->GetPosDist()->SetPosDisType("Volume");
			fCDEXGPS->GetCurrentSource()->GetPosDist()->SetPosDisShape("Cylinder");
			fCDEXGPS->GetCurrentSource()->GetPosDist()->SetRadius(Radius * 2);
			fCDEXGPS->GetCurrentSource()->GetPosDist()->SetHalfZ(Length / 2);
			fCDEXGPS->GetCurrentSource()->GetPosDist()->SetCentreCoords(CentCoord);
			fCDEXGPS->GetCurrentSource()->GetPosDist()->ConfineSourceToVolume("av_1_impr_" + std::to_string(ImprintID) + "_logicOuterShell_pv_0");
		}

		if (SrcType == "SingleInnerShell") {
			G4ThreeVector CentCoord;
			G4double zCoord;
			if (ImprintID % 2 == 0) {
				zCoord = -(ImprintID / 2 - 1) * 65 - 32.5;
				CentCoord = G4ThreeVector(0, 0, zCoord * mm);
			}
			if (ImprintID % 2 == 1) {
				zCoord = (ImprintID - 1) / 2 * 65 + 32.5;
				CentCoord = G4ThreeVector(0, 0, zCoord * mm);
			}
			G4double Radius = fDetCons->GetPENShellRadius();
			G4double Length = fDetCons->GetPENShellLength();

			if (G4RunManager::GetRunManager()->GetRunManagerType() == 1) {
				G4cout << "==========================Primary Info==========================" << G4endl;
				G4cout << "Sample Region Radius: " << Radius << G4endl;
				G4cout << "Sample Region Length: " << Length << G4endl;
				G4cout << "================================================================" << G4endl;
			}

			fCDEXGPS->GetCurrentSource()->GetEneDist()->SetEnergyDisType("Mono");
			fCDEXGPS->GetCurrentSource()->GetAngDist()->SetAngDistType("iso");
			fCDEXGPS->GetCurrentSource()->GetPosDist()->SetPosDisType("Volume");
			fCDEXGPS->GetCurrentSource()->GetPosDist()->SetPosDisShape("Cylinder");
			fCDEXGPS->GetCurrentSource()->GetPosDist()->SetRadius(Radius);
			fCDEXGPS->GetCurrentSource()->GetPosDist()->SetHalfZ(Length / 2);
			fCDEXGPS->GetCurrentSource()->GetPosDist()->SetCentreCoords(CentCoord);
			fCDEXGPS->GetCurrentSource()->GetPosDist()->ConfineSourceToVolume("av_1_impr_" + std::to_string(ImprintID) + "_logicInnerShell_pv_1");
		}

		else if (SrcType == "StringBoxCrystal") {
			G4double Radius = fDetCons->GetPENShellRadius();
			G4double Length = fDetCons->GetPENShellLength() * 22;

			if (G4RunManager::GetRunManager()->GetRunManagerType() == 1) {
				G4cout << "==========================Primary Info==========================" << G4endl;
				G4cout << "Sample Region Radius: " << Radius << G4endl;
				G4cout << "Sample Region Length: " << Length << G4endl;
				G4cout << "================================================================" << G4endl;
			}

			fCDEXGPS->GetCurrentSource()->GetEneDist()->SetEnergyDisType("Mono");
			fCDEXGPS->GetCurrentSource()->GetAngDist()->SetAngDistType("iso");
			fCDEXGPS->GetCurrentSource()->GetPosDist()->SetPosDisType("Volume");
			fCDEXGPS->GetCurrentSource()->GetPosDist()->SetPosDisShape("Cylinder");
			fCDEXGPS->GetCurrentSource()->GetPosDist()->SetRadius(Radius * 2);
			fCDEXGPS->GetCurrentSource()->GetPosDist()->SetHalfZ(Length / 2);
			fCDEXGPS->GetCurrentSource()->GetPosDist()->ConfineSourceToVolume("StringBoxCrystal");
		}
		else if (SrcType == "Wire") {
		G4double Radius = fDetCons->GetWireRadius();
		G4double Length = fDetCons->GetPENShellLength() * 22;
		G4ThreeVector WirePos = fDetCons->GetWirePos();
		if (G4RunManager::GetRunManager()->GetRunManagerType() == 1) {
			G4cout << "==========================Primary Info==========================" << G4endl;
			G4cout << "Sample Region Radius: " << Radius << G4endl;
			G4cout << "Sample Region Length: " << Length << G4endl;
			G4cout << "================================================================" << G4endl;
		}
		fCDEXGPS->GetCurrentSource()->GetEneDist()->SetEnergyDisType("Mono");
		fCDEXGPS->GetCurrentSource()->GetAngDist()->SetAngDistType("iso");
		fCDEXGPS->GetCurrentSource()->GetPosDist()->SetPosDisType("Volume");
		fCDEXGPS->GetCurrentSource()->GetPosDist()->SetPosDisShape("Cylinder");
		fCDEXGPS->GetCurrentSource()->GetPosDist()->SetRadius(Radius * 2);
		fCDEXGPS->GetCurrentSource()->GetPosDist()->SetHalfZ(Length / 2);
		fCDEXGPS->GetCurrentSource()->GetPosDist()->SetCentreCoords(WirePos);
		fCDEXGPS->GetCurrentSource()->GetPosDist()->ConfineSourceToVolume("Wire");
		}

		else if (SrcType == "SingleWire") {
		G4ThreeVector WirePosIni = fDetCons->GetWirePos();
		G4ThreeVector CentCoord;
		G4double zCoord;
		if (ImprintID % 2 == 0) {
			zCoord = -(ImprintID / 2 - 1) * 65 - 32.5;
			CentCoord = G4ThreeVector(0, 0, zCoord * mm) + WirePosIni;
		}
		if (ImprintID % 2 == 1) {
			zCoord = (ImprintID - 1) / 2 * 65 + 32.5;
			CentCoord = G4ThreeVector(0, 0, zCoord * mm) + WirePosIni;
		}
		G4double Radius = fDetCons->GetWireRadius();
		G4double Length = fDetCons->GetPENShellLength();

		if (G4RunManager::GetRunManager()->GetRunManagerType() == 1) {
			G4cout << "==========================Primary Info==========================" << G4endl;
			G4cout << "Sample Region Radius: " << Radius << G4endl;
			G4cout << "Sample Region Length: " << Length << G4endl;
			G4cout << "================================================================" << G4endl;
		}

		fCDEXGPS->GetCurrentSource()->GetEneDist()->SetEnergyDisType("Mono");
		fCDEXGPS->GetCurrentSource()->GetAngDist()->SetAngDistType("iso");
		fCDEXGPS->GetCurrentSource()->GetPosDist()->SetPosDisType("Volume");
		fCDEXGPS->GetCurrentSource()->GetPosDist()->SetPosDisShape("Cylinder");
		fCDEXGPS->GetCurrentSource()->GetPosDist()->SetRadius(Radius * 2);
		fCDEXGPS->GetCurrentSource()->GetPosDist()->SetHalfZ(Length / 2);
		fCDEXGPS->GetCurrentSource()->GetPosDist()->SetCentreCoords(CentCoord);
		fCDEXGPS->GetCurrentSource()->GetPosDist()->ConfineSourceToVolume("Wire");
		}

		else if (SrcType == "SingleASIC") {
		G4ThreeVector CentCoord;
		G4double zCoord;
		if (ImprintID % 2 == 0) {
			zCoord = -(ImprintID / 2 - 1) * 65 - 32.5 - 25 * mm;
			CentCoord = G4ThreeVector(0, 0, zCoord * mm);
		}
		if (ImprintID % 2 == 1) {
			zCoord = (ImprintID - 1) / 2 * 65 + 32.5 - 25 * mm;
			CentCoord = G4ThreeVector(0, 0, zCoord * mm);
		}
		G4double Radius = 3 * cm;
		G4double Length = fDetCons->GetASICThickness();

		if (G4RunManager::GetRunManager()->GetRunManagerType() == 1) {
			G4cout << "==========================Primary Info==========================" << G4endl;
			G4cout << "Sample Region Radius: " << Radius << G4endl;
			G4cout << "Sample Region Length: " << Length << G4endl;
			G4cout << "================================================================" << G4endl;
		}

		fCDEXGPS->GetCurrentSource()->GetEneDist()->SetEnergyDisType("Mono");
		fCDEXGPS->GetCurrentSource()->GetAngDist()->SetAngDistType("iso");
		fCDEXGPS->GetCurrentSource()->GetPosDist()->SetPosDisType("Volume");
		fCDEXGPS->GetCurrentSource()->GetPosDist()->SetPosDisShape("Cylinder");
		fCDEXGPS->GetCurrentSource()->GetPosDist()->SetRadius(Radius);
		fCDEXGPS->GetCurrentSource()->GetPosDist()->SetHalfZ(Length);
		fCDEXGPS->GetCurrentSource()->GetPosDist()->SetCentreCoords(CentCoord);
		fCDEXGPS->GetCurrentSource()->GetPosDist()->ConfineSourceToVolume("ASIC");
		}
		else {
			G4cout << "Error: Src type not found! Using Geant4 default settings." << G4endl;
		}
	}
	else if (mode == "CDEXBucket") {
		if (SrcType == "Bucket") {
			G4ThreeVector CentCoord;
			G4double zCoord;
			G4double BucketHeight = fDetCons->GetBucketHeight();
			G4double BucketThickness = fDetCons->GetBucketThickness();
			G4double BucketRadius = fDetCons->GetBucketOuterRadius();

			if (G4RunManager::GetRunManager()->GetRunManagerType() == 1) {
				G4cout << "==========================Primary Info==========================" << G4endl;
				G4cout << "Sample Region Radius: " << BucketRadius << G4endl;
				G4cout << "Sample Region Height: " << BucketHeight << G4endl;
				G4cout << "================================================================" << G4endl;
			}
			
			fCDEXGPS->GetCurrentSource()->GetEneDist()->SetEnergyDisType("Mono");
			fCDEXGPS->GetCurrentSource()->GetAngDist()->SetAngDistType("iso");
			fCDEXGPS->GetCurrentSource()->GetPosDist()->SetPosDisType("Volume");
			fCDEXGPS->GetCurrentSource()->GetPosDist()->SetPosDisShape("Cylinder");
			fCDEXGPS->GetCurrentSource()->GetPosDist()->SetRadius(BucketThickness * 1.1);
			fCDEXGPS->GetCurrentSource()->GetPosDist()->SetHalfZ(BucketHeight / 2);
			fCDEXGPS->GetCurrentSource()->GetPosDist()->ConfineSourceToVolume("Bucket");
		}

		else if (SrcType == "Wire") {
			G4double Radius = fDetCons->GetWireRadius();
			G4double Length = fDetCons->GetSmallestUnitHeight() * fDetCons->GetUnitNumber();
			G4ThreeVector WirePos = fDetCons->GetWirePos();
			if (G4RunManager::GetRunManager()->GetRunManagerType() == 1) {
				G4cout << "==========================Primary Info==========================" << G4endl;
				G4cout << "Sample Region Radius: " << Radius << G4endl;
				G4cout << "Sample Region Length: " << Length << G4endl;
				G4cout << "================================================================" << G4endl;
			}
			fCDEXGPS->GetCurrentSource()->GetEneDist()->SetEnergyDisType("Mono");
			fCDEXGPS->GetCurrentSource()->GetAngDist()->SetAngDistType("iso");
			fCDEXGPS->GetCurrentSource()->GetPosDist()->SetPosDisType("Volume");
			fCDEXGPS->GetCurrentSource()->GetPosDist()->SetPosDisShape("Cylinder");
			fCDEXGPS->GetCurrentSource()->GetPosDist()->SetRadius(Radius);
			fCDEXGPS->GetCurrentSource()->GetPosDist()->SetHalfZ(Length / 2);
			fCDEXGPS->GetCurrentSource()->GetPosDist()->SetCentreCoords(WirePos);
			fCDEXGPS->GetCurrentSource()->GetPosDist()->ConfineSourceToVolume("Wire");
		}

		else if (SrcType == "SingleASIC") {
			G4ThreeVector CentCoord;
			G4double zCoord;
			if (ImprintID % 2 == 0) {
				zCoord = -(ImprintID / 2 - 1) * fDetCons->GetSmallestUnitHeight() - fDetCons->GetSmallestUnitHeight()/2;
				CentCoord = G4ThreeVector(0, 0, zCoord * mm);
			}
			if (ImprintID % 2 == 1) {
				zCoord = (ImprintID - 1) / 2 * fDetCons->GetSmallestUnitHeight() + fDetCons->GetSmallestUnitHeight() / 2;
				CentCoord = G4ThreeVector(0, 0, zCoord * mm);
			}
			G4double Radius = 3 * cm;
			G4double Length = fDetCons->GetSmallestUnitHeight();

			if (G4RunManager::GetRunManager()->GetRunManagerType() == 1) {
				G4cout << "==========================Primary Info==========================" << G4endl;
				G4cout << "Sample Region Radius: " << Radius << G4endl;
				G4cout << "Sample Region Length: " << Length << G4endl;
				G4cout << "================================================================" << G4endl;
			}

			fCDEXGPS->GetCurrentSource()->GetEneDist()->SetEnergyDisType("Mono");
			fCDEXGPS->GetCurrentSource()->GetAngDist()->SetAngDistType("iso");
			fCDEXGPS->GetCurrentSource()->GetPosDist()->SetPosDisType("Volume");
			fCDEXGPS->GetCurrentSource()->GetPosDist()->SetPosDisShape("Cylinder");
			fCDEXGPS->GetCurrentSource()->GetPosDist()->SetRadius(Radius);
			fCDEXGPS->GetCurrentSource()->GetPosDist()->SetHalfZ(Length/2);
			fCDEXGPS->GetCurrentSource()->GetPosDist()->SetCentreCoords(CentCoord);
			fCDEXGPS->GetCurrentSource()->GetPosDist()->ConfineSourceToVolume("ASIC");
		}
		else {
			G4cout << "Error: Src type not found! Using Geant4 default settings." << G4endl;
		}
	}
	else {
		G4cout << "Error: Mode not found! Using Geant4 default settings." << G4endl;
	}


    fCDEXGPS->GeneratePrimaryVertex(anEvent);


	if (anEvent->GetEventID() == 0) {
		PrimaryE = fCDEXGPS->GetCurrentSource()->GetParticleEnergy();
		PrimaryName = fCDEXGPS->GetCurrentSource()->GetParticleDefinition()->GetParticleName();
	}
}