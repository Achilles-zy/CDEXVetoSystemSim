// CDEXSteppingAction.cc

#include "CDEXSteppingAction.hh"

#include "CDEXDetectorConstruction.hh"
#include "CDEXEventAction.hh"
#include "CDEXRunAction.hh"

#include "G4Track.hh"
#include "G4SteppingManager.hh"

#include "G4ThreeVector.hh"
#include "G4Step.hh"
#include "globals.hh"
#include "G4RunManager.hh"
#include "G4HadronicProcessType.hh"
#include "G4EmProcessSubType.hh"
#include "G4SystemOfUnits.hh"

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

CDEXSteppingAction::CDEXSteppingAction(
				     CDEXEventAction* evt,CDEXRunAction* run)
:CDEXEvent(evt),CDEXRun(run),EnableAcc(false)
{ }
//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

void CDEXSteppingAction::UserSteppingAction(const G4Step* aStep)
{
	EnableAcc = CDEXRun->GetAccelerate();
	//G4cout << "EnableAcc: " << CDEXRun->GetAccelerate() << G4endl;
	auto volume = aStep->GetPostStepPoint()->GetTouchableHandle()->GetVolume();
	auto touchable = aStep->GetPostStepPoint()->GetTouchableHandle();
	auto particle_name = aStep->GetTrack()->GetParticleDefinition()->GetParticleName();

	auto edep = aStep->GetTotalEnergyDeposit();
	G4int TrackID = aStep->GetTrack()->GetTrackID();
	//**********************For Acceleration**********************//
	if (EnableAcc == true) {
		if (CDEXEvent->GetAcceletateStatus() == false && TrackID % 200 == 0) {
			SignalSiPMCount = 0;
			ContainerSignalSiPMCount = 0;
			G4int rownb = CDEXEvent->GetRowNb();
			G4int columnnb = CDEXEvent->GetColumnNb();
			for (int i = 0; i < rownb; i++)
			{
				for (int j = 0; j < columnnb; j++)
				{
					if (CDEXEvent->GetSiPMSignalCount(i, j) > 0) {
						SignalSiPMCount++;
					}
				}
			}

			G4int crownb = CDEXEvent->GetContainerRowNb();
			G4int ccolumnnb = CDEXEvent->GetContainerColumnNb();
			for (int i = 0; i < crownb; i++)
			{
				for (int j = 0; j < ccolumnnb; j++)
				{
					if (CDEXEvent->GetContainerSiPMSignalCount(i, j) > 0) {
						ContainerSignalSiPMCount++;
					}
				}
			}
			G4int totalcnt = SignalSiPMCount + ContainerSignalSiPMCount;
			if (totalcnt >= 2) {
				CDEXEvent->SetAccelerate(true);
			}
		}

		if (CDEXEvent->GetAcceletateStatus() == true && particle_name == "opticalphoton") {
			aStep->GetTrack()->SetTrackStatus(fStopAndKill);
		}
		//if (CDEXEvent->GetTotalSiPMPhotonCnt() > 2 && particle_name == "opticalphoton") {
		//	aStep->GetTrack()->SetTrackStatus(fStopAndKill);
		//}
	}

	//***********************************************************//

	const CDEXDetectorConstruction* detectorConstruction
		= static_cast<const CDEXDetectorConstruction*>
		(G4RunManager::GetRunManager()->GetUserDetectorConstruction());

	if (volume == detectorConstruction->GetBulk() && particle_name != "opticalphoton") {
		//CDEXEvent->BulkTrue();
		CDEXEvent->AddBulkEnergy(edep);
	}
	
	//if (volume == detectorConstruction->GetEnv() && particle_name == "opticalphoton") {
	//	CDEXEvent->DetectableTrue();
	//}

	if (particle_name == "opticalphoton") {
		CDEXEvent->DetectableTrue();
	}

	//G4cout << aStep->GetPostStepPoint()->GetPosition() << G4endl;
	for (int i = 0; i < 4; i++) {
		if (volume == detectorConstruction->GetSiPM(i) && detectorConstruction->GetSiPM(i) != nullptr && particle_name == "opticalphoton" ) {
			aStep->GetTrack()->SetTrackStatus(fStopAndKill);
			CDEXEvent->DetectableTrue();
			G4int GetCopyNumber0 = touchable->GetCopyNumber(0);
			G4int GetCopyNumber1 = touchable->GetCopyNumber(1);
			G4int GetCopyNumber2 = touchable->GetCopyNumber(2);
			//G4cout << "CpNb0=" << GetCopyNumber0 << G4endl;
			//G4cout << "CpNb1=" << GetCopyNumber1 << G4endl;
			//G4cout << "CpNb2=" << GetCopyNumber2 << G4endl;
			//G4cout << "EndStep"<< G4endl;
			G4double Energy = aStep->GetPostStepPoint()->GetKineticEnergy() / (1 * eV);
			G4double WaveLength = 1242 / Energy;//nm
			G4double SiPMEff = GetEfficiency(WaveLength);
			//G4cout << "Wavelength = " << WaveLength << " nm" << G4endl;	
			//G4cout << "Efficiency = " << SiPMEff  << G4endl;
			G4double rnd = G4UniformRand();
			if (rnd < SiPMEff) {
				CDEXEvent->AddToSiPMSignal(GetCopyNumber1, GetCopyNumber0);
			}
			CDEXEvent->AddToSiPM(GetCopyNumber1, GetCopyNumber0);
			CDEXEvent->CountTotalSiPMPhoton(1);
		}
	}

	for (int j = 0; j < 5; j++) {
		if (volume == detectorConstruction->GetContainerSiPM(j) && detectorConstruction->GetContainerSiPM(j) != nullptr && particle_name == "opticalphoton") {
			aStep->GetTrack()->SetTrackStatus(fStopAndKill);

			G4int GetCopyNumber0 = touchable->GetCopyNumber(0);
			G4int GetCopyNumber1 = touchable->GetCopyNumber(1);
			G4double Energy = aStep->GetPostStepPoint()->GetKineticEnergy() / (1 * eV);
			G4double WaveLength = 1242 / Energy;//nm
			G4double SiPMEff = GetEfficiency(WaveLength);
			G4double rnd = G4UniformRand();
			if (rnd < SiPMEff) {
				CDEXEvent->AddToContainerSiPMSignal(GetCopyNumber1, GetCopyNumber0);
			}
			CDEXEvent->AddToContainerSiPM(GetCopyNumber1, GetCopyNumber0);
			CDEXEvent->CountTotalSiPMPhoton(1);
		}
	}

	//G4int processtype = aStep->GetPostStepPoint()->GetProcessDefinedStep()->GetProcessSubType();
	//G4int creatorprocess = aStep->GetTrack()->GetCreatorProcess()->GetProcessSubType();
	//G4int parentID = aStep->GetTrack()->GetParentID();

	//if (parentID == 1 ) {
	//	if (processtype == fScintillation || processtype == fRadioactiveDecay) {
	//		//aStep->GetTrack()->SetTrackStatus(fKillTrackAndSecondaries);
	//	}
	//}

	//if (processtype == fRadioactiveDecay) {
	//	//G4cout << "Parent ID =" <<parentID << G4endl;
	//}

	//if (parentID == 1) {
	//	//G4cout << "Parent ID =" << parentID << "Particle =" << aStep->GetTrack()->GetParticleDefinition()->GetParticleName() << "Process =" << aStep->GetPostStepPoint()->GetProcessDefinedStep()->GetProcessName() <<aStep->GetTotalEnergyDeposit()<< G4endl;
	//}
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

G4double CDEXSteppingAction::GetEfficiency(G4double wavelength) {
	G4double eff;
	if (wavelength > 400 && wavelength <= 900) {
		eff = 0.5 - (wavelength - 400) / 1000;
	}
	else if (wavelength > 100 && wavelength <= 400) {
		eff = (wavelength - 100) * 5 / 3000;
	}
	else {
		eff = 0;
	}
	return eff;
}