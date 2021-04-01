//CDEXTrackingAction.cc

#include "CDEXTrackingAction.hh"
#include "CDEXDetectorConstruction.hh"
#include "CDEXEventAction.hh"

#include "G4Event.hh"
#include "G4Track.hh"
#include "G4EventManager.hh"
#include "G4TrajectoryContainer.hh"
#include "G4Trajectory.hh"
#include "G4VVisManager.hh"
#include "G4ios.hh"
#include "G4RunManager.hh"
#include "G4SystemOfUnits.hh"

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......
 
CDEXTrackingAction::CDEXTrackingAction(CDEXEventAction* evt, CDEXDetectorConstruction* cons)
    :CDEXEvent(evt),CDEXCons(cons)
{
    EdepTrack = 0;
    TrackPos = G4ThreeVector(0, 0, 0);
}
//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......
 
CDEXTrackingAction::~CDEXTrackingAction()
{}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......
 
void CDEXTrackingAction::PreUserTrackingAction(const G4Track* track)
{
    EdepTrack = 0;
    TrackPos = G4ThreeVector(0, 0, 0);
    
    G4double trackTime = track->GetGlobalTime();
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......
 
void CDEXTrackingAction::PostUserTrackingAction(const G4Track* trk)
{
    const CDEXDetectorConstruction* detectorConstruction
        = static_cast<const CDEXDetectorConstruction*>
        (G4RunManager::GetRunManager()->GetUserDetectorConstruction());

	auto volume = trk->GetVolume();
	G4LogicalVolume* logicvolume;
	if (volume) {
		logicvolume = volume->GetLogicalVolume();
		//G4cout << volume << volume->GetName() << G4endl;
		//G4cout << logicvolume << logicvolume->GetName() << G4endl;
		//G4cout << G4endl;
	}

	G4String ParticleName = trk->GetParticleDefinition()->GetParticleName();
	//G4cout << ParticleName << G4endl;
	G4int ParticleType = -1;
	G4bool ifFastDeposition = false;
	if (ParticleName == "opticalphoton") {
		ParticleType = 0;
	}
	else if (ParticleName == "gamma") {
		ParticleType = 1;
	}
	else if (ParticleName == "e-") {
		ParticleType = 2;
		ifFastDeposition = true;
	}
	else if (ParticleName == "e+") {
		ParticleType = 3;
		ifFastDeposition = true;
	}
	else if (ParticleName == "alpha") {
		ParticleType = 4;
		ifFastDeposition = true;
	}
	else {
		ParticleType = 5;
	}
	G4int ID = trk->GetTrackID();
	G4String CreatorProcessName = "Default";
	if (ID > 1) {
		CreatorProcessName = trk->GetCreatorProcess()->GetProcessName();
	}

	G4int CreatorProcessType = -1;
	if (CreatorProcessName == "RadioactiveDecay") {
		CreatorProcessType = 0;
	}
	else if (CreatorProcessName == "conv") {
		CreatorProcessType = 1;
	}
	else if (CreatorProcessName == "phot") {
		CreatorProcessType = 2;
	}
	else if (CreatorProcessName == "compt") {
		CreatorProcessType = 3;
	}
	else {
		CreatorProcessType = 4;
	}
	//Record Gamma Photon-electric Process

	if (volume && logicvolume != detectorConstruction->GetLogicBulk() && logicvolume != detectorConstruction->GetLogicBEGe() && EdepTrack > 1 * eV && ifFastDeposition == true) {
		CDEXEvent->RecordStepInfo(ParticleType, CreatorProcessType, TrackPos.getX(), TrackPos.getY(), TrackPos.getZ(), EdepTrack);
	}

	G4String Mode = CDEXCons->GetMode();
	if (volume && logicvolume == detectorConstruction->GetArgonVolume(Mode) && EdepTrack > 1 * eV && ParticleType != 0) {
		CDEXEvent->DetectableTrue();
		if (ifFastDeposition == true) {
			//Record Gamma Photon-electric Process
			CDEXEvent->RecordStepInfoInScintillator(ParticleType, CreatorProcessType, TrackPos.getX(), TrackPos.getY(), TrackPos.getZ(), EdepTrack);
		}
	}
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

void CDEXTrackingAction::AddEdepTrack(G4double edepstep) {
    EdepTrack += edepstep;
}

void CDEXTrackingAction::RecordTrackPos(G4ThreeVector trackpos) {
    TrackPos = trackpos;
}