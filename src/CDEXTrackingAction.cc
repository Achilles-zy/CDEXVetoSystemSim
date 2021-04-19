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

	G4double Ekin = trk->GetKineticEnergy();
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
	
	ParticleType = GetParticleIntType(ParticleName);

	G4bool ifFastDeposition = false;
	if (ParticleName == "e-" || ParticleName == "e+" || ParticleName == "alpha") {
		ifFastDeposition = true;
	}

	G4int ID = trk->GetTrackID();
	G4String CreatorProcessName = "Default";
	if (ID > 1) {
		CreatorProcessName = trk->GetCreatorProcess()->GetProcessName();
	}

	G4int ParentTrackID = trk->GetParentID();
	G4int TrackID = trk->GetTrackID();
	//if (trk->GetTrackID() == 1) {
	//	G4cout << ParticleName << G4endl;
	//	G4cout << TrackID << G4endl;
	//	G4cout << ParentTrackID << G4endl;
	//	G4cout << trk->GetCreatorProcess()->GetProcessName() << G4endl;

	//	G4cout << G4endl;
	//}

	//if (ParticleName == "gamma") {
	//	G4cout << ParticleName << G4endl;
	//	G4cout << TrackID << G4endl;
	//	G4cout << ParentTrackID << G4endl;
	//	G4cout << trk->GetCreatorProcess()->GetProcessName() << G4endl;

	//	G4cout << G4endl;
	//}
	
	G4int CreatorProcessType = -1;
	CreatorProcessType = GetCreatorProcessIntType(CreatorProcessName);

	//if (volume && logicvolume != detectorConstruction->GetLogicBulk() && logicvolume != detectorConstruction->GetLogicBEGe() && EdepTrack > 1 * eV && ifFastDeposition == true) {
	//	CDEXEvent->RecordStepInfo(ParticleType, CreatorProcessType, TrackPos.getX(), TrackPos.getY(), TrackPos.getZ(), EdepTrack);
	//	//G4cout << EdepTrack << G4endl;
	//	//G4cout << "Recorded" << G4endl;
	//}

	//G4String Mode = CDEXCons->GetMode();
	//if (volume && logicvolume == detectorConstruction->GetArgonVolume(Mode) && EdepTrack > 1 * eV && ParticleType != 0) {
	//	CDEXEvent->DetectableTrue();
	//}

	//if (volume && logicvolume == detectorConstruction->GetArgonVolume(Mode) && EdepTrack > 1 * eV && ifFastDeposition == true) {
	//	//CDEXEvent->DetectableTrue();
	//	//Record Gamma Photon-electric Process
	//	CDEXEvent->RecordStepInfoInScintillator(ParticleType, CreatorProcessType, TrackPos.getX(), TrackPos.getY(), TrackPos.getZ(), EdepTrack);
	//	//G4cout << "Recorded" << G4endl;
	//}

}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

void CDEXTrackingAction::AddEdepTrack(G4double edepstep) {
    EdepTrack += edepstep;
}

void CDEXTrackingAction::RecordTrackPos(G4ThreeVector trackpos) {
    TrackPos = trackpos;
}