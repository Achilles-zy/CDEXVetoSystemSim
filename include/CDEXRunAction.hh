//CDEXRunAction.hh

#ifndef CDEXRunAction_h
#define CDEXRunAction_h

#include "G4UserRunAction.hh"
#include "globals.hh"
#include "G4AccumulableManager.hh"
#include "G4Accumulable.hh"
#include "CDEXPrimaryGeneratorAction.hh"
#include "CDEXDetectorConstruction.hh"

class G4Run;
class CDEXPrimaryGeneratorAction;
class CDEXDetectorConstruction;
class CDEXRunMessenger;

class CDEXRunAction : public G4UserRunAction
{
    public:
    CDEXRunAction(CDEXPrimaryGeneratorAction*, CDEXDetectorConstruction*);
    ~CDEXRunAction();

    void BeginOfRunAction(const G4Run*);
    void EndOfRunAction(const G4Run*);

	void CDEXOutput(const G4Run* aRun);
	void CDEXArrayOutput(const G4Run* aRun);

	void RefreshOutput(G4bool b) {
		ifRefresh = b;
	}
	void SetAccelerate(G4bool b) {
		ifAccelerate = b;
	}
	G4bool GetAccelerate() {
		return ifAccelerate;
	}
	void CountSiPMEvent() { 
		SiPMEventCount += 1;
	};


	inline 	void CountBulkEvent() {
		BulkEventCount += 1;
	};
	inline	void CountVetoEvent() {
		VetoEventCount += 1;
	}
	inline	void CountVetoEvent_0() {
		VetoEventCount_0 += 1;
	}
	inline	void CountVetoEvent_1() {
		VetoEventCount_1 += 1;
	}
	inline	void CountVetoEvent_2() {
		VetoEventCount_2 += 1;
	}
	inline	void CountVetoEvent_3() {
		VetoEventCount_3 += 1;
	}
	inline	void CountVetoEvent_4() {
		VetoEventCount_4 += 1;
	}
	inline 	void CountDetectableEvent() {
		DetectableEventCount += 1;
	}
	inline 	void CountVetoPossibleEvent() {
		VetoPossibleEvtCount += 1;
	}


	inline 	void CountROIVetoEvent() {
		VetoEventCount += 1;
	}
	inline 	void CountROIVetoPossibleEvent() {
		VetoPossibleEvtCount += 1;
	}
	inline 	void CountROIEvent() {
		ROIEventCount += 1;
	}

private:
	//number of events that generate signals in bulk
	//G4int BulkEventCount;
	//number of events that generate signals in SiPMs
	//G4int SiPMEventCount;
	//G4int VetoEventCount;
	CDEXPrimaryGeneratorAction* fPrimaryGenerator;
	CDEXDetectorConstruction* fDetCons;
	G4Accumulable<G4int> SiPMEventCount;

	G4Accumulable<G4int> VetoEventCount;
	G4Accumulable<G4int> VetoEventCount_0;
	G4Accumulable<G4int> VetoEventCount_1;
	G4Accumulable<G4int> VetoEventCount_2;
	G4Accumulable<G4int> VetoEventCount_3;
	G4Accumulable<G4int> VetoEventCount_4;
	G4Accumulable<G4int> BulkEventCount;
	G4Accumulable<G4int> DetectableEventCount;
	G4Accumulable<G4int> VetoPossibleEvtCount;

	G4Accumulable<G4int> ROIVetoEventCount;

	G4Accumulable<G4int> ROIEventCount;
	G4Accumulable<G4int> ROIVetoPossibleEvtCount;


	CDEXRunMessenger* fRunMessenger;
	G4String filename;
	G4String txtname;
	G4bool ifRefresh;
	G4bool ifAccelerate;
	G4int runID;
};

#endif