//CDEXRunAction.cc
 
#include "CDEXRunAction.hh"
#include "CDEXRunMessenger.hh"
#include "G4Run.hh"
#include "G4RunManager.hh"
#include "G4SystemOfUnits.hh"
#include <time.h>
#include "g4root.hh"

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......
 
CDEXRunAction::CDEXRunAction(CDEXPrimaryGeneratorAction* gen, CDEXDetectorConstruction* det):
	SiPMEventCount(0),
	BulkEventCount(0),
	VetoEventCount(0),
	VetoEventCount_0(0),
	VetoEventCount_1(0),
	VetoEventCount_2(0),
	VetoEventCount_3(0),
	VetoEventCount_4(0),
	DetectableEventCount(0),
	VetoPossibleEvtCount(0),
	ROIEventCount(0),
	ROIVetoEventCount(0),
	ROIVetoPossibleEvtCount(0),
	ifRefresh(false),
	ifAccelerate(true),
	runID(0)
{
  fRunMessenger = new CDEXRunMessenger(this);
  auto analysisManager = G4AnalysisManager::Instance();
  analysisManager->SetNtupleMerging(true);
  analysisManager->SetVerboseLevel(1);
  analysisManager->CreateH1("edepBulk", "Edep in Bulk", 200, 0 * keV, 20 * keV);
  analysisManager->CreateH1("PhotonCount", "Photon Count", 1000, 0, 3000);
  analysisManager->CreateH1("PhotonCount2", "Photon Count", 100, 0, 100);
  analysisManager->CreateH1("VetoPhotonCount", "Photon Count", 100, 0, 3000);
  //analysisManager->CreateH1("EdepBulk1", "Edep in Bulk", 200, 0 * keV, 2500 * keV);
  //analysisManager->CreateH1("EdepBulk2", "Edep in Bulk", 200, 2000 * keV, 2100 * keV);

  analysisManager->CreateNtuple("RunRes", "EventCount");
  analysisManager->CreateNtupleIColumn(0, "VetoCount");
  analysisManager->CreateNtupleIColumn(0, "SiPMCount");
  analysisManager->CreateNtupleIColumn(0, "BulkCount");
  analysisManager->CreateNtupleIColumn(0, "DetectableCount");
  analysisManager->CreateNtupleIColumn(0, "VetoPossibleCount");
  analysisManager->CreateNtupleIColumn(0, "ROIEventCount");
  analysisManager->CreateNtupleIColumn(0, "ROIVetoEventCount");
  analysisManager->CreateNtupleIColumn(0, "ROIVetoPossibleEvtCount");
  analysisManager->FinishNtuple(0);

  analysisManager->CreateNtuple("PhotonCount", "PhotonCoutnt");
  for (int i = 0; i < 30; i++) {
	  G4String Ntuplename = "SiPMPhotonCount" + std::to_string(i);
	  analysisManager->CreateNtupleIColumn(1, Ntuplename);
  }
  analysisManager->FinishNtuple(1);

  analysisManager->CreateNtuple("EdepBulk", "Edep in Bulk");
  analysisManager->CreateNtupleDColumn(2, "Edep");
  analysisManager->FinishNtuple(2);

  G4AccumulableManager* accumulableManager = G4AccumulableManager::Instance();
  accumulableManager->RegisterAccumulable(SiPMEventCount);
  accumulableManager->RegisterAccumulable(VetoEventCount);
  accumulableManager->RegisterAccumulable(VetoEventCount_0);
  accumulableManager->RegisterAccumulable(VetoEventCount_1);
  accumulableManager->RegisterAccumulable(VetoEventCount_2);
  accumulableManager->RegisterAccumulable(VetoEventCount_3);
  accumulableManager->RegisterAccumulable(VetoEventCount_4);
  accumulableManager->RegisterAccumulable(BulkEventCount);
  accumulableManager->RegisterAccumulable(DetectableEventCount);
  accumulableManager->RegisterAccumulable(VetoPossibleEvtCount);
  accumulableManager->RegisterAccumulable(ROIEventCount);
  accumulableManager->RegisterAccumulable(ROIVetoEventCount);
  accumulableManager->RegisterAccumulable(ROIVetoPossibleEvtCount);
  fDetCons = det;
  fPrimaryGenerator = gen;
  filename = "Simulation Result";
  txtname = "Simulation Result";

}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......
 
CDEXRunAction::~CDEXRunAction()
{
  delete G4AnalysisManager::Instance();
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......
 
void CDEXRunAction::BeginOfRunAction(const G4Run* aRun)
{
  auto analysisManager = G4AnalysisManager::Instance();

  G4int RunID = aRun->GetRunID();
  G4String fileName;
  if (fDetCons->GetMode() == "Unit") {
	  if (fPrimaryGenerator->GetSrcType() == "PENShell") {
		  //fileName = "Unit_" + fDetCons->GetConfine() + "_" + std::to_string(RunID);
		  fileName = "Unit_" + std::to_string(fDetCons->GetPENPropertiesID()) + "_" + fPrimaryGenerator->GetSrcType() + "_" + fDetCons->GetReflectorType();
		  filename = fileName;
		  txtname = "Unit_" + std::to_string(fDetCons->GetPENPropertiesID()) + "_" + fPrimaryGenerator->GetSrcType() + "_" + fDetCons->GetReflectorType();
	  }
	  else if (fPrimaryGenerator->GetSrcType() == "InnerShellSurface") {
		  //fileName = "Unit_" + fDetCons->GetConfine() + "_" + std::to_string(RunID);
		  fileName = "Unit_" + std::to_string(fDetCons->GetPENPropertiesID()) + "_" + fPrimaryGenerator->GetSrcType() + "_" + fDetCons->GetReflectorType();
		  filename = fileName;
		  txtname = "Unit_" + std::to_string(fDetCons->GetPENPropertiesID()) + "_" + fPrimaryGenerator->GetSrcType() + "_" + fDetCons->GetReflectorType();
	  }
	  else if (fPrimaryGenerator->GetSrcType() == "Wire") {
		  //fileName = "Unit_" + fDetCons->GetConfine() + "_" + fDetCons->GetWireType() + "_" + std::to_string(RunID);
		  fileName = "Unit_" + std::to_string(fDetCons->GetPENPropertiesID()) + "_" + fPrimaryGenerator->GetSrcType() + "_" + fDetCons->GetReflectorType() + "_" + fDetCons->GetWireType();
		  filename = fileName;
		  txtname = "Unit_" + std::to_string(fDetCons->GetPENPropertiesID()) + "_" + fPrimaryGenerator->GetSrcType() + "_" + fDetCons->GetReflectorType() + fDetCons->GetWireType();
	  }
  }
  else if (fDetCons->GetMode() == "Array-1") {
	  if (fPrimaryGenerator->GetSrcType() == "PENShell") {
		  //fileName = "Unit_" + fDetCons->GetConfine() + "_" + std::to_string(RunID);
		  fileName = "Array-1_" + std::to_string(fDetCons->GetPENPropertiesID()) + "_" + fPrimaryGenerator->GetSrcType() + "_" + fDetCons->GetReflectorType();
		  filename = fileName;
		  txtname = "Array-1_" + std::to_string(fDetCons->GetPENPropertiesID()) + "_" + fPrimaryGenerator->GetSrcType() + "_" + fDetCons->GetReflectorType();
	  }
	  else if (fPrimaryGenerator->GetSrcType() == "InnerShellSurface") {
		  //fileName = "Unit_" + fDetCons->GetConfine() + "_" + std::to_string(RunID);
		  fileName = "Array-1_" + std::to_string(fDetCons->GetPENPropertiesID()) + "_" + fPrimaryGenerator->GetSrcType() + "_" + fDetCons->GetReflectorType();
		  filename = fileName;
		  txtname = "Array-1_" + std::to_string(fDetCons->GetPENPropertiesID()) + "_" + fPrimaryGenerator->GetSrcType() + "_" + fDetCons->GetReflectorType();
	  }
	  else if (fPrimaryGenerator->GetSrcType() == "Wire") {
		  //fileName = "Unit_" + fDetCons->GetConfine() + "_" + fDetCons->GetWireType() + "_" + std::to_string(RunID);
		  fileName = "Array-1_" + std::to_string(fDetCons->GetPENPropertiesID()) + "_" + fPrimaryGenerator->GetSrcType() + "_" + fDetCons->GetReflectorType() + "_" + fDetCons->GetWireType();
		  filename = fileName;
		  txtname = "Array-1_" + std::to_string(fDetCons->GetPENPropertiesID()) + "_" + fPrimaryGenerator->GetSrcType() + "_" + fDetCons->GetReflectorType() + fDetCons->GetWireType();
	  }
  }
  else if (fDetCons->GetMode() == "CDEXBucket") {
	  if (fPrimaryGenerator->GetSrcType() == "Bucket") {
		  //fileName = "Unit_" + fDetCons->GetConfine() + "_" + std::to_string(RunID);
		  fileName = "CDEXBucket_" + fPrimaryGenerator->GetSrcType();
		  filename = fileName;
		  txtname = "CDEXBucket_" + fPrimaryGenerator->GetSrcType();
	  }
	  else if (fPrimaryGenerator->GetSrcType() == "Wire") {
		  //fileName = "Unit_" + fDetCons->GetConfine() + "_" + std::to_string(RunID);
		  fileName = "CDEXBucket_"+ fPrimaryGenerator->GetSrcType();
		  filename = fileName;
		  txtname = "CDEXBucket_"+ fPrimaryGenerator->GetSrcType();
	  }
	  else if (fPrimaryGenerator->GetSrcType() == "SingleASIC") {
		  //fileName = "Unit_" + fDetCons->GetConfine() + "_" + fDetCons->GetWireType() + "_" + std::to_string(RunID);
		  fileName = "CDEXBucket_" + fPrimaryGenerator->GetSrcType() + "_" + fDetCons->GetWireType();
		  filename = fileName;
		  txtname = "CDEXBucket_" + fPrimaryGenerator->GetSrcType() + "_" + fDetCons->GetWireType();
	  }
  }

  //analysisManager->OpenFile(fileName);
  analysisManager->OpenFile(filename);

  G4AccumulableManager* accumulableManager = G4AccumulableManager::Instance();
  accumulableManager->Reset();
  if (G4RunManager::GetRunManager()->GetRunManagerType() == 1) {
	  G4cout << "Run started." << G4endl;
  }
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......
 
void CDEXRunAction::EndOfRunAction(const G4Run* aRun)
{
  G4AccumulableManager* accumulableManager = G4AccumulableManager::Instance();
  accumulableManager->Merge();

  auto analysisManager = G4AnalysisManager::Instance();
  if (G4RunManager::GetRunManager()->GetRunManagerType() != 1) {
	  analysisManager->FillNtupleIColumn(0, 0, VetoEventCount.GetValue());
	  analysisManager->FillNtupleIColumn(0, 1, SiPMEventCount.GetValue());
	  analysisManager->FillNtupleIColumn(0, 2, BulkEventCount.GetValue());
	  analysisManager->FillNtupleIColumn(0, 3, DetectableEventCount.GetValue());
	  analysisManager->FillNtupleIColumn(0, 4, VetoPossibleEvtCount.GetValue());
	  analysisManager->FillNtupleIColumn(0, 5, ROIEventCount.GetValue());
	  analysisManager->FillNtupleIColumn(0, 6, ROIVetoEventCount.GetValue());
	  analysisManager->FillNtupleIColumn(0, 7, ROIVetoPossibleEvtCount.GetValue());
	  analysisManager->AddNtupleRow(0);
  }

  if (fDetCons->GetMode() == "Unit") {
	  CDEXOutput(aRun);
  }
  else if (fDetCons->GetMode() == "Array-1") {
	  CDEXArrayOutput(aRun);
  }
  else if (fDetCons->GetMode() == "SArUnit") {
	  CDEXOutput(aRun);
  }
  else if (fDetCons->GetMode() == "CDEXBucket") {
	  CDEXOutput(aRun);
  }
  else {
	  G4cout << "ERRRO! Mode does not exsist, nothing to output!" << G4endl;
  }

  ifRefresh = false;
  G4String fileName0 = filename + ".root";
  G4String fileName1 = filename + "_" + fPrimaryGenerator->GetPrimaryName()+ ".root";
  //analysisManager->SetFileName(filename);
  analysisManager->Write();
  analysisManager->CloseFile();
  if (G4RunManager::GetRunManager()->GetRunManagerType() == 1) {
	  std::rename(fileName0, fileName1);
  }
  G4cout << "End" << G4endl;
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

void CDEXRunAction::CDEXOutput(const G4Run* aRun) {
	if (G4RunManager::GetRunManager()->GetRunManagerType() == 1) {
		G4cout << "End Running ..." << G4endl;
		G4cout << "==========================Run Summary==========================" << G4endl;
		G4cout << G4endl;
		G4cout << G4endl;
		G4cout << "VetoEventCount =" << VetoEventCount.GetValue() << G4endl;
		G4cout << "VetoEventCount_1 =" << VetoEventCount_1.GetValue() << G4endl;
		G4cout << "VetoEventCount_2 =" << VetoEventCount_2.GetValue() << G4endl;
		G4cout << "VetoEventCount_3 =" << VetoEventCount_3.GetValue() << G4endl;
		G4cout << "VetoEventCount_4 =" << VetoEventCount_4.GetValue() << G4endl;
		G4cout << "SiPMEventCount =" << SiPMEventCount.GetValue() << G4endl;
		G4cout << "BulkEventCount =" << BulkEventCount.GetValue() << G4endl;
		G4cout << "DetectableEventCount =" << DetectableEventCount.GetValue() << G4endl;
		G4cout << "VetoPossibleEventCount =" << VetoPossibleEvtCount.GetValue() << G4endl;
		G4cout << "ROIEventCount =" << ROIEventCount.GetValue() << G4endl;
		G4cout << "ROIVetoEventCount =" << ROIVetoEventCount.GetValue() << G4endl;
		G4cout << "ROIVetoPossibleEventCount =" << ROIVetoPossibleEvtCount.GetValue() << G4endl;
		G4cout << G4endl;
		G4cout << G4endl;
		G4cout << "===============================================================" << G4endl;
	
		std::ofstream output;
		if (aRun->GetRunID() == 0) {
			output.open(txtname + ".txt", std::ios::ate);
			if (fPrimaryGenerator->GetSrcType() == "PENShell") {
				output 
					<< "Source Distribution:\t"<< "PENShell" << G4endl;
			}
			else if (fPrimaryGenerator->GetSrcType() == "Wire") {
				output
					<< "Source Distribution:\t" << "Wire" << G4endl
					<< "Wire Type:\t" << fDetCons->GetWireType() << G4endl
				    << "Reflector Type:\t" << fDetCons->GetReflectorType() << G4endl;
			}
			output
				//<< "Confine Info:\t" << fDetCons->GetConfine() << G4endl
				<< "Simulation Result:" << G4endl;
		}
		else
		{
			output.open(txtname + ".txt", std::ios::app);
		}
		G4double eff;
		if (BulkEventCount.GetValue() != 0) {
			eff = VetoEventCount.GetValue() / BulkEventCount.GetValue();
		}
		else
		{
			eff = 0;
		}

		G4double ROIeff;
		if (ROIEventCount.GetValue() != 0) {
			ROIeff = ROIVetoEventCount.GetValue() / ROIEventCount.GetValue();
		}
		else
		{
			ROIeff = 0;
		}

		if (aRun->GetRunID() == 0) {
			output
				<< std::setw(10) << std::left << "Run ID" << '\t'
				<< std::setw(40) << std::left << "Number of Event" << '\t'
				<< std::setw(40) << std::left << "Primary Particle" << '\t'
				<< std::setw(40) << std::left << "Primary Energy(MeV)" << '\t'
				<< std::setw(40) << std::left << "Reflector Type" << '\t'
				<< std::setw(40) << std::left << "SiPMEvent" << '\t'
				<< std::setw(40) << std::left << "VetoEvent" << '\t'
				<< std::setw(40) << std::left << "VetoEvent_1" << '\t'
				<< std::setw(40) << std::left << "VetoEvent_2" << '\t'
				<< std::setw(40) << std::left << "VetoEvent_3" << '\t'
				<< std::setw(40) << std::left << "VetoEvent_4" << '\t'
				<< std::setw(40) << std::left << "BulkEvent" << '\t'
				<< std::setw(40) << std::left << "DetectableEvent" << '\t'
				<< std::setw(40) << std::left << "VetoPossibleEvent" << '\t'
				<< std::setw(40) << std::left << "VetoEfficiency" << '\t'
				<< std::setw(40) << std::left << "ROIEvent" << '\t'
				<< std::setw(40) << std::left << "ROIVetoEvent" << '\t'
				<< std::setw(40) << std::left << "ROIVetoPossibleEvent" << '\t'
				<< std::setw(40) << std::left << "ROIVetoEfficiency" << '\t' << G4endl;
		}

		output
			<< std::setw(10) << std::left << aRun->GetRunID() << '\t'
			<< std::setw(40) << std::left << aRun->GetNumberOfEvent() << '\t'
			<< std::setw(40) << std::left << fPrimaryGenerator->GetPrimaryName() << '\t'
			<< std::setw(40) << std::left << std::setiosflags(std::ios::fixed) << std::setprecision(2) << fPrimaryGenerator->GetPrimaryE() << '\t'
			<< std::setw(40) << std::left << fDetCons->GetReflectorType() << '\t'
			<< std::setw(40) << std::left << SiPMEventCount.GetValue() << '\t'
			<< std::setw(40) << std::left << VetoEventCount.GetValue() << '\t'
			<< std::setw(40) << std::left << VetoEventCount_1.GetValue() << '\t'
			<< std::setw(40) << std::left << VetoEventCount_2.GetValue() << '\t'
			<< std::setw(40) << std::left << VetoEventCount_3.GetValue() << '\t'
			<< std::setw(40) << std::left << VetoEventCount_4.GetValue() << '\t'
			<< std::setw(40) << std::left << BulkEventCount.GetValue() << '\t'
			<< std::setw(40) << std::left << DetectableEventCount.GetValue() << '\t'
			<< std::setw(40) << std::left << VetoPossibleEvtCount.GetValue() << '\t'
			<< std::setw(40) << std::left << eff << '\t'
			<< std::setw(40) << std::left << ROIEventCount.GetValue() << '\t'
			<< std::setw(40) << std::left << ROIVetoEventCount.GetValue() << '\t'
			<< std::setw(40) << std::left << ROIVetoPossibleEvtCount.GetValue() << '\t'
			<< std::setw(40) << std::left << ROIeff << G4endl;
		output.close();
		//std::DecimalFormat df1 = new DecimalFormat("0.0");
	}
}

void CDEXRunAction::CDEXArrayOutput(const G4Run* aRun) {
	if (G4RunManager::GetRunManager()->GetRunManagerType() == 1) {
		G4cout << "End Running ..." << G4endl;
		G4cout << "==========================Run Summary==========================" << G4endl;
		G4cout << G4endl;
		G4cout << G4endl;
		G4cout << "VetoEventCount =" << VetoEventCount.GetValue() << G4endl;
		G4cout << "VetoEventCount_0 =" << VetoEventCount_0.GetValue() << G4endl;
		G4cout << "VetoEventCount_1 =" << VetoEventCount_1.GetValue() << G4endl;
		G4cout << "VetoEventCount_2 =" << VetoEventCount_2.GetValue() << G4endl;
		G4cout << "VetoEventCount_3 =" << VetoEventCount_3.GetValue() << G4endl;
		G4cout << "VetoEventCount_4 =" << VetoEventCount_4.GetValue() << G4endl;
		G4cout << "SiPMEventCount =" << SiPMEventCount.GetValue() << G4endl;
		G4cout << "BulkEventCount =" << BulkEventCount.GetValue() << G4endl;
		G4cout << "DetectableEventCount =" << DetectableEventCount.GetValue() << G4endl;
		G4cout << "VetoPossibleEventCount =" << VetoPossibleEvtCount.GetValue() << G4endl;
		G4cout << "ROIEventCount =" << ROIEventCount.GetValue() << G4endl;
		G4cout << "ROIVetoEventCount =" << ROIVetoEventCount.GetValue() << G4endl;
		G4cout << "ROIVetoPossibleEventCount =" << ROIVetoPossibleEvtCount.GetValue() << G4endl;
		G4cout << G4endl;
		G4cout << G4endl;
		G4cout << "===============================================================" << G4endl;

		std::ofstream output;
		if (aRun->GetRunID() == 0) {
			output.open(txtname + ".txt", std::ios::ate);
			if (fPrimaryGenerator->GetSrcType() == "PENShell") {
				output
					<< "Source Distribution:\t" << "PENShell" << G4endl;
			}
			else if (fPrimaryGenerator->GetSrcType() == "Wire") {
				output
					<< "Source Distribution:\t" << "Wire" << G4endl
					<< "Wire Type:\t" << fDetCons->GetWireType() << G4endl
					<< "Reflector Type:\t" << fDetCons->GetReflectorType() << G4endl;
			}
			output
				//<< "Confine Info:\t" << fDetCons->GetConfine() << G4endl
				<< "Simulation Result:" << G4endl;
		}
		else
		{
			output.open(txtname + ".txt", std::ios::app);
		}
		G4double eff;
		if (BulkEventCount.GetValue() != 0) {
			eff = VetoEventCount.GetValue() / BulkEventCount.GetValue();
		}
		else
		{
			eff = 0;
		}

		G4double ROIeff;
		if (ROIEventCount.GetValue() != 0) {
			ROIeff = ROIVetoEventCount.GetValue() / ROIEventCount.GetValue();
		}
		else
		{
			ROIeff = 0;
		}

		if (aRun->GetRunID() == 0) {
			output
				<< std::setw(10) << std::left << "Run ID" << '\t'
				<< std::setw(40) << std::left << "Number of Event" << '\t'
				<< std::setw(40) << std::left << "Primary Particle" << '\t'
				<< std::setw(40) << std::left << "Primary Energy(MeV)" << '\t'
				<< std::setw(40) << std::left << "SiPMEvent" << '\t'
				<< std::setw(40) << std::left << "VetoEvent" << '\t'
				<< std::setw(40) << std::left << "VetoEvent_0" << '\t'
				<< std::setw(40) << std::left << "VetoEvent_1" << '\t'
				<< std::setw(40) << std::left << "VetoEvent_2" << '\t'
				<< std::setw(40) << std::left << "VetoEvent_3" << '\t'
				<< std::setw(40) << std::left << "VetoEvent_4" << '\t'
				<< std::setw(40) << std::left << "BulkEvent" << '\t'
				<< std::setw(40) << std::left << "DetectableEvent" << '\t'
				<< std::setw(40) << std::left << "VetoPossibleEvent" << '\t'
				<< std::setw(40) << std::left << "VetoEfficiency" << '\t'
				<< std::setw(40) << std::left << "ROIEvent" << '\t'
				<< std::setw(40) << std::left << "ROIVetoEvent" << '\t'
				<< std::setw(40) << std::left << "ROIVetoPossibleEvent" << '\t'
				<< std::setw(40) << std::left << "ROIVetoEfficiency" << '\t' << G4endl;
		}

		output
			<< std::setw(10) << std::left << aRun->GetRunID() << '\t'
			<< std::setw(40) << std::left << aRun->GetNumberOfEvent() << '\t'
			<< std::setw(40) << std::left << fPrimaryGenerator->GetPrimaryName() << '\t'
			<< std::setw(40) << std::left << std::setiosflags(std::ios::fixed) << std::setprecision(2) << fPrimaryGenerator->GetPrimaryE() << '\t'
			<< std::setw(40) << std::left << SiPMEventCount.GetValue() << '\t'
			<< std::setw(40) << std::left << VetoEventCount.GetValue() << '\t'
			<< std::setw(40) << std::left << VetoEventCount_0.GetValue() << '\t'
			<< std::setw(40) << std::left << VetoEventCount_1.GetValue() << '\t'
			<< std::setw(40) << std::left << VetoEventCount_2.GetValue() << '\t'
			<< std::setw(40) << std::left << VetoEventCount_3.GetValue() << '\t'
			<< std::setw(40) << std::left << VetoEventCount_4.GetValue() << '\t'
			<< std::setw(40) << std::left << BulkEventCount.GetValue() << '\t'
			<< std::setw(40) << std::left << DetectableEventCount.GetValue() << '\t'
			<< std::setw(40) << std::left << VetoPossibleEvtCount.GetValue() << '\t'
			<< std::setw(40) << std::left << eff << '\t'
			<< std::setw(40) << std::left << ROIEventCount.GetValue() << '\t'
			<< std::setw(40) << std::left << ROIVetoEventCount.GetValue() << '\t'
			<< std::setw(40) << std::left << ROIVetoPossibleEvtCount.GetValue() << '\t'
			<< std::setw(40) << std::left << ROIeff << G4endl;
		output.close();
		//std::DecimalFormat df1 = new DecimalFormat("0.0");
	}
}


