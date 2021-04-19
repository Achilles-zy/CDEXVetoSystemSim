//CDEXEventAction.cc
#include <string.h>
#include "CDEXEventAction.hh"

#include "G4Event.hh"
#include "G4EventManager.hh"
#include "G4TrajectoryContainer.hh"
#include "G4Trajectory.hh"
#include "G4VVisManager.hh"
#include "G4ios.hh"
#include "g4root.hh"
#include "G4SystemOfUnits.hh"
#include "G4RunManager.hh"
//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

CDEXEventAction::CDEXEventAction(CDEXRunAction* runaction)
	: edepBulk(0.),
	SiPMPhotonCount{ 0 },
	SiPMSignalCount{ 0 },
	ContainerSiPMPhotonCount{ 0 },
	ContainerSiPMSignalCount{ 0 },
	Total(0),
	ID(0),
	EscapedPhotonCount(0),
	ifSiPM(false),
	ifBulk(false),
	ifDetectable(false),
	ifAccelerate(false),
	run(runaction)
	//ResultFile("Distribution_Results_NTuple.root","RECREATE"),
	//Distribution_Results("Distribution_Results","Distribution_Results")
{
	PhotonCut_0 = 1;
	PhotonCut_1 = 3;
	PhotonCut_2 = 5;
	PhotonCut_3 = 7;
	PhotonCut_4 = 10;
	SignalSiPMCount = 0;
	DepositeInfo;
	DepositeInfoInScintillator;
	ContainerSignalSiPMCount = 0;
	SignalSiPMCount_0 = 0;
	SignalSiPMCount_1 = 0;
	SignalSiPMCount_2 = 0;
	SignalSiPMCount_3 = 0;
	SignalSiPMCount_4 = 0;
	EnergyThreshold = 160 * eV;
	DepositeID = 0;
	RowNb = sizeof(SiPMPhotonCount) / sizeof(SiPMPhotonCount[0]);
	ColumnNb = sizeof(SiPMPhotonCount[0]) / sizeof(SiPMPhotonCount[0][0]);
	ContainerRowNb = sizeof(ContainerSiPMPhotonCount) / sizeof(ContainerSiPMPhotonCount[0]);
	ContainerColumnNb = sizeof(ContainerSiPMPhotonCount[0]) / sizeof(ContainerSiPMPhotonCount[0][0]);
}
//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

CDEXEventAction::~CDEXEventAction()
{
	//Distribution_Results.Write();
	//ResultFile.Close();
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

void CDEXEventAction::BeginOfEventAction(const G4Event* evt)
{
	edepBulk = 0;
	memset(SiPMPhotonCount, 0, sizeof(SiPMPhotonCount));
	memset(ContainerSiPMPhotonCount, 0, sizeof(ContainerSiPMPhotonCount));
	memset(SiPMSignalCount, 0, sizeof(SiPMPhotonCount));
	memset(ContainerSiPMSignalCount, 0, sizeof(ContainerSiPMPhotonCount));
	Total = 0;
	TotalSiPMPhotonCount = 0;
	SignalSiPMCount = 0;
	ContainerSignalSiPMCount = 0;
	SignalSiPMCount_0 = 0;
	SignalSiPMCount_1 = 0;
	SignalSiPMCount_2 = 0;
	SignalSiPMCount_3 = 0;
	SignalSiPMCount_4 = 0;
	MinSignalSiPMCount = 1;
	EscapedPhotonCount = 0;
	ifSiPM = false;
	ifBulk = false;
	ifDetectable = false;
	ifAccelerate = false;

	DepositeInfo.clear();
	DepositeInfoInScintillator.clear();

	EdepInfo.clear();
	TempPosList.clear();
	EdepInfoInScintillator.clear();
	TempPosListInScintillator.clear();
	//G4cout << evt->GetEventID() << G4endl;
	// G4cout<<ID<<G4endl;
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

void CDEXEventAction::EndOfEventAction(const G4Event* evt)
{
	auto analysisManager = G4AnalysisManager::Instance();
	analysisManager->FillH1(0, edepBulk);
	analysisManager->FillH1(1, TotalSiPMPhotonCount);
	analysisManager->FillH1(2, TotalSiPMPhotonCount);
	analysisManager->FillNtupleDColumn(2, 0, edepBulk);
	analysisManager->AddNtupleRow(2);
	//analysisManager->FillNtupleDColumn(2, 1, edepBulk);

	//G4int rownb = sizeof(SiPMPhotonCount) / sizeof(SiPMPhotonCount[0]);
	//G4int columnnb = sizeof(SiPMPhotonCount[0]) / sizeof(SiPMPhotonCount[0][0]);

	for (int i = 0; i < RowNb; i++)
	{
		for (int j = 0; j < ColumnNb; j++)
		{
			//G4cout << SiPMSignalCount[i][j] << " ";
			G4int NtupleColumnID = i * ColumnNb + j;
			//analysisManager->FillNtupleIColumn(1, NtupleColumnID, TotalSiPMPhotonCount);
			if (SiPMPhotonCount[i][j] >= PhotonCut_0) {
				SignalSiPMCount_0++;
			}
			if (SiPMPhotonCount[i][j] >= PhotonCut_1) {
				SignalSiPMCount_1++;
			}
			if (SiPMPhotonCount[i][j] >= PhotonCut_2) {
				SignalSiPMCount_2++;
			}
			if (SiPMPhotonCount[i][j] >= PhotonCut_3) {
				SignalSiPMCount_3++;
			}
			if (SiPMPhotonCount[i][j] >= PhotonCut_4) {
				SignalSiPMCount_4++;
			}
			if (SiPMSignalCount[i][j] > 0) {
				SignalSiPMCount++;
			}
		}
		//G4cout<<G4endl;
	}
	//G4cout << "Matrix End" << G4endl;


	for (int i = 0; i < ContainerRowNb; i++)
	{
		for (int j = 0; j < ContainerColumnNb; j++)
		{
			//G4cout << SiPMPhotonCount[i][j] << " ";
			G4int NtupleColumnID = i * ContainerColumnNb + j;
			//analysisManager->FillNtupleIColumn(1, NtupleColumnID, TotalSiPMPhotonCount);
			if (ContainerSiPMSignalCount[i][j] > 0) {
				ContainerSignalSiPMCount++;
			}
		}
		//G4cout<<G4endl;
	}



	analysisManager->FillNtupleIColumn(1, 0, TotalSiPMPhotonCount);
	analysisManager->AddNtupleRow(1);
	//G4cout << SiPMPhotonCount << G4endl;

	G4int TotalSignalSiPMCount = ContainerSignalSiPMCount + SignalSiPMCount;
	if (edepBulk > EnergyThreshold && TotalSignalSiPMCount >= MinSignalSiPMCount) {
		//G4cout << "SignalSiPMCount0Ture" << G4endl;
		run->CountVetoEvent();
	}

	if (edepBulk > EnergyThreshold && SignalSiPMCount_0 >= MinSignalSiPMCount) {
		//G4cout << "SignalSiPMCount0Ture" << G4endl;
		run->CountVetoEvent_0();
	}
	if (edepBulk > EnergyThreshold && SignalSiPMCount_1 >= MinSignalSiPMCount) {
		//G4cout << "SignalSiPMCount1Ture" << G4endl;
		run->CountVetoEvent_1();
	}
	if (edepBulk > EnergyThreshold && SignalSiPMCount_2 >= MinSignalSiPMCount) {
		//G4cout << "SignalSiPMCount2Ture" << G4endl;
		run->CountVetoEvent_2();
	}
	if (edepBulk > EnergyThreshold && SignalSiPMCount_3 >= MinSignalSiPMCount) {
		//G4cout << "SignalSiPMCount3Ture" << G4endl;
		run->CountVetoEvent_3();
	}

	if (edepBulk > EnergyThreshold && SignalSiPMCount_4 >= MinSignalSiPMCount) {
		//G4cout << "SignalSiPMCount3Ture" << G4endl;
		run->CountVetoEvent_4();
	}

	if (TotalSignalSiPMCount > 0) {
	//for (int i = 0; i < RowNb; i++)
	//{
	//	for (int j = 0; j < ColumnNb; j++)
	//	{
	//		G4cout << SiPMSignalCount[i][j] << " ";
	//	}
	//	G4cout << G4endl;
	//}
		//G4cout << "TotalSignalSiPMCount = " << TotalSignalSiPMCount << G4endl;
		//G4cout << "Matrix End" << G4endl;
		run->CountSiPMEvent();
	}

	if (edepBulk > EnergyThreshold) {
		run->CountBulkEvent();
		analysisManager->FillH1(3, TotalSiPMPhotonCount);
	}

	if (ifDetectable == true) {
		run->CountDetectableEvent();
	}

	if (ifDetectable == true && edepBulk > EnergyThreshold) {
		run->CountVetoPossibleEvent();
	}

	if (edepBulk > 2000 * keV && edepBulk < 2100 * keV) {
		run->CountROIEvent();

		if (TotalSignalSiPMCount >= MinSignalSiPMCount) {
			run->CountROIVetoEvent();
		}

		if (ifDetectable == true) {
			run->CountROIVetoPossibleEvent();
		}
	}

	G4int evtID = evt->GetEventID();

	if (evtID % 5000 == 0) {
		G4cout << evtID << G4endl;
	}

	if (EdepInfo.empty() == false) {
		for (G4int i = 0; i < EdepInfo.size(); i++) {
			analysisManager->FillNtupleIColumn(3, 0, EdepInfo[i][0]);
			analysisManager->FillNtupleIColumn(3, 1, EdepInfo[i][1]);
			analysisManager->FillNtupleDColumn(3, 2, EdepInfo[i][2]);
			analysisManager->FillNtupleDColumn(3, 3, EdepInfo[i][3]);
			analysisManager->FillNtupleDColumn(3, 4, EdepInfo[i][4]);
			analysisManager->FillNtupleDColumn(3, 5, EdepInfo[i][5]);
			analysisManager->AddNtupleRow(3);

			if (edepBulk > 160 * eV) {
				analysisManager->FillNtupleIColumn(4, 0, EdepInfo[i][0]);
				analysisManager->FillNtupleIColumn(4, 1, EdepInfo[i][1]);
				analysisManager->FillNtupleDColumn(4, 2, EdepInfo[i][2]);
				analysisManager->FillNtupleDColumn(4, 3, EdepInfo[i][3]);
				analysisManager->FillNtupleDColumn(4, 4, EdepInfo[i][4]);
				analysisManager->FillNtupleDColumn(4, 5, EdepInfo[i][5]);
				analysisManager->AddNtupleRow(4);
			}
		}
	}

	if (EdepInfoInScintillator.empty() == false) {
		for (G4int i = 0; i < EdepInfoInScintillator.size(); i++) {
			analysisManager->FillNtupleIColumn(5, 0, EdepInfoInScintillator[i][0]);
			analysisManager->FillNtupleIColumn(5, 1, EdepInfoInScintillator[i][1]);
			analysisManager->FillNtupleDColumn(5, 2, EdepInfoInScintillator[i][2]);
			analysisManager->FillNtupleDColumn(5, 3, EdepInfoInScintillator[i][3]);
			analysisManager->FillNtupleDColumn(5, 4, EdepInfoInScintillator[i][4]);
			analysisManager->FillNtupleDColumn(5, 5, EdepInfoInScintillator[i][5]);
			analysisManager->AddNtupleRow(5);

			if (edepBulk > 160 * eV) {
				analysisManager->FillNtupleIColumn(6, 0, EdepInfoInScintillator[i][0]);
				analysisManager->FillNtupleIColumn(6, 1, EdepInfoInScintillator[i][1]);
				analysisManager->FillNtupleDColumn(6, 2, EdepInfoInScintillator[i][2]);
				analysisManager->FillNtupleDColumn(6, 3, EdepInfoInScintillator[i][3]);
				analysisManager->FillNtupleDColumn(6, 4, EdepInfoInScintillator[i][4]);
				analysisManager->FillNtupleDColumn(6, 5, EdepInfoInScintillator[i][5]);
				analysisManager->AddNtupleRow(6);
			}
		}
	}

	//if (DepositeInfoInScintillator.empty() == false) {
	//	for (G4int i = 0; i < DepositeInfoInScintillator.size(); i++) {
	//		analysisManager->FillNtupleIColumn(5, 0, DepositeInfoInScintillator[i][0]);
	//		analysisManager->FillNtupleIColumn(5, 1, DepositeInfoInScintillator[i][1]);
	//		if (edepBulk > 160 * eV) {
	//			analysisManager->FillNtupleIColumn(6, 0, DepositeInfoInScintillator[i][0]);
	//			analysisManager->FillNtupleIColumn(6, 1, DepositeInfoInScintillator[i][1]);
	//			analysisManager->AddNtupleRow(6);
	//		}
	//		//G4cout << DepositeInfo[i][0] << " " ;
	//		for (G4int j = 2; j < DepositeInfoInScintillator[0].size(); j++) {
	//			analysisManager->FillNtupleDColumn(5, j, DepositeInfoInScintillator[i][j]);
	//			if (edepBulk > 160 * eV) {
	//				analysisManager->FillNtupleDColumn(6, j, DepositeInfoInScintillator[i][j]);
	//				analysisManager->AddNtupleRow(6);
	//			}
	//			//G4cout << DepositeInfo[i][j] << " " ;
	//		}
	//		analysisManager->AddNtupleRow(5);

	//		//G4cout << G4endl;
	//	}

	//}
	//G4cout << "Size=" << DepositeInfo.size() << G4endl;
	//G4cout << "Size2=" << DepositeInfo[0].size() << G4endl;
	//ID++;

}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

void CDEXEventAction::AddBulkEnergy(G4double de)
{
	edepBulk += de;
}

void CDEXEventAction::AddToSiPM(G4int i, G4int j) {
	Total++;
	SiPMPhotonCount[i][j]++;
}

void CDEXEventAction::AddToSiPMSignal(G4int i, G4int j) {
	SiPMSignalCount[i][j]++;
}

void CDEXEventAction::AddToContainerSiPM(G4int i, G4int j) {
	Total++;
	ContainerSiPMPhotonCount[i][j]++;
}

void CDEXEventAction::AddToContainerSiPMSignal(G4int i, G4int j) {
	ContainerSiPMSignalCount[i][j]++;
}

void CDEXEventAction::RecordStepInfo(G4int particletype, G4int creatorprocess, G4double posx, G4double posy, G4double posz, G4double edep) {
	std::vector<G4double> StepInfo;
	StepInfo.push_back(particletype);
	StepInfo.push_back(creatorprocess);
	StepInfo.push_back(posx);
	StepInfo.push_back(posy);
	StepInfo.push_back(posz);
	StepInfo.push_back(edep);

	DepositeInfo.push_back(StepInfo);
	StepInfo.clear();
}

void CDEXEventAction::RecordStepInfoInScintillator(G4int particletype, G4int creatorprocess, G4double posx, G4double posy, G4double posz, G4double edep) {
	std::vector<G4double> StepInfo;
	StepInfo.push_back(particletype);
	StepInfo.push_back(creatorprocess);
	StepInfo.push_back(posx);
	StepInfo.push_back(posy);
	StepInfo.push_back(posz);
	StepInfo.push_back(edep);

	DepositeInfoInScintillator.push_back(StepInfo);
	StepInfo.clear();
}

void CDEXEventAction::RecordEdepInfo(G4int particletype, G4int creatorprocess, G4double posx, G4double posy, G4double posz, G4double edep) {
	std::vector<G4double> StepInfo;

	StepInfo.push_back(particletype);
	StepInfo.push_back(creatorprocess);
	StepInfo.push_back(posx);
	StepInfo.push_back(posy);
	StepInfo.push_back(posz);
	StepInfo.push_back(edep);

	std::vector<G4double> Pos;
	Pos.push_back(posx);
	Pos.push_back(posy);
	Pos.push_back(posz);

	if (TempPosList.empty()) {
		TempPosList.push_back(Pos);
		EdepInfo.push_back(StepInfo);
	}
	else
	{
		G4bool ifNewPoint = true;
		G4int PointID;
		for (G4int i = 0; i < TempPosList.size(); i++) {
			G4double dist = GetDistance(TempPosList[i][0], TempPosList[i][1], TempPosList[i][2], Pos[0], Pos[1], Pos[2]);
			if (dist < 5 * mm) {
				ifNewPoint = false;
				PointID = i;
			}
		}
		if (ifNewPoint == true) {
			TempPosList.push_back(Pos);
			EdepInfo.push_back(StepInfo);
		}
		else
		{
			EdepInfo[PointID][5] += edep;
			if (particletype == 1) {
				EdepInfo[PointID][0] == 1;
			}
		}
	}

	StepInfo.clear();
	Pos.clear();
}

void CDEXEventAction::RecordEdepInfoInScintillator(G4int particletype, G4int creatorprocess, G4double posx, G4double posy, G4double posz, G4double edep) {
	
	std::vector<G4double> StepInfo;
	StepInfo.push_back(particletype);
	StepInfo.push_back(creatorprocess);
	StepInfo.push_back(posx);
	StepInfo.push_back(posy);
	StepInfo.push_back(posz);
	StepInfo.push_back(edep);

	std::vector<G4double> Pos;
	Pos.push_back(posx);
	Pos.push_back(posy);
	Pos.push_back(posz);

	if (TempPosListInScintillator.empty()) {
		TempPosListInScintillator.push_back(Pos);
		EdepInfoInScintillator.push_back(StepInfo);
	}
	else
	{
		G4bool ifNewPoint = true;
		G4int PointID;
		for (G4int i = 0; i < TempPosListInScintillator.size(); i++) {
			G4double dist = GetDistance(TempPosListInScintillator[i][0], TempPosListInScintillator[i][1], TempPosListInScintillator[i][2], Pos[0], Pos[1], Pos[2]);
			if (dist < 5 * mm) {
				ifNewPoint = false;
				PointID = i;
			}
		}
		if (ifNewPoint == true) {
			TempPosListInScintillator.push_back(Pos);
			EdepInfoInScintillator.push_back(StepInfo);
		}
		else
		{
			EdepInfoInScintillator[PointID][5] += edep;
			if (particletype == 1) {
				EdepInfoInScintillator[PointID][0] == 1;
			}
		}
	}

	StepInfo.clear();
	Pos.clear();
}

G4double CDEXEventAction::GetDistance(G4double x0, G4double y0, G4double z0, G4double x1, G4double y1, G4double z1) {
	G4double res;
	res = std::sqrt((x0 - x1) * (x0 - x1) + (y0 - y1) * (y0 - y1) + (z0 - z1) * (z0 - z1));
	return res;
}
