extern "C" {
	extern struct{
		int pdg, evt, bch, evtcas;
		double x1, y1, z1, px, py, pz, ee, wi;
	} root_;
}

#include <iostream>
#include <TGraph.h>
#include <TCanvas.h>
#include <TFile.h>
#include <TTree.h>

using namespace std;

TFile* file = NULL;
TTree* tree = NULL;

int EvtPdgcode;
int EvtEvt;
int EvtBch;
int EvtCas;
double EvtX;
double EvtY;
double EvtZ;
double EvtPX;
double EvtPY;
double EvtPZ;
double EvtE;
double EvtWght;

bool set_up = false;

extern "C"{
	void rootlink_(){
		if (set_up == false){
			file = new TFile("output.root", "recreate");
			tree = new TTree("gROOTracker", "PHITS event tree");

			tree->Branch("pdgcode", &EvtPdgcode, "pdgcode/I");
			tree->Branch("event", &EvtEvt, "event/I");
			tree->Branch("bunch", &EvtBch, "bunch/I");
			tree->Branch("cascade", &EvtCas, "cascade/I");
			tree->Branch("x", &EvtX, "x/D");
			tree->Branch("y", &EvtY, "y/D");
			tree->Branch("z", &EvtZ, "z/D");
			tree->Branch("px", &EvtPX, "px/D");
			tree->Branch("py", &EvtPY, "py/D");
			tree->Branch("pz", &EvtPZ, "pz/D");
			tree->Branch("energy", &EvtE, "energy/D");
			tree->Branch("weight", &EvtWght, "weight/D");
			
			set_up = true;
		}

		EvtPdgcode = root_.pdg;
		EvtEvt = root_.evt;
		EvtBch = root_.bch;
		EvtCas = root_.evtcas;
		EvtX = root_.x1;
		EvtY = root_.y1;
		EvtZ = root_.z1;
		EvtPX = root_.px;
		EvtPY = root_.py;
		EvtPZ = root_.pz;
		EvtE = root_.ee;
		EvtWght = root_.wi;
	
		tree->Fill();
		tree->AutoSave();
	}
}
