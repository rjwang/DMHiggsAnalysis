//
//  PhysicsEvent.h
//
//
//  Created by RENJIE WANG on 2/17/15.
//
//

#ifndef ____PhysicsEvent__
#define ____PhysicsEvent__

#include <stdio.h>
#include <vector>
#include "TH2F.h"

#include "DMHiggsAnalysis/DataEvtSummaryHandler.h"

enum PhysicsObjects   { MET=0, JET=1, TOP=6, ELECTRON=11, MUON=13, TAU=15, GLUON=21, PHOTON=22, Z=23, W=24};



class PhysicsObject : public LorentzVector {
public :
    PhysicsObject(LorentzVector vec, int id_):
        LorentzVector(vec), id(id_) { }
    int id;
};



//
class PhysicsObject_Photon : public LorentzVector {
public :
    PhysicsObject_Photon(LorentzVector vec, int id_):
        LorentzVector(vec), id(id_) { }
    int id;
};



//
class PhysicsObject_Electron : public LorentzVector {
public :
    PhysicsObject_Electron(LorentzVector vec, int id_):
        LorentzVector(vec), id(id_) { }
    int id;
};


//
class PhysicsObject_Muon : public LorentzVector {
public :
    PhysicsObject_Muon(LorentzVector vec, int id_):
        LorentzVector(vec), id(id_) { }
    int id;
};



//
class PhysicsObject_Jet : public LorentzVector {
public :
    PhysicsObject_Jet(LorentzVector vec):
        LorentzVector(vec) { }
};


typedef std::vector<PhysicsObject>        	PhysicsObjectCollection;
typedef std::vector<PhysicsObject_Electron> 	PhysicsObjectElectronCollection;
typedef std::vector<PhysicsObject_Muon>   	PhysicsObjectMuonCollection;
typedef std::vector<PhysicsObject_Photon> 	PhysicsObjectPhotonCollection;
typedef std::vector<PhysicsObject_Jet>    	PhysicsObjectJetCollection;

//
struct PhysicsEvent_t {
    int run,event,lumi;
    int nvtx;

    PhysicsObjectElectronCollection electrons;
    PhysicsObjectMuonCollection muons;
    PhysicsObjectPhotonCollection photons;
    PhysicsObjectJetCollection jets;
    LorentzVector met;
    LorentzVector met_hv;

    float pv_hard_z;
    float pv_diphot_z;

    PhysicsObjectCollection genneutrinos,genleptons,genWIMPs,genGravitons,genPhotons;
    PhysicsObjectCollection genjets;
};



//
PhysicsEvent_t getPhysicsEventFrom(DataEvtSummary_t &ev);

float getSFfrom2DHist(double xval, double yval, TH2F* h_);
float getSFfrom1DHist(double xval, TH1F* h_);




#endif /* defined(____PhysicsEvent__) */
