//
//  PhysicsEvent.cpp
//
//
//  Created by RENJIE WANG on 2/17/15.
//
//

#include "DMHiggsAnalysis/PhysicsEvent.h"

using namespace std;


//
PhysicsEvent_t getPhysicsEventFrom(DataEvtSummary_t &ev)
{
    PhysicsEvent_t phys;
//
    phys.run=ev.RunNumber;
    phys.event=ev.EventNumber;
    phys.lumi=ev.LumiBlock;
    phys.nvtx = ev.NPV;


    size_t npho(0);
    for(Int_t i=0; i<ev.nPhotons; i++) {
        LorentzVector P4(ev.photon_Px[i]/1000.,ev.photon_Py[i]/1000.,ev.photon_Pz[i]/1000.,ev.photon_E[i]/1000.);
        phys.photons.push_back( PhysicsObject_Photon(P4, 22) );
        npho++;
    }


    size_t nele(0);
    for(Int_t i=0; i<ev.nElectrons; i++) {
        LorentzVector P4(ev.electron_Px[i]/1000.,ev.electron_Py[i]/1000.,ev.electron_Pz[i]/1000.,ev.electron_E[i]/1000.);
        phys.electrons.push_back( PhysicsObject_Electron(P4, 11*ev.electron_charge[i]) );
        nele++;
    }


    size_t nmuon(0);
    for(Int_t i=0; i<ev.nMuons; i++) {
        LorentzVector P4(ev.muon_Px[i]/1000.,ev.muon_Py[i]/1000.,ev.muon_Pz[i]/1000.,ev.muon_E[i]/1000.);
        phys.muons.push_back( PhysicsObject_Muon(P4, 13*ev.muon_charge[i]) );
        nmuon++;
    }


    size_t njet(0);
    for(Int_t i=0; i<ev.nJets; i++) {
        LorentzVector P4(ev.jet_Px[i]/1000.,ev.jet_Py[i]/1000.,ev.jet_Pz[i]/1000.,ev.jet_E[i]/1000.);
        phys.jets.push_back( PhysicsObject_Jet(P4) );
        njet++;
    }


    // MET
    phys.met 	 = LorentzVector( ev.met*cos(ev.phi_met)/1000., ev.met*sin(ev.phi_met)/1000., 0, ev.met/1000. );


    return phys;
}


float getSFfrom2DHist(double xval, double yval, TH2F* h_)
{

    if(h_==NULL) {
        cout << "[getSFfrom2DHist]: empty hist! " << endl;
        return 1;
    }
    int xbins = h_->GetXaxis()->GetNbins();
    if(xval > h_->GetXaxis()->GetBinUpEdge(xbins)    ) xval = h_->GetXaxis()->GetBinUpEdge(xbins);
    if(xval < h_->GetXaxis()->GetBinLowEdge(1)       ) xval = h_->GetXaxis()->GetBinLowEdge(1);

    int ybins = h_->GetYaxis()->GetNbins();
    if(yval > h_->GetYaxis()->GetBinUpEdge(ybins)    ) yval = h_->GetYaxis()->GetBinUpEdge(ybins);
    if(yval < h_->GetYaxis()->GetBinLowEdge(1)       ) yval = h_->GetYaxis()->GetBinLowEdge(1);

    int binx = h_->GetXaxis()->FindBin(xval);
    int biny = h_->GetYaxis()->FindBin(yval);
    float sf_ = h_->GetBinContent(binx,biny);

    if(sf_==0.) return 1.;
    else return sf_;
}

float getSFfrom1DHist(double xval, TH1F* h_)
{

    if(h_==NULL) {
        cout << "[getSFfrom1DHist]: empty hist! " << endl;
        return 1;
    }
    int xbins = h_->GetXaxis()->GetNbins();
    if(xval > h_->GetXaxis()->GetBinUpEdge(xbins)    ) xval = h_->GetXaxis()->GetBinUpEdge(xbins);
    if(xval < h_->GetXaxis()->GetBinLowEdge(1)       ) xval = h_->GetXaxis()->GetBinLowEdge(1);


    int binx = h_->GetXaxis()->FindBin(xval);
    float sf_ = h_->GetBinContent(binx);

    if(sf_==0.) return 1.;
    else return sf_;
}




