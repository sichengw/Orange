// C++
#include <iostream>
#include <stdio.h>  
#include <vector>
#include <math.h>
#include <iomanip>

// ROOT
#include "TBenchmark.h"
#include "TCanvas.h"
#include "TChain.h"
#include "TDirectory.h"
#include "TFile.h"
#include "TROOT.h"
#include "TTreeCache.h"
#include "TLegend.h"
#include "TRandom1.h"
#include "TPad.h"
#include "TH2F.h"
#include "TH1.h"
#include "THStack.h"
#include "TMath.h"
#include "Math/VectorUtil.h"
#include "Math/Vector4D.h"
#include "TLorentzVector.h"

// CMS2
#include "CMS2.h"
#include "muonSelections.h"
#include "susySelections.h"
#include "ssSelections.h"
#include "eventSelections.h"
#include "MT2/MT2.h"
#include "MT2/MT2Utility.h"


// Good run list

#include "/home/users/jgran/CMSSW_5_3_2_patch4_V05-03-23/src/CMS2/NtupleMacros/Tools/goodrun.cc"

// My includes
#include "myheader.h"
#include "NSBABY.h"
typedef ROOT::Math::LorentzVector<ROOT::Math::PxPyPzE4D<float> > LorentzVec;

using namespace std;
//using namespace tas;
using namespace ROOT::Math;
using namespace nsbaby;


int ScanChain( TChain* chain, char* sampleType = "", bool ismc = true, bool fast = true, int nEvents = -1, string skimFilePrefix = "test") {

  // Benchmark
  TBenchmark *bmark = new TBenchmark();
  bmark->Start("benchmark");

  const int FIRTBIN = 0;
  const int LASTBIN = 200;
  const int BINNUM = 50;
  const int INVMLASTBIN = 360;
  const int METLASTBIN = 300;
  const int MTLASTBIN = 200;

  // TFile* f1 = new TFile("./hists/hists_mcW_sl_1.root");
  TFile* f1 = new TFile("./hists/hists_odata_40cmetg1j0b30lpt90HTTTI_recW_12.root");
  TH1F* RWPT = (TH1F*)f1->Get("wpt"); 
  // f1->Close();

////////////////////////////////////////////////////////////////////////////////////////////
  TH1F* MT2_hist  = new TH1F("MT2","MT2 distribution", BINNUM, FIRTBIN, LASTBIN);
  TH1F* MT2_el    = new TH1F("MT2e","MT2 distribution for e  events", BINNUM, FIRTBIN, LASTBIN);
  TH1F* MT2_mu    = new TH1F("MT2m","MT2 distribution for mu events", BINNUM, FIRTBIN, LASTBIN);
  TH1F* h2_nozero = new TH1F("hnz","MT2 without 0 bin", BINNUM, 1, LASTBIN);
  TH1F* MT2_Cut50 = new TH1F("h50","MT2 without 0 bin with met cut 50", BINNUM, 1, LASTBIN);
        
  TH1F* MT_hist  = new TH1F("MT",  "MT distribution", BINNUM, FIRTBIN, MTLASTBIN+1);
  TH1F* MTc_hist = new TH1F("MTc", "MT distribution", BINNUM, FIRTBIN, MTLASTBIN+1);
  TH1F* MT_el    = new TH1F("MTe", "MT distribution for e  events", BINNUM, FIRTBIN, MTLASTBIN+1);
  TH1F* MT_mu    = new TH1F("MTm", "MT distribution for mu events", BINNUM, FIRTBIN, MTLASTBIN+1);
        
  TH1F* JetMult_a = new TH1F("jeta", "Jet mutiplicity all", 13, 0, 13);
  TH1F* Jeta_el   = new TH1F("jae",  "Jet mutiplicity for e  events", 13, 0, 13);
  TH1F* Jeta_mu   = new TH1F("jam",  "Jet mutiplicity for mu events", 13, 0, 13);
  TH1F* JetMult_b = new TH1F("jetb", "b-Jet mutiplicity", 13, 0, 13);
  TH1F* Jetb_el   = new TH1F("jbe",  "b-Jet mutiplicity for e  events", 13, 0, 13);
  TH1F* Jetb_mu   = new TH1F("jbm",  "b-Jet mutiplicity for mu events", 13, 0, 13);
        
  TH1F* Met_all = new TH1F("met",  "MET for all events",  60, 0, METLASTBIN+1);
  TH1F* Met_el  = new TH1F("mete",  "MET for e  events",  60, 0, METLASTBIN+1);
  TH1F* Met_mu  = new TH1F("metm",  "MET for mu events",  60, 0, METLASTBIN+1);
  TH1F* FMet    = new TH1F("fmet",  "Real Met + Fake #nu",60, 0, METLASTBIN+1);
        
  TH1F* MetPhi    = new TH1F("metphi",  "MetPhi ",		     90, -3.3, 3.3);
  TH1F* MetPhi_el = new TH1F("metphie",  "Corrected MetPhi of elec", 90, -3.3, 3.3);
  TH1F* MetPhi_mu = new TH1F("metphim",  "Corrected MetPhi of muon", 90, -3.3, 3.3);
  TH1F* MetPhiCor = new TH1F("cmetphi", "Corrected MetPhi",       90, -3.3, 3.3);
  TH1F* CMetPhi_e = new TH1F("cmpe",  "Corrected MetPhi of elec", 90, -3.3, 3.3);
  TH1F* CMetPhi_m = new TH1F("cmpm",  "Corrected MetPhi of muon", 90, -3.3, 3.3);
  TH1F* FMetPhi  = new TH1F("fmetphi", "Metphi + fake nu",         90, -3.3, 3.3);
  TH1F* Metx    = new TH1F("metx", "Met x component",           120, -240, 240);
  TH1F* Mety    = new TH1F("mety", "Met y component",           120, -240, 240);
  TH1F* MetCorx = new TH1F("cmetx","Corrected Met x component", 120, -240, 240);
  TH1F* MetCory = new TH1F("cmety","Corrected Met y component", 120, -240, 240);
  TH1F* DiffPhi  = new TH1F("dphi",  "Difference in phi for lep and met",  60, 0, 3.3);
  TH1F* DiffPhic = new TH1F("dphic", "Difference in cphi for lep and met", 60, 0, 3.3);
  TH1F* DiffPhi_e  = new TH1F("dpe",  "Difference in phi for ele and met", 60, 0, 3.3);
  TH1F* DiffPhi_m  = new TH1F("dpm",  "Difference in phi for mu and met",  60, 0, 3.3);
        
  TH1F* InvM_lep = new TH1F("invm","Dileptons InvMass all",   70, 0, INVMLASTBIN+1);
  TH1F* FW_P     = new TH1F("fwp",  "Momentum of the fake W ",  80, 0, 1200);
  TH1F* W_P      = new TH1F("rwp",  "W_P value by (lpt + met)/cos#theta ",  80, 0, 1601);
  TH1F* W_P_el   = new TH1F("rwpe",  "W_P value by (lpt + met)/cos#theta elec",   80, 0, 1601);
  TH1F* W_P_mu   = new TH1F("rwpm",  "W_P value by (lpt + met)/cos#theta muon",  80, 0, 1601);
  TH1F* W_PT     = new TH1F("rwpt","Pt value of lpt + met ",  70, 0, 350);
  TH1F* W_PT_el  = new TH1F("wpte","Pt value of lpt + met for ele",   70, 0, 350);
  TH1F* W_PT_mu  = new TH1F("wptm","Pt value of lpt + met for muon",  70, 0, 350);
  TH1F* Lep_PT = new TH1F("lpt","Pt value of leptons ",    70, 0, 350);
  TH1F* Mu_PT  = new TH1F("mpt","Pt value of muons ",      70, 0, 350);
  TH1F* El_PT  = new TH1F("ept","Pt value of electron ",   70, 0, 350);
  TH1F* FLep_PT = new TH1F("flpt","Pt value of the fake leptons ",    70, 0, 350);
  TH1F* TEST  = new TH1F("test","Testing histo for an interesting variable",    70, 0, 350);
  TH1F* TEST3 = new TH1F("test3","Testing histo 3 for an interesting variable", 70, 0, 350);
  TH1F* TEST2 = new TH1F("test2","Testing histo 2 for an interesting variable", 70, 0, 3.2);
  TH1F* FLepEta = new TH1F("flEta","Eta Distribution of fake leps",  90, -2.7, 2.7);
  TH1F* LepEta = new TH1F("lepEta","Eta Distribution of real leps",  90, -2.7, 2.7);
  TH1F* ElEta  = new TH1F("elEta" ,"Eta Distribution of electrons",  90, -2.7, 2.7);
  TH1F* MuEta  = new TH1F("muEta" ,"Eta Distribution of muons",      90, -2.7, 2.7);
  TH1F* LepPhi = new TH1F("lepPhi","Phi Distribution of real leps",  90, -3.4, 3.4);
  TH1F* ElPhi  = new TH1F("elPhi" ,"Phi Distribution of electrons",  90, -3.4, 3.4);
  TH1F* MuPhi  = new TH1F("muPhi" ,"Phi Distribution of muons",      90, -3.4, 3.4);

  TH1F* JPt1     = new TH1F("jpt1"  ,"Highest Jet Pt all",	       70, 0, 350);
  TH1F* JPt1_el  = new TH1F("jpt1e" ,"Highest Jet Pt for elec",        70, 0, 350);
  TH1F* JPt1_mu  = new TH1F("jpt1m" ,"Highest Jet Pt for muon",        70, 0, 350);
  TH1F* JPt2     = new TH1F("jpt2"  ,"2nd Highest Jet Pt all",	       70, 0, 350);
  TH1F* JPt2_el  = new TH1F("jpt2e" ,"2nd Highest Jet Pt for elec",    70, 0, 350);
  TH1F* JPt2_mu  = new TH1F("jpt2m" ,"2nd Highest Jet Pt for muon",    70, 0, 350);
  TH1F* Jdphi    = new TH1F("jdp"   ,"delta Phi of the highest pt 2 jets all",   60, 0, 3.3);
  TH1F* Jdphi_el = new TH1F("jdpe"  ,"delta Phi of the highest pt 2 jets elec",  60, 0, 3.3);
  TH1F* Jdphi_mu = new TH1F("jdpm"  ,"delta Phi of the highest pt 2 jets muon",  60, 0, 3.3);
  TH1F* SJPt     = new TH1F("sjpt"  ,"Pt of vector sum of highest 2 jets all",    70, 0, 350);
  TH1F* SJPt_el  = new TH1F("sjpte" ,"Pt of vector sum of highest 2 jets elec",   70, 0, 350);
  TH1F* SJPt_mu  = new TH1F("sjptm" ,"Pt of vector sum of highest 2 jets muon",   70, 0, 350);
  TH1F* JetHT    = new TH1F("jht"  ,"Sum of Pt of all jets",   60, 0, 600);
  TH1F* JetHT_el = new TH1F("jhte" ,"HT for e  events",   60, 0, 600);
  TH1F* JetHT_mu = new TH1F("jhtm" ,"HT for mu events",   60, 0, 600);

  TH2F* EtaPhi    = new TH2F("etaphi" , "Phi vs Eta Distribution of real leps",  90, -2.7, 2.7, 90, -3.4, 3.4);
  TH2F* EtaPhi_el = new TH2F("etaphie", "Phi vs Eta Distribution of electrons",  90, -2.7, 2.7, 90, -3.4, 3.4);
  TH2F* EtaPhi_mu = new TH2F("etaphim", "Phi vs Eta Distribution of muons",      90, -2.7, 2.7, 90, -3.4, 3.4);

  MT_hist   ->Sumw2(); 
  MT_el     ->Sumw2();
  MT_mu     ->Sumw2();
  MTc_hist  ->Sumw2(); 
  MT2_hist  ->Sumw2(); 
  MT2_el    ->Sumw2();
  MT2_mu    ->Sumw2();
  h2_nozero ->Sumw2();

  Met_all   ->Sumw2();
  Met_el    ->Sumw2();
  Met_mu    ->Sumw2();

  W_PT      ->Sumw2();
  W_PT_el   ->Sumw2();
  W_PT_mu   ->Sumw2();

  InvM_lep  ->Sumw2();
  JetMult_a ->Sumw2();
  JetMult_b ->Sumw2();
  Jeta_el   ->Sumw2();
  Jeta_mu   ->Sumw2();
  FMet      ->Sumw2();


///////////////////////////////////////////////////////////////////////////////////////////

  int file_count = 0;

  // int less_jets = 0;
  // int nGoodEvents = 0;
  // float bTagDiscriminant = 0.244; 

  // Loop over events to Analyze
  unsigned int nEventsTotal = 0;
  unsigned int nEventsChain = chain->GetEntries();
  if( nEvents >= 0 ) nEventsChain = nEvents;
  TObjArray *listOfFiles = chain->GetListOfFiles();
  TIter fileIter(listOfFiles);
  TFile *currentFile = 0;

  unsigned int nDuplicates = 0;
  // Set Good Run List
  // set_goodrun_file("/home/users/jgran/analysis/sswh/fakes/json/final_19p49fb_cms2.txt");

  // File Loop
  while ( (currentFile = (TFile*)fileIter.Next()) ) {
    
    // Get File Content
    TFile *file = new TFile( currentFile->GetTitle() );
    TTree *tree = (TTree*)file->Get("tree");
    if(fast) TTreeCache::SetLearnEntries(10);
    if(fast) tree->SetCacheSize(128*1024*1024);
    // cms2.Init(tree);
    baby.Init(tree);  
    // Loop over Events in current file
    if( nEventsTotal >= nEventsChain ) continue;
    unsigned int nEventsTree = tree->GetEntriesFast();
    for( unsigned int event = 0; event < nEventsTree; ++event) {
    
      // Get Event Content
      if( nEventsTotal >= nEventsChain ) continue;
      if(fast) tree->LoadTree(event);
      // cms2.GetEntry(event);
      baby.GetEntry(event);  
      ++nEventsTotal;

      if(isRealData()){
        DorkyEventIdentifier id = {runNumber(), eventNumber(), lumiBlock()};
        if (is_duplicate(id)){
	nDuplicates++;
          continue;
        }
      }
    
      // Progress	
      CMS2::progress( nEventsTotal, nEventsChain );

      // Select Good Runs

      Int_t n_jets = 0;
      Int_t n_bTag = 0;
 	
      // Applying cuts
      if( lr_p4().pt() < 30  )			 continue;
      if( fabs(lr_p4().eta()) > 2.1 )		 continue;
      if( abs(lr_id()) == 11 && fabs(lr_p4().eta()) > 1.1 )		 continue;

      if(isRealData()){
      	if( abs(lr_id()) == 11 && !ele27wp80() ) continue;
      	if( abs(lr_id()) == 13 && !isomu24()   ) continue;
      }

      if( njets() < 2  )			 continue;
      if( nbTag() != 0 )			 continue;

      if( trackingProblemVeto() )                continue;    
      if( tauVeto() )                            continue;
      if( isoTrackVeto() )                       continue;


      //// Do the MetPhi correction
      float metx = met() * cos( metPhi() );
      float mety = met() * sin( metPhi() );

      float shiftx = 0.;
      float shifty = 0.;

      shiftx = (! isRealData()) ?  (0.1166 + 0.0200*nvtxs()) : (0.2661 + 0.3217*nvtxs());
      shifty = (! isRealData()) ?  (0.2764 - 0.1280*nvtxs()) : (-0.2251 - 0.1747*nvtxs());
      
      metx -= shiftx;
      mety -= shifty;

      float cmet = sqrt( metx*metx + mety*mety ); // cmet = corrected met
      float cmetphi = atan2( mety , metx );

      float metu = cmet;			// metu = met to be used in the following calculations
      float metphiu = cmetphi;

      if(metu < 40 )    continue;

      /// HT cut 
      float sum_jetPt = 0;
      for(unsigned int c = 0; c < jets_p4().size(); c++) {
      	float _jetPt = jets_p4().at(c).pt() * jets_p4Correction().at(c);
      	// jet pt times correction to jets
      	if(_jetPt < 30) continue;       // disgard those with small pt

      	if(fabs(jets_p4().at(c).eta()) > 2.5) continue;   // Jet eta check

      	float dr_lr = VectorUtil::DeltaR(jets_p4().at(c), lr_p4() );
      	// delta R is the distance in the eta-phi plane
      	if(dr_lr < 0.4) continue;      // disgard small delta R 
	
	sum_jetPt += _jetPt;

      }//pfjets_p4().size()

      if(sum_jetPt < 90) continue;
      
      /// all selection 

      float scale1fb = ( isRealData() )? 1		// : scale_1fb();  
	: (abs(lr_id()) == 13) ? getTriggerEfficiency_HLT_IsoMu24(lr_p4().pt(), lr_p4().eta()) * scale_1fb()
	: getTriggerEfficiency_HLT_Ele27_WP80(lr_p4().pt(), lr_p4().eta()) * scale_1fb();

      
      Fill1F(Met_all, metu, scale1fb);
      Fill1F(MetPhi, metphiu, scale1fb);
      Fill1F(JetHT, sum_jetPt, scale1fb);

      //// Later ///    GenRandomW:   //////			
      TRandom1 rand; 
      rand.SetSeed(0);

      float M_W = 80.2;
      
      ///////////////// fake w initialization  ///////////////////
      float fakeW_eta = -4;	
      float fakeW_theta = -4;
      float fakeW_phi = -4;
      float fakeW_pt  = -1;
      Float_t fakeW_px = 0, fakeW_py = 0 , fakeW_pz = 0;

      ///////////////// fake lepton initialization  ///////////////////
      Float_t fakeLep_px = 0, fakeLep_py = 0 , fakeLep_pz = 0;
      float fakeLep_p  = M_W/2;		// Later can apply a width here
      float fakeLep_eta = -4;
      float fakeLep_theta = -4;
      float fakeLep_phi = -4;

      LorentzVec lf_p4;		               // lf_p4 == fakeLep_p4;
      TLorentzVector fakeLepP4, fakeNeuP4;     // fakeLep_p4 in TLorentzVector

      do{
	//fakeW_p = 30; // WP->GetRandom();
	fakeW_pt = RWPT->GetRandom();          // Get the W_pt value from histogram
        Fill1F(TEST, fakeW_pt, scale1fb);

	fakeW_phi = 2*3.141562653589793*(rand.Rndm()-0.5);  // random W_phi, -pi to pi, flat 
	fakeW_eta = 6*rand.Rndm() - 3 ;			    // eta range from -3 to 3
	fakeW_theta = 2*atan(exp(-fakeW_eta));		    // theta determined by eta
      
        // Fill1F(TEST2, fakeW_theta, scale1fb);
	
	fakeW_px = fakeW_pt * cos(fakeW_phi);		    // simple calculation 
	fakeW_py = fakeW_pt * sin(fakeW_phi);
	fakeW_pz = fakeW_pt / tan(fakeW_theta);

	float fakeW_p2 = fakeW_pt*fakeW_pt + fakeW_pz*fakeW_pz;  // momentum square of the fake W

	Fill1F(FW_P, sqrt(fakeW_p2), scale1fb);             // momentum of the fake W

	/// Generate W decay at it's rest frame
	float r1 = rand.Rndm();
	float r2 = rand.Rndm();

	fakeLep_theta = acos(2*r2-1);
	fakeLep_eta = - log(0.5*fakeLep_theta);
	fakeLep_phi = 2*3.14156265359*(r1-0.5);

	fakeLep_px = fakeLep_p* sin(fakeLep_theta)* cos(fakeLep_phi);
	fakeLep_py = fakeLep_p* sin(fakeLep_theta)* sin(fakeLep_phi);
	fakeLep_pz = fakeLep_p* cos(fakeLep_theta);

	// float v2 = 1/((M_W/(fakeW_p))*(M_W/(fakeW_p))+1);
	// float gamma = 1/(sqrt(1-v2));
	float gammaM = sqrt(fakeW_p2 + M_W*M_W); // gammaM == gamma times Mass of W
	// if((gammaM/ M_W)<1)  cout << "Hey! you have a bug!!\n";

	// Initialize lep and nu at W rest frame
	fakeLepP4.SetPxPyPzE( fakeLep_px,  fakeLep_py,  fakeLep_pz, fakeLep_p);
	fakeNeuP4.SetPxPyPzE(-fakeLep_px, -fakeLep_py, -fakeLep_pz, fakeLep_p);

	// Boost from W rest frame to Lab frame as W_p / gamma*Mass == v 
	fakeLepP4.Boost(fakeW_px / gammaM, fakeW_py / gammaM, fakeW_pz / gammaM );
	fakeNeuP4.Boost(fakeW_px / gammaM, fakeW_py / gammaM, fakeW_pz / gammaM );

      } while(fakeLepP4.Eta() > 2.1 && fakeLepP4.Pt() > 30);
      
      // cout << "Sanity check: M_W = (fakeLepP4+fakeNeuP4).M() = " << (fakeLepP4+fakeNeuP4).M() << "  p_W = (fakeLepP4 + fakeNeuP4).P() = " << (fakeLepP4+fakeNeuP4).P() << endl;
      // cout << "original vs boosted: " 
      // 	   << "vector: " << fakeW_px / M_W << "  " << fakeW_py / M_W << "  " <<  fakeW_pz / M_W << endl
      //      << "org fakeLep: " << fakeLep_px << "  " <<  fakeLep_py << "  " << fakeLep_pz << endl
      // 	   << "boosted fakeLep: " << fakeLepP4.Px() << "  " << fakeLepP4.Py() << "  " << fakeLepP4.Pz() << endl << endl;

      lf_p4.SetPxPyPzE(fakeLepP4.Px(), fakeLepP4.Py(), fakeLepP4.Pz(), fakeLepP4.P());
      // lf_p4.SetPxPyPzE(fakeLep_px, fakeLep_py, fakeLep_pz, fakeLep_p);

      float nux =  fakeNeuP4.Px();
      float nuy =  fakeNeuP4.Py();

      float fakeW_recWpx = nux + fakeLepP4.Px();
      float fakeW_recWpy = nuy + fakeLepP4.Py();

      float fakeW_recWpt = sqrt(fakeW_recWpx*fakeW_recWpx + fakeW_recWpy*fakeW_recWpy);

      Fill1F(TEST3, fakeW_recWpt, scale1fb);

      float M_ll = (lr_p4() + lf_p4).M();
      Fill1F(InvM_lep, M_ll , scale1fb);     

      Fill1F(FLep_PT, lf_p4.pt(), scale1fb);

      float diffPhi = fabs(lr_p4().Phi() - metPhi());
      if(diffPhi >  3.1415926) diffPhi = 6.2831853 - diffPhi;

      Fill1F(DiffPhi, diffPhi , scale1fb);

      // Plot PT of real W
      float wr_px = metx + lr_p4().Px();
      float wr_py = mety + lr_p4().Py();
      float wr_pt = sqrt(wr_px* wr_px  + wr_py*wr_py);

      // bool under80 = wr_pt < 80 ;
      // bool under65 = wr_pt < 65 ;
      // bool under50 = wr_pt < 50 ;

      // if(under80&& !isRealData())    scale1fb = (abs(lr_id()) == 11) ? scale1fb*1.1 : scale1fb*0.9;
      // if(under65&& !isRealData())    scale1fb = (abs(lr_id()) == 11) ? scale1fb*0.9 : scale1fb*0.9;
      // if(under50&& !isRealData())    scale1fb = (abs(lr_id()) == 11) ? scale1fb*0.8 : scale1fb*0.89;
      
      Fill1F(W_PT, wr_pt, scale1fb);

      // if( cmet < 50 ) continue;
      float diffPhic = fabs(lr_p4().phi() - cmetphi);      
      if(diffPhic >  3.1415926) diffPhic = 6.2831853 - diffPhic;
      Fill1F(DiffPhic, diffPhic , scale1fb);

      float fmetx = metx + nux;	                 // fmet == met after adding fake neutrinos
      float fmety = mety + nuy;

      float fmet = sqrt( fmetx*fmetx + fmety*fmety );
      float fmetphi = atan2( fmety , fmetx );

      Fill1F(FMet, fmet, scale1fb);
      Fill1F(FMetPhi, fmetphi, scale1fb);


      Fill1F(Lep_PT, lr_p4().pt(), scale1fb);

      if(njets() == 0) cerr << "nJets Branch 0 !! \n";
      // if(njets() > 9)
      // 	cout << "nJets = " << njets() << " at eventNumber: " << eventNumber() << ", runNumber: "
      // 	     << runNumber()  << ", lumiBlock: " << lumiBlock() << endl;
 
      // Fill the jet mutiplicity to the histogram
      Fill1F(JetMult_a, njets(), scale1fb);
      Fill1F(JetMult_b, nbTag(), scale1fb);

      float mt2 = MT2(fmet, fmetphi, lr_p4() , lf_p4);

      Fill1F(MT2_hist, mt2, scale1fb); 
      h2_nozero->Fill(mt2, scale1fb); 

      float mt  = sqrt(2*lr_p4().Pt()* metu*(1 - cos(lr_p4().Phi() - metphiu)));
      float mtc = sqrt(2*lr_p4().Pt()* cmet *(1 - cos(lr_p4().Phi() - cmetphi )));

      Fill1F(MT_hist,  mt,  scale1fb); 
      Fill1F(MTc_hist, mtc, scale1fb); 

      float wr_cosTheta = mtc/80.2;			// assume W mass at 80, calculate cosTheta from M_T/80.2
      float wr_p = wr_pt / wr_cosTheta;			// calculate W momentum from W_pt by cosTheta

      Fill1F(W_P, wr_p, scale1fb);
      
      Fill1F(LepEta, lr_p4().eta(), scale1fb);
      Fill1F(LepPhi, lr_p4().phi(), scale1fb);
      Fill2F(EtaPhi, lr_p4().eta(), lr_p4().phi(), scale1fb);

      Fill1F(MetPhiCor, cmetphi, scale1fb);

      Fill1F(JPt1, jets_p4().at(j1idx()).Pt(), scale1fb);
      Fill1F(JPt2, jets_p4().at(j2idx()).Pt(), scale1fb);

      float jet_dphi = fabs(jets_p4().at(j1idx()).phi() - jets_p4().at(j2idx()).phi());      
      if(jet_dphi >  3.1415926)  jet_dphi = 6.2831853 - jet_dphi;

      Fill1F(Jdphi, jet_dphi, scale1fb);

      float vecsumPt_h2jet = (jets_p4().at(j1idx()) + jets_p4().at(j2idx())).Pt();

      Fill1F(SJPt, vecsumPt_h2jet, scale1fb);

      if(abs(lr_id()) == 11 ){
	Fill1F(Jeta_el, njets(), scale1fb); 
	Fill1F(Jetb_el, nbTag(), scale1fb); 
	Fill1F(MT2_el, mt2, scale1fb); 
	Fill1F(MT_el, mt, scale1fb); 
	Fill1F(Met_el, metu, scale1fb); 
	Fill1F(MetPhi_el, metphiu, scale1fb); 
	Fill1F(ElEta, lr_p4().eta(), scale1fb);
	Fill1F(ElPhi, lr_p4().phi(), scale1fb);
	Fill1F(El_PT, lr_p4().pt(), scale1fb);
	Fill1F(W_PT_el, wr_pt, scale1fb);
	Fill1F(W_P_el, wr_p, scale1fb);
	Fill1F(DiffPhi_e, diffPhic , scale1fb);
	Fill2F(EtaPhi_el, lr_p4().eta(), lr_p4().phi(), scale1fb);
	Fill1F(JPt1_el, jets_p4().at(j1idx()).Pt(), scale1fb);
	Fill1F(JPt2_el, jets_p4().at(j2idx()).Pt(), scale1fb);
	Fill1F(Jdphi_el, jet_dphi, scale1fb);
	Fill1F(SJPt_el, vecsumPt_h2jet, scale1fb);
	Fill1F(JetHT_el, sum_jetPt, scale1fb);

      }
      if(abs(lr_id()) == 13 ){
	Fill1F(Jeta_mu, njets(), scale1fb);
	Fill1F(Jetb_mu, nbTag(), scale1fb);
	Fill1F(MT2_mu, mt2, scale1fb); 
	Fill1F(MT_mu, mt, scale1fb); 
	Fill1F(Met_mu, metu, scale1fb); 
	Fill1F(MetPhi_mu, metphiu, scale1fb); 
	Fill1F(MuEta, lr_p4().eta(), scale1fb);
	Fill1F(MuPhi, lr_p4().phi(), scale1fb);
	Fill1F(Mu_PT, lr_p4().pt(), scale1fb);
	Fill1F(W_PT_mu, wr_pt, scale1fb);
	Fill1F(W_P_mu, wr_p, scale1fb);
	Fill1F(DiffPhi_m, diffPhic , scale1fb);
	Fill2F(EtaPhi_mu, lr_p4().eta(), lr_p4().phi(), scale1fb);
	Fill1F(JPt1_mu, jets_p4().at(j1idx()).Pt(), scale1fb);
	Fill1F(JPt2_mu, jets_p4().at(j2idx()).Pt(), scale1fb);
	Fill1F(Jdphi_mu, jet_dphi, scale1fb);
	Fill1F(SJPt_mu, vecsumPt_h2jet, scale1fb);
	Fill1F(JetHT_mu, sum_jetPt, scale1fb);

      }

    } //loop over events in the current file

    file_count++;

    // Clean Up
    delete tree;
    file->Close();
    delete file;
  } //file loop
  
  if ( nEventsChain != nEventsTotal ) {
    cout << Form( "ERROR: number of events from files (%d) is not equal to total number of events (%d)", nEventsChain, nEventsTotal ) << endl;
  }

  cout << "\nNumber of total Events: " << nEventsTotal 
       //<< " at 20GeV cut, " << nttbarEvents-metCut30 
       <<endl <<endl;
  cout << "For the data samples: "
       << "nDuplicates: " << nDuplicates
       //<< "\n# notGoodRun: " << notGoodRun 
       //<< "\n# goodRun: " << goodRun
       //<< "\n# noGoodVtx: " << noGoodVtx
       << endl;

  char* suffix = "40cmetg1j0b30lpt90HTTTI_recW_13";
  TFile* fout = new TFile(Form("./hists/hists%s%s.root",sampleType, suffix),"RECREATE");

  MT2_hist ->Write();
  MT2_el   ->Write();
  MT2_mu   ->Write();
  MT_hist  ->Write();
  MTc_hist ->Write();
  MT_el    ->Write();
  MT_mu    ->Write();
  MT2_Cut50->Write();
  h2_nozero->Write();
  Met_all  ->Write();
  Met_el   ->Write();
  Met_mu   ->Write();
  MetPhi   ->Write();
  MetPhi_el->Write();
  MetPhi_mu->Write();
  MetPhiCor->Write();
  FMet     ->Write();
  FMetPhi  ->Write();
  InvM_lep ->Write();
  JetMult_a->Write();
  JetMult_b->Write();
  Jeta_el  ->Write();
  Jeta_mu  ->Write();
  Jetb_el  ->Write();
  Jetb_mu  ->Write();
  InvM_lep ->Write();
  FW_P     ->Write();
  W_P      ->Write();
  W_PT     ->Write();
  Lep_PT   ->Write();
  FLep_PT  ->Write();
  Mu_PT    ->Write();
  El_PT    ->Write();
  W_PT_el  ->Write();
  W_PT_mu  ->Write();
  W_P_el   ->Write();
  W_P_mu   ->Write();
  LepEta   ->Write();
  MuEta    ->Write();
  ElEta    ->Write();
  DiffPhi  ->Write();
  DiffPhi_e->Write();
  DiffPhi_m->Write();
  DiffPhic ->Write();
  EtaPhi   ->Write();
  EtaPhi_el->Write();
  EtaPhi_mu->Write();
  LepPhi   ->Write();
  ElPhi    ->Write();
  MuPhi    ->Write();
  JPt1     ->Write();
  JPt1_el  ->Write();
  JPt1_mu  ->Write();
  JPt2     ->Write();
  JPt2_el  ->Write();
  JPt2_mu  ->Write();
  Jdphi    ->Write();
  Jdphi_el ->Write();
  Jdphi_mu ->Write();
  SJPt     ->Write();
  SJPt_el  ->Write();
  SJPt_mu  ->Write();
  JetHT    ->Write();
  JetHT_el ->Write();
  JetHT_mu ->Write();

  TEST	   ->Write();
  TEST2    ->Write();
  TEST3    ->Write();


  fout->Close();

  					   
  // return
  bmark->Stop("benchmark");
  cout << endl;
  cout << nEventsTotal << "Events Processed" << endl;
  cout << "------------------------------" << endl;
  cout << "CPU  Time:   " << Form( "%.01f", bmark->GetCpuTime("benchmark")  ) << endl;
  cout << "Real Time:   " << Form( "%.01f", bmark->GetRealTime("benchmark") ) << endl;
  cout << endl;
  delete bmark;

  return 0;
 }

