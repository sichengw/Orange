#include "TStyle.h"
#include "TFile.h"
#include "dataMCplotMaker.h"
#include "TH2F.h"

//int mergeHist(char* datafile, char* mcfile)
int mergeOrange_2()
{
  // Get File 
  TFile* data = new TFile("/home/users/sicheng/build/hists/hists_odata_40cmetg1j0b30lpt90HTTTI_recW_13.root");
  TFile* mcwj = new TFile("/home/users/sicheng/build/hists/hists_omcwj_40cmetg1j0b30lpt90HTTTI_recW_13.root");
  TFile* mctt = new TFile("/home/users/sicheng/build/hists/hists_omctt_40cmetg1j0b30lpt90HTTTI_recW_13.root");

  // Param Defining
  float lumi = 19.5;		// for AB files

  // Hists Getting
  
  TH1F* data_MT2   = (TH1F*) data->Get("MT2"); 
  TH1F* mcwj_MT2   = (TH1F*) mcwj->Get("MT2"); 
  TH1F* mctt_MT2   = (TH1F*) mctt->Get("MT2"); 

  TH1F* data_MT2e  = (TH1F*) data->Get("MT2e"); 
  TH1F* mcwj_MT2e  = (TH1F*) mcwj->Get("MT2e"); 
  TH1F* mctt_MT2e  = (TH1F*) mctt->Get("MT2e"); 

  TH1F* data_MT2m  = (TH1F*) data->Get("MT2m"); 
  TH1F* mcwj_MT2m  = (TH1F*) mcwj->Get("MT2m"); 
  TH1F* mctt_MT2m  = (TH1F*) mctt->Get("MT2m"); 

  
  TH1F* data_hnz   = (TH1F*) data->Get("hnz"); 
  TH1F* mcwj_hnz   = (TH1F*) mcwj->Get("hnz"); 
  TH1F* mctt_hnz   = (TH1F*) mctt->Get("hnz"); 

  TH1F* data_jeta  = (TH1F*) data->Get("jeta"); 
  TH1F* mcwj_jeta  = (TH1F*) mcwj->Get("jeta"); 
  TH1F* mctt_jeta  = (TH1F*) mctt->Get("jeta"); 
		   
  TH1F* data_jetb  = (TH1F*) data->Get("jetb"); 
  TH1F* mcwj_jetb  = (TH1F*) mcwj->Get("jetb"); 
  TH1F* mctt_jetb  = (TH1F*) mctt->Get("jetb"); 

  TH1F* data_jae   = (TH1F*) data->Get("jae"); 
  TH1F* mcwj_jae   = (TH1F*) mcwj->Get("jae"); 
  TH1F* mctt_jae   = (TH1F*) mctt->Get("jae"); 
  TH1F* data_jbe   = (TH1F*) data->Get("jbe"); 
  TH1F* mcwj_jbe   = (TH1F*) mcwj->Get("jbe"); 
  TH1F* mctt_jbe   = (TH1F*) mctt->Get("jbe"); 
  TH1F* data_jam   = (TH1F*) data->Get("jam"); 
  TH1F* mcwj_jam   = (TH1F*) mcwj->Get("jam"); 
  TH1F* mctt_jam   = (TH1F*) mctt->Get("jam"); 
  TH1F* data_jbm   = (TH1F*) data->Get("jbm"); 
  TH1F* mcwj_jbm   = (TH1F*) mcwj->Get("jbm"); 
  TH1F* mctt_jbm   = (TH1F*) mctt->Get("jbm"); 


  TH1F* data_Mll   = (TH1F*) data->Get("invm"); 
  TH1F* mcwj_Mll   = (TH1F*) mcwj->Get("invm"); 
  TH1F* mctt_Mll   = (TH1F*) mctt->Get("invm"); 

  TH1F* data_cmet  = (TH1F*) data->Get("met"); 
  TH1F* mcwj_cmet  = (TH1F*) mcwj->Get("met"); 
  TH1F* mctt_cmet  = (TH1F*) mctt->Get("met"); 

  TH1F* data_mete  = (TH1F*) data->Get("mete"); 
  TH1F* mcwj_mete  = (TH1F*) mcwj->Get("mete"); 
  TH1F* mctt_mete  = (TH1F*) mctt->Get("mete"); 
  TH1F* data_metm  = (TH1F*) data->Get("metm"); 
  TH1F* mcwj_metm  = (TH1F*) mcwj->Get("metm"); 
  TH1F* mctt_metm  = (TH1F*) mctt->Get("metm"); 

  TH1F* data_fmet  = (TH1F*) data->Get("fmet"); 
  TH1F* mcwj_fmet  = (TH1F*) mcwj->Get("fmet"); 
  TH1F* mctt_fmet  = (TH1F*) mctt->Get("fmet"); 

  TH1F* data_fmep  = (TH1F*) data->Get("fmetphi"); 
  TH1F* mcwj_fmep  = (TH1F*) mcwj->Get("fmetphi"); 
  TH1F* mctt_fmep  = (TH1F*) mctt->Get("fmetphi"); 

  TH1F* data_mpe   = (TH1F*) data->Get("metphie"); 
  TH1F* mcwj_mpe   = (TH1F*) mcwj->Get("metphie"); 
  TH1F* mctt_mpe   = (TH1F*) mctt->Get("metphie"); 
  TH1F* data_mpm   = (TH1F*) data->Get("metphim"); 
  TH1F* mcwj_mpm   = (TH1F*) mcwj->Get("metphim"); 
  TH1F* mctt_mpm   = (TH1F*) mctt->Get("metphim"); 

  TH1F* data_mt    = (TH1F*) data->Get("MT"); 
  TH1F* mcwj_mt    = (TH1F*) mcwj->Get("MT"); 
  TH1F* mctt_mt    = (TH1F*) mctt->Get("MT"); 
  TH1F* data_mtc   = (TH1F*) data->Get("MTc"); 
  TH1F* mcwj_mtc   = (TH1F*) mcwj->Get("MTc"); 
  TH1F* mctt_mtc   = (TH1F*) mctt->Get("MTc"); 

  TH1F* data_mte   = (TH1F*) data->Get("MTe"); 
  TH1F* mcwj_mte   = (TH1F*) mcwj->Get("MTe"); 
  TH1F* mctt_mte   = (TH1F*) mctt->Get("MTe"); 
  TH1F* data_mtm   = (TH1F*) data->Get("MTm"); 
  TH1F* mcwj_mtm   = (TH1F*) mcwj->Get("MTm"); 
  TH1F* mctt_mtm   = (TH1F*) mctt->Get("MTm"); 

  TH1F* data_wpt   = (TH1F*) data->Get("rwpt"); 
  TH1F* mcwj_wpt   = (TH1F*) mcwj->Get("rwpt"); 
  TH1F* mctt_wpt   = (TH1F*) mctt->Get("rwpt"); 
  TH1F* data_wpte  = (TH1F*) data->Get("wpte"); 
  TH1F* mcwj_wpte  = (TH1F*) mcwj->Get("wpte"); 
  TH1F* mctt_wpte  = (TH1F*) mctt->Get("wpte"); 
  TH1F* data_wptm  = (TH1F*) data->Get("wptm"); 
  TH1F* mcwj_wptm  = (TH1F*) mcwj->Get("wptm"); 
  TH1F* mctt_wptm  = (TH1F*) mctt->Get("wptm"); 
  
  TH1F* data_fwp   = (TH1F*) data->Get("fwp"); 
  TH1F* mcwj_fwp   = (TH1F*) mcwj->Get("fwp"); 
  TH1F* mctt_fwp   = (TH1F*) mctt->Get("fwp"); 
  TH1F* data_rwpe  = (TH1F*) data->Get("rwpe"); 
  TH1F* mcwj_rwpe  = (TH1F*) mcwj->Get("rwpe"); 
  TH1F* mctt_rwpe  = (TH1F*) mctt->Get("rwpe"); 
  TH1F* data_rwpm  = (TH1F*) data->Get("rwpm"); 
  TH1F* mcwj_rwpm  = (TH1F*) mcwj->Get("rwpm"); 
  TH1F* mctt_rwpm  = (TH1F*) mctt->Get("rwpm"); 
  
  TH1F* data_lpt   = (TH1F*) data->Get("lpt"); 
  TH1F* mcwj_lpt   = (TH1F*) mcwj->Get("lpt"); 
  TH1F* mctt_lpt   = (TH1F*) mctt->Get("lpt"); 
  TH1F* data_mpt   = (TH1F*) data->Get("mpt"); 
  TH1F* mcwj_mpt   = (TH1F*) mcwj->Get("mpt"); 
  TH1F* mctt_mpt   = (TH1F*) mctt->Get("mpt"); 
  TH1F* data_ept   = (TH1F*) data->Get("ept"); 
  TH1F* mcwj_ept   = (TH1F*) mcwj->Get("ept"); 
  TH1F* mctt_ept   = (TH1F*) mctt->Get("ept"); 
  TH1F* data_flpt  = (TH1F*) data->Get("flpt"); 
  TH1F* mcwj_flpt  = (TH1F*) mcwj->Get("flpt"); 
  TH1F* mctt_flpt  = (TH1F*) mctt->Get("flpt"); 

  TH1F* data_dphi  = (TH1F*) data->Get("dphi"); 
  TH1F* mcwj_dphi  = (TH1F*) mcwj->Get("dphi"); 
  TH1F* mctt_dphi  = (TH1F*) mctt->Get("dphi"); 
  TH1F* data_dpc   = (TH1F*) data->Get("dphic"); 
  TH1F* mcwj_dpc   = (TH1F*) mcwj->Get("dphic"); 
  TH1F* mctt_dpc   = (TH1F*) mctt->Get("dphic"); 
  TH1F* data_dpe   = (TH1F*) data->Get("dpe"); 
  TH1F* mcwj_dpe   = (TH1F*) mcwj->Get("dpe"); 
  TH1F* mctt_dpe   = (TH1F*) mctt->Get("dpe"); 
  TH1F* data_dpm   = (TH1F*) data->Get("dpm"); 
  TH1F* mcwj_dpm   = (TH1F*) mcwj->Get("dpm"); 
  TH1F* mctt_dpm   = (TH1F*) mctt->Get("dpm"); 

  
  TH1F* data_test  = (TH1F*) data->Get("test"); 
  TH1F* mcwj_test  = (TH1F*) mcwj->Get("test"); 
  TH1F* mctt_test  = (TH1F*) mctt->Get("test"); 
  TH1F* data_tes3  = (TH1F*) data->Get("test3"); 
  TH1F* mcwj_tes3  = (TH1F*) mcwj->Get("test3"); 
  TH1F* mctt_tes3  = (TH1F*) mctt->Get("test3"); 

  TH1F* data_etam  = (TH1F*) data->Get("muEta"); 
  TH1F* data_etae  = (TH1F*) data->Get("elEta"); 
  TH1F* mcwj_etam  = (TH1F*) mcwj->Get("muEta"); 
  TH1F* mcwj_etae  = (TH1F*) mcwj->Get("elEta"); 
  TH1F* mctt_etam  = (TH1F*) mctt->Get("muEta"); 
  TH1F* mctt_etae  = (TH1F*) mctt->Get("elEta"); 
  TH1F* data_phim  = (TH1F*) data->Get("muPhi"); 
  TH1F* data_phie  = (TH1F*) data->Get("elPhi"); 
  TH1F* mcwj_phim  = (TH1F*) mcwj->Get("muPhi"); 
  TH1F* mcwj_phie  = (TH1F*) mcwj->Get("elPhi"); 
  TH1F* mctt_phim  = (TH1F*) mctt->Get("muPhi"); 
  TH1F* mctt_phie  = (TH1F*) mctt->Get("elPhi"); 


  TH1F* data_jht   = (TH1F*) data->Get("jht"); 
  TH1F* mcwj_jht   = (TH1F*) mcwj->Get("jht"); 
  TH1F* mctt_jht   = (TH1F*) mctt->Get("jht"); 
  TH1F* data_jhte  = (TH1F*) data->Get("jhte"); 
  TH1F* mcwj_jhte  = (TH1F*) mcwj->Get("jhte"); 
  TH1F* mctt_jhte  = (TH1F*) mctt->Get("jhte"); 
  TH1F* data_jhtm  = (TH1F*) data->Get("jhtm"); 
  TH1F* mcwj_jhtm  = (TH1F*) mcwj->Get("jhtm"); 
  TH1F* mctt_jhtm  = (TH1F*) mctt->Get("jhtm"); 


  TH1F* data_jdpe  = (TH1F*) data->Get("jdpe"); 
  TH1F* mcwj_jdpe  = (TH1F*) mcwj->Get("jdpe"); 
  TH1F* mctt_jdpe  = (TH1F*) mctt->Get("jdpe"); 
  TH1F* data_jdpm  = (TH1F*) data->Get("jdpm"); 
  TH1F* mcwj_jdpm  = (TH1F*) mcwj->Get("jdpm"); 
  TH1F* mctt_jdpm  = (TH1F*) mctt->Get("jdpm"); 
  
  TH1F* data_jpt1e  = (TH1F*) data->Get("jpt1e"); 
  TH1F* mcwj_jpt1e  = (TH1F*) mcwj->Get("jpt1e"); 
  TH1F* mctt_jpt1e  = (TH1F*) mctt->Get("jpt1e"); 
  TH1F* data_jpt1m  = (TH1F*) data->Get("jpt1m"); 
  TH1F* mcwj_jpt1m  = (TH1F*) mcwj->Get("jpt1m"); 
  TH1F* mctt_jpt1m  = (TH1F*) mctt->Get("jpt1m"); 

  TH1F* data_jpt2e  = (TH1F*) data->Get("jpt2e"); 
  TH1F* mcwj_jpt2e  = (TH1F*) mcwj->Get("jpt2e"); 
  TH1F* mctt_jpt2e  = (TH1F*) mctt->Get("jpt2e"); 
  TH1F* data_jpt2m  = (TH1F*) data->Get("jpt2m"); 
  TH1F* mcwj_jpt2m  = (TH1F*) mcwj->Get("jpt2m"); 
  TH1F* mctt_jpt2m  = (TH1F*) mctt->Get("jpt2m"); 

  TH1F* data_sjpte  = (TH1F*) data->Get("sjpte"); 
  TH1F* mcwj_sjpte  = (TH1F*) mcwj->Get("sjpte"); 
  TH1F* mctt_sjpte  = (TH1F*) mctt->Get("sjpte"); 
  TH1F* data_sjptm  = (TH1F*) data->Get("sjptm"); 
  TH1F* mcwj_sjptm  = (TH1F*) mcwj->Get("sjptm"); 
  TH1F* mctt_sjptm  = (TH1F*) mctt->Get("sjptm"); 


  // MC hists Scaling to lumi
  mctt_MT2 ->Scale(lumi);        mcwj_MT2 ->Scale(lumi);
  mctt_MT2e->Scale(lumi);        mcwj_MT2e->Scale(lumi);
  mctt_MT2m->Scale(lumi);        mcwj_MT2m->Scale(lumi);
  mctt_hnz ->Scale(lumi);        mcwj_hnz ->Scale(lumi);
  mctt_jeta->Scale(lumi);        mcwj_jeta->Scale(lumi);
  mctt_jetb->Scale(lumi);        mcwj_jetb->Scale(lumi);
  mctt_jae ->Scale(lumi);        mcwj_jae ->Scale(lumi);
  mctt_jam ->Scale(lumi);        mcwj_jam ->Scale(lumi);
  mctt_jbe ->Scale(lumi);        mcwj_jbe ->Scale(lumi);
  mctt_jbm ->Scale(lumi);        mcwj_jbm ->Scale(lumi);
  mctt_Mll ->Scale(lumi);        mcwj_Mll ->Scale(lumi);
  mctt_cmet->Scale(lumi);        mcwj_cmet->Scale(lumi);
  mctt_mete->Scale(lumi);        mcwj_mete->Scale(lumi);
  mctt_metm->Scale(lumi);        mcwj_metm->Scale(lumi);
  mctt_mpe ->Scale(lumi);        mcwj_mpe ->Scale(lumi);
  mctt_mpm ->Scale(lumi);        mcwj_mpm ->Scale(lumi);
  mctt_fmep->Scale(lumi);        mcwj_fmep->Scale(lumi);
  mctt_fmet->Scale(lumi);        mcwj_fmet->Scale(lumi);
  mctt_mt  ->Scale(lumi);        mcwj_mt  ->Scale(lumi);
  mctt_mtc ->Scale(lumi);        mcwj_mtc ->Scale(lumi);
  mctt_mte ->Scale(lumi);        mcwj_mte ->Scale(lumi);
  mctt_mtm ->Scale(lumi);        mcwj_mtm ->Scale(lumi);
  mctt_fwp ->Scale(lumi);        mcwj_fwp ->Scale(lumi);
  mctt_rwpe->Scale(lumi);        mcwj_rwpe->Scale(lumi);
  mctt_rwpm->Scale(lumi);        mcwj_rwpm->Scale(lumi);
  mctt_dphi->Scale(lumi);        mcwj_dphi->Scale(lumi);
  mctt_dpc ->Scale(lumi);        mcwj_dpc ->Scale(lumi);
  mctt_dpe ->Scale(lumi);        mcwj_dpe ->Scale(lumi);
  mctt_dpm ->Scale(lumi);        mcwj_dpm ->Scale(lumi);
  mctt_wpt ->Scale(lumi);        mcwj_wpt ->Scale(lumi);
  mctt_wpte->Scale(lumi);        mcwj_wpte->Scale(lumi);
  mctt_wptm->Scale(lumi);        mcwj_wptm->Scale(lumi);
  mctt_lpt ->Scale(lumi);        mcwj_lpt ->Scale(lumi);
  mctt_flpt->Scale(lumi);        mcwj_flpt->Scale(lumi);
  mctt_mpt ->Scale(lumi);        mcwj_mpt ->Scale(lumi);
  mctt_ept ->Scale(lumi);        mcwj_ept ->Scale(lumi);
  mctt_etam->Scale(lumi);        mcwj_etam->Scale(lumi);
  mctt_etae->Scale(lumi);        mcwj_etae->Scale(lumi);
  mctt_phie->Scale(lumi);        mcwj_phie->Scale(lumi);
  mctt_phim->Scale(lumi);        mcwj_phim->Scale(lumi);
  mctt_test->Scale(lumi);        mcwj_test->Scale(lumi);
  mctt_tes3->Scale(lumi);        mcwj_tes3->Scale(lumi);
  
  mctt_jht ->Scale(lumi);        mcwj_jht ->Scale(lumi);
  mctt_jhte->Scale(lumi);        mcwj_jhte->Scale(lumi);
  mctt_jhtm->Scale(lumi);        mcwj_jhtm->Scale(lumi);
  mctt_jdpe->Scale(lumi);        mcwj_jdpe->Scale(lumi);
  mctt_jdpm->Scale(lumi);        mcwj_jdpm->Scale(lumi);
  mctt_jpt1e->Scale(lumi);        mcwj_jpt1e->Scale(lumi);
  mctt_jpt1m->Scale(lumi);        mcwj_jpt1m->Scale(lumi);
  mctt_jpt2e->Scale(lumi);        mcwj_jpt2e->Scale(lumi);
  mctt_jpt2m->Scale(lumi);        mcwj_jpt2m->Scale(lumi);
  mctt_sjpte->Scale(lumi);        mcwj_sjpte->Scale(lumi);
  mctt_sjptm->Scale(lumi);        mcwj_sjptm->Scale(lumi);


  // Optimize for MT2 0-70
  int bin70  = data_MT2->FindBin(70 );
  float scalee = (data_MT2e->Integral(0,bin70) - mctt_MT2e->Integral(0,bin70)) / mcwj_MT2e->Integral(0,bin70);
  mcwj_MT2e->Scale(scalee);
  mcwj_mte ->Scale(scalee);
  mcwj_mete->Scale(scalee);
  mcwj_ept ->Scale(scalee);
  mcwj_jae ->Scale(scalee);
  mcwj_mpe ->Scale(scalee);
  mcwj_dpe ->Scale(scalee);
  mcwj_etae->Scale(scalee);
  mcwj_jhte->Scale(scalee);
  mcwj_wpte->Scale(scalee);
  mcwj_sjpte->Scale(scalee);

  float scalem = (data_MT2m->Integral(0,bin70) - mctt_MT2m->Integral(0,bin70)) / mcwj_MT2m->Integral(0,bin70);
  mcwj_MT2m->Scale(scalem);
  mcwj_mtm ->Scale(scalem);
  mcwj_metm->Scale(scalem);
  mcwj_mpt ->Scale(scalem);
  mcwj_jam ->Scale(scalem);
  mcwj_mpm ->Scale(scalem);
  mcwj_dpm ->Scale(scalem);
  mcwj_etam->Scale(scalem);
  mcwj_jhtm->Scale(scalem);
  mcwj_wptm->Scale(scalem);
  mcwj_sjptm->Scale(scalem);
  
  // MC Stacking
  
  std::vector<TH1F*> hmcs_MT2 ;
  std::vector<TH1F*> hmcs_MT2e;
  std::vector<TH1F*> hmcs_MT2m;
  std::vector<TH1F*> hmcs_hnz ;
  std::vector<TH1F*> hmcs_jeta;
  std::vector<TH1F*> hmcs_jae ;
  std::vector<TH1F*> hmcs_jam ;
  std::vector<TH1F*> hmcs_jetb;
  std::vector<TH1F*> hmcs_jbe ;
  std::vector<TH1F*> hmcs_jbm ;
  std::vector<TH1F*> hmcs_Mll ;
  std::vector<TH1F*> hmcs_cmet;
  std::vector<TH1F*> hmcs_mete;
  std::vector<TH1F*> hmcs_metm;
  std::vector<TH1F*> hmcs_mpe ;
  std::vector<TH1F*> hmcs_mpm ;
  std::vector<TH1F*> hmcs_fmep;
  std::vector<TH1F*> hmcs_fmet;
  std::vector<TH1F*> hmcs_mt  ;
  std::vector<TH1F*> hmcs_mtc ;
  std::vector<TH1F*> hmcs_mte ;
  std::vector<TH1F*> hmcs_mtm ;
  std::vector<TH1F*> hmcs_fwp ;
  std::vector<TH1F*> hmcs_rwpe;
  std::vector<TH1F*> hmcs_rwpm;
  std::vector<TH1F*> hmcs_dphi;
  std::vector<TH1F*> hmcs_dpc ;
  std::vector<TH1F*> hmcs_dpe ;
  std::vector<TH1F*> hmcs_dpm ;
  std::vector<TH1F*> hmcs_wpt ;
  std::vector<TH1F*> hmcs_wpte;
  std::vector<TH1F*> hmcs_wptm;
  std::vector<TH1F*> hmcs_lpt ;
  std::vector<TH1F*> hmcs_flpt;
  std::vector<TH1F*> hmcs_mpt ;
  std::vector<TH1F*> hmcs_ept ;
  std::vector<TH1F*> hmcs_etam;
  std::vector<TH1F*> hmcs_etae;
  std::vector<TH1F*> hmcs_phie;
  std::vector<TH1F*> hmcs_phim;
  std::vector<TH1F*> hmcs_test;
  std::vector<TH1F*> hmcs_tes3;

  std::vector<TH1F*> hmcs_jht ;
  std::vector<TH1F*> hmcs_jhte;
  std::vector<TH1F*> hmcs_jhtm;
  std::vector<TH1F*> hmcs_jdpe;
  std::vector<TH1F*> hmcs_jdpm;
  std::vector<TH1F*> hmcs_jpt1e;
  std::vector<TH1F*> hmcs_jpt1m;
  std::vector<TH1F*> hmcs_jpt2e;
  std::vector<TH1F*> hmcs_jpt2m;
  std::vector<TH1F*> hmcs_sjpte;
  std::vector<TH1F*> hmcs_sjptm;
  
  hmcs_MT2 .push_back(mctt_MT2 );	  hmcs_MT2 .push_back(mcwj_MT2 );
  hmcs_MT2e.push_back(mctt_MT2e);	  hmcs_MT2e.push_back(mcwj_MT2e);
  hmcs_MT2m.push_back(mctt_MT2m);	  hmcs_MT2m.push_back(mcwj_MT2m);
  hmcs_hnz .push_back(mctt_hnz );	  hmcs_hnz .push_back(mcwj_hnz );
  hmcs_jeta.push_back(mctt_jeta);	  hmcs_jeta.push_back(mcwj_jeta);
  hmcs_jae .push_back(mctt_jae );	  hmcs_jae .push_back(mcwj_jae );
  hmcs_jam .push_back(mctt_jam );	  hmcs_jam .push_back(mcwj_jam );
  hmcs_jetb.push_back(mctt_jetb);	  hmcs_jetb.push_back(mcwj_jetb);
  hmcs_jbe .push_back(mctt_jbe );	  hmcs_jbe .push_back(mcwj_jbe );
  hmcs_jbm .push_back(mctt_jbm );	  hmcs_jbm .push_back(mcwj_jbm );
  hmcs_Mll .push_back(mctt_Mll );	  hmcs_Mll .push_back(mcwj_Mll );
  hmcs_cmet.push_back(mctt_cmet);	  hmcs_cmet.push_back(mcwj_cmet);
  hmcs_mete.push_back(mctt_mete);	  hmcs_mete.push_back(mcwj_mete);
  hmcs_metm.push_back(mctt_metm);	  hmcs_metm.push_back(mcwj_metm);
  hmcs_mpe .push_back(mctt_mpe );	  hmcs_mpe .push_back(mcwj_mpe );
  hmcs_mpm .push_back(mctt_mpm );	  hmcs_mpm .push_back(mcwj_mpm );
  hmcs_fmep.push_back(mctt_fmep);	  hmcs_fmep.push_back(mcwj_fmep);
  hmcs_fmet.push_back(mctt_fmet);	  hmcs_fmet.push_back(mcwj_fmet);
  hmcs_mt  .push_back(mctt_mt  );	  hmcs_mt  .push_back(mcwj_mt  );
  hmcs_mtc .push_back(mctt_mtc );	  hmcs_mtc .push_back(mcwj_mtc );
  hmcs_mte .push_back(mctt_mte );	  hmcs_mte .push_back(mcwj_mte );
  hmcs_mtm .push_back(mctt_mtm );	  hmcs_mtm .push_back(mcwj_mtm );
  hmcs_fwp .push_back(mctt_fwp );	  hmcs_fwp .push_back(mcwj_fwp );
  hmcs_rwpe.push_back(mctt_rwpe);	  hmcs_rwpe.push_back(mcwj_rwpe);
  hmcs_rwpm.push_back(mctt_rwpm);	  hmcs_rwpm.push_back(mcwj_rwpm);
  hmcs_wpt .push_back(mctt_wpt );	  hmcs_wpt .push_back(mcwj_wpt );
  hmcs_wpte.push_back(mctt_wpte);	  hmcs_wpte.push_back(mcwj_wpte);
  hmcs_wptm.push_back(mctt_wptm);	  hmcs_wptm.push_back(mcwj_wptm);
  hmcs_dphi.push_back(mctt_dphi);	  hmcs_dphi.push_back(mcwj_dphi);
  hmcs_dpc .push_back(mctt_dpc );	  hmcs_dpc .push_back(mcwj_dpc );
  hmcs_dpe .push_back(mctt_dpe );	  hmcs_dpe .push_back(mcwj_dpe );
  hmcs_dpm .push_back(mctt_dpm );	  hmcs_dpm .push_back(mcwj_dpm );
  hmcs_lpt .push_back(mctt_lpt );	  hmcs_lpt .push_back(mcwj_lpt );
  hmcs_flpt.push_back(mctt_flpt);	  hmcs_flpt.push_back(mcwj_flpt);
  hmcs_mpt .push_back(mctt_mpt );	  hmcs_mpt .push_back(mcwj_mpt );
  hmcs_ept .push_back(mctt_ept );	  hmcs_ept .push_back(mcwj_ept );
  hmcs_etae.push_back(mctt_etae);	  hmcs_etae.push_back(mcwj_etae);
  hmcs_etam.push_back(mctt_etam);	  hmcs_etam.push_back(mcwj_etam);
  hmcs_phie.push_back(mctt_phie);	  hmcs_phie.push_back(mcwj_phie);
  hmcs_phim.push_back(mctt_phim);	  hmcs_phim.push_back(mcwj_phim);
  hmcs_test.push_back(mctt_test);	  hmcs_test.push_back(mcwj_test);
  hmcs_tes3.push_back(mctt_tes3);	  hmcs_tes3.push_back(mcwj_tes3);

  hmcs_jht .push_back(mctt_jht );	  hmcs_jht .push_back(mcwj_jht );
  hmcs_jhte.push_back(mctt_jhte);	  hmcs_jhte.push_back(mcwj_jhte);
  hmcs_jhtm.push_back(mctt_jhtm);	  hmcs_jhtm.push_back(mcwj_jhtm);
  hmcs_jdpe.push_back(mctt_jdpe);	  hmcs_jdpe.push_back(mcwj_jdpe);
  hmcs_jdpm.push_back(mctt_jdpm);	  hmcs_jdpm.push_back(mcwj_jdpm);
  hmcs_jpt1e.push_back(mctt_jpt1e);	  hmcs_jpt1e.push_back(mcwj_jpt1e);
  hmcs_jpt1m.push_back(mctt_jpt1m);	  hmcs_jpt1m.push_back(mcwj_jpt1m);
  hmcs_jpt2e.push_back(mctt_jpt2e);	  hmcs_jpt2e.push_back(mcwj_jpt2e);
  hmcs_jpt2m.push_back(mctt_jpt2m);	  hmcs_jpt2m.push_back(mcwj_jpt2m);
  hmcs_sjpte.push_back(mctt_sjpte);	  hmcs_sjpte.push_back(mcwj_sjpte);
  hmcs_sjptm.push_back(mctt_sjptm);	  hmcs_sjptm.push_back(mcwj_sjptm);

  
  // Titles setting
  std::vector<char*> titmcs;
  titmcs.push_back("ttbar");	  
  titmcs.push_back("Wjets");

  TH1F* empty = new TH1F("empty","empty", 60, 0, 60);

  // Hists Drawing
  // TCanvas* c1 = new TCanvas;

  char* suffix  = "40cmetge2j0b30lpt90HTTTIMT2MT2opt_recW_13";
  char* suffix2 = "cmet>40,lpt>30,>=2j,0b,HT>90,recW,MT2opt";
  char* opts = "--noOverflow --noDivisionLabel --outputName ";

  dataMCplotMaker(data_MT2 , hmcs_MT2,  titmcs, (char*)data_MT2 ->GetTitle(), suffix2, Form("--xAxisLabel MT2 %s MT2_%s",    opts, suffix) );
  dataMCplotMaker(data_MT2e, hmcs_MT2e, titmcs, (char*)data_MT2e->GetTitle(), suffix2, Form("--xAxisLabel MT2 %s MT2_el_%s", opts, suffix) );
  dataMCplotMaker(data_MT2m, hmcs_MT2m, titmcs, (char*)data_MT2m->GetTitle(), suffix2, Form("--xAxisLabel MT2 %s MT2_mu_%s", opts, suffix) );
  
  // dataMCplotMaker(data_hnz , hmcs_hnz,  titmcs, (char*)data_hnz ->GetTitle(), suffix2, Form("--isLinear --xAxisLabel MT2 %s hnz_%s", opts, suffix) );
  // dataMCplotMaker(data_mt  , hmcs_mt ,  titmcs, (char*)data_mt  ->GetTitle(), suffix2, Form("--xAxisLabel M_T %s MT_%s", opts, suffix) );
  // dataMCplotMaker(data_mtc , hmcs_mtc,  titmcs, (char*)data_mtc ->GetTitle(), suffix2, Form("--xAxisLabel M_T %s MTc_%s", opts, suffix) );
  // dataMCplotMaker(data_mte , hmcs_mte,  titmcs, (char*)data_mte ->GetTitle(), suffix2, Form("--xAxisLabel M_T %s MTe_%s", opts, suffix) );
  // dataMCplotMaker(data_mtm , hmcs_mtm,  titmcs, (char*)data_mtm ->GetTitle(), suffix2, Form("--xAxisLabel M_T %s MTm_%s", opts, suffix) );
  // dataMCplotMaker(data_jeta, hmcs_jeta, titmcs, (char*)data_jeta->GetTitle(), suffix2, Form("--xAxisLabel nJets --noXaxisUnit %s jeta_%s", opts, suffix) );
  // dataMCplotMaker(data_jae , hmcs_jae,  titmcs, (char*)data_jae ->GetTitle(), suffix2, Form("--xAxisLabel nJets --noXaxisUnit %s jae_%s", opts, suffix) );
  // dataMCplotMaker(data_jam , hmcs_jam,  titmcs, (char*)data_jam ->GetTitle(), suffix2, Form("--xAxisLabel nJets --noXaxisUnit %s jam_%s", opts, suffix) );
  // dataMCplotMaker(data_cmet, hmcs_cmet, titmcs, (char*)data_cmet->GetTitle(), suffix2, Form("-xAxisLabel MET %s cmet_%s", opts, suffix) );
  // dataMCplotMaker(data_mete, hmcs_mete, titmcs, (char*)data_mete->GetTitle(), suffix2, Form("-xAxisLabel MET %s mete_%s", opts, suffix) );
  // dataMCplotMaker(data_metm, hmcs_metm, titmcs, (char*)data_metm->GetTitle(), suffix2, Form("-xAxisLabel MET %s metm_%s", opts, suffix) );
  // dataMCplotMaker(data_mpe , hmcs_mpe , titmcs, (char*)data_mpe ->GetTitle(), suffix2, Form("--isLinear -xAxisLabel MET %s mpe_%s", opts, suffix) );
  // dataMCplotMaker(data_mpm , hmcs_mpm , titmcs, (char*)data_mpm ->GetTitle(), suffix2, Form("--isLinear -xAxisLabel MET %s mpm_%s", opts, suffix) );
  // dataMCplotMaker(data_fwp,  hmcs_fwp,  titmcs, (char*)data_fwp ->GetTitle(), suffix2, Form("-xAxisLabel p^{W} %s fwp_%s", opts, suffix) );
  // dataMCplotMaker(data_rwpe, hmcs_rwpe, titmcs, (char*)data_rwpe->GetTitle(), suffix2, Form("--isLinear -xAxisLabel W_{p} %s rwpe_%s", opts, suffix) );
  // dataMCplotMaker(data_rwpm, hmcs_rwpm, titmcs, (char*)data_rwpm->GetTitle(), suffix2, Form("--isLinear -xAxisLabel W_{p} %s rwpm_%s", opts, suffix) );
  // dataMCplotMaker(data_wpt , hmcs_wpt,  titmcs, (char*)data_wpt ->GetTitle(), suffix2, Form("--isLinear -xAxisLabel W_{pt} %s wpt_l_%s", opts, suffix) );
  // dataMCplotMaker(data_wpte, hmcs_wpte, titmcs, (char*)data_wpte->GetTitle(), suffix2, Form("--isLinear -xAxisLabel W_{pt} %s wpte_%s", opts, suffix) );
  // dataMCplotMaker(data_wptm, hmcs_wptm, titmcs, (char*)data_wptm->GetTitle(), suffix2, Form("--isLinear -xAxisLabel W_{pt} %s wptm_%s", opts, suffix) );
  // dataMCplotMaker(data_dphi, hmcs_dphi, titmcs, (char*)data_dphi->GetTitle(), suffix2, Form("--isLinear -xAxisLabel W_{pt} %s dphi_%s", opts, suffix) );
  // dataMCplotMaker(data_dpc , hmcs_dpc , titmcs, (char*)data_dpc ->GetTitle(), suffix2, Form("--isLinear -xAxisLabel W_{pt} %s dphic_%s", opts, suffix) );
  // dataMCplotMaker(data_dpe , hmcs_dpe , titmcs, (char*)data_dpe ->GetTitle(), suffix2, Form("--isLinear -xAxisLabel W_{pt} %s dpe_%s", opts, suffix) );
  // dataMCplotMaker(data_dpm , hmcs_dpm , titmcs, (char*)data_dpm ->GetTitle(), suffix2, Form("--isLinear -xAxisLabel W_{pt} %s dpm_%s", opts, suffix) );
  // dataMCplotMaker(data_lpt , hmcs_lpt , titmcs, (char*)data_lpt ->GetTitle(), suffix2, Form("-xAxisLabel l_{pt} %s lpt_%s", opts, suffix) );
  // dataMCplotMaker(data_flpt, hmcs_flpt, titmcs, (char*)data_flpt->GetTitle(), suffix2, Form("-xAxisLabel p_{T}^{fLep} %s flpt_%s", opts, suffix) );
  // dataMCplotMaker(data_ept , hmcs_ept , titmcs, (char*)data_ept ->GetTitle(), suffix2, Form("-xAxisLabel p_T^e %s ept_%s", opts, suffix) );
  // dataMCplotMaker(data_mpt , hmcs_mpt , titmcs, (char*)data_mpt ->GetTitle(), suffix2, Form("-xAxisLabel p_T^#mu %s mpt_%s", opts, suffix) );
  // dataMCplotMaker(data_etae, hmcs_etae, titmcs, (char*)data_etae->GetTitle(), suffix2, Form("-xAxisLabel Eta --isLinear --noXaxisUnit %s etae_l_%s", opts, suffix) );
  // dataMCplotMaker(data_etam, hmcs_etam, titmcs, (char*)data_etam->GetTitle(), suffix2, Form("-xAxisLabel Eta --isLinear --noXaxisUnit %s etam_l_%s", opts, suffix) );
  // dataMCplotMaker(data_phie, hmcs_phie, titmcs, (char*)data_phie->GetTitle(), suffix2, Form("-xAxisLabel Phi --isLinear --noXaxisUnit %s phie_l_%s", opts, suffix) );
  // dataMCplotMaker(data_phim, hmcs_phim, titmcs, (char*)data_phim->GetTitle(), suffix2, Form("-xAxisLabel Phi --isLinear --noXaxisUnit %s phim_l_%s", opts, suffix) );
  // dataMCplotMaker(data_jht , hmcs_jht , titmcs, (char*)data_jht ->GetTitle(), suffix2, Form("-xAxisLabel Jdp --isLinear --noXaxisUnit %s jht_l_%s", opts, suffix) );
  // dataMCplotMaker(data_jhte, hmcs_jhte, titmcs, (char*)data_jhte->GetTitle(), suffix2, Form("-xAxisLabel Jdp --isLinear --noXaxisUnit %s jhte_l_%s", opts, suffix) );
  // dataMCplotMaker(data_jhtm, hmcs_jhtm, titmcs, (char*)data_jhtm->GetTitle(), suffix2, Form("-xAxisLabel Jdp --isLinear --noXaxisUnit %s jhtm_l_%s", opts, suffix) );
  // dataMCplotMaker(data_jdpe, hmcs_jdpe, titmcs, (char*)data_jdpe->GetTitle(), suffix2, Form("-xAxisLabel Jdp --isLinear --noXaxisUnit %s jdpe_l_%s", opts, suffix) );
  // dataMCplotMaker(data_jdpm, hmcs_jdpm, titmcs, (char*)data_jdpm->GetTitle(), suffix2, Form("-xAxisLabel Jdp --isLinear --noXaxisUnit %s jdpm_l_%s", opts, suffix) );
  // dataMCplotMaker(data_jpt1e,hmcs_jpt1e,titmcs, (char*)data_jpt1e->GetTitle(),suffix2, Form("-xAxisLabel Jpt1 --isLinear --noXaxisUnit %s jpt1e_l_%s", opts, suffix) );
  // dataMCplotMaker(data_jpt1m,hmcs_jpt1m,titmcs, (char*)data_jpt1m->GetTitle(),suffix2, Form("-xAxisLabel Jpt1 --isLinear --noXaxisUnit %s jpt1m_l_%s", opts, suffix) );
  // dataMCplotMaker(data_jpt2e,hmcs_jpt2e,titmcs, (char*)data_jpt2e->GetTitle(),suffix2, Form("-xAxisLabel Jpt2 --isLinear --noXaxisUnit %s jpt2e_l_%s", opts, suffix) );
  // dataMCplotMaker(data_jpt2m,hmcs_jpt2m,titmcs, (char*)data_jpt2m->GetTitle(),suffix2, Form("-xAxisLabel Jpt2 --isLinear --noXaxisUnit %s jpt2m_l_%s", opts, suffix) );
  // dataMCplotMaker(data_sjpte,hmcs_sjpte,titmcs, (char*)data_sjpte->GetTitle(),suffix2, Form("-xAxisLabel Sjpt --isLinear --noXaxisUnit %s sjpte_l_%s", opts, suffix) );
  // dataMCplotMaker(data_sjptm,hmcs_sjptm,titmcs, (char*)data_sjptm->GetTitle(),suffix2, Form("-xAxisLabel Sjpt --isLinear --noXaxisUnit %s sjptm_l_%s", opts, suffix) );

  // std::vector<TH1F*> singdata;  
  // singdata.push_back(data_test);
  // std::vector<TH1F*> singdata2;  
  // singdata2.push_back(data_tes3);
  // std::vector<char*> singtitle;  
  // singtitle.push_back("data");

  // dataMCplotMaker(empty, singdata, singtitle, "fake W_PT", "", Form("--isLinear -xAxisLabel p_{T}^{W} %s fwpt_%s", opts, suffix) );
  // dataMCplotMaker(empty, singdata2, singtitle, "fake reconstruced W_PT", "", Form("--isLinear -xAxisLabel p_{T}^{W} %s rfwpt_%s", opts, suffix) );
 
  // dataMCplotMaker(empty, hmcs_test, titmcs, "fake W_PT", suffix, Form("--isLinear -xAxisLabel p_{T}^{W} %s fwpt_%s", opts, suffix) );
  // dataMCplotMaker(empty, hmcs_tes3, titmcs, "fake W #theta", suffix, Form("--isLinear -xAxisLabel #eta_{W} %s fwtheta_%s", opts, suffix) );
  // dataMCplotMaker(data_fmet, hmcs_fmet, titmcs, (char*)data_fmet->GetTitle(), suffix2, Form("-xAxisLabel MET %s fmet_%s", opts, suffix) );

  // dataMCplotMaker(data_fmep, hmcs_fmep, titmcs, (char*)data_fmep->GetTitle(), suffix2, Form("--isLinear -xAxisLabel #phi MET %s fmep_%s", opts, suffix) );
  // dataMCplotMaker(data_Mll , hmcs_Mll , titmcs, (char*)data_Mll ->GetTitle(), suffix2, Form("-xAxisLabel #phi MET %s Mll_%s", opts, suffix) );


  cout << "MT2-:"  << data_MT2->Integral() << "  " << (hmcs_MT2[0]->Integral() + hmcs_MT2[1]->Integral()) << "  ("
       << data_MT2->Integral()/(hmcs_MT2[0]->Integral() + hmcs_MT2[1]->Integral()) << ") "<< endl;
  cout << "MT2e:"  << data_MT2e->Integral() << "  " << (hmcs_MT2e[0]->Integral() + hmcs_MT2e[1]->Integral()) << "  ("
       << data_MT2e->Integral()/(hmcs_MT2e[0]->Integral() + hmcs_MT2e[1]->Integral()) << ") "<< endl;
  cout << "MT2m:"  << data_MT2m->Integral() << "  " << (hmcs_MT2m[0]->Integral() + hmcs_MT2m[1]->Integral()) << "  ("
       << data_MT2m->Integral()/(hmcs_MT2m[0]->Integral() + hmcs_MT2m[1]->Integral()) << ") "<< endl;
  cout << "mtc-:"  << data_mtc->Integral() << "  " << (hmcs_mtc[0]->Integral() + hmcs_mtc[1]->Integral()) << endl;
  cout << "mtm-:"  << data_mtm->Integral() << "  " << (hmcs_mtm[0]->Integral() + hmcs_mtm[1]->Integral()) << endl;
  cout << "metm:"  << data_metm->Integral() << "  " << (hmcs_metm[0]->Integral() + hmcs_metm[1]->Integral()) << endl;
  cout << "jeta:"  << data_jeta->Integral() << "  " << (hmcs_jeta[0]->Integral() + hmcs_jeta[1]->Integral()) << endl;
  cout << "wpt-:"  << data_wpt->Integral() << "  " << (hmcs_wpt[0]->Integral() + hmcs_wpt[1]->Integral()) << endl;
  cout << "lpt-:"  << data_lpt->Integral() << "  " << (hmcs_lpt[0]->Integral() + hmcs_lpt[1]->Integral()) << endl;
  cout << "mpt-:"  << data_mpt->Integral() << "  " << (hmcs_mpt[0]->Integral() + hmcs_mpt[1]->Integral()) << endl;
  cout << "etam:"  << data_etam->Integral() << "  " << (hmcs_etam[0]->Integral() + hmcs_etam[1]->Integral()) << endl;
  cout << "etae:"  << data_etae->Integral() << "  " << (hmcs_etae[0]->Integral() + hmcs_etae[1]->Integral()) << endl;
  //cout << data_test->Integral() << "  " << (hmcs_test[0]->Integral() + hmcs_test[1]->Integral()) << endl;

  int bin80  = data_MT2->FindBin(80 );
  int bin90  = data_MT2->FindBin(90 );
  int bin100 = data_MT2->FindBin(100);
  int bin110 = data_MT2->FindBin(110);
  int bin120 = data_MT2->FindBin(120);
  int bin140 = data_MT2->FindBin(140);
  int lastbin = data_MT2->GetSize();

  cout << "MT2   80-90: "  << data_MT2->Integral(bin80 , bin90 )  << "  "
       << (hmcs_MT2[0]->Integral(bin80 , bin90 ) + hmcs_MT2[1]->Integral(bin80 , bin90 ))<< endl; 
  cout << "MT2  90-100: "  << data_MT2->Integral(bin90 +1, bin100) << "   "
       << (hmcs_MT2[0]->Integral(bin90 +1, bin100) + hmcs_MT2[1]->Integral(bin90 +1, bin100))<< endl; 
  cout << "MT2 100-110: "  << data_MT2->Integral(bin100+1, bin110) << "   "
       << (hmcs_MT2[0]->Integral(bin100+1, bin110) + hmcs_MT2[1]->Integral(bin100+1, bin110))<< endl; 
  cout << "MT2 110-120: "  << data_MT2->Integral(bin110+1, bin120) << "    "
       << (hmcs_MT2[0]->Integral(bin110+1, bin120) + hmcs_MT2[1]->Integral(bin110+1, bin120))<< endl; 
  cout << "MT2 over120: "  << data_MT2->Integral(bin120+1, lastbin) << "    "
       << (hmcs_MT2[0]->Integral(bin120+1, lastbin) + hmcs_MT2[1]->Integral(bin120+1, lastbin))<< endl; 

  cout << "MT2 >= 100: "  << data_MT2->Integral(bin100, lastbin) << "   "
       << (hmcs_MT2[0]->Integral(bin100, lastbin) + hmcs_MT2[1]->Integral(bin100, lastbin))<< endl; 
  cout << "MT2 >= 120: "  << data_MT2->Integral(bin120, lastbin) << "   "
       << (hmcs_MT2[0]->Integral(bin120, lastbin) + hmcs_MT2[1]->Integral(bin120, lastbin))<< endl; 
  cout << "MT2 >= 140: "  << data_MT2->Integral(bin140, lastbin) << "   "
       << (hmcs_MT2[0]->Integral(bin140, lastbin) + hmcs_MT2[1]->Integral(bin140, lastbin))<< endl; 

  cout << "MT2e >= 100: "  << data_MT2e->Integral(bin100, lastbin) << "   "
       << (hmcs_MT2e[0]->Integral(bin100, lastbin) + hmcs_MT2e[1]->Integral(bin100, lastbin))<< endl; 
  cout << "MT2e >= 120: "  << data_MT2e->Integral(bin120, lastbin) << "   "
       << (hmcs_MT2e[0]->Integral(bin120, lastbin) + hmcs_MT2e[1]->Integral(bin120, lastbin))<< endl; 
  cout << "MT2e >= 140: "  << data_MT2e->Integral(bin140, lastbin) << "   "
       << (hmcs_MT2e[0]->Integral(bin140, lastbin) + hmcs_MT2e[1]->Integral(bin140, lastbin))<< endl; 

  cout << "MT2m >= 100: "  << data_MT2m->Integral(bin100, lastbin) << "   "
       << (hmcs_MT2m[0]->Integral(bin100, lastbin) + hmcs_MT2m[1]->Integral(bin100, lastbin))<< endl; 
  cout << "MT2m >= 120: "  << data_MT2m->Integral(bin120, lastbin) << "   "
       << (hmcs_MT2m[0]->Integral(bin120, lastbin) + hmcs_MT2m[1]->Integral(bin120, lastbin))<< endl; 
  cout << "MT2m >= 140: "  << data_MT2m->Integral(bin140, lastbin) << "   "
       << (hmcs_MT2m[0]->Integral(bin140, lastbin) + hmcs_MT2m[1]->Integral(bin140, lastbin))<< endl; 

  // TH1F* test = (TH1F*)data_MT2m->Clone("testing");
  // test->Add(hmcs_MT2m[0], -1);
  // test->Add(hmcs_MT2m[1], -1);
  // test->Divide(data_MT2m);
  
  // test->Draw("");

  // dataMCplotMaker(data_jetb, hmcs_jetb, titmcs, "# bTag distribution", suffix,		Form("%s  jetb_40met0b30w_4" );
  // dataMCplotMaker(data_jbe, hmcs_jbe, titmcs, "# bTag distribution for e events ", suffix,   Form("%s  jbe_40met0b30w_4" );
  // dataMCplotMaker(data_jbm, hmcs_jbm, titmcs, "# bTag distribution for mu events", suffix,   Form("%s  jbm_40met0b30w_4" );
    
  return 0;
}


  // h7_da->Rebin(2);
  // h7_mc->Rebin(2);
  // gStyle->SetErrorX(0);
  // h7_da->SetMarkerStyle(20);

