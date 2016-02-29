//____________________________________________________________________________
/*!

\program gtestRewghtPiKAbs

\brief   A program to run reweighting of number of nucleons ejected when a 
         pion or kaon is absorbed in a final-state interaction.

\syntax  gtestRewghtPiKAbs -f filename [-n nev]

         where 
         [] is an optional argument
         -f specifies a GENIE event file (GHEP format)
         -n specifies the number of events to process (default: all)

\author  Nick Grant <n.grant.3 \at warwick.ac.uk>
         University of Warwick

\created 26 February 2016

\cpright Copyright (c) 2003-2015, GENIE Neutrino MC Generator Collaboration
         For the full text of the license visit http://copyright.genie-mc.org
         or see $GENIE/LICENSE
*/
//____________________________________________________________________________


#include <string>

#include <TSystem.h>
#include <TFile.h>
#include <TTree.h>
#include <TH1.h>
#include <TGraph.h>
#include <TSpline.h>
#include <fstream>

#include "EVGCore/EventRecord.h"
#include "HadronTransport/INukeUtils.h"
#include "Messenger/Messenger.h"
#include "Ntuple/NtpMCFormat.h"
#include "Ntuple/NtpMCTreeHeader.h"
#include "Ntuple/NtpMCEventRecord.h"
#include "ReWeight/GReWeightI.h"
#include "ReWeight/GSystSet.h"
#include "ReWeight/GReWeight.h"
#include "ReWeight/GReWeightNuXSecCCQE.h"
#include "ReWeight/GReWeightNuXSecCCQEvec.h"
#include "ReWeight/GReWeightNuXSecCCRES.h"
#include "ReWeight/GReWeightNuXSecNCRES.h"
#include "ReWeight/GReWeightNuXSecDIS.h"
#include "ReWeight/GReWeightNuXSecCOH.h"
#include "ReWeight/GReWeightNonResonanceBkg.h"
#include "ReWeight/GReWeightFGM.h"
#include "ReWeight/GReWeightDISNuclMod.h"
#include "ReWeight/GReWeightResonanceDecay.h"
#include "ReWeight/GReWeightFZone.h"
#include "ReWeight/GReWeightINuke.h"
#include "ReWeight/GReWeightAGKY.h"
#include "Utils/CmdLnArgParser.h"

using std::string;

using namespace genie;
using namespace genie::rew;

void GetCommandLineArgs (int argc, char ** argv);

int    gOptNEvt;
string gOptInpFilename;

//___________________________________________________________________
int main(int argc, char ** argv)
{

  GetCommandLineArgs (argc, argv);

  //units of tweaks are ns and nd themselves (they are not any kind of error or "sigma")
  double ns_tweak = -1.0;
  double nd_tweak = 1.0;

  std::string splineFile = "PiKAbsSplines.root";

  // open the ROOT file and get the TTree & its header
  TTree *           tree = 0;
  NtpMCTreeHeader * thdr = 0;

  TFile file(gOptInpFilename.c_str(),"READ");

  tree = dynamic_cast <TTree *>           ( file.Get("gtree")  );
  thdr = dynamic_cast <NtpMCTreeHeader *> ( file.Get("header") );

  if(!tree) return 1;

  //LOG("test", pNOTICE) << "Input tree header: " << *thdr;

  int nev = (gOptNEvt > 0) ?
        TMath::Min(gOptNEvt, (int)tree->GetEntries()) :
        (int) tree->GetEntries();

  //LOG("test", pNOTICE) << "Will process " << nev << " events";

  NtpMCEventRecord * mcrec = 0;
  tree->SetBranchAddress("gmcrec", &mcrec);

  double nsnd_weight_mean = 1.0;

  if(fabs(ns_tweak) <= 1e-6 && fabs(nd_tweak) <= 1e-6)
    {
      LOG("test", pNOTICE) << "No tweak of either ns or nd specified"; 
      return 1;
    }

  if(std::ifstream(splineFile.c_str()).good())
    remove(splineFile.c_str());

  genie::utils::intranuke::MakePiKAbsSplines(gOptInpFilename.c_str(), nev, ns_tweak, nd_tweak);

  TFile *splFile = TFile::Open(splineFile.c_str());
  splFile->cd();

  std::map<double, double> ns_weights;
  std::map<double, double> nd_weights;

  std::map<double, double>::iterator weightmapit;

  double x, y, weight;

  double nevents_ns = 0.0;
  double nevents_nd = 0.0;
  double ns_weight_mean = 0.0;
  double nd_weight_mean = 0.0;

  double ns_min = -99999.9;
  double nd_min = -99999.9;
  double ns_max = -99999.9; 
  double nd_max = -99999.9;;

  std::map<double, double> ns_with_zero_events, nd_with_zero_events;

  if(fabs(ns_tweak) > 1e-6)
    {
     TGraph *gNs = (TGraph*)splFile->Get("nsGraph");
     TSpline3 *nsSpl = (TSpline3*)splFile->Get("nsSpline");
     TSpline3 *nsSplTwk = (TSpline3*)splFile->Get("nsSplineTweaked");

     for(int p=0; p<gNs->GetN(); p++)
       {
	 gNs->GetPoint(p, x, y);

	 if(y > 0.0)
	   {
	   if(ns_min == -99999.9)
	     ns_min = x;
	   ns_max = x;
	   }

	 if(nsSplTwk->Eval(x) < 0.01 || nsSpl->Eval(x) < 0.01)
	   weight = 0.0;
	 else
	   weight = nsSplTwk->Eval(x) / nsSpl->Eval(x);
	 ns_weights.insert(std::make_pair(x, weight));
         nevents_ns += y;
         ns_weight_mean += weight * y;
       }

     //check whether there are any values of ns with zero events that are within the range of values with non-zero events
     for(int p=0; p<gNs->GetN(); p++)
       {
         gNs->GetPoint(p, x, y);
         if(y == 0 && x > ns_min && x < ns_max)
	   ns_with_zero_events.insert(std::make_pair(x, TMath::Max(0.0, nsSplTwk->Eval(x))));
       }
     
     //add one fake "event" to mean for any value of ns with zero events that is within the range of values with non-zero events
     for(weightmapit=ns_with_zero_events.begin(); weightmapit!=ns_with_zero_events.end(); ++weightmapit)
       {
	 nevents_ns++;
	 ns_weight_mean += weightmapit->second;
       }
     
     ns_weight_mean /= nevents_ns;

    }//end of if(fabs(ns_tweak) > 1e-6)

  if(fabs(nd_tweak) > 1e-6)
    {
    TGraph *gNd = (TGraph*)splFile->Get("ndGraph");
    TSpline3 *ndSpl = (TSpline3*)splFile->Get("ndSpline");
    TSpline3 *ndSplTwk = (TSpline3*)splFile->Get("ndSplineTweaked"); 

    for(int p=0; p<gNd->GetN(); p++)
      {
	gNd->GetPoint(p, x, y);

	if(y > 0.0)
	  {
	  if(nd_min == -99999.9)
	    nd_min = x;
	  nd_max = x;
	  }

	if(ndSplTwk->Eval(x) < 0.01 || ndSpl->Eval(x) < 0.01)
	  weight = 0.0;
	else
	  weight = ndSplTwk->Eval(x) / ndSpl->Eval(x);
	nd_weights.insert(std::make_pair(x, weight));
	nevents_nd += y;
	nd_weight_mean += weight * y;
      }

    //check whether there are any values of nd with zero events that are within the range of values with non-zero events                                                                                 
    for(int p=0; p<gNd->GetN(); p++)
      {
	gNd->GetPoint(p, x, y);
	if(y == 0 && x > nd_min && x < nd_max)
	  nd_with_zero_events.insert(std::make_pair(x, TMath::Max(0.0, ndSplTwk->Eval(x))));
      }    

    //add one fake "event" to mean for any value of nd with zero events that is within the range of values with non-zero events
    for(weightmapit=nd_with_zero_events.begin(); weightmapit!=nd_with_zero_events.end(); ++weightmapit)
      {
	nevents_nd++;
        nd_weight_mean += weightmapit->second;
      }

    nd_weight_mean /= nevents_nd;

    }//end of if(fabs(nd_tweak) > 1e-6)

  //calculate mean of weights when both ns and nd are tweaked
  if(fabs(ns_tweak) > 1e-6 && fabs(nd_tweak) > 1e-6)
    nsnd_weight_mean = genie::utils::intranuke::CalcPiKAbsWeightMean(gOptInpFilename.c_str(), nev, ns_with_zero_events, nd_with_zero_events);

  int nPiKAbsEvents = 0;

  double ns_weight, nd_weight, event_weight;
  double weight_mean = 0.0;

  //
  // Event loop
  //

  for(int i = 0; i < nev; i++) {
    tree->GetEntry(i);

    EventRecord & event = *(mcrec->event);

    if(i%1000 == 0)    
      std::cout<<"Doing event "<<i<<std::endl;

    int nPiKAbs = genie::utils::intranuke::GetNPiKAbs(event);

    if(nPiKAbs == 1)
      {
	nPiKAbsEvents++;

	ns_weight = 0.0;
        nd_weight = 0.0;

	if(fabs(ns_tweak) > 1e-6 && ns_weights.count(genie::utils::intranuke::GetNs(event)))
	  ns_weight = ns_weights.find(genie::utils::intranuke::GetNs(event))->second;

	if(fabs(nd_tweak) > 1e-6 && nd_weights.count(genie::utils::intranuke::GetNd(event)))
	  nd_weight = nd_weights.find(genie::utils::intranuke::GetNd(event))->second;
	
	if(fabs(ns_tweak) > 1e-6 && fabs(nd_tweak) > 1e-6)
	  event_weight = ns_weight * nd_weight / nsnd_weight_mean;
	else if(fabs(ns_tweak) > 1e-6 && fabs(nd_tweak) <= 1e-6)
          event_weight = ns_weight / ns_weight_mean;
	else if(fabs(ns_tweak) <= 1e-6 && fabs(nd_tweak) > 1e-6)
	  event_weight = nd_weight / nd_weight_mean;
	else
	  event_weight = 0.0;
	
	weight_mean += event_weight; 
      }
    
    mcrec->Clear();
  }

  //add one fake "event" for any value of ns with zero events that is within the range of values with non-zero events
  //value of nd is not defined for this "event", assume it has weight 1.0
  for(weightmapit=ns_with_zero_events.begin(); weightmapit!=ns_with_zero_events.end(); ++weightmapit)
    {
      nPiKAbsEvents++;
      if(fabs(ns_tweak) > 1e-6 && fabs(nd_tweak) > 1e-6)
	weight_mean += weightmapit->second / nsnd_weight_mean;
      else if(fabs(ns_tweak) > 1e-6 && fabs(nd_tweak) <= 1e-6)
        weight_mean += weightmapit->second / ns_weight_mean;
    }
  //add one fake "event" for any value of nd with zero events that is within the range of values with non-zero events
  //value of ns is not defined for this "event", assume it has weight 1.0
  for(weightmapit=nd_with_zero_events.begin(); weightmapit!=nd_with_zero_events.end(); ++weightmapit)
    {
      nPiKAbsEvents++;
      if(fabs(ns_tweak) > 1e-6 && fabs(nd_tweak) > 1e-6)
	weight_mean += weightmapit->second / nsnd_weight_mean;
      else if(fabs(ns_tweak) <= 1e-6 && fabs(nd_tweak) > 1e-6)
        weight_mean += weightmapit->second / nd_weight_mean;
    }
  
  file.Close();

  LOG("test", pNOTICE) << "Reweighting of number of nucleons ejected when a pion or kaon is absorbed in a FSI has been done for " << nev << " events from the file " << gOptInpFilename.c_str()<<".";

  LOG("test", pNOTICE) << "ns was tweaked by " << ns_tweak << " and nd by " << nd_tweak <<". The mean of the weights is " << weight_mean / nPiKAbsEvents<< ".";

  LOG("test", pNOTICE)  << "Done!";

  return 0;
}
//___________________________________________________________________
void GetCommandLineArgs(int argc, char ** argv)
{
  LOG("test", pINFO) << "*** Parsing command line arguments";

  CmdLnArgParser parser(argc,argv);

  // get GENIE event sample
  if( parser.OptionExists('f') ) {  
    LOG("test", pINFO) << "Reading event sample filename";
    gOptInpFilename = parser.ArgAsString('f');
  } else {
    LOG("test", pFATAL) 
      << "Unspecified input filename - Exiting";
    exit(1);
  }

  // number of events:
  if( parser.OptionExists('n') ) {  
    LOG("test", pINFO) << "Reading number of events to analyze";
    gOptNEvt = parser.ArgAsInt('n');
  } else {
    LOG("test", pINFO)
       << "Unspecified number of events to analyze - Use all";
    gOptNEvt = -1;
  }
}
//_________________________________________________________________________________
