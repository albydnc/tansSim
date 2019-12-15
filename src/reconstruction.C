// @(#)root/reconstruction
// Author: Alberto Perro 15/11/19
#include "lib/detector.h"
#include "lib/collision.h"
#include "lib/hit.h"
#include "lib/tracklet.h"
#include "TClonesArray.h"
#include "TStopwatch.h"
#include "TObject.h"
#include "TCanvas.h"
#include "TFile.h"
#include "TTree.h"
#include "TString.h"
#include "TGraphAsymmErrors.h"
#include <Riostream.h>
#include <vector>

const double minTol = 0.1, minDTheta = 0.01;
int notRecon = 0;

double findMax(TH1I& data, const double& tol);
double nearAvg(const std::vector<double>& zRecon, const double& zTrend, const double& tol);
double findZ(TClonesArray* L1Hits,TClonesArray* L2Hits,const double& maxTheta,const double& tol, bool debug);

int reconstruction(TString fileName = "simulation.root", uint32_t nevents = 0, bool showSimData = true){
  TFile* sourceFile = new TFile(fileName);
	if (!sourceFile) return 1;
	//histogram declarations
	TCanvas* cResAll = new TCanvas("cResAll","Risoluzione",10,10,1200,800);
	TH1I* hResAll = new TH1I("hResAll","Risoluzione inclusva", 200,-20,20);

	TCanvas* cEffVsMult = new TCanvas("cEffVsMult","Efficiency vs Multiplicity",10,10,1200,800);
	TGraphAsymmErrors* hEffVsMult;
	TH1I* hHitVsMult = new TH1I("hHitVsMult", "Totali;Molteplicita';# eventi", 20,0,50);
	TH1I* hHitReconVsMult = new TH1I("hHitReconVsMult", "Ricostruiti;Molteplicita';# eventi", 20,0,50);

	TCanvas* cEffVsMult1S = new TCanvas("cEffVsMult1S","Efficienza (entro 1 sigma)",10,10,1200,800);
	TGraphAsymmErrors* hEffVsMult1S;
	TH1I* hHitVsMult1S = new TH1I("hHitVsMult1S", "Totali;Molteplicita';# eventi", 20,0,50);
	TH1I* hHitReconVsMult1S = new TH1I("hHitReconVsMult1S", "Ricostruiti;Molteplicita';# eventi", 20,0,50);

	TCanvas* cEffVsZ = new TCanvas("cEffVsZ","cEffVsZ",10,10,1200,800);
	TGraphAsymmErrors* hEffVsZ;
	TH1I* hHitVsZ = new TH1I("hHitVsZ", "Totali;Z;# eventi",13,-195,195);
	TH1I* hHitReconVsZ = new TH1I("hHitReconVsZ", "Ricostruiti;Z;# eventi",13,-195,195);

	TTree* tree = (TTree*)(sourceFile->Get("hits"));
	if(!tree) return 2;
  TClonesArray * hitsL1 = new TClonesArray("hit",50);
  TClonesArray * hitsL2 = new TClonesArray("hit",50);
  double ptc[4];
  tree->GetBranch("Vertex")->SetAddress(&ptc[0]);
	tree->GetBranch("L1hit")->SetAddress(&hitsL1);
	tree->GetBranch("L2hit")->SetAddress(&hitsL2);
	tree->SetBranchStatus("*",1); //Activates reading of all branches
	//tree->Print();
  TStopwatch watch;
  uint32_t nEntries = tree->GetEntries();
	uint32_t toProcess = (nevents == 0 ? nEntries : min(nEntries, nevents));
	printf("Starting reconstruction on %u events.\n",toProcess);
	for(uint32_t e = 0; e < toProcess; e++){
		tree->GetEvent(e);
		//printf(" SIZE L1: %d, L2: %d\n",hitsL1->GetEntriesFast(),hitsL2->GetEntriesFast());
		//hitsL1->Print();
    hHitVsMult->Fill(ptc[3]); //hit vs molteplicita
		if(ptc[2] < 53 && ptc[2] > -53) hHitVsMult1S->Fill(ptc[3]); // hit sotto 1 sigma di Z
		hHitVsZ->Fill(ptc[2]); //hit al variare di z

		double zRecon = findZ(hitsL1, hitsL2, minDTheta, minTol, false);
		if(zRecon <-9999.){
			notRecon++;
		}
		else{
			double dZ = zRecon - ptc[2];
			hResAll->Fill(dZ);
			hHitReconVsMult->Fill(ptc[3]);
			if(ptc[2] < 53 && ptc[2] > -53) hHitReconVsMult1S->Fill(ptc[3]);
			hHitReconVsZ->Fill(ptc[2]);
		}
		hitsL1->Clear("C");
		hitsL2->Clear("C");
		printf("\rReconstruction Progress: %u%%",(uint8_t)(e*100/toProcess+1));
  }
	printf("\nNot reconstructed: %d\n",notRecon);
	printf("Elapsed time: %.2f s",watch.CpuTime());
	cResAll->cd();
	hResAll->Draw();
	cEffVsMult->cd();
	hEffVsMult = new TGraphAsymmErrors(hHitReconVsMult, hHitVsMult);
	hEffVsMult->SetTitle("Efficienza;Molteplicita';Ricostruiti/Totali");
	hEffVsMult->SetMinimum(0.5);
	hEffVsMult->Draw();
	cEffVsMult1S->cd();
	hEffVsMult1S = new TGraphAsymmErrors(hHitReconVsMult1S, hHitVsMult1S);
	hEffVsMult1S->SetTitle("Efficienza;Molteplicita' 1 sigma;Ricostruiti/Totali");
	hEffVsMult1S->SetMinimum(0.5);
	hEffVsMult1S->Draw();
	cEffVsZ->cd();
	hEffVsZ = new TGraphAsymmErrors(hHitReconVsZ, hHitVsZ);
	hEffVsZ->SetTitle("Efficienza;Z (mm);Ricostruiti/Totali");
	hEffVsZ->SetMinimum(0.5);
	hEffVsZ->Draw();
	//sourceFile->Close();
	return 0;
}

double findZ(TClonesArray* L1Hits, TClonesArray* L2Hits,const double& maxTheta,const double& tol, bool debug){
	TIter l1iterator(L1Hits);
	TIter l2iterator(L2Hits);
	hit * hitL1;
	hit * hitL2;
	std::vector<double> zRecon;
	TH1I hFindTrend = TH1I("hFindTrend", "findTrend",1200,-200,200);
	while((hitL1=(hit*)l1iterator())){
		while((hitL2=(hit*)l2iterator())){
			if(hitL1->deltaTheta(hitL1,hitL2) < maxTheta){
				tracklet tracker(hitL1, hitL2);
				double zvtx = tracker.findZVertex();
				hFindTrend.Fill(zvtx);
				zRecon.push_back(zvtx);
			}
		}
	}
	if(zRecon.size() == 0) return -99999.;
	double zTrend = findMax(hFindTrend, tol);
	return nearAvg(zRecon, zTrend, tol);
}

double nearAvg(const std::vector<double>& zRecon, const double& zTrend, const double& tol){
	if(zRecon.size() == 1) return zRecon.at(0);
	double temp = 0.;
	double zAvg = 0;
	for(uint32_t k = 0; k < zRecon.size(); k++){
		if(fabs(zRecon.at(k) - zTrend) < tol){
			temp += zRecon.at(k);
			zAvg++;
		}
	}
	return temp/zAvg;
}

double findMax(TH1I& data, const double& tol){
	int binMax = data.GetMaximumBin();
	double zTrend = data.GetXaxis()->GetBinCenter(binMax);

	data.GetXaxis()->SetRange(binMax+1,1200);
	int binMax2 = data.GetMaximumBin();
	double zTrend2 = data.GetXaxis()->GetBinCenter(binMax2);

	if(data.GetBinContent(binMax2) == data.GetBinContent(binMax)){
		if(fabs(zTrend - zTrend2) < tol)
			return ((zTrend + zTrend2)/2.);
		else
			return (fabs(zTrend) < fabs(zTrend2) ? zTrend : zTrend2 );
	}
	else return zTrend;
}
