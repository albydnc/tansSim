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
#include "TGraphErrors.h"
#include <Riostream.h>
#include <vector>

int notRecon = 0;
// methods declarations
double findMax(TH1I &data, const double &tol);
double nearAvg(const std::vector<double> &zRecon, const double &zTrend, const double &tol);
double findZ(TClonesArray *L1Hits, TClonesArray *L2Hits, const double &maxTheta, const double &tol, bool debug);

int nRec = 0;
//------------------------RECONSTRUCTION METHOD -------------------------
int reconstruction(TString fileName = "simulation.root", uint32_t nevents = 0, const double minTol = 1,
                   const double minDTheta = 0.01)
{
   // loading source file
   TFile *sourceFile = new TFile(fileName);
   if (!sourceFile) return 1;
   // histogram and canvas declarations
   TCanvas *     cResAll = new TCanvas("cResAll", "Risoluzione", 10, 10, 1200, 800);
   TH1I *        hResAll = new TH1I("hResAll", "Risoluzione inclusva;Z [mm];#eventi", 200, -2, 2);
   double        sigmaZ[13], ssigmaZ[13], xZ[13], sxZ[13];
   TCanvas *     cResVsZtrue = new TCanvas("cResVsZtrue", "Risoluzione vs Z", 10, 10, 1200, 800);
   TGraphErrors *hResVsZtrue;
   TH1I *        hdZ[13];
   for (uint8_t i = 0; i < 13; i++) {
      xZ[i]  = -180 + i * 30;
      hdZ[i] = new TH1I(Form("dZ%u", i), Form("dZ%u", i), 100, -5, 5);
      sxZ[i] = 15;
   }
   double        sigmaMult[14], ssigmaMult[14], xMult[14], sxMult[14];
   TCanvas *     cResVsMult = new TCanvas("cResVsMult", "Risoluzione vs Mult", 10, 10, 1200, 800);
   TGraphErrors *hResVsMult;
   TH1I *        hdMult[14];
   for (uint8_t i = 0; i < 14; i++) {
      xMult[i]  = 2 + i * 4;
      hdMult[i] = new TH1I(Form("dMult%u", i), Form("dMult%u", i), 100, -5, 5);
      sxMult[i] = 2;
   }
   TCanvas *          cEffVsMult = new TCanvas("cEffVsMult", "Efficiency vs Multiplicity", 10, 10, 1200, 800);
   TGraphAsymmErrors *hEffVsMult;
   TH1I *             hHitVsMult      = new TH1I("hHitVsMult", "Totali;Molteplicita';# eventi", 20, 0, 50);
   TH1I *             hHitReconVsMult = new TH1I("hHitReconVsMult", "Ricostruiti;Molteplicita';# eventi", 20, 0, 50);

   TCanvas *          cEffVsMult1S = new TCanvas("cEffVsMult1S", "Efficienza (entro 1 sigma)", 10, 10, 1200, 800);
   TGraphAsymmErrors *hEffVsMult1S;
   TH1I *             hHitVsMult1S = new TH1I("hHitVsMult1S", "Totali;Molteplicita';# eventi", 20, 0, 50);
   TH1I *hHitReconVsMult1S         = new TH1I("hHitReconVsMult1S", "Ricostruiti;Molteplicita';# eventi", 20, 0, 50);

   TCanvas *          cEffVsZ = new TCanvas("cEffVsZ", "cEffVsZ", 10, 10, 1200, 800);
   TGraphAsymmErrors *hEffVsZ;
   TH1I *             hHitVsZ      = new TH1I("hHitVsZ", "Totali;Z;# eventi", 13, -195, 195);
   TH1I *             hHitReconVsZ = new TH1I("hHitReconVsZ", "Ricostruiti;Z;# eventi", 13, -195, 195);
   // getting simulation data from TTree
   TTree *tree = (TTree *)(sourceFile->Get("hits"));
   if (!tree) return 2;
   TClonesArray *hitsL1 = new TClonesArray("hit", 100);
   TClonesArray *hitsL2 = new TClonesArray("hit", 100);
   double        ptc[4];
   tree->GetBranch("Vertex")->SetAddress(&ptc[0]);
   tree->GetBranch("L1hit")->SetAddress(&hitsL1);
   tree->GetBranch("L2hit")->SetAddress(&hitsL2);
   tree->SetBranchStatus("*", 1); // Activates reading of all branches
   TStopwatch watch;
   uint32_t   nEntries = tree->GetEntries();
   // get event number to process
   uint32_t toProcess = (nevents == 0 ? nEntries : min(nEntries, nevents));
   printf("Starting reconstruction on %u events.\n", toProcess);
   //------------------------Reconstruction Algorithm
   for (uint32_t e = 0; e < toProcess; e++) {
      //read event
      tree->GetEvent(e);
      hHitVsMult->Fill(ptc[3]);                                    // hit vs molteplicita
      if (ptc[2] < 53 && ptc[2] > -53) hHitVsMult1S->Fill(ptc[3]); // hit sotto 1 sigma di Z
      hHitVsZ->Fill(ptc[2]);                                       // hit al variare di z
      // reconstruct z
      double zRecon = findZ(hitsL1, hitsL2, minDTheta, minTol, false);
      if (zRecon < -90000.) {
         notRecon++;
      } else {
         // calculate deviation
         double dZ = zRecon - ptc[2];
         hResAll->Fill(dZ);
         hHitReconVsMult->Fill(ptc[3]);
         if (ptc[2] < 53 && ptc[2] > -53) hHitReconVsMult1S->Fill(ptc[3]);
         hHitReconVsZ->Fill(ptc[2]);
         for (uint8_t i = 0; i < 13; i++) {
            if (ptc[2] - xZ[i] < 15 && ptc[2] - xZ[i] > -15) {
               hdZ[i]->Fill(dZ);
               break;
            }
         }
         for (uint8_t i = 0; i < 14; i++) {
            if (ptc[3] - xMult[i] <= 2 && ptc[3] - xMult[i] >= -2) {
               hdMult[i]->Fill(dZ);
               break;
            }
         }
      }
      hitsL1->Clear("C");
      hitsL2->Clear("C");
      printf("\rReconstruction Progress: %u%%", (uint8_t)(e * 100 / toProcess + 1));
   }
   printf("\nNot reconstructed: %d\n", notRecon);
   printf("\nNot reconstructed 2: %d\n", nRec);
   printf("Elapsed time: %.2f s", watch.CpuTime());
   // draw histograms
   for (uint8_t i = 0; i < 13; i++) {
      sigmaZ[i]  = hdZ[i]->GetStdDev();
      ssigmaZ[i] = hdZ[i]->GetStdDevError();
   }
   cResVsZtrue->cd();
   hResVsZtrue = new TGraphErrors(13, xZ, sigmaZ, sxZ, ssigmaZ);
   hResVsZtrue->SetTitle("Risoluzione vs Z;Z [mm];Risoluzione [mm]");
   hResVsZtrue->SetMinimum(0);
   hResVsZtrue->Draw();

   for (uint8_t i = 0; i < 14; i++) {
      sigmaMult[i]  = hdMult[i]->GetStdDev();
      ssigmaMult[i] = hdMult[i]->GetStdDevError();
   }
   cResVsMult->cd();
   hResVsMult = new TGraphErrors(13, xMult, sigmaMult, sxMult, ssigmaMult);
   hResVsMult->SetTitle("Risoluzione vs Molteplicita';Molteplicita' ;Risoluzione [mm]");
   hResVsMult->SetMinimum(0);
   hResVsMult->Draw();
   cResAll->cd();
   hResAll->Draw();
   cEffVsMult->cd();
   hEffVsMult = new TGraphAsymmErrors(hHitReconVsMult, hHitVsMult);
   hEffVsMult->SetTitle("Efficienza;Molteplicita';Ricostruiti/Totali");
   hEffVsMult->SetMinimum(0);
   hEffVsMult->Draw();
   cEffVsMult1S->cd();
   hEffVsMult1S = new TGraphAsymmErrors(hHitReconVsMult1S, hHitVsMult1S);
   hEffVsMult1S->SetTitle("Efficienza;Molteplicita' 1 sigma;Ricostruiti/Totali");
   hEffVsMult1S->SetMinimum(0);
   hEffVsMult1S->Draw();
   cEffVsZ->cd();
   hEffVsZ = new TGraphAsymmErrors(hHitReconVsZ, hHitVsZ);
   hEffVsZ->SetTitle("Efficienza;Z (mm);Ricostruiti/Totali");
   hEffVsZ->SetMinimum(0);
   hEffVsZ->Draw();
   return 0;
}

// finding z method
double findZ(TClonesArray *L1Hits, TClonesArray *L2Hits, const double &maxTheta, const double &tol, bool debug)
{
   TIter               l1iterator(L1Hits);
   TIter               l2iterator(L2Hits);
   hit *               hitL1;
   hit *               hitL2;
   std::vector<double> zRecon;
   TH1I                hFindTrend = TH1I("hFindTrend", "findTrend", 200, -200, 200);
   // get hits from array
   l1iterator.Reset();
   l2iterator.Reset();
   while ((hitL1 = (hit *)l1iterator.Next())) {
      while ((hitL2 = (hit *)l2iterator.Next())) {
         // max theta filtering
         if (hitL1->deltaTheta(hitL1, hitL2) < maxTheta) {
            // finding vertex via tracklets
            tracklet tracker(hitL1, hitL2, 40., 70.);
            double   zvtx = tracker.findZVertex();
            hFindTrend.Fill(zvtx);
            zRecon.push_back(zvtx);
         }
      }
      l2iterator.Reset();
   }

   if (zRecon.size() == 0) {
      nRec++;
      return -99999.;
   }
   // finding mode value
   double zTrend = findMax(hFindTrend, tol);
   // averaging in tolerance
   return nearAvg(zRecon, zTrend, tol);
}
// averaging value in a tolerance span method
double nearAvg(const std::vector<double> &zRecon, const double &zTrend, const double &tol)
{
   if (zRecon.size() == 1) return zRecon.at(0);
   double temp = 0.;
   int    zAvg = 0;
   for (uint32_t k = 0; k < zRecon.size(); k++) {
      if (fabs(zRecon.at(k) - zTrend) < tol) {
         temp += zRecon.at(k);
         zAvg++;
      }
   }
   if (zAvg > 0) {
      return temp / zAvg;
   } else {
      return -99999.;
   }
}
// finding mode method
double findMax(TH1I &data, const double &tol)
{
   int    binMax = data.GetMaximumBin();
   double zTrend = data.GetXaxis()->GetBinCenter(binMax);
   data.GetXaxis()->SetRange(binMax + 1, 200);
   int    binMax2 = data.GetMaximumBin();
   double zTrend2 = data.GetXaxis()->GetBinCenter(binMax2);

   if (data.GetBinContent(binMax2) == data.GetBinContent(binMax)) {
      if (fabs(zTrend - zTrend2) < tol) {
         return ((zTrend + zTrend2) / 2.);
      } else {
         return (fabs(zTrend) < fabs(zTrend2) ? zTrend : zTrend2);
      }
   } else {
      return zTrend;
   }
}
