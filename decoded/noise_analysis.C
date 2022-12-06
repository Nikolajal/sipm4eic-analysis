#include "analysis_utils.h"
#include "style.h"

void
noise_analysis(std::vector<std::string> filenames)
{

  style();
  
  /** histograms **/
  TH1F *hCoarse[32];
  for (int i = 0; i < 32; ++i) {
    std::string name = "hCoarse_" + std::to_string(i);
    hCoarse[i] = new TH1F(name.c_str(), "", 32768, 0., 32768);
  }

  /** loop over input file list **/
  for (const auto filename : filenames) {
    
    /** open file **/
    std::cout << " --- opening decoded file: " << filename << std::endl; 
    auto fin = TFile::Open(filename.c_str());
    if (!fin || !fin->IsOpen()) continue;
    
    /** retrieve tree and link it **/
    data_t data;
    auto tin = (TTree *)fin->Get("alcor");
    auto nev = tin->GetEntries();
    std::cout << " --- found " << nev << " entries in tree " << std::endl;
    tin->SetBranchAddress("fifo", &data.fifo);
    tin->SetBranchAddress("type", &data.type);
    tin->SetBranchAddress("counter", &data.counter);
    tin->SetBranchAddress("column", &data.column);
    tin->SetBranchAddress("pixel", &data.pixel);
    tin->SetBranchAddress("tdc", &data.tdc);
    tin->SetBranchAddress("rollover", &data.rollover);
    tin->SetBranchAddress("coarse", &data.coarse);
    tin->SetBranchAddress("fine", &data.fine);

    /** loop over events in tree **/
    int spill_id = 0;
    for (int iev = 0; iev < nev; ++iev) {
      tin->GetEntry(iev);
      /** only ALCOR hits **/
      if (data.type == kEndSpill) {
        spill_id++;
        continue;
      }
      if (data.type != kAlcorHit) continue;
      
      auto doch = get_dochannel(data.pixel, data.column);
      hCoarse[doch]->Fill(data.coarse);
      
    }
    fin->Close();
  }
  
  auto cNoiseAnalysis = new TCanvas("cNoiseAnalysis", "cNoiseAnalysis", 1600, 800);
  cNoiseAnalysis->Divide(8, 4, 0., 0.);
  for (int i = 0; i < 32; ++i) {
    int cx = i / 4;
    int cy = 3 - i % 4;
    auto hframe = cNoiseAnalysis->cd(cx + 8 * cy + 1)->DrawFrame(0., 0.5, 2000., 50000.);
    if (i == 0) hframe->SetTitle(";coarse (clock cycles);counts");
    hCoarse[i]->Rebin(8);
    hCoarse[i]->Draw("same");
    cNoiseAnalysis->cd(cx + 8 * cy + 1)->SetLogy();    
  }
  cNoiseAnalysis->SaveAs("noise_analysis.png");
  
}
