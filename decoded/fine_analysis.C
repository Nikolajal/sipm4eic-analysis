/**
   example of usage:

   root [0] .L fine_analysis.C 
   root [1] TH2 *hFine_all = fine_fill({"alcdaq.fifo_0.root", "alcdaq.fifo_1.root"});
   --- opening decoded file: alcdaq.fifo_0.root
   --- found 2151601 entries in tree 
   --- opening decoded file: alcdaq.fifo_1.root
   --- found 1875405 entries in tree 
   root [2] hFine_all->Draw("colz");
   root [3] TH1 *hFine_0 = fine_histo(hFine_all, 0);
   root [4] hFine_0->Draw();

 **/

int get_index(int fifo, int pixel, int column, int tdc);
TH2 *fine_fill(std::vector<std::string> filenames);
TH1 *fine_histo(TH2 *hFine_all, int index);

struct data_t {
  int fifo;
  int type;
  int counter;
  int column;
  int pixel;
  int tdc;
  int rollover;
  int coarse;
  int fine;
};

enum type_t {
  kAlcorHit = 1,
  kTriggerTag = 9,
  kStartSpill = 7,
  kEndSpill = 15
};

/*******************************************************************************/

int
get_index(int fifo, int pixel, int column, int tdc)
{
  /**
     4 TDCs for each pixel
     4 pixel for each column
     8 columns for each chip
  **/
  int chip = fifo / 4;
  int index = tdc + 4 * pixel + 16 * column + 128 * chip;
  return index;
}

/*******************************************************************************/

TH2 *
fine_fill(std::vector<std::string> filenames)
{

  /** histogram with fine distribution **/
  auto hFine_all = new TH2F("hFine_all", ";index;fine", 768, 0, 768, 512, 0, 512);

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
    for (int iev = 0; iev < nev; ++iev) {
      tin->GetEntry(iev);
      /** only ALCOR hits **/
      if (data.type != kAlcorHit) continue;
      int index = get_index(data.fifo, data.pixel, data.column, data.tdc);
      hFine_all->Fill(index, data.fine);
    }
    fin->Close();
  }
    
  return hFine_all;
}

/*******************************************************************************/

TH1 *
fine_histo(TH2 *hFine_all, int index)
{
  std::string name = "hFine_" + std::to_string(index);
  TH1 *hFine = hFine_all->ProjectionY(name.c_str(), index + 1, index + 1);
  return hFine;
}

/*******************************************************************************/
