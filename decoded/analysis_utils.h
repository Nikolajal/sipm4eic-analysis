#pragma once
//!
string
get_environment_variable
( string search_key ) {
    extern char **environ;
    int iTer = 0;
    string result = "./";
    while(environ[iTer]) {
        string current_string = environ[iTer++];
        auto key_position   = current_string.find(search_key);
        auto key_length     = search_key.size();
        if ( key_position < current_string.size() && key_position == 0 && current_string.substr(key_length,1) == "=" ) {
            result = current_string.substr(key_length+key_position+1);
        }
    }
    return result;
}
//!
namespace analysis_utils {
//!  Data Structures
//!  --- ALCOR output structure
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
//!
//! --- Datatype Structure                     -------------------------------------------------------------------------------------------------------
enum type_t {
  kAlcorHit = 1,
  kTriggerTag = 9,
  kStartSpill = 7,
  kEndSpill = 15
};
//!
//! Analysis types                                -------------------------------------------------------------------------------------------------------
typedef std::vector<data_t>             channel_data_t;
typedef std::map<int, channel_data_t>   chip_data_t;
typedef std::map<int, chip_data_t>      frame_data_t;
typedef std::map<int, frame_data_t>     spill_data_t;
typedef std::map<int, spill_data_t>     framed_data_t;
typedef std::tuple<TString,int>         global_identifier_t;
//!
//! Variables                                        ---------------------------------------------------------------------------------------------------------
//! --- General
const int       frame_size      = 1024; // [clock cycles]
const Int_t     kGlobalIndexRange = 768;
const Int_t     kGlobalIndexTiming_Start = 512;
const Int_t     kGlobalIndexTiming_End = kGlobalIndexRange;
std::map<std::string, int> _next_spill;
//! --- Finetune analysis
const Int_t     kFineRange  = 512;
const TString   intput_rawdata_directory                        = TString(get_environment_variable("ALCOR_DATA_DIR"))+TString("/");
const TString   intput_rawdata_decoded_file_dir                 = intput_rawdata_directory + TString("./%s/decoded/");
const TString   intput_rawdata_decoded_file                     = intput_rawdata_decoded_file_dir + TString("./alcdaq.fifo_%i.root");
const TString   output_preprocess_directory                     = TString(get_environment_variable("ALCOR_WORK_DIR"))+TString("/")+TString("PreProcessedData/");
//! --- Convertion
//! ---    --- Coarse
const double coarse_to_s  = 3.1250000e-09;
const double coarse_to_us = 3.1250000e-03;
const double coarse_to_ns = 3.1250000;
//! --- --- Rollover
const int rollover_to_coarse = 32768.;
const double rollover_to_s  = 0.0001024;
const double rollover_to_us = 102.4;
const double rollover_to_ns = 102400.;
//! --- --- Electronics-Detector conversions
int eo2do[32] = {22, 20, 18, 16, 24, 26, 28, 30, 25, 27, 29, 31, 23, 21, 19, 17, 9, 11, 13, 15, 7, 5, 3, 1, 6, 4, 2, 0, 8, 10, 12, 14};
//! --- Filenames
std::vector<std::string> data_filenames = {
  "alcdaq.fifo_0.root", "alcdaq.fifo_1.root", "alcdaq.fifo_2.root", "alcdaq.fifo_3.root",
  "alcdaq.fifo_4.root", "alcdaq.fifo_5.root", "alcdaq.fifo_6.root", "alcdaq.fifo_7.root",
  "alcdaq.fifo_8.root", "alcdaq.fifo_9.root", "alcdaq.fifo_10.root", "alcdaq.fifo_11.root",
  "alcdaq.fifo_12.root", "alcdaq.fifo_13.root", "alcdaq.fifo_14.root", "alcdaq.fifo_15.root"
};
std::vector<std::string> timing_filenames = {
  "alcdaq.fifo_16.root", "alcdaq.fifo_17.root", "alcdaq.fifo_18.root", "alcdaq.fifo_19.root",
  "alcdaq.fifo_20.root", "alcdaq.fifo_21.root", "alcdaq.fifo_22.root", "alcdaq.fifo_23.root"
};
std::vector<std::string> extra_filenames = {
    "alcdaq.fifo_24.root"
};
std::vector<std::string> all_filenames = {"alcdaq.fifo_0.root", "alcdaq.fifo_1.root", "alcdaq.fifo_2.root", "alcdaq.fifo_3.root",
    "alcdaq.fifo_4.root", "alcdaq.fifo_5.root", "alcdaq.fifo_6.root", "alcdaq.fifo_7.root",
    "alcdaq.fifo_8.root", "alcdaq.fifo_9.root", "alcdaq.fifo_10.root", "alcdaq.fifo_11.root",
    "alcdaq.fifo_12.root", "alcdaq.fifo_13.root", "alcdaq.fifo_14.root", "alcdaq.fifo_15.root",
    "alcdaq.fifo_16.root", "alcdaq.fifo_17.root", "alcdaq.fifo_18.root", "alcdaq.fifo_19.root",
    "alcdaq.fifo_20.root", "alcdaq.fifo_21.root", "alcdaq.fifo_22.root", "alcdaq.fifo_23.root",
    "alcdaq.fifo_24.root"};
//! --- Verbose levels
bool kVerboseInfo       = false;
bool kVerboseWarning    = true;
bool kVerboseError      = true;
bool kVerboseFatal      = true;
//!
//! Functions                                       ---------------------------------------------------------------------------------------------------------
//! --- Declaration
//! --- --- Getter & Setters
std::pair<int, int>     get_index                       ( int pixel, int column );
int                     get_eochannel                   ( int pixel, int column );
int                     get_global_index                ( int fifo, int pixel, int column, int tdc, bool kUseFIFO = true );
std::tuple<Int_t,Int_t,Int_t,Int_t>
                        get_full_info                   ( Int_t iIndex );
void                    set_verbose_level               ( int kVerboseLevel );
//! --- --- Data handling
bool                    sort_data                       ( data_t i, data_t j );
void                    load_tree                       ( TTree* kInputTree, data_t &kReadData );
//bool                    populate_framed_data            ( framed_data_t &framed_data, std::string dirname, std::vector<std::string> filenames = all_filenames, int frame_size = 1024 );
//! --- Implementation
//! --- --- Getter & Setters
//!
int
get_dochannel
 ( int pixel, int column ) {
  int eoch = pixel + 4 * column;
  int doch = eo2do[eoch];
  return doch;
}
//!
std::pair<int, int>
get_index
 ( int pixel, int column ) {
  int doch = get_dochannel(pixel, column);
  int ix = doch / 4;
  int iy = doch % 4;
  return {ix, iy};
}
//!
int
get_global_index
( int fifo, int pixel, int column, int tdc, bool kUseFIFO = true ) {
    //! 4 TDCs for each pixel
    //! 4 pixel for each column
    //! 8 columns for each chip
    Bool_t kSkipCorrupted = false;
    if ( fifo   < 0 ) kSkipCorrupted = true;
    if ( pixel  < 0 ) kSkipCorrupted = true;
    if ( column < 0 ) kSkipCorrupted = true;
    if ( tdc    < 0 ) kSkipCorrupted = true;
    if ( kSkipCorrupted ) {
        std::cout << "[ERROR] Invalid negative value given, failed to determine global index!" << std::endl;
        return -1;
    }
    int chip = kUseFIFO ? fifo / 4 : fifo;
    int index = tdc + 4 * pixel + 16 * column + 128 * chip;
    return index;
}
//!
std::tuple<Int_t,Int_t,Int_t,Int_t>
get_full_info
( Int_t iIndex ) {
    //! Result { Chip, Pixel, Column, TDC }
    std::tuple<Int_t,Int_t,Int_t,Int_t> kResult;
    get<0>(kResult) = (iIndex/128)%6;
    get<1>(kResult) = (iIndex/4)%4;
    get<2>(kResult) = (iIndex/16)%8;
    get<3>(kResult) = iIndex%4;
    return kResult;
}
//!
bool
sort_data
 ( data_t i, data_t j ) {
  auto itime = i.coarse + i.rollover * rollover_to_coarse;
  auto jtime = j.coarse + j.rollover * rollover_to_coarse;
  return (itime < jtime);
}
//!
void
load_tree
( TTree* kInputTree, data_t &kReadData ) {
    if ( !kInputTree ) { std::cout << "[ERROR] analysis_utils::load_tree for data_t: Invalid Tree given" << std::endl; return; }
    kInputTree->SetBranchAddress(   "fifo",     &kReadData.fifo     );
    kInputTree->SetBranchAddress(   "type",     &kReadData.type     );
    kInputTree->SetBranchAddress(   "counter",  &kReadData.counter  );
    kInputTree->SetBranchAddress(   "column",   &kReadData.column   );
    kInputTree->SetBranchAddress(   "pixel",    &kReadData.pixel    );
    kInputTree->SetBranchAddress(   "tdc",      &kReadData.tdc      );
    kInputTree->SetBranchAddress(   "rollover", &kReadData.rollover );
    kInputTree->SetBranchAddress(   "coarse",   &kReadData.coarse   );
    kInputTree->SetBranchAddress(   "fine",     &kReadData.fine     );
}
//!
bool
populate_framed_data(framed_data_t &framed_data, std::string dirname, std::vector<std::string> filenames = all_filenames, int frame_size = 1024)
{
  /** start clean **/
  bool has_data = false;
  framed_data.clear();
  
  /** loop over input file list **/
  for (const auto filename : filenames) {

    /** reset spill pointer **/
    if (!_next_spill.count(filename)) _next_spill[filename] = 0;
    
    const auto pathname = dirname + "/" + filename;
    
    /** open file **/
    std::cout << " --- opening decoded file: " << pathname << std::endl;
    auto fin = TFile::Open(pathname.c_str());
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
    int spill = 0;
    for (int iev = _next_spill[filename]; iev < nev; ++iev) {
      tin->GetEntry(iev);

      /** start of spill **/
      if (data.type == kStartSpill) {
        has_data = true;
        continue;
      }

      /** ALCOR hit **/
      if (data.type == kAlcorHit) {
        auto chip = data.fifo / 4;
        auto doch = get_dochannel(data.pixel, data.column);
        auto coarse = data.coarse + data.rollover * rollover_to_coarse;
        auto frame = coarse / frame_size;
        framed_data[spill][frame][chip][doch].push_back(data);
      }
      
      /** end of spill **/
      if (data.type == kEndSpill) {
        /** sort channel data **/
        auto &frames = framed_data[spill];
        for (auto &frame_data : frames) {
          auto &chips = frame_data.second;
          for (auto &chip_data : chips) {
            auto &channels = chip_data.second;
            for (auto &channel_data : channels) {
              auto &hits = channel_data.second;
              std::sort(hits.begin(), hits.end(), sort_data);
            }
          }
        }
        spill++;
        _next_spill[filename] = iev + 1;
        break;
      }
      
    } /** end of loop over events **/
    
    fin->Close();
    
  } /** end of loop over input files **/

    //! Reset spill counter if end of data
    if ( !has_data ) _next_spill.clear();
  return has_data;
}
//!
void
set_verbose_level
( int kVerboseLevel = 4 ) {
    switch ( kVerboseLevel ) {
        case 0:
            kVerboseInfo       = false;
            kVerboseWarning    = false;
            kVerboseError      = false;
            kVerboseFatal      = false;
            break;
            
        case 1:
            kVerboseInfo       = false;
            kVerboseWarning    = false;
            kVerboseError      = false;
            kVerboseFatal      = true;
            break;
            
        case 2:
            kVerboseInfo       = false;
            kVerboseWarning    = false;
            kVerboseError      = true;
            kVerboseFatal      = true;
            break;
            
        case 3:
            kVerboseInfo       = false;
            kVerboseWarning    = true;
            kVerboseError      = true;
            kVerboseFatal      = true;
            break;
            
        case 4:
            kVerboseInfo       = true;
            kVerboseWarning    = true;
            kVerboseError      = true;
            kVerboseFatal      = true;
            break;
        
        default:
            kVerboseInfo       = false;
            kVerboseWarning    = true;
            kVerboseError      = true;
            kVerboseFatal      = true;
            break;
    }
}
//!
} //! namespace analysis_utils
