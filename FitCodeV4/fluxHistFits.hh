#include <string>
#include <cstdlib>

#include "TROOT.h"

class TCanvas;
class TFile;
class TH1;
class TLegend;

class HistObject {

public:

    // Constructor
    HistObject(TH1* theHist);

    // Function value for given bin x
    Double_t getValue(Double_t x, Double_t param);

    TH1* theHist_;
    
private:

//    TH1* theHist_;
    Double_t xMin_;
    Double_t xMax_;
    Double_t dx_;
    Int_t N_;

};

class fluxHistFits {

public:

    // Constructor
    fluxHistFits(std::string FDFileName, std::string FDHistName, std::vector<std::string> listOfNDFileNames, std::string NDHistName);

    // Destructor
    ~fluxHistFits();

    // Run function called by "main"
    void run();

protected:

    // Initialisation
    void init();

    void fitHisto(std::vector<HistObject>& histVector);

private:

    // Input fitted file name
    std::string inFileName_;  

    // Fitted file hist name
    std::string FDHistName_;

    // Input fitted file
    TFile* inFile_;

    // Input fitting file names
    std::vector<std::string> fitFileNames_;

    // Fitting files hist name(s)
    std::string NDHistName_;

    // Drawing canvas
    TCanvas* theCanvas_;

};

class TheFitFunction {

public:

    // Constructor
    TheFitFunction(std::vector<HistObject>& histVect) {
	histos_ = histVect;
    }

    // the operator
    Double_t operator() (Double_t* xArr, Double_t* par);

    // Vector of flux histograms
    std::vector<HistObject> histos_;

};
