#include "fluxHistFits.hh"

#include "TCanvas.h"
#include "TFile.h"
#include "TF1.h"
#include "TH1.h"
#include "TAxis.h"
#include "TStyle.h"
#include "TLegend.h"

#include <cstdlib>
#include <iostream>
#include <fstream>

fluxHistFits::fluxHistFits(std::string FDFileName, std::string FDHistName, std::vector<std::string> listOfNDFileNames, std::string NDHistName) :
    inFileName_(FDFileName),
    FDHistName_(FDHistName),
    inFile_(0),
    fitFileNames_(listOfNDFileNames),
    NDHistName_(NDHistName),
    theCanvas_(0)
{
    this->init();
}

fluxHistFits::~fluxHistFits() 
{

    // Delete drawing canvas
    delete theCanvas_;

    // Close file
    inFile_->Close();

}

void fluxHistFits::init()
{
    // Initialisation
    std::cout<<"Opening file to fit, "<<inFileName_<<std::endl;
    inFile_ = new TFile(inFileName_.c_str(), "read");

    // Define drawing canvas
    theCanvas_ = new TCanvas("theCanvas", "", 900, 700);
    // Plain drawing style
    gROOT->SetStyle("Plain");
    // Remove title, entries, mean, std dev from legend
    gStyle->SetOptStat(0);
    // Fit display options
//    gStyle->SetOptFit(0111);
    gStyle->SetOptFit(0000);
    // Force loaded histograms to follow style
    gROOT->ForceStyle();

}


void fluxHistFits::run() 
{
    //Define vector of histogram objects
    std::vector<HistObject> histos;
    int nFiles = fitFileNames_.size();
    std::cout<<"A total of "<<nFiles<<" files are being accessed for fit parameters"<<std::endl;

    for (int i = 0; i < nFiles; i++) {

	TFile* fitFile = new TFile(fitFileNames_[i].c_str(), "read");
	
	if (fitFile) {

	    // Get the histogram
	    TH1F* theHist = dynamic_cast<TH1F*>(fitFile->Get(NDHistName_.c_str()));
	    if (theHist) {

		HistObject hObj(theHist);
		histos.push_back(hObj);
	    } else {
		std::cout<<"Error, can't find histogram "<<NDHistName_<<" for file "<<fitFileNames_[i]<<std::endl;
	    }

	} else {
	    std::cout<<"Error, file "<<fitFileNames_[i]<<" is null"<<std::endl;
	}
	
    } 
    
    std::cout<<"HistObject vector size = "<<histos.size()<<std::endl;

    this->fitHisto(histos);

}

void fluxHistFits::fitHisto(std::vector<HistObject>& histVector)
{

    // Get the histogram we want to fit
    TH1D* theHist = dynamic_cast<TH1D*>(inFile_->Get(FDHistName_.c_str()));
    if (!theHist) {
	std::cout<<"Error, can't find histogram "<<FDHistName_<<" for file "<<inFileName_<<std::endl;
	return;
    }

    int nHist = histVector.size();
    std::cout<<"nHist = "<<nHist<<std::endl;
    
    TheFitFunction fObj(histVector);

    Double_t xMin(0.75), xMax(7.75);
    Int_t nParams = histVector.size();

    TF1* theFun = new TF1("theFun", fObj, xMin, xMax, nParams);
    // Initialise parameters
    for (int i = 0; i < nParams; i++) {
	theFun->SetParameter(i, 1.0);
	theFun->SetParName(i, (fitFileNames_[i].substr(15,6)+" ").c_str());
    }

    theHist->Fit("theFun","R");
    std::string CanvasName = "xRange=" + std::to_string(xMin).substr(0,4) + "-" + std::to_string(xMax).substr(0,3) + "_" + std::to_string(nHist) + "ParamFit";
    theCanvas_->Print((CanvasName+".png").c_str());


    histVector[0].theHist_->Scale(theFun->GetParameter(0)*pow(10,-20));
    histVector[0].theHist_->SetAxisRange(-60.0E-12, 60.0E-12, "Y");

    TH1* lincomb = dynamic_cast<TH1*>(histVector[0].theHist_->Clone());
    auto legend = new TLegend(0.7, 0.6, 0.985, 0.995);
    legend->AddEntry(histVector[0].theHist_, (theFun->GetParName(0) + std::to_string(theFun->GetParameter(0))).c_str() );

    for (Int_t i = 1; i < nParams; i++) {
	histVector[i].theHist_->Scale(theFun->GetParameter(i)*pow(10,-20));
	lincomb->Add(histVector[i].theHist_);

	legend->SetHeader("Horn Fluxes", "C");
	std::string legendname = theFun->GetParName(i) + std::to_string(theFun->GetParameter(i));
	legend->AddEntry(histVector[i].theHist_, legendname.c_str() );
    }
    theHist->Draw("PLC");
    lincomb->SetLineColor(kViolet-1);
    lincomb->Draw("L SAME PMC"); 
    theCanvas_->Print((CanvasName+"lincomb.png").c_str());

    histVector[0].theHist_->Draw("PLC PMC");
    for (Int_t i = 1; i < nParams; i++) {
	histVector[i].theHist_->Draw("SAME PLC PMC");
    }
/**/    theHist->Draw("L SAME");
    lincomb->Draw("L SAME PMC"); 
    legend->AddEntry(lincomb, "FD Osc Flux Fit" );
    legend->Draw();
//    theCanvas_->Update(); 
    theCanvas_->Print((CanvasName+"weights.png").c_str());

    lincomb->Add(theHist,-1);
    lincomb->Draw("PLC PMC"); 
    theCanvas_->Print((CanvasName+"residuals.png").c_str());
    lincomb->Add(theHist);
    lincomb->Divide(theHist);
    theHist->Divide(theHist);
    lincomb->Draw("PLC PMC"); 
    theHist->Draw("SAME PLC");
    theCanvas_->Print((CanvasName+"normed.png").c_str());

}

Double_t TheFitFunction::operator()(Double_t* xArr, Double_t* par) 
{

    double result(0.0);
    
    double x = xArr[0];
//    std::cout<<"x from operator is.. "<<x<<std::endl;

    int nHist = histos_.size();
    for (int i = 0; i < nHist; i++) {

	HistObject hObj = histos_[i];
	result += hObj.getValue(x, par[i]);

    }

    return result*pow(10,-8);//scaling for readability/fitting
}

HistObject::HistObject(TH1* theHist) :
    theHist_(theHist),
    xMin_(0.0),
    xMax_(0.0),
    dx_(0.0),
    N_(0)
{

    if (theHist) {
	TAxis* xAxis = theHist->GetXaxis();
	xMin_ = xAxis->GetXmin();
	xMax_ = xAxis->GetXmax();
	N_ = xAxis->GetNbins();
	dx_ = (xMax_ - xMin_)/(N_*1.0);
    }

}

Double_t HistObject::getValue(Double_t x, Double_t param) {

    Double_t val(0.0), val1(0.0), val2(0.0), Result(0.0);
//    Double_t xMinFit(-1.0), xMaxFit(10.0);
    
    if (theHist_) {
//     if (xMinFit <= x && x <= xMaxFit) {
	// Find the bin centre for the given x
	Double_t binC = Double_t(((x - xMin_)/dx_) + 0.5);
	if ( binC == floor(binC) ) {
	    val = theHist_->GetBinContent(binC);
	    Result = param*val;
	}
	else {
	    val1 = theHist_->GetBinContent(floor(binC));
	    val2 = theHist_->GetBinContent(ceil(binC));
//	std::cout<<"param*val1 = "<<param*val1<<std::endl;
//	std::cout<<"param*val2 ="<<param*val2<<std::endl;
	Result = ((param*val1)*(ceil(binC)-binC) + (param*val2)*(binC-floor(binC)));
//	std::cout<<Result<<std::endl;
	}
//    }
   }

    return Result*pow(10,-12);//correction factor for difference in histogram Y axis scale

}

// Main program
int main(const int argc, const char** argv) 
{
    // Define input files:
    // ./fitProj fileName
    std::vector<std::string> NDFileN;
    std::string FDFileN;

    // First filename is location of oscillated FD histogram
    FDFileN = "Histo5/fluxhistN2934FD.root";
//    FDFileN = "Histo2/fluxhist0N2934FD.root";
//    FDFileN = "exampleHistos/fluxhist1504.root";
    std::cout<<"Fitting to "<<FDFileN<<std::endl;

    // Following files give location of fitting histograms
    NDFileN.push_back("Histo5/fluxhist0N3001.root");
    NDFileN.push_back("Histo5/fluxhist0N3002.root");
    //NDFileN.push_back("Histo5/fluxhist0N3003.root");
    NDFileN.push_back("Histo5/fluxhist0N3004.root");
//    NDFileN.push_back("Histo5/fluxhist0N3005.root");
//    NDFileN.push_back("Histo5/fluxhist0N1501.root");
    NDFileN.push_back("Histo5/fluxhist0N0752.root");
//    NDFileN.push_back("Histo5/fluxhist0N1002.root");//
    NDFileN.push_back("Histo5/fluxhist0N1252.root");
    NDFileN.push_back("Histo5/fluxhist0N2002.root");
    //NDFileN.push_back("Histo5/fluxhist0N2502.root");//
    //NDFileN.push_back("Histo5/fluxhist0N1503.root");//
//    NDFileN.push_back("Histo5/fluxhist0N2003.root");//
    NDFileN.push_back("Histo5/fluxhist0N1254.root");
    NDFileN.push_back("Histo5/fluxhist0N1505.root");
    //NDFileN.push_back("Histo5/fluxhist0N2505.root");
//    NDFileN.push_back("Histo5/fluxhist0N3505.root");
    //NDFileN.push_back("Histo5/fluxhist0A0752.root");
    NDFileN.push_back("Histo5/fluxhist0A3002.root");
//    NDFileN.push_back("Histo5/fluxhist1N3004.root");
//    NDFileN.push_back("Histo5/fluxhist2N3004.root");
//    NDFileN.push_back("Histo5/fluxhist3N3004.root");
//    NDFileN.push_back("Histo5/fluxhist4N3004.root");
//    NDFileN.push_back("Histo5/fluxhist5N3004.root");
//    NDFileN.push_back("Histo5/fluxhist6N3004.root");
    
    int nNDFiles = NDFileN.size();
    for (int i=0; i<nNDFiles; i++) {
	std::cout<<"Fitting with "<<NDFileN[i]<<std::endl;
    }

    std::string NDHistN = "numufluxatDUNEND;1";
    std::string FDHistN = "numu_fluxosc_forplots;1";

    std::cout<<"Using ND Hists "<<NDHistN<<std::endl;
    std::cout<<"Using FD Hist "<<FDHistN<<std::endl;

    // Declare object
    fluxHistFits a(FDFileN, FDHistN, NDFileN, NDHistN);
    // Run code
    a.run();

    // Return zero for OK
    return 0;

}
