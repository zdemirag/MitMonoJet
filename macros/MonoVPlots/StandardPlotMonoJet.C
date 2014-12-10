#include<vector>

#if !defined (__CINT__) || defined (__MAKECINT__)
#include "THStack.h"
#include "TGaxis.h"
#include "TH1D.h"
#include "TLatex.h"
#include "TPad.h"
#include "TCanvas.h"
#include "TAxis.h"
#include "TLegend.h"
#include "TFrame.h"
#include "TExec.h"
#include <iostream>
#include "TPaveText.h"
#endif

enum samp {iHiggs,iOther,iDiboson,iTop,iWjets,iZjets,iGjets,ifromQ,ifromB,ifromG,ifromW,nSamples};

float xPos[nSamples+1] = {0.21,0.21,0.21,0.21,0.21,0.45,0.45,0.45,0.45,0.45,0.45,0.45}; 
float yOff[nSamples+1] = {0,1,2,3,4,0,1,2,3,4,5,6};

const Float_t _tsize   = 0.035;
const Float_t _xoffset = 0.20;
const Float_t _yoffset = 0.05;

//------------------------------------------------------------------------------
// GetMaximumIncludingErrors
//------------------------------------------------------------------------------
Float_t GetMaximumIncludingErrors(TH1D* h)
{
    Float_t maxWithErrors = 0;

    for (Int_t i=1; i<=h->GetNbinsX(); i++) {

        Float_t binHeight = h->GetBinContent(i) + h->GetBinError(i);

        if (binHeight > maxWithErrors) maxWithErrors = binHeight;
    }

    return maxWithErrors;
}


//------------------------------------------------------------------------------
// AxisFonts
//------------------------------------------------------------------------------
void AxisFonts(TAxis*  axis,
        TString coordinate,
        TString title)
{
    axis->SetLabelFont  (   42);
    axis->SetLabelOffset(0.015);
    axis->SetLabelSize  (0.050);
    axis->SetNdivisions (  505);
    axis->SetTitleFont  (   42);
    axis->SetTitleOffset(  1.5);
    axis->SetTitleSize  (0.050);

    if (coordinate == "y") axis->SetTitleOffset(1.6);

    axis->SetTitle(title);
}

//------------------------------------------------------------------------------
// THStackAxisFonts
//------------------------------------------------------------------------------
void THStackAxisFonts(THStack* h,
        TString  coordinate,
        TString  title)
{
    TAxis* axis = NULL;

    if (coordinate.Contains("x")) axis = h->GetHistogram()->GetXaxis();
    if (coordinate.Contains("y")) axis = h->GetHistogram()->GetYaxis();

    AxisFonts(axis, coordinate, title);
}

//------------------------------------------------------------------------------
// DrawLegend
//------------------------------------------------------------------------------
void DrawLegend(Float_t x1,
        Float_t y1,
        TH1D*   hist,
        TString label,
        TString option)
{
    TLegend* legend = new TLegend(x1,
            y1,
            x1 + _xoffset,
            y1 + _yoffset);

    legend->SetBorderSize(     0);
    legend->SetFillColor (     0);
    legend->SetTextAlign (    12);
    legend->SetTextFont  (    42);
    legend->SetTextSize  (_tsize);

    legend->AddEntry(hist, label.Data(), option.Data());

    legend->Draw();
}


class StandardPlot {

    public: 
        StandardPlot() { _hist.resize(nSamples,0); _isBlinded = false; _data = 0; _breakdown = false; _mass = 0;_isHWWOverlaid = false; _nsel = 0; _isVarBins = false;}
        void setBlinded  (bool b)                   { _isBlinded       = b;} 
        void setMCHist   (const samp &s, TH1D * h)  { _hist[s]       = h;  } 
        void setDataHist  (TH1D * h)                { _data          = h;  } 
        void setHWWOverlaid(bool b)                 { _isHWWOverlaid = b;  }
        void setNsel  (Int_t x)                     { _nsel          = x;  } 
        void setVarBins  (bool b)                   { _isVarBins     = b;  } 
 
        TH1D* getDataHist() { return _data; }

        void setMass(const int &m) {_mass=m;}

        TH1* DrawAndRebinTo(const int &rebinTo) {

            if(rebinTo == 0) return Draw();
            int rebin = 0, nbins = 0;
            for (int i=0; i<nSamples; i++) {

                // in case the user doesn't set it
                if( !_hist[i] ) continue;
                nbins = _hist[i]->GetNbinsX();
            }
            if (nbins == 0) return Draw();

            rebin = nbins / rebinTo;
            while(nbins % rebin != 0) rebin--;
            return Draw(rebin);

        }

        TH1* Draw(const int &rebin=1) {

            Color_t _sampleColor[nSamples];
            _sampleColor[iTop  ]   = kRed+1;
            _sampleColor[iWjets]   = kAzure-2;
            _sampleColor[iZjets]   = kSpring-7;
            _sampleColor[iDiboson] = kGray+1;
            _sampleColor[iGjets  ] = kBlue-1;
            _sampleColor[ifromW]   = kAzure-9;
            _sampleColor[ifromQ]   = kSpring+4;
            _sampleColor[ifromB]   = kYellow-7;
            _sampleColor[ifromG]   = kOrange-2;
            _sampleColor[iOther]   = kOrange+2;
            _sampleColor[iHiggs]   = kOrange-2;

            THStack* hstack = new THStack();
            TH1D* hSum = (TH1D*)_data->Clone();
            hSum->Rebin(rebin);
            hSum->Scale(0.0);
            TAxis *xa;
	    
            for (int i=0; i<nSamples; i++) {

                // in case the user doesn't set it
                if( !_hist[i] ) continue;
     	    	bool modifyXAxis = false;
         		if(modifyXAxis == true){
                    xa =_hist[i]->GetXaxis();
   	                xa->SetLabelSize(0.001);
	                xa->SetLabelColor(4);
  	    	        for(Int_t k=1;k<=_hist[i]->GetNbinsX();++k){
                        xa->SetBinLabel(1 ,"m_{ll}< 250&m_{jj}< 750");
                        xa->SetBinLabel(2 ,"m_{ll}>=250&m_{jj}< 750");
                        xa->SetBinLabel(3 ,"m_{ll}< 250&m_{jj}>=750");
                        xa->SetBinLabel(4 ,"m_{ll}>=250&m_{jj}>=750");
                        xa->SetRangeUser(1,4);
                    }
     		}
            _hist[i]->Rebin(rebin);
            _hist[i]->SetLineColor(_sampleColor[i]);
            _hist[i]->SetFillColor(_sampleColor[i]);
    
            if(i == iHiggs) continue;
    
            _hist[i]->SetFillStyle(1001);
    
            hstack->Add(_hist[i]);
            hSum->Add(_hist[i]);
            }

            if(_data) _data->Rebin(rebin);
            if(_data) _data->SetLineColor  (kBlack);
            if(_data) _data->SetMarkerStyle(kFullCircle);
            hstack->Draw("hist");

            if(_hist[iHiggs] && _hist[iHiggs]->GetSumOfWeights() > 0) {
              _hist[iHiggs]->SetFillStyle(1);
              _hist[iHiggs]->SetLineWidth(3);
              _hist[iHiggs]->Draw("hist,same");
            }

	    bool plotSystErrorBars = true;
	    if(plotSystErrorBars) {
  	      TGraphAsymmErrors * gsyst = new TGraphAsymmErrors(hSum);
              for (int i = 0; i < gsyst->GetN(); ++i) {
                double systBck = 0.06;
                gsyst->SetPointEYlow (i, sqrt(hSum->GetBinError(i+1)*hSum->GetBinError(i+1)+hSum->GetBinContent(i+1)*hSum->GetBinContent(i+1)*systBck*systBck));
                gsyst->SetPointEYhigh(i, sqrt(hSum->GetBinError(i+1)*hSum->GetBinError(i+1)+hSum->GetBinContent(i+1)*hSum->GetBinContent(i+1)*systBck*systBck));
	      }
              gsyst->SetFillColor(12);
              gsyst->SetFillStyle(3345);
              gsyst->SetMarkerSize(0);
              gsyst->SetLineWidth(0);
              gsyst->SetLineColor(kWhite);
	      gsyst->Draw("E2same");
              //TExec *setex1 = new TExec("setex1","gStyle->SetErrorX(0)");
              //setex1->Draw();
	    }
	    else { // plot statistical errors
  	      TGraphAsymmErrors * gsyst = new TGraphAsymmErrors(hSum);
              for (int i = 0; i < gsyst->GetN(); ++i) {
                gsyst->SetPointEYlow (i, hSum->GetBinError(i+1));
                gsyst->SetPointEYhigh(i, hSum->GetBinError(i+1));
	      }
              gsyst->SetFillColor(12);
              gsyst->SetFillStyle(3345);
              gsyst->SetMarkerSize(0);
              gsyst->SetLineWidth(0);
              gsyst->SetLineColor(kWhite);
	      gsyst->Draw("E2same");
	    }

            if(_data && _data->GetSumOfWeights() > 0 && !(_isBlinded)) {
	      bool plotCorrectErrorBars = true;
 	      if(plotCorrectErrorBars == true) {
  		TGraphAsymmErrors * g = new TGraphAsymmErrors(_data);
  		for (int i = 0; i < g->GetN(); ++i) {
               double errScaling = 1.;
               if (_isVarBins) 
                 errScaling = _data->GetBinWidth(i+1);
 	     	   double N = g->GetY()[i] * errScaling;
  	     	   double alpha=(1-0.6827);
  	     	   double L =  (N==0) ? 0  : (ROOT::Math::gamma_quantile(alpha/2,N,1.));
  	     	   double U =  (N==0) ?  ( ROOT::Math::gamma_quantile_c(alpha,N+1,1.) ) :
  	     	      ( ROOT::Math::gamma_quantile_c(alpha/2,N+1,1.) );
  	     	   g->SetPointEYlow(i,double(N/errScaling)-L/errScaling);
		   if(N >= 0)
  	     	     g->SetPointEYhigh(i, U/errScaling-double(N/errScaling));
		   else
		     g->SetPointEYhigh(i, 0.0);
  		}
  		g->Draw("P");
	      }
	      else {
	        if (!_isBlinded) _data->Draw("ep,same");
	      }
            }
	    
            Float_t theMax = hstack->GetMaximum();
            Float_t theMin = hstack->GetMinimum();

            if (_hist[iTop]) {
                if (_hist[iTop]->GetMaximum() > theMax) theMax = _hist[iTop]->GetMaximum();
                if (_hist[iTop]->GetMinimum() < theMin) theMin = _hist[iTop]->GetMinimum();
            }

            if (_data) {

                Float_t dataMax = GetMaximumIncludingErrors(_data);

                if (dataMax > theMax) theMax = dataMax;
            }

            if (gPad->GetLogy()) {
                float maxFactor = 300;
                float minLimit = 0.5;
                if (_isVarBins) {
                  maxFactor = 30;
                  minLimit = 0.005;
                }
            	hstack->SetMaximum(maxFactor * theMax);
            	hstack->SetMinimum(TMath::Max(0.9 * theMin,minLimit));
            } else {
              hstack->SetMaximum(1.55 * theMax);
            }

            if(_breakdown) {
                THStackAxisFonts(hstack, "x", _xLabel.Data());
                THStackAxisFonts(hstack, "y", "Events / bin");
                hstack->GetHistogram()->LabelsOption("v");
            } else {
                THStackAxisFonts(hstack, "x", TString::Format("%s [%s]",_xLabel.Data(),_units.Data()));
                if(_units.Sizeof() == 1) {
                    THStackAxisFonts(hstack, "x", _xLabel.Data());
                    THStackAxisFonts(hstack, "y", "Events / bin");
                } else {
                    THStackAxisFonts(hstack, "x", TString::Format("%s [%s]",_xLabel.Data(),_units.Data()));
                    if(_data->GetBinWidth(0) < 1) THStackAxisFonts(hstack, "y", TString::Format("Events / %.1f %s", _data->GetBinWidth(0),_units.Data()));
		    else                          THStackAxisFonts(hstack, "y", TString::Format("Events / %.0f %s", _data->GetBinWidth(0),_units.Data()));
                }
            }

            size_t j=0;
            //TString higgsLabel = " HWW";
            //higgsLabel.Form("W^{#pm}W^{#pm}jj");

            TString fromW = " W #rightarrow j";
            TString fromQ = " u/d/s #rightarrow j";
            TString fromB = " b/c #rightarrow j";
            TString fromG = " g #rightarrow j";
            if     (_nsel == 1) {
	      fromW = " Top (W #rightarrow j)";
	      fromQ = " Top (q #rightarrow j)";
	      fromB = " Top (b/c #rightarrow j)";
	      fromG = " Top (g #rightarrow j)";
	    }
            else if(_nsel == 2) {
	      fromW = " V+jets (W #rightarrow j)";
	      fromQ = " V+jets (q #rightarrow j)";
	      fromB = " V+jets (b/c #rightarrow j)";
	      fromG = " V+jets (g #rightarrow j)";
	    }

            if(                   _data->GetSumOfWeights() > 0)         { DrawLegend(xPos[j], 0.84 - yOff[j]*_yoffset, _data          , " Data",     	         "lp"); j++; }
            if(_hist[iZjets]   &&_hist[iZjets]->GetSumOfWeights() > 0)  { DrawLegend(xPos[j], 0.84 - yOff[j]*_yoffset, _hist[iZjets]  , " Z+jets",   	         "f" ); j++; }
            if(_hist[iWjets]   &&_hist[iWjets]->GetSumOfWeights() > 0)  { DrawLegend(xPos[j], 0.84 - yOff[j]*_yoffset, _hist[iWjets]  , " W+jets",   	         "f" ); j++; }
            if(_hist[iTop  ]   &&_hist[iTop  ]->GetSumOfWeights() > 0)  { DrawLegend(xPos[j], 0.84 - yOff[j]*_yoffset, _hist[iTop  ]  , " Top",      	         "f" ); j++; }
            if(_hist[iDiboson] &&_hist[iDiboson]->GetSumOfWeights() > 0){ DrawLegend(xPos[j], 0.84 - yOff[j]*_yoffset, _hist[iDiboson], " Diboson",       	     "f" ); j++; }
            if(_hist[iGjets]   &&_hist[iGjets]->GetSumOfWeights() > 0)  { DrawLegend(xPos[j], 0.84 - yOff[j]*_yoffset, _hist[iGjets]  , " #gamma+jets",      	 "f" ); j++; }
            if(_hist[ifromW]   &&_hist[ifromW]->GetSumOfWeights() > 0)  { DrawLegend(xPos[j], 0.84 - yOff[j]*_yoffset, _hist[ifromW]  , Form("%s",fromW.Data()), "f" ); j++; }
            if(_hist[ifromQ]   &&_hist[ifromQ]->GetSumOfWeights() > 0)  { DrawLegend(xPos[j], 0.84 - yOff[j]*_yoffset, _hist[ifromQ]  , Form("%s",fromQ.Data()), "f" ); j++; }
            if(_hist[ifromB]   &&_hist[ifromB]->GetSumOfWeights() > 0)  { DrawLegend(xPos[j], 0.84 - yOff[j]*_yoffset, _hist[ifromB]  , Form("%s",fromB.Data()), "f" ); j++; }
            if(_hist[ifromG]   &&_hist[ifromG]->GetSumOfWeights() > 0)  { DrawLegend(xPos[j], 0.84 - yOff[j]*_yoffset, _hist[ifromG]  , Form("%s",fromG.Data()), "f" ); j++; }
            if(_hist[iOther]   &&_hist[iOther]->GetSumOfWeights() > 0)  { DrawLegend(xPos[j], 0.84 - yOff[j]*_yoffset, _hist[iOther]  , " Other bkgs.",          "f" ); j++; }
            if(_hist[iHiggs]   &&_hist[iHiggs]->GetSumOfWeights() > 0)  { DrawLegend(xPos[j], 0.84 - yOff[j]*_yoffset, _hist[iHiggs]  , "H_{125}#rightarrowinv.","f" ); j++; }

            return hstack->GetHistogram();
        }
        void setLabel(const TString &s) { _xLabel = s; }
        void setUnits(const TString &s) { _units = s; }
        void setBreakdown(const bool &b = true) { _breakdown = b; }
        void setLumiLabel (const std::string &s) {
            _lumiLabel = new TLatex (0.95, 0.93, TString (s));
            _lumiLabel->SetNDC ();
            _lumiLabel->SetTextAlign (30);
            _lumiLabel->SetTextFont (42);
            _lumiLabel->SetTextSize (_tsize);
        }

    private: 
        std::vector<TH1D*> _hist;
        TH1D* _data;

        //MWL
        TLatex * _lumiLabel;     //PG label with the centre of mass energy and lumi info
        TString  _xLabel;
        TString  _units;
        bool     _isBlinded;
        bool     _breakdown;
        int      _mass;
        int      _nsel;
        bool    _isHWWOverlaid;
        bool    _isVarBins;

};
