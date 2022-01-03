#include "Plots.h"
#include "SimpleFits/FitSoftware/interface/Logger.h"

#include "THStack.h"
#include "TLegend.h"
#include "TStyle.h"
#include "TLatex.h"
#include "TFile.h"
#include "TGaxis.h"
#include "TView.h"
#include "TMath.h"
#include <cmath>
#include <iostream>

int Plots::plotLabel = 0;
TString Plots::File_;
std::vector<TString> Plots::HistogramNames_;

Plots::Plots() :
		doscale(false), dooneprofile(false) {

}

Plots::~Plots() {
}

void Plots::Set_Plot_Type(TString style, TString label) {
	style.ToLower();
	label.ToLower();
	if (label.Contains("internal"))
		plotLabel = cmsInternal;
	if (label.Contains("private"))
		plotLabel = cmsPrivate;
	if (label.Contains("preliminary"))
		plotLabel = cmsPreliminary;
	if (label.Contains("public"))
		plotLabel = cmsPublic;
	if (style.Contains("style1")) {
		CMSStyle1();
	} else if (style.Contains("style2")) {
		CMSStyle2();
	} else {
		CMSStyle2();
	}
}

void Plots::Plot1D(std::vector<TH1D> histo, std::vector<int> colour, std::vector<TString> legend) {
	std::vector<std::vector<TH1D> > histos;
	histos.push_back(histo);
	Plot1D(histos, colour, legend);
}

void Plots::SaveHistograms(TString File, std::vector<TString> HistogramNames) {
	File_ = File;
	HistogramNames_ = HistogramNames;
}

void Plots::Plot1D(std::vector<std::vector<TH1D> > histo, std::vector<int> colour, std::vector<TString> legend) {

	Logger(Logger::Verbose) << "Create 1D plots" << std::endl;
	TCanvas canv("canv", "canv", 200, 10, 750, 750);
	canv.Update();
	canv.Clear();
	canv.Update();
	TLatex latex;
	latex.SetTextSize(0.03);
	latex.SetNDC();

	// loop over all plots
	for (unsigned int i_plot = 0; i_plot < histo.size(); i_plot++) {
		Logger(Logger::Debug) << "N(plots)= " << histo.size() << ", i_plot= " << i_plot << ", color.size()= " << colour.size() << ", N(histos)=  " << histo.at(i_plot).size() << std::endl;

		// check that histos and colors exist
		if (histo.at(i_plot).size() > 0 && histo.at(i_plot).size() == colour.size()) {

			// === basic settings ===
			// legend
			TLegend leg(0.25, 0.72, 0.9, 0.925);
			leg.SetBorderSize(0);
			leg.SetFillStyle(4000);
			leg.SetFillColor(0);
			leg.SetTextSize(0.03);
			leg.SetMargin(0.15);
			leg.SetNColumns(3);
			leg.SetColumnSeparation(0.05);
			// create the histogram containing sum of all other histograms
			TH1D Total("Total", "Total", histo.at(i_plot).at(0).GetNbinsX(), histo.at(i_plot).at(0).GetXaxis()->GetXmin(), histo.at(i_plot).at(0).GetXaxis()->GetXmax());
			Total.SetFillStyle(3005);
			Total.SetFillColor(1);
			Total.SetLineColor(18);
			Total.SetMarkerColor(1);
			Total.SetMarkerSize(0.001);
			Total.Sumw2();
			// create stack of histograms
			TString N = "MCStack";
			N += i_plot;
			N += histo.at(i_plot).at(0).GetName();
			THStack MCHistoStack(N, histo.at(i_plot).at(0).GetName());
			// container for the sum of integrals
			double MC_Integral(0);
			double Data_Integral(0);
			// flags to check if histogram should be combined with previous histgram
			// histograms are combined if they have the same name (i.e. legend)
			bool isNewHist = true;
			TString legOfPreviousHist = "";
			// first histogram is expected to be data, and drawn as black dots
			histo.at(i_plot).at(0).SetLineColor(colour.at(0));
			histo.at(i_plot).at(0).SetMarkerColor(colour.at(0));

			// === settings for individual histograms ===
			Logger(Logger::Debug) << "Setup plot " << histo.at(i_plot).at(0).GetTitle() << std::endl;

			// set colours and legends of histograms
			for (int i_hist = histo.at(i_plot).size() - 1; i_hist >= 0; i_hist--) {
				// don't draw histogram if colour is set to -999
				if( colour.at(i_hist) == -999 ) continue;

				// set colors according to values in "colour" vector
				histo.at(i_plot).at(i_hist).SetLineColor(colour.at(i_hist));
				histo.at(i_plot).at(i_hist).SetFillColor(colour.at(i_hist));

				// combine histograms if they have the same legend
				if (i_hist + 1 < histo.at(i_plot).size()) { // first element is always a new histogram
					if (legend.at(i_hist) != legOfPreviousHist) {
						legOfPreviousHist = legend.at(i_hist);
						isNewHist = true;
					}
				}
				// draw black line of histogram only if it should NOT be combined with previous histogram, and if it contains some entries
				if (isNewHist && histo.at(i_plot).at(i_hist).Integral() > 0) {
					histo.at(i_plot).at(i_hist).SetLineColor(1);
					isNewHist = false;
				}
			}

			// calculate data and MC integrals and add legends entries
			double combinedIntegral(0); // sum of integrals of combined histograms
			for (unsigned int i_hist = 0; i_hist < histo.at(i_plot).size(); i_hist++) {
				if( colour.at(i_hist) == -999 ) continue; // don't draw if colour is set to -999

				combinedIntegral += histo.at(i_plot).at(i_hist).Integral();
				// data histogram
				if (i_hist == 0) {
					Data_Integral += histo.at(i_plot).at(0).Integral();
					leg.AddEntry(&histo.at(i_plot).at(i_hist), legend.at(i_hist), "pe");
					combinedIntegral = 0;
				}
				// MC histograms
				else {
					if ( ( (colour.size() <= i_hist +1) // this is the last histogram
							|| (colour.at(i_hist) != colour.at(i_hist + 1)) ) // OR this is an individual histogram
							&& combinedIntegral >= 0) { // AND this histogram has content
						histo.at(i_plot).at(i_hist).SetLineColor(1);
						leg.AddEntry(&histo.at(i_plot).at(i_hist), legend.at(i_hist), "F");
						combinedIntegral = 0;
					}
					if (histo.at(i_plot).at(i_hist).Integral() > 0) {
						// this histogram has content, and may or may not be combined with the next one
						MCHistoStack.Add(&histo.at(i_plot).at(i_hist));
						Total.Add(&histo.at(i_plot).at(i_hist));
					}
					MC_Integral += histo.at(i_plot).at(i_hist).Integral();
					Logger(Logger::Debug) << histo.at(i_plot).at(i_hist).GetTitle() << " " << histo.at(i_plot).at(i_hist).Integral() << " i_hist= " << i_hist << std::endl;
				}
			}

			// scale MC to data for shape comparisons
			// currently the doscale variable is hardcoded to FALSE in the constructor
			if (doscale) {
				for (unsigned int i_hist = 0; i_hist < histo.at(i_plot).size(); i_hist++) {
					if( colour.at(i_hist) == -999 ) continue; // don't draw if colour is set to -999
					if (i_hist != 0 && Data_Integral > 0 && MC_Integral > 0)
						histo.at(i_plot).at(i_hist).Scale(Data_Integral / MC_Integral);
				}
				if (Data_Integral > 0 && MC_Integral > 0)
					Total.Scale(Data_Integral / MC_Integral);
			}

			// === draw the plot ===
			Logger(Logger::Verbose) << "Start drawing ... " << std::endl;

			// loop to create linear (linlog = 0) and logarithmic (linlog = 1) plots
			for (unsigned int linlog = 0; linlog < 2; linlog++) {
				// linear y-axis for linlog = 0
				canv.Clear();
				canv.SetLogy(0);
				if (linlog == 1) {
					// logarithmic y-axis for linlog = 1
					canv.Clear();
					canv.SetLogy(1);
				}

				// draw in case there is data
				if (Data_Integral > 0) {
					// set y-axis range based on data histogram
					double max = histo.at(i_plot).at(0).GetBinContent(histo.at(i_plot).at(0).GetMaximumBin());
					if (linlog == 1) {
						histo.at(i_plot).at(0).SetMinimum(0.01);
						histo.at(i_plot).at(0).SetMaximum(100 * max);
					} else {
						histo.at(i_plot).at(0).SetMinimum(0);
						histo.at(i_plot).at(0).SetMaximum(1.7 * max);
					}
					// draw data histogram
					histo.at(i_plot).at(0).Draw("E");
					// draw MC stack
					MCHistoStack.Draw("Histsame");
					// additionally draw sum of all MC samples
					Total.Draw("E2same");
					// draw data again (to bring it to the top)
					histo.at(i_plot).at(0).Draw("Esame");
					histo.at(i_plot).at(0).Draw("axissame");

					// add Kolmogorov test result to plot
					TString hTitle = histo.at(i_plot).at(0).GetName();
					if (hTitle.Contains("KSTest")) {
						double KS = histo.at(i_plot).at(0).KolmogorovTest(&Total);
						hTitle = histo.at(i_plot).at(0).GetTitle();
						hTitle += " KS=";
						hTitle += KS;
						latex.DrawLatex(0.225, 0.95, hTitle);
					}

					Logger(Logger::Debug) << "Integrals of Data: " << histo.at(i_plot).at(0).Integral() << ", MC: " << Total.Integral() << std::endl;
				}

				// draw in case this is MC only
				else {
					MCHistoStack.Draw("Hist");
					double max = Total.GetBinContent(Total.GetMaximumBin());
					if (linlog == 1) {
						histo.at(i_plot).at(0).SetMinimum(0.01);
						if (max != 0) histo.at(i_plot).at(0).SetMaximum(100 * max);
						else histo.at(i_plot).at(0).SetMaximum(100);
					} else {
						histo.at(i_plot).at(0).SetMinimum(0);
						histo.at(i_plot).at(0).SetMaximum(1.7 * max);
					}
					histo.at(i_plot).at(0).Draw("E");
					MCHistoStack.Draw("Histsame");
					MCHistoStack.Draw("axissame");
					Total.Draw("E2same");
				}

				// draw the legend
				leg.Draw();
				canv.Update();

				// give plot a proper name
				TString filename = histo.at(i_plot).at(0).GetName();
				if (linlog == 1) {
					filename += "_log";
				}
				filename += "_index_";
				filename += i_plot;
				// filename including folder structure and type
				TString EPSName = "EPS/";
				EPSName += filename;
				EPSName += ".eps";

				// save plot only if it is linear OR ( logarithmic AND ( Nminus plot OR Accumdist plot) )
				if (linlog == 0 || !filename.Contains("Nminus") || filename.Contains("Accumdist")) {
					canv.Print(EPSName);
				}

				// save some histograms in extra root file
				// this code does nothing unless method "SaveHistograms" was called before
				if (linlog == 0) {
					for (unsigned int a = 0; a < HistogramNames_.size(); a++) {
						if (histo.at(i_plot).at(0).GetName() == HistogramNames_.at(a)) {
							TFile f(File_ + "_EXTRA.root", "RECREATE");
							TString n = HistogramNames_.at(a);
							n.ReplaceAll("Data", "");
							unsigned int s = histo.at(i_plot).size() - 1;
							histo.at(i_plot).at(0).Write((histo.at(i_plot).at(0)).GetName());
							histo.at(i_plot).at(s).Write(n + "sig");
							for (unsigned int t = 2; t < s; t++) {
								histo.at(i_plot).at(2).Add(&histo.at(i_plot).at(t));
							}
							histo.at(i_plot).at(2).Write(n + "background");
							f.Close();
						}
					}
				}
			}
		}
		else{
			Logger(Logger::Error) << "List of histograms empty OR the list of colors is not of equal length than the list of histograms" << std::endl;
		}
	}
	Logger(Logger::Verbose) << "done. " << std::endl;
}

void Plots::Plot1DSignificance(std::vector<TH1D> histo, bool gt, bool lt, std::vector<int> colour, std::vector<TString> legend) {
	// plots significance: assumes the histograms have been normalized
	Logger(Logger::Verbose) << "Start" << std::endl;
	gStyle->SetPadTopMargin(0.05);
	gStyle->SetPadBottomMargin(0.20);
	gStyle->SetPadLeftMargin(0.175);
	gStyle->SetPadRightMargin(0.05);
	TCanvas c("c", "c", 200, 10, 750, 750);
	c.Update();
	c.Clear();
	c.Update();

	histo.at(0).Reset();
	if (histo.size() > 0) {
		for (int i_bin = 1; i_bin <= histo.at(0).GetNbinsX(); i_bin++) {
			if (gt) {
				double sig = histo.at(histo.size() - 1).Integral(i_bin, histo.at(0).GetNbinsX() + 1);
				double bkg(0);
				for (unsigned int i = 1; i < histo.size() - 1; i++) {
					bkg += histo.at(i).Integral(i_bin, histo.at(0).GetNbinsX() + 1);
				}
				if (/*1+sig/bkg >0 && */(sig + bkg) != 0) {
					histo.at(0).SetBinContent(i_bin, sig / sqrt(sig + bkg)); //(b,sqrt(2)*sqrt((sig+bkg)*TMath::Log(1+sig/bkg)-sig));//(b,sig/sqrt(sig+bkg));
				} else {
					histo.at(0).SetBinContent(i_bin, 0);
				}
			}
			if (lt) {
				double sig = histo.at(histo.size() - 1).Integral(0, i_bin);
				double bkg(0);
				for (unsigned int i = 1; i < histo.size() - 1; i++) {
					bkg += histo.at(i).Integral(0, i_bin);
				}
				if (/*1+sig/bkg >0 && */sig + bkg != 0) {
					histo.at(0).SetBinContent(i_bin, sig / sqrt(sig + bkg)); //(b,sqrt(2)*sqrt((sig+bkg)*TMath::Log(1+sig/bkg)-sig));
				} else {
					histo.at(0).SetBinContent(i_bin, 0);
				}
			}
		}

		histo.at(0).SetYTitle("Signal/#sqrt{Signal+Bkg}"); //("#sqrt{2((S+B)ln(1+S/B)-S)}");//("Signal/#sqrt{Signal+Bkg}");
		histo.at(0).Draw("hist");
		histo.at(0).Draw("axissame");

		TString name = histo.at(0).GetName();
		name += "_sig";
		if (lt)
			name += "lt";
		if (gt)
			name += "gt";
		name += "_index_";
		name += 0;
		TString EPSName = "EPS/";
		EPSName += name;
		EPSName += ".eps";

		if (!name.Contains("Accumdist")) {
			c.Print(EPSName);
			int bmax = histo.at(0).GetMaximumBin();
			if (gt)
				bmax += 1;
			Logger(Logger::Debug) << "Opt. Sig. " << EPSName << " Sig: " << histo.at(0).GetMaximum() << " bin: " << bmax << std::endl;
		}
	}

	Logger(Logger::Verbose) << "done" << std::endl;
}

void Plots::Plot1Dsigtobkg(std::vector<TH1D> histo, bool gt, bool lt, std::vector<int> colour, std::vector<TString> legend) {
	// plots significance: assumes the histograms have been normilized
	Logger(Logger::Verbose) << "Start" << std::endl;
	gStyle->SetPadTopMargin(0.05);
	gStyle->SetPadBottomMargin(0.20);
	gStyle->SetPadLeftMargin(0.175);
	gStyle->SetPadRightMargin(0.05);
	TCanvas c("c", "c", 200, 10, 750, 750);
	c.Update();
	c.Clear();
	c.Update();

	histo.at(0).Reset();
	if (histo.size() > 2) {
		for (int i_bin = 1; i_bin <= histo.at(0).GetNbinsX(); i_bin++) {
			if (gt) {
				double sig = histo.at(histo.size() - 1).Integral(i_bin, histo.at(0).GetNbinsX() + 1);
				double bkg(0);
				for (unsigned int i_hist = 1; i_hist < histo.size() - 1; i_hist++) {
					if( colour.at(i_hist) == -999 ) continue; // don't draw if colour is set to -999
					bkg += histo.at(i_hist).Integral(i_bin, histo.at(0).GetNbinsX() + 1);
				}
				if ((bkg + sig) != 0) {
					histo.at(0).SetBinContent(i_bin, sig / (sig + bkg));
				} else {
					histo.at(0).SetBinContent(i_bin, 0);
				}
			}
			if (lt) {
				double sig = histo.at(histo.size() - 1).Integral(0, i_bin);
				double bkg(0);
				for (unsigned int i_hist = 1; i_hist < histo.size() - 1; i_hist++) {
					if( colour.at(i_hist) == -999 ) continue; // don't draw if colour is set to -999
					bkg += histo.at(i_hist).Integral(0, i_bin);
				}
				if (sig + bkg != 0) {
					histo.at(0).SetBinContent(i_bin, sig / (sig + bkg));
				} else {
					histo.at(0).SetBinContent(i_bin, 0);
				}
			}
		}

		for (int i_bin = 1; i_bin <= histo.at(0).GetNbinsX(); i_bin++) {
			double total = histo.at(histo.size() - 1).Integral();
			if (gt) {
				double sig = histo.at(histo.size() - 1).Integral(i_bin, histo.at(0).GetNbinsX() + 1);
				if (total != 0) {
					histo.at(1).SetBinContent(i_bin, sig / total);
				} else {
					histo.at(1).SetBinContent(i_bin, 0);
				}
			}
			if (lt) {
				double sig = histo.at(histo.size() - 1).Integral(0, i_bin);
				if (total != 0) {
					histo.at(1).SetBinContent(i_bin, sig / total);
				} else {
					histo.at(1).SetBinContent(i_bin, 0);
				}
			}
		}

		histo.at(0).SetYTitle("Purity (Black)/Efficiency (Red)");
		histo.at(0).SetMinimum(0.0);
		histo.at(0).SetMaximum(1.0);
		histo.at(0).Draw("hist");
		histo.at(0).Draw("axissame");
		histo.at(1).SetLineColor(2);
		histo.at(1).Draw("histsame");

		TString name = histo.at(0).GetName();
		name += "_sigtobkg";
		if (lt)
			name += "lt";
		if (gt)
			name += "gt";
		name += "_index_";
		name += 0;
		TString EPSName = "EPS/";
		EPSName += name;
		EPSName += ".eps";
		if (!name.Contains("Accumdist")) {
			c.Print(EPSName);
		}
	}
	Logger(Logger::Verbose) << "Plots::Plot1Dsigtobkg done" << std::endl;
}

void Plots::Plot1D_DataMC_Compare(std::vector<TH1D> histo, std::vector<int> colour, std::vector<TString> legend) {
	Logger(Logger::Verbose) << "Start" << std::endl;
	gStyle->SetPadTopMargin(0.05);
	gStyle->SetPadBottomMargin(0.20);
	gStyle->SetPadLeftMargin(0.175);
	gStyle->SetPadRightMargin(0.05);
	TCanvas c("c", "c", 200, 10, 750, 750);
	c.Clear();

	if (histo.size() > 4) {
		for (unsigned int i_hist = 0; i_hist < histo.size(); i_hist++) {
			if (i_hist > 1) {
				if( colour.at(i_hist) == -999 ) continue; // don't draw if colour is set to -999
				histo.at(1).Add(&histo.at(i_hist), 1);
			}
		}
		if (doscale) {
			Logger(Logger::Debug) << "Comparison Integral Scale:" << histo.at(0).Integral() << " " << histo.at(1).Integral() << std::endl;
			histo.at(1).Scale(histo.at(0).Integral() / histo.at(1).Integral());
		}
		for (int b = 1; b <= histo.at(0).GetNbinsX(); b++) {
			double data = histo.at(0).GetBinContent(b);
			double dataerror = histo.at(0).GetBinError(b);
			double mc = histo.at(1).GetBinContent(b);
			double mcerror = histo.at(1).GetBinError(b);
			if (mc != 0) {
				histo.at(0).SetBinContent(b, data / mc);
				histo.at(0).SetBinError(b, dataerror / mc);
				histo.at(2).SetBinContent(b, 1);
				histo.at(2).SetBinError(b, mcerror / mc);

			} else {
				histo.at(0).SetBinContent(b, 0);
				histo.at(0).SetBinError(b, 0);
				histo.at(2).SetBinContent(b, 0);
				histo.at(2).SetBinError(b, 2);
			}
			histo.at(1).SetBinContent(b, data - mc);
			histo.at(1).SetBinError(b, dataerror);
			histo.at(3).SetBinContent(b, 0);
			histo.at(3).SetBinError(b, mcerror);

		}

		histo.at(0).SetYTitle("Data/MC");
		histo.at(0).SetMinimum(0.0);
		histo.at(0).SetMaximum(2.0);
		histo.at(0).Draw("E");
		histo.at(0).Draw("axissame");
		histo.at(2).SetFillColor(5);
		histo.at(2).SetLineColor(5);
		histo.at(2).SetMarkerColor(5);
		histo.at(2).SetMarkerSize(0.0001);
		histo.at(2).Draw("E2same");
		histo.at(0).Draw("Esame");
		TString name = histo.at(0).GetName();
		name += "_sig";
		name += "_index_";
		name += 0;
		TString EPSName = "EPS/";
		EPSName += name;
		EPSName += ".eps";
		if (!name.Contains("Accumdist")) {
			c.Print(EPSName);
		}
		///////////////////////////////////////////////////////////////////////////////////////////////////////
		histo.at(1).SetYTitle("Data-MC");
		histo.at(1).Draw("E");
		histo.at(1).Draw("axissame");
		histo.at(3).SetFillColor(5);
		histo.at(3).SetLineColor(5);
		histo.at(3).SetMarkerColor(5);
		histo.at(3).SetMarkerSize(0.0001);
		histo.at(3).Draw("E2same");
		histo.at(1).Draw("Esame");

		name = histo.at(0).GetName();
		name += "_sigtobkg";
		name += "_index_";
		name += 0;
		EPSName = "EPS/";
		EPSName += name;
		EPSName += ".eps";
		if (!name.Contains("Accumdist")) {
			c.Print(EPSName);
		}

	}
	Logger(Logger::Verbose) << "Finished" << std::endl;
}

void Plots::Plot2D(std::vector<TH2D> histo, std::vector<int> colour, std::vector<TString> legend) {
	Logger(Logger::Verbose) << "Create 2D plots" << std::endl;
	gStyle->SetPadTopMargin(0.10);
	gStyle->SetPadBottomMargin(0.22);
	gStyle->SetPadLeftMargin(0.180);
	gStyle->SetPadRightMargin(0.175);
//	gStyle->SetPalette(53,0);
	Logger(Logger::Debug) << "Starting Plot2D" << histo.size() << " " << colour.size() << " " << legend.size() << std::endl;
	if (histo.size() > 0) {
		std::vector<int> types;
		TString LegType = "";
		double integral = 0;
		for (unsigned int i_hist = 0; i_hist < histo.size(); i_hist++) {
			if( colour.at(i_hist) == -999 ) continue; // don't draw if colour is set to -999
			Logger(Logger::Debug) << "Legend " << legend.at(i_hist) << " Integral " << histo.at(i_hist).Integral() << std::endl;
			if (LegType != legend.at(i_hist)) {
				LegType = legend.at(i_hist);
				integral = 0;
				for (unsigned int j = i_hist; j < histo.size(); j++) {
					if (LegType == legend.at(j))
						integral += histo.at(j).Integral();
				}
				if (integral != 0) {
					types.push_back((i_hist != 0) ? i_hist + 1 : i_hist); // +1 as later on the TotalMC histo is added
					Logger(Logger::Debug) << "adding one type to types list" << std::endl;
				}
				integral = 0;
			}
		}

		LegType = "";
		std::vector<TH2D>::iterator it;
		it = histo.begin();
		TString tempname = histo.at(0).GetName();
		tempname += "TotalMC";
		histo.insert(it + 1, *((TH2D*) histo.at(0).Clone(tempname.Data())));
		if (types.size() > 0)
			types.insert(types.begin() + 1, 1);
		std::vector<TString>::iterator itleg = legend.begin();
		legend.insert(itleg + 1, "Total MC");
		std::vector<int>::iterator itcol = colour.begin();
		colour.insert(itcol + 1, 55);
		histo.at(1).Reset();
		for (unsigned int i_hist = 2; i_hist < histo.size(); i_hist++) {
			if( colour.at(i_hist) == -999 ) continue; // don't draw if colour is set to -999
			histo.at(1).Add(&histo.at(i_hist));
		}
		Logger(Logger::Debug) << "n(types) = " << types.size() << std::endl;
		LegType = "null";

		for (unsigned int i_hist = 0; i_hist < histo.size(); i_hist++) {
			if( colour.at(i_hist) == -999 ) continue; // don't draw if colour is set to -999
			Logger(Logger::Debug) << "histo i" << i_hist << std::endl;
			if (LegType != legend.at(i_hist)) {
				LegType = legend.at(i_hist);
				for (unsigned int j = i_hist + 1; j < histo.size(); j++) {
					Logger(Logger::Debug) << "histo j" << j << std::endl;
					if (legend.at(j) == legend.at(i_hist)) {
						histo.at(i_hist).Add(&histo.at(j));
						Logger(Logger::Debug) << "Added histo " << j << " to histo " << i_hist << std::endl;
					}
				}
				Logger(Logger::Debug) << "legend " << i_hist << " " << histo.at(i_hist).Integral() << std::endl;
				if (histo.at(i_hist).Integral() != 0) {
					Logger(Logger::Debug) << "legend " << i_hist << histo.at(i_hist).Integral() << std::endl;
				}
			}
		}

		// loop to allow splitting of TH2D histograms on multiple canvases
		unsigned int nPlotted(0);
		while (nPlotted < types.size()) {
			Logger(Logger::Debug) << "nPlotted = " << nPlotted << ", ntypes = " << types.size() << std::endl;
			// split canvas for TH2D histograms
			TCanvas c("c", "c", 200, 10, 750, 750);
			c.Clear();

			if (types.size() - nPlotted == 1) {
				Logger(Logger::Debug) << "Canvas: Divide(1,1)" << std::endl;
				c.Divide(1, 1);
			} else if (types.size() - nPlotted == 2 || types.size() - nPlotted == 3 || types.size() - nPlotted == 4) {
				Logger(Logger::Debug) << "Canvas: Divide(2,2)" << std::endl;
				c.Divide(2, 2);
			} else {
				Logger(Logger::Debug) << "Canvas: Divide(2,3)" << std::endl;
				// plot only up to 6 TH2D in one canvas
				c.Divide(2, 3);
			}

			for (unsigned int i = nPlotted; (i < nPlotted + 6) && (i < types.size()); i++) {
				int type = types.at(i);
				Logger(Logger::Debug) << "loop: " << i << ", type = " << types.at(i) << " legend.at(type) = " << legend.at(type) << std::endl;
				c.cd(i - nPlotted + 1);
				Logger(Logger::Debug) << "Draw histogram " << type << " into subcanvas " << i - nPlotted + 1 << " with Integral = " << histo.at(type).Integral() << std::endl;
				histo.at(type).Draw("colz");
				c.Update();
				TLatex l;
				l.SetNDC();
				l.DrawLatex(0.25, 0.90, legend.at(type));
			}

			TString name = histo.at(0).GetName();
			TString EPSName = "EPS/";
			c.Update();
			name = histo.at(0).GetName();
			name += "_index_";
			name += 0;
			EPSName = "EPS/";
			EPSName += name;
			EPSName += nPlotted;
			EPSName += ".eps";
			c.Print(EPSName);

			nPlotted += 6;
		}
	}
}

void Plots::Plot3D(std::vector<TH3F> histo, std::vector<int> colour, std::vector<TString> legend) {
	Logger(Logger::Verbose) << "Create 3D plots" << std::endl;
	gStyle->SetPadTopMargin(0.10);
	gStyle->SetPadBottomMargin(0.22);
	gStyle->SetPadLeftMargin(0.180);
	gStyle->SetPadRightMargin(0.175);
	gStyle->SetTitleXOffset(2.4);
	gStyle->SetTitleYOffset(2.4);

	TCanvas c("c", "c", 200, 10, 750, 750);
	if (histo.size() > 0) {
		Logger(Logger::Debug) << "Starting Plot3D" << histo.size() << " " << colour.size() << " " << legend.size() << std::endl;
		c.Clear();
		c.Divide(2, 2);
		TLegend leg(0.1, 0.4, 0.9, 0.9);
		leg.SetBorderSize(0);
		leg.SetFillStyle(4000);
		leg.SetFillColor(0);
		leg.SetTextSize(0.03);
		leg.SetMargin(0.15);
		leg.SetNColumns(1);
		leg.SetColumnSeparation(0.05);
		for (unsigned int p = 0; p < 4; p++) {
			for (unsigned int i_hist = 0; i_hist < histo.size(); i_hist++) {
				if( colour.at(i_hist) == -999 ) continue; // don't draw if colour is set to -999
				histo.at(i_hist).SetMarkerSize(0.25);
				histo.at(i_hist).SetMarkerColor(colour.at(i_hist));
				histo.at(i_hist).GetXaxis()->SetTitleOffset(2.0);
				histo.at(i_hist).GetYaxis()->SetTitleOffset(2.0);
				histo.at(i_hist).GetZaxis()->SetTitleOffset(2.0);
				if (p == 0)
					gPad->SetPhi(30);
				if (p == 1)
					gPad->SetPhi(150);
				if (p == 2)
					gPad->SetPhi(240);
				if (colour.at(i_hist) == 10)
					histo.at(i_hist).SetMarkerColor(16);
				if (colour.at(i_hist) == 414)
					histo.at(i_hist).SetMarkerColor(3);
				if (p != 3) {
					c.cd(p + 1);

					if (i_hist == 0) {
						histo.at(i_hist).Draw();
					} else {
						histo.at(i_hist).Draw("same");
					}
				} else {
					leg.AddEntry(&histo.at(i_hist), legend.at(i_hist), "pe");
				}
			}
		}
		c.cd(4);
		leg.Draw();
		TString name = histo.at(0).GetName();
		TString EPSName = "EPS/";
		c.Update();
		name = histo.at(0).GetName();
		name += "_index_";
		name += 0;
		EPSName = "EPS/";
		EPSName += name;
		EPSName += ".eps";
		c.Print(EPSName);
		Logger(Logger::Debug) << "Histo 3D " << EPSName << std::endl;
	}
}

void Plots::CMSStyle1() {
	Logger(Logger::Info) << "Configuring Plots::CMSStyle1()" << std::endl;
	// use plain black on white colors
	int icol = 0; // WHITE
	gStyle->SetFrameBorderMode(icol);
	gStyle->SetFrameFillColor(icol);
	gStyle->SetCanvasBorderMode(icol);
	gStyle->SetCanvasColor(icol);
	gStyle->SetPadBorderMode(icol);
	gStyle->SetPadColor(icol);
	gStyle->SetStatColor(icol);
	//gStyle->SetFillColor(icol); // don't use: white fill color for *all* objects
	// set the paper & margin sizes
	gStyle->SetPaperSize(20, 26);

	// set margin sizes
	gStyle->SetPadTopMargin(0.05);
	gStyle->SetPadRightMargin(0.05);
	gStyle->SetPadBottomMargin(0.16);
	gStyle->SetPadLeftMargin(0.16);

	// set title offsets (for axis label)
	gStyle->SetTitleXOffset(1.4);
	gStyle->SetTitleYOffset(1.4);

	// use large fonts
	//int font=72; // Helvetica italics
	int font = 42; // Helvetica
	Double_t tsize = 0.05;
	gStyle->SetTextFont(font);

	gStyle->SetTextSize(tsize);
	gStyle->SetLabelFont(font, "x");
	gStyle->SetTitleFont(font, "x");
	gStyle->SetLabelFont(font, "y");
	gStyle->SetTitleFont(font, "y");
	gStyle->SetLabelFont(font, "z");
	gStyle->SetTitleFont(font, "z");

	gStyle->SetLabelSize(tsize, "x");
	gStyle->SetTitleSize(tsize, "x");
	gStyle->SetLabelSize(tsize, "y");
	gStyle->SetTitleSize(tsize, "y");
	gStyle->SetLabelSize(tsize, "z");
	gStyle->SetTitleSize(tsize, "z");

	// use bold lines and markers
	gStyle->SetMarkerStyle(20);
	gStyle->SetMarkerSize(1.2);
	gStyle->SetHistLineWidth(2);
	gStyle->SetLineStyleString(2, "[12 12]"); // postscript dashes

	// get rid of X error bars
	//gStyle->SetErrorX(0.001);
	// get rid of error bar caps
	gStyle->SetEndErrorSize(0.);

	// do not display any of the standard histogram decorations
	gStyle->SetOptTitle(0);
	//gStyle->SetOptStat(1111);
	gStyle->SetOptStat(0);
	//gStyle->SetOptFit(1111);
	gStyle->SetOptFit(0);
	// put tick marks on top and RHS of plots
	gStyle->SetPadTickX(1);
	gStyle->SetPadTickY(1);
	gStyle->SetPalette(1);
	gROOT->ForceStyle(true);

}

void Plots::CMSStyle2() {
	Logger(Logger::Info) << "Configuring Plots::CMSStyle2()" << std::endl;
	// use plain black on white colors
	int icol = 0; // WHITE
	gStyle->SetFrameBorderMode(icol);
	gStyle->SetFrameFillColor(icol);
	gStyle->SetCanvasBorderMode(icol);
	gStyle->SetCanvasColor(icol);
	gStyle->SetPadBorderMode(icol);
	gStyle->SetPadColor(icol);
	gStyle->SetStatColor(icol);
	//gStyle->SetFillColor(icol); // don't use: white fill color for *all* objects

	// set the paper & margin sizes
	gStyle->SetPaperSize(20, 26);

	// set margin sizes
	gStyle->SetPadTopMargin(0.05);
	gStyle->SetPadRightMargin(0.05);
	gStyle->SetPadBottomMargin(0.16);
	gStyle->SetPadLeftMargin(0.16);

	// set title offsets (for axis label)
	gStyle->SetTitleXOffset(1.4);
	gStyle->SetTitleYOffset(1.4);

	// use large fonts
	//int font=72; // Helvetica italics
	int font = 42; // Helvetica
	Double_t tsize = 0.05;
	gStyle->SetTextFont(font);

	gStyle->SetTextSize(tsize);
	gStyle->SetLabelFont(font, "x");
	gStyle->SetTitleFont(font, "x");
	gStyle->SetLabelFont(font, "y");
	gStyle->SetTitleFont(font, "y");
	gStyle->SetLabelFont(font, "z");
	gStyle->SetTitleFont(font, "z");

	gStyle->SetLabelSize(tsize, "x");
	gStyle->SetTitleSize(tsize, "x");
	gStyle->SetLabelSize(tsize, "y");
	gStyle->SetTitleSize(tsize, "y");
	gStyle->SetLabelSize(tsize, "z");
	gStyle->SetTitleSize(tsize, "z");

	// use bold lines and markers
	gStyle->SetMarkerStyle(20);
	gStyle->SetMarkerSize(1.2);
	gStyle->SetHistLineWidth(2);
	gStyle->SetLineStyleString(2, "[12 12]"); // postscript dashes

	// get rid of X error bars
	//gStyle->SetErrorX(0.001);
	// get rid of error bar caps
	gStyle->SetEndErrorSize(0.);

	// do not display any of the standard histogram decorations
	gStyle->SetOptTitle(0);
	//gStyle->SetOptStat(1111);
	gStyle->SetOptStat(0);
	//gStyle->SetOptFit(1111);
	gStyle->SetOptFit(0);

	// put tick marks on top and RHS of plots
	gStyle->SetPadTickX(1);
	gStyle->SetPadTickY(1);
	gStyle->SetPalette(1);
	gROOT->ForceStyle(true);

}

void Plots::CMSLabel(Double_t x, Double_t y, Color_t color) {
	TLatex l; //l.SetTextAlign(12); l.SetTextSize(tsize);
	l.SetNDC();
	l.SetTextFont(72);
	l.SetTextColor(color);

	double delx = 0.115 * 696 * gPad->GetWh() / (472 * gPad->GetWw());

	l.DrawLatex(x, y, "CMS");
	TLatex p;
	p.SetNDC();
	p.SetTextFont(42);
	p.SetTextColor(color);
	if (plotLabel == cmsPrivate) {
		p.DrawLatex(x + delx, y, "Private");
	}
	if (plotLabel == cmsPreliminary) {
		p.DrawLatex(x + delx, y, "Preliminary");
	}
}

