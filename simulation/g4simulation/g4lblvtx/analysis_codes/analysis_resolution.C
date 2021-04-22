/*
   Winston DeGraw (wdegraw@lbl.gov)
   Rey Cruz-Torres (reynier@lbl.gov)
   */

// Forward-declaring functions
void prettyTH1F( TH1F * h1 , int color , int marker , float min , float max );
// ============================================================================================================================================
void analysis_resolution(){
	// -------------------------------------------------------------
	// Some important parameters
	const TString partic = "mu-";	// particle to be studied
	const float Bfield = 3.0;	// [T] Magnetic field
	bool use_widths = true;
	bool update_tab = true;
	// -------------------------
	// Binning
	float eta_bin[] = {0.,0.5,1.0,1.5,2.0,2.5,3.0,3.5,4.0};		const int size_eta_bin = sizeof(eta_bin)/sizeof(*eta_bin);
	float mom_bin[] = {0.,1.,2.,3.,4.,5.,10.,15.,20.,25.,30.};	const int size_mom_bin = sizeof(mom_bin)/sizeof(*mom_bin);
	// -------------------------
	string Bfield_str = Form("%.1f",Bfield);
	Bfield_str.replace(Bfield_str.find("."), sizeof(".") - 1, "_");
	TString tab_name = Form("sigma_eta_%i_p_%i_",size_eta_bin-1,size_mom_bin-1) + partic + Form("_B%.1fT.txt",Bfield);
	// -------------------------------------------------------------
	// Some settings
	TH1::SetDefaultSumw2();
	TH2::SetDefaultSumw2();
	gStyle -> SetOptStat(0);	
	// -------------------------------------------------------------
	// Loading all the needed info from the root file
	TFile * F = new TFile("../data/files_"+Bfield_str+"_T/total_"+partic+"_B_"+Form("%.1f",Bfield)+"T_FastTrackingEval.root");
	TTree * T = (TTree*) F -> Get("tracks");
	float gpx, gpy, gpz, px, py, pz;
	T -> SetBranchAddress("gpx",&gpx);
	T -> SetBranchAddress("gpy",&gpy);
	T -> SetBranchAddress("gpz",&gpz);
	T -> SetBranchAddress("px" ,&px );
	T -> SetBranchAddress("py" ,&py );
	T -> SetBranchAddress("pz" ,&pz );
	int nEntries = T -> GetEntries();
	// -------------------------------------------------------------
	fstream tab;
	float approx_sig_dpp[size_eta_bin-1][size_mom_bin-1] = {{0}};
	float approx_sig_dth[size_eta_bin-1][size_mom_bin-1] = {{0}};
	float approx_sig_dph[size_eta_bin-1][size_mom_bin-1] = {{0}};
	TString temp_str;
	if(use_widths){
		tab.open("tables/"+tab_name);
		if(!tab){cout << "Could not find file '" << tab_name << "'" << endl; use_widths = false; update_tab = true;}
		else{
			cout << "Loading parameters from file '" << tab_name << "'" << endl;
			for(int et = 0 ; et < size_eta_bin-1 ; et++){ for(int p = 0 ; p < size_mom_bin-1 ; p++){tab >> approx_sig_dpp[et][p];}}	//tab >> temp_str;
			for(int et = 0 ; et < size_eta_bin-1 ; et++){ for(int p = 0 ; p < size_mom_bin-1 ; p++){tab >> approx_sig_dth[et][p];}}	//tab >> temp_str;
			for(int et = 0 ; et < size_eta_bin-1 ; et++){ for(int p = 0 ; p < size_mom_bin-1 ; p++){tab >> approx_sig_dph[et][p];}}
		}
		tab.close();
	}

	float approx_sig_dpp_3_0[size_eta_bin-1][size_mom_bin-1] = {0}; float approx_sig_dpp_1_2[size_eta_bin-1][size_mom_bin-1] = {0};
	float approx_sig_dth_3_0[size_eta_bin-1][size_mom_bin-1] = {0}; float approx_sig_dth_1_2[size_eta_bin-1][size_mom_bin-1] = {0};
	float approx_sig_dph_3_0[size_eta_bin-1][size_mom_bin-1] = {0}; float approx_sig_dph_1_2[size_eta_bin-1][size_mom_bin-1] = {0};

	for(int et = 0 ; et < size_eta_bin-1 ; et++){
		for(int p = 0 ; p < size_mom_bin-1 ; p++){
			approx_sig_dpp_3_0[et][p] = 3.0*approx_sig_dpp[et][p]; approx_sig_dpp_1_2[et][p] = 1.2*approx_sig_dpp[et][p];
			approx_sig_dth_3_0[et][p] = 3.0*approx_sig_dth[et][p]; approx_sig_dth_1_2[et][p] = 1.2*approx_sig_dth[et][p];
			approx_sig_dph_3_0[et][p] = 3.0*approx_sig_dph[et][p]; approx_sig_dph_1_2[et][p] = 1.2*approx_sig_dph[et][p];
		}
	}
	// -------------------------------------------------------------
	// Defining histograms
	TH1F *** h1_dpp_p_et_bins = new TH1F**[size_eta_bin-1];	// delta p / p vs. p in eta bins
	TH1F *** h1_dth_p_et_bins = new TH1F**[size_eta_bin-1];	// delta theta vs. p in eta bins
	TH1F *** h1_dph_p_et_bins = new TH1F**[size_eta_bin-1];	// delta phi   vs. p in eta bins
	for(int et = 0 ; et < size_eta_bin-1 ; et++){
		h1_dpp_p_et_bins[et] = new TH1F*[size_mom_bin-1];
		h1_dth_p_et_bins[et] = new TH1F*[size_mom_bin-1];
		h1_dph_p_et_bins[et] = new TH1F*[size_mom_bin-1];
		for(int p = 0 ; p < size_mom_bin-1 ; p++){
			if(use_widths){
				h1_dpp_p_et_bins[et][p] = new TH1F(Form("h1_dpp_p_et_bins_%i_%i",et,p),";dp/p;Counts"         ,100,-approx_sig_dpp_3_0[et][p],approx_sig_dpp_3_0[et][p]);
				h1_dth_p_et_bins[et][p] = new TH1F(Form("h1_dth_p_et_bins_%i_%i",et,p),";d#theta [rad];Counts",100,-approx_sig_dth_3_0[et][p],approx_sig_dth_3_0[et][p]);
				h1_dph_p_et_bins[et][p] = new TH1F(Form("h1_dph_p_et_bins_%i_%i",et,p),";d#phi [rad];Counts"  ,100,-approx_sig_dph_3_0[et][p],approx_sig_dph_3_0[et][p]);
			}
			else{
				h1_dpp_p_et_bins[et][p] = new TH1F(Form("h1_dpp_p_et_bins_%i_%i",et,p),";dp/p;Counts"         ,100,-0.08  ,0.08  );
				h1_dth_p_et_bins[et][p] = new TH1F(Form("h1_dth_p_et_bins_%i_%i",et,p),";d#theta [rad];Counts",100,-0.0014,0.0014);
				h1_dph_p_et_bins[et][p] = new TH1F(Form("h1_dph_p_et_bins_%i_%i",et,p),";d#phi [rad];Counts"  ,100,-0.04  ,0.04  );
			}

			h1_dpp_p_et_bins[et][p] -> SetTitle(Form("%.1f < |#eta| < %.1f, %.1f < p < %.1f GeV/c",eta_bin[et],eta_bin[et],mom_bin[p],mom_bin[p]));
			h1_dth_p_et_bins[et][p] -> SetTitle(Form("%.1f < |#eta| < %.1f, %.1f < p < %.1f GeV/c",eta_bin[et],eta_bin[et],mom_bin[p],mom_bin[p]));
			h1_dph_p_et_bins[et][p] -> SetTitle(Form("%.1f < |#eta| < %.1f, %.1f < p < %.1f GeV/c",eta_bin[et],eta_bin[et],mom_bin[p],mom_bin[p]));
		}
	}
	// -------------------------------------------------------------	
	int color[] = {1,2,62,8,95,52,6,28,209,92,15};
	int marker[] = {20,21,23,24,25,26,27,28,29,30};

	TH1F ** h1_dpp_v_p_et_bins = new TH1F*[size_eta_bin-1];
	TH1F ** h1_dth_v_p_et_bins = new TH1F*[size_eta_bin-1];
	TH1F ** h1_dph_v_p_et_bins = new TH1F*[size_eta_bin-1];

	for(int et = 0 ; et < size_eta_bin-1 ; et++){
		h1_dpp_v_p_et_bins[et] = new TH1F(Form("h1_dpp_v_p_et_bins_%i",et),";p [GeV/c];dp/p [%]"      ,size_mom_bin-1,mom_bin);	prettyTH1F( h1_dpp_v_p_et_bins[et] , color[et] , marker[et] , 0. , 10. );
		h1_dth_v_p_et_bins[et] = new TH1F(Form("h1_dth_v_p_et_bins_%i",et),";p [GeV/c];d#theta [mrad]",size_mom_bin-1,mom_bin);	prettyTH1F( h1_dth_v_p_et_bins[et] , color[et] , marker[et] , 0. , 1.  );
		h1_dph_v_p_et_bins[et] = new TH1F(Form("h1_dph_v_p_et_bins_%i",et),";p [GeV/c];d#phi [mrad]"  ,size_mom_bin-1,mom_bin);	prettyTH1F( h1_dph_v_p_et_bins[et] , color[et] , marker[et] , 0. , 25. );
	}

	TH1F ** h1_dpp_v_et_p_bins = new TH1F*[size_mom_bin-1];
	TH1F ** h1_dth_v_et_p_bins = new TH1F*[size_mom_bin-1];
	TH1F ** h1_dph_v_et_p_bins = new TH1F*[size_mom_bin-1];

	for(int p = 0 ; p < size_mom_bin-1 ; p++){
		h1_dpp_v_et_p_bins[p] = new TH1F(Form("h1_dpp_v_et_p_bins_%i",p),";#eta;dp/p [%]"      ,size_eta_bin-1,eta_bin);	prettyTH1F( h1_dpp_v_et_p_bins[p] , color[p] , marker[p] , 0. , 10. );
		h1_dth_v_et_p_bins[p] = new TH1F(Form("h1_dth_v_et_p_bins_%i",p),";#eta;d#theta [mrad]",size_eta_bin-1,eta_bin);	prettyTH1F( h1_dth_v_et_p_bins[p] , color[p] , marker[p] , 0. , 1.  );
		h1_dph_v_et_p_bins[p] = new TH1F(Form("h1_dph_v_et_p_bins_%i",p),";#eta;d#phi [mrad]"  ,size_eta_bin-1,eta_bin);	prettyTH1F( h1_dph_v_et_p_bins[p] , color[p] , marker[p] , 0. , 25. );
	}

	// -------------------------------------------------------------
	// Declaring other useful variables and functions
	float width_dpp[size_eta_bin-1][size_mom_bin-1] = {{0}};
	float error_dpp[size_eta_bin-1][size_mom_bin-1] = {{0}};
	float width_dth[size_eta_bin-1][size_mom_bin-1] = {{0}};
	float error_dth[size_eta_bin-1][size_mom_bin-1] = {{0}};
	float width_dph[size_eta_bin-1][size_mom_bin-1] = {{0}};
	float error_dph[size_eta_bin-1][size_mom_bin-1] = {{0}};

	TF1 *** f_gaus_dpp = new TF1**[size_eta_bin-1];
	TF1 *** f_gaus_dth = new TF1**[size_eta_bin-1];
	TF1 *** f_gaus_dph = new TF1**[size_eta_bin-1];

	for(int et = 0 ; et < size_eta_bin-1 ; et++){
		f_gaus_dpp[et] = new TF1*[size_mom_bin-1];
		f_gaus_dth[et] = new TF1*[size_mom_bin-1];
		f_gaus_dph[et] = new TF1*[size_mom_bin-1];

		for(int p = 0 ; p < size_mom_bin-1 ; p++){
			if(use_widths){
				f_gaus_dpp[et][p] = new TF1(Form("f_gaus_dpp_%i_%i",et,p),"gaus",-approx_sig_dpp_1_2[et][p],approx_sig_dpp_1_2[et][p]);
				f_gaus_dth[et][p] = new TF1(Form("f_gaus_dth_%i_%i",et,p),"gaus",-approx_sig_dth_1_2[et][p],approx_sig_dth_1_2[et][p]);
				f_gaus_dph[et][p] = new TF1(Form("f_gaus_dph_%i_%i",et,p),"gaus",-approx_sig_dph_1_2[et][p],approx_sig_dph_1_2[et][p]);
			}
			else{
				f_gaus_dpp[et][p] = new TF1(Form("f_gaus_dpp_%i_%i",et,p),"gaus",-0.007 ,0.007 );
				f_gaus_dth[et][p] = new TF1(Form("f_gaus_dth_%i_%i",et,p),"gaus",-0.0007,0.0007);
				f_gaus_dph[et][p] = new TF1(Form("f_gaus_dph_%i_%i",et,p),"gaus",-0.02  ,0.02  );
			}
		}
	}
	// -------------------------------------------------------------
	// Loop over entries of the tree
	for(int ev = 0 ; ev < nEntries ; ev++){
		T -> GetEntry(ev);
		if(ev%1000000==0) cout << "Looping over entry " << ev << " out of " << nEntries << endl;

		// Calculating some variables
		float gtheta = TMath::ACos(gpz/sqrt(gpx*gpx+gpy*gpy+gpz*gpz));	
		float theta = TMath::ACos(pz/sqrt(px*px+py*py+pz*pz));
		float dth = theta - gtheta;

		float geta = -TMath::Log(TMath::Tan(gtheta/2.));

		float p_reco = sqrt(px*px+py*py+pz*pz);
		float p_truth = sqrt(gpx*gpx+gpy*gpy+gpz*gpz);
		float dp_p = (p_reco-p_truth)/p_truth;

		float gphi = TMath::ATan(gpy/gpx);
		float phi = TMath::ATan(py/px);
		float dph = phi - gphi;

		// Filling histograms
		for(int et = 0 ; et < size_eta_bin-1 ; et++){
			if( abs(geta) >  eta_bin[et] &&  abs(geta) <= eta_bin[et+1] ){
				for(int p = 0 ; p < size_mom_bin-1 ; p++){
					if( p_truth > mom_bin[p] && p_truth <= mom_bin[p+1] ){
						h1_dpp_p_et_bins[et][p] -> Fill( dp_p );
						h1_dth_p_et_bins[et][p] -> Fill( dth  );
						h1_dph_p_et_bins[et][p] -> Fill( dph  );
					}	
				}
			}
		}
	}
	// -------------------------------------------------------------
	// Doing fits
	TCanvas ** c_fits_p  = new TCanvas*[size_eta_bin-1];
	TCanvas ** c_fits_th = new TCanvas*[size_eta_bin-1];
	TCanvas ** c_fits_ph = new TCanvas*[size_eta_bin-1];

	for(int et = 0 ; et < size_eta_bin-1 ; et++){
		c_fits_p [et] = new TCanvas(Form("c_fits_p_%i" ,et),Form("dp/p  , %.1f<eta<%.1f",eta_bin[et],eta_bin[et+1]),1000,800);	c_fits_p [et] -> Divide(3,2);
		c_fits_th[et] = new TCanvas(Form("c_fits_th_%i",et),Form("dtheta, %.1f<eta<%.1f",eta_bin[et],eta_bin[et+1]),1000,800);	c_fits_th[et] -> Divide(3,2);
		c_fits_ph[et] = new TCanvas(Form("c_fits_ph_%i",et),Form("dphi  , %.1f<eta<%.1f",eta_bin[et],eta_bin[et+1]),1000,800);	c_fits_ph[et] -> Divide(3,2);

		for(int p = 0 ; p < size_mom_bin-1 ; p++){
			// Momentum resolutions
			c_fits_p [et] -> cd(p+1);
			h1_dpp_p_et_bins[et][p] -> Draw();	h1_dpp_p_et_bins[et][p] -> Fit(Form("f_gaus_dpp_%i_%i",et,p),"R");
			width_dpp[et][p] = f_gaus_dpp[et][p] -> GetParameter(2);
			error_dpp[et][p] = (f_gaus_dpp[et][p] -> GetParError(2))*(f_gaus_dpp[et][p] -> GetChisquare())/(f_gaus_dpp[et][p] -> GetNDF());

			// Theta resolution
			c_fits_th[et] -> cd(p+1);
			h1_dth_p_et_bins[et][p] -> Draw();	h1_dth_p_et_bins[et][p] -> Fit(Form("f_gaus_dth_%i_%i",et,p),"R");
			width_dth[et][p] = f_gaus_dth[et][p] -> GetParameter(2);
			error_dth[et][p] = (f_gaus_dth[et][p] -> GetParError(2))*(f_gaus_dth[et][p] -> GetChisquare())/(f_gaus_dth[et][p] -> GetNDF());

			// Phi resolution
			c_fits_ph[et] -> cd(p+1);
			h1_dph_p_et_bins[et][p] -> Draw();	h1_dph_p_et_bins[et][p] -> Fit(Form("f_gaus_dph_%i_%i",et,p),"R");
			width_dph[et][p] = f_gaus_dph[et][p] -> GetParameter(2);
			error_dph[et][p] = (f_gaus_dph[et][p] -> GetParError(2))*(f_gaus_dph[et][p] -> GetChisquare())/(f_gaus_dph[et][p] -> GetNDF());

			// ----

			h1_dpp_v_p_et_bins[et] -> SetBinContent(p +1,width_dpp[et][p]*100. );
			h1_dth_v_p_et_bins[et] -> SetBinContent(p +1,width_dth[et][p]*1000.);
			h1_dph_v_p_et_bins[et] -> SetBinContent(p +1,width_dph[et][p]*1000.);
			h1_dpp_v_et_p_bins[ p] -> SetBinContent(et+1,width_dpp[et][p]*100. );
			h1_dth_v_et_p_bins[ p] -> SetBinContent(et+1,width_dth[et][p]*1000.);
			h1_dph_v_et_p_bins[ p] -> SetBinContent(et+1,width_dph[et][p]*1000.);

			h1_dpp_v_p_et_bins[et] -> SetBinError  (p +1,error_dpp[et][p]*100. );
			h1_dth_v_p_et_bins[et] -> SetBinError  (p +1,error_dth[et][p]*1000.);
			h1_dph_v_p_et_bins[et] -> SetBinError  (p +1,error_dph[et][p]*1000.);
			h1_dpp_v_et_p_bins[ p] -> SetBinError  (et+1,error_dpp[et][p]*100. );
			h1_dth_v_et_p_bins[ p] -> SetBinError  (et+1,error_dth[et][p]*1000.);
			h1_dph_v_et_p_bins[ p] -> SetBinError  (et+1,error_dph[et][p]*1000.);
		}
	}

	// -------------------------------------------------------------
	// Updating table with width values
	ofstream updated_tab;
	if(update_tab){
		updated_tab.open("tables/"+tab_name);
		for(int et = 0 ; et < size_eta_bin-1 ; et++){
			for(int p = 0 ; p < size_mom_bin-1 ; p++){
				updated_tab << width_dpp[et][p];
				if(p == size_mom_bin-2) updated_tab << "\n";
				else updated_tab << "\t";
			}
		}
		updated_tab << "\n";
		for(int et = 0 ; et < size_eta_bin-1 ; et++){
			for(int p = 0 ; p < size_mom_bin-1 ; p++){
				updated_tab << width_dth[et][p];
				if(p == size_mom_bin-2) updated_tab << "\n";
				else updated_tab << "\t";
			}
		}
		updated_tab << "\n";
		for(int et = 0 ; et < size_eta_bin-1 ; et++){
			for(int p = 0 ; p < size_mom_bin-1 ; p++){
				updated_tab << width_dph[et][p];
				if(p == size_mom_bin-2) updated_tab << "\n";
				else updated_tab << "\t";
			}
		}
		updated_tab.close();
	}

	// -------------------------------------------------------------
	// Plotting histograms
	TCanvas * c1 = new TCanvas("c1","c1",1300,900);
	c1 -> Divide(3,2);

	c1 -> cd(1);
	h1_dpp_v_p_et_bins[0] -> Draw();
	for(int et = 0 ; et < size_eta_bin-1 ; et++) h1_dpp_v_p_et_bins[et] -> Draw("same");
	c1 -> cd(2);
	h1_dth_v_p_et_bins[0] -> Draw();
	for(int et = 0 ; et < size_eta_bin-1 ; et++) h1_dth_v_p_et_bins[et] -> Draw("same");
	c1 -> cd(3);
	h1_dph_v_p_et_bins[0] -> Draw();
	for(int et = 0 ; et < size_eta_bin-1 ; et++) h1_dph_v_p_et_bins[et] -> Draw("same");
	c1 -> cd(4);
	h1_dpp_v_et_p_bins[0] -> Draw();
	for(int p = 0 ; p < size_mom_bin-1 ; p++) h1_dpp_v_et_p_bins[p] -> Draw("same");
	c1 -> cd(5);
	h1_dth_v_et_p_bins[0] -> Draw();
	for(int p = 0 ; p < size_mom_bin-1 ; p++) h1_dth_v_et_p_bins[p] -> Draw("same");
	c1 -> cd(6);
	h1_dph_v_et_p_bins[0] -> Draw();
	for(int p = 0 ; p < size_mom_bin-1 ; p++) h1_dph_v_et_p_bins[p] -> Draw("same");
	TLegend * leg1 = new TLegend(0.50,0.6,0.70,0.85);
	leg1 -> SetLineColor(0);
	for(int et = 0 ; et < size_eta_bin-1 ; et++) leg1 -> AddEntry(h1_dph_v_p_et_bins[et],Form("%.1f < |#eta| < %.1f",eta_bin[et],eta_bin[et+1]));
	c1 -> cd(3);
	leg1 -> Draw("same");
	TLegend * leg2 = new TLegend(0.20,0.6,0.40,0.85);
	leg2 -> SetLineColor(0);
	for(int p = 0 ; p < size_mom_bin-1 ; p++) leg2 -> AddEntry(h1_dph_v_et_p_bins[p],Form("%.1f < p < %.1f GeV/c",mom_bin[p],mom_bin[p+1]));
	c1 -> cd(6);
	leg2 -> Draw("same");

	// -------------------------------------------------------------
	// Saving fits to pdf
	TString out_pdf_name = "fits_AllS_"+partic+"_B_"+Form("%.1f",Bfield)+"T.pdf";

	for(int et = 0 ; et < size_eta_bin-1 ; et++){
		TString fname = out_pdf_name;
		if(et == 0) fname+="(";
		else if(et == size_eta_bin-2) fname+=")";
		c_fits_p[et] -> Print(fname);
	}

	// -------------------------------------------------------------
	// Saving histograms
	TFile * Fout = new TFile("histos_"+partic+"_B_"+Form("%.1f",Bfield)+"T.root","recreate");
	for(int et = 0 ; et < size_eta_bin-1 ; et++){
		h1_dpp_v_p_et_bins[et] -> Write(Form("h1_dpp_v_p_et_bins_%i",et));
		h1_dth_v_p_et_bins[et] -> Write(Form("h1_dth_v_p_et_bins_%i",et));
		h1_dph_v_p_et_bins[et] -> Write(Form("h1_dph_v_p_et_bins_%i",et));
	}
	for(int p = 0 ; p < size_mom_bin-1 ; p++){
		h1_dpp_v_et_p_bins[p] -> Write(Form("h1_dpp_v_et_p_bins_%i",p));
		h1_dth_v_et_p_bins[p] -> Write(Form("h1_dth_v_et_p_bins_%i",p));
		h1_dph_v_et_p_bins[p] -> Write(Form("h1_dph_v_et_p_bins_%i",p));
	}
	c1 -> Write("c1");
	Fout -> Close();

}
// ============================================================================================================================================
void prettyTH1F( TH1F * h1 , int color , int marker , float min , float max ){
	h1 -> SetLineWidth(2);
	h1 -> SetLineColor(color);
	h1 -> SetMarkerStyle(marker);
	h1 -> SetMarkerColor(color);

	h1 -> SetMinimum(min);
	h1 -> SetMaximum(max);

	h1 -> GetXaxis() -> CenterTitle();
	h1 -> GetXaxis() -> SetNdivisions(107); // to draw less tick marks
	h1 -> GetYaxis() -> CenterTitle();
	h1 -> GetYaxis() -> SetNdivisions(107); // to draw less tick marks

	h1 -> SetMinimum(0.001);
}

