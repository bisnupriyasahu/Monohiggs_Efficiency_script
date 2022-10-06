void mutau_analyzer::selections(float weight, int shift, string uncObject)
{
  check_unc = false; // set true for printing unc pt, values
  //std::cout<<"coming here 1"<<std::endl;
  double event_weight = weight;
  TLorentzVector metP4, event_met_p4;
  if (shift > 0)
    unc_shift = "up";
  else if (shift < 0)
    unc_shift = "down";
  else
    unc_shift = "nominal";
  shift_index = shift;
  selected_systematic = uncObject;
  // cout<<" selected_systematic = "<< selected_systematic << " shift = "<< shift<<endl;
  std::vector<int> event_mu, event_tau;
  event_mu.clear();
  event_tau.clear();
  muCand.clear();
  tauCand.clear();
  aisrtauCand.clear();
  jetCand.clear();
  for (int i = 0; i < nJet; i++)
    jetPt->at(i) = orginal_jetPt[i];
 
  int mu_index = -1;
  int tau_index = -1;
  int gentau1_index = -1;
  int gentau2_index = -1;
  int genmatchedmu_idx = -1;
  int genmatchedtau_idx = -1;

  if (is_MC)
    event_weight = weight;
  else
    event_weight = 1.0;

  vector<int> gentaucand, genmucand, genhiggscand;
  gentaucand.clear();
  genmucand.clear();
  genhiggscand.clear();

  int pdgid_tocheck = 25;
  if(found_DYjet_sample)
    pdgid_tocheck =23;
  for(int i=0; i<nMC; i++)
  {
   if(abs(mcPID->at(i)) == pdgid_tocheck)
   {
      //std::cout<<"coming inside gen and pdgid_tocheck "<<pdgid_tocheck<<std::endl;
   
      //std::cout<<"coming inside gen higgs"<<pdgid_tocheck<<std::endl;
	    //if(abs(mcDaughterPID->at(i)) == 13){
      // if (abs(mcDaughterPID->at(i)) == 13)
      //cout<< i << " " << mcPID->at(i) << " " << " " << mcDaughterPID->at(i) <<endl;
      //cout<< i <<  " " << mcDaughterPID->at(i) <<endl;
      genhiggscand.push_back(i);
      
    }
  }
  //std::cout<<"coming here 1 "<<std::endl;
  for(int i=0; i<nMC; i++){
    //if(abs(mcMotherPID->at(i))==pdgid_tocheck && mcPID->at(i)==15  )
    if(abs(mcPID->at(i)) == 15  && abs(mcMotherPID->at(i)) == pdgid_tocheck && abs(mcDaughterPID->at(i)) == 13)
      {
	      genmucand.push_back(i);
	      gentau1_index = i;
	      gentau1.SetPtEtaPhiE(mcPt->at(i) , mcEta->at(i) , mcPhi->at(i)  , mcE->at(i));
      }
      //std::cout<<"coming here 2 "<<std::endl;
      if(abs(mcPID->at(i)) == 15 &&  abs(mcMotherPID->at(i)) == pdgid_tocheck && abs(mcDaughterPID->at(i)) != 13 )
      {
	      gentaucand.push_back(i);
	      gentau2_index = i;
	      gentau2.SetPtEtaPhiE(mcPt->at(i) , mcEta->at(i) , mcPhi->at(i)  , mcE->at(i));
      }
    }
    TLorentzVector genmatch_mup4, genmatch_taup4, genmatchboost_mup4, genmatchboost_taup4 ;
    TLorentzVector my_muP4, my_tauP4, my_metP4, myboost_muP4, myboost_tauP4, myboost_metP4;
    double gen_TauPt, gen_subleadingtauPt, gen_taudeltaR, gen_tauhiggsPt, gen2_TauPt, gen2_subleadingtauPt, gen2_taudeltaR, gen2_tauhiggsPt;
    metP4.SetPtEtaPhiE(pfMET, 0, pfMETPhi, pfMET);
    event_met_p4 = metP4;
    if(gentau1_index >= 0 and gentau2_index >= 0)
    {
      //std::cout<<"coming here 3 "<<std::endl;
      double genTauPt = max( gentau1.Pt(), gentau2.Pt() );
      double gensubleadingtauPt = min( gentau1.Pt(), gentau2.Pt() );
      double gentau_deltaR = gentau1.DeltaR(gentau2);
      double gentau_higgsPt = (gentau1 + gentau2).Pt();
      
      
      //if( gentau1.Pt()>20 && gentau2.Pt()>20  && abs(gentau1.Eta())<2.3 && abs(gentau2.Eta())<2.3 )
      if(abs(gentau1.Eta())<2.3 && abs(gentau2.Eta())<2.3 )
	    {
        
	      //std::cout<<"gen index "<<gentau1_index<<"2nd index : "<<gentau2_index<<std::endl;
	      plotFill("gentauPt_raw_0", genTauPt, 970, 30, 1000, event_weight);
	      plotFill("gensubleadingtauPt_raw_0", gensubleadingtauPt, 970, 30, 1000, event_weight);
	      plotFill("genHiggsPt_raw_0", gentau_higgsPt, 970, 30, 1000, event_weight);
	      plotFill("gendeltaR_raw_0", gentau_deltaR, 600, 0, 6, event_weight);
	      if(gentau1.Pt()>20 && gentau2.Pt()>20)
		    {
		  
		      if (nMu>0 && nTau > 0)
		      {
		        pair<int, int> selected_indices = get_index();
		        make_iso_plot = false;
		        mu_index = selected_indices.first;
		        tau_index = selected_indices.second;
		        if (mu_index >= 0 && tau_index >= 0)
		        {
		          //std::cout<<"mu index "<<mu_index<<"2nd index : "<<tau_index<<std::endl;
              //MuIndex = mu_index;
		          //TauIndex = tau_index;
		          my_muP4.SetPtEtaPhiE(muPt->at(mu_index), muEta->at(mu_index),muPhi->at(mu_index), muE->at(mu_index));
		          my_tauP4.SetPtEtaPhiE(tau_Pt->at(tau_index), tau_Eta->at(tau_index), tau_Phi->at(tau_index), tau_Energy->at(tau_index));
		          //my_metP4.SetPtEtaPhiE(pfMET ,0,pfMETPhi,pfMET);
		          //genmatching loop:
		          if(my_muP4.DeltaR(gentau1)<0.1 && my_tauP4.DeltaR(gentau2)<0.1)
		        	{
		        	  genmatchedmu_idx = gentau1_index;
		        	  genmatchedtau_idx = gentau2_index;
		        	}
		          else if(my_muP4.DeltaR(gentau2)<0.1 && my_tauP4.DeltaR(gentau1)<0.1)
		        	{
		        	  genmatchedmu_idx = gentau2_index;
		        	  genmatchedtau_idx = gentau1_index;
		        	}
		          if(genmatchedmu_idx >= 0 && genmatchedtau_idx >= 0)
              {
			          std::cout<<"coming here 1"<<std::endl;
			          genmatch_mup4.SetPtEtaPhiE(mcPt->at(genmatchedmu_idx) , mcEta->at(genmatchedmu_idx) , mcPhi->at(genmatchedmu_idx)  , mcE->at(genmatchedmu_idx));
			          genmatch_taup4.SetPtEtaPhiE(mcPt->at(genmatchedtau_idx) , mcEta->at(genmatchedtau_idx) , mcPhi->at(genmatchedtau_idx)  , mcE->at(genmatchedtau_idx));
			          gen_TauPt = genmatch_mup4.Pt();
			          gen_subleadingtauPt = genmatch_taup4.Pt();
			          gen_taudeltaR = genmatch_mup4.DeltaR(genmatch_taup4);
			          gen_tauhiggsPt = (genmatch_mup4 + genmatch_taup4).Pt();
			          cout<<__LINE__<<endl; // 

			          plotFill("gentauPt_num_1", gen_TauPt, 970, 30, 1000, event_weight);
			          plotFill("gensubleadingtauPt_num_1", gen_subleadingtauPt, 970, 30, 1000, event_weight);
			          plotFill("genHiggsPt_num_1", gen_tauhiggsPt, 970, 30, 1000, event_weight);
			          plotFill("gendeltaR_num_1", gen_taudeltaR, 600, 0, 6, event_weight);
			          cout<<__LINE__<<endl; // 			    
			          jetCand.clear();
			          jetCand = getJetCand(mu_index, tau_index);
			          my_njets = jetCand.size();
			          plot_resolved_taus(genmatch_mup4, genmatch_taup4, tau_index, "4", event_weight);
			        }
			      }
		      }
		    } 
      }
    plot_boosted = false;
    mu_index = -1;
    tau_index = -1;
    genmatchedmu_idx = -1;
    genmatchedtau_idx = -1;
    if( abs(gentau1.Eta())<2.3 && abs(gentau2.Eta())<2.3 && gentau1.Pt()>20 && gentau2.Pt()>20)
	  {
	    cout<<__LINE__<<endl; // 
	    plotFill("gentauPt_boostedraw_1", genTauPt, 970, 30, 1000, event_weight);
	    plotFill("gensubleadingtauPt_boostedraw_1", gensubleadingtauPt, 970, 30, 1000, event_weight);
	    plotFill("genHiggsPt_boostedraw_1", gentau_higgsPt, 970, 30, 1000, event_weight);
	    plotFill("gendeltaR_boostedraw_1", gentau_deltaR, 600, 0, 6, event_weight);
	    cout<<__LINE__<<endl; // 	  
	    if (nMu> 0 && nBoostedTau>0 )
	    {
	      pair<int, int> selected_indices = get_index_2();
	      make_iso_plot = false;
	      mu_index = selected_indices.first;
	      tau_index = selected_indices.second;
	      //Reco Mu and tau
	      cout<<__LINE__<<endl; // 
	      if (mu_index >= 0 && tau_index >= 0)
	      {
	        //MuIndex = mu_index;
	        //TauIndex = tau_index;
	        myboost_muP4.SetPtEtaPhiE(muPt->at(mu_index), muEta->at(mu_index),muPhi->at(mu_index), muE->at(mu_index));
	        myboost_tauP4.SetPtEtaPhiE(boostedTauPt->at(tau_index), boostedTauEta->at(tau_index), boostedTauPhi->at(tau_index), boostedTauEnergy->at(tau_index));
	        //myboost_metP4.SetPtEtaPhiE(pfMET ,0,pfMETPhi,pfMET);
	        //gen matching loop
	        cout<<__LINE__<<endl; // 	
	        if(myboost_muP4.DeltaR(gentau1)<0.1 && myboost_tauP4.DeltaR(gentau2)<0.1)
	          {
	            cout<<__LINE__<<endl; //
	            genmatchedmu_idx = gentau1_index;
	            genmatchedtau_idx = gentau2_index;
	          }
	        else if(myboost_muP4.DeltaR(gentau2)<0.1 && myboost_tauP4.DeltaR(gentau1)<0.1)
	          {
	            cout<<__LINE__<<endl; //
	            genmatchedmu_idx = gentau2_index;
	            genmatchedtau_idx = gentau1_index;
	          }
	        if (genmatchedmu_idx >= 0 && genmatchedtau_idx >=0)
	        {
	          cout<<__LINE__<<endl; //
	          genmatchboost_mup4.SetPtEtaPhiE(mcPt->at(genmatchedmu_idx) , mcEta->at(genmatchedmu_idx) , mcPhi->at(genmatchedmu_idx)  , mcE->at(genmatchedmu_idx));
	          genmatchboost_taup4.SetPtEtaPhiE(mcPt->at(genmatchedtau_idx) , mcEta->at(genmatchedtau_idx) , mcPhi->at(genmatchedtau_idx)  , mcE->at(genmatchedtau_idx));
	          gen2_TauPt = genmatchboost_mup4.Pt();
	          gen2_subleadingtauPt = genmatchboost_taup4.Pt();
	          gen2_taudeltaR = genmatchboost_mup4.DeltaR(genmatchboost_taup4);
	          gen2_tauhiggsPt = (genmatchboost_mup4 + genmatchboost_taup4).Pt();
	          cout<<__LINE__<<endl; //
	          plotFill("gentauPt_boostedraw_2", gen2_TauPt, 970, 30, 1000, event_weight);
            cout<<__LINE__<<endl; //
	          plotFill("gensubleadingtauPt_boostedraw_2", gen2_subleadingtauPt, 970, 30, 1000, event_weight);
	          cout<<__LINE__<<endl; //
            plotFill("genHiggsPt_boostedraw_2", gen2_tauhiggsPt, 970, 30, 1000, event_weight);
	          cout<<__LINE__<<endl; //
            plotFill("gendeltaR_boostedraw_2", gen2_taudeltaR, 600, 0, 6, event_weight);
            cout<<__LINE__<<endl; //
	          plot_boosted_taus(genmatchboost_mup4, genmatchboost_taup4, tau_index , "4", event_weight);
		        cout<<__LINE__<<endl; //
		      }
          cout<<__LINE__<<endl; //
	      }	  
         cout<<__LINE__<<endl; //
	    }
       cout<<__LINE__<<endl; //
    }
     cout<<__LINE__<<endl; //
  }
   cout<<__LINE__<<endl; //
}
void mutau_analyzer::plot_resolved_taus(TLorentzVector muP4, TLorentzVector tauP4, int tauindex, string hnumber, double event_weight)
{
cout<<__LINE__<<endl; //  
  make_plot("raw", hnumber, muP4, tauP4, tauindex, event_weight);
cout<<__LINE__<<endl; //
  if (tau_byVVVLooseDeepTau2017v2p1VSjet->at(tauindex) == 1)
  {
    make_plot("deepVVVLoose", hnumber, muP4, tauP4, tauindex, event_weight);
  }
  cout<<__LINE__<<endl; //
  if (tau_byVVLooseDeepTau2017v2p1VSjet->at(tauindex) == 1)
  {
    make_plot("deepVVLoose", hnumber, muP4, tauP4, tauindex, event_weight);
  }
  cout<<__LINE__<<endl; //
  if (tau_byVLooseDeepTau2017v2p1VSjet->at(tauindex) == 1)
  {
    make_plot("deepVLoose", hnumber, muP4, tauP4, tauindex, event_weight);
  }
  cout<<__LINE__<<endl;

  if (tau_byLooseDeepTau2017v2p1VSjet->at(tauindex) == 1)
  {
    make_plot("deepLoose", hnumber, muP4, tauP4, tauindex, event_weight);
  }
  cout<<__LINE__<<endl;
  if (tau_byMediumDeepTau2017v2p1VSjet->at(tauindex) == 1)
  {
    make_plot("deepMedium", hnumber, muP4, tauP4, tauindex, event_weight);
  }
  cout<<__LINE__<<endl;
  if (tau_byTightDeepTau2017v2p1VSjet->at(tauindex) == 1)
  {
    make_plot("deepTight", hnumber, muP4, tauP4, tauindex, event_weight);
  }
  cout<<__LINE__<<endl;
  if (tau_byVTightDeepTau2017v2p1VSjet->at(tauindex) == 1)
  {
    make_plot("deepVTight", hnumber, muP4, tauP4, tauindex, event_weight);
  }
  cout<<__LINE__<<endl;
  if (tau_byVVTightDeepTau2017v2p1VSjet->at(tauindex) == 1)
  {
    make_plot("deepVVTight", hnumber, muP4, tauP4, tauindex, event_weight);
  }
  cout<<__LINE__<<endl;

  if (tau_IDbits->at(tauindex) >> 12 & 1 == 1)
  {
    make_plot("2017VVLoose", hnumber, muP4, tauP4, tauindex, event_weight);
  }
  cout<<__LINE__<<endl;
  if (tau_IDbits->at(tauindex) >> 13 & 1 == 1)
  {
    make_plot("2017VLoose", hnumber, muP4, tauP4, tauindex, event_weight);
  }
  cout<<__LINE__<<endl;
  if (tau_IDbits->at(tauindex) >> 14 & 1 == 1)
  {
    make_plot("2017Loose", hnumber, muP4, tauP4, tauindex, event_weight);
  }
  cout<<__LINE__<<endl;
  if (tau_IDbits->at(tauindex) >> 15 & 1 == 1)
  {
    make_plot("2017Medium", hnumber, muP4, tauP4, tauindex, event_weight);
  }
  cout<<__LINE__<<endl;
  if (tau_IDbits->at(tauindex) >> 16 & 1 == 1)
  {
    make_plot("2017Tight", hnumber, muP4, tauP4, tauindex, event_weight);
  }
  cout<<__LINE__<<endl;
  if (tau_IDbits->at(tauindex) >> 17 & 1 == 1)
  {
    make_plot("2017VTight", hnumber, muP4, tauP4, tauindex, event_weight);
  }
  cout<<__LINE__<<endl;
  if (tau_IDbits->at(tauindex) >> 18 & 1 == 1)
  {
    make_plot("2017VVTight", hnumber, muP4, tauP4, tauindex, event_weight);
  }
  if (tau_IDbits->at(tauindex) >> 20 & 1 == 1)
  {
    make_plot("2016VVLoose", hnumber, muP4, tauP4, tauindex, event_weight);
  }
  if (tau_IDbits->at(tauindex) >> 21 & 1 == 1)
  {
    make_plot("2016VLoose", hnumber, muP4, tauP4, tauindex, event_weight);
  }

  if (tau_IDbits->at(tauindex) >> 22 & 1 == 1)
  {
    make_plot("2016Loose", hnumber, muP4, tauP4, tauindex, event_weight);;
  }
  if (tau_IDbits->at(tauindex) >> 23 & 1 == 1)
  {
    make_plot("2016Medium", hnumber, muP4, tauP4, tauindex, event_weight);
  }
  if (tau_IDbits->at(tauindex) >> 24 & 1 == 1)
  {
    make_plot("2016Tight", hnumber, muP4, tauP4, tauindex, event_weight);
  }
  if (tau_IDbits->at(tauindex) >> 25 & 1 == 1)
  {
    make_plot("2016VTight", hnumber, muP4, tauP4, tauindex, event_weight);
  }
  if (tau_IDbits->at(tauindex) >> 26 & 1 == 1)
  {
    make_plot("2016VVTight", hnumber, muP4, tauP4, tauindex, event_weight);
  }

}

void mutau_analyzer::plot_boosted_taus(TLorentzVector muP4, TLorentzVector tauP4, int tauindex, string hnumber, double event_weight)
{
cout<<__LINE__<<endl; //
  make_plot("boostedraw", hnumber, muP4, tauP4, tauindex, event_weight);
  // cout<<__LINE__<<endl;
  if (boostedTauByVLooseIsolationMVArun2v1DBoldDMwLTNew->at(tauindex) == 1)
  {
    make_plot("boostedVLoose", hnumber, muP4, tauP4, tauindex, event_weight);
  }
  // cout<<__LINE__<<endl;
  if (boostedTauByLooseIsolationMVArun2v1DBoldDMwLTNew->at(tauindex) == 1)
  {
    make_plot("boostedLoose", hnumber, muP4, tauP4, tauindex, event_weight);
  }
  // cout<<__LINE__<<endl;
  if (boostedTauByMediumIsolationMVArun2v1DBoldDMwLTNew->at(tauindex) == 1)
  {
    make_plot("boostedMedium", hnumber, muP4, tauP4, tauindex, event_weight);
  }
  if (boostedTauByTightIsolationMVArun2v1DBoldDMwLTNew->at(tauindex) == 1)
  {
    make_plot("boostedTight", hnumber, muP4, tauP4, tauindex, event_weight);
  }
  if (boostedTauByVTightIsolationMVArun2v1DBoldDMwLTNew->at(tauindex) == 1)
  {
    make_plot("boostedVTight", hnumber, muP4, tauP4, tauindex, event_weight);
  }
}

pair<int, int> mutau_analyzer::get_index()
{
cout<<__LINE__<<endl; //
  int mu_index = -1;
  int tau_index = -1;
  float relMuIso, relMuIso_old, relMuIso_v3;

  if (nMu > 0 && nTau>0)
  {
    for (int i = 0; i < nMu; i++)
	  {
	  for (int j = 0; j < nTau; j++)
	    {
	      int iMu = i;
	      int iTau = j;
	      TLorentzVector mu_p4, tau_p4;
	      mu_p4.SetPtEtaPhiE(muPt->at(iMu), muEta->at(iMu), muPhi->at(iMu), muE->at(iMu));
	      tau_p4.SetPtEtaPhiE(tau_Pt->at(iTau), tau_Eta->at(iTau), tau_Phi->at(iTau), tau_Energy->at(iTau));
	      float mu_tau_dr = mu_p4.DeltaR(tau_p4);
	      double higgspt = (mu_p4 + tau_p4).Pt();
	      relMuIso = (muPFChIso->at(iMu) + max(muPFNeuIso->at(iMu) + muPFPhoIso->at(iMu) - 0.5 * muPFPUIso->at(iMu), 0.0) - tau_Pt->at(iTau)) / (muPt->at(iMu));
	      bool pass_bjet_veto = ((bJet_medium(iMu, iTau).size() == 0) && (bJet_loose(iMu, iTau).size() < 2));
	      bool pass3rdLeptonVeto = (passDiMuonVeto(iMu) == true && eVetoZTTp001dxyz(iMu, iTau) && mVetoZTTp001dxyz(iMu, iTau));
	      
	      double dr_aa = mu_p4.DeltaR(gentau1);
	      double dr_ab = mu_p4.DeltaR(gentau2);
	      double dr_ba = tau_p4.DeltaR(gentau1);
	      double dr_bb = tau_p4.DeltaR(gentau2);


	      if (  muPt->at(iMu) > 20 
		    && fabs(muEta->at(iMu)) < 2.4
		    && fabs(muDz->at(iMu)) < 0.2 
		    && fabs(muD0->at(iMu)) < 0.045
		    && muCharge->at(iMu) * tau_Charge->at(iTau) < 0 
		    //&& muIDbit->at(iMu) >> 1 & 1 == 1 // muon id
		    && tau_Pt->at(iTau) > 30 
		    && fabs(tau_Eta->at(iTau)) < 2.3
		    && tau_LeadChargedHadron_dz->at(iTau) < 0.2
		    //&& relMuIso < 0.15               // muon relatice isolation
		    //&& tau_byMediumDeepTau2017v2p1VSjet->at(iTau)==1
		    //&& tau_IDbits->at(iTau)>>1&1==1
		    && ( tau_DecayMode->at(iTau)>=0)
		    //&& ( (dr_aa < 0.1 && dr_bb < 0.1) ||  ( dr_ab < 0.1 &&  dr_ba < 0.1) )
		    //&& ( dr_aa < 0.1 || dr_ab < 0.1 || dr_ba < 0.1 || dr_bb < 0.1 )
		    )
  		  {
  		    return make_pair(iMu, iTau);
  		    /////////// okay we found the mu-au pair, exit these loops
  		  }
	    }
	  }
  }
  return make_pair(-1, -1);
}

pair<int, int> mutau_analyzer::get_index_2()
{
cout<<__LINE__<<endl; //
  int mu_index = -1;
  int tau_index = -1;
  
  if (nMu > 0 && nTau>0)
  {
    for (int i = 0; i < nMu; i++)
	  {
	    for (int j = 0; j < nBoostedTau; j++)
	    {
	      int iMu = i;
	      int iTau = j;
	      TLorentzVector mu_p4, tau_p4;
	      mu_p4.SetPtEtaPhiE(muPt->at(iMu), muEta->at(iMu), muPhi->at(iMu), muE->at(iMu));
	      tau_p4.SetPtEtaPhiE(boostedTauPt->at(iTau), boostedTauEta->at(iTau), boostedTauPhi->at(iTau), boostedTauEnergy->at(iTau));
	      
	      float mu_tau_dr = mu_p4.DeltaR(tau_p4);
	      double higgspt = (mu_p4 + tau_p4).Pt();

	      double dr_aa = mu_p4.DeltaR(gentau1);
	      double dr_ab = mu_p4.DeltaR(gentau2);
	      double dr_ba = tau_p4.DeltaR(gentau1);
	      double dr_bb = tau_p4.DeltaR(gentau2);

	      if ( muPt->at(iMu) > 20 
		       && fabs(muEta->at(iMu)) < 2.4
		       && fabs(muDz->at(iMu)) < 0.2 
		       && fabs(muD0->at(iMu)) < 0.045
		       && muCharge->at(iMu) * boostedTauCharge->at(iTau) < 0 
		       //&& muIDbit->at(iMu) >> 1 & 1 == 1 // muon id
		       && tau_p4.Pt() > 30
		       && fabs(tau_p4.Eta()) < 2.3
		       && boostedTauZImpact->at(iTau) < 0.2
		       //&& muIDbit->at(iMu) >> 1 & 1 == 1 // muon id
		       && boostedTaupfTausDiscriminationByDecayModeFindingNewDMs->at(iTau) > 0.5 
		       //&& ( (dr_aa < 0.1 && dr_bb < 0.1) ||  ( dr_ab < 0.1 &&  dr_ba < 0.1) )
		       //&& ( dr_aa < 0.1 || dr_ab < 0.1 || dr_ba < 0.1 || dr_bb < 0.1 )
		     )
		    {
		      return make_pair(iMu, iTau);
		    }
	    }
	  }
  }
  return make_pair(-1, -1);
}
void mutau_analyzer::make_plot(string idtype, string hnumber, TLorentzVector muP4, TLorentzVector tauP4, int tauindex, double event_weight)
{
  double MuPt = muP4.Pt();
  double subleadingtauPt = tauP4.Pt();
  double deltaR = muP4.DeltaR(tauP4);
  double higgsPt = (muP4 + tauP4).Pt();
  cout<<__LINE__<<endl; //
  plotFill("muPt_"+ idtype + "_" + hnumber, MuPt, 970, 30, 1000, event_weight);
  plotFill("tauPt_"+ idtype + "_" + hnumber, subleadingtauPt , 970, 30, 1000, event_weight);
  plotFill("deltaR_"+ idtype + "_" + hnumber, deltaR, 600, 0, 6, event_weight);
  plotFill("HiggsPt_"+ idtype + "_" + hnumber, higgsPt, 970, 30, 1000, event_weight);
  /*
    plotFill("genmuPt_"+ idtype + "_" + hnumber, genmuPt, 970, 30, 1000, event_weight);
    plotFill("gentauPt_"+ idtype + "_" + hnumber, gentauPt, 970, 30, 1000, event_weight);
    plotFill("genHiggsPt_"+ idtype + "_" + hnumber, gentau_higgsPt, 970, 30, 1000, event_weight);
    plotFill("gendeltaR_"+ idtype + "_" + hnumber, gentau_deltaR, 600, 0, 6, event_weight);
  */  
}
