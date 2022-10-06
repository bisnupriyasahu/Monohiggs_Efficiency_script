import ROOT as R
from ROOT import TCanvas, TGraph
from ROOT import gROOT
from math import sin
from array import array
from itertools import product

R.gStyle.SetFrameLineWidth(1)
R.gStyle.SetLineWidth(2)
R.gStyle.SetOptStat(0)
R.gStyle.SetEndErrorSize(5)

def make_legend():
  output = R.TLegend(0.85, 0.2, 0.99, 0.50, "", "brNDC")
  #output = ROOT.TLegend(0.2, 0.1, 0.47, 0.65, "", "brNDC")
  output.SetLineWidth(1)
  output.SetLineStyle(1)
  output.SetFillStyle(1001) #0
  output.SetFillColor(0)
  output.SetBorderSize(1)
  output.SetTextFont(42)
  output.SetNColumns(1)                                                                                                                                                                                                                                              
                                       
  return output

def add_lumi():
	lowX=0.550
	lowY=0.82
	lumi  = R.TPaveText(lowX, lowY+0.06, lowX+0.40, lowY+0.15, "NDC")
	lumi.SetBorderSize(   0 )
	lumi.SetFillStyle(    0 )
	lumi.SetTextAlign(   32 )#12
	lumi.SetTextColor(    1 )
	lumi.SetTextSize(0.03)
	lumi.SetTextFont (   42 )
	lumiProcessed="41.52"
	channel_ = 'mutau'
	if channel_=="combined":
		lumi.AddText("4 channels combined 2017, "+lumiProcessed+" fb^{-1} (13 TeV)")
	if channel_=="mutau":
		lumi.AddText("#mu#tau_{h} 2017, "+lumiProcessed+" fb^{-1} (13 TeV)")
	if channel_=="etau":
		lumi.AddText("e#tau_{h} 2017, "+lumiProcessed+" fb^{-1} (13 TeV)")
	if channel_=="tautau":
		lumi.AddText("#tau_{h}#tau_{h} 2017, "+lumiProcessed+" fb^{-1} (13 TeV)")
	if channel_=="emu":
		lumi.AddText("e#mu 2017, "+lumiProcessed+" fb^{-1} (13 TeV)")
	return lumi

#mzp_map = {'2_0310_v10' : 'MZp_1000_MChi_1' , '7_0310_v10': 'MZp_100_MChi_1', '10_0310_v10': 'MZp_1500_MChi_1',  'DY': 'DYJetsToLL'}
mzp_map = {'2_0610_v3' : 'MZp_1000_MChi_1' }
canv=R.TCanvas("canvas","",0,0,1300,1200)
canv.cd()

vars = {'HiggsPt', 'muPt', 'deltaR', 'tauPt'}
#selections = ['boostedraw', 'num']
#idtypes = { 'deep' : ['VVVLoose', 'Loose', 'Medium', 'Tight', 'VVTight'] ,
#'boosted' : ['VLoose', 'Loose', 'Medium', 'Tight', 'VTight'] ,
#            '2017'    : ['VLoose' , 'Loose', 'Medium', 'Tight', 'VTight'] ,
#           '2016'    : ['VLoose' , 'Loose', 'Medium', 'Tight', 'VTight']
#  }
#new_binning = array('d', [30, 50, 70, 90, 110, 130, 150, 170, 190, 210, 230, 250, 270, 300, 350, 400, 450,  500, 600, 700,800,900, 1000])


#parameters = list(product(vars, selections, idtypes))
#parameters = list(product(vars, idtypes))

for idx in mzp_map:
	mzp_idx = idx
	file_name = 'output_0610/Zpbaryonic2017_'+mzp_idx+'.root'
	if idx=='DY':
		file_name = 'output_0610/DYJetsToLL_M-50_TuneCP5_0610_v3.root'
	infile = R.TFile(file_name, 'update')

	for var in vars:
		signal = mzp_map[mzp_idx]
                
		raw_name = 'gen'+var+'_raw_0'
 
   
                print raw_name
		hist_raw = infile.Get(raw_name)
                #		idlist = idtypes[idtype]
                #case 1
		#hist_1 = infile.Get('gen'+var+'_num_1')
		hist_2 = infile.Get('gen'+var+'_boostedraw_2')
		
                #tmpHist_num1.Divide(tmpHist_den)
                #tmpHist_num2.Divide(tmpHist_den)

		if var != 'deltaR':
			#tmpHist_num1 = hist_1.Rebin(22, 'gen'+var+'_num_1', new_binning)
			tmpHist_num2 = hist_2.Rebin(22, 'gen'+var+'_boostedraw_2')
                        tmpHist_den = hist_raw.Rebin(22, raw_name )  
		else:
			dr_bins = 5
			#tmpHist_num1 = hist_1.Rebin(dr_bins, 'gen'+var+'_num_1')
			tmpHist_num2 = hist_2.Rebin(dr_bins, 'gen'+var+'_boostedraw_2')
                        tmpHist_den = hist_raw.Rebin(dr_bins, raw_name)  		

		
                #tmpHist_num1.Divide(tmpHist_den)
                tmpHist_num2.Divide(tmpHist_den)
                #print"num1 {}/ denum {}".format(tmpHist_num1.GetName(),tmpHist_den.GetName())
                print"num2 {}/ denum {}".format(tmpHist_num2.GetName(),tmpHist_den.GetName())

                if var=='deltaR' :		
                  tmpHist_num2.GetXaxis().SetRangeUser(0, 4.0) 

		#if var=='deltaR' :		
                # tmpHist_num1.GetXaxis().SetRangeUser(0, 4.0) 


		tmpHist_num2.GetXaxis().SetTitle(var)

		# if var=='tauPt':
		# 	var = 'muonPt'
		# elif var=='subleadingtauPt':
		# 	var = 'tauPt'
		tmpHist_num2.SetTitle(var+ ' efficiency case1 booosted tau isolation '+signal)
		tmpHist_num2.SetMaximum(1.2)
		tmpHist_num2.SetMinimum(0.0)

		#tmpHist_num2.SetMarkerStyle(20)
		tmpHist_num2.SetMarkerStyle(21)
		
		#tmpHist_num2.SetMarkerSize(3)
		tmpHist_num2.SetMarkerSize(3)

		#tmpHist_num2.SetMarkerColor(2)
		tmpHist_num2.SetMarkerColor(3)
		
		#tmpHist_num2.SetLineColor(2)
		tmpHist_num2.SetLineColor(3)
		

		nbinsX = tmpHist_num2.GetNbinsX()

		tmpHist_num2.Draw("e1")
		#tmpHist_num2.Draw("e1same")
		

		legende=make_legend()
		legende.AddEntry(tmpHist_num2,"boosted eficiency","lp")
		#legende.AddEntry(tmpHist_num2,"boosted efficiency","lp")

		legende.Draw("same")
		l1=add_lumi()
		l1.Draw("same")

		canv.SaveAs('plots_case1/eff_2017mutau_Boosted'+signal+'_'+var+'_'+'.png')

