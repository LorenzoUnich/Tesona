#!/usr/bin/env python
#the same as the other, just a copy to do another run 

"""
arca point source analysis 

Run the pseudo-experiments for candiates nr 1-10
 arca_binned -c 1...10 

Combine all the availalbe pseudo-experiment outputs in one BinnedPointSourceAnalysis,
writing pickle file with the combined analysis.
 araca_binned -a

Produce results, from a combined analysis:
 arca_binnned -s (or -u for unblined)

 options:
 -a : agragate/combine the candidates searches       default=False
 -A : ANTARES sample                                 default=No
 -b : batch                                          default=True
 -B : band for the noise modelling                   default="all"
 -c : run a list of candidates                       default=none
 -d : debug                                          default=False
 -F : fluxexpression                                 default=1
 -W : fluxtype                                       default=Expr
 -D : differential limits                            default=False
 -M : minimum energy                                 default=1
 -m : maximum energy                                 default=8
 -f : ARCA21 period                                  default=ARCA21_first
 -h : this help message and exit                     default=False       
 -i : Python interative mode (prompt when done)      default=False  
 -n : output name                                    default=Source
 -p : output path                                    default=standard
 -q : only first period                              default=True
 -s : produce results                                default=True
 -t : path of the root of this project               default=../
 -L : make Likelihood curve                          default=True
 -T : make Data scramble                             default=False
 -w : write webpage                                  default=False
"""

import aa, ROOT, os, sys, pickle, copy, glob
from math import *
#from array import array
import array
import libROOTPythonizations
sys.modules['libROOTPythonizations3_7'] = libROOTPythonizations


from ana.search import BinnedPointSourceAnalysis, SearchPeriod, CandidateList, DataSample, BinnedPointSourceSearch

from util.webfile import JsROOTFile
from background_model import ArcaBackgroundModel
from data_sample import ArcaDataSample

# channels
from ROOT.defs import track, shower
from ROOT import TApplication
# mcfiletypes
from ROOT.defs import anutauCCshowerdecay, anutauCCmuondecay, anutauNC, anumuCC, anumuNC, anueCC, anueNC, anueGLRES
from ROOT.defs import muonMuon, nueCC, nueNC, numuCC, numuNC, nutauCCshowerdecay, nutauCCmuondecay, nutauNC
from iminuit import Minuit
import multiprocessing
import pandas as pd
#multiprocessing.set_start_method('spawn')

import time
import math

L = "anutauCCshowerdecay anutauCCmuondecay anutauNC anumuCC anumuNC anueCC anueNC muonMuon nueCC nueNC numuCC numuNC nutauCCshowerdecay nutauCCmuondecay nutauNC"

allmcfiletypes = [ getattr(ROOT.defs, x) for x in L.split() ]

options = aa.Options( __doc__, sys.argv[1:])

ROOT.stringutil.verbosity( ROOT.stringutil.debug )

ROOT.gSystem.Load("${ROOTSYS}/lib/libMinuit.so")
ROOT.gSystem.Load("${ROOTSYS}/lib/libMinuit2.so")
ROOT.gStyle.SetOptStat(0)
ROOT.TH1.AddDirectory( False ); 
ROOT.gROOT.SetBatch( options.b )

if options.p != 'standard':
    outputpath = options.p
else:
    outputpath = options.t + "/ana_pickles"
os.system(f"mkdir -p {outputpath}")
#print(ROOT.Det())

# The following trick is needed to read old (existing) picklefiles.
# It loads the streamerinfo needed to read them.# # next time, let's store the evt's in a root file.

dum = ROOT.TFile("./N2022_PS_selection/datav6.2.jchain.aashower.dst.merged_9635_10005_pre_upm01_antinoise_upaam01.root")

#Differential sensitivity

energy_min = options.M
energy_max = options.m

#@staticmethod
catalouge_name = "catalogue_MeerKAT_final.csv"

def get_gauss_expr_string(params):
    """
    Gives as a result the string of the gaussian for the input parameters.
    params: list of 9 paramters [norm1, mean1, sigma1, norm2, mean2, sigma2, norm3, mean3, sigma3]
    """
    return   f"{params[0]}*exp(-0.5*((x - {params[1]})/{params[2]})**2) + " f"{params[3]}*exp(-0.5*((x - {params[4]})/{params[5]})**2) + " f"{params[6]}*exp(-0.5*((x - {params[7]})/{params[8]})**2)"
    import numpy as np
from math import exp, sqrt, pi

def triple_gaussian(x, norm1, mean1, sigma1,
                         norm2, mean2, sigma2,
                         norm3, mean3, sigma3):
    def gauss(n, m, s):
        return n * exp(-0.5 * ((x - m)/s)**2) / (s * sqrt(2*pi))
    return gauss(norm1, mean1, sigma1) + gauss(norm2, mean2, sigma2) + gauss(norm3, mean3, sigma3)
def load( picklefile, s_min, s_max ) :
    'workaround for loading old pickle files'
    import __main__
    __main__.DataSample = DataSample
    datasample = pickle.load( open(picklefile,"rb") )  
    datasample.__class__ = ArcaDataSample # you happen to pull this shit while I'm in a transitional phase
    datasample.dir_z_cut = -0.1
    datasample.compute_coords( scramble = True ) 
    datasample.fill_hist(s_min, s_max)
    return datasample

def weight_source(x):
    #x =sys.argv[1]
    
    with open("new_S_radio_weights.txt", "r") as file:
        cont = file.read()
        values = cont.split()
        weight = float(values[x])

    return weight


def my_graph_flux():

        f = open("./file_modelling_Seyferts/file_flux_"+str(options.c)+".txt")
        loge_s, flux_s = array.array('d'), array.array('d')

        for line in f:
            if line.startswith("#") : continue
            energy, flux = map( float, line.split() )
            loge_s.append(log10(energy))
            flux_s.append(flux)

        graph_flux = ROOT.TGraph( len(loge_s), loge_s, flux_s)

        return graph_flux


class ArcaBinnedPointSourceAnalysis ( BinnedPointSourceAnalysis ):

    """
    This is where we define all things that are specific to this analysis:
    datasets, irfs, etc.
    
    """
    def __init__( self , name ) :
        
        self.name = name
        BinnedPointSourceAnalysis.__init__(self)

        self.candlist = CandidateList(options.t + '/inputdata/' + catalouge_name )
        print(self.candlist.pretty_table())

        self.picklepath  = outputpath
        #self.fluxexpr    = options.F 
 
        if options.D==True:
 
            self.fluxexpr    = "("+str(weight_source(int(options.c)))+"2*0.5*1e-4*(x>10**"+str(energy_min)+" && x<10**"+str(energy_max)+")*x**-2)" #there was a 5 as a normalization. now there is a two 
            #self.fluxexpr    = "(2e4*0.5*1e-4*(x>10**"+str(energy_min)+" && x<10**"+str(energy_max)+")*x**-2)"
            
        elif options.D==False: #PER ORA CONSIDERA QUESTO CAPRA!

           #self.fluxexpr    = "("+str(weight_source(int(options.c)))+"*1e-4*x**-2)" #Seyfert spectrum -2 weighted  Walid non è sicuro che fatto re può mettere.
           #hai messo il fattore per due rispetto all'iniziale
           #non è detta che con il peso che tu metti che 6 vada bene (potrebbe essere 10!!) allora metti il pes e poi metti un fattore avanti a tutte elsorgenti.
           #quindi metti un fattoe che rinormalizzi
           #stackling ha flag L e quella la lancio per tutto il catalogo
            #self.fluxexpr    = "("+str(weight_source(int(options.c)))+"*1e4*0.5*(1e4)*(1.5332e-11)*((x/100)**-2.8))" #SBG weighted
            #self.fluxexpr="(0.5*(1e4)*(1.5332e-11)*((x/100)**-2.0))" # SBG weighted equal
           self.fluxexpr="("+str(weight_source(int(options.c)))+"*1e-4*x**-2)" #CAMBIATO IL FLUSSO 
#            self.fluxexpr="10*0.5*0.5*1e-4*x**-2"

        self.flag_flux= options.W

        if self.flag_flux == "Expr":
            print("Power-law expression")
            #self.fluxexpr    = "("+str(weight_source(int(options.c)))+"*0.5*1e-4*x**-2)"
        elif self.flag_flux == "Graph":
            print("Graph flux")
            self.fluxgraph = my_graph_flux()

        det_ARCA = ROOT.Det()
        det_ARCA.longitude = 0.278819;
        det_ARCA.latitude = 0.633407
        det_ARCA.meridian_convergence_angle = 0.0100733

        det_ANTARES = ROOT.Det()
        deg = 180./pi
        det_ANTARES.set_lonlat( (6+9.942/60.)/deg, (42+47.935/60.)/deg )

        # define the 'periods'
        def addperiod( name, irf, dataset, bg_formula,**kwargs ) :

            common = { 'chan'         : track,
                       'channels'     : [track],
                       'mcfiletypes'  : [numuCC, anumuCC, muonMuon ], 
                       'det'          : ROOT.Det(),
                       'binning'      : (50,5,14,1,8), 
                       'syst_acc'     : 0.3,
                       'syst_PSF'     : 0.5,
                       'psf_mergebins': 4 } 

            common.update( kwargs )

            
            ds = load ( dataset, self.s_min, self.s_max)

            bg = ArcaBackgroundModel( name,ds.hist, bg_formula  )

            self.add_period( name            = name, 
                             veto_runs       = ds.missing_runs(),
                             irf             = irf,
                             datasample      = ds, 
                             backgroundmodel = bg,
                             **common )
        
        '''
        addperiod( "ARCA6good","inputdata/2022428_14_58_DETECTORESPONSE_ARCA6_75.0_zen_CutLevel1_nikhef_b40_pre_upm01_antinoise_upaam01.root","inputdata/DS_P75.0_pre_upm01_antinoise_upaam01.pickle","604.979*exp(-0.5*((x-4.08273)/0.528)*((x-4.08273)/0.528))+54.6157*exp(-0.5*((x-2.75205)/1.17)*((x-2.75205)/1.17))", signal_efficiency = 0.9, det = det_ARCA)
        '''
        
        '''
        addperiod("ARCA21_prova","inputdata/detres_ARCA21_zen_b_40_dst_v4_dynamic_bdttrk_v2_bdtcasc_v8_cut_track_an_bdt095_casc_bdt095.root","inputdata/DS_ARCA21_dst_v4_goldsilver_dynamic_bdttrk_v2_cut_track_an_bdt095_nosparks.pickle","199.98539655808148*exp(-0.5*pow((x-3.0333892447931126)/0.8311531517955035,2))+126.6992776946563*exp(-0.5*pow((x-3.717497074409196)/0.4081747170709585,2))+72.77426215481778*exp(-0.5*pow((x-1.3889354346167384)/0.27031229042785315,2))", mcfiletypes = allmcfiletypes, det = det_ARCA)
        '''
        '''
        addperiod ("ARCA21_prova","inputdata/202475_4_54_35_DETECTORESPONSE_ARCA21_133_tmva_tracks_v10.ALL_zen_b40_pre_upm01_antinoise_upaam01.root","inputdata/DS_P133.all_tmva_tracksC_v9_pre_upm01_antinoise_upaam01.pickle","9900.0*exp(-0.5*pow((x-3.75)/0.5,2))+2000*exp(-0.5*pow((x-2.8)/0.8,2))", mcfiletypes = allmcfiletypes, det = det_ARCA, syst_acc = 0.30, syst_PSF = 0.50 )
        '''
        '''
        addperiod ("ARCA21_prova","inputdata/202475_4_55_12_DETECTORESPONSE_ARCA21_133_tmva_tracks_v10.ALL_zen_b40_pre_upm01_antinoise_len200.root","inputdata/DS_P133.all_tmva_tracks_v9_pre_upm01_antinoise_len200_TestShift.pickle","1033.9462443148202*exp(-0.5*pow((x-3.593382403792754)/0.5532062653963777,2))+9.750099808580014e-07*exp(-0.5*pow((x-4.292137178498489)/3.3293596715027034,2))+316.17663638080796*exp(-0.5*pow((x-2.2435144672644527)/1.1909797511602955,2))", mcfiletypes = allmcfiletypes, det = det_ARCA)

        '''
        print(options.B)
        band =  options.B
        functions_sigmoid  = ['(93.34446907532332+x)/(1+exp(2.660452255180891*x-6.999999290909457)) + (394.9487337400987/(sigma*sqrt(2*pi)))*exp(-power((x-3.471539726708467)/0.6398306223834892,2)/2)','(84.5232943419101+x)/(1+exp(1.901442374723312*x-6.99999999477764)) + (209.9999999592336/(sigma*sqrt(2*pi)))*exp(-power((x-3.539552559551672)/0.5358492201586051,2)/2)',
							  '(46.007713104863726+x)/(1+exp(1.9592611936935997*x-6.999958002696216)) + (208.27389300135758/(sigma*sqrt(2*pi)))*exp(-power((x-3.5800856999053847)/0.5648488242921944,2)/2)',
							  '(21.157315491834506+x)/(1+exp(1.7074160387113217*x-6.91814452892212)) + (201.35812737863176/(sigma*sqrt(2*pi)))*exp(-power((x-3.533225039671027)/0.6067112941056997,2)/2)']
        if band==  "all":
         self.s_min, self.s_max = -1, 1   # all declination!
         
         #all neutrino flavours
         func = "1033.9462443148202*exp(-0.5*pow((x-3.593382403792754)/0.5532062653963777,2))+9.750099808580014e-07*exp(-0.5*pow((x-4.292137178498489)/3.3293596715027034,2))+316.17663638080796*exp(-0.5*pow((x-2.2435144672644527)/1.1909797511602955,2))"
         print("No bands used! Fit of the background without cuts")
        elif band == "adaptive":
          candum = options.c          
          df = pd.read_csv("/sps/km3net/users/lunich/binned/arca-ps-aart_update_bands/inputdata/"+ catalouge_name, sep = "\t", names = ["name", "type", "idk ", "declination", "extended", "galacitc"])
          print(df.declination)
          dec =  df.declination[int(candum)]
          dec_rad = math.radians(dec)
          s_dec = sin(dec_rad)
          print("Dec is: ",dec_rad)
          print("Sin(dec) is: ",sin(dec_rad))
          self.s_min,self.s_max = s_dec-0.2, s_dec+0.2
          #min,max=s_dec-0.2, s_dec+0.2
          if s_dec<-0.8:
            self.s_min = -1
            self.s_max = -0.6
          if s_dec>0.6:
            self.s_max = 0.8
            self.s_min = 0.4
          ds = load("inputdata/DS_P133.all_tmva_tracks_v10_pre_upm01_antinoise_bdt095.pickle", self.s_min, self.s_max)
          #params = fit_histogram(ds)
          #func = get_gauss_expr_string(params)  
          func = ""
        elif band == "auto":
          candum = options.c          
          df =  df = pd.read_csv("/sps/km3net/users/lunich/binned/arca-ps-aart_update_bands/inputdata/"+ catalouge_name, sep = "\t", names = ["name", "type", "idk ", "declination", "extended", "galacitc"])
          print(df.declination)
          dec =  df.declination[int(candum)]
          dec_rad = math.radians(dec)
          s_dec = sin(dec_rad)
          print("Dec is: ",dec_rad)
          print("Sin(dec) is: ",sin(dec_rad))
          if s_dec <=- 0.6:
              
              func=functions_sigmoid[0]

          elif s_dec >- 0.6 and s_dec <= -0.2:
              func=functions_sigmoid[1]
          elif s_dec >-0.2 and s_dec <= 0.2:
              func=functions_sigmoid[2]

          elif s_dec >0.2:
              func=functions_sigmoid[3]
          else: 
              print("Something was wrong in the placement of this source in a band. I will use the traditional fit with no bands")
              func = "1033.9462443148202*exp(-0.5*pow((x-3.593382403792754)/0.5532062653963777,2))+9.750099808580014e-07*exp(-0.5*pow((x-4.292137178498489)/3.3293596715027034,2))+316.17663638080796*exp(-0.5*pow((x-2.2435144672644527)/1.1909797511602955,2))" 

        elif band == "1006":
               self.s_min, self.s_max = -1,-0.6
               gauss_params = gauss_params_dictionary[0]
               func=get_gauss_expr_string(gauss_params)
        elif band == "0602":
               self.s_min, self.s_max = -0.6,-0.2
               gauss_params = gauss_params_dictionary[1]
               func=get_gauss_expr_string(gauss_params)
        elif band == "0202":
               self.s_min, self.s_max = -0.2,0.2
               gauss_params = gauss_params_dictionary[2]
               func=get_gauss_expr_string(gauss_params)
        elif band == "0206":
               self.s_min, self.s_max = 0.2,0.6
               gauss_params = gauss_params_dictionary[3]
               func=get_gauss_expr_string(gauss_params)

        elif band == "hist":
               self.s_min, self.s_max = -1,0.8
               func="hist"
        else:
         print("!!!!ATTENTION!!!!")
         print( "I did NOT recognise the band parameter you gave me. This means that I am going to go with a band wide as all the declinations")
         print("!!!ATTENTION!!!")
         func = "1033.9462443148202*exp(-0.5*pow((x-3.593382403792754)/0.5532062653963777,2))+9.750099808580014e-07*exp(-0.5*pow((x-4.292137178498489)/3.3293596715027034,2))+316.17663638080796*exp(-0.5*pow((x-2.2435144672644527)/1.1909797511602955,2))" 

        addperiod ("ARCA21_track","inputdata/202475_5_53_29_DETECTORESPONSE_ARCA21_133_tmva_tracks_v10.ALL_zen_b40_pre_upm01_antinoise_bdt095.root","inputdata/DS_P133.all_tmva_tracks_v10_pre_upm01_antinoise_bdt095.pickle",func,mcfiletypes = allmcfiletypes, det = det_ARCA,syst_acc = 0.3, syst_PSF = 0.5)
        
        width = "from sin "+str(self.s_min)+" to sin " + str(self.s_max )
        print("band", width, " with background function", func)
             
        if not options.q :

            '''
            addperiod( "ARCA6good","inputdata/2022428_14_58_DETECTORESPONSE_ARCA6_75.0_zen_CutLevel1_nikhef_b40_pre_upm01_antinoise_upaam01.root","inputdata/DS_P75.0_pre_upm01_antinoise_upaam01.pickle","604.979*exp(-0.5*((x-4.08273)/0.528)*((x-4.08273)/0.528))+54.6157*exp(-0.5*((x-2.75205)/1.17)*((x-2.75205)/1.17))", signal_efficiency = 0.9, det = det_ARCA)
            '''

            addperiod( "ARCA6bad",
                        "inputdata/2022428_14_56_DETECTORESPONSE_ARCA6_75.1_zen_CutLevel1_nikhef_b40_pre_upm01_antinoise_upaam01.root",
                        "inputdata/DS_P75.1_pre_upm01_antinoise_upaam01.pickle",
                        "335.326*exp(-0.5*((x-4.11981)/0.522495)*((x-4.11981)/0.522495))+30.925*exp(-0.5*((x-2.93386)/1.11772)*((x-2.93386)/1.11772))",
                        signal_efficiency = 0.9, det = det_ARCA )

            addperiod( "ARCA8",
                    "inputdata/2023113_16_58_DETECTORESPONSE_ARCA8_94.ALL_zen_CutLevel1_nikhef_b40_pre_upm01_antinoise_upaam01.root",
                    "inputdata/DS_P94.ALL_pre_upm01_antinoise_upaam01.pickle",
                       "394.963*exp(-0.5*((x-2.598)/1.0524)*((x-2.598)/1.0524))+3791.88*exp(-0.5*((x-3.89447)/0.503468)*((x-3.89447)/0.503468))", lifetime_data=18346400.0, det = det_ARCA)

#copia arca 8 al posto di quella fuori. e basta così.             
            if options.f=="ARCA21_first":
	 
                # first period of ARCA21 and 
                addperiod( "ARCA19_track","inputdata/detres_ARCA19_zen_b_40_dst_v1_dynamic_bdttrk_v1_bdtcasc_v1_cut_track_an_bdt095_casc_bdt095.root","inputdata/DS_ARCA19_dst_v1_goldsilver_dynamic_bdttrk_v1_cut_track_an_bdt095.pickle","72.99999976151094*exp(-0.5*pow((x-3.5817783977388475)/0.43656084263473616,2))+82.9999993518388*exp(-0.5*pow((x-2.9439618440730775)/0.8564133538587623,2))+29.9999993600465*exp(-0.5*pow((x-1.5600816090021006)/0.30468350260232924,2))",
                        mcfiletypes = allmcfiletypes, det = det_ARCA )

                addperiod ("ARCA21","inputdata/detres_ARCA21_zen_b_40_dst_v4_dynamic_bdttrk_v2_bdtcasc_v8_cut_track_an_bdt095_casc_bdt095.root","inputdata/DS_ARCA21_dst_v4_goldsilver_dynamic_bdttrk_v2_cut_track_an_bdt095_nosparks.pickle","199.98539655808148*exp(-0.5*pow((x-3.0333892447931126)/0.8311531517955035,2))+126.6992776946563*exp(-0.5*pow((x-3.717497074409196)/0.4081747170709585,2))+72.77426215481778*exp(-0.5*pow((x-1.3889354346167384)/0.27031229042785315,2))",
                    mcfiletypes = allmcfiletypes, det = det_ARCA )


              
                addperiod( "ARCA19_track",
                        "inputdata/detres_ARCA19_zen_b_40_dst_v0_bdttrk_v0_bdtcasc_v0_cut_track_an_bdt095_casc_bdt_Edependent.root",
                        "inputdata/DS_ARCA19_dst_v0_dynamic_bdttrk_v0_bdtcasc_v0_cut_track_an_bdt095_casc_bdt_Edependent.pickle",
                        "", det = det_ARCA)
            
            
                addperiod ("ARCA21_track",
                    "inputdata/detres_ARCA21_zen_b_40_dst_v3_dynamic_bdttrk_v1_bdtcasc_v7_cut_track_an_bdt095_casc_bdt_Edependent.root",
                    "inputdata/DS_ARCA21_dst_v3_dynamic_dttrk_v1_bdtcasc_v7_cut_track_an_bdt095_casc_bdt_Edependent.pickle",
                  "", det = det_ARCA )
                
            elif options.f=="ARCA21_all_period_upaam01":


                addperiod( "ARCA19_track","inputdata/202475_4_30_20_DETECTORESPONSE_ARCA19_116_tmva_tracks_v5.ALL_zen_b40_pre_upm01_antinoise_upaam01.root","inputdata/DS_P116.all_tmva_tracks_v9_pre_upm01_antinoise_upaam01_TestShift.pickle","1500*exp(-0.5*pow((x-3.75)/0.46,2))+300.0*exp(-0.5*pow((x-2.8)/1.12,2))", mcfiletypes = allmcfiletypes, lifetime_data=4181454.0,det = det_ARCA )
                addperiod ("ARCA21_track","inputdata/202475_4_54_35_DETECTORESPONSE_ARCA21_133_tmva_tracks_v10.ALL_zen_b40_pre_upm01_antinoise_upaam01.root","inputdata/DS_P133.all_tmva_tracks_v9_pre_upm01_antinoise_upaam01_TestShift.pickle","8841.421661086519*exp(-0.5*pow((x-3.735473879066079)/-0.5255379644753796,2))+3.5762317816844837e-09*exp(-0.5*pow((x-6.908495391664351)/-117.77721701970478,2))+1336.4352913608418*exp(-0.5*pow((x-2.3678892278672343)/-1.1440334344037448,2))", mcfiletypes = allmcfiletypes,lifetime_data= 24831400.0 ,det = det_ARCA )

            elif options.f=="ARCA21_all_period_len200":
            
                addperiod( "ARCA19_track","inputdata/202475_4_31_34_DETECTORESPONSE_ARCA19_116_tmva_tracks_v5.ALL_zen_b40_pre_upm01_antinoise_len200.root","inputdata/DS_P116.all_tmva_tracks_v9_pre_upm01_antinoise_len200_TestShift.pickle","55.5*exp(-0.5*pow((x-2.34)/1.00164,2))+177.17*exp(-0.5*pow((x-3.553)/0.4595,2))+30.63*exp(-0.5*pow((x-3.509)/0.275,2))", mcfiletypes = allmcfiletypes,lifetime_data=4181454.0,det = det_ARCA )

                addperiod ("ARCA21_track","inputdata/202475_4_54_35_DETECTORESPONSE_ARCA21_133_tmva_tracks_v10.ALL_zen_b40_pre_upm01_antinoise_upaam01.root","inputdata/DS_P133.all_tmva_tracks_v9_pre_upm01_antinoise_len200_TestShift.pickle","1245.5129563979392*exp(-0.5*pow((x-3.525172964630171)/0.5016495499860752,2))+1.4037087818605416e-07*exp(-0.5*pow((x-6.775223646445769)/2.8817975300126677,2))+380.57875293744286*exp(-0.5*pow((x-2.2803219008040454)/1.0579923893616108,2))", mcfiletypes = allmcfiletypes, lifetime_data= 24831400.0,det = det_ARCA )

            elif options.f=="ARCA_21_all_period_bdt095":

                addperiod( "ARCA19_track","inputdata/202475_4_37_33_DETECTORESPONSE_ARCA19_116_tmva_tracks_v5.ALL_zen_b40_pre_upm01_antinoise_bdt095.root","inputdata/DS_P116.all_tmva_tracks_v9_pre_upm01_antinoise_bdt095_TestShift.pickle","96.62811065789774*exp(-0.5*pow((x-3.5661666227811875)/0.3812943029755955,2))+70.03654742954299*exp(-0.5*pow((x-2.9999997985326217)/0.9516548043397174,2))+18.222917969960314*exp(-0.5*pow((x-0.3366699529083293)/0.7310122642257921,2))", mcfiletypes = allmcfiletypes,lifetime_data=4181454.0,det = det_ARCA )

                addperiod ("ARCA21_track","inputdata/202475_5_53_29_DETECTORESPONSE_ARCA21_133_tmva_tracks_v10.ALL_zen_b40_pre_upm01_antinoise_bdt095.root","inputdata/DS_P133.all_tmva_tracks_v9_pre_upm01_antinoise_bdt095_TestShift.pickle","1033.9462443148202*exp(-0.5*pow((x-3.593382403792754)/0.5532062653963777,2))+9.750099808580014e-07*exp(-0.5*pow((x-4.292137178498489)/3.3293596715027034,2))+316.17663638080796*exp(-0.5*pow((x-2.2435144672644527)/1.1909797511602955,2))", mcfiletypes = allmcfiletypes,lifetime_data= 24831400.0,det = det_ARCA )

            elif options.f=="ARCA_21_test_energy":

                addperiod ("ARCA21_track","inputdata/202475_5_53_29_DETECTORESPONSE_ARCA21_133_tmva_tracks_v10.ALL_zen_b40_pre_upm01_antinoise_bdt095.root","inputdata/DS_P133.all_tmva_tracks_v9_pre_upm01_antinoise_bdt095_TestShift.pickle","1033.9462443148202*exp(-0.5*pow((x-3.593382403792754)/0.5532062653963777,2))+9.750099808580014e-07*exp(-0.5*pow((x-4.292137178498489)/3.3293596715027034,2))+316.17663638080796*exp(-0.5*pow((x-2.2435144672644527)/1.1909797511602955,2))", mcfiletypes = allmcfiletypes,lifetime_data= 24831400.0,det = det_ARCA )

                
            if options.A=="ANTARES_2020":
        
                addperiod( "ANTAREStrack","inputdata/2023_4_17_21_2_DETECTORESPONSE_ANTARES_zen_CutLevel1_local_b40_watm_ANTARES.root","inputdata/DS_Adata_all4tracks_2020.pickle","1378.12*exp(-0.5*((x-2.71)/0.143)*((x-2.71)/0.143))+536.64*exp(-0.5*((x-3.062)/0.244)*((x-3.062)/0.244))+219.03*exp(-0.5*((x-3.505)/0.414)*((x-3.505)/0.414))", det = det_ANTARES) 

            elif options.A=="ANTARES_2022":

                addperiod( "ANTAREStrack","inputdata/2024_5_21_15_4_DETECTORESPONSE_ANTARES_TRACK_HONDA_zen_CutLevel1_local_b40_watm_honda_nutau.root","inputdata/DS_Adata_all4tracks.pickle","1666.81*exp(-0.5*pow((x-2.72)/0.141,2))+630.34*exp(-0.5*pow((x-3.071)/0.2489,2))+ 264.90*exp(-0.5*pow((x-3.505)/0.419,2))", det = det_ANTARES, binning = (20,10, 12, 0, 12), syst_acc = 0.15, syst_PSF = 0.15)
    

            elif options.A=="ANTARES_track_shower_2022":

                #addperiod("ANTARESshower","inputdata/2024_1_26_16_12_DETECTORESPONSE_ANTARESSH_zen_CutLevel1_local_b40_honda_tEreco.root","inputdata/DS_Adata_all4showers_Ereco.pickle","6.30e+02*TMath::Landau(x,1.01,8.5e-02)",det=det_ANTARES,syst_acc = 0.15, syst_PSF = 0.15, binning = (20,10, 12, 0, 12))
                
                addperiod("ANTARES_shower","inputdata/2024_1_26_16_12_DETECTORESPONSE_ANTARESSH_zen_CutLevel1_local_b40_honda_tEreco.root","inputdata/DS_Adata_all4showers_Ereco.pickle","5.03568e+01*exp(-0.5*pow((x-3.62471)/0.327182,2.))+5.05382*exp(-0.5*pow((x-4.41389)/0.728542,2.))",det=det_ANTARES,syst_acc = 0.15, syst_PSF = 0.15, binning = (20,10, 12, 0, 12))
                
                addperiod( "ANTAREStrack","inputdata/2024_5_21_15_4_DETECTORESPONSE_ANTARES_TRACK_HONDA_zen_CutLevel1_local_b40_watm_honda_nutau.root","inputdata/DS_Adata_all4tracks.pickle","1666.81*exp(-0.5*pow((x-2.72)/0.141,2))+630.34*exp(-0.5*pow((x-3.071)/0.2489,2))+ 264.90*exp(-0.5*pow((x-3.505)/0.419,2))", det = det_ANTARES,syst_acc = 0.15, syst_PSF = 0.15, binning = (20,10, 12, 0, 12))


            elif options.A=="ANTARES_check":

                addperiod( "ANTAREStrack","inputdata/2025_4_8_10_37_DETECTORESPONSE_ANTARES_TRACK_HONDA_zen_CutLevel1_local_b40_watm_honda_nutau.root","inputdata/DS_Adata_all4tracks.pickle","1666.81*exp(-0.5*pow((x-2.72)/0.141,2))+630.34*exp(-0.5*pow((x-3.071)/0.2489,2))+ 264.90*exp(-0.5*pow((x-3.505)/0.419,2))", det = det_ANTARES,syst_acc = 0.15, syst_PSF = 0.15, binning = (20,10, 12, 2, 8))

                
            elif options.A=="No":
                print("No ANTARES data")
            
        #---------------
        # number of PEs
        #---------------

        print("Periods are: ",self.periods)
        
        self.npe_h0       = 80000  # number of PEs to run
        self.npe_h1       = 10000
        self.npe_sigcheck = 1
        self.npe_nsigs    = list( aa.frange( 0, 10, 1.0 ) )
        self.make_likelihood_curve = options.L #lascia così per ora)
        self.make_data_scramble = options.T     #idem
##################################################################################################
#                                           MAIN                                                 #
##################################################################################################
if __name__ == "__main__" :


    #options = aa.Options( __doc__, sys.argv[1:])

    start_time = time.time()
    print("--- Code start: %s seconds ---" % (time.time() - start_time)) 
    print(aa.__file__)
    #multiprocessing.set_start_method('spawn')
    name = options.n

    if options.c != "none" : # run PE's for a single candidates

        seq = options.c

        if '...' in seq : 
            cands = range( *map(int, seq.split("...")))
        else :
            cands = map(int, options.c.split(","))

        for candnum in cands:
            print(candnum)
#	CONTROLLA QUA ATTENTOI AAAAAAAAAAAAA
            ana = ArcaBinnedPointSourceAnalysis( name = name +"_cand_"+str(candnum))
            print(options.p)
            ana.picklepath = options.p  #same of the  seccond
            ana.run_cand_search( ana.candlist.sources[candnum] )
            for k,v in ana.cand_searches.items():
                #print("sono in cycle k",k)
                print(v.summary_table())
            print( ana.limit_table_cands() )
            print("Checks on p-value",ana.post_p(0.1))
            ana.skymap(ana.picklepath+"/skymap_"+str(candnum)+".png")
            #ana.save()
            ana.write_rootfile() #update path.  Usciranno due file root, uno sulle pdf e uno su istogramma
	#ogni volta cambia il path ogni volta per avere in cartelle diverse . Oppure automatizzalo
#                                               destination= options.p
            ana.write_webpage(destination=options.p)

        print("--- Code finished: %s seconds ---" % (time.time() - start_time))
        sys.exit()

    if options.d : # debug 
        ana = pickle.load( open( outputpath +"/"+ options.n +"_ALL.pickle", 'rb') )
        s = list(ana.cand_searches.values())[0]

        ddres = s.periods[0].detres
        
        for flav,v in ddres.response_map :  # 1st period, muons->tracks
            print(k,v)
            for channel,detresponse in v :
                print (flav,channel,detresponse )
        
    if options.a : # aggregate/combine
        ana = ArcaBinnedPointSourceAnalysis( name = options.n +"_ALL" )
        L = glob.glob( ana.picklepath +"/"+ options.n +"_cand_*.pickle" )
        print ( options.n +"_cand_*.pickle" )
        print(L)
        for f in L :
            print ("loading ",f )
            a = pickle.load( open( f,'rb') )
            print(a)
            ana.cand_searches.update( a.cand_searches )
        print ( ana.limit_table_cands() )
        ana.save()
        sys.exit()

    if options.s :
        ana = pickle.load( open( outputpath +"/"+ options.n +"_ALL.pickle", 'rb') )
        ana.scrabmle = options.s
        ana.init_datasamples()
        ana.produce_results()
        print ( ana.limit_table_cands() )
        ana.save()
        #sys.exit()

    if options.w :
        #ana = pickle.load( open( outputpath +"/"+options.n +"_ALL.pickle", 'rb') ) 
        ana.write_webpage( destination  =  "t61@login.nikhef.nl:~/public_html/arcapnt2/" )
