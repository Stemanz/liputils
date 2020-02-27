# -*- coding: utf-8 -*-
# Author: Manzini Stefano; stefano.manzini@gmail.com
# 240220

import re
import pandas as pd
import numpy as np
from numpy import average

class Lipid:
    """
    A class for lipids. Use like lipid = Lipid("CE 12:0")

    It defaults to: pmol / μg (solid samples) or μM (plasma), which is basically
    the same (μM → pmol / μl). The volume is not considered, so the total number
    of molecules is computed per normalized volume unit.
    
    Parameters
    ==========
    
    amount: <int> amount of picomoles per microgram or picomoles per microliter
        defaults to 0

    unit: <str> Choose from: "moles", "millimoles", "micromoles", "nanomoles",
        "picomoles", "femtomoles", "attomoles", "zeptomoles". The exact number of molecules
        can be accessed via the attribute self.molecules. defaults to "picomoles"
    """
    def __init__(self, lipid, *, amount = 0, unit="picomoles",
                 **kwargs):
    
        # some useful info
        # ================
        mole = 6.022140857e+23
        picomole = 6.022140857e+11
        
        unit_conversion = {
            "moles":      6.022140857e+23,
            "millimoles": 6.022140857e+20,
            "micromoles": 6.022140857e+17,
            "nanomoles":  6.022140857e+14,
            "picomoles":  6.022140857e+11,
            "femtomoles": 6.022140857e+8,
            "attomoles":  6.022140857e+5,
            "zeptomoles": 6.022140857e+2
        }

        double_bond_mass_reduction = 2.0155 # (-CH2-CH2- to -C=C-); also Sigma
        mass_difference_per_additional_C_atom = 14.0156 # (-CH2- to -CH2-CH2-)
        cholesterol = 386.3549   # Sigma-Aldrich: 386.65
        glycerol = 92.09    # Sigma-Aldrich
        C12_0 = 200.32  # Sigma-aldrich, 12:0 dodecanoic acid
        C13_0 = 214.34  # Sigma-aldrich, 13:0 tridecanoic acid
        C14_0 = 231.39  # Sigma-aldrich, 14:0 myristic acid
        C15_0 = 242.40  # Sigma-aldrich, 15:0 pentadecanoic acid
        C16_0 = 256.42  # Sigma-aldrich, 16:0 palmitic acid
        C16_1 = 254.41  # Sigma-aldrich, 16:1 palmitoleic acid
        C17_0 = 270.45  # Sigma-aldrich, 17:0 margaric acid
        C18_0 = 284.48  # Sigma-aldrich, 18:0 stearic acid
        C18_1 = 282.46  # Sigma-aldrich, 18:1 vaccenic acid
        C18_2 = 280.45  # Sigma-aldrich, 18:2 linolenic acid
        C18_3 = 278.4345    # predicted
        C19_0 = 298.50  # Sigma-aldrich, 19:0 nonadecanoic acid
        C20_0 = 312.53  # Sigma-aldrich, 20:0 arachidic acid
        C20_1 = 310.51  # Sigma-aldrich, 20:1 gondoic acid
        C20_2 = 308.4945    # predicted
        C20_3 = 306.479 # predicted
        C21_0 = 326.56  # Sigma-aldrich, 21:0 heneicosanoic acid
        C22_0 = 340.58  # Sigma-aldrich, 22:0 behenic acid
        C23_0 = 354.5956    # predicted
        C24_0 = 368.6112   # predicted
        C22_1 = 338.5645    # predicted
        C22_6 = 328.4869    # predicted
        C26_0 = 396.6423    # predicted
        C27_0 = 410.6579    # predicted (does this exist?)
        C26_1 = 394.6268    # predicted
        C28_1 = 422.658 #predicted (does this exist?)
        # longer saturated fatty acids for total TAG content
        H = 14.0156 / 14    # (-CH2- / 14)
        H2O = 18.02 # Sigma-Aldrich
        PEP_backbone = 162.5847 # predicted
        Gb3_backbone = 484.7905 # predicted
        Glc_GalCer_backbone = 160.6848  # predicted
        
        self.amount = amount    # in the specified unit

        self.molecules = amount * unit_conversion.get(unit, 0)  # number of molecules
        
        self.name = lipid
        
        lipid_weights = {   # source: http://www.lipidmaps.org/
            "CE 12:0" : 568.5219,
            "CE 13:0" : 582.5375,   # average between CE 12:0 and CE 14:0
            "CE 14:0" : 596.5532,
            "CE 14:1" : 594.5376,
            "CE 15:0" : 610.5689,
            "CE 15:1" : 608.5532,
            "CE 16:0" : 624.5845,
            "CE 16:1" : 622.5689,
            "CE 17:0" : 638.6002,
            "CE 17:1" : 636.5845,
            "CE 18:0" : 652.6158,
            "CE 18:1" : 650.6002,
            "CE 18:2" : 648.5845,
            "CE 18:3" : 646.5689,
            "CE 19:0" : 666.6315,
            "CE 19:1" : 664.616,    # predicted
            "CE 19:2" : 662.6005,   # predicted
            "CE 20:0" : 680.6471,
            "CE 20:1" : 678.6315,
            "CE 20:2" : 676.6158,
            "CE 20:3" : 674.6002,
            "CE 20:4" : 672.5845,
            "CE 20:5" : 670.5689,
            "CE 21:0" : 694.6627,   # predicted
            "CE 21:1" : 692.6472,   # predicted
            "CE 22:0" : 708.6784,
            "CE 22:1" : 706.6628,
            "CE 22:2" : 704.6471,
            "CE 22:3" : 702.6315,
            "CE 22:4" : 700.6158,
            "CE 22:5" : 698.6002,
            "CE 22:6" : 696.5845,
            "CE 24:0" : 736.7096,   # predicted
            "CE 24:1" : 734.6941,
            "CE 24:2" : 732.6786,   # predicted
            "CE 24:3" : 730.6631,   # predicted
            "CE 24:4" : 728.6476,   # predicted
            "CE 24:5" : 726.6321,   # predicted
            "CE 24:6" : 724.6166,   # predicted
            "CE 26:0" : 764.7408,   # predicted
            "CE 26:1" : 762.7253,   # predicted
            "CE 26:2" : 760.7098,   # predicted
            "Cer(d18:0/16:0)" : 539.5277,
            "Cer(d18:0/18:0)" : 567.5590,
            "Cer(d18:0/20:0)" : 595.5903,
            "Cer(d18:0/22:0)" : 785.6745,
            "Cer(d18:0/24:0)" : 651.6529,
            "Cer(d18:0/24:1)" : 649.6373,
            "Cer(d18:1/16:0)" : 537.5121,
            "Cer(d18:1/18:0)" : 565.5434,
            "Cer(d18:1/18:1)" : 563.5277,
            "Cer(d18:1/20:0)" : 593.5747,
            "Cer(d18:1/22:0)" : 649.6371,   # predicted
            "Cer(d18:1/22:1)" : 619.5904,   # predicted
            "Cer(d18:1/23:0)" : 635.6216,
            "Cer(d18:1/24:0)" : 649.6373,
            "Cer(d18:1/24:1)" : 647.6216,
            "Cer(d18:1/26:0)" : 677.6686,
            "Cer(d18:1/26:1)" : 675.6529,
            "DAG 14:0/14:0" : 518.8299, # glycerol - 2*H20 + 2*C14_0
            "DAG 14:0/16:0" : 543.86,   # glycerol - 2*H20 + C14_0 + C16_0
            "DAG 14:0/18:1" : 569.9,    # glycerol - 2*H20 + C14_0 + C18_1
            "DAG 14:0/18:2" : 567.89,   # glycerol - 2*H20 + C14_0 + C18_2
            "DAG 16:0/16:0" : 568.89,   # glycerol - 2*H20 + 2*C16_0
            "DAG 16:0/16:1" : 566.88,   # glycerol - 2*H20 + C16_0 + C16_1
            "DAG 16:0/18:1" : 594.93,   # glycerol - 2*H20 + C16_0 + C18_1
            "DAG 16:0/18:2" : 592.92,   # glycerol - 2*H20 + C16_0 + C18_2
            "DAG 16:0/20:2" : 620.9645, # glycerol - 2*H20 + C16_0 + C20_2
            "DAG 16:0/20:3" : 618.949,  # glycerol - 2*H20 + C16_0 + C20_3
            "DAG 16:0/22:6" : 640.9569, # glycerol - 2*H20 + C16_0 + C22_6
            "DAG 16:1/16:1" : 564.87,   # glycerol - 2*H20 + 2*C16_1
            "DAG 16:1/18:0" : 594.94,   # glycerol - 2*H20 + C16_1 + C18_0
            "DAG 16:1/18:2" : 590.91,   # glycerol - 2*H20 + C16_1 + C18_2
            "DAG 18:0/18:1" : 622.99,   # glycerol - 2*H20 + C18_0 + C18_1
            "DAG 18:0/18:2" : 620.98,   # glycerol - 2*H20 + C18_0 + C18_2
            "DAG 18:0/20:3" : 647.009,  # glycerol - 2*H20 + C18_0 + C20_3
            "DAG 18:1/18:1" : 620.9699, # glycerol - 2*H20 + 2*C18_1
            "DAG 18:1/18:2" : 618.96,   # glycerol - 2*H20 + C18_1 + C18_2
            "DAG 18:1/20:1" : 649.02,   # glycerol - 2*H20 + C18_1 + C20_1
            "DAG 18:1/20:3" : 644.989,  # glycerol - 2*H20 + C18_1 + C20_3
            "DAG 18:1/22:6" : 666.9969, # glycerol - 2*H20 + C18_1 + C22_6
            "DAG 18:2/18:2" : 616.9499, # glycerol - 2*H20 + 2*C18_2
            "DAG 18:2/18:3" : 614.9345, # glycerol - 2*H20 + C18_2 + C18_3
            "DAG 18:2/20:1" : 647.01,   # glycerol - 2*H20 + C18_2 + C20_1
            "DAG 18:2/22:6" : 664.9869, # glycerol - 2*H20 + C18_2 + C22_6
            "FC" : 386.3549,
            "Gb3(d18:1/16:0)" : 1023.670565,  # http://www.swisslipids.org/#/entity/SLM:000487911/?chemInfo
            "Gb3(d18:1/18:0)" : 1051.7305,  # Gb3_backbone + C18_1 + C18_0
            "Gb3(d18:1/20:0)" : 1079.7804,  # Gb3_backbone + C18_1 + C18_0
            "Gb3(d18:1/22:0)" : 1107.8305,  # Gb3_backbone + C18_1 + C22_0
            "Gb3(d18:1/22:1)" : 1105.815,  # Gb3_backbone + C18_1 + C22_1
            "Gb3(d18:1/23:0)" : 1121.8461,  # Gb3_backbone + C18_1 + C23_0
            "Gb3(d18:1/24:0)" : 1135.8617,  # Gb3_backbone + C18_1 + C24_0
            "Gb3(d18:1/24:1)" : 1133.8462,  # predicted
            "Gb3(d18:1/26:1)" : 1161.8773,  # Gb3_backbone + C18_1 + C26_1
            "Glc/GalCer(d18:1/16:0)" : 699.5649,
            "Glc/GalCer(d18:1/18:0)" : 727.5962,
            "Glc/GalCer(d18:1/20:0)" : 755.6275,
            "Glc/GalCer(d18:1/22:0)" : 783.6588,
            "Glc/GalCer(d18:1/22:1)" : 781.6433,    # predicted
            "Glc/GalCer(d18:1/23:0)" : 797.6745,
            "Glc/GalCer(d18:1/24:0)" : 811.6901,
            "Glc/GalCer(d18:1/24:1)" : 809.6745,
            "Glc/GalCer(d18:1/26:0)" : 839.7214,
            "Glc/GalCer(d18:1/26:1)" : 837.7058,
            "LacCer(d18:1/16:0)" : 861.6177,
            "LacCer(d18:1/18:0)" : 889.6490,
            "LacCer(d18:1/20:0)" : 917.6803,
            "LacCer(d18:1/22:0)" : 945.7116,
            "LacCer(d18:1/22:1)" : 973.7428,    # predicted
            "LacCer(d18:1/23:0)" : 959.7271,    # predited
            "LacCer(d18:1/24:0)" : 973.7429,
            "LacCer(d18:1/24:1)" : 971.7273,
            "LacCer(d18:1/26:0)" : 1001.7742,
            "LacCer(d18:1/26:1)" : 999.7586,
            "LPC 14:0" : 467.577,   # MW, https://avantilipids.com/product/855575/
            "LPC 16:0" : 495.630,
            "LPC 16:1" : 493.6145,  # predicted
            "LPC 17:1" : 507.641,   # https://avantilipids.com/product/lm1601/
            "LPC 18:0" : 523.683,   # https://avantilipids.com/product/855775/
            "LPC 18:1" : 521.667,   # https://avantilipids.com/product/845875/
            "LPC 18:2" : 519.6515,  # predicted
            "LPC 18:3" : 517.636,   # predicted
            "LPC 19:0" : 537.710,   # https://avantilipids.com/product/855776/
            "LPC 20:0" : 551.736,   # https://avantilipids.com/product/855777/
            "LPC 20:1" : 549.7205,  # predicted
            "LPC 20:2" : 547.705,   # predicted
            "LPC 20:3" : 545.6895,  # predicted
            "LPC 20:4" : 543.674,   # predicted
            "LPC 20:5" : 541.6585,  # predicted
            "LPC 22:0" : 902.358,   # https://avantilipids.com/product/850371/
            "LPC 22:1" : 900.3425,  # predicted
            "LPC 22:5" : 892.2805,  # predicted
            "LPC 22:6" : 890.265,   # predicted
            "LPE 16:0" : 453.550,   # https://avantilipids.com/product/856705/
            "LPE 18:0" : 481.603,   # https://avantilipids.com/product/856715/
            "LPE 18:1" : 479.5875,  # predicted
            "LPE 18:2" : 477.572,   # predicted
            "LPI 18:1" : 615.691,   # https://avantilipids.com/product/850100/
            "PA 16:0/16:0" : 648.4730,
            "PA 16:0/18:1" : 674.4887,
            "PA 18:1/18:1" : 700.5043,
            "PC 14:0/14:0" : 677.4996,
            "PC 14:0/16:0" : 705.5309,
            "PC 14:0/16:1" : 703.5152,
            "PC 14:0/18:1" : 731.5465,
            "PC 14:0/18:2" : 729.5309,
            "PC 16:0/16:0" : 733.5622,
            "PC 16:0/16:1" : 731.5465,
            "PC 16:0/17:1" : 745.5622,
            "PC 16:0/18:0" : 761.5935,
            "PC 16:0/18:1" : 759.5778,
            "PC 16:0/18:2" : 757.5622,
            "PC 16:0/18:3" : 755.5465,
            "PC 16:0/20:0" : 789.6248,
            "PC 16:0/20:1" : 787.6091,
            "PC 16:0/20:2" : 785.5935,
            "PC 16:0/20:3" : 783.5778,
            "PC 16:0/20:4" : 781.5622,
            "PC 16:0/20:5" : 779.5465,
            "PC 16:0/22:4" : 809.5935,
            "PC 16:0/22:5" : 807.5778,
            "PC 16:0/22:6" : 805.5622,
            "PC 16:1/16:1" : 729.5309,
            "PC 16:1/17:0" : 745.5622,
            "PC 16:1/17:1" : 743.5465,
            "PC 16:1/18:0" : 759.5778,
            "PC 16:1/18:1" : 757.5622,
            "PC 16:1/18:2" : 755.5465,
            "PC 16:1/18:3" : 753.5309,
            "PC 16:1/20:1" : 785.5935,
            "PC 16:1/20:3" : 781.5622,
            "PC 16:1/20:4" : 779.5465,
            "PC 16:1/20:5" : 777.5309,
            "PC 16:1/22:6" : 803.5465,
            "PC 17:0/18:1" : 773.5935,
            "PC 17:0/18:2" : 771.5778,
            "PC 17:0/20:3" : 797.5935,
            "PC 17:0/20:4" : 795.5778,
            "PC 17:0/22:6" : 819.5778,
            "PC 17:1/18:0" : 773.5935,
            "PC 17:1/18:1" : 771.5778,
            "PC 17:1/18:2" : 769.5622,
            "PC 18:0/18:1" : 787.6091,
            "PC 18:0/18:2" : 785.5935,
            "PC 18:0/18:3" : 783.5778,
            "PC 18:0/20:1" : 815.6404,
            "PC 18:0/20:2" : 813.6248, #!
            "PC 18:0/20:3" : 811.6091,
            "PC 18:0/20:4" : 809.5935,
            "PC 18:0/20:5" : 807.5778,
            "PC 18:0/22:4" : 837.6248,
            "PC 18:0/22:6" : 833.5935,
            "PC 18:1/18:1" : 785.5935,
            "PC 18:1/18:2" : 783.5778,
            "PC 18:1/18:3" : 781.5622,
            "PC 18:1/20:0" : 815.6404,
            "PC 18:1/20:1" : 813.6248, #!
            "PC 18:1/20:2" : 811.6091,
            "PC 18:1/20:3" : 809.5935,
            "PC 18:1/20:4" : 807.5778,
            "PC 18:1/20:5" : 805.5622,
            "PC 18:1/22:5" : 833.5935,
            "PC 18:1/22:6" : 831.5778,
            "PC 18:2/18:2" : 781.5622,
            "PC 18:2/18:3" : 779.5465,
            "PC 18:2/20:0" : 813.6248,
            "PC 18:2/20:1" : 811.6091,
            "PC 18:2/20:2" : 809.5935,
            "PC 18:2/20:3" : 807.5778,
            "PC 18:2/20:4" : 805.5622,
            "PC 18:2/20:5" : 803.5465,
            "PC 18:2/22:6" : 829.5622,
            "PC 20:1/20:4" : 835.6091,
            "PC 20:3/20:3" : 833.5935,
            "PC 20:3/20:4" : 831.5778,
            "PC 20:4/20:4" : 829.5622,
            "PC O-18:0/16:0" : 747.6142,
            "PC P-16:0/16:1 (PC O-16:1/16:1)" : 715.5516,
            "PC P-16:0/18:2 (PC O-16:1/18:2)" : 741.5672,
            "PC P-18:0/16:0 (PC O-18:1/16:0)" : 745.5985,
            "PC P-18:0/18:2 (PC O-18:1/18:2)" : 769.5985,
            "PC P-18:0/20:3 (PC O-18:1/20:3)" : 795.6142,
            "PC P-18:0/20:4 (PC O-18:1/20:4)" : 793.5985,
            "PC P-18:0/22:5 (PC O-18:1/22:5)" : 819.6139,   # predicted
            "PC P-18:0/22:6 (PC O-18:1/22:6)" : 817.5985,
            "PE 14:0/18:1" : 689.4996,
            "PE 16:0/16:0" : 691.5152,
            "PE 16:0/16:1" : 689.4996,
            "PE 16:0/17:1" : 703.5152,
            "PE 16:0/18:1" : 717.5309,
            "PE 16:0/18:2" : 715.5152,
            "PE 16:0/18:3" : 713.4996,
            "PE 16:0/20:2" : 743.5465,
            "PE 16:0/20:3" : 741.5309,
            "PE 16:0/20:4" : 739.5152,
            "PE 16:0/20:5" : 737.4996,
            "PE 16:0/22:4" : 767.5465,
            "PE 16:1/16:1" : 687.4839,
            "PE 16:1/18:0" : 717.5309,
            "PE 16:1/18:1" : 757.5622,
            "PE 16:1/18:2" : 670.4574,
            "PE 16:1/20:3" : 781.5622,
            "PE 16:1/20:4" : 779.5465,
            "PE 16:1/22:6" : 803.5465,
            "PE 17:0/18:1" : 731.5465,
            "PE 17:0/18:2" : 729.5309,
            "PE 17:0/20:4" : 753.5309,
            "PE 17:1/18:1" : 729.5309,
            "PE 17:1/18:2" : 727.5152,
            "PE 17:1/20:4" : 751.5152,
            "PE 18:0/18:1" : 745.5622,
            "PE 18:0/18:2" : 743.5465,
            "PE 18:0/20:1" : 773.5935,
            "PE 18:0/20:2" : 771.5778,
            "PE 18:0/20:3" : 769.5622,
            "PE 18:0/20:4" : 767.5465,
            "PE 18:0/20:5" : 765.5309,
            "PE 18:0/22:3" : 797.5934,  # predicted
            "PE 18:0/22:4" : 795.5778,
            "PE 18:0/22:6" : 791.5465,
            "PE 18:1/18:1" : 743.5465,
            "PE 18:1/18:2" : 741.5309,
            "PE 18:1/18:3" : 739.5152,
            "PE 18:1/20:1" : 771.5778,
            "PE 18:1/20:2" : 769.5622,
            "PE 18:1/20:3" : 767.5465,
            "PE 18:1/20:4" : 765.5309,
            "PE 18:1/20:5" : 763.5152,
            "PE 18:1/22:4" : 793.5622,
            "PE 18:1/22:5" : 791.5467,  # predicted
            "PE 18:1/22:6" : 789.5309,
            "PE 18:2/18:2" : 739.5152,
            "PE 18:2/20:4" : 763.5152,
            "PE 20:0/20:4" : 795.5778,
            "PE 20:1/20:1" : 799.6091,
            "PE 20:1/20:4" : 793.5622,
            "PE 22:6/22:6" : 835.5152,
            "PE O-18:0/16:0" : 705.5672,
            "PE O-18:0/18:1" : 731.5829,
            "PE O-18:0/20:4" : 753.5672,
            "PE O-20:0/18:2" : 757.5985,
            "PE P-16:0/18:1 (PE O-16:1/18:1)" : 701.5359,
            "PE P-16:0/18:2 (PE O-16:1/18:2)" : 699.5203,
            "PE P-16:0/20:3 (PE O-16:1/20:3)" : 725.5359,
            "PE P-16:0/20:4 (PE O-16:1/20:4)" : 723.5203,
            "PE P-16:0/20:5 (PE O-16:1/20:5)" : 721.5046,
            "PE P-16:0/22:6 (PE O-16:1/22:6)" : 747.5203,
            "PE P-18:0/16:0 (PE O-18:1/16:0)" : 703.5516,
            "PE P-18:0/18:1 (PE O-18:1/18:1)" : 729.5672,
            "PE P-18:0/18:2 (PE O-18:1/18:2)" : 727.5516,
            "PE P-18:0/20:3 (PE O-18:1/20:3)" : 753.5672,
            "PE P-18:0/20:4 (PE O-18:1/20:4)" : 751.5516,
            "PE P-18:0/22:5 (PE O-18:1/22:5)" : 777.5671,   # predicted
            "PE P-18:0/22:6 (PE O-18:1/22:6)" : 775.5516,
            "PE P-18:1/16:0" : 701.4647,    # predicted
            "PE P-18:1/18:1" : 727.5047,    # predicted
            "PE P-18:1/18:2" : 725.4947,    # predicted
            "PE P-18:1/20:4" : 749.5081,    # predicted
            "PE P-18:1/22:6" : 773.5316,    # predicted
            "PE P-20:0/16:0 (PE O-20:1/16:0)" : 731.5346,   # predicted
            "PE P-20:0/20:4 (PE O-20:1/20:4)" : 779.5826,   # predicted
            "PG 16:0/18:1" : 748.5254,
            "PG 16:1/18:1" : 746.5098,
            "PG 17:1/18:1" : 760.5254,
            "PG 18:0/18:1" : 776.5567,
            "PG 18:1/18:1" : 774.5411,
            "PG 18:1/18:2" : 772.5254,
            "PG 18:1/20:1" : 802.5724,
            "PG 18:1/20:2" : 800.5567,
            "PG 18:1/20:3" : 798.5411,
            "PG 18:1/22:4" : 824.5567,
            "PI 16:0/20:3" : 860.5415,
            "PI 16:0/20:4" : 858.5258,
            "PI 18:0/20:3" : 888.5728,
            "PI 18:0/20:4" : 886.5571,
            "PI 18:1/18:1" : 862.5568,  # predicted
            "PI 18:1/20:3" : 886.5571,
            "PI 18:1/20:4" : 884.5415,
            "PS 16:0/16:0" : 735.5050,
            "PS 16:0/18:0" : 763.5363,
            "PS 16:0/18:1" : 761.5207,
            "PS 16:0/20:4" : 783.5050,
            "PS 16:1/18:0" : 761.5207,
            "PS 17:0/18:0" : 777.5520,
            "PS 18:0/18:0" : 791.5676,
            "PS 18:0/18:1" : 789.5520,
            "PS 18:0/18:2" : 787.5363,
            "PS 18:0/20:1" : 817.5833,
            "PS 18:0/20:3" : 813.5520,
            "PS 18:0/20:4" : 811.5363,
            "PS 18:0/22:4" : 839.5676,
            "PS 18:0/22:5" : 837.5521,  # predicted
            "PS 18:1/18:1" : 787.5363,
            "PS 18:1/20:4" : 809.5207,
            "PS 20:4/20:4" : 831.5050,
            "PS 22:6/22:6" : 879.5050,
            "SM (d18:0/13:0) (d18:1/12:0-OH)" : 662.5363,
            "SM (d18:0/16:0) (d18:1/15:0-OH)" : 704.5832,
            "SM (d18:0/17:0) (d18:1/16:0-OH)" : 718.5989,
            "SM (d18:0/19:0) (d18:1/18:0-OH)" : 746.6301,   # predicted
            "SM (d18:0/21:0) (d18:1/20:0-OH)" : 774.6613,   # predicted
            "SM (d18:0/22:0) (d18:1/21:0-OH)" : 788.6771,
            "SM (d18:0/23:0) (d18:1/22:0-OH)" : 802.6927,   # predicted
            "SM (d18:0/24:0) (d18:1/23:0-OH)" : 816.7084,
            "SM (d18:1/14:0) (d18:1/13:1-OH)" : 674.5363,
            "SM (d18:1/14:1) (d18:1/13:2-OH)" : 672.5208,   # predicted
            "SM (d18:1/15:0) (d18:1/14:1-OH)" : 688.5519,
            "SM (d18:1/16:0) (d18:1/15:1-OH)" : 702.5676,
            "SM (d18:1/16:1) (d18:1/15:2-OH)" : 700.5519,
            "SM (d18:1/17:0) (d18:1/16:1-OH)" : 716.5832,
            "SM (d18:1/18:0)" : 730.5989,
            "SM (d18:1/18:1)" : 728.5832,
            "SM (d18:1/18:2)" : 726.5677,   # predicted
            "SM (d18:1/21:0) (d18:1/20:1-OH)" : 772.6458,
            "SM (d18:1/22:0) (d18:1/21:1-OH)" : 786.6615,
            "SM (d18:1/23:0) (d18:1/22:1-OH)" : 800.6771,
            "SM (d18:1/23:1) (d18:1/22:2-OH)" : 798.6616,   # predicted
            "SM (d18:1/24:0) (d18:1/23:1-OH)" : 814.6928,
            "SM (d18:1/24:1) (d18:1/23:2-OH)" : 812.6771,
            "SM (d18:1/24:2)" : 810.6616,   # predicted
            "SM (d18:1/26:0) (d18:1/25:1-OH)" : 842.7241,
            "SM (d18:1/26:1) (d18:1/25:2-OH)" : 840.7084,
            "TAG 48:0 total (16:0/16:0/16:0)" : 807.29,  # glycerol - 3*H2O + 3*C16_0
            "TAG 48:1 total (14:0/16:0/18:1)(16:0/16:0/16:1)" : 805.2745,   # predicted
            "TAG 48:2 total (14:0/16:0/18:2)(14:0/16:1/18:1)(16:0/16:1/16:1)" : 803.259,    # predicted
            "TAG 48:3 total (16:1/16:1/16:1)" : 801.2435,   # predicted
            "TAG 50:0 total (16:0/16:0/18:0)" : 835.35, # glycerol - 3*H2O + 2*C16_0 + C18_0
            "TAG 50:1 total (14:0/18:0/18:1)(16:0/16:0/18:1)" : 833.3345,   # predicted
            "TAG 50:2 total (14:0/18:1/18:1)(16:0/16:0/18:2)(16:0/16:1/18:1)" : 831.319,    #predicted
            "TAG 50:3 total (14:0/18:1/18:2)(16:0/16:1/18:2)(16:1/16:1/18:1)" : 829.3035,   # predicted
            "TAG 50:4 total (16:1/16:1/18:2)" : 827.288,    # predicted
            "TAG 52:0 total (16:0/18:0/18:0)" : 863.41, #glycerol - 3*H2O + C16_0 + 2*C18_0
            "TAG 52:1 total (16:0/18:0/18:1)" : 861.3945,   # predicted
            "TAG 52:2 total (16:0/18:0/18:2)(16:0/18:1/18:1)(16:1/18:0/18:1)" : 859.379,    # predicted
            "TAG 52:3 total (16:0/18:1/18:2)(16:1/18:1/18:1)" : 857.3635,   # predicted
            "TAG 52:4 total (16:0/18:1/18:3)(16:0/18:2/18:2)(16:1/18:1/18:2)" : 855.348,    # predicted
            "TAG 52:5 total (16:0/18:2/18:3)(16:1/18:2/18:2)" : 853.3325,   # predicted
            "TAG 54:1 total (18:0/18:0/18:1)" : 889.4545,   # predicted
            "TAG 54:2 total (16:0/18:1/20:1)(18:0/18:1/18:1)" : 887.439,    # predicted
            "TAG 54:3 total (18:0/18:1/18:2)(18:1/18:1/18:1)" : 885.4235,   # predicted
            "TAG 54:4 total (18:1/18:1/18:2)" : 883.408,    # predicted
            "TAG 54:5 total (18:1/18:2/18:2)" : 881.3925,   # predicted
            "TAG 54:6 total (18:2/18:2/18:2)" : 879.377,    # predicted
            "TAG 56:3 total (18:1/18:1/20:1)" : 913.4599,   # glycerol - 3*H2O + 2*C18_1 + C20_1
            "TAG 56:7 total (16:0/18:1/22:6)" : 905.3978,   # predicted
            "TAG 56:8 total (16:0/18:2/22:6)" : 903.3823   # predicted
        }
        
        self.mass = lipid_weights.get(lipid, "Unrecognized lipid")

    # ========================
    def lipid_class(self):
    # ========================
        """
        Attribution of given lipid to a class.
        """
        string = self.name
        
        # when working with whole classes, there's no PC or PE followed by digits.
        if string in ("PC", "PE"):
            return string
    
        try:
            match = re.search(r"(CE|Cer|DAG|FC|Gb3|Glc|Lac|LPC|LPE|LPI|PA|PC \d|PC O|PC P|PE \d|PE P|PE O|PG|PI|PS|SM|TAG)", string)
            matched = match.group(1)
        
            # grouping PC \d and PE \d classes
            if re.search(r"PC \d", matched):
                return "PC"
            elif re.search(r"PE \d", matched):
                return "PE"
            else:
                return matched
        
        except:
            fatty_acids = [
                "17:1", "13:2", "14:1", "15:0", "20:3", "13:1",
                "20:1", "19:0", "17:0", "22:2", "19:1", "14:0",
                "22:1", "21:1", "16:1", "18:1", "18:2", "20:4",
                "18:0", "22:3", "22:4", "23:0", "16:0", "22:6",
                "19:2", "23:1", "25:1", "18:3", "12:0", "22:5",
                "20:5", "20:2", "20:0", "24:1", "15:1", "23:2",
                "24:0", "26:0", "15:2", "26:1", "21:0", "22:0",
                "24:5", "25:2", "24:2"
            ]
            if string in fatty_acids:
                return string
            else:
                return "unknown lipid type"
    
    # ===========================================
    def residues(self, *, drop_ambiguous = False):
    # ===========================================
        """
        Reads something like PG 18:1/22:4
        Spits out: a <tuple> something like (["18:1", "22:4"], 1)
        Where <int> n hints at how many times each residue should be divided
        before taking note of each residue
        (this is for counting ambiguous lipids).
    
        Options
        =======
        If drop_ambiguous is set to True, ambiguous lipids are not eveluated
        and an empty <list> is returned instead (and also a 0 dividend).
        For example: TAG 52:4 total (16:0/18:1/18:3)(16:0/18:2/18:2)
        returns (["16:0", "18":"1", "18":"3", "16:0", "18":"2", "18":"2"], 2)
        if drop_ambiguous = False; else it returns ([], 0)
        """
        lipid = self.name
        lipid_residue = r"([\d]*:[\d]*)"  # what we want
        match = re.findall(lipid_residue, lipid)
    
        lipid_class = self.lipid_class()
    
        ambiguity = False # lipids are assumed unambiguous unless proven otherwise
    
        def pack_returntuple():
            if drop_ambiguous == True and ambiguity == True:
                return ([], 0)
            else:
                if lipid_class == "TAG":
                    return (match[1:], dividend)
                else:
                    return (match, dividend)
    
        if lipid_class == "TAG":
            if len(lipid) == 31:
                dividend = 1
                return pack_returntuple()
            elif 31 < len(lipid) < 63:
                ambiguity = True
                dividend = 2
                return pack_returntuple()
            else:
                ambiguity = True
                dividend = 3
                return pack_returntuple()
        elif lipid_class == "SM":
            if len(lipid) == 15:
                dividend = 1
                return pack_returntuple()
            else:
                ambiguity = True
                dividend = 2
                return pack_returntuple()
        elif lipid_class == "PC P":
            ambiguity = True
            dividend = 2
            return pack_returntuple()
        elif lipid_class == "PE P":
            if len(lipid) == 14:
                dividend = 1
                return pack_returntuple()
            else:
                ambiguity = True
                dividend = 2
                return pack_returntuple()
        else:
            dividend = 1
            return pack_returntuple()


def make_residues_table(dataframe, *, drop_ambiguous=False, name="residues_table",
                        replace_nan=0, cleanup=True, absolute_amount=False,
                        unwanted=["total", "fc", "tc"], **kwargs):

    """ takes a pandas DataFrame as input, and outputs a pandas DataFrame that
    contains individual residues as index, and their amount for every sample/column.

    Parameters
    ==========

    dataframe: a pandas dataframe of data. Lipid names as index, and samples as columns
        (just unlike sklearn wants it, but as you might get it from Tableau software
        tables. Just dataframe.T your table - that would just do the trick).

    drop_ambiguous: <bool> don't take isobars into consideration. Defaults to False. If True,
        each residue is divided by its uncertainty.

    name: <str> a tag that gets attached to the returned dataframe, so you can use it
        to save it afterwards. The tag is found in the .name attribute.

    replace_nan: <object> the object you would like to replace your missing values with.
        It can be set to False, but I would suggest against what.

    cleanup: <bool> Whether to perform a cleanup of unwanted lipids that can be present
        in the index. Unwanted strings are read from the 'unwanted' parameter. Defaults
        to True
    
    absolute_amount <bool> Wheter to count the individual number of residues, rather to
        sticking to the same units found in the original table. Defaults to False

    unwanted: <list> <set> <tuple> Strings that must be removed from the lipid index.
        Defaults to ["total", "fc", "tc"]

    returns:
    ========

    pandas DataFrame

    """
    
    if not isinstance(dataframe, pd.core.frame.DataFrame):
        raise TypeError("Input table must be a pandas DataFrame")

    df = dataframe._get_numeric_data().copy()

    if replace_nan:
        df = df.replace(to_replace=np.nan, value=replace_nan)

    if cleanup:
        # removing some common unwanted lipids from the index,
        # like total counts and free cholesterol
        newindex = [x for x in df.index if not x.lower() in unwanted]
        df = df.reindex(newindex)


    # beginning to transform data
    # ===========================
    masterlist = []

    # We will slice the dataframe, creating a serie for each sample (c). The index
    # of the series will hold the individual lipids (like "CE 16:1"), joined
    # to their values in that sample (c). Yes, I could have called it sample.
    # For each sample, a dictionary (mm) will be constructed with simple fatty
    # acids as keys, and their amount as values.

    for c in df.columns:
        serie = df[c]

        # each column now becomes a series with the index == df.index, and the
        # corresponding values

        mm = {}
        for lipid, picomoles in zip(serie.index, serie):

            try:
                lip = Lipid(lipid, amount = picomoles, **kwargs)
            except TypeError:
                print(f"Something went wrong with lipid {lipid} and amount {picomoles}")
                raise

            residues, coeff = lip.residues(drop_ambiguous=drop_ambiguous)
            if len(residues) < 0:
                continue
            else:
                for residue in residues:
                    mm.setdefault(residue, 0)

                    if absolute_amount:
                        mm[residue] += lip.molecules / coeff # exact number of residues
                    else:
                        mm[residue] += lip.amount / coeff # units of residue 

        masterlist.append(mm)


    dfinal = pd.concat([pd.DataFrame.from_dict(x, orient="index") for x in masterlist], axis=1)
    dfinal.columns = df.columns

    dfinal.name = name

    return dfinal
