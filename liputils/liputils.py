# -*- coding: utf-8 -*-
# Author: Manzini Stefano; stefano.manzini@gmail.com
# 040320

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

    Returns:
    ========

    Lipid object
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

        
        self.amount = amount    # in the specified unit

        self.molecules = amount * unit_conversion.get(unit, 0)  # number of molecules
        
        self.name = lipid
        
        # self.mass argument was removed in 0.13. I'd prefer to stick with consistent
        # features that are actually used. Ready to bring it back for all RefMet
        # fatty-acids contaning compounds if requested by the users.


    def lipid_class(self):
        """
        Attribution of given lipid to a class. This function tries to extract lipid
        identifiers that fall outside RefMet official naming scheme.

        It would be best to translate your particular naming scheme with the aid of
        RefMet's online translator, which is available at:

        https://www.metabolomicsworkbench.org/databases/refmet/name_to_refmet_form.php

        If this suits your naming scheme but fails to get lipids right, report the issue.
        """
        stringlike = self.name
        
        # when working with whole classes, there's no PC or PE followed by digits.
        if stringlike in ("PC", "PE"):
            return stringlike
    
        try:
            match = re.search(r"(CE|Cer|DAG|FC|Gb3|Glc|Lac|LPC|LPE|LPI|PA|PC \d|PC O|PC P|PE \d|PE P|PE O|PG|PI|PS|SM|TAG)", stringlike)
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
            if stringlike in fatty_acids:
                return stringlike
            else:
                return "unknown lipid type"


    def refmet_class(self):
        """
        Attribution of given lipid to a class. This function extracts the lipid
        class from the official RefMet nomenclature.
        """

        stringlike = self.name

        # 
        if stringlike.startswith("1,"):
            return "DG"
    
        coord = stringlike.find("(")
        if coord != -1:
            refclass = stringlike[:coord]
        else:
            refclass = stringlike

        if len(refclass) == 0:
            return "unknown lipid type"

        # let's give a chance to lipids with annotations with a whitespace
        if refclass.endswith(" "):
            refclass = refclass[:-1]
            if len(refclass) ==  0:
                return "unknown lipid type"

        return refclass


    def refmet_fullclass(self):
        """
        Attribution of given lipid to a class. This function extracts the lipid
        class (full name) from the official RefMet nomenclature.
        """

        stringlike = self.name

        lipids_book = {
            "CAR": "Fatty esters",
            "CoA": "Fatty esters",
            "FAHFA": "Fatty esters",
            "DG": "Diradylglycerols",
            "DGDG": "Glycosyldiradylglycerols",
            "MGDG": "Glycosyldiradylglycerols",
            "MG": "Monoradylglycerols",
            "MeDAG": "Triradylglycerols",
            "TG": "Triradylglycerols",
            "CL": "Cardiolipins",
            "CDP-DG": "CDP-Glycerols",
            "LPA": "Glycerophosphates",
            "PA": "Glycerophosphates",
            "LPC": "Glycerophosphocholines",
            "PC": "Glycerophosphocholines",
            "LPE": "Glycerophosphoethanolamines",
            "PE": "Glycerophosphoethanolamines",
            "LPG": "Glycerophosphoglycerols",
            "LPGP": "Glycerophosphoglycerols",
            "PG": "Glycerophosphoglycerols",
            "PIP2": "Glycerophosphoinositol phosphates",
            "PIP3": "Glycerophosphoinositol phosphates",
            "LPI": "Glycerophosphoinositols",
            "LPIP": "Glycerophosphoinositols",
            "LPIP2": "Glycerophosphoinositols",
            "LPIP3": "Glycerophosphoinositols",
            "PI": "Glycerophosphoinositols",
            "LPS": "Glycerophosphoserines",
            "PS": "Glycerophosphoserines",
            "CPA": "Other Glycerophospholipids",
            "1-DeoxyCer": "Ceramides",
            "Cer": "Ceramides",
            "CerP": "Ceramides",
            "PI-Cer": "Ceramides",
            "GD1": "Glycosphingolipids",
            "GD2": "Glycosphingolipids",
            "GD3": "Glycosphingolipids",
            "GM1": "Glycosphingolipids",
            "GM2": "Glycosphingolipids",
            "GM3": "Glycosphingolipids",
            "GM4": "Glycosphingolipids",
            "GT1a": "Glycosphingolipids",
            "GT1b": "Glycosphingolipids",
            "GT2": "Glycosphingolipids",
            "GT3": "Glycosphingolipids",
            "GlcAbeta-Cer": "Glycosphingolipids",
            "GalCer": "Glycosphingolipids",
            "GlcCer": "Glycosphingolipids",
            "HexCer": "Glycosphingolipids",
            "Hex2Cer": "Glycosphingolipids",
            "LacCer": "Glycosphingolipids",
            "Sulfatide": "Glycosphingolipids",
            "PE-Cer": "Phosphosphingolipids",
            "unknown lipid type": "unknown lipid type",
            "Iso": "Sphingoid bases",
            "SM": "Sphingomyelins",
            "Campesterol ester": "Sterol esters",
            "CE": "Sterol esters",
            "Sitosterol Ester": "Sterol esters",
            "Stigmasterol Ester": "Sterol esters",
            "FA": "Fatty acid", # example: FA(27:1), Heptacosadienoic acid
        }

        return lipids_book.get(self.refmet_class(), "unknown lipid type")



    def residues(self, *, drop_ambiguous=False):
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
        
        # TODO: we could regex the dividend int instead of guessing from
        # string length.
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


    def refmet_residues(self, *, drop_ambiguous=False):
        """
        Reads something like PE(22:6(4Z,7Z,10Z,13Z,16Z,19Z)_22:6(4Z,7Z,10Z,13Z,16Z,19Z))
        Spits out: a <tuple> something like (['22:6', '22:6'], 1)
        Where <int> n hints at how many times each residue should be divided
        before taking note of each residue
        (this is for counting ambiguous lipids).
    
        Options
        =======
        If drop_ambiguous is set to True, ambiguous lipids are not eveluated
        and an empty <list> is returned instead (and also a 0 dividend).
        For example: PG(P-22:1/18:1)/PG(O-22:2/18:1)
        returns (['22:1', '18:1', '22:2', '18:1'], 2)
        if drop_ambiguous = False; else it returns ([], 0)
        """
        lipid = self.name


        # == preprocessing ==

        def preprocess(stringlike):
            """ Taming some troublesome lipid names.
            """
            stringlike = stringlike.replace("(OH)", "")
            stringlike = stringlike.replace("(9E)", "") # nasty DG(18:1(9E)/0:0/0:0/18:1(9E))
                                                        # trying to save some regex time here

            #.replace("(3OH)", "") not dangerous
            #.replace("(2OH)", "") absent
            # (*Me*) series: not dangerous TODO: verify
            #lipid = lipid.replace("(Me)", "").replace("(2Me)", "")\
            #.replace("(3Me)", "").replace("(Me2)", "")\
            #.replace("(2Me,6Me)", "") taken care of
            # lipid = lipid.replace("(2COOH)", "") # not dangerous

            if "," in stringlike:
                if stringlike.startswith("1,"):
                    return stringlike
                else:
                    stringlike = stringlike.replace(",", "")
                    pattern = r"\(.[^:]*?\)"
                    match = re.findall(pattern, stringlike)
                    if len(match) == 0:
                        return stringlike
                    else:
                        for m in match:
                            stringlike = stringlike.replace(m, "")
    
            return stringlike


        def sloppy_lipid_identifier(stringlike, fatty_amides=True, fatty_alcohols=True):
            """Not very sophisticated. Get things done. ~90,000+ queries/sec on old i5-3427U
            for a failed search (must go through all of them)
            
            Special, alkinyl-fatty acids (with C≡C insaturations) are NOT decoded.
            l
            Other compounds could get decoded, even if decoding them is **not the primary
            goal** of this function. To allow to try and decode other fatty compounds,
            set target molecule to True (default) . It is possible to exclude recognition
            for these molecules by setting the appropriate class to False.
            """
            stringlike = stringlike.lower()
            
            if fatty_amides is False:
                if stringlike.endswith("mide"):
                    return ([], 0)
                
            if fatty_alcohols is False:
                if stringlike.endswith("ol"):
                    return ([], 0)
            
            if "acetic acid" in stringlike:
                return (["2:0"], 1)
            if "acet" in stringlike: # no need to worry: acetates catched with "-ate" catcher
                return (["2:0"], 1)
            elif "propion" in stringlike:
                return (["3:0"], 1)
            elif "propan" in stringlike:
                return (["3:0"], 1)
            elif "propen" in stringlike:
                return (["3:1"], 1)    
            elif "acrylic" in stringlike:
                return (["3:1"], 1)
            elif "butyr" in stringlike:
                return (["4:0"], 1)
            elif "butan" in stringlike:
                return (["4:0"], 1)
            elif "buten" in stringlike:
                return (["4:1"], 1)
            elif "pentan" in stringlike:
                return (["5:0"], 1)
            elif "penten" in stringlike:
                return (["5:1"], 1)
            elif "valer" in stringlike: # valeric acid
                return (["5:0"], 1)
            elif "levulin" in stringlike: # levulinic acid (KETO-acid)
                return (["5:0"], 1)
            elif "hexan" in stringlike:
                return (["6:0"], 1)
            elif "hexen" in stringlike:
                return (["6:1"], 1)
            elif "hexadien" in stringlike:
                return (["6:2"], 1)
            elif "heptan" in stringlike:
                return (["7:0"], 1)
            elif "enanthi" in stringlike: # enanthic acid
                return (["7:0"], 1)
            elif "heptylic" in stringlike: # enanthic acid
                return (["7:0"], 1)
            elif "hepten" in stringlike:
                return (["7:1"], 1)
            elif "heptadien" in stringlike:
                return (["7:2"], 1)
            elif "heptatrien" in stringlike:
                return (["7:3"], 1)
            elif "nonan" in stringlike:
                return (["9:0"], 1)
            elif "nonen" in stringlike:
                return (["9:1"], 1)
            elif "nonadien" in stringlike:
                return (["9:2"], 1)
            elif "pelargo" in stringlike: # pelargonic acid
                return (["9:0"], 1)
            elif "capric" in stringlike:
                return (["10:0"], 1)
            elif "caproyl" in stringlike:
                return (["10:0"], 1)
            elif "caproleic" in stringlike:
                return (["10:1"], 1)
            elif "caproleyl" in stringlike:
                return (["10:1"], 1)
            elif "undecan" in stringlike:
                return (["11:0"], 1)
            elif "undecyl" in stringlike:
                return (["11:0"], 1)
            elif "undecen" in stringlike:
                return (["11:1"], 1)
            elif "undecadien" in stringlike:
                return (["11:2"], 1)
            elif "undecatrien" in stringlike:
                return (["11:3"], 1)
            elif "dodecan" in stringlike:
                return (["12:0"], 1)
            elif "lauric" in stringlike:
                return (["12:0"], 1)
            elif "lauryl" in stringlike: # lauric acid
                return (["12:0"], 1)
            elif "lauroyl" in stringlike: # lauric acid
                return (["12:0"], 1)
            elif "dodecen" in stringlike:
                return (["12:1"], 1)
            elif "lauroleyc" in stringlike:
                return (["12:1"], 1)
            elif "lauroleyl" in stringlike:
                return (["12:1"], 1)
            elif "dodecadien" in stringlike:
                return (["12:2"], 1)    
            elif "dodecatrien" in stringlike:
                return (["12:3"], 1)
            elif "dodecatetraen" in stringlike:
                return (["12:4"], 1)
            elif "dodecapentaen" in stringlike:
                return (["12:5"], 1)
            elif "dodecapenten" in stringlike:
                return (["12:5"], 1)
            elif "tridecan" in stringlike:
                return (["13:0"], 1)
            elif "tridecylic" in stringlike:
                return (["13:0"], 1)
            elif "tridecen" in stringlike:
                return (["13:1"], 1)
            elif "tridecadien" in stringlike:
                return (["13:2"], 1)
            elif "tridecatrien" in stringlike:
                return (["13:3"], 1)
            elif "tetradecan" in stringlike:
                return (["14:0"], 1)
            elif "myristic" in stringlike:
                return (["14:0"], 1)
            elif "myristoyl" in stringlike:
                return (["14:0"], 1)
            elif "tetradecen" in stringlike:
                return (["14:1"], 1)
            elif "myristoleic" in stringlike:
                return (["14:1"], 1)
            elif "myristoleyl" in stringlike:
                return (["14:1"], 1)
            elif "tetradecadien" in stringlike:
                return (["14:2"], 1)
            elif "tetradecatrien" in stringlike:
                return (["14:3"], 1)
            elif "tetradecatetraen" in stringlike:
                return (["14:4"], 1)
            elif "tetradecapentaen" in stringlike:
                return (["14:5"], 1)
            elif "tetradecapenten" in stringlike:
                return (["14:5"], 1)
            elif "hexadecan" in stringlike:
                return (["16:0"], 1)
            elif "hexadecan" in stringlike:
                return (["16:1"], 1)
            elif "pentadecan" in stringlike:
                return (["15:0"], 1)
            elif "pentadecyl" in stringlike:
                return (["15:0"], 1)
            elif "pentadecen" in stringlike:
                return (["15:1"], 1)
            elif "pentadecadien" in stringlike:
                return (["15:2"], 1)
            elif "pentadecatrien" in stringlike:
                return (["15:3"], 1)
            elif "pentadecatetraen" in stringlike:
                return (["15:4"], 1)
            elif "pentadecapentaen" in stringlike:
                return (["15:5"], 1)
            elif "pentadecapenten" in stringlike:
                return (["15:5"], 1)
            elif "hexadecan" in stringlike:
                return (["16:0"], 1)    
            elif "palmit" in stringlike: # palmitic acid
                return (["16:0"], 1)
            elif "hexadecen" in stringlike:
                return (["16:1"], 1)
            elif "hexadecadien" in stringlike:
                return (["16:2"], 1)
            elif "hexadecatrien" in stringlike:
                return (["16:3"], 1)
            elif "hexadecatetraen" in stringlike:
                return (["16:4"], 1)
            elif "hexadecapentaen" in stringlike:
                return (["16:5"], 1)
            elif "hexadecapenten" in stringlike:
                return (["16:5"], 1)
            elif "hexadecahexaen" in stringlike:
                return (["16:6"], 1)
            elif "hexadecahexen" in stringlike:
                return (["16:6"], 1)
            elif "hexadecaheptaen" in stringlike:
                return (["16:7"], 1)
            elif "hexadecahepten" in stringlike:
                return (["16:7"], 1)
            elif "heptadecan" in stringlike:
                return (["17:0"], 1)
            elif "margar" in stringlike: # margaric acid
                return (["17:0"], 1)
            elif "heptadecen" in stringlike:
                return (["17:1"], 1)
            elif "heptadecadien" in stringlike:
                return (["17:2"], 1)
            elif "heptadecatrien" in stringlike:
                return (["17:3"], 1)
            elif "heptadecatetraen" in stringlike:
                return (["17:4"], 1)
            elif "heptadecapentaen" in stringlike:
                return (["17:5"], 1)
            elif "heptadecapenten" in stringlike:
                return (["17:5"], 1)
            elif "heptadecahexaen" in stringlike:
                return (["17:6"], 1)
            elif "heptadecahexen" in stringlike:
                return (["17:6"], 1)
            elif "heptadecaheptaen" in stringlike:
                return (["17:7"], 1)
            elif "heptadecahepten" in stringlike:
                return (["17:7"], 1)
            elif "octadecan" in stringlike:
                return (["18:0"], 1)
            elif "octadecen" in stringlike:
                return (["18:1"], 1)
            elif "vaccen" in stringlike: # vaccenic acid
                return (["18:1"], 1)
            elif "octadecadien" in stringlike:
                return (["18:2"], 1)
            elif "linolelaid" in stringlike: # linolelaidic acid
                return (["18:2"], 1)
            elif "octadecatrien" in stringlike:
                return (["18:3"], 1)
            elif "linolen" in stringlike: # linolenic acid
                return (["18:3"], 1)
            elif "columbin" in stringlike: # columbinic acid (5E,9E,12E)-octadeca-5,9,12-trienoic acid
                return (["18:3"], 1)
            elif "stearidon" in stringlike: # stearidoic acid (6Z,9Z,12Z,15Z)-octadeca-6,9,12,15-tetraenoic acid
                return (["18:3"], 1)
            elif "octadecatrien" in stringlike: # undifferentiated
                return (["18:3"], 1)
            elif "parinar" in stringlike: # parinaric acid
                return (["18:4"], 1)
            elif "octadecatetraen" in stringlike: # undifferentiated
                return (["18:4"], 1)
            elif "octadecapentaen" in stringlike: # undifferentiated
                return (["18:5"], 1)
            elif "octadecapenten" in stringlike: # undifferentiated
                return (["18:5"], 1)
            elif "octadecahexaen" in stringlike: # undifferentiated
                return (["18:6"], 1)
            elif "octadecahexen" in stringlike: # undifferentiated
                return (["18:6"], 1)
            elif "octadecaheptaen" in stringlike: # undifferentiated
                return (["18:7"], 1)
            elif "octadecahepten" in stringlike: # undifferentiated
                return (["18:7"], 1)
            elif "nonadecan" in stringlike:
                return (["19:0"], 1)
            elif "nonadecylic" in stringlike:
                return (["19:0"], 1)
            elif "nonadecen" in stringlike:
                return (["19:1"], 1)
            elif "nonadecadien" in stringlike:
                return (["19:2"], 1)
            elif "nonadecatrien" in stringlike:
                return (["19:3"], 1)
            elif "nonadecatetraen" in stringlike:
                return (["19:4"], 1)
            elif "nonadecapentaen" in stringlike:
                return (["19:5"], 1)
            elif "nonadecapenten" in stringlike:
                return (["19:5"], 1)
            elif "nonadecahexaen" in stringlike:
                return (["19:6"], 1)
            elif "nonadecahexen" in stringlike:
                return (["19:6"], 1)
            elif "nonadecaheptaen" in stringlike:
                return (["19:7"], 1)
            elif "nonadecahepten" in stringlike:
                return (["19:7"], 1)
            elif "arachidic" in stringlike:
                return (["20:0"], 1)
            elif "arachida" in stringlike:
                return (["20:0"], 1)
            elif "arachidyl" in stringlike:
                return (["20:0"], 1)
            elif "eicosenoic" in stringlike:
                return (["20:1"], 1)
            elif "eicosenyl" in stringlike:
                return (["20:1"], 1)
            elif "gadoleic" in stringlike:
                return (["20:1"], 1)
            elif "gadoleyl" in stringlike:
                return (["20:1"], 1)
            elif "eicosadienoic" in stringlike:
                return (["20:2"], 1)
            elif "eicosadienoyl" in stringlike:
                return (["20:2"], 1)
            elif "eicosatrienoic" in stringlike:
                return (["20:3"], 1)
            elif "eicosatrienoyl" in stringlike:
                return (["20:3"], 1)
            elif "mead acid" in stringlike:
                return (["20:3"], 1)
            elif "eicosatetraen" in stringlike:
                return (["20:4"], 1)
            elif "arachidon" in stringlike:
                return (["20:4"], 1)
            elif "eicosapentenoic" in stringlike:
                return (["20:5"], 1)
            elif "eicosapentenoyl" in stringlike:
                return (["20:5"], 1)
            elif "eicosapentaenoic" in stringlike:
                return (["20:5"], 1)
            elif "eicosapentaenoyl" in stringlike:
                return (["20:5"], 1)
            elif "timnodon" in stringlike: # timnodonic acid
                return (["20:5"], 1)
            elif stringlike == "epa":
                return (["20:5"], 1)
            elif "heneicosan" in stringlike:
                return (["21:0"], 1)
            elif "heneicosen" in stringlike:
                return (["21:1"], 1)    
            elif "heneicosaen" in stringlike:
                return (["21:1"], 1)
            elif "heneicosadien" in stringlike:
                return (["21:2"], 1)
            elif "heneicosatrien" in stringlike:
                return (["21:3"], 1)
            elif "heneicosatetraen" in stringlike:
                return (["21:4"], 1)
            elif "heneicosapenten" in stringlike:
                return (["21:5"], 1)
            elif "heneicosapentaen" in stringlike:
                return (["21:5"], 1)
            elif "heneicosahexen" in stringlike: # heneicosahexaenoic acid
                return (["21:6"], 1)
            elif "heneicosahexaen" in stringlike: # heneicosahexaenoic acid
                return (["21:6"], 1)
            elif "docosanoic" in stringlike:
                return (["22:1"], 1)
            elif "docosanoyl" in stringlike:
                return (["22:1"], 1)
            elif "behen" in stringlike:
                return (["22:0"], 1)
            elif "erucic" in stringlike: # (Z)-docos-13-enoic acid
                return (["22:1"], 1)
            elif "erucyl" in stringlike: # (Z)-docos-13-enoic acid
                return (["22:1"], 1)
            elif "brassidic" in stringlike: # (E)-docos-13-enoic acid
                return (["22:1"], 1)
            elif "brassidoyl" in stringlike: # (E)-docos-13-enoic acid
                return (["22:1"], 1)
            elif "docosenoic" in stringlike:
                return (["22:1"], 1)
            elif "docosenoyl" in stringlike:
                return (["22:1"], 1)
            elif "docosadieno" in stringlike:
                return (["22:2"], 1)
            elif "brassic" in stringlike: #brassic acid
                return (["22:2"], 1)
            elif "docosatrien" in stringlike:
                return (["22:3"], 1)
            elif "docosatetraenoic" in stringlike:
                return (["22:4"], 1)
            elif "docosatetraenoyl" in stringlike:
                return (["22:4"], 1)
            elif "docosapentenoic" in stringlike:
                return (["22:5"], 1)
            elif "docosapentenoyl" in stringlike:
                return (["22:5"], 1)
            elif "docosapentaenoic" in stringlike:
                return (["22:5"], 1)
            elif "docosapentaenoyl" in stringlike:
                return (["22:5"], 1)
            elif "clupanodon" in stringlike: #clupanodonic acid
                return (["22:5"], 1)
            elif stringlike == "dpa":
                return (["22:5"], 1)
            elif "docosahexen" in stringlike:
                return (["22:6"], 1)
            elif "docosahexaen" in stringlike:
                return (["22:6"], 1)
            elif stringlike == "dha":
                return (["22:6"], 1)
            elif "tricosan" in stringlike:
                return (["23:0"], 1)
            elif "tricosen" in stringlike:
                return (["23:1"], 1)
            elif "tricosadien" in stringlike:
                return (["23:2"], 1)
            elif "tricosatrien" in stringlike:
                return (["23:3"], 1)
            elif "tricosatetraen" in stringlike:
                return (["23:4"], 1)
            elif "tricosapentaen" in stringlike:
                return (["23:5"], 1)
            elif "tricosapenten" in stringlike:
                return (["23:5"], 1)
            elif "tricosahexaen" in stringlike:
                return (["23:6"], 1)
            elif "tricosahexen" in stringlike:
                return (["23:6"], 1)
            elif "tricosaheptaen" in stringlike:
                return (["23:7"], 1)
            elif "tricosahepten" in stringlike:
                return (["23:7"], 1)
            elif "tetracosanoic" in stringlike:
                return (["24:0"], 1)
            elif "tetracosanoyl" in stringlike:
                return (["24:0"], 1)
            elif "lignocer" in stringlike:
                return (["24:0"], 1)
            elif "tetracosanoic" in stringlike:
                return (["24:0"], 1)
            elif "tetracosanoyl" in stringlike:
                return (["24:0"], 1)
            elif "tetracosen" in stringlike:
                return (["24:1"], 1)
            elif "nervon" in stringlike:
                return (["24:1"], 1)
            elif "tetracosadien" in stringlike:
                return (["24:2"], 1)
            elif "tetracosatrien" in stringlike:
                return (["24:3"], 1)
            elif "tetracosatetraen" in stringlike:
                return (["24:4"], 1)
            elif "tetracosapentaen" in stringlike:
                return (["24:5"], 1)
            elif "tetracosapenten" in stringlike:
                return (["24:5"], 1)
            elif "tetracosahexaen" in stringlike:
                return (["24:6"], 1)
            elif "tetracosahexen" in stringlike:
                return (["24:6"], 1)
            elif stringlike == "tha":
                return (["24:6"], 1)
            elif "tetracosaheptaen" in stringlike:
                return (["24:7"], 1)
            elif "tetracosahepten" in stringlike:
                return (["24:7"], 1)
            elif "pentacosan" in stringlike:
                return (["25:0"], 1)
            elif "pentacosen" in stringlike:
                return (["25:1"], 1)
            elif "pentacosadien" in stringlike:
                return (["25:2"], 1)
            elif "pentacosatrien" in stringlike:
                return (["25:3"], 1)
            elif "pentacosatetraen" in stringlike:
                return (["25:4"], 1)
            elif "pentacosapentaen" in stringlike:
                return (["25:5"], 1)
            elif "pentacosapenten" in stringlike:
                return (["25:5"], 1)
            elif "pentacosahexaen" in stringlike:
                return (["25:6"], 1)
            elif "pentacosahexen" in stringlike:
                return (["25:6"], 1)
            elif "pentacosaheptaen" in stringlike:
                return (["25:7"], 1)
            elif "pentacosahepten" in stringlike:
                return (["25:7"], 1)
            elif "hexacosan" in stringlike:
                return (["26:0"], 1)
            elif "cerot" in stringlike: # cerotic acid
                return (["26:0"], 1)
            elif "hexacosen" in stringlike:
                return (["26:1"], 1)
            elif "hexacosadien" in stringlike:
                return (["26:2"], 1)
            elif "hexacosatrien" in stringlike:
                return (["26:3"], 1)
            elif "hexacosatetraen" in stringlike:
                return (["26:4"], 1)
            elif "hexacosapentaen" in stringlike:
                return (["26:5"], 1)
            elif "hexacosapenten" in stringlike:
                return (["26:5"], 1)
            elif "hexacosahexaen" in stringlike:
                return (["26:6"], 1)
            elif "hexacosahexen" in stringlike:
                return (["26:6"], 1)
            elif "hexacosaheptaen" in stringlike:
                return (["26:7"], 1)
            elif "hexacosahepten" in stringlike:
                return (["26:7"], 1)
            elif "heptacosanoic" in stringlike:
                return (["27:0"], 1)
            elif "heptacosanoyl" in stringlike:
                return (["27:0"], 1)
            elif "carboceric" in stringlike:
                return (["27:0"], 1)
            elif "heptacosen" in stringlike:
                return (["27:1"], 1)
            elif "heptacosadien" in stringlike:
                return (["27:2"], 1)
            elif "heptacosatrien" in stringlike:
                return (["27:3"], 1)
            elif "heptacosatetraen" in stringlike:
                return (["27:4"], 1)
            elif "heptacosatetren" in stringlike:
                return (["27:4"], 1)
            elif "heptacosapentaen" in stringlike:
                return (["27:5"], 1)
            elif "heptacosapenten" in stringlike:
                return (["27:5"], 1)
            elif "heptacosahexaen" in stringlike:
                return (["27:6"], 1)
            elif "heptacosahexten" in stringlike:
                return (["27:6"], 1)
            elif "heptacosaheptaen" in stringlike:
                return (["27:7"], 1)
            elif "heptacosahepten" in stringlike:
                return (["27:7"], 1)
            elif "octacosan" in stringlike:
                return (["28:0"], 1)
            elif "montan" in stringlike: # montanic acid
                return (["28:0"], 1)
            elif "octacosen" in stringlike:
                return (["28:1"], 1)
            elif "octacosadien" in stringlike:
                return (["28:2"], 1)
            elif "octacosatrien" in stringlike:
                return (["28:3"], 1)
            elif "octacosatetraen" in stringlike:
                return (["28:4"], 1)
            elif "octacosapentaen" in stringlike:
                return (["28:5"], 1)
            elif "octacosapenten" in stringlike:
                return (["28:5"], 1)
            elif "octacosahexaen" in stringlike:
                return (["28:6"], 1)
            elif "octacosahexen" in stringlike:
                return (["28:6"], 1)
            elif "octacosaheptaen" in stringlike:
                return (["28:7"], 1)
            elif "octacosahepten" in stringlike:
                return (["28:7"], 1)
            elif "octacosaoctaen" in stringlike:
                return (["28:8"], 1)
            elif "nonacosan" in stringlike:
                return (["29:0"], 1)
            elif "nonacosen" in stringlike:
                return (["29:1"], 1)
            elif "nonacosadien" in stringlike:
                return (["29:2"], 1)
            elif "nonacosatrien" in stringlike:
                return (["29:3"], 1)
            elif "nonacosatetraen" in stringlike:
                return (["29:4"], 1)
            elif "nonacosapentaen" in stringlike:
                return (["29:5"], 1)
            elif "nonacosapenten" in stringlike:
                return (["29:5"], 1)
            elif "nonacosahexaen" in stringlike:
                return (["29:6"], 1)
            elif "nonacosahexen" in stringlike:
                return (["29:6"], 1)
            elif "nonacosaheptaen" in stringlike:
                return (["29:7"], 1)
            elif "nonacosahepten" in stringlike:
                return (["29:7"], 1)
            elif "triacontan" in stringlike:
                return (["30:0"], 1)
            elif "melissic" in stringlike:
                return (["30:0"], 1)
            elif "triaconten" in stringlike:
                return (["30:1"], 1)
            elif "triacontadien" in stringlike:
                return (["30:2"], 1)
            elif "triacontatrien" in stringlike:
                return (["30:3"], 1)
            elif "triacontatetraen" in stringlike:
                return (["30:4"], 1)
            elif "triacontapentaen" in stringlike:
                return (["30:5"], 1)
            elif "triacontapenten" in stringlike:
                return (["30:5"], 1)
            elif "triacontahexaen" in stringlike:
                return (["30:6"], 1)
            elif "triacontahexen" in stringlike:
                return (["30:6"], 1)
            elif "triacontaheptaen" in stringlike:
                return (["30:7"], 1)
            elif "triacontahepten" in stringlike:
                return (["30:7"], 1)
            elif "hentriacontan" in stringlike:
                return (["31:0"], 1)
            elif "hentriaconten" in stringlike:
                return (["31:1"], 1)
            elif "hentriacontaen" in stringlike:
                return (["31:1"], 1)
            elif "hentriacontadien" in stringlike: # hentriacontadienoic acid
                return (["31:2"], 1)
            elif "hentriacontatrien" in stringlike: # hentriacontadienoic acid
                return (["31:3"], 1)
            elif "hentriacontatetraen" in stringlike:
                return (["31:4"], 1)
            elif "hentriacontapentaen" in stringlike:
                return (["31:5"], 1)
            elif "hentriacontapenten" in stringlike:
                return (["31:5"], 1)
            elif "hentriacontahexaen" in stringlike:
                return (["31:6"], 1)
            elif "hentriacontahexen" in stringlike:
                return (["31:6"], 1)
            elif "hentriacontaheptaen" in stringlike:
                return (["31:7"], 1)
            elif "hentriacontahepten" in stringlike:
                return (["31:7"], 1)
            elif "lacceroic" in stringlike: # dotriacontanoic acid
                return (["32:0"], 1)
            elif "dotriacontan" in stringlike: # dotriacontanoic acid
                return (["32:0"], 1)
            elif "dotriaconten" in stringlike: # dotriacontanoic acid
                return (["32:1"], 1)
            elif "dotriacontadien" in stringlike:
                return (["32:2"], 1)
            elif "dotriacontatrien" in stringlike:
                return (["32:3"], 1)
            elif "dotriacontatetraen" in stringlike:
                return (["32:4"], 1)
            elif "dotriacontapentaen" in stringlike:
                return (["32:5"], 1)
            elif "dotriacontapenten" in stringlike:
                return (["32:5"], 1)
            elif "dotriacontahexaen" in stringlike:
                return (["32:6"], 1)
            elif "dotriacontahexen" in stringlike:
                return (["32:6"], 1)
            elif "dotriacontaheptaen" in stringlike:
                return (["32:7"], 1)
            elif "dotriacontahepten" in stringlike:
                return (["32:7"], 1)
            elif "tritriacontan" in stringlike:
                return (["33:0"], 1)
            elif "psyllic" in stringlike:
                return (["33:0"], 1)
            # no FA(33:1) !
            elif "tritriacontadien" in stringlike:
                return (["33:2"], 1)
            elif "tritriacontatrien" in stringlike:
                return (["33:3"], 1)
            elif "tritriacontatetraen" in stringlike:
                return (["33:4"], 1)
            elif "tritriacontapentaen" in stringlike:
                return (["33:5"], 1)
            elif "tritriacontapenten" in stringlike:
                return (["33:5"], 1)
            elif "tritriacontahmevaexaen" in stringlike:
                return (["33:6"], 1)
            elif "tritriacontahexen" in stringlike:
                return (["33:6"], 1)
            elif "tritriacontaheptaen" in stringlike:
                return (["33:7"], 1)
            elif "tritriacontahepten" in stringlike:
                return (["33:7"], 1)
            elif "tetratriacontan" in stringlike:
                return (["34:0"], 1)
            elif "gheddic" in stringlike: # Gheddic acid
                return (["34:0"], 1)
            elif "tetratriaconten" in stringlike:
                return (["34:1"], 1)
            elif "tetratriacontadien" in stringlike:
                return (["34:2"], 1)
            elif "tetratriacontatrien" in stringlike:
                return (["34:3"], 1)
            elif "tetratriacontatetraen" in stringlike:
                return (["34:4"], 1)
            elif "tetratriacontapentaen" in stringlike:
                return (["34:5"], 1)
            elif "tetratriacontapenten" in stringlike:
                return (["34:5"], 1)
            elif "tetratriacontahexaen" in stringlike:
                return (["34:6"], 1)
            elif "tetratriacontahexen" in stringlike:
                return (["34:6"], 1)
            elif "tetratriacontaheptaen" in stringlike:
                return (["34:7"], 1)
            elif "tetratriacontahepten" in stringlike:
                return (["34:7"], 1)
            elif "ceroplastic" in stringlike:
                return (["35:0"], 1)
            elif "pentatriacontan" in stringlike:
                return (["35:0"], 1)
            elif "hexatriacontylic" in stringlike:
                return (["36:0"], 1)
            elif "hexatriacontan" in stringlike:
                return (["36:0"], 1)
            elif "octatriacontan" in stringlike:
                return (["38:0"], 1)
            elif "hexatetracontan" in stringlike:
                return (["46:0"], 1)

            # special

            elif "caproic" in stringlike: # LAST
                return (["6:0"], 1)
            elif "octan" in stringlike: # LsteAST (must search other *octanoi)
                return (["8:0"], 1)
            elif "caprylic" in stringlike: # caprylic acid
                return (["8:0"], 1)
            elif "capryloyl" in stringlike: # caprylic acid
                return (["8:0"], 1)
            elif "octen" in stringlike: # LAST
                return (["8:1"], 1)
            elif "octaen" in stringlike: # LAST
                return (["8:1"], 1)
            elif "octadien" in stringlike: # LAST
                return (["8:2"], 1)
            elif "decan" in stringlike: # LAST (must search other *decan)
                return (["10:0"], 1)
            elif "decen" in stringlike: # LAST
                return (["10:1"], 1)
            elif "decadien" in stringlike: # LAST
                return (["10:2"], 1)
            elif "decatrien" in stringlike: # LAST
                return (["10:3"], 1)
            elif "decatetraen" in stringlike: # LAST
                return (["10:4"], 1)
            elif "stear" in stringlike: #stearic acid LAST
                return (["18:0"], 1)
            elif "olei" in stringlike: # C18:1 cis(n9) LAST
                return (["18:1"], 1)
            elif "olea" in stringlike: # C18:1 cis(n9) LAST
                return (["18:1"], 1)
            elif "oleo" in stringlike: # C18:1 cis(n9) LAST
                return (["18:1"], 1)
            elif "oley" in stringlike: # C18:1 cis(n9) LAST
                return (["18:1"], 1)
            #elif "oleic" in stringlike: # C18:1 cis(n9) LAST
            #    return (["18:1"], 1)
            #elif "oleyl" in stringlike: # C18:1 cis(n9) LAST
            #    return (["18:1"], 1)
            #elif "oleoyl" in stringlike: # C18:1 cis(n9) LAST
            #    return (["18:1"], 1)
            elif "elaid" in stringlike: # C18:1 trans(n9)
                return (["18:1"], 1)
            elif "linol" in stringlike: # linoleic acid LAST
                return (["18:2"], 1)
            elif "eicosan" in stringlike: # LAST
                return (["20:0"], 1)
            elif "docosan" in stringlike:
                return (["22:0"], 1) # LAST
            else:
                return ([], 0)


        lipid = preprocess(lipid)

        # Check if there are lipid residues in the lipid name.
        # If not, we still try to see whether it is a known, named fatty acid.
        
        lipid_residue_pattern = r"([\d]*:[\d]*)" # we don't care about _ or / separators
        residues = re.findall(lipid_residue_pattern, lipid)
        
        if len(residues) == 0: # [] no fatty acid residues
                               # trying to recognize residues in the name

            if lipid.endswith("ate"):
                compound_lipid = lipid.split(" ")
                
                if len(compound_lipid) != 2:
                    return ([], 0)
                else:
                    allresidues = []
                    for chunk in compound_lipid:
                        res, div = sloppy_lipid_identifier(chunk)
                        if div != 0:
                            allresidues.extend(res)
                
                if len(allresidues) != 0:
                    return (allresidues, 1)
                else:
                    return ([], 0)

            else:
                return sloppy_lipid_identifier(lipid)
        
        else:
            # _ means that all lipid residues are present. Example: TG(11:0_18:4_20:5)
            # / is used with mixed meaning. It both indicates uncertainty between
            # mass isomers, and it looks like it is used ambiguously with the same
            # meaning as _ within single entities in parentheses. Examples:
            #
            # TG(23:0/23:0/23:0) mass: 1101.06499
            # TG(23:0_23:0_24:0) mass: 1115.08064
            # PE(P-18:0/25:0)/PE(O-18:1/25:0)
            # TG(22:5/22:5/22:5) mass: 1028.78329
            # TG(22:5_22:5_22:5) mass: 1028.78329
            #
            # The number of mass isomers lies in the number of parenthesized
            # entities.
            parentheses_pattern = r"\(.*?\)"
            parenthesized_entities = re.findall(parentheses_pattern, lipid)

            dividend = len(parenthesized_entities)
            if dividend > 1:
                ambiguity = True
            else:
                ambiguity = False

            lipid_residue_pattern = r"([\d]*:[\d]*)" # we don't care about _ or / separators
            lipid_residues = re.findall(lipid_residue_pattern, lipid)

            if drop_ambiguous == True and ambiguity == True:
                return ([], 0)
            else:
                return (lipid_residues, dividend)       



def make_residues_table(dataframe, *, drop_ambiguous=False, name="residues_table",
                        replace_nan=0, cleanup=True, absolute_amount=False,
                        unwanted=["total", "fc", "tc"], liptype="refmet", **kwargs):

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

    liptype: Defaults to RefMet nomenclature. It can be set to None to use the alternative
        method of residue extraction, .residues(), based on another nomenclature that
        is used in some commercial lipidomics services.

    returns:
    ========

    pandas DataFrame

    """
    
    if not isinstance(dataframe, pd.core.frame.DataFrame):
        raise TypeError("Input table must be a pandas DataFrame")

    if liptype not in ("refmet", None):
        raise TypeError("Unrecognized lipids data type. Try 'refmet' or None")

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

            if liptype == "refmet":
                residues, coeff = lip.refmet_residues(drop_ambiguous=drop_ambiguous)
            else:
                residues, coeff = lip.residues(drop_ambiguous=drop_ambiguous)
            
            if len(residues) < 0:
                continue
            else:
                for residue in residues:
                    mm.setdefault(residue, 0)

                    try:
                        if absolute_amount:
                            mm[residue] += lip.molecules / coeff # exact number of residues
                        else:
                            mm[residue] += lip.amount / coeff # units of residue
                    except ZeroDivisionError:
                        print("Trying to use non-Refmet names? Try using liptype=None instead.")
                        raise

        masterlist.append(mm)


    dfinal = pd.concat([pd.DataFrame.from_dict(x, orient="index") for x in masterlist], axis=1)
    dfinal.columns = df.columns

    dfinal.name = name

    return dfinal
