import unittest
from liputils import Lipid, make_residues_table
import pandas as pd
import liputils

class LiputilsGeneralTests(unittest.TestCase):

    # liputils.Lipid - RefMet

    def test_refmet_decode_simple(self):
        self.l = Lipid("TG(18:4_20:4_27:0)")
        self.assertEqual(self.l.refmet_residues(), (["18:4", "20:4", "27:0"], 1))
        self.assertEqual(self.l.refmet_class(), "TG")
        self.assertEqual(self.l.refmet_fullclass(), "Triradylglycerols")

    def test_refmet_decode_complex(self):
        self.l = Lipid("PS(P-16:1/17:0)/PS(O-16:2/17:0)")
        self.assertEqual(self.l.refmet_residues(), (['16:1', '17:0', '16:2', '17:0'], 2))
        self.assertEqual(self.l.refmet_class(), "PS")
        self.assertEqual(self.l.refmet_fullclass(), "Glycerophosphoserines")

    def test_refmet_decode_composite(self):
        self.l = Lipid("Oleyl arachidonate")
        self.assertEqual(self.l.refmet_residues(), (["18:1", "20:4"], 1))
        self.assertEqual(self.l.refmet_class(), "Oleyl arachidonate")
        self.assertEqual(self.l.refmet_fullclass(), "unknown lipid type")

    def test_refmet_decode_nolipid(self):
        self.l = Lipid("Uranium phosphate")
        self.assertEqual(self.l.refmet_residues(), ([], 0))
        self.assertEqual(self.l.refmet_class(), "Uranium phosphate")
        self.assertEqual(self.l.refmet_fullclass(), "unknown lipid type")

    def test_refmet_decode_nolipid(self):
        self.l = Lipid("Uranium phosphate")
        self.assertEqual(self.l.refmet_residues(), ([], 0))
        self.assertEqual(self.l.refmet_class(), "Uranium phosphate")
        self.assertEqual(self.l.refmet_fullclass(), "unknown lipid type")

        self.assertEqual(self.l.residues(), ([], 1)) # intended
        self.assertEqual(self.l.lipid_class(), "unknown lipid type")

    def test_decode_zero_length_string(self):
        self.l = Lipid("")
        self.assertEqual(self.l.refmet_residues(), ([], 0))
        self.assertEqual(self.l.refmet_class(), "unknown lipid type")
        self.assertEqual(self.l.refmet_fullclass(), "unknown lipid type")

     # liputils.Lipid - other nomenclature

    def test_decode_simple(self):
        self.l = Lipid("CE 12:0")
        self.assertEqual(self.l.residues(), (['12:0'], 1))
        self.assertEqual(self.l.refmet_residues(), (['12:0'], 0)) # intended
        self.assertEqual(self.l.lipid_class(), "CE")
        self.assertEqual(self.l.refmet_fullclass(), "unknown lipid type")

    def test_decode_complex(self):
        self.l = Lipid("TAG 52:4 total (16:0/18:1/18:3)(16:0/18:2/18:2)")
        self.assertEqual(self.l.residues(), (['16:0', '18:1', '18:3', '16:0', '18:2', '18:2'], 2))
        self.assertEqual(
            self.l.refmet_residues(), (['52:4', '16:0', '18:1', '18:3', '16:0', '18:2', '18:2'], 2)
            ) # intended
        self.assertEqual(self.l.lipid_class(), "TAG")
        self.assertEqual(self.l.refmet_class(), "TAG 52:4 total")
        self.assertEqual(self.l.refmet_fullclass(), "unknown lipid type")

    # liputils.make_residues_table

    def test_refmet_tablemaking(self):
        RM_input = pd.read_csv("refmet_compliant_table.csv", sep="\t", index_col=0)
        RM_output = pd.read_csv("refmet_compliant_table_residues.csv", sep="\t", index_col=0)
        RM_test = make_residues_table(RM_input)

        for x,y in zip(RM_output.index, RM_test.index):
            self.assertEqual(x, y)
        for x,y in zip(RM_output.columns, RM_test.columns):
            self.assertEqual(x, y)
        for x,y in zip(RM_output[[RM_output.columns[0]]], RM_test[[RM_test.columns[0]]]):
            self.assertEqual(x, y)
        for x,y in zip(RM_output[[RM_output.columns[1]]], RM_test[[RM_test.columns[1]]]):
            self.assertEqual(x, y)

    def test_generic_type_tablemaking(self):
        RM_input = pd.read_csv("generic_lipid_table.csv", sep="\t", index_col=0)
        RM_output = pd.read_csv("generic_lipid_table_residues.csv", sep="\t", index_col=0)
        RM_test = make_residues_table(RM_input, liptype=None)

        for x,y in zip(RM_output.index, RM_test.index):
            self.assertEqual(x, y)
        for x,y in zip(RM_output.columns, RM_test.columns):
            self.assertEqual(x, y)
        for x,y in zip(RM_output[[RM_output.columns[0]]], RM_test[[RM_test.columns[0]]]):
            self.assertEqual(x, y)
        for x,y in zip(RM_output[[RM_output.columns[1]]], RM_test[[RM_test.columns[1]]]):
            self.assertEqual(x, y)


if __name__ == "__main__":
    print(f"Testing liputils version {liputils.__version__}")
    print(f"{liputils.__file__}")
    unittest.main()