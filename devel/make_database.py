import requests
import pandas as pd
import re

refmet_url = "https://www.metabolomicsworkbench.org/databases/refmet/refmet_latest.xlsx"
refmet_file = "refmet.xlsx"


def fetch_database(url):
    req = requests.get(url, allow_redirects=False)
    
    with open(refmet_file, "wb") as o:
        o.write(req.content)


def process_database(filename, wanted_cols = ["name", "exactmass", "main_class"],
                     to_translate=["Fatty esters", "Fatty amides", "Fatty acids"],
    ):
    """This processes an Excel RefMet list of compounds, by removing
    all residues that do not contain any fatty acid or similar residues.
    
    The resulting table is then saved as a tab-delimited .csv text.

    It also saves another table of named compounds that need to be
    manually translated. The intended use of this is to periodically
    check new entries against updated RefMet versions.
    """

    database = pd.read_excel(filename)

    # == step 1 ==

    df = database.copy()
    df = df[wanted_cols]

    #choosing lipid names that seemingly contain fatty acids residues
    new_index = df["name"].str.extract(r"(\d*:\d*)")
    new_index = new_index.dropna()
    df = df.reindex(new_index.index)
    
    #removing some troublesome classes by hand
    df = df[df["main_class"] != "Disaccharides"]
    df = df[df["main_class"] != "Steroids"]
    
    df.to_csv(filename[:-5]+".csv", sep="\t")

    # == step 2 ==

    database = database[database["main_class"].isin(to_translate)]
    database.to_csv(filename[:-5]+"_to_translate.csv", sep="\t")


def extract_fatty_classes(filename):
    """Reads fatty classes from RefMet list of compounds, returns set of classes
    """
    df = pd.read_excel(filename)
    return set([x for x in df["main_class"] if "fatt" in x.lower()])




if __name__ == "__main__":

    print("Downloading latest RefMet database..")
    fetch_database(refmet_url)
    print("Processing RefMet database..")
    process_database(refmet_file)