import glob
import xml.etree.ElementTree as ET

import numpy as np
import pandas as pd
import seaborn as sns


class RNASequences:
    def __init__(self, data_dir="Data/"):
        filenames = glob.glob("*.txt", root_dir=data_dir)
        dfs = []

        for filename in filenames:
            df = pd.read_csv(
                data_dir + filename, sep="\t", names=[filename[:10]], skiprows=1
            ).T
            dfs.append(df)

        self.__rna_counts = pd.concat(dfs)

        self.__annotations = pd.DataFrame(
            columns=("Subject ID", "Sample Group", "CNS Subregion")
        )
        tree = ET.parse(open("Data/GSE124439_family.xml"))
        root = tree.getroot()
        namespace = {"ns": "http://www.ncbi.nlm.nih.gov/geo/info/MINiML"}

        for sample in root.findall("ns:Sample", namespace):
            sample_id = sample.attrib["iid"]
            for channel in sample.iterfind(".//ns:Channel", namespace):
                self.__annotations.loc[sample_id] = [
                    channel[4].text.strip(),
                    channel[5].text.strip(),
                    channel[3].text.strip(),
                ]

        self.__annotations = self.__annotations
        self.__check_annotations()

    def __check_annotations(self):
        assert self.__rna_counts.index.difference(self.__annotations.index).empty

    # def _ipython_display_(self):
    #     display(self.__rna_counts)
    #     display(self.__annotations)

    def get_counts(self):
        return self.__rna_counts

    def get_annotations(self):
        return self.__annotations

    def get_count(self, item):
        if isinstance(item, str):
            return self.__rna_counts[item]
        else:
            return self.__rna_counts.at[item[0], item[1]]

    def get_annotation(self, item):
        if isinstance(item, str):
            return self.__annotations[item]
        else:
            return self.__annotations.at[item[0], item[1]]

    def get_sample(self, sample=None):
        match sample:
            case "ALS":
                return self.__annotations.loc[
                    self.__annotations["Sample Group"] == "ALS Spectrum MND"
                ]
            case "Control":
                return self.__annotations.loc[
                    self.__annotations["Sample Group"] == "Non-Neurological Control"
                ]
            case "Other":
                return self.__annotations.loc[
                    self.__annotations["Sample Group"] == "Other Neurological Disorders"
                ]
            case _:
                return self.__annotations

    def get_sample_count(self, sample=None):
        return self.__rna_counts.loc[self.get_sample(sample).index]

    def mean(self, sample=None):
        return (
            self.get_sample_count(sample).mean().to_frame().rename(columns={0: "Means"})
        )

    def median(self, sample=None):
        return (
            self.get_sample_count(sample)
            .median()
            .to_frame()
            .rename(columns={0: "Medians"})
        )

    def std(self, sample=None):
        return (
            self.get_sample_count(sample)
            .std()
            .to_frame()
            .rename(columns={0: "Standard Deviations"})
        )
