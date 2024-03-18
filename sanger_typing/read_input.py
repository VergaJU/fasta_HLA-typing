from Bio import SeqIO
import Bio.Seq
import os
import re

# start class_import_ab1, it takes a file_path as input. __init__ take the inputs,
# then the function load_ab1 load the ab1 file and return the sequence and the quality score
class import_ab1:
    def __init__(self, folder_path: str) -> None:
        self.folder_path = folder_path

    def list_ab1_files(self) -> list[str]:
        """
        List all the .ab1 files in the folder
        Args:
            folder_path (str): /path/to/folder

        Returns:
            list[str]: list of .ab1 files in the folder
        """
        ab1_files = []
        for file in os.listdir(self.folder_path):
            if file.endswith(".ab1"):
                ab1_files.append(file)
        return ab1_files
    
    def load_ab1(self, filename: str) -> dict[Bio.Seq.Seq,  list[int]]:
        """
        Load the sequence and the quality score from the ab1 file.
        Args:
            filename (str): filename.ab1

        Returns:
            dict[str,  list[int]]: dictionay with the sequence and the quality score
        """

        with open(filename, "rb") as handle:
            for record in SeqIO.parse(handle, "abi"):
                sequence = record.seq
                quality_score = record.letter_annotations["phred_quality"]
        return {"seq":sequence, "phred":quality_score} 


    def check_reverse(self, exon_direction: str, ab1: dict) -> dict[Bio.Seq.Seq,  list[int]]:
        """
        Reverse complement the sequence if the exon_direction contains an R
        Args:
            exon_direction (str): exon_direction
            ab1 (dict): dictionary with sequence and quality score

        Returns:
            dict[str,  list[int]]: dictionary with the sequence and the quality score (reverse complemented if needed)
        """
        if exon_direction.endswith('r') or exon_direction.endswith('R'):
            ab1["seq"] = ab1["seq"].reverse_complement()
            ab1["phred"] = ab1["phred"][::-1]
        else:
            pass
        return ab1  
    
    def get_exon_direction(self, filename: str) -> str:
        """
        Given a file name with this pattern:
        sample_name_gene_exon_direction.ab1
        return exon direction

        Args:
            filename (str): filename.ab1

        Returns:
            str: exon_direction
        """        
        pattern = re.compile(r"\d[fFrR]")
        list_strings = filename.split("_")
        list_strings[-1] = list_strings[-1].split(".")[0]
        matches = [string for string in list_strings if re.match(pattern, string)]
        if matches:
            return matches[0]
        else:
            return ValueError("Filename does not match the pattern: \
                                sample_exondirection.ab1 \nexample: sample_A_2F.ab1")
    
    
    def load_ab1_files(self) -> dict[str,  dict[Bio.Seq.Seq, list[int]]]:
        ab1_files = self.list_ab1_files()
        ab1_dict = {}
        for file in ab1_files:
            exon_direction = self.get_exon_direction(filename=file)
            file_path = self.folder_path + file
            ab1 = self.load_ab1(filename=file_path)
            ab1 = self.check_reverse(exon_direction=exon_direction, ab1=ab1)
            ab1_dict[exon_direction] = ab1
        return ab1_dict
    


if __name__=="__main__":
    pass