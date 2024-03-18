from Bio.Align import PairwiseAligner
import Bio.Seq



class get_exons:
    def __init__(self, ab1: dict[str,  list[int]]) -> None:
        self.ab1 = ab1


    def align_local(self, seqA: Bio.Seq.Seq, seqB: Bio.Seq.Seq) -> Bio.Align.PairwiseAlignments:
        """
        Given 2 sequences, it align them and return the start and end of the overlapping area.
        Matches: 2, Mismatches: -1, Gap open: -2, Gap extend: -0.1. trying to force the sequences tostay compact
        Args:
            seqA (Bio.Seq.Seq): sequence A
            seqB (Bio.Seq.Seq): sequence B

        Returns:
           tuple(int, int): start and end of the overlapping area
        """
        aligner = PairwiseAligner()
        aligner.mode = 'local'
        aligner.match_score = 2
        aligner.mismatch_score = -1
        aligner.open_gap_score = -10
        aligner.extend_gap_score = -0.1
        alignments = aligner.align(seqA, seqB)
        return alignments
    

    def get_indices(self, exon: str) -> dict[int,np.ndarray[int, int]]:
        """
        Given the exon number, it run the alignment and return the indexes for the 2 sequences:
        Args:
            exon (str): exon number

        Returns:
            dict[int, Bio.Align.PairwiseAlignments]: dictionary with the exon number as key and the alignment object as value
        """
        if exon != str:
            try:
                exon = str(exon)
            except ValueError:
                raise ValueError(f"Exon must be a string or an integer, not {type(exon)}")
        indices = self.align_local(seqA=self.ab1[exon + "F"]["seq"], seqB=self.ab1[exon +"R"]["seq"]).__getitem__(0).indices.T
        return {exon:indices}
    

    def check_index(self, index: int, strand= str) ->  tuple[str,  int]:
        """
        Given the index, it checks which bases and scores are assigned to that position
        Args:
            index (int): position to retrieve base and score
            strand (str): forward or reverse primers: (F or R)

        Returns:
            
        """
        if index == -1:
            seq_base = "-"
            seq_score = 0
        else:
            seq_base = self.ab1[index + strand]["seq"][index]
            seq_score = self.ab1[index + strand]["phred"][index]
        return seq_base, seq_score


    def get_base(self, indices = numpy.ndarray) -> dict[str,  str]:
        """
        Given the exon number and the indices, it return the base in that position with better score
        Args:
            indices (numpy.ndarray): array with the indices of the position for seq1 and seq2
            exon (str): exon number
        Returns:
            dict[str,  str]: dictionary with the exon number as key and the base as value
        """
        seq1_base, seq1_score = self.check_index(indices[0], "F")
        seq2_base, seq2_score = self.check_index(indices[1], "R")

        if seq1_base == seq2_base:
            return seq1_base
        else:
            if seq1_score > seq2_score:
                return seq1_base
            else:
                return seq2_base    


# TODO: finish exon sequence and test it
   def exon_sequence(self, indices=numpy.ndarray, exon=str) -> dict[str,  str]:
        """
        Given the exon number and the indices, it return the consensus sequence  based on the quality scores
        Args:
           indices (dict): dictionary with the exon number as key and the indexes as value

        Returns:
            dict[str,  str]: dictionary with the exon number as key and the sequence as value
        """
        exon_seq = self.ab1[exon + "F"]["seq"][indices[0]:indices[1]]
        exon_phred = self.ab1[exon + "F"]["phred"][indices[0]:indices[1]]
        return {exon:{"seq":exon_seq, "phred":exon_phred}}


if __name__=="__main__":
    pass