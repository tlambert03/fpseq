FPseq
=====

This is a small package used mostly for the FPbase project, and contains functions and methods for dealing with fluorescent protein sequences, alignments, mutations and HGVS-formatted mutation strings.

The `fpseq.FPseq` class pretty much just a subclass of the [scikit-bio](https://github.com/biocore/scikit-bio) [Protein](https://github.com/biocore/scikit-bio/blob/master/skbio/sequence/_protein.py#L17) class; however, the skbio code is duplicated here rather than imported, because I was having issues with numpy version dependencies with skbio that were causing problems for my automatic deployments on Heroku, so I cherry picked just what I needed, and the `FPSeq` class is the result.

The `align.ParasailAlignment` class is just a thin wrapper around the results of a [parasail-python](https://github.com/jeffdaily/parasail-python) alignment object, and us usually instantiated as the result of the `align.align_seqs` function.  It provides methods for converting the CIGAR from a global alignment of two protein sequences into a `MutationSet`.

`mutations.Mutation` is an object representing a single mutation operation (substitution, deletion, insertion, extension, etc...) in HGVS format and `mutations.MutationSet` is a sets of `Mutations`, such as would be required to mutate one fluorescent protein into another.  `MutationSets` can be applied to `FPSeqs` or sequence strings to generate a new `FPseq` sequence (with error checking to make sure that the `MutationSet` is actually valid given the parent sequence to which it is being applied).  `MutationSets` can be calculated from two protein sequences using, for example, `FPSeq.mutations_to()`.  MutationSets have many of the same methods as python sets (add, remove, intersection, union, contains, etc...)

`fpseq.from_fpbase` is a convenience function for retrieving a fluorescent protein sequence from FPbase, using the "slug" of the protein (i.e. the name of the protein that appears in the URL at FPbase after www.fpbase.org/protein/).

### Example Usage

```python
In [1]: from fpseq import from_fpbase

In [2]: avGFP = from_fpbase('avgfp')  # retrieve sequence from FPbase.org

In [3]: avGFP # sequences from FPbase.org
Out[3]: 
Protein
MSKGEELFTG VVPILVELDG DVNGHKFSVS GEGEGDATYG KLTLKFICTT
GKLPVPWPTL VTTFSYGVQC FSRYPDHMKQ HDFFKSAMPE GYVQERTIFF
KDDGNYKTRA EVKFEGDTLV NRIELKGIDF KEDGNILGHK LEYNYNSHNV
YIMADKQKNG IKVNFKIRHN IEDGSVQLAD HYQQNTPIGD GPVLLPDNHY
LSTQSALSKD PNEKRDHMVL LEFVTAAGIT HGMDELYK

In [4]: EGFP = from_fpbase('egfp')

In [6]: avGFP.mutations_to(EGFP)  # calculate HGVS mutation string
Out[6]: <MutationSet: M1_S2insV/F64L/S65T/H231L>

In [7]: mEGFP = from_fpbase('megfp')

In [8]: EGFP.mutations_to(mEGFP)
Out[8]: <MutationSet: A207K>

# A207K does not match the literature, because of V1a...
# use reference parameter to enforce position numbering relative to avGFP 
In [9]: EGFP.mutations_to(mEGFP, reference=avGFP)
Out[9]: <MutationSet: A206K>

In [10]: mCherry = from_fpbase('mcherry')

# attempt to apply the ‘mCherry2’ mutation from Shen et al. (2017)
In [11]: newseq = mCherry.mutate('K92N/K138C/K139R/S147T/N196D/T202L')
SequenceMismatch: Mutation K138C does not align with the parent seq: PSD>G<PVM.. But a match was found 5 positions away: K97N/K143C/K144R/S152T/N201D/T207L

# use correct_offset to apply a shift to the mutation set, if a match is found
In [12]: newseq, offset = mCherry.mutate('K92N/K138C/K139R/S147T/N196D/T202L', correct_offset=True)
UserWarning: An offset of 5 amino acids was detected between the sequence and the mutation set, and automatically corrected

In [13]: newseq == from_fpbase('mcherry2') # sequence equivalence checks
Out[13]: True
```