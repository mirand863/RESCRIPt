# ----------------------------------------------------------------------------
# Copyright (c) 2019-2023, QIIME 2 development team.
#
# Distributed under the terms of the Modified BSD License.
#
# The full license is in the file LICENSE, distributed with this software.
# ----------------------------------------------------------------------------


import qiime2
from qiime2.plugin.testing import TestPluginBase
from q2_types.feature_data import DNAIterator
from q2_types.multiplexed_sequences import FASTQIterator
from qiime2.plugins import rescript
from rescript.orient import _add_optional_parameters

import_data = qiime2.Artifact.import_data


class TestOrientSeqs(TestPluginBase):
    package = 'rescript.tests'

    def setUp(self):
        super().setUp()
        self.fasta_seqs = \
            import_data('FeatureData[Sequence]',
                        self.get_data_path('mixed-orientations.fasta'))
        self.fastq_seqs = \
            import_data('MultiplexedSingleEndBarcodeInSequence',
                        self.get_data_path('mixed-orientations.fastq.gz'))
        self.ref = import_data(
            'FeatureData[Sequence]', self.get_data_path('derep-test.fasta'))

        self.rc = \
            import_data('FeatureData[Sequence]',
                        self.get_data_path('mixed-orientations-rc.fasta')
                        )

    def test_reorient_fasta_default(self):
        # TODO: fix FASTA tests
        # this test checks that expected IDs AND reoriented seqs are returned
        reoriented, unmatched, = rescript.actions.orient_seqs(
            sequences=self.fasta_seqs, reference_sequences=self.ref)
        # all input seqs should pass, minus the two junk seqs. Those that pass
        # are copies of the ref seqs with some mismatches/gaps inserted
        # so the oriented seq ids should match the ref ids.
        reoriented_seqs = {seq.metadata['id']: str(seq)
                           for seq in reoriented.view(DNAIterator)}
        exp_reoriented_seqs = {seq.metadata['id']: str(seq)
                               for seq in self.ref.view(DNAIterator)}
        self.assertEqual(reoriented_seqs.keys(), exp_reoriented_seqs.keys())
        # the oriented seqs should also match the ref seqs, except for those
        # with some mismatches inserted...
        exp_reoriented_seqs['A1'] = (
            'AAAAAAGCAGCCGCGGTAATACGGAGGGTGCGAGCGTTGTCCGGAATCACTGGGCGTAAAGGGCG'
            'CGTAGGTGGCCCGTTAAGTGGCTGGTGAAATCCCGGGGCTCAACTCCGGGGCTGCCGGTCAGACT'
            'GGCGAGCTAGAGCACGGTAGGGGCAGATGGAATTCCCGGTGTAGCGGTGGAATGCGTAGATATCG'
            'GGAAGAATACCAGTGGCGAAGGCGTTCTGCTGGGCCGTTGCTGACACTGAGGCGCGACAGCGTGG'
            'GGAGCAAACAGGATTAGATACCCTGGTAGTC')
        exp_reoriented_seqs['A2'] = (
            'AGCAGCCGCGGTAATACGGAGGGTGCGAGCGTTGTCCGGAATCACTGGGCGTAAAGGGCGCGTAG'
            'GTGGCCCGTTAAGTGGCTGGTGAAATCCCGGGGCTCAACTCCGGGGCTGCCGGTCAGACTGGCGA'
            'GGAAGAGCACGGTAGGGGCAGATGGAATTCCCGGTGTAGCGGTGGAATGCGTAGATATCGGGAAG'
            'AATACCAGTGGCGAAGGCGTTCTGCTGGGCCGTTGCTGACACTGAGGCGCGACAGCGTGGGGAGC'
            'AAACAGGATTAGATACCCTGGTAGTC')
        exp_reoriented_seqs['A3'] = exp_reoriented_seqs['A3'][:-21] + 'T' * 21
        self.assertEqual(
            set(reoriented_seqs.values()), set(exp_reoriented_seqs.values()))
        # make sure the junk sequences don't match anything and are discarded
        unmatched_ids = {
            seq.metadata['id'] for seq in unmatched.view(DNAIterator)}
        exp_unmatched_ids = {'JUNK', 'MOREJUNK'}
        self.assertEqual(unmatched_ids, exp_unmatched_ids)

    def test_reorient_fasta_no_ref(self):
        reoriented, unmatched = rescript.actions.orient_seqs(
            sequences=self.fasta_seqs, reference_sequences=None,
            )
        unmatched_ids = {seq.metadata['id']
                         for seq in unmatched.view(DNAIterator)}
        self.assertEqual(unmatched_ids, set([]))
        exp_seqs = [seq for seq in self.rc.view(DNAIterator)]
        test_seqs = [seq for seq in reoriented.view(DNAIterator)]
        for exp, test in zip(*(exp_seqs, test_seqs)):
            self.assertEqual(str(exp), str(test))
            self.assertEqual(exp.metadata['id'], test.metadata['id'])

    def test_reorient_fastq_default(self):
        # this test checks that expected IDs AND reoriented seqs are returned
        reoriented, unmatched, = rescript.actions.orient_seqs(
            sequences=self.fastq_seqs, reference_sequences=self.ref)
        # all input seqs should pass, minus the two junk seqs. Those that pass
        # are copies of the ref seqs with some mismatches/gaps inserted
        # so the oriented seq ids should match the ref ids.
        reoriented_seqs = {seq.metadata['id']: str(seq)
                           for seq in reoriented.view(FASTQIterator)}
        exp_reoriented_seqs = {seq.metadata['id']: str(seq)
                               for seq in self.ref.view(DNAIterator)}
        self.assertEqual(reoriented_seqs.keys(), exp_reoriented_seqs.keys())
        # the oriented seqs should also match the ref seqs, except for those
        # with some mismatches inserted...
        exp_reoriented_seqs['A1'] = (
            'AAAAAAGCAGCCGCGGTAATACGGAGGGTGCGAGCGTTGTCCGGAATCACTGGGCGTAAAGGGCG'
            'CGTAGGTGGCCCGTTAAGTGGCTGGTGAAATCCCGGGGCTCAACTCCGGGGCTGCCGGTCAGACT'
            'GGCGAGCTAGAGCACGGTAGGGGCAGATGGAATTCCCGGTGTAGCGGTGGAATGCGTAGATATCG'
            'GGAAGAATACCAGTGGCGAAGGCGTTCTGCTGGGCCGTTGCTGACACTGAGGCGCGACAGCGTGG'
            'GGAGCAAACAGGATTAGATACCCTGGTAGTC')
        exp_reoriented_seqs['A2'] = (
            'AGCAGCCGCGGTAATACGGAGGGTGCGAGCGTTGTCCGGAATCACTGGGCGTAAAGGGCGCGTAG'
            'GTGGCCCGTTAAGTGGCTGGTGAAATCCCGGGGCTCAACTCCGGGGCTGCCGGTCAGACTGGCGA'
            'GGAAGAGCACGGTAGGGGCAGATGGAATTCCCGGTGTAGCGGTGGAATGCGTAGATATCGGGAAG'
            'AATACCAGTGGCGAAGGCGTTCTGCTGGGCCGTTGCTGACACTGAGGCGCGACAGCGTGGGGAGC'
            'AAACAGGATTAGATACCCTGGTAGTC')
        exp_reoriented_seqs['A3'] = exp_reoriented_seqs['A3'][:-21] + 'T' * 21
        self.assertEqual(
            set(reoriented_seqs.values()), set(exp_reoriented_seqs.values()))
        # make sure the junk sequences don't match anything and are discarded
        unmatched_ids = {
            seq.metadata['id'] for seq in unmatched.view(FASTQIterator)}
        exp_unmatched_ids = {'JUNK', 'MOREJUNK'}
        self.assertEqual(unmatched_ids, exp_unmatched_ids)


# Needed to edit the file:
# /home/fabio/mambaforge/envs/rescript/lib/python3.8/site-packages/q2_types/multiplexed_sequences/_transformer.py
# # ----------------------------------------------------------------------------
# # Copyright (c) 2016-2023, QIIME 2 development team.
# #
# # Distributed under the terms of the Modified BSD License.
# #
# # The full license is in the file LICENSE, distributed with this software.
# # ----------------------------------------------------------------------------
# import collections.abc
#
# import skbio.io
#
# from . import (MultiplexedFastaQualDirFmt,
#                MultiplexedSingleEndBarcodeInSequenceDirFmt)
#
# from ..plugin_setup import plugin
#
#
# @plugin.register_transformer
# def _1(df: MultiplexedFastaQualDirFmt) -> \
#      MultiplexedSingleEndBarcodeInSequenceDirFmt:
#     seqs = open(df.sequences.path_maker())
#     qual = open(df.quality.path_maker())
#
#     result = MultiplexedSingleEndBarcodeInSequenceDirFmt()
#
#     with open(result.path / 'forward.fastq.gz', 'wb') as fh:
#         for seq in skbio.io.read(seqs, qual=qual, format='fasta',
#                                  verify=False):
#             seq.write(fh, format='fastq', variant='illumina1.8',
#                       compression='gzip')
#
#     return result
#
#
# class NucleicAcidIterator(collections.abc.Iterable):
#     def __init__(self, generator):
#         self.generator = generator
#
#     def __iter__(self):
#         yield from self.generator
#
#
# class FASTQIterator(NucleicAcidIterator):
#     pass
#
#
# # DNA Transformers
# @plugin.register_transformer
# def _2(ff: MultiplexedSingleEndBarcodeInSequenceDirFmt) -> FASTQIterator:
#     generator = _read_from_fastq(str(ff) + "/forward.fastq.gz", skbio.DNA)
#     return FASTQIterator(generator)
#
#
# def _read_from_fastq(path, constructor=skbio.DNA, lowercase=False):
#     return skbio.read(
#     	path, format='fastq',
#     	constructor=constructor,
#     	lowercase=lowercase,
#     	phred_offset=33,
#     )

    # def test_reorient_fastq_no_ref(self):
    #     reoriented, unmatched = rescript.actions.orient_seqs(
    #         sequences=self.fasta_seqs, reference_sequences=None,
    #         )
    #     unmatched_ids = {seq.metadata['id']
    #                      for seq in unmatched.view(DNAIterator)}
    #     self.assertEqual(unmatched_ids, set([]))
    #     exp_seqs = [seq for seq in self.rc.view(DNAIterator)]
    #     test_seqs = [seq for seq in reoriented.view(DNAIterator)]
    #     for exp, test in zip(*(exp_seqs, test_seqs)):
    #         self.assertEqual(str(exp), str(test))
    #         self.assertEqual(exp.metadata['id'], test.metadata['id'])

    def test_add_optional_parameters(self):
        expected = [
            "--dbmask", "dust",
            "--relabel", "new_id",
            "--relabel_keep",
            "--relabel_md5",
            "--relabel_self",
            "--relabel_sha1",
            "--sizein",
            "--sizeout",
        ]
        result = []
        _add_optional_parameters(
            result,
            dbmask="dust",
            relabel="new_id",
            relabel_keep=True,
            relabel_md5=True,
            relabel_self=True,
            relabel_sha1=True,
            sizein=True,
            sizeout=True,
        )
        self.assertEqual(result, expected)
