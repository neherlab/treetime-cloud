#!/usr/bin/env python
"""
This module defines functions used to facilitate the Beast run from other scripts
and to parse the Beast output results.
"""
import pandas
import xml.etree.ElementTree as XML
import utility_functions_general as gen_utils
import os, sys
import subprocess
from Bio import AlignIO, Phylo
import StringIO
from external_binaries import BEAST_BIN
import treetime

def read_beast_log(logfile, nearest_leaf_date, take_last_lines=500):
    """
    Reads the log file, produced by Beast. Note the column names are defined in
    the Beast template file from the resources. If another template used, make
    sure the  column names are the same.

    Args:
     - logfile(str): filename to be parsed

     - nearest_leaf_date(float): the date of the youngest leaf in the tree. This
     data is needed since Beast only outputs the height of the tree in years.
     Additional conversion to the absolute date is needed.

     - take_last_lines(int): Use the last lines from the log to estimate the convergence
     and stability of the simulations. If the log file has less lines, no output is
     possible. Hence, this is additional control of the Beast simulation finished.

    Returns:

     - LogFile parsed as pandas dataframe is there is enough data lines
     (specified by take_last_lines). None otherwise.
    """

    with open(logfile) as inlog:
        ss = inlog.readlines()
    if len(ss) < take_last_lines:
        return None
    df = pandas.read_csv(logfile, delimiter='\t', skiprows=3)
    df['treeModel.rootHeight']  = nearest_leaf_date  - df['treeModel.rootHeight']
    return df.iloc[-take_last_lines:]

def create_beast_xml(tree, aln, dates, log_file, template_file):
    """
    Take template XML configuration and create a valid Beast configuration.
    Basically, the function adds the initial tree, alignment and leaf dates to
    the corresponding sections of the template XML.

    Args:

     - tree(str or Biopython tree): initial tree for beast simulations. Tree
     object of path to the tree file

     - aln(str or Biopython multiple seq alignment object): alignment to set
     sequences to the tree leaves.

     - dates(dict): dictionary of the leaf dates.

     - log_file(str): path to save the BEAST configuration and log files. Only
     basename is needed. The following files will be produced after the successful
     BEAST run:
        - <log_file>.log.txt - the log of the BEAST operation, which represents
        the main results of the BEAST run. Can be opened with BEAST Tracer program
        - <log_file>.config.xml - the configuration file to run BEAST
        - <log_file>.trees.txt - tree samples taken every 1e6 iterations.

     - template_file(str):  path to the template XML to be used to produce
     config.

    Returns:

     - config as ElementTree.Xml object
    """

    def _set_taxa_dates(xml_root, tree, dates):

        def _create_taxon(name, date=None):

            """
            create XML structure, which holds taxon info: id, date, etc
            Args:

             - date: numeric date (e.g. 2013.7243) as float or string

             - name: name of the sequence/tree leaf. The name should match exactly

            Returns:
             - xml_taxon: Taxon data as XML structure ready to be plugged to the BEAST xml
            """
            xml_taxon = XML.Element('taxon')
            xml_taxon.attrib = {"id" : name}
            if date is not None:
                xml_date = XML.Element("date")
                xml_date.attrib = {"value": str(date), "direction" : "forwards", "units":"years"}
                xml_taxon.append(xml_date)
            return xml_taxon

        xml_taxa = xml_root.find('taxa')
        for leaf in tree.get_terminals():
            name = leaf.name
            if name not in dates:
                date=None
            else:
                date = dates[name]
            xml_taxa.append(_create_taxon(name, date))

    def _set_aln(xml_root, aln, tree=None):
        xml_aln = xml_root.find('alignment')

        if tree is not None:
            leaf_names = [k.name for k in tree.get_terminals()]
        else:
            leaf_names = None

        for seq in aln:
            if leaf_names is not None and seq.name not in leaf_names:
                continue

            xml_taxon = XML.Element("taxon")
            xml_taxon.attrib = {"idref" : seq.name}

            xml_seq = XML.Element("sequence")
            xml_seq.append(xml_taxon)
            xml_seq.text = str(seq.seq)

            xml_aln.append(xml_seq)

    def _set_newick(xml_root, tree):
        xml_nwk = xml_root.find('newick')
        st_io = StringIO.StringIO()
        Phylo.write(tree, st_io, 'newick')
        xml_nwk.text = st_io.getvalue()

    def _set_log_output(xml_root, log_file):

        xml_filelog = XML.Element("log")

        xml_filelog.attrib = {"id" : "filelog",
                "fileName" : log_file + ".log.txt",
                "overwrite" : "true",
                "logEvery": "10000"}

        posterior = XML.Element("posterior")
        posterior.attrib = {"idref" : "posterior"}
        xml_filelog.append(posterior)

        prior = XML.Element("prior")
        prior.attrib = {"idref" : "prior"}
        xml_filelog.append(prior)

        likelihood = XML.Element("likelihood")
        likelihood.attrib = {"idref" : "likelihood"}
        xml_filelog.append(likelihood)

        parameter = XML.Element("parameter")
        parameter.attrib = {"idref" : "treeModel.rootHeight"}
        xml_filelog.append(parameter)

        parameter = XML.Element("parameter")
        parameter.attrib = {"idref" : "constant.popSize"}
        xml_filelog.append(parameter)

        parameter = XML.Element("parameter")
        parameter.attrib = {"idref" : "CP1+2.kappa"}
        xml_filelog.append(parameter)

        parameter = XML.Element("parameter")
        parameter.attrib = {"idref" : "CP3.kappa"}
        xml_filelog.append(parameter)

        parameter = XML.Element("parameter")
        parameter.attrib = {"idref" : "CP1+2.frequencies"}
        xml_filelog.append(parameter)

        parameter = XML.Element("parameter")
        parameter.attrib = {"idref" : "CP3.frequencies"}
        xml_filelog.append(parameter)

        compoundParameter = XML.Element("compoundParameter")
        compoundParameter.attrib = {"idref" : "allMus"}
        xml_filelog.append(compoundParameter)

        parameter = XML.Element("parameter")
        parameter.attrib = {"idref" : "clock.rate"}
        xml_filelog.append(parameter)

        treeLikelihood = XML.Element("treeLikelihood")
        treeLikelihood.attrib = {"idref" : "CP1+2.treeLikelihood"}
        xml_filelog.append(treeLikelihood)

        treeLikelihood = XML.Element("treeLikelihood")
        treeLikelihood.attrib = {"idref" : "CP3.treeLikelihood"}
        xml_filelog.append(treeLikelihood)

        xml_root.append(xml_filelog)

        xml_logtree = XML.Element("logTree")
        xml_logtree.attrib = {"id" : "treeFileLog",
                    "logEvery" : "1000000",
                    "nexusFormat" : "true",
                    "fileName" : log_file + ".trees.txt",
                    "sortTranslationTable": "true"}

        treeModel = XML.Element("treeModel")
        treeModel.attrib = {"idref" : "treeModel"}
        xml_logtree.append(treeModel)

        posterior = XML.Element("posterior")
        posterior.attrib = {"idref" : "posterior"}
        xml_logtree.append(posterior)

        xml_root.append(xml_logtree)

        return

    # prepare input data
    if isinstance(tree, str):
        tree = Phylo.read(tree, 'newick')

    if isinstance(aln, str):
        aln = AlignIO.read(aln, 'fasta')

    # read template
    xml = XML.parse(template_file)
    xml_root = xml.getroot()

    # set data to the template
    _set_taxa_dates(xml_root, tree, dates)
    _set_aln(xml_root, aln, tree)
    _set_newick(xml_root, tree)

    _set_log_output(xml_root.find("mcmc"), log_file)

    return xml

def run_beast(tree, aln, dates, out_filename_prefix, template_file, log_post_process = None):
    """
    Run Beast for the specified tree, alignmentm, dates. It first prouces the
    Beast template using the specified data, and then calls Beast binary in a
    subprocess for the configuration produced.

    Args:

     - tree(str or Biopython tree): initial tree for beast simulations. Tree
     object of path to the tree file

     - aln(str or Biopython multiple seq alignment object): alignment to set
     sequences to the tree leaves.

     - dates(dict): dictionary of the leaf dates.

     - out_file_prefix(str): path to save the resulting configuration. Beast
     results and configuration will be saved to different files, which names will
     share the same specified prefix. Configuration will be '<prefix>.config.txt'
     the results will be saved under '<prefix>.log.txt' and '<prefix>.trees.txt'

     - template_file(str):  path to the template XML to be used to produce
     config.

    Returns:
     - None

    """

    config_filename = out_filename_prefix + ".config.xml"
    config_xml = create_beast_xml(tree, aln, dates, out_filename_prefix, template_file)
    config_xml.write(config_filename)
    call = ["java", "-jar", BEAST_BIN, "-beagle_off", "-overwrite",  config_filename]
    subprocess.call(call)

    #  process log, save the data to a pivot table
    if log_post_process is not None:
        print ("BEAST log post-processing...")
        # get log file
        log_file = out_filename_prefix + ".log.txt"
        # processing the log file using the external callable:
        log_post_process(log_file)


if __name__ == '__main__':
    pass
