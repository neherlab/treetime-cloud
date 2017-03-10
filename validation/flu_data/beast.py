import numpy as np
from Bio import AlignIO, Phylo
from Bio.Align import  MultipleSeqAlignment
import random
import subprocess
import datetime
import os, copy
import matplotlib.pyplot as plt
from scipy.stats import linregress
from collections import Counter
import xml.etree.ElementTree as XML
import StringIO
import pandas

plt.ion()
plt.show()

def date_from_seq_name(name):
    def str2date_time(instr):
        """
        Convert input string to datetime object.
        Args:
         - instr (str): input string. Accepts one of the formats:
         {MM.DD.YYYY, MM.YYYY, MM/DD/YYYY, MM/YYYY, YYYY}.

        Returns:
         - date (datetime.datetime): parsed date object. If the parsing failed,
         None is returned
        """

        instr = instr.replace('/', '.')
        # import ipdb; ipdb.set_trace()
        try:
            date = datetime.datetime.strptime(instr, "%m.%d.%Y")
        except ValueError:
            date = None
        if date is not None:
            return date

        try:
            date = datetime.datetime.strptime(instr, "%m.%Y")
        except ValueError:
            date = None

        if date is not None:
            return date

        try:
            date = datetime.datetime.strptime(instr, "%Y")
        except ValueError:
            date = None

        return date

    date = str2date_time(name.split('|')[2].strip())

    return date.year + (date - datetime.datetime(date.year, 1, 1)).days / 365.25

def remove_polytomies(tree):
    """
    Scan tree, and remove the polytomies (if any) by random merging
    Returns: tree without polytomies
    """
    for clade in tree.find_clades():
        #if not hasattr(clade, "name") or clade.name is None:
        #    clade.name = "None"
        if len(clade.clades) < 3:
            continue
        clades = clade.clades
        while len(clades) > 2:
            c1 = clades.pop()
            c2 = clades.pop()
            new_node = Phylo.BaseTree.Clade()
            new_node.name="None"
            new_node.branch_length = 1e-5
            new_node.clades = [c1,c2]
            clades.append(new_node)
        clade.clades = clades

    return tree

def dates_from_flu_tree(tree):

    if isinstance(tree, str):
        tree = Phylo.read(tree, 'newick')

    dates = {k.name:date_from_seq_name(k.name) for k in tree.get_terminals()
                if date_from_seq_name(k.name) is not None}
    return dates

def internal_regress(myTree):
    resarr = []
    for node in myTree.tree.get_nonterminals():
        try:
            resarr.append((node.numdate, node.dist2root))
        except:
            continue
    resarr = np.array(resarr)
    if resarr.shape[0] == 0:
        return 0.
    else:
        return linregress(resarr[:, 0], resarr[:, 1]).rvalue**2

def create_beast_xml(tree, aln, dates, log_file, template_file='../beast/template.xml'):

    def _set_taxa_dates(xml_root, tree, dates):

        def _create_taxon(name, date):

            """
            create XML structure, which holds taxon info: id, date, etc
            Args:

             - date: numeric date (e.g. 2013.7243) as float or string

             - name: name of the sequence/tree leaf. The name should match exactly

            Returns:
             - xml_taxon: Taxon data as XML structure ready to be plugged to the BEAST xml
            """

            xml_date = XML.Element("date")
            xml_date.attrib = {"value": str(date), "direction" : "forward", "units":"years"}

            xml_taxon = XML.Element('taxon')
            xml_taxon.attrib = {"id" : name}
            xml_taxon.append(xml_date)

            return xml_taxon

        xml_taxa = xml_root.find('taxa')
        for leaf in tree.get_terminals():
            name = leaf.name
            if name not in dates:
                continue
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
                "fileName" : log_file + ".log.xml",
                "overwrite" : "true",
                "logEvery": "1000"}

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

        #parameter = XML.Element("parameter")
        #parameter.attrib = {"idref" : "kappa"}
        #xml_filelog.append(parameter)

        parameter = XML.Element("parameter")
        parameter.attrib = {"idref" : "JC69.frequencies"}
        xml_filelog.append(parameter)

        #parameter = XML.Element("parameter")
        #parameter.attrib = {"idref" : "siteModel.alpha"}
        #xml_filelog.append(parameter)

        parameter = XML.Element("parameter")
        parameter.attrib = {"idref" : "ucld.mean"}
        xml_filelog.append(parameter)

        parameter = XML.Element("parameter")
        parameter.attrib = {"idref" : "ucld.stdev"}
        xml_filelog.append(parameter)

        rateStatistic = XML.Element("rateStatistic")
        rateStatistic.attrib = {"idref" : "meanRate"}
        xml_filelog.append(rateStatistic)

        rateStatistic = XML.Element("rateStatistic")
        rateStatistic.attrib = {"idref" : "coefficientOfVariation"}
        xml_filelog.append(rateStatistic)

        rateCovarianceStatistic = XML.Element("rateCovarianceStatistic")
        rateCovarianceStatistic.attrib = {"idref" : "covariance"}
        xml_filelog.append(rateCovarianceStatistic)

        treeLikelihood = XML.Element("treeLikelihood")
        treeLikelihood.attrib = {"idref" : "treeLikelihood"}
        xml_filelog.append(treeLikelihood)

        coalescentLikelihood = XML.Element("coalescentLikelihood")
        coalescentLikelihood.attrib = {"idref" : "coalescent"}
        xml_filelog.append(coalescentLikelihood)

        tmrcaStatistic = XML.Element("tmrcaStatistic")
        tmrcaStatistic.attrib = {"idref" : "Tmrca_stat"}
        xml_filelog.append(tmrcaStatistic)

        xml_root.append(xml_filelog)

        xml_logtree = XML.Element("logTree")
        xml_logtree.attrib = {"id" : "treeFileLog",
                    "logEvery" : "10000",
                    "nexusFormat" : "true",
                    "fileName" : log_file + ".trees.txt",
                    "sortTranslationTable": "true"}


        treeModel = XML.Element("treeModel")
        treeModel.attrib = {"idref" : "treeModel"}
        xml_logtree.append(treeModel)

        posterior = XML.Element("posterior")
        posterior.attrib = {"idref" : "posterior"}
        xml_logtree.append(posterior)

        discretizedBranchRates = XML.Element("discretizedBranchRates")
        discretizedBranchRates.attrib = {"idref" :"branchRates"}
        trait = XML.Element("trait")
        trait.attrib = {"name" : "rate",
                        "tag" : "rate"}
        trait.append(discretizedBranchRates)
        xml_logtree.append(trait)

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

def run_beast(subtree_file, out_dir, aln="./H3N2_HA_1980_2015_NA.fasta"):

    BEAST_BIN = "/ebio/ag-neher/share/programs/bundles/BEASTv1.8.4/lib/beast.jar"

    print ("Running beast for tree: " + subtree_file)

    if not os.path.exists(out_dir):
        os.makedirs(out_dir)

    tt = Phylo.read(subtree_file, 'newick')
    aln = AlignIO.read(aln, 'fasta')
    dates = dates_from_flu_tree(tt)

    #import ipdb; ipdb.set_trace();
    basename = os.path.join(out_dir, os.path.split(subtree_file)[-1][:-4])
    config = create_beast_xml(tt, aln, dates, basename)

    config_name = basename + ".config.xml"
    print (basename, config_name)
    config.write(config_name)

    call = ["java", "-jar", BEAST_BIN, "-beagle_off", "-overwrite",  config_name]
    subprocess.call(call)

if __name__ == '__main__':

    pass




