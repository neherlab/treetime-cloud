import pandas
import xml.etree.ElementTree as XML
import utility_functions_general as gen_utils
from external_binaries import  *
import os, sys
import subprocess
from Bio import AlignIO, Phylo
import StringIO
def read_beast_log(logfile, nearest_leaf_date, take_last_lines=500):
    """
    """

    with open(logfile) as inlog:
        ss = inlog.readlines()
    if len(ss) < 100:
        return None
    df = pandas.read_csv(logfile, delimiter='\t', skiprows=3)
    df['treeModel.rootHeight']  = nearest_leaf_date  - df['treeModel.rootHeight']
    return df.iloc[-take_last_lines:]

def create_beast_xml(tree, aln, dates, log_file, template_file):

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

def run_beast(tree, aln, dates, out_filename_prefix, template_file):

    config_filename = out_filename_prefix + ".config.xml"
    config_xml = create_beast_xml(tree, aln, dates, out_filename_prefix, template_file)
    config_xml.write(config_filename)
    call = ["java", "-jar", BEAST_BIN, "-beagle_off", "-overwrite",  config_filename]
    subprocess.call(call)

