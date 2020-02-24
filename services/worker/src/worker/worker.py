import json
import os

from typing import Any, Dict, NamedTuple, Optional, Tuple, Union

from Bio import Phylo, AlignIO, Align, Seq, SeqRecord
import numpy as np
import pandas

from treetime import TreeTime
from treetime.utils import numeric_date


class TreetimeInputFilepaths(NamedTuple):
  DATES: str
  FASTA: str
  NWK: str


class TreetimeOutputFilepaths(NamedTuple):
  NWK: Optional[str]
  NWK_GENERATED: str
  NEX: str
  FASTA: str
  GTR: str
  META_CSV: str
  # MOLECULAR_CLOCK: str
  TREE_JSON: str
  # LIKELIHOODS_JSON: str
  CONFIG_JSON: str


class TreetimeRelaxedClockConfig(NamedTuple):
  slack: float
  coupling: float


class TreetimeConfig(NamedTuple):
  input_filenames: TreetimeInputFilepaths
  output_filenames: TreetimeOutputFilepaths
  output_zip_filename: str
  generate_tree: bool
  gtr: str
  do_marginal: bool
  resolve_polytomies: bool
  max_iter: int
  relaxed_clock: Union[bool, TreetimeRelaxedClockConfig]
  root: Optional[str] = "best"
  slope: Optional[float] = None
  coalescent_prior: Optional[float] = None


def isfloat(value):
  try:
    float(value)
    return True
  except ValueError:
    return False


def write_json(data: Dict, filepath: str):
  with open(filepath, "w") as f:
    json.dump(data, f, indent=2)


def read_metadata_from_file(infile: str, log: str) -> Tuple[Any, Any]:
  """
    @brief      Reads a metadata from file or handle.
    @param      self    The object
    @param      infile  The input file or file handle
    """
  try:
    # read the metadata file into pandas dataframe.
    df = pandas.read_csv(infile, index_col=0, sep=r"\s*,\s*", engine="python")
    # check the metadata has strain names in the first column
    # look for the column containing sampling dates
    # We assume that the dates might be given either in human-readable format
    # (e.g. ISO dates), or be already converted to the numeric format.
    if "name" not in df.index.name.lower():
      print(
          "Cannot read metadata: first column should contain the names of the strains",
          file=log,
      )
      return
    potential_date_columns = []
    potential_numdate_columns = []

    # Scan the dataframe columns and find ones which likely to store the dates
    for ci, col in enumerate(df.columns):
      d = df.iloc[0, ci]
      if type(d) == str and d[0] in ['"', "'"] and d[-1] in ['"', "'"]:
        for i, tmp_d in enumerate(df.iloc[:, ci]):
          df.iloc[i, ci] = tmp_d.strip(d[0])
      if "date" in col.lower():
        if isfloat(df.iloc[0, ci]):
          # avoid date parsing when can be parsed as float
          potential_numdate_columns.append((ci, col))
        else:
          # otherwise add as potential date column
          potential_date_columns.append((ci, col))

    # if a potential numeric date column was found, use it (use the first, if there are more than one)
    if len(potential_numdate_columns) >= 1:
      name = potential_numdate_columns[0][1]
      # Use this column as numdate_given
      dates = df[name].to_dict()
      for k, val in dates.items():
        try:
          dates[k] = float(val)
        except ValueError:
          dates[k] = None

    elif len(potential_date_columns) >= 1:
      # try to parse the csv file with dates in the idx column:
      idx = potential_date_columns[0][0]
      name = potential_date_columns[0][1]
      # NOTE as the 0th column is the index, we should parse the dates
      # for the column idx + 1
      df = pandas.read_csv(
          infile,
          index_col=0,
          sep=r"\s*,\s*",
          parse_dates=[1 + idx],
          engine="python",
      )
      dates = {k: numeric_date(df.loc[k, name]) for k in df.index}
      df.loc[:, name] = map(lambda x: str(x.date()), df.loc[:, name])
    else:
      print(
          "Metadata file has no column which looks like a sampling date!",
          file=log,
      )
    metadata = df.to_dict(orient="index")
    for k, val in metadata.items():
      if type(k) == str and k[0] in ["'", '"'] and k[-1] in ["'", '"']:
        metadata[k.strip(k[0])] = val
        dates[k.strip(k[0])] = dates[k]
    return dates, metadata
  except:
    print("Cannot read the metadata file. Exception caught!", file=log)
    raise
    return {}, {}


def generate_tree(input_aln_filename: str, output_tree_filename: str) -> None:

  call = [
      "fasttree", "-nt", "-quiet", input_aln_filename, " > ",
      output_tree_filename
  ]

  cmdline = " ".join(call)

  res = os.system(cmdline)
  if res != 0:
    raise RuntimeError(
        "fasttree: Exception caught while building the phylogenetic tree")


def layout(tt: TreeTime):
  """Add clade, xvalue, yvalue, mutation and trunk attributes to all nodes in tree"""
  tt.tree.root.branch_length = 0.001
  clade = 0
  yvalue = 0
  for node in tt.tree.find_clades(order="preorder"):
    # set mutations
    if node.up is not None:
      node.muts = ", ".join([
          node.up.sequence[p] + str(p) + node.sequence[p]
          for p in np.where(node.up.sequence != node.sequence)[0]
      ])

    # set sequences
    node.strseq = "".join(node.sequence)

    # set clade No
    node.clade = clade
    clade += 1
    if node.up is not None:  # try:
      # Set xValue, tValue, yValue
      node.xvalue = node.up.xvalue + node.mutation_length
      node.tvalue = node.numdate - tt.tree.root.numdate
    else:
      node.xvalue = 0.0
      node.tvalue = 0.0
    if node.is_terminal():
      node.yvalue = yvalue
      yvalue += 1
    # check numdate
    if not hasattr(node, "numdate"):
      node.numdate = 0.0
  for node in tt.tree.get_nonterminals(order="postorder"):
    node.yvalue = np.mean([x.yvalue for x in node.clades])


def node_metadata(node, metadata, config):

  if (hasattr(node, "branch_length_interpolator") and
      node.branch_length_interpolator is not None):
    gamma = node.branch_length_interpolator.gamma
  else:
    gamma = 1.0

  meta = []
  # if node has user-provided metadata, append it
  if node.name in metadata:
    node_meta = metadata[node.name]
    meta += [{"name": k, "value": node_meta[k]} for k in node_meta]

  # append numdates to the metadata
  if hasattr(node, "numdate"):
    meta.append({"name": "numdate", "value": node.numdate})

  # append deviation of the branch length from average
  stretch = np.min([
      node.branch_length / (node.mutation_length / gamma + 0.000001),
      2.0,
  ])
  meta.append({"name": "Branch length stretch", "value": stretch})

  if config.relaxed_clock:
    # append mutation rate deviation from average
    meta.append({"name": "Local substitution rate", "value": gamma})
    # else:
    #     meta.append({"name": "Relaxed mutation rate", "value":1.0})

  return meta


def node_to_json(node, meta, config):
  """
  convert node data to dictionary
  """

  tree_json = {}
  str_attr = ["clade", "strain", "date", "muts", "strseq"]
  num_attr = ["xvalue", "yvalue", "tvalue", "numdate"]

  if hasattr(node, "name"):
    tree_json["strain"] = node.name
    tree_json["name"] = node.name

  for prop in str_attr:
    if hasattr(node, prop):
      tree_json[prop] = node.__getattribute__(prop)
  for prop in num_attr:
    if hasattr(node, prop):
      try:
        tree_json[prop] = round(node.__getattribute__(prop), 5)
      except:
        logger(
            "cannot round:",
            node.__getattribute__(prop),
            "assigned as is",
            1,
            warn=True,
        )
        tree_json[prop] = node.__getattribute__(prop)

  if node.clades:  # node is internal
    tree_json["children"] = []
    for ch in node.clades:
      tree_json["children"].append(node_to_json(ch, meta, config))

  tree_json["metadata"] = node_metadata(node, meta, config)

  return tree_json


def likelihoods_to_json(tt: TreeTime):
  """
  Save the likelihoods fro the node positions in JSON format
  """
  out_dic = {}
  for node in tt.tree.find_clades():
    x, y = _distribution_to_human_readable(tt, node.marginal_pos_LH)
    arr = [{"x": f, "y": b} for f, b in zip(x, y)]
    out_dic[node.name] = arr

  return out_dic


def _distribution_to_human_readable(tt: TreeTime, dist: Any):

  if dist.is_delta:
    date = tt.date2dist.to_numdate(dist.peak_pos)
    return (
        [date - 0.5, date - 1e-10, date, date + 1e-10, date + 0.5],
        [0, 0, 1.0, 0, 0],
    )

  peak_pos = dist.peak_pos
  fwhm = dist.fwhm
  raw_x = dist.x[(dist.x > peak_pos - 3 * fwhm) &
                 (dist.x < peak_pos + 3 * fwhm)]
  dates_x = np.array(map(tt.date2dist.to_numdate, raw_x))
  y = dist.prob_relative(raw_x)
  return dates_x, y


def decorate(tt: TreeTime):
  for node in tt.tree.find_clades():
    if node.up is None:
      continue

    node.confidence = None

    if len(node.mutations) > 0:
      node.comment = (
          '&mutations="' +
          ",".join([a + str(pos) + d for (a, pos, d) in node.mutations]) + '"')


def save_alignment(tt: TreeTime, config: TreetimeConfig):
  records = [
      SeqRecord.SeqRecord(
          Seq.Seq("".join(n.sequence)),
          id=n.name,
          name=n.name,
          description="",
      ) for n in tt.tree.find_clades()
  ]

  aln = Align.MultipleSeqAlignment(records)

  with open(config.output_filenames.FASTA, "w") as ofile:
    AlignIO.write(aln, ofile, "fasta")


def save_metadata_to_csv(tt: TreeTime, meta: Any, config: TreetimeConfig):
  meta = {
      node: node_metadata(node, meta, config) for node in tt.tree.find_clades()
  }

  rows_dic = {
      key.name: {e["name"]: e["value"] for e in meta[key]
                } for key in meta.keys()
  }

  df = pandas.DataFrame(rows_dic).T
  df.to_csv(config.output_filenames.META_CSV)


def save_molecular_clock_to_csv(tt: TreeTime, config: TreetimeConfig):
  # save molecular clock in normal format
  mclock = np.array([
      (tip.dist2root, tip.numdate_given)
      for tip in tt.tree.get_terminals()
      if hasattr(tip, "dist2root") and hasattr(tip, "numdate_given")
  ])
  np.savetxt(
      config.output_filenames.MOLECULAR_CLOCK,
      mclock,
      delimiter=",",
      header="Distance_to_root,Sampling_date",
  )


def save_gtr(tt: TreeTime, config: TreetimeConfig):
  with open(config.output_filenames.GTR, "w") as outf:
    outf.write(str(tt.gtr))


def run_treetime(config: TreetimeConfig):

  if config.generate_tree:
    generate_tree(config.input_filenames.FASTA,
                  config.output_filenames.NWK_GENERATED)

  input_nwk = config.output_filenames.NWK_GENERATED if config.generate_tree else config.input_filenames.NWK

  tree = Phylo.read(input_nwk, "newick")
  aln = AlignIO.read(config.input_filenames.FASTA, "fasta")
  dates, meta = read_metadata_from_file(config.input_filenames.DATES,
                                        '.cache/log.txt')

  # TODO: in case of config.gtr == "infer" shall we default to "jc" here ?
  # TreeTime.run() has additional "infer_gtr" option.
  # Why same thing in two places?
  tt = TreeTime(dates=dates, tree=tree, aln=aln, gtr="jc")

  tt.run(root=config.root,
         infer_gtr=config.gtr == "infer",
         relaxed_clock=config.relaxed_clock,
         resolve_polytomies=config.resolve_polytomies,
         max_iter=config.max_iter,
         Tc=config.coalescent_prior,
         fixed_slope=config.slope,
         do_marginal=config.do_marginal)

  config_json = config._asdict()
  write_json(config_json, config.output_filenames.CONFIG_JSON)

  layout(tt)

  tree_json = node_to_json(tt.tree.root, meta, config)
  write_json(tree_json, config.output_filenames.TREE_JSON)

  # likelihoods_json = likelihoods_to_json(tt)
  # write_json(likelihoods_json, config.output_filenames.LIKELIHOODS_JSON)

  Phylo.write(tt.tree, config.output_filenames.NWK, "newick")

  decorate(tt)

  Phylo.write(tt.tree, config.output_filenames.NEX, "nexus")

  save_alignment(tt, config)

  save_metadata_to_csv(tt, meta, config)

  # TODO: Do we need this? It was commented out in the original treetime-web.
  # save_molecular_clock_to_csv(tt, config)

  save_gtr(tt, config)
