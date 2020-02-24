import os

from typing import Callable, List

from filestore import FileStore

from worker import (TreetimeConfig, TreetimeInputFilepaths,
                    TreetimeOutputFilepaths, make_zip, run_treetime)

from .types import Task


class TaskConsumer(Callable):
  """ Consumes Task objects """

  def __init__(self, file_store: FileStore):
    self._file_store = file_store

  def __call__(self, task: Task) -> None:
    """
    Implements Callable interface, so that we can pass this object as a callback
    """
    self._consume_task(task)

  def _consume_task(self, task: Task) -> None:
    """
    Consumes a task by running a useful job on it
    """
    task_id = task.task_id
    input_filenames = task.input_filenames

    data_dir = ".data"
    task_dir = os.path.join(data_dir, task_id)
    input_dir = os.path.join(task_dir, "input")
    output_dir = os.path.join(task_dir, "output")
    zip_dir = os.path.join(task_dir, "zip")

    inputs = TreetimeInputFilepaths(
        NWK=os.path.join(input_dir, input_filenames["NWK"]),
        FASTA=os.path.join(input_dir, input_filenames["FASTA"]),
        DATES=os.path.join(input_dir, input_filenames["DATES"]),
    )

    outputs = TreetimeOutputFilepaths(
        NWK=os.path.join(output_dir, f"{task_id}.nwk"),
        NWK_GENERATED=os.path.join(output_dir, f"{task_id}.generated.nwk"),
        NEX=os.path.join(output_dir, f"{task_id}.nex"),
        FASTA=os.path.join(output_dir, f"{task_id}.fasta"),
        GTR=os.path.join(output_dir, f"{task_id}.gtr.txt"),
        META_CSV=os.path.join(output_dir, f"{task_id}.metadata.csv"),
        # MOLECULAR_CLOCK=os.path.join(output_dir, f"{task_id}.molecular_clock.csv"),
        TREE_JSON=os.path.join(output_dir, f"{task_id}.tree.json"),
        # LIKELIHOODS_JSON=os.path.join(output_dir, f"{task_id}.likelihoods.json"),
        CONFIG_JSON=os.path.join(output_dir, f"{task_id}.config.json"),
    )

    output_zip_filename = os.path.join(zip_dir, f"{task_id}.zip")

    config = TreetimeConfig(
        input_filenames=inputs,
        output_filenames=outputs,
        output_zip_filename=output_zip_filename,
        generate_tree=True,
        gtr="jc",
        root="best",
        do_marginal=False,
        resolve_polytomies=False,
        slope=None,
        coalescent_prior=None,
        max_iter=2,
        relaxed_clock=False,
    )

    os.makedirs(input_dir, exist_ok=True)
    os.makedirs(output_dir, exist_ok=True)
    os.makedirs(zip_dir, exist_ok=True)
    for local_filepath in config.input_filenames:
      self._file_store.download_file(task_id, 'input', local_filepath)

    run_treetime(config)

    nwk_filepath = (
        outputs.NWK_GENERATED
        if config.generate_tree else outputs.NWK_GENERATED)
    final_output_filepaths: List[str] = [
        nwk_filepath,
        outputs.NEX,
        outputs.FASTA,
        outputs.GTR,
        outputs.META_CSV,
        outputs.TREE_JSON,
        outputs.CONFIG_JSON,
    ]

    nwk_input_filepath: List[str] = [] if config.generate_tree else [inputs.NWK]
    final_zip_filepaths = final_output_filepaths + nwk_input_filepath + [
        inputs.DATES, inputs.FASTA
    ]

    for local_filepath in final_output_filepaths:
      self._file_store.upload_file(task_id, 'output', local_filepath)

    make_zip(
        final_zip_filepaths, config.output_zip_filename, relative_to=data_dir)
    self._file_store.upload_file(task_id, 'zip', config.output_zip_filename)
