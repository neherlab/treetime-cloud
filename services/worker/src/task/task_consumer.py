from filestore import FileStore

from .types import Task


class TaskConsumer:
  """ Consumes Task objects """

  def __init__(self, file_store: FileStore):
    self._file_store = file_store

  def __call__(self, task: Task) -> None:
    """
    Implements call interface, so that we can pass this object as a callblack
    """
    self._consume_task(task)

  def _consume_task(self, task: Task) -> None:
    """
    Consumes a task by running a useful job on it
    """
    task_id = task.task_id
    input_filenames = task.input_filenames

    for _, filename in input_filenames.items():
      self._file_store.download_input_file(task_id, filename)

    # TODO: run the actual algorithm
    # root = ".data/"
    # cfg = tree_time_config.treetime_webconfig
    # ttw = TreeTimeWeb(root, cfg)
    # ttw.run()

    # TODO: upload output files, not inputs
    for _, filename in input_filenames.items():
      self._file_store.upload_output_file(task_id, filename)
