import os

from typing import Iterable, Optional

from zipfile import ZipFile


def make_zip(input_filepaths: Iterable[str],
             output_filepath: str,
             relative_to: Optional[str] = None):
  with ZipFile(output_filepath, 'w') as zipfile:
    for filepath in input_filepaths:
      dest_filepath = os.path.relpath(filepath,
                                      relative_to) if relative_to else None
      zipfile.write(filepath, dest_filepath)
