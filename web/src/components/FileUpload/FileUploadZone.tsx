import path from 'path'

import React, { useCallback } from 'react'
import { connect } from 'react-redux'

import { useDropzone } from 'react-dropzone'

import { State } from '../../state/reducer'
import { selectTaskId } from '../../state/task/task.selectors'
import { TaskId } from '../../state/task/task.types'
import {
  triggerUploadFiles,
  UploadFilesPayload,
} from '../../state/upload/upload.actions'
import { selectFiles } from '../../state/upload/upload.selectors'
import { FileType } from '../../state/upload/upload.types'

/* Adds relevant files to a Map to be dispatched */
function reduceDroppedFiles(files: Map<FileType, File>, file: File) {
  const type = fileExtToType(path.extname(file.name))
  if (type) {
    files.set(type, file)
  }
  return files
}

/* Converts file extension to file FileType enum */
function fileExtToType(ext: string) {
  const extMap = new Map<string, FileType>(
    Object.entries({
      '.nwk': FileType.NWK,
      '.fasta': FileType.FASTA,
      '.csv': FileType.DATES,
      '.json': FileType.CONFIG,
    }),
  )
  return extMap.get(ext)
}

export interface FileUploadZoneProps {
  files: Map<FileType, File>
  taskId?: TaskId
  triggerFileUpload(payload: UploadFilesPayload): void
}

function FileUploadZone({
  files,
  taskId,
  triggerFileUpload,
}: FileUploadZoneProps) {
  const onDrop = useCallback(
    (droppedFiles: File[]) => {
      if (taskId) {
        const files = droppedFiles.reduce(reduceDroppedFiles, new Map())
        triggerFileUpload({ files, taskId })
      }
    },
    [taskId, triggerFileUpload],
  )

  const { getRootProps, getInputProps, isDragActive } = useDropzone({ onDrop })

  return (
    <div {...getRootProps()}>
      <p>{`TaskID: ${taskId}`}</p>
      <input type="file" {...getInputProps()} />
      {taskId &&
        (isDragActive ? (
          <p>{'Drop the files here ...'}</p>
        ) : (
          <p>{`Drag 'n' drop some files here, or click to select files`}</p>
        ))}
      <ul>
        {[...files.values()].map(({ name }: File) => (
          <li key={name}>{name}</li>
        ))}
      </ul>
    </div>
  )
}

const mapStateToProps = (state: State) => ({
  files: selectFiles(state),
  taskId: selectTaskId(state),
})

const mapDispatchToProps = { triggerFileUpload: triggerUploadFiles }

export default connect(mapStateToProps, mapDispatchToProps)(FileUploadZone)
