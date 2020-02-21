import path from 'path'

import React, { useCallback } from 'react'
import { connect } from 'react-redux'

import { useDropzone } from 'react-dropzone'

import { State } from '../../state/reducer'
import { PostTaskPayload, triggerPostTask } from '../../state/task/task.actions'
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

/* Converts file extension to FileType enum */
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
  triggerUploadFiles(payload: UploadFilesPayload): void
  triggerPostTask(payload: PostTaskPayload): void
}

function FileUploadZone({
  files,
  taskId,
  triggerUploadFiles,
  triggerPostTask,
}: FileUploadZoneProps) {
  const onDrop = useCallback(
    (droppedFiles: File[]) => {
      if (taskId) {
        const files = droppedFiles.reduce(reduceDroppedFiles, new Map())
        triggerUploadFiles({ files, taskId })
      }
    },
    [taskId, triggerUploadFiles],
  )

  const { getRootProps, getInputProps, isDragActive } = useDropzone({ onDrop })

  const canRun =
    Boolean(taskId) &&
    Boolean(files.get(FileType.DATES)) &&
    Boolean(files.get(FileType.FASTA)) &&
    Boolean(files.get(FileType.NWK))

  return (
    <div>
      <div>
        <p>{`TaskID: ${taskId}`}</p>
      </div>
      <div {...getRootProps()}>
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

      <div>
        <button
          type="submit"
          disabled={!canRun}
          onClick={() => taskId && triggerPostTask({ task: { taskId } })}
        >
          Run
        </button>
      </div>
    </div>
  )
}

const mapStateToProps = (state: State) => ({
  files: selectFiles(state),
  taskId: selectTaskId(state),
})

const mapDispatchToProps = {
  triggerUploadFiles,
  triggerPostTask,
}

export default connect(mapStateToProps, mapDispatchToProps)(FileUploadZone)
