import React, { useCallback } from 'react'

import { connect } from 'react-redux'
import { Button } from 'reactstrap'
import styled from 'styled-components'
import { useTranslation } from 'react-i18next'

import { Uploader } from 'src/components/Main/Uploader'

import { FileType } from 'src/state/upload/upload.types'
import type { TaskId } from 'src/state/task/task.types'
import { triggerUploadFiles, UploadFilesPayload } from 'src/state/upload/upload.actions'
import { PostTaskPayload, postTaskTrigger } from 'src/state/task/task.actions'
import { State } from 'src/state/reducer'
import { selectFiles } from 'src/state/upload/upload.selectors'
import { selectTaskId } from 'src/state/task/task.selectors'
import path from 'path'

const ButtonRun = styled(Button)`
  width: 200px;
  height: 80px;
  margin-left: auto;
`

/* Adds relevant files to a Map to be dispatched */
export function reduceDroppedFiles(files: Map<FileType, File>, file: File) {
  const type = fileExtToType(path.extname(file.name))
  if (type) {
    files.set(type, file)
  }
  return files
}

/* Converts file extension to FileType enum */
export function fileExtToType(ext: string) {
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

export interface MainPageProps {
  files: Map<FileType, File>
  taskId?: TaskId
  triggerUploadFiles(payload: UploadFilesPayload): void
  triggerPostTask(payload: PostTaskPayload): void
}

const mapStateToProps = (state: State) => ({
  files: selectFiles(state),
  taskId: selectTaskId(state),
})

const mapDispatchToProps = {
  triggerUploadFiles,
  triggerPostTask: postTaskTrigger,
}

export const MainPage = connect(mapStateToProps, mapDispatchToProps)(MainPageDisconnected)

export function MainPageDisconnected({ files, taskId, triggerUploadFiles, triggerPostTask }: MainPageProps) {
  const { t } = useTranslation()

  const onUpload = useCallback(
    (droppedFiles: File[]) => {
      if (taskId) {
        const files = droppedFiles.reduce(reduceDroppedFiles, new Map())
        triggerUploadFiles({ files, taskId })
      }
    },
    [taskId, triggerUploadFiles],
  )

  const canRun =
    Boolean(taskId) &&
    Boolean(files.get(FileType.DATES)) &&
    Boolean(files.get(FileType.FASTA)) &&
    Boolean(files.get(FileType.NWK))

  return (
    <div>
      <Uploader onUpload={onUpload} />
      <div>
        <p>{`TaskID: ${taskId ?? ''}`}</p>
      </div>
      <ul>
        {[...files.values()].map(({ name }: File) => (
          <li key={name}>{name}</li>
        ))}
      </ul>
      <div className="d-flex">
        <div className="ml-auto">
          <ButtonRun type="button" disabled={!canRun} onClick={() => taskId && triggerPostTask({ task: { taskId } })}>
            {t('Run')}
          </ButtonRun>
        </div>
      </div>
    </div>
  )
}
