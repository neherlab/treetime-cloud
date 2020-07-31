import actionCreatorFactory from 'typescript-fsa'

import { TaskId } from '../task/task.types'
import { FileType } from './upload.types'

export interface RemoveFilesPayload {
  files: Map<FileType, File>
}

export interface UploadFilesPayload {
  files: Map<FileType, File>
  taskId: TaskId
}

export interface FileUploadAsyncResult {
  taskId: TaskId
}

export interface UploadError {
  error: Error
}

const action = actionCreatorFactory('Upload')

export const removeFiles = action<RemoveFilesPayload>('removeFiles')

export const triggerUploadFiles = action<UploadFilesPayload>('triggerUploadFiles')

export const uploadFilesAsync = action.async<UploadFilesPayload, FileUploadAsyncResult, UploadError>('uploadFilesAsync')
