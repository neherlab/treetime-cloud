import actionCreatorFactory from 'typescript-fsa'

import { FileType } from './upload.types'

export interface RemoveFilesPayload {
  files: Map<FileType, File>
}

export interface UploadFilesPayload {
  files: Map<FileType, File>
}

export interface FileUploadAsyncPayload {} // eslint-disable-line @typescript-eslint/no-empty-interface

export interface FileUploadAsyncResult {} // eslint-disable-line @typescript-eslint/no-empty-interface

export interface UploadError {
  error: Error
}

const action = actionCreatorFactory('UPLOAD')

export const removeFiles = action<RemoveFilesPayload>('REMOVE_FILES')

export const triggerUploadFiles = action<UploadFilesPayload>(
  'UPLOAD_FILES_TRIGGER',
)

export const uploadFilesAsync = action.async<
  UploadFilesPayload,
  FileUploadAsyncResult,
  UploadError
>('UPLOAD_FILES_ASYNC')
