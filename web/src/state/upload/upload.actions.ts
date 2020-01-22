import actionCreatorFactory from 'typescript-fsa'

import { FileType } from './upload.types'

export interface AddFilesPayload {
  files: Map<FileType, string>
}

export interface RemoveFilesPayload {
  files: Map<FileType, string>
}
export interface UploadError {
  error: Error
}

const action = actionCreatorFactory('UPLOAD')

export const addFiles = action<AddFilesPayload>('ADD_FILES')

export const removeFiles = action<RemoveFilesPayload>('REMOVE_FILES')
