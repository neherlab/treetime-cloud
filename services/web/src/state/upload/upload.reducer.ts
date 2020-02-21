import { reducerWithInitialState } from 'typescript-fsa-reducers'

import immerCase from '../util/fsaImmerReducer'

import { removeFiles, triggerUploadFiles } from './upload.actions'
import { FileType } from './upload.types'

export interface UploadState {
  files: Map<FileType, File>
  pending: number
  error: Error | null
}

export const uploadDefaultState: UploadState = {
  files: new Map<FileType, File>(),
  pending: 0, // number of async updates in-flight
  error: null,
}

export const uploadReducer = reducerWithInitialState(uploadDefaultState)
  .withHandling(
    immerCase(triggerUploadFiles, (draft, payload) => {
      // Merge new files into state
      draft.files = new Map([...draft.files, ...payload.files])
    }),
  )

  .withHandling(
    immerCase(removeFiles, (draft, payload) => {
      if (!draft.files) {
        return
      }

      // Remove files by type from state
      payload.files.forEach((_, fileType) => draft.files.delete(fileType))
    }),
  )
