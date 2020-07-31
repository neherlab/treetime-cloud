import { reducerWithInitialState } from 'typescript-fsa-reducers'

import immerCase from 'src/state/util/fsaImmerReducer'

import { removeFiles, triggerUploadFiles } from './upload.actions'
import { uploadDefaultState } from './upload.state'

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
