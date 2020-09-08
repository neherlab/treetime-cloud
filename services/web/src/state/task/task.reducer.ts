import { reducerWithInitialState } from 'typescript-fsa-reducers'

import immerCase from 'src/state/util/fsaImmerReducer'

import { getTaskIdAsync } from './task.actions'
import { taskDefaultState } from './task.state'

export const taskReducer = reducerWithInitialState(taskDefaultState)
  .withHandling(
    immerCase(getTaskIdAsync.started, (draft, payload) => {
      draft.pending += 1
    }),
  )

  .withHandling(
    immerCase(getTaskIdAsync.done, (draft, payload) => {
      draft.pending -= 1
      draft.task = { id: payload.result.taskId }
    }),
  )

  .withHandling(
    immerCase(getTaskIdAsync.failed, (draft, payload) => {
      draft.pending -= 1
      draft.error = payload.error.error
      draft.task = undefined
    }),
  )
