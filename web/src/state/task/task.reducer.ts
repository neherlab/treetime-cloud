import { reducerWithInitialState } from 'typescript-fsa-reducers'

import immerCase from '../util/fsaImmerReducer'

import { getTaskIdAsync } from './task.actions'

import { Task } from './task.types'

export interface TaskState {
  task?: Task
  pending: number
  error: Error | null
}

export const taskDefaultState: TaskState = {
  task: undefined,
  pending: 0, // number of async updates in-flight
  error: null,
}

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
