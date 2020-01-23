import actionCreatorFactory from 'typescript-fsa'

import { TaskId } from './task.types'

export type GetTaskIdPayload = undefined

export interface GetTaskIdAsyncResult {
  taskId: TaskId
}

export interface TaskError {
  error: Error
}

const action = actionCreatorFactory('UPLOAD')

export const triggerGetTaskId = action('GET_TASK_ID/TRIGGER')

export const getTaskIdAsync = action.async<
  void,
  GetTaskIdAsyncResult,
  TaskError
>('GET_TASK_ID_ASYNC')
