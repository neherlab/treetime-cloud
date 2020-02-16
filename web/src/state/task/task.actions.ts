import actionCreatorFactory from 'typescript-fsa'

import { TaskId } from './task.types'

export interface GetTaskIdAsyncResult {
  taskId: TaskId
}

export interface PostTaskPayload {
  task: {
    taskId: TaskId
  }
}

export interface PostTaskResult {
  taskId: TaskId
}

export interface TaskError {
  error: Error
}

const action = actionCreatorFactory('TASK')

export const triggerGetTaskId = action('GET_TASK_ID_TRIGGER')

export const getTaskIdAsync = action.async<
  void,
  GetTaskIdAsyncResult,
  TaskError
>('GET_TASK_ID_ASYNC')

export const triggerPostTask = action<PostTaskPayload>('POST_TASK_TRIGGER')

export const asyncPostTask = action.async<
  PostTaskPayload,
  PostTaskResult,
  TaskError
>('POST_TASK_ASYNC')
