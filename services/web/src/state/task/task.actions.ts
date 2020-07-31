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

const action = actionCreatorFactory('Task')

export const getTaskIdTrigger = action('getTaskIdTrigger')

export const getTaskIdAsync = action.async<void, GetTaskIdAsyncResult, TaskError>('getTaskIdAsync')

export const postTaskTrigger = action<PostTaskPayload>('postTaskTrigger')

export const postTaskAsync = action.async<PostTaskPayload, PostTaskResult, TaskError>('postTaskAsync')
