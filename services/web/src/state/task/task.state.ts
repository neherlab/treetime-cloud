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
