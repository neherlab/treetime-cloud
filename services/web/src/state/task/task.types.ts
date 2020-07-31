import { Id } from 'src/helpers/types'

export type TaskId = Id<'TaskId'>

export interface Task {
  id: TaskId
}
