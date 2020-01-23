import { State } from '../reducer'

export const selectTaskId = (state: State) => state.task.task?.id
