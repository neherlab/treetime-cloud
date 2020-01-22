import { State } from '../reducer'

export const selectFiles = (state: State) => state.upload.files
