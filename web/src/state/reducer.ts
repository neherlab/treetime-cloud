import { connectRouter, RouterState } from 'connected-react-router'
import { History } from 'history'
import { combineReducers } from 'redux'

import { errorReducer, ErrorState } from './error/error.reducer'
import { uploadReducer, UploadState } from './upload/upload.reducer'

const rootReducer = (history: History) =>
  combineReducers({
    upload: uploadReducer,
    error: errorReducer,
    router: connectRouter(history),
  })

export interface State {
  upload: UploadState
  error: ErrorState
  router: RouterState
}

export default rootReducer
